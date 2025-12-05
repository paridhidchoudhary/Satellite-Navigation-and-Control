% =======================================================================
%         NMPC with Receding Horizon for Spacecraft Control
% =======================================================================

import casadi.*

% ----------------------------- PARAMETERS ------------------------------
T_horizon = 480;       % [s] prediction horizon
N = 160;               % prediction steps
dt = T_horizon / N;    % 3 s time step

T_sim = 1200;          % [s] total simulation time (longer than horizon)
N_sim = T_sim / dt;    % simulation steps

% Setup NMPC problem (same as before)
fdyn  = f_dyn();
fcost = f_cost();
nx = 21; nu = 8;

% Build NLP and create a separate integrator function
X = SX.sym('X', nx, N+1);
U = SX.sym('U', nu, N);

J = 0; g = {}; lbg = []; ubg = [];
Wr = 1e-2 * eye(nu);

for k = 1:N
    J = J + fcost(X(:,k), U(:,k));

    J = J + 100 * (U(:,k)' * U(:,k));
    
    if k > 1
        du = U(:,k) - U(:,k-1);
        J = J + du' * Wr * du;
    end
    
    % RK-4 integration
    xk = X(:,k); uk = U(:,k);
    k1 = fdyn(xk, uk);
    k2 = fdyn(xk + 0.5*dt*k1, uk);
    k3 = fdyn(xk + 0.5*dt*k2, uk);
    k4 = fdyn(xk + dt*k3, uk);
    x_next = xk + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    
    g{end+1} = X(:,k+1) - x_next;
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];
end

% Terminal cost
P = 1e1 * eye(6);
rhoT = X(13:15,end);
rhodT = X(16:18,end);
J = J + [rhoT; rhodT]' * P * [rhoT; rhodT];

% Create solver
w = [X(:); U(:)];
g = vertcat(g{:});
nlp = struct('x', w, 'f', J, 'g', g);
opts = struct('ipopt', struct('print_level', 0, 'max_iter', 500));
solver = nlpsol('solver', 'ipopt', nlp, opts);

% Create a separate integrator function for simulation
x_sym = SX.sym('x', nx);
u_sym = SX.sym('u', nu);
k1 = fdyn(x_sym, u_sym);
k2 = fdyn(x_sym + 0.5*dt*k1, u_sym);
k3 = fdyn(x_sym + 0.5*dt*k2, u_sym);
k4 = fdyn(x_sym + dt*k3, u_sym);
x_next_sym = x_sym + dt/6*(k1 + 2*k2 + 2*k3 + k4);
integrator = Function('integrator', {x_sym, u_sym}, {x_next_sym});

% Bounds
u_max = 15;
lbx = -inf(size(w)); ubx = inf(size(w));
u_start = nx*(N+1) + 1;
lbx(u_start:end) = 0;
ubx(u_start:end) = u_max;

% ----------------------- RECEDING HORIZON SIMULATION -------------------
% Initial condition
x0 = zeros(nx,1);
x0([1 5 9]) = 1;
x0(13:15) = [-1.1394; 113.0847; -7.6204];
x0(16:18) = [0.0081; 0.0039; 1.3490];

% Storage for simulation results
X_sim = zeros(nx, N_sim+1);
U_sim = zeros(nu, N_sim);
X_sim(:,1) = x0;

% Initial guess
X_guess = repmat(x0, 1, N+1);
U_guess = 0.5 + 0.5*rand(nu, N);

fprintf('Running receding horizon NMPC...\n');
for i = 1:N_sim
    if mod(i,50) == 0
        fprintf('Step %d/%d\n', i, N_sim);
    end
    
    % Set initial condition constraint
    x_current = X_sim(:,i);
    
    % Update initial condition in bounds
    lbx(1:nx) = x_current;
    ubx(1:nx) = x_current;
    
    % Solve NMPC
    sol = solver('x0', [X_guess(:); U_guess(:)], ...
                 'lbx', lbx, 'ubx', ubx, ...
                 'lbg', lbg, 'ubg', ubg);
    
    % Extract solution
    w_opt = full(sol.x);
    X_opt = reshape(w_opt(1:nx*(N+1)), nx, N+1);
    U_opt = reshape(w_opt(nx*(N+1)+1:end), nu, N);
    
    % Apply first control input
    U_sim(:,i) = U_opt(:,1);
    
    % Simulate one step forward using the integrator function
    X_sim(:,i+1) = full(integrator(x_current, U_opt(:,1)));
    
    % Warm start for next iteration
    X_guess = [X_opt(:,2:end), X_opt(:,end)]; % Shift solution
    U_guess = [U_opt(:,2:end), U_opt(:,end)]; % Shift solution
end

% ----------------------------- PLOTTING --------------------------------
t_sim = 0:dt:T_sim;
t_control_sim = 0:dt:(T_sim-dt);

% States
figure; 
subplot(3,1,1);
plot(t_sim, X_sim(13:15,:));
title('\rho – Relative Position Error'); ylabel('[m]');
legend('\rho_x', '\rho_y', '\rho_z');

subplot(3,1,2);
plot(t_sim, X_sim(16:18,:));
title('\dot{\rho} – Relative Velocity Error'); ylabel('[m/s]');
legend('\dot{\rho}_x', '\dot{\rho}_y', '\dot{\rho}_z');

subplot(3,1,3);
plot(t_sim, X_sim(10:12,:));
title('\omega – Angular Velocity Error'); xlabel('Time [s]'); ylabel('[rad/s]');
legend('\omega_x', '\omega_y', '\omega_z');

% Controls
figure; hold on;
for k = 1:8
    plot(t_control_sim, U_sim(k,:));
end
hold off; grid on;
xlabel('Time [s]'); ylabel('Force [N]');
title('Thrust per Thruster');
legend(arrayfun(@(k) sprintf('f_%d',k), 1:8, 'Uni', 0));