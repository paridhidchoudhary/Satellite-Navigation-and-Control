% =======================================================================
%         NMPC Simulation Loop for Coupled Spacecraft Control
%         (rev-B – scalar 8-jet model, finer grid, thrust bounds)
% =======================================================================

import casadi.*

% ----------------------------- HORIZON ---------------------------------
T  = 960;          % [s] total horizon (unchanged)
N  = 320;          % •finer grid: dt = 3 s 
dt = T / N;        % 6 s control period

% -------------------------- DYNAMICS & COST ----------------------------
fdyn  = f_dyn();   % (new scalar-jet dynamics)
fcost = f_cost();  % (new cost with thrust penalty)

% ---------------------------- DIMENSIONS -------------------------------
nx = 21;           % state length (unchanged)
nu = 8;            % •8 scalar thrusters (was 24)•

% ------------------------ DECISION VARIABLES ---------------------------
X = SX.sym('X', nx, N+1);   % states
U = SX.sym('U', nu, N);     % controls

% ------------------------- OBJECTIVE & CONSTRAINTS ---------------------
J   = 0;
g   = {};   lbg = [];   ubg = [];

Wr  = 1e-2 * eye(nu);       % Δu regularisation (resized to 8×8)

for k = 1:N
    % stage cost
    J = J + fcost(X(:,k), U(:,k));

    % Δu penalty
    if k > 1
        du = U(:,k) - U(:,k-1);
        J  = J + du' * Wr * du;
    end

    % RK-4 integration
    xk = X(:,k);   uk = U(:,k);
    k1 = fdyn(xk,                 uk);
    k2 = fdyn(xk + 0.5*dt*k1,     uk);
    k3 = fdyn(xk + 0.5*dt*k2,     uk);
    k4 = fdyn(xk +     dt*k3,     uk);
    x_next = xk + dt/6*(k1 + 2*k2 + 2*k3 + k4);

    % dynamics equality
    g{end+1} = X(:,k+1) - x_next;
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];
end

% terminal cost on (ρ, ρ̇)
P = 1e4 * eye(6);
rhoT    = X(13:15,end);
rhodT   = X(16:18,end);
J = J + [rhoT; rhodT]' * P * [rhoT; rhodT];

% --------------------------- NLP STRUCTURE -----------------------------
w = [X(:); U(:)];
g = vertcat(g{:});
nlp = struct('x', w, 'f', J, 'g', g);
opts = struct('ipopt', struct('print_level', 0, 'max_iter', 500));
solver = nlpsol('solver', 'ipopt', nlp, opts);

% ----------------------- INITIAL GUESS ---------------------------------
x0 = zeros(nx,1);
x0([1 5 9]) = 1;                 % identity rotation
x0(13:15)   = [-1.1394; 113.0847; -7.6204];
x0(16:18)   = [ 0.0081;   0.0039;  1.3490];

X_init = repmat(x0, 1, N+1);
U_init = 0.5 + 0.5*rand(nu, N);        % small bias helps IPOPT

% ------------------------ VARIABLE BOUNDS ------------------------------
u_max = 15;                       % [N] per-jet limit
lbx   = -inf(size(w));
ubx   =  inf(size(w));

u_start = nx*(N+1) + 1;           % first control index in w
lbx(u_start:end) = 0;             % 0 ≤ u ≤ u_max
ubx(u_start:end) = u_max;

% ------------------------------- SOLVE ---------------------------------
sol   = solver('x0', [X_init(:); U_init(:)], ...
               'lbx', lbx, 'ubx', ubx, ...
               'lbg', lbg, 'ubg', ubg);

% ---------------------------- UNPACK -----------------------------------
w_opt = full(sol.x);
X_opt = reshape(w_opt(1:nx*(N+1)), nx, N+1);
U_opt = reshape(w_opt(nx*(N+1)+1:end), nu, N);

% -------------------------- TIME VECTORS -------------------------------
t_full    = 0 : dt : T;
t_control = 0 : dt : T-dt;

% ----------------------------- PLOTS -----------------------------------
% ρ (position error)
figure; subplot(3,1,1);
plot(t_full, X_opt(13:15,:));
title('\rho – Relative Position Error'); ylabel('[m]');
legend('\rho_x', '\rho_y', '\rho_z');

% ρ̇ (velocity error)
subplot(3,1,2);
plot(t_full, X_opt(16:18,:));
title('\dot{\rho} – Relative Velocity Error'); ylabel('[m/s]');
legend('\dot{\rho}_x', '\dot{\rho}_y', '\dot{\rho}_z');

% ω (angular velocity error)
subplot(3,1,3);
plot(t_full, X_opt(10:12,:));
title('\omega – Angular Velocity Error'); xlabel('Time [s]'); ylabel('[rad/s]');
legend('\omega_x', '\omega_y', '\omega_z');

% 3-D trajectory
figure;
plot3(X_opt(13,:), X_opt(14,:), X_opt(15,:));
hold on; plot3(0,0,0,'r*','MarkerSize',10); hold off;
grid on; xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Chaser Trajectory in Target Frame'); legend('path','target');

% Thrust magnitudes per jet
figure; hold on;
for k = 1:8
    plot(t_control, U_opt(k,:));
end
hold off; grid on;
xlabel('Time [s]'); ylabel('Force [N]');
title('Thrust per Thruster'); 
legend(arrayfun(@(k) sprintf('f_%d',k), 1:8, 'Uni', 0));
