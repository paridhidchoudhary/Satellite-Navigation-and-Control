% NMPC Simulation Loop for Coupled Spacecraft Control
% Assumes f_dyn() and f_cost() are already defined

import casadi.*

% Time horizon and steps
T = 960; N = 16; dt = T/N;  % Increased simulation time

% Instantiate dynamics and cost functions
fdyn = f_dyn();
fcost = f_cost();

% State and control dimensions
nx = 21; % 9 (Rct) + 3 (omega) + 3 (rho) + 3 (rho_dot) + 3 (xi)
nu = 24; % 8 thrusters * 3D force vectors

% Symbolic decision variables
X = SX.sym('X', nx, N+1);
U = SX.sym('U', nu, N);

% Initialize objective and constraints
J = 0; g = {}; lbg = []; ubg = [];

% Cost penalty for control rate change
Wr = 1e-4 * eye(nu);

% Thruster magnitude limit (N)
fmax = 10;  % tighten to force the optimiser to share load

% Loop over horizon
for k = 1:N
    % Stage cost
    J = J + fcost(X(:,k), U(:,k));

    % Control rate penalty (for k > 1)
    if k > 1
        delta_u = U(:,k) - U(:,k-1);
        J = J + delta_u' * Wr * delta_u;
    end

    % -----------------------------
    % Add thruster‑magnitude constraints
    % -----------------------------
    for th = 1:8
        fk = U(3*(th-1)+1:3*th, k);
        g{end+1}  = norm_2(fk) - fmax;  % must be <= 0
        lbg       = [lbg; -inf];
        ubg       = [ubg;  0];
    end

    % RK4 integration instead of Euler
    xk = X(:,k);
    uk = U(:,k);
    k1 = fdyn(xk, uk);
    k2 = fdyn(xk + 0.5*dt*k1, uk);
    k3 = fdyn(xk + 0.5*dt*k2, uk);
    k4 = fdyn(xk + dt*k3, uk);
    x_next = xk + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    % Continuity constraint
    g{end+1} = X(:,k+1) - x_next;
    lbg = [lbg; zeros(nx,1)];
    ubg = [ubg; zeros(nx,1)];
end
    

% Terminal cost (only on rho and rho_dot)
P = 1000 * eye(6);
rho_terminal = X(13:15, end);
rhodot_terminal = X(16:18, end);
J = J + [rho_terminal; rhodot_terminal]' * P * [rho_terminal; rhodot_terminal];

% Flatten variables for solver
w = [reshape(X, nx*(N+1), 1); reshape(U, nu*N, 1)];
g = vertcat(g{:});

% Define NLP
nlp = struct('x', w, 'f', J, 'g', g);
opts = struct('ipopt', struct('print_level', 0, 'max_iter', 500));
solver = nlpsol('solver', 'ipopt', nlp, opts);

% Initial conditions (from paper values)
x0 = zeros(nx,1);
x0(1) = 1; x0(5) = 1; x0(9) = 1;  % Rc/t = identity rotation
x0(10:12) = [0; 0; 0];            % Angular velocity error
x0(13:15) = [-1.1394; 113.0847; -7.6204];  % Relative position (pc/pt)
x0(16:18) = [0.0081; 0.0039; 1.3490];       % Relative velocity
x0(19:21) = [0; 0; 0];           % Integral of position error

% Define target position (final goal, origin in target frame)
target_position = [0; 0; 0];  % End effector should reach surface rock

X_init = repmat(x0, 1, N+1);
U_init = zeros(nu, N);

% Solve NLP
sol = solver('x0', [X_init(:); U_init(:)], 'lbg', lbg, 'ubg', ubg);

% Extract optimal trajectories
w_opt = full(sol.x);
X_opt = reshape(w_opt(1:nx*(N+1)), nx, N+1);
U_opt = reshape(w_opt(nx*(N+1)+1:end), nu, N);

% Time vectors
t_full = 0:dt:T;
t_control = 0:dt:T-dt;

% Extract variables
Rct_traj = X_opt(1:9, :);
omega_traj = X_opt(10:12, :);
rho_traj = X_opt(13:15, :);
rho_dot_traj = X_opt(16:18, :);
xi_traj = X_opt(19:21, :);

% Plot relative position error (rho)
figure;
subplot(3,1,1);
plot(t_full, rho_traj);
title('Relative Position Error (\rho)');
xlabel('Time [s]'); ylabel('\rho [m]');
legend('\rho_x', '\rho_y', '\rho_z');

% Plot velocity error (rho_dot)
subplot(3,1,2);
plot(t_full, rho_dot_traj);
title('Relative Velocity Error (\rho dot)');
xlabel('Time [s]'); ylabel('\dot{\rho} [m/s]');
legend('\dot{\rho}_x', '\dot{\rho}_y', '\dot{\rho}_z');

% Plot angular velocity error (omega)
subplot(3,1,3);
plot(t_full, omega_traj);
title('Angular Velocity Error (\omega)');
xlabel('Time [s]'); ylabel('\omega [rad/s]');
legend('\omega_x', '\omega_y', '\omega_z');

% Plot trajectory in 3D (separate figure)
figure;
plot3(rho_traj(1,:), rho_traj(2,:), rho_traj(3,:), '-o');
hold on;
plot3(rho_traj(1,1), rho_traj(2,1), rho_traj(3,1), 'bo', 'MarkerSize', 10, 'DisplayName', 'Start');
plot3(target_position(1), target_position(2), target_position(3), 'r*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Target');
hold off;
grid on;
title('Chaser Trajectory in Target Frame');
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
legend;

% Plot thrust magnitudes of all thrusters
figure;
hold on;
for k = 1:8
    fk = U_opt(3*(k-1)+1:3*k, :);
    plot(t_control, vecnorm(fk, 2, 1));
end
hold off;
title('Thrust Magnitudes per Thruster');
xlabel('Time [s]'); ylabel('Force [N]');
legend('f_1','f_2','f_3','f_4','f_5','f_6','f_7','f_8');
