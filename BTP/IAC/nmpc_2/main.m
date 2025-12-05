%% Nonlinear Model Predictive Control for Spacecraft Relative Motion
% Implementation based on "Nonlinear Model Predictive Control of Coupled 
% Rotational-Translational Spacecraft Relative Motion" by Malladi et al.

clear; clc; close all;

%% Spacecraft Parameters
params = struct();
params.mc = 4000;  % Chaser mass [kg]
params.Jc = diag([6.0833, 1.5, 6.0833]) * 1e3;  % Chaser inertia [kg*m^2]
params.Jt = diag([1000, 800, 1200]);  % Target inertia [kg*m^2]
params.mu = 3.986e14;  % Earth gravitational parameter [m^3/s^2]
params.fmax = 20;  % Maximum thrust per thruster [N]

% Feature point locations
params.rpt_t = [1.1404; 3.3462; 5.8907];  % Target feature point [m]
params.rpc_c = [0; 0; -1.75];  % Chaser feature point [m]

% Thruster positions (8 thrusters at corners of spacecraft)
L = 1.5; W = 4; H = 1.5;  % Spacecraft dimensions
params.thruster_pos = [
    -L/2, -W/2, -H/2;
    -L/2, -W/2,  H/2;
    -L/2,  W/2, -H/2;
    -L/2,  W/2,  H/2;
     L/2, -W/2, -H/2;
     L/2, -W/2,  H/2;
     L/2,  W/2, -H/2;
     L/2,  W/2,  H/2
]';

%% Initial Conditions
% Initial relative position and velocity
rc_t_0 = [0.001016052; 116.430912621; 0.020321028];  % [m]
vc_t_0 = [0.020321028; 0; -0.000117733];  % [m/s]

% Initial orientations
Rt_e_0 = eye(3);
Rc_e_0 = diag([-1, -1, 1]);
Rc_t_0 = Rt_e_0' * Rc_e_0;

% Initial angular velocities
omega_t_0 = 0.0046 * [1; 1; 1];  % [rad/s]
omega_c_0 = [0; 0; 0];  % [rad/s]
omega_c_t_0 = Rc_t_0 * omega_c_0 - omega_t_0;

% Initial feature point relative state
rho_c_t_0 = rc_t_0 + Rc_t_0 * params.rpc_c - params.rpt_t;
rho_dot_c_t_0 = [0.0081; 0.0039; 1.3490];  % [m/s]

% Integral state
xi_0 = zeros(3, 1);

% Initial state vector
x0 = [reshape(Rc_t_0, 9, 1); omega_c_t_0; rho_c_t_0; rho_dot_c_t_0; xi_0];

%% MPC Parameters
mpc_params = struct();
mpc_params.N = 8;  % Prediction horizon
mpc_params.T = 480;  % Total prediction time [s]
mpc_params.Ts = mpc_params.T / mpc_params.N;  % Sampling time [s]

% Weight matrices
mpc_params.QR = 70;
mpc_params.Qomega = 5 * eye(3);
mpc_params.Qrho = 0.01 * eye(3);
mpc_params.Qxi = 1e-6 * eye(3);
mpc_params.Wrho_dot = 0.01 * eye(3);
mpc_params.Wtau = 9e-13 * eye(3);
mpc_params.Wf = eye(24);

% Terminal weight
mpc_params.P = eye(6);

%% Constraint Parameters
% Line-of-sight cone constraint (π/4 rad)
los_angle = pi/4;
mpc_params.A_los = [
    1, 0, -1/tan(los_angle);
   -1, 0, -1/tan(los_angle);
    0, 1, -1/tan(los_angle);
    0,-1, -1/tan(los_angle)
];
mpc_params.b_los = zeros(4, 1);

% Thruster gimbal constraints (π/4 rad polyhedral constraint)
gimbal_angle = pi/4;
cos_angle = cos(gimbal_angle);
mpc_params.A_gimbal = [
    1, 0, -1/tan(gimbal_angle);
   -1, 0, -1/tan(gimbal_angle);
    0, 1, -1/tan(gimbal_angle);
    0,-1, -1/tan(gimbal_angle)
];
mpc_params.b_gimbal = zeros(4, 1);

%% Simulation Parameters
sim_time = 800;  % Total simulation time [s]
dt = 1;  % Integration time step [s]
time_vec = 0:dt:sim_time;
num_steps = length(time_vec);

%% Initialize Storage Arrays
state_history = zeros(length(x0), num_steps);
control_history = zeros(24, num_steps);
cost_history = zeros(num_steps, 1);

state_history(:, 1) = x0;

%% Main Simulation Loop
fprintf('Starting spacecraft NMPC simulation...\n');
fprintf('Progress: ');

for k = 1:num_steps-1
    % Progress indicator
    if mod(k, floor(num_steps/20)) == 0
        fprintf('█');
    end
    
    current_state = state_history(:, k);
    current_time = time_vec(k);
    
    % Solve MPC optimization problem
    [u_opt, cost] = solve_mpc(current_state, params, mpc_params);
    
    % Store results
    control_history(:, k) = u_opt;
    cost_history(k) = cost;
    
    % Integrate dynamics
    next_state = integrate_dynamics(current_state, u_opt, dt, params);
    state_history(:, k+1) = next_state;
end

fprintf(' Complete!\n');

%% Post-processing and Visualization
process_results(time_vec, state_history, control_history, cost_history, params);

%% Function Definitions

function [u_opt, cost] = solve_mpc(x_current, params, mpc_params)
    % Solve the MPC optimization problem using fmincon
    
    % Number of control inputs (8 thrusters × 3 components each)
    nu = 24;
    
    % Decision variables: control inputs over prediction horizon
    u_vars = zeros(nu * mpc_params.N, 1);
    
    % Set up optimization options
    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'interior-point', ...
        'MaxIterations', 100, ...
        'ConstraintTolerance', 1e-6, ...
        'OptimalityTolerance', 1e-6);
    
    % Bounds on control inputs
    lb = repmat(-params.fmax * ones(nu, 1), mpc_params.N, 1);
    ub = repmat(params.fmax * ones(nu, 1), mpc_params.N, 1);
    
    % Solve optimization
    try
        [u_opt_vec, cost] = fmincon(@(u) mpc_cost_function(u, x_current, params, mpc_params), ...
            u_vars, [], [], [], [], lb, ub, ...
            @(u) mpc_constraints(u, x_current, params, mpc_params), options);
        
        % Extract first control input
        u_opt = u_opt_vec(1:nu);
    catch
        warning('MPC optimization failed, using zero control');
        u_opt = zeros(nu, 1);
        cost = inf;
    end
end

function cost = mpc_cost_function(u_vec, x0, params, mpc_params)
    % Compute the MPC cost function
    
    nu = 24;
    N = mpc_params.N;
    Ts = mpc_params.Ts;
    
    cost = 0;
    x = x0;
    
    % Stage costs
    for i = 1:N
        u = u_vec((i-1)*nu + 1 : i*nu);
        
        % Extract state components
        Rc_t = reshape(x(1:9), 3, 3);
        omega_c_t = x(10:12);
        rho_c_t = x(13:15);
        rho_dot_c_t = x(16:18);
        xi = x(19:21);
        
        % Compute torques from thrust forces
        tau_c = compute_torques(u, params);
        
        % Stage cost
        cost_R = mpc_params.QR * trace(eye(3) - Rc_t);
        cost_omega = omega_c_t' * mpc_params.Qomega * omega_c_t;
        cost_rho = rho_c_t' * mpc_params.Qrho * rho_c_t;
        cost_rho_dot = rho_dot_c_t' * mpc_params.Wrho_dot * rho_dot_c_t;
        cost_xi = xi' * mpc_params.Qxi * xi;
        cost_tau = tau_c' * mpc_params.Wtau * tau_c;
        cost_f = u' * mpc_params.Wf * u;
        
        stage_cost = cost_R + cost_omega + cost_rho + cost_rho_dot + cost_xi + cost_tau + cost_f;
        cost = cost + stage_cost;
        
        % Propagate dynamics
        if i < N
            x = integrate_dynamics(x, u, Ts, params);
        end
    end
    
    % Terminal cost
    rho_terminal = x(13:15);
    rho_dot_terminal = x(16:18);
    terminal_state = [rho_terminal; rho_dot_terminal];
    cost = cost + terminal_state' * mpc_params.P * terminal_state;
end

function [c, ceq] = mpc_constraints(u_vec, x0, params, mpc_params)
    % MPC constraints
    
    nu = 24;
    N = mpc_params.N;
    Ts = mpc_params.Ts;
    
    c = [];
    ceq = [];
    x = x0;
    
    for i = 1:N
        u = u_vec((i-1)*nu + 1 : i*nu);
        
        % Thruster constraints
        for j = 1:8
            f_thruster = u((j-1)*3 + 1 : j*3);
            
            % Magnitude constraint
            c = [c; norm(f_thruster) - params.fmax];
            
            % Gimbal constraints (simplified as magnitude bounds)
            % In practice, this would be the polyhedral constraint
        end
        
        % Propagate state
        x = integrate_dynamics(x, u, Ts, params);
        
        % Line-of-sight constraint
        rho_c_t = x(13:15);
        if norm(rho_c_t) > 1e-6  % Avoid division by zero
            los_violations = mpc_params.A_los * rho_c_t - mpc_params.b_los;
            c = [c; los_violations];
        end
    end
end

function x_next = integrate_dynamics(x, u, dt, params)
    % Integrate spacecraft dynamics using Runge-Kutta 4th order
    
    k1 = spacecraft_dynamics(x, u, params);
    k2 = spacecraft_dynamics(x + dt/2 * k1, u, params);
    k3 = spacecraft_dynamics(x + dt/2 * k2, u, params);
    k4 = spacecraft_dynamics(x + dt * k3, u, params);
    
    x_next = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
    
    % Ensure rotation matrix remains orthogonal
    Rc_t = reshape(x_next(1:9), 3, 3);
    [U, ~, V] = svd(Rc_t);
    Rc_t = U * V';
    x_next(1:9) = reshape(Rc_t, 9, 1);
end

function dx = spacecraft_dynamics(x, u, params)
    % Spacecraft relative motion dynamics
    
    % Extract state components
    Rc_t = reshape(x(1:9), 3, 3);
    omega_c_t = x(10:12);
    rho_c_t = x(13:15);
    rho_dot_c_t = x(16:18);
    xi = x(19:21);
    
    % Compute forces and torques
    f_total = sum(reshape(u, 3, 8), 2);
    tau_c = compute_torques(u, params);
    
    % Attitude kinematics
    dRc_t = Rc_t * skew_symmetric(omega_c_t);
    
    % Attitude dynamics (simplified - assumes zero target angular acceleration)
    Jc_inv = inv(params.Jc);
    domega_c_t = Jc_inv * (-cross(params.Jc * omega_c_t, omega_c_t) + tau_c);
    
    % Translational dynamics (simplified orbital mechanics)
    % This is a simplified version - full implementation would include
    % all coupling terms from equation (9) in the paper
    drho_c_t = rho_dot_c_t;
    
    % Simplified relative acceleration
    a_control = f_total / params.mc;
    drho_dot_c_t = Rc_t * a_control + cross(omega_c_t, cross(omega_c_t, Rc_t * params.rpc_c));
    
    % Integral action
    dxi = rho_c_t;
    
    % Combine derivatives
    dx = [reshape(dRc_t, 9, 1); domega_c_t; drho_c_t; drho_dot_c_t; dxi];
end

function tau = compute_torques(u_vec, params)
    % Compute torques from thruster forces
    
    tau = zeros(3, 1);
    for i = 1:8
        f_thruster = u_vec((i-1)*3 + 1 : i*3);
        r_thruster = params.thruster_pos(:, i);
        tau = tau + cross(r_thruster, f_thruster);
    end
end

function S = skew_symmetric(v)
    % Create skew-symmetric matrix from vector
    S = [0, -v(3), v(2);
         v(3), 0, -v(1);
         -v(2), v(1), 0];
end

function process_results(time_vec, state_history, control_history, cost_history, params)
    % Process and visualize simulation results
    
    fprintf('Processing results and generating plots...\n');
    
    % Extract state components
    num_steps = length(time_vec);
    attitude_error = zeros(num_steps, 1);
    omega_error = zeros(3, num_steps);
    position_error = zeros(3, num_steps);
    velocity_error = zeros(3, num_steps);
    
    for k = 1:num_steps
        Rc_t = reshape(state_history(1:9, k), 3, 3);
        attitude_error(k) = acos((trace(Rc_t) - 1) / 2);
        omega_error(:, k) = state_history(10:12, k);
        position_error(:, k) = state_history(13:15, k);
        velocity_error(:, k) = state_history(16:18, k);
    end
    
    % Create figure with subplots
    figure('Position', [100, 100, 1200, 800]);
    
    % Attitude error
    subplot(2, 3, 1);
    plot(time_vec, rad2deg(attitude_error), 'b-', 'LineWidth', 2);
    grid on;
    xlabel('Time [s]');
    ylabel('Attitude Error [deg]');
    title('Relative Attitude Error');
    
    % Angular velocity error
    subplot(2, 3, 2);
    plot(time_vec, omega_error(1, :), 'r-', 'LineWidth', 2); hold on;
    plot(time_vec, omega_error(2, :), 'g-', 'LineWidth', 2);
    plot(time_vec, omega_error(3, :), 'b-', 'LineWidth', 2);
    grid on;
    xlabel('Time [s]');
    ylabel('Angular Velocity Error [rad/s]');
    title('Relative Angular Velocity Error');
    legend('\omega_x', '\omega_y', '\omega_z', 'Location', 'best');
    
    % Position error
    subplot(2, 3, 3);
    plot(time_vec, position_error(1, :), 'r-', 'LineWidth', 2); hold on;
    plot(time_vec, position_error(2, :), 'g-', 'LineWidth', 2);
    plot(time_vec, position_error(3, :), 'b-', 'LineWidth', 2);
    grid on;
    xlabel('Time [s]');
    ylabel('Position Error [m]');
    title('Relative Position Error');
    legend('X', 'Y', 'Z', 'Location', 'best');
    
    % Velocity error
    subplot(2, 3, 4);
    plot(time_vec, velocity_error(1, :), 'r-', 'LineWidth', 2); hold on;
    plot(time_vec, velocity_error(2, :), 'g-', 'LineWidth', 2);
    plot(time_vec, velocity_error(3, :), 'b-', 'LineWidth', 2);
    grid on;
    xlabel('Time [s]');
    ylabel('Velocity Error [m/s]');
    title('Relative Velocity Error');
    legend('V_x', 'V_y', 'V_z', 'Location', 'best');
    
    % Control effort
    subplot(2, 3, 5);
    control_magnitude = sqrt(sum(reshape(control_history.^2, 3, 8, []), 1));
    control_magnitude = squeeze(control_magnitude);
    plot(time_vec(1:end-1), sum(control_magnitude, 1), 'k-', 'LineWidth', 2);
    grid on;
    xlabel('Time [s]');
    ylabel('Total Thrust Magnitude [N]');
    title('Control Effort');
    
    % Cost evolution
    subplot(2, 3, 6);
    plot(time_vec(1:end-1), cost_history(1:end-1), 'm-', 'LineWidth', 2);
    grid on;
    xlabel('Time [s]');
    ylabel('Cost Function Value');
    title('MPC Cost Evolution');
    
    % 3D trajectory plot
    figure('Position', [200, 200, 800, 600]);
    plot3(position_error(1, :), position_error(2, :), position_error(3, :), 'b-', 'LineWidth', 2);
    hold on;
    plot3(position_error(1, 1), position_error(2, 1), position_error(3, 1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot3(position_error(1, end), position_error(2, end), position_error(3, end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    grid on;
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    title('3D Relative Trajectory');
    legend('Trajectory', 'Start', 'End', 'Location', 'best');
    axis equal;
    
    % Display final errors
    fprintf('\n=== Final Performance Metrics ===\n');
    fprintf('Final attitude error: %.4f deg\n', rad2deg(attitude_error(end)));
    fprintf('Final position error magnitude: %.4f m\n', norm(position_error(:, end)));
    fprintf('Final velocity error magnitude: %.4f m/s\n', norm(velocity_error(:, end)));
    fprintf('Average control effort: %.2f N\n', mean(sum(control_magnitude, 1)));
    fprintf('Simulation completed successfully!\n');
end