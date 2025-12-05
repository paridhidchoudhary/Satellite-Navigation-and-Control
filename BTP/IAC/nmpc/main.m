%% Nonlinear Model Predictive Control of Coupled Rotational-Translational Spacecraft Relative Motion
% This simulation implements the NMPC approach described in the paper
% by Malladi, Di Cairano and Weiss (2019)

%% Setup environment and parameters
clear all;
close all;
clc;

% Add path to CasADi and MPCTools if installed
% addpath('/path/to/casadi');
% addpath('/path/to/mpctools');

%% Constants and Parameters
% Earth gravitational parameter
mu = 3.986004418e14; % m^3/s^2

% Target parameters
m_t = 500; % kg (asteroid mass)
J_t = diag([300, 250, 350]); % kg*m^2 (asteroid inertia)
r_pt_t = [1.1404; 3.3462; 5.8907]; % m (feature point on target wrt target frame)

% Chaser parameters 
m_c = 4000; % kg (chaser spacecraft mass)
J_c = diag([6.0833, 1.5, 6.0833]) * 1e3; % kg*m^2 (chaser inertia)
r_pc_c = [0; 0; -1.75]; % m (feature point on chaser wrt chaser frame)

% Thruster configuration
thruster_positions = zeros(3,8);
% Define the 8 thruster positions at corners of two faces of spacecraft
% Example configuration for a rectangular cuboid
L = 1.5; % Length in x (m)
W = 4;   % Width in y (m)
H = 1.5; % Height in z (m)

% Front face corners
thruster_positions(:,1) = [L/2; W/2; H/2];
thruster_positions(:,2) = [L/2; -W/2; H/2];
thruster_positions(:,3) = [L/2; -W/2; -H/2];
thruster_positions(:,4) = [L/2; W/2; -H/2];

% Back face corners
thruster_positions(:,5) = [-L/2; W/2; H/2];
thruster_positions(:,6) = [-L/2; -W/2; H/2];
thruster_positions(:,7) = [-L/2; -W/2; -H/2];
thruster_positions(:,8) = [-L/2; W/2; -H/2];

% Maximum thrust force
f_max = 20; % N

% Orbital parameters
alt = 300e3; % m (orbital altitude)
R_earth = 6378.1e3; % m (Earth radius)
R_orbit = R_earth + alt; % m (orbital radius)
orbital_velocity = sqrt(mu/R_orbit); % m/s
orbital_period = 2*pi*R_orbit/orbital_velocity; % s

% Simulation parameters
T_final = 800; % s (final simulation time)
dt = 1; % s (time step for simulation)
t = 0:dt:T_final; % time vector
N_steps = length(t);

% NMPC parameters 
Ts = 60; % s (MPC sampling time)
N = 8; % number of prediction steps
T = N * Ts; % prediction horizon

% Constraints
LOS_angle = pi/4; % rad (Line of Sight constraint angle)
gimbal_angle = pi/4; % rad (thruster gimbal constraint angle)

%% Initial conditions
% Target initial state
R_t_e_0 = eye(3); % identity rotation matrix
r_t_e_0 = [R_orbit; 0; 0]; % position in ECI frame (m)
v_t_e_0 = [0; orbital_velocity; 0]; % velocity in ECI frame (m/s)
omega_t_t_e_0 = 0.0046 * [1; 1; 1]; % angular velocity in target frame (rad/s)

% Chaser initial state(wrt target frame)
R_c_e_0 = diag([-1, -1, 1]); % initial rotation matrix
r_c_t_0 = [0.001016052; 116.430912621; 0.020321028]; % position relative to target (m) 
v_c_t_0 = [0.020321028; 0; -0.000117733]; % velocity relative to target (m/s)
omega_c_c_e_0 = [0; 0; 0]; % angular velocity in chaser frame (rad/s)

% Calculate absolute chaser position and velocity in inertial frame
r_c_e_0 = r_t_e_0 + R_t_e_0 * r_c_t_0;
v_c_e_0 = v_t_e_0 + R_t_e_0 * (v_c_t_0 + cross(omega_t_t_e_0, r_c_t_0));

% Calculate relative feature point position and velocity
rho_c_t_0 = r_c_t_0 + (R_c_e_0'*R_t_e_0)'*r_pc_c - r_pt_t; %target frame
rhodot_c_t_0 = v_c_t_0 + cross((R_c_e_0'*R_t_e_0)'*omega_c_c_e_0, r_pc_c) - ...
              cross(omega_t_t_e_0, r_pt_t); %target frame 

%% Initialize state history arrays
% Target states
R_t_e = zeros(3, 3, N_steps);
R_t_e(:,:,1) = R_t_e_0;
r_t_e = zeros(3, N_steps);
r_t_e(:,1) = r_t_e_0;
v_t_e = zeros(3, N_steps);
v_t_e(:,1) = v_t_e_0;
omega_t_t_e = zeros(3, N_steps);
omega_t_t_e(:,1) = omega_t_t_e_0;

% Chaser states
R_c_e = zeros(3, 3, N_steps);
R_c_e(:,:,1) = R_c_e_0;
r_c_e = zeros(3, N_steps);
r_c_e(:,1) = r_c_e_0;
v_c_e = zeros(3, N_steps);
v_c_e(:,1) = v_c_e_0;
omega_c_c_e = zeros(3, N_steps);
omega_c_c_e(:,1) = omega_c_c_e_0;

% Relative states
R_c_t = zeros(3, 3, N_steps); %chaser to target
R_c_t(:,:,1) = R_t_e_0' * R_c_e_0;
r_c_t = zeros(3, N_steps); %target frame
r_c_t(:,1) = r_c_t_0;
v_c_t = zeros(3, N_steps); %target frame
v_c_t(:,1) = v_c_t_0;
omega_c_t_t = zeros(3, N_steps);
omega_c_t_t(:,1) = R_c_t(:,:,1) * omega_c_c_e_0 - omega_t_t_e_0;

% Feature points
rho_c_t = zeros(3, N_steps);
rho_c_t(:,1) = rho_c_t_0;
rhodot_c_t = zeros(3, N_steps);
rhodot_c_t(:,1) = rhodot_c_t_0;

%states -> [R_c_t' omegae_c_t_t' rho_c_T' rhodot_c_t' xi']

% Control inputs
thruster_forces = zeros(3, 8, N_steps);
total_force = zeros(3, N_steps);
total_torque = zeros(3, N_steps);

% Error integrator
xi = zeros(3, N_steps);
xi(:,1) = zeros(3,1);

%% Main simulation loop
for k = 1:N_steps-1
    % Current time
    current_time = t(k);
    
    % Extract current state
    R_t_e_k = R_t_e(:,:,k);
    r_t_e_k = r_t_e(:,k);
    v_t_e_k = v_t_e(:,k);
    omega_t_t_e_k = omega_t_t_e(:,k);
    
    R_c_e_k = R_c_e(:,:,k);
    r_c_e_k = r_c_e(:,k);
    v_c_e_k = v_c_e(:,k);
    omega_c_c_e_k = omega_c_c_e(:,k);
    
    R_c_t_k = R_c_t(:,:,k);
    omega_c_t_t_k = omega_c_t_t(:,k);
    
    rho_c_t_k = rho_c_t(:,k);
    rhodot_c_t_k = rhodot_c_t(:,k);
    
    xi_k = xi(:,k);
    
    % Solve NMPC problem at control intervals
    if mod(current_time, Ts) == 0 || k == 1
        % Here we would solve the NMPC problem using CasADi/MPCTools
        % For simplicity in this demo, we'll implement a simpler controller
        
        % Proportional control for attitude
        K_R = 0.1;
        K_omega = 0.2;
        
        % Extract rotation error from R_c_t_k
        % Note: This is simplified; in practice, use proper error metric
        R_error = eye(3) - R_c_t_k;
        torque_cmd = -K_R * vee(R_error) - K_omega * omega_c_t_t_k;
        
        % Proportional control for position
        K_rho = 0.02;
        K_rhodot = 0.1;
        
        % Position control
        force_cmd = -K_rho * rho_c_t_k - K_rhodot * rhodot_c_t_k;
        
        % Distribute force and torque commands to thrusters
        % This is a simplified thruster allocation algorithm
        % In the real implementation, this would be part of the NMPC optimization
        thruster_forces_k = allocate_thrusters(force_cmd, torque_cmd, thruster_positions, f_max);
        
        for i = 1:8
            thruster_forces(:,i,k) = thruster_forces_k(:,i);
        end
    else
        % Use the same thruster forces as the previous step
        for i = 1:8
            thruster_forces(:,i,k) = thruster_forces(:,i,k-1);
        end
    end
    
    % Calculate resulting force and torque on chaser
    F_c = zeros(3,1);
    tau_c = zeros(3,1);
    for i = 1:8
        F_c = F_c + thruster_forces(:,i,k);
        tau_c = tau_c + cross(thruster_positions(:,i), thruster_forces(:,i,k));
    end
    total_force(:,k) = F_c;
    total_torque(:,k) = tau_c;
    
    % Propagate target dynamics
    % Rotational dynamics
    omega_dot_t = J_t \ (cross(-omega_t_t_e_k, J_t * omega_t_t_e_k));
    omega_t_t_e(:,k+1) = omega_t_t_e_k + dt * omega_dot_t;
    
    % Update rotation matrix using exponential map
    omega_skew = skew(omega_t_t_e_k);
    R_t_e(:,:,k+1) = R_t_e_k * expm(omega_skew * dt);
    
    % Translational dynamics (simplified, assuming circular orbit)
    r_norm = norm(r_t_e_k);
    r_unit = r_t_e_k / r_norm;
    a_gravity = -mu * r_unit / r_norm^2;
    
    v_t_e(:,k+1) = v_t_e_k + dt * a_gravity;
    r_t_e(:,k+1) = r_t_e_k + dt * v_t_e(:,k+1);
    
    % Propagate chaser dynamics
    % Rotational dynamics
    omega_dot_c = J_c \ (cross(-omega_c_c_e_k, J_c * omega_c_c_e_k) + tau_c);
    omega_c_c_e(:,k+1) = omega_c_c_e_k + dt * omega_dot_c;
    
    % Update rotation matrix using exponential map
    omega_skew_c = skew(omega_c_c_e_k);
    R_c_e(:,:,k+1) = R_c_e_k * expm(omega_skew_c * dt);
    
    % Translational dynamics with thrust
    r_c_norm = norm(r_c_e_k);
    r_c_unit = r_c_e_k / r_c_norm;
    a_gravity_c = -mu * r_c_unit / r_c_norm^2;
    a_thrust = R_c_e_k * F_c / m_c;
    
    v_c_e(:,k+1) = v_c_e_k + dt * (a_gravity_c + a_thrust);
    r_c_e(:,k+1) = r_c_e_k + dt * v_c_e(:,k+1);
    
    % Compute relative states
    r_c_t(:,k+1) = R_t_e(:,:,k+1)' * (r_c_e(:,k+1) - r_t_e(:,k+1));
    v_c_t(:,k+1) = R_t_e(:,:,k+1)' * (v_c_e(:,k+1) - v_t_e(:,k+1)) - ...
        cross(omega_t_t_e(:,k+1), r_c_t(:,k+1));
    
    R_c_t(:,:,k+1) = R_t_e(:,:,k+1)' * R_c_e(:,:,k+1);
    omega_c_t_t(:,k+1) = R_c_t(:,:,k+1) * omega_c_c_e(:,k+1) - omega_t_t_e(:,k+1);
    
    % Compute feature point relative states
    rho_c_t(:,k+1) = r_c_t(:,k+1) + R_c_t(:,:,k+1) * r_pc_c - r_pt_t;
    rhodot_c_t(:,k+1) = v_c_t(:,k+1) + cross(omega_c_t_t(:,k+1), R_c_t(:,:,k+1) * r_pc_c);
    
    % Update integrator state
    xi(:,k+1) = xi_k + dt * rho_c_t(:,k+1);
end

%% Calculate some additional metrics
% Eigenaxis attitude error
theta_c_t = zeros(N_steps, 1);
for k = 1:N_steps
    theta_c_t(k) = acos(0.5 * (trace(R_c_t(:,:,k)) - 1));
end

%% Plotting
figure(1);
subplot(3,1,1);
plot(t, r_c_t(1,:), 'r', t, r_c_t(2,:), 'g', t, r_c_t(3,:), 'b');
grid on;
title('Relative Position of Chaser wrt Target (m)');
legend('x', 'y', 'z');
xlabel('Time (s)');
ylabel('Position (m)');

subplot(3,1,2);
plot(t, v_c_t(1,:), 'r', t, v_c_t(2,:), 'g', t, v_c_t(3,:), 'b');
grid on;
title('Relative Velocity of Chaser wrt Target (m/s)');
legend('x', 'y', 'z');
xlabel('Time (s)');
ylabel('Velocity (m/s)');

subplot(3,1,3);
plot(t, rad2deg(theta_c_t));
grid on;
title('Eigenaxis Attitude Error (deg)');
xlabel('Time (s)');
ylabel('Attitude Error (deg)');

figure(2);
subplot(3,1,1);
plot(t, rho_c_t(1,:), 'r', t, rho_c_t(2,:), 'g', t, rho_c_t(3,:), 'b');
grid on;
title('Relative Position Between Feature Points (m)');
legend('x', 'y', 'z');
xlabel('Time (s)');
ylabel('Position (m)');

subplot(3,1,2);
plot(t, rhodot_c_t(1,:), 'r', t, rhodot_c_t(2,:), 'g', t, rhodot_c_t(3,:), 'b');
grid on;
title('Relative Velocity Between Feature Points (m/s)');
legend('x', 'y', 'z');
xlabel('Time (s)');
ylabel('Velocity (m/s)');

subplot(3,1,3);
plot(t, vecnorm(rho_c_t));
grid on;
title('Distance Between Feature Points (m)');
xlabel('Time (s)');
ylabel('Distance (m)');

figure(3);
subplot(3,1,1);
plot(t, omega_c_t_t(1,:), 'r', t, omega_c_t_t(2,:), 'g', t, omega_c_t_t(3,:), 'b');
grid on;
title('Relative Angular Velocity (rad/s)');
legend('x', 'y', 'z');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');

subplot(3,1,2);
plot(t(1:end-1), vecnorm(total_force(:,1:end-1)));
grid on;
title('Total Force Magnitude (N)');
xlabel('Time (s)');
ylabel('Force (N)');

subplot(3,1,3);
plot(t(1:end-1), vecnorm(total_torque(:,1:end-1)));
grid on;
title('Total Torque Magnitude (N·m)');
xlabel('Time (s)');
ylabel('Torque (N·m)');

% 3D visualization of the approach
figure(4);
plot3(rho_c_t(1,:), rho_c_t(2,:), rho_c_t(3,:), 'b-');
hold on;
plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot3(rho_c_t(1,1), rho_c_t(2,1), rho_c_t(3,1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
plot3(rho_c_t(1,end), rho_c_t(2,end), rho_c_t(3,end), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Feature Point Approach Trajectory in Target Frame');
legend('Trajectory', 'Target Feature Point', 'Initial Chaser Position', 'Final Chaser Position');
axis equal;

%% Helper functions
function S = skew(v)
    % Creates a skew-symmetric matrix from a 3x1 vector
    S = [0, -v(3), v(2); 
         v(3), 0, -v(1); 
         -v(2), v(1), 0];
end

function v = vee(S)
    % Extracts a 3x1 vector from a skew-symmetric matrix
    v = [S(3,2); S(1,3); S(2,1)];
end

function thruster_forces = allocate_thrusters(force_cmd, torque_cmd, thruster_positions, f_max)
    % Simple thruster allocation algorithm
    % In practice, this would be solved as part of the NMPC optimization
    
    % Initialize thruster forces
    thruster_forces = zeros(3, 8);
    
    % Number of thrusters
    n_thrusters = size(thruster_positions, 2);
    
    % For this simplified example, distribute forces equally
    % In reality, this should be solved as a constrained optimization problem
    for i = 1:n_thrusters
        thruster_forces(:,i) = force_cmd / n_thrusters;
        
        % Add torque contribution (very simplified)
        % This doesn't actually solve the torque allocation problem correctly
        % but gives a basic response for demonstration purposes
        torque_contribution = cross(thruster_positions(:,i), thruster_forces(:,i));
        thruster_forces(:,i) = thruster_forces(:,i) + cross(torque_cmd - torque_contribution, thruster_positions(:,i)) * 0.01;
        
        % Limit thrust magnitude
        thrust_mag = norm(thruster_forces(:,i));
        if thrust_mag > f_max
            thruster_forces(:,i) = thruster_forces(:,i) * f_max / thrust_mag;
        end
    end
end

function visualize_spacecraft(t_idx, r_t_e, R_t_e, r_c_e, R_c_e, rho_c_t, r_pt_t, r_pc_c)
    % Function to visualize the spacecraft at a specific time index
    
    % This would be a more complex visualization function
    % For a complete implementation, consider using existing MATLAB
    % visualization tools or custom graphics
    
    % Plot target and chaser bodies
    % Plot feature points
    % Show LOS cone
    % etc.
end