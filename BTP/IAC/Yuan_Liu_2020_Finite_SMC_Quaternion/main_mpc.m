clearvars;
clear global;
close all;
clc;

% Parameters
param;
import casadi.*

%% Initial Conditions
qv_t0 = [0;0;0];
omega_t0 = [0.01;0.01;0.01];
r_t0 = [1; 1; 1] * 7.078e8;
v_t0 = [2; 3; -2] * 1e4;
Xt0 = [qv_t0; omega_t0; r_t0; v_t0];

qv0 = [sin(deg2rad(19.9984)/2);
       sin(deg2rad(-9.9987)/2);
       sin(deg2rad(15.0050)/2)];
omega0 = [0;0;0];
r0 = [10;10;-10];
v0 = [0;0;0];

% Initial State Vector
X0 = [qv0; omega0; r0; v0];

%% Initial trajectory generation
disp('generating initial trajectory');

% CONSISTENT TIMESTEP - this is crucial!
dt = 0.1;  % Use same timestep throughout
t = 150;   % simulation time
tspan = (0:dt:t)';

% Generate target trajectory
[target_state,q_t,omega_t] = rk4t(@EOM_Target,tspan,Xt0);

% Run simulation with consistent timestep
[X,u,s] = rk4(@EOM_RvD_fixed,tspan,X0, target_state);

disp('done');

%% Plots
% Relative position
figure(1);
subplot(2,2,1);
plot(tspan, X(:,7), 'r', 'LineWidth',1.5); hold on;
plot(tspan, X(:,8), 'g', 'LineWidth',1.5);
plot(tspan, X(:,9), 'b', 'LineWidth',1.5); grid on;
title('Relative Position'); xlabel('Time [s]'); ylabel('Position Error [m]');
legend('x','y','z');

% Relative velocity
subplot(2,2,2);
plot(tspan, X(:,10), 'r', 'LineWidth',1.5); hold on;
plot(tspan, X(:,11), 'g', 'LineWidth',1.5);
plot(tspan, X(:,12), 'b', 'LineWidth',1.5); grid on;
title('Relative Velocity'); xlabel('Time [s]'); ylabel('Velocity Error [m/s]');
legend('x','y','z');

% Quaternion vector part
subplot(2,2,3);
plot(tspan, X(:,1), 'r', 'LineWidth',1.5); hold on;
plot(tspan, X(:,2), 'g', 'LineWidth',1.5);
plot(tspan, X(:,3), 'b', 'LineWidth',1.5); grid on;
title('Quaternion Vector Part Error'); xlabel('Time [s]'); ylabel('q_v');
legend('q_1','q_2','q_3');

% Relative angular velocity
subplot(2,2,4);
plot(tspan, X(:,4), 'r', 'LineWidth',1.5); hold on;
plot(tspan, X(:,5), 'g', 'LineWidth',1.5);
plot(tspan, X(:,6), 'b', 'LineWidth',1.5); grid on;
title('Relative Angular Velocity'); xlabel('Time [s]'); ylabel('Angular Velocity [rad/s]');
legend('x','y','z');

% Control plots
figure(2);
subplot(2,2,1);
plot(tspan, u(:,1), 'r', 'LineWidth',1.5); hold on;
plot(tspan, u(:,2), 'g', 'LineWidth',1.5);
plot(tspan, u(:,3), 'b', 'LineWidth',1.5); grid on;
title('Control Torque \tau_r'); xlabel('Time [s]'); ylabel('\tau_r [Nm]');
legend('\tau_{r1}','\tau_{r2}','\tau_{r3}');

subplot(2,2,2);
plot(tspan, u(:,4), 'r', 'LineWidth',1.5); hold on;
plot(tspan, u(:,5), 'g', 'LineWidth',1.5);
plot(tspan, u(:,6), 'b', 'LineWidth',1.5); grid on;
title('Control Force \tau_t'); xlabel('Time [s]'); ylabel('\tau_t [N]');
legend('\tau_{t1}','\tau_{t2}','\tau_{t3}');

subplot(2,2,3);
plot(tspan, s(:,1), 'r', 'LineWidth',1.5); hold on;
plot(tspan, s(:,2), 'g', 'LineWidth',1.5);
plot(tspan, s(:,3), 'b', 'LineWidth',1.5); grid on;
title('Sliding variable for Rotation'); xlabel('Time [s]'); ylabel('Sliding Variable (S)');
legend('x','y','z');

subplot(2,2,4);
plot(tspan, s(:,4), 'r', 'LineWidth',1.5); hold on;
plot(tspan, s(:,5), 'g', 'LineWidth',1.5);
plot(tspan, s(:,6), 'b', 'LineWidth',1.5); grid on;
title('Sliding variable for Translation'); xlabel('Time [s]'); ylabel('Sliding Variable (S)');
legend('x','y','z');

% Fixed EOM_RvD function
function [X_dot,u,s] = EOM_RvD_fixed(t,X, target_state_traj)
    %including parameters
    param;
    
    %% Extract relative states from vector
    qv(1:3,1)        = X(1:3,1);   % quaternion vector
    omega(1:3,1)     = X(4:6,1);   % angular velocity error
    r(1:3,1)         = X(7:9,1);   % postion
    v(1:3,1)         = X(10:12,1); % velocity
    
    %% Extract target states from vector
    dt = 0.1;  % MATCH YOUR SIMULATION TIMESTEP
    idx = min(size(target_state_traj,1), round(t/dt) + 1);  % safeguard last point
    target_state = target_state_traj(idx, :).';
    
    qv_t = target_state(1:3);
    omega_t = target_state(4:6);
    r_t = target_state(7:9);
    v_t = target_state(10:12);
    
    %% Target satellite state
    L_t = J_t*omega_t;  % target angular momentum
    q0_t = sqrt(max(1 - norm(qv_t)^2, 0));    % scalar component of target quaternion vector
    R_t = rotation_quaternion(q0_t,qv_t);   % rotation matrix of target satellite
    
    r_pt = r_t + sigma_t;      %  position of target docking port 
    v_pt = v_t + skew(omega_t)*sigma_t;   % velocity of target docking port
    
    r_tb = R_t * r_t;   % target position in inertial frame
    v_tb = R_t * v_t;   % target velocity in inertial frame
    
    %% Relative Quaternion vector
    q0 = sqrt(max(1 - norm(qv)^2, 0));
    q_cross = skew(qv);  
    
    %% Relative dynamics term
    R = rotation_quaternion(q0,qv);               % relative rotation matrix of target w.r.t chaser satellite
    omega_t_d = inv(J_t)*cross(L_t,omega_t);      % target angular acceleration
    v_t_d = cross(v_t,omega_t);                   % target translational acceleration
    
    Cr = J_c*skew(R*omega_t) + skew(R*omega_t)*J_c - skew(J_c*(omega + R*omega_t));
    Ct = skew(omega + R*omega_t);
    
    nr = skew(R*omega_t)*J_c*R*omega_t + J_c*R*omega_t_d;
    nt = skew(R*omega_t)*R*v_t + R*v_t_d + skew(omega)*R*skew(sigma_t)*omega_t - R*skew(sigma_t)*omega_t_d;
    
    %% Chaser Satellite states
    R_c = R_t * transpose(R);    % rotation matrix of chaser satellite
    q_c = rotm2quat(R_c);        % chaser satellite quaternion 4x1 vector
    omega_c = omega + R*omega_t;   % chaser angular velocity
    r_c = r + R*r_pt;   % chaser position in its own body frame
    v_c = v + R*v_pt;   % chaser velocity in its own body frame
    
    r_cb = R_c * r_c;   % chaser position in inertial frame
    v_cb = R_c * v_c;   % chaser velocity in inertial frame
    
    %% External force and torque acting on chaser and target satellite in chaser frame
    F_g_t = m_t * compute_gravity_acceleration(R*r_t, r_tb);      % gravity force of target satellite in chaser body frame
    F_g_c = m_c * compute_gravity_acceleration(r_c, r_cb);      % gravity force of chaser satellite in chaser body frame
    
    F_j2_t = m_t * transpose(R) * computeJ2Acceleration(r_tb);    % j2 force of target satellite in chaser frame
    F_j2_c = m_c * transpose(R_c) * computeJ2Acceleration(r_cb);    % j2 force on chaser satellite in chaser frame
    
    tau_gg_t = J_t * computeTauGG(R*r_t, r_tb, J_t);   % gravity gradient torque on target satellite in chaser frame  
    tau_gg_c = J_c * computeTauGG(r_c, r_cb, J_c);   % gravity gradient torque on chaser satellite in chaser frame  
    
    % MPC control
    [u, ~] = mpc_controller(t, X, target_state_traj);
    tau_mpc = u(1:3);
    F_mpc = u(4:6);
    
    % Sliding surfaces for logging
    s1 = omega + k1*qv + k2*tanh(qv);
    s2 = v     + k3*r  + k4*tanh(r);
    s  = [s1; s2];
    
    %% Calculation
    X_dot(1:3,1) = 0.5 * (q0 * eye(3) + q_cross) * omega;  % qv_dot
    X_dot(4:6,1) = inv(J_c) * (-Cr*omega - nr + tau_mpc + tau_gg_c - tau_gg_t);    % omega_dot
    X_dot(7:9,1) = v - cross(omega, r);  % r_dot
    X_dot(10:12,1) = (1/m_c) * (-m_c*Ct*v - m_c*nt + F_mpc + F_g_c - F_g_t + F_j2_c - F_j2_t);   % v_dot
    
end