% NMPC for Coupled Rotational-Translational Spacecraft Motion
clear all; close all; clc;
import casadi.*
opti = Opti();

%% Constants
Jc = diag([6.0833, 1.5, 6.0833]) * 1e3;
Jt = diag([300, 250, 350]);
m_c = 4000;
fmax = 20;
r_pc_c = [0; 0; -1.75];
r_pt_t = [1.1404; 3.3462; 5.8907];

thruster_positions = [... 
    [1.5/2; 4/2; 1.5/2], [1.5/2; -4/2; 1.5/2], [1.5/2; -4/2; -1.5/2], [1.5/2; 4/2; -1.5/2], ...
    [-1.5/2; 4/2; 1.5/2], [-1.5/2; -4/2; 1.5/2], [-1.5/2; -4/2; -1.5/2], [-1.5/2; 4/2; -1.5/2]
];

%% NMPC Setup
omega_cce_param = opti.parameter(3,1);
omega_tte_param = opti.parameter(3,1);
thruster_param = opti.parameter(3,8);

% %% NMPC Problem
N = 8; Ts = 60;
 
%% Simulation Setup
T_final = 800; dt = 1; t = 0:dt:T_final; N_steps = length(t);
R_c_t = repmat(eye(3),1,1,N_steps); rho_c_t = zeros(3,N_steps); rhodot_c_t = zeros(3,N_steps);
omega_c_t_t = zeros(3,N_steps); omega_c_c_e = zeros(3,N_steps); omega_t_t_e = repmat(0.0046*[1;1;1],1,N_steps);
xi = zeros(3,N_steps); thruster_forces = zeros(3,8,N_steps);
rho_c_t(:,1) = [-1.1394; 113.0847; -7.6204]; rhodot_c_t(:,1) = [0.0081; 0.0039; 1.3490];

% Initial orientations
R_t_e_0 = eye(3);
R_c_e_0 = diag([-1, -1, 1]);
R_c_t_0 = R_t_e_0' * R_c_e_0;

% Initial angular velocities
omega_c_c_e_0 = [0; 0; 0];              % Given
omega_t_t_e_0 = 0.0046 * [1; 1; 1];     % Given

% Compute and assign relative angular velocity in target frame
omega_c_t_t_0 = R_c_t_0 * omega_c_c_e_0 - omega_t_t_e_0;
omega_c_t_t(:,1) = omega_c_t_t_0;

% Also store these for use later
omega_c_c_e(:,1) = omega_c_c_e_0;
omega_t_t_e(:,1) = omega_t_t_e_0;

% Set R_c_t(:,:,1) correctly too
R_c_t(:,:,1) = R_c_t_0;


%% Simulation Loop
for k = 1:N_steps-1
    disp('omega_c_t_t(:,1):'); disp(omega_c_t_t(:,k));
    disp('rho_c_t(:,1):'); disp(rho_c_t(:,k));
    disp('rhodot_c_t(:,1):'); disp(rhodot_c_t(:,k));

    Rct_k = R_c_t(:,:,k);
    xk = [reshape(Rct_k,9,1); omega_c_t_t(:,k); rho_c_t(:,k); rhodot_c_t(:,k); xi(:,k)];
    % opti.set_initial(X, repmat(xk, 1, N+1));
    % opti.set_initial(U, zeros(24, N));

    if mod(t(k), Ts) == 0 || k == 1
        % Bundle constants for convenience
        params = struct('Jc', Jc, 'Jt', Jt, 'm_c', m_c, ...
                    'fmax', fmax, 'r_pc_c', r_pc_c, ...
                    'N', N, 'Ts', Ts);

        u_opt = solve_nmpc_step(xk, omega_c_c_e(:,k), omega_t_t_e(:,k), thruster_positions, params);
    end

    fps_k = reshape(u_opt, 3, 8); thruster_forces(:,:,k) = fps_k;
    F_c = sum(fps_k, 2); Tau_c = sum(cross(thruster_positions, fps_k, 1), 2);
    omega_ct = omega_c_t_t(:,k);
    R_dot = Rct_k * [0 -omega_ct(3) omega_ct(2); omega_ct(3) 0 -omega_ct(1); -omega_ct(2) omega_ct(1) 0];
    R_c_t(:,:,k+1) = Rct_k + dt * R_dot;
    omega_dot = Rct_k*Jc\(cross(Jc*omega_c_c_e(:,k), omega_c_c_e(:,k)) + Tau_c) - ...
                Jt\(cross(Jt*omega_t_t_e(:,k), omega_t_t_e(:,k))) + ...
                cross(omega_ct, Rct_k * omega_c_c_e(:,k));
    omega_c_t_t(:,k+1) = omega_ct + dt * omega_dot;
    rot_pc = Rct_k * r_pc_c;
    rho_ddot = F_c/m_c + cross(omega_dot, rot_pc) + cross(omega_ct, cross(omega_ct, rot_pc));
    rhodot_c_t(:,k+1) = rhodot_c_t(:,k) + dt * rho_ddot;
    rho_c_t(:,k+1) = rho_c_t(:,k) + dt * rhodot_c_t(:,k);
    xi(:,k+1) = xi(:,k) + dt * rho_c_t(:,k+1);
end


%% Plots
figure;
subplot(3,1,1); plot(t, rho_c_t); title('Relative Position'); ylabel('m'); grid on;
subplot(3,1,2); plot(t, rhodot_c_t); title('Relative Velocity'); ylabel('m/s'); grid on;
subplot(3,1,3); plot(t, vecnorm(rho_c_t)); title('Distance'); ylabel('m'); grid on;

% figure; plot3(rho_c_t(1,:), rho_c_t(2,:), rho_c_t(3,:), 'b'); grid on;
% xlabel('X'); ylabel('Y'); zlabel('Z'); title('3D Trajectory');
figure;
plot3(rho_c_t(1,:), rho_c_t(2,:), rho_c_t(3,:), 'b-', 'LineWidth', 1.5); hold on;

% Start position (initial)
plot3(rho_c_t(1,1), rho_c_t(2,1), rho_c_t(3,1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
text(rho_c_t(1,1), rho_c_t(2,1), rho_c_t(3,1), ' Start', 'Color', 'green');

% Target feature point (origin)
plot3(0, 0, 0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
text(0, 0, 0, ' Target', 'Color', 'red');

xlabel('X'); ylabel('Y'); zlabel('Z'); grid on; axis equal;
title('3D Trajectory of Chaser Feature Point');
legend('Chaser Trajectory', 'Initial Position', 'Target Position');


theta = zeros(1, N_steps);
for k = 1:N_steps
    theta(k) = acos(0.5 * (trace(R_c_t(:,:,k)) - 1));
end
figure; plot(t, rad2deg(theta)); title('Attitude Error'); ylabel('deg'); grid on;

disp("Final feature point position error (m):");
disp(rho_c_t(:,end));

disp("Final velocity error (m/s):");
disp(rhodot_c_t(:,end));

F_net = squeeze(sum(thruster_forces,2));  % 3 x N_steps
F_mag = vecnorm(F_net, 2, 1);             % 1 x N_steps
disp('f_mag size:');
disp(size(F_mag));
disp('t size:');
disp(size(t));
plot(t, F_mag);                           % both are 1 x N_steps
title('Total Applied Force Magnitude (N)');
xlabel('Time (s)'); ylabel('Force (N)'); grid on;


