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


%% NMPC Problem
N = 8; Ts = 60;
X = opti.variable(21, N+1);
U = opti.variable(24, N);
x0_param = opti.parameter(21,1);
omega_cce_param = opti.parameter(3,1);
omega_tte_param = opti.parameter(3,1);
thruster_param = opti.parameter(3,8);

% Initialize optimizer with feasible values
opti.set_initial(X, zeros(21,N+1));
opti.set_initial(U, zeros(24,N));

total_cost = 0;

for k = 1:N
    xk = X(:,k); uk = U(:,k);
    % Extract state components:
    R_ct = reshape(xk(1:9), 3, 3);
    omega_ct = xk(10:12);
    rho_ct = xk(13:15);
    rhodot_ct = xk(16:18);
    xi = xk(19:21);

    % Controls
    fps = reshape(uk, 3, 8);
    
    % Ensure inputs are bounded to avoid numerical issues
    for j = 1:8
        opti.subject_to(norm_2(fps(:,j)) <= fmax);
    end
    
    % Now build expressions:
    F_total = sum(fps, 2);
    Tau_total = zeros(3,1);
    for i = 1:8
        Tau_total = Tau_total + cross_casadi(thruster_param(:,i), fps(:,i));
    end

    % Create skew-symmetric matrix properly
    omega_skew = [0, -omega_ct(3), omega_ct(2);
                  omega_ct(3), 0, -omega_ct(1);
                  -omega_ct(2), omega_ct(1), 0];

    R_dot = R_ct * omega_skew;

    % Fix the rotation matrix dynamics
    % Add orthogonality constraint for R_ct
    for i = 1:3
        for j = 1:3
            if i == j
                % Diagonal elements should be close to 1 or -1
                opti.subject_to(R_ct(i,j)^2 <= 1.001);
                opti.subject_to(R_ct(i,j)^2 >= 0.999);
            else
                % Off-diagonal elements should be small
                opti.subject_to(R_ct(i,j) <= 1.001);
                opti.subject_to(R_ct(i,j) >= -1.001);
            end
        end
    end

    % Properly handle angular velocity dynamics
    omega_dot = inv(Jc) * (cross_casadi(Jc*omega_cce_param, omega_cce_param) + Tau_total) - ...
                inv(Jt) * (cross_casadi(Jt*omega_tte_param, omega_tte_param)) + ...
                cross_casadi(omega_ct, R_ct * omega_cce_param);

    rot_pc = R_ct * r_pc_c;

    % Translational dynamics
    rho_ddot = F_total/m_c + cross_casadi(omega_dot, rot_pc) + cross_casadi(omega_ct, cross_casadi(omega_ct, rot_pc));

    xi_dot = rho_ct;

    % State derivative
    x_dot = [reshape(R_dot,9,1); omega_dot; rhodot_ct; rho_ddot; xi_dot];
    
    % Define cost function
    QR = 70*eye(3); 
    Qomega = 5*eye(3); 
    Qrho = 10*eye(3); 
    Qxi = 1e-3*eye(3); 
    Wtau = 9e13*eye(3); 
    Wf = eye(24);
    
    % Cost function components
    F_cost = trace(QR*(eye(3) - R_ct)) + ...
             omega_ct'*Qomega*omega_ct + ...
             rho_ct'*Qrho*rho_ct + ...
             xi'*Qxi*xi + ...
             Tau_total'*Wtau*Tau_total + ...
             uk'*Wf*uk;
    
    P = blkdiag(Qrho, Qrho);
    
    % Euler integration for dynamics
    x_next = xk + Ts * x_dot;
    
    % Add dynamics constraint
    opti.subject_to(X(:,k+1) == x_next);
    
    % Accumulate stage cost
    total_cost = total_cost + F_cost;
end

% Terminal cost
rho_terminal = X(13:15,N+1);
rhodot_terminal = X(16:18,N+1);
E_cost = [rho_terminal; rhodot_terminal]' * blkdiag(Qrho, Qrho) * [rho_terminal; rhodot_terminal];
total_cost = total_cost + E_cost;

% Set objective
opti.minimize(total_cost);

% Initial constraint
opti.subject_to(X(:,1) == x0_param);

% Solver options
opts = struct;
opts.ipopt.print_level = 5;
opts.print_time = false;
opts.ipopt.max_iter = 500;
opts.ipopt.tol = 1e-4;  % Relaxed tolerance
opts.ipopt.acceptable_tol = 1e-3;  % More relaxed acceptable tolerance
opts.ipopt.hessian_approximation = 'limited-memory';  % Use BFGS approximation instead of exact Hessian
opti.solver('ipopt', opts);

%% Simulation Setup
T_final = 800; dt = 1; t = 0:dt:T_final; N_steps = length(t);
R_c_t = repmat(eye(3),1,1,N_steps); rho_c_t = zeros(3,N_steps); rhodot_c_t = zeros(3,N_steps);
omega_c_t_t = zeros(3,N_steps); omega_c_c_e = zeros(3,N_steps); omega_t_t_e = repmat(0.0046*[1;1;1],1,N_steps);
xi = zeros(3,N_steps); thruster_forces = zeros(3,8,N_steps);
rho_c_t(:,1) = [-1.1394; 113.0847; -7.6204]; rhodot_c_t(:,1) = [0.0081; 0.0039; 1.3490];

% Initial orientations
R_t_e_0 = eye(3);
R_c_e_0 = diag([-1, -1, 1]);
R_c_t_0 = R_c_e_0 * R_t_e_0';  % Corrected multiplication order

% Initial angular velocities
omega_c_c_e_0 = [0; 0; 0];
omega_t_t_e_0 = 0.0046 * [1; 1; 1];

% Compute and assign relative angular velocity in target frame
omega_c_t_t_0 = R_c_t_0 * omega_c_c_e_0 - omega_t_t_e_0;
omega_c_t_t(:,1) = omega_c_t_t_0;

% Also store these for use later
omega_c_c_e(:,1) = omega_c_c_e_0;
omega_t_t_e(:,1) = omega_t_t_e_0;

% Set R_c_t(:,:,1) correctly
R_c_t(:,:,1) = R_c_t_0;

% Verify initial rotation matrix is valid
if abs(det(R_c_t_0) - 1) > 1e-8
    warning('Initial rotation matrix is not valid: det = %f', det(R_c_t_0));
    % Fix it by SVD
    [U, ~, V] = svd(R_c_t_0);
    R_c_t(:,:,1) = U * V';
end

%% Simulation Loop
for k = 1:N_steps-1
    % Print state info for debugging
    if mod(k, 100) == 1
        fprintf('Step %d/%d\n', k, N_steps-1);
        disp('omega_c_t_t:'); disp(omega_c_t_t(:,k));
        disp('rho_c_t:'); disp(rho_c_t(:,k));
        disp('rhodot_c_t:'); disp(rhodot_c_t(:,k));
    end

    Rct_k = R_c_t(:,:,k);
    
    % Ensure rotation matrix stays orthogonal
    if mod(k, 50) == 0
        [U, ~, V] = svd(Rct_k);
        Rct_k = U * V';
        R_c_t(:,:,k) = Rct_k;
    end
    
    xk = [reshape(Rct_k,9,1); omega_c_t_t(:,k); rho_c_t(:,k); rhodot_c_t(:,k); xi(:,k)];
    
    % Solve MPC problem at specified intervals or first step
    if mod(t(k), Ts) == 0 || k == 1
        opti.set_value(x0_param, xk);
        opti.set_value(omega_cce_param, omega_c_c_e(:,k));
        opti.set_value(omega_tte_param, omega_t_t_e(:,k));
        opti.set_value(thruster_param, thruster_positions);
        
        try
            sol = opti.solve();
            u_opt = sol.value(U(:,1));
            fprintf('✅ NMPC solved at step %d\n', k);
            
            % For debugging: check solution quality
            x_opt = sol.value(X);
            fprintf('Final position error: [%.4f, %.4f, %.4f]\n', x_opt(13:15,end));
            
        catch ME
            fprintf('❌ NMPC failed at step %d: %s\n', k, ME.message);
            
            % Use zero control or previous control as fallback
            if k > 1
                u_opt = reshape(thruster_forces(:,:,k-1),[],1);
            else
                u_opt = zeros(24,1);
            end
        end
    end
    
    % Extract and apply control
    fps_k = reshape(u_opt, 3, 8); 
    thruster_forces(:,:,k) = fps_k;
    
    % Compute dynamics for simulation
    F_c = sum(fps_k, 2); 
    Tau_c = zeros(3,1);
    for j = 1:8
        Tau_c = Tau_c + cross(thruster_positions(:,j), fps_k(:,j));
    end
    
    omega_ct = omega_c_t_t(:,k);
    omega_skew = [0, -omega_ct(3), omega_ct(2); 
                  omega_ct(3), 0, -omega_ct(1); 
                  -omega_ct(2), omega_ct(1), 0];
    
    % Rotation dynamics
    R_dot = Rct_k * omega_skew;
    R_c_t(:,:,k+1) = Rct_k + dt * R_dot;
    
    % Angular velocity dynamics
    omega_dot = Rct_k * (Jc \ (cross(Jc*omega_c_c_e(:,k), omega_c_c_e(:,k)) + Tau_c)) - ...
                Jt \ (cross(Jt*omega_t_t_e(:,k), omega_t_t_e(:,k))) + ...
                cross(omega_ct, Rct_k * omega_c_c_e(:,k));
    
    omega_c_t_t(:,k+1) = omega_ct + dt * omega_dot;
    
    % Position dynamics
    rot_pc = Rct_k * r_pc_c;
    rho_ddot = F_c/m_c + cross(omega_dot, rot_pc) + cross(omega_ct, cross(omega_ct, rot_pc));
    rhodot_c_t(:,k+1) = rhodot_c_t(:,k) + dt * rho_ddot;
    rho_c_t(:,k+1) = rho_c_t(:,k) + dt * rhodot_c_t(:,k);
    xi(:,k+1) = xi(:,k) + dt * rho_c_t(:,k+1);
    
    % Ensure rotation matrix stays orthogonal after integration
    [U, ~, V] = svd(R_c_t(:,:,k+1));
    R_c_t(:,:,k+1) = U * V';
end

%% Plots
figure;
subplot(3,1,1); 
plot(t, rho_c_t); 
title('Relative Position'); 
ylabel('m'); 
legend('X', 'Y', 'Z');
grid on;

subplot(3,1,2); 
plot(t, rhodot_c_t); 
title('Relative Velocity'); 
ylabel('m/s'); 
legend('Vx', 'Vy', 'Vz');
grid on;

subplot(3,1,3); 
plot(t, vecnorm(rho_c_t)); 
title('Distance'); 
ylabel('m'); 
xlabel('Time (s)');
grid on;

figure;
plot3(rho_c_t(1,:), rho_c_t(2,:), rho_c_t(3,:), 'b-', 'LineWidth', 1.5); 
hold on;

% Start position (initial)
plot3(rho_c_t(1,1), rho_c_t(2,1), rho_c_t(3,1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
text(rho_c_t(1,1), rho_c_t(2,1), rho_c_t(3,1), ' Start', 'Color', 'green');

% Target feature point (origin)
plot3(0, 0, 0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
text(0, 0, 0, ' Target', 'Color', 'red');

xlabel('X'); ylabel('Y'); zlabel('Z'); grid on; axis equal;
title('3D Trajectory of Chaser Feature Point');
legend('Chaser Trajectory', 'Initial Position', 'Target Position');

% Attitude error plot
theta = zeros(1, N_steps);
for k = 1:N_steps
    theta(k) = acos(min(1, max(-1, 0.5 * (trace(R_c_t(:,:,k)) - 1))));
end
figure; 
plot(t, rad2deg(theta)); 
title('Attitude Error'); 
ylabel('deg'); 
xlabel('Time (s)');
grid on;

% Force magnitude plot
figure;
F_net = squeeze(sum(thruster_forces,2));  % 3 x N_steps
F_mag = vecnorm(F_net, 2, 1);             % 1 x N_steps
plot(t(1:end-1), F_mag);
title('Total Applied Force Magnitude (N)');
xlabel('Time (s)'); 
ylabel('Force (N)'); 
grid on;

% Final results
disp("Final feature point position error (m):");
disp(rho_c_t(:,end));

disp("Final velocity error (m/s):");
disp(rhodot_c_t(:,end));

%% Define cross product function for CasADi
% Using a proper function for cross product to avoid issues
function c = cross_casadi(a, b)
    c = [a(2)*b(3) - a(3)*b(2);
         a(3)*b(1) - a(1)*b(3);
         a(1)*b(2) - a(2)*b(1)];
end
