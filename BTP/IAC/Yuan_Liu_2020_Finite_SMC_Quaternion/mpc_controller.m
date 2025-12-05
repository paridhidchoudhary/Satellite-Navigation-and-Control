% 
function [u, solver_stats] = mpc_controller(t, X, target_state_traj)
% NMPC for 6-DOF rendezvous using full target state propagation

% USER PARAMETERS
H        = 2.0;      % prediction horizon [s]
dt_ctrl  = 0.2;      % controller step [s]
Q_att    =  20 * eye(3);
Q_omega  =  10 * eye(3);
Q_pos    = 500 * eye(3);
Q_vel    =   5 * eye(3);
R_tau    =   1 * eye(3);
R_F      =  10 * eye(3);
u_max    = [5; 5; 5; 20; 20; 20];

import casadi.*

persistent opti Xv Uv X0 TargetHorizon N dt last_dt prev_U

if isempty(opti) || abs(dt_ctrl - last_dt) > 1e-9
    N       = round(H / dt_ctrl);
    dt      = dt_ctrl;
    last_dt = dt;

    opti  = Opti();

    Xv    = opti.variable(12, N+1);      % [qv; ω; r; v] trajectory
    Uv    = opti.variable(6, N);         % [τ;  F] trajectory
    X0    = opti.parameter(12, 1);       % current state
    TargetHorizon = opti.parameter(12, N+1); % full target state trajectory

    Q  = blkdiag(Q_att, Q_omega, Q_pos, Q_vel);
    R  = blkdiag(R_tau, R_F);

    J = 0;
    for k = 1:N
        err = Xv(:,k);  % desired relative state is zero
        J   = J + err.'*Q*err + Uv(:,k).' * R * Uv(:,k);

        % Use full target state for each horizon step
        f = @(x, u, target) dynamics_rel(x, u, target);
        Xv(:,k+1) = Xv(:,k) + dt * f(Xv(:,k), Uv(:,k), TargetHorizon(:,k));
    end
    errN = Xv(:,N+1);
    J    = J + errN.' * Q * errN;

    opti.minimize(J);
    opti.subject_to(Xv(:,1) == X0);
    opti.subject_to(-u_max <= Uv <= u_max);

    opti.solver('ipopt', struct('print_time', false), ...
                         struct('print_level', 0, 'sb', 'yes'));
end

% Prepare full target state for each horizon step
dt_ref = 0.01;  % must match your simulation timestep
k0 = min(size(target_state_traj, 1), round(t / dt_ref) + 1);
target_horizon = zeros(12, N+1);
for k = 0:N
    kk = min(size(target_state_traj, 1), k0 + k);
    target_horizon(:, k+1) = target_state_traj(kk, :).';
end

opti.set_value(X0, X);
opti.set_value(TargetHorizon, target_horizon);

if ~isempty(prev_U)
    opti.set_initial(Uv, prev_U);
end

sol  = opti.solve();
Uopt = sol.value(Uv);
u    = Uopt(:,1);
prev_U = Uopt;

if nargout > 1
    solver_stats = sol.stats();
end
end


% function [u, solver_stats] = mpc_controller(t, X, target_state_traj)
% % NMPC for 6-DOF rendezvous using full target state propagation
% % USER PARAMETERS
% H = 2.0; % prediction horizon [s]
% dt_ctrl = 0.1; % controller step [s] - MATCH YOUR SIMULATION TIMESTEP
% Q_att = 20 * eye(3);
% Q_omega = 10 * eye(3);
% Q_pos = 500 * eye(3);
% Q_vel = 5 * eye(3);
% R_tau = 1 * eye(3);
% R_F = 10 * eye(3);
% u_max = [5; 5; 5; 20; 20; 20];
% 
% import casadi.*
% persistent opti Xv Uv X0 TargetHorizon N dt last_dt prev_U
% 
% if isempty(opti) || abs(dt_ctrl - last_dt) > 1e-9
%     N = round(H / dt_ctrl);
%     dt = dt_ctrl;
%     last_dt = dt;
% 
%     opti = Opti();
%     Xv = opti.variable(12, N+1); % [qv; ω; r; v] trajectory
%     Uv = opti.variable(6, N); % [τ; F] trajectory
%     X0 = opti.parameter(12, 1); % current state
%     TargetHorizon = opti.parameter(12, N+1); % full target state trajectory
% 
%     Q = blkdiag(Q_att, Q_omega, Q_pos, Q_vel);
%     R = blkdiag(R_tau, R_F);
% 
%     J = 0;
%     for k = 1:N
%         err = Xv(:,k); % desired relative state is zero
%         J = J + err.'*Q*err + Uv(:,k).' * R * Uv(:,k);
% 
%         % Use consistent dynamics function
%         f = @(x, u, target) dynamics_rel_fixed(x, u, target);
%         Xv(:,k+1) = Xv(:,k) + dt * f(Xv(:,k), Uv(:,k), TargetHorizon(:,k));
%     end
% 
%     errN = Xv(:,N+1);
%     J = J + errN.' * Q * errN;
% 
%     opti.minimize(J);
%     opti.subject_to(Xv(:,1) == X0);
%     opti.subject_to(-u_max <= Uv <= u_max);
% 
%     opti.solver('ipopt', struct('print_time', false), ...
%         struct('print_level', 0, 'sb', 'yes'));
% end
% 
% % Prepare full target state for each horizon step
% dt_ref = 0.1; % MATCH YOUR SIMULATION TIMESTEP
% k0 = min(size(target_state_traj, 1), round(t / dt_ref) + 1);
% target_horizon = zeros(12, N+1);
% for k = 0:N
%     kk = min(size(target_state_traj, 1), k0 + k);
%     target_horizon(:, k+1) = target_state_traj(kk, :).';
% end
% 
% opti.set_value(X0, X);
% opti.set_value(TargetHorizon, target_horizon);
% 
% if ~isempty(prev_U)
%     opti.set_initial(Uv, prev_U);
% end
% 
% sol = opti.solve();
% Uopt = sol.value(Uv);
% u = Uopt(:,1);
% prev_U = Uopt;
% 
% if nargout > 1
%     solver_stats = sol.stats();
% end
% end
% 
% function xdot = dynamics_rel_fixed(x, u, target)
% % Fixed dynamics function that matches EOM_RvD exactly
% % x = [qv; ω; r; v] relative state (12x1)
% % u = [τ; F] control input (6x1)
% % target = 12x1 vector, target state at this horizon step
% 
% param; % brings J_c, m_c, mu_earth, etc. into scope
% 
% % Unpack states and inputs
% qv = x(1:3); omega = x(4:6);
% r = x(7:9); v = x(10:12);
% tau_mpc = u(1:3); F_mpc = u(4:6);
% 
% % Extract target state
% qv_t = target(1:3);
% omega_t = target(4:6);
% r_t = target(7:9);
% v_t = target(10:12);
% 
% %% Target satellite state (matching EOM_RvD)
% L_t = J_t*omega_t;  % target angular momentum
% q0_t = sqrt(max(1 - norm(qv_t)^2, 0));    % scalar component of target quaternion vector
% R_t = rotation_quaternion(q0_t,qv_t);   % rotation matrix of target satellite
% 
% r_pt = r_t + sigma_t;      %  position of target docking port 
% v_pt = v_t + skew(omega_t)*sigma_t;   % velocity of target docking port
% 
% r_tb = R_t * r_t;   % target position in inertial frame
% v_tb = R_t * v_t;   % target velocity in inertial frame
% 
% %% Relative Quaternion vector
% q0 = sqrt(max(1 - norm(qv)^2, 0));
% q_cross = skew(qv);  
% 
% %% Relative dynamics term (matching EOM_RvD)
% R = rotation_quaternion(q0,qv);               % relative rotation matrix of target w.r.t chaser satellite
% omega_t_d = inv(J_t)*cross(L_t,omega_t);      % target angular acceleration
% v_t_d = cross(v_t,omega_t);                   % target translational acceleration
% 
% Cr = J_c*skew(R*omega_t) + skew(R*omega_t)*J_c - skew(J_c*(omega + R*omega_t));
% Ct = skew(omega + R*omega_t);
% 
% nr = skew(R*omega_t)*J_c*R*omega_t + J_c*R*omega_t_d;
% nt = skew(R*omega_t)*R*v_t + R*v_t_d + skew(omega)*R*skew(sigma_t)*omega_t - R*skew(sigma_t)*omega_t_d;
% 
% %% Chaser Satellite states
% R_c = R_t * transpose(R);    % rotation matrix of chaser satellite
% omega_c = omega + R*omega_t;   % chaser angular velocity
% r_c = r + R*r_pt;   % chaser position in its own body frame
% v_c = v + R*v_pt;   % chaser velocity in its own body frame
% 
% r_cb = R_c * r_c;   % chaser position in inertial frame
% v_cb = R_c * v_c;   % chaser velocity in inertial frame
% 
% %% External force and torque (matching EOM_RvD)
% F_g_t = m_t * compute_gravity_acceleration(R*r_t, r_tb);      % gravity force of target satellite in chaser body frame
% F_g_c = m_c * compute_gravity_acceleration(r_c, r_cb);      % gravity force of chaser satellite in chaser body frame
% 
% F_j2_t = m_t * transpose(R) * computeJ2Acceleration(r_tb);    % j2 force of target satellite in chaser frame
% F_j2_c = m_c * transpose(R_c) * computeJ2Acceleration(r_cb);    % j2 force on chaser satellite in chaser frame
% 
% tau_gg_t = J_t * computeTauGG(R*r_t, r_tb, J_t);   % gravity gradient torque on target satellite in chaser frame  
% tau_gg_c = J_c * computeTauGG(r_c, r_cb, J_c);   % gravity gradient torque on chaser satellite in chaser frame  
% 
% %% Calculation (matching EOM_RvD exactly)
% qv_dot = 0.5 * (q0 * eye(3) + q_cross) * omega;  % qv_dot
% 
% omega_dot = inv(J_c) * (-Cr*omega - nr + tau_mpc + tau_gg_c - tau_gg_t);    % omega_dot
% 
% r_dot = v - cross(omega, r);  % r_dot
% 
% v_dot = (1/m_c) * (-m_c*Ct*v - m_c*nt + F_mpc + F_g_c - F_g_t + F_j2_c - F_j2_t);   % v_dot
% 
% xdot = [qv_dot; omega_dot; r_dot; v_dot];
% end