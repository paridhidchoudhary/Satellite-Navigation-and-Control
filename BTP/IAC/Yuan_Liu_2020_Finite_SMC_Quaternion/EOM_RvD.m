function [X_dot,u,s] = EOM_RvD(t,X, target_state_traj)

%including parameters
param;

%% Extract relative states from vector

qv(1:3,1)        = X(1:3,1);   % quaternion vector
omega(1:3,1)     = X(4:6,1);   % angular velocity error
r(1:3,1)         = X(7:9,1);   % postion
v(1:3,1)         = X(10:12,1); % velocity




%% Extract target states from vector
dt = 0.01;  % must match your simulation timestep
idx = min(size(target_state_traj,1), round(t/dt) + 1);  % safeguard last point
target_state = target_state_traj(idx, :).';

qv_t = target_state(1:3);
omega_t = target_state(4:6);
r_t = target_state(7:9);
v_t = target_state(10:12);





%% Target satellite state
L_t = J_t*omega_t;  % target angular momentum
q0_t = sqrt(1 - norm(qv_t)^2);    % scalar component of target quaternion vector
R_t = rotation_quaternion(q0_t,qv_t);   % rotation matrix of target satellite

r_pt = r_t + sigma_t;      %  position of target docking port 
v_pt = v_t + skew(omega_t)*sigma_t;   % velocity of target docking port


r_tb = R_t * r_t;   % target position in inertial frame
v_tb = R_t * v_t;   % target velocity in inertial frame




%% Relative Quaternion vector
q0 = sqrt(1 - norm(qv)^2);
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



% %% Sliding surfaces
% s1 = omega + k1 * qv + k2 * tanh(qv);
% s2 = v + k3 * r + k4 * tanh(r);
% 
% s = [s1 ; s2];





% %% Control law (finite time SMC)
% epsilon_r = -k5 * tanh(s1) - k6 * s1;
% epsilon_t = -k7 * tanh(s2) - k8 * s2;
% 
% epsilon_r = min(max(epsilon_r,-2),2);
% epsilon_t = min(max(epsilon_t,-10),10);
% 
% u = [epsilon_r ; epsilon_t];

% --------- MPC control --------------------------------------------------
[u, ~] = mpc_controller(t, X, target_state_traj);
tau_mpc = u(1:3)
F_mpc = u(4:6)


% (Optional) keep sliding surfaces for logging only
s1 = omega + k1*qv + k2*tanh(qv);
s2 = v     + k3*r  + k4*tanh(r);
s  = [s1; s2];




%% Calculation----------------------------------------------------------------------------
X_dot(1:3,1) = 0.5 * (q0 * eye(3) + q_cross) * omega;  % qv_dot

X_dot(4:6,1) = inv(J_c) * (-Cr*omega - nr + tau_mpc + tau_gg_c - tau_gg_t);    % omega_dot

X_dot(7:9,1) = v - cross(omega, r);  % r_dot

X_dot(10:12,1) = (1/m_c) * (-m_c*Ct*v - m_c*nt + F_mpc + F_g_c - F_g_t + F_j2_c - F_j2_t);   % v_dot

%disp(t);

end