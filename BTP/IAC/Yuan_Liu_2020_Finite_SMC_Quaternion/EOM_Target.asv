function [X_dot,qv_t,omega_t] = EOM_Target(t,X)

%including parameters
param;

%% Extract states from vector

qv_t(1:3,1)        = X(1:3,1);   % quaternion vector
omega_t(1:3,1)     = X(4:6,1);   % angular velocity error
r_t(1:3,1)         = X(7:9,1);   % postion
v_t(1:3,1)         = X(10:12,1); % velocity


L_t = J_t*omega_t; % Target angular momentum


%% Target control torques / forces (assumed zero here, as in paper)
tau_c_t = [0;0;0]; % no target control torque
F_c_t = [0;0;0];   % no target control force




%% Quaternion kinematics
q0_t = sqrt(1 - norm(qv_t)^2);
q_cross = skew(qv_t);


R_t = rotation_quaternion(q0_t,qv_t);   % rotation matrix of target satellite


r_tb = R_t * r_t;   % target position in inertial frame
v_tb = R_t * v_t;   % target velocity in inertial frame



%% External Force and Torque acting on target satellite
F_g_t = m_t * compute_gravity_acceleration(r_t, r_tb);      % gravity force in target body frame
F_j2_t = m_t * transpose(R_t) * computeJ2Acceleration(r_tb);    % j2 force in target body frame

tau_gg_t = J_t * computeTauGG(r_t, r_tb, J_t);   % gravity gradient torque in target body frame  


%% Calculation----------------------------------------------------------------------------
X_dot(1:3,1) = 0.5 * (q0_t * eye(3) + q_cross) * omega_t;  % qv_dot

X_dot(4:6,1) = inv(J_t) * (cross(L_t,omega_t) + tau_c_t + tau_gg_t);    % omega_dot

X_dot(7:9,1) = v_t - cross(omega_t, r_t);  % r_dot

X_dot(10:12,1) = (1/m_t) * (m_t*cross(v_t,omega_t) + F_c_t + F_g_t + F_j2_t);   % v_dot

disp(t); 


end