function f_dyn = f_dyn()
import casadi.*

% Constants
Jc = diag([6.0833, 1.5, 6.0833]) * 1e3;  % Inertia of chaser
mc = 4000;                              % Mass of chaser
r_pc_c = [0; 0; -1.75];                 % Effector position in chaser frame
r_pt_t = [1.1404; 3.3462; 5.8907];      % Rock position in target frame

% Define symbolic variables
Rct = SX.sym('Rct', 3, 3);        % Relative rotation matrix
omega = SX.sym('omega', 3);       % Angular velocity error
rho = SX.sym('rho', 3);           % Position error (non-center-of-mass)
rho_dot = SX.sym('rho_dot', 3);   % Velocity of error
xi = SX.sym('xi', 3);             % Integral of position error
x = [reshape(Rct, 9, 1); omega; rho; rho_dot; xi];

fps = SX.sym('fps', 3*8);         % Control inputs: 8 thrusters × 3D
u = fps;

% Thruster parameters
r_tk = getThrusterPositions();    % [3 x 8]
R_kc = repmat(eye(3), [1, 1, 8]); % Assume R_k/c = Identity for simplicity

% Compute total force and torque from thrusts
fc = SX.zeros(3,1);
tau_c = SX.zeros(3,1);
for k = 1:8
    fk = fps(3*(k-1)+1:3*k);           % Force from thruster k
    fc = fc + fk;                      % Total force
    tau_c = tau_c + cross(r_tk(:,k), fk); % Torque
end

% Dynamics
skew = @(v) [ 0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0 ];
Rct_dot = Rct * skew(omega);                % Relative attitude kinematics
omega_dot = Jc \ (tau_c - skew(omega)*Jc*omega); % Attitude dynamics

% Translational dynamics
accel_rot = cross(omega, cross(omega, Rct * r_pc_c)) + ...
            cross(omega_dot, Rct * r_pc_c);
rho_ddot = fc/mc + accel_rot;  % Kinematically coupled translational acc.

xi_dot = rho;

% Return time derivative of state
x_dot = [reshape(Rct_dot, 9, 1); omega_dot; rho_dot; rho_ddot; xi_dot];
f_dyn = Function('f_dyn', {x, u}, {x_dot});
end

function r_tk = getThrusterPositions()
% Returns 3x8 matrix of thruster positions in chaser frame
% Assuming cuboid (1.5m × 4m × 1.5m), thrusters at corners of top/bottom faces

% Coordinates of corners (simplified symmetric model)
x = [0.75, -0.75];
y = [2.0, -2.0];
z = [0.75, -0.75];

r_tk = [];
for i = 1:2
    for j = 1:2
        for k = 1:2
            r_tk = [r_tk, [x(i); y(j); z(k)]];
        end
    end
end
end
