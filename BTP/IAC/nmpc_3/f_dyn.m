function f_dyn = f_dyn()
import casadi.*

% -------------------------------------------------------------------------
%                           CONSTANTS  (unchanged)
% -------------------------------------------------------------------------
Jc = diag([6.0833, 1.5, 6.0833]) * 1e3;
mc = 4000;
r_pc_c = [0; 0; -1.75];
r_pt_t = [1.1404; 3.3462; 5.8907];

% -------------------------------------------------------------------------
%                           SYMBOLIC STATES  (unchanged)
% -------------------------------------------------------------------------
Rct     = SX.sym('Rct',3,3);
omega   = SX.sym('omega',3);
rho     = SX.sym('rho',3);
rho_dot = SX.sym('rho_dot',3);
xi      = SX.sym('xi',3);
x       = [reshape(Rct,9,1); omega; rho; rho_dot; xi];

% -------------------------------------------------------------------------
%                           ** NEW CONTROL INPUT **
% -------------------------------------------------------------------------
u  = SX.sym('u',8);        % 8 scalar thruster magnitudes  (was 24-vector)
[B, r_tk] = getThrusterConfig();   % 3×8 unit dirs  & 3×8 positions

% -------------------------------------------------------------------------
%                           FORCE & TORQUE
% -------------------------------------------------------------------------
Fc = B * u;                                 % total force
Tc = SX.zeros(3,1);
for k = 1:8
    Tc = Tc + cross(r_tk(:,k), B(:,k)*u(k));   % moment of thruster k
end

% -------------------------------------------------------------------------
%                           KINEMATICS & DYNAMICS  (unchanged)
% -------------------------------------------------------------------------
skew = @(v)[ 0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0 ];
Rct_dot  = Rct * skew(omega);
omega_dot = Jc \ (Tc - skew(omega)*Jc*omega);
accel_rot = cross(omega, cross(omega, Rct * r_pc_c)) + ...
            cross(omega_dot, Rct * r_pc_c);
rho_ddot = Fc/mc + accel_rot;
xi_dot   = rho;

x_dot = [reshape(Rct_dot,9,1); omega_dot; rho_dot; rho_ddot; xi_dot];
disp('fc (symbolic):'); disp(Fc);
f_dyn = Function('f_dyn',{x,u},{x_dot});

end

