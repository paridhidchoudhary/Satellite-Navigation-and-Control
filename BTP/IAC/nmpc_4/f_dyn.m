function f_dyn = f_dyn()
% Dynamics in the target frame
% -------------------------------------------------------------------------
%   x  = [ vec(R_ct) ;           %  9  rotation (target → chaser)
%          ω_tc/t     ;          %  3  angular‑rate error  (in Ft)
%          ρ          ;          %  3  relative position  (pc wrt pt, in Ft)
%          ρ̇         ;          %  3  relative velocity  (in Ft)
%          ξ ];                  %  3  integral of ρ
%   u  = [ f_t1 ; … ; f_t8 ]     % 24  thruster forces (each 3×1), in Fc
%   p  = ω_t                     %  3  target angular‑rate (constant or known)
% -------------------------------------------------------------------------

import casadi.*

% ---------- constants ----------
Jc   = diag([6.0833 1.5 6.0833])*1e3;     % chaser inertia  
mc   = 4000;                              % chaser mass   
r_pc_c = [0; 0; -1.75];                   % COM → end‑effector in Fc
[r_tk, j_nom]   = getThrusterConfig();        % 3×8 thruster moment arms

% ---------- symbolic variables ----------
Rct   = SX.sym('Rct',3,3);                % rotation  Ft → Fc
omega_e = SX.sym('omega_e',3);            % ω_tc/t   (error in Ft)
rho     = SX.sym('rho',3);                % position error
rho_dot = SX.sym('rho_dot',3);            % velocity error
xi      = SX.sym('xi',3);                 % integral state
x    = [reshape(Rct,9,1); omega_e; rho; rho_dot; xi];

u    = SX.sym('u',24);                    % 8×3 forces stacked
omega_t = SX.sym('omega_t',3);            % target angular‑rate 

% ---------- unpack thruster forces ----------
f_tk = reshape(u,3,8);                   

% ---------- total force & torque on chaser ----------
Fc   = sum(f_tk,2);                     
Tc   = SX.zeros(3,1);
for k = 1:8
    Tc = Tc + cross(r_tk(:,k), f_tk(:,k));
end

% ---------- recover chaser angular rate in Fc ----------
% ω_c = R_ctᵀ (ω_e + ω_t)   (paper Eq. 4 rearranged)
omega_c = Rct'*(omega_e + omega_t);

% ---------- chaser rotational dynamics  (Euler eq.) ----------
omega_c_dot = Jc \ ( Tc - cross(omega_c, Jc*omega_c) );

% ---------- error‑rate dynamics  (paper Eq. 4,) ----------
omega_e_dot = Rct*omega_c_dot + cross(omega_e, Rct*omega_c);

% ---------- rotation kinematics  (paper Eq. 3) ----------
skew = @(v)[   0   -v(3)  v(2);
             v(3)   0    -v(1);
            -v(2)  v(1)   0 ];
Rct_dot = Rct * skew(omega_e);

% ---------- translational dynamics  ----------
accel_rot = cross(omega_c, cross(omega_c, Rct*r_pc_c)) + ...
            cross(omega_c_dot, Rct*r_pc_c);
rho_ddot = Fc/mc + accel_rot;   
xi_dot   = rho;

% ---------- state derivative ----------
x_dot = [reshape(Rct_dot,9,1);
         omega_e_dot;
         rho_dot;
         rho_ddot;
         xi_dot];

% ---------- CasADi function ----------
f_dyn = Function('f_dyn',{x,u,omega_t},{x_dot});

end
