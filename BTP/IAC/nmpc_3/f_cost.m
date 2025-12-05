% function f_cost = f_cost()
% import casadi.*
% 
% % ------------------ SYMBOLIC VARIABLES ------------------
% Rct     = SX.sym('Rct', 3, 3);
% omega   = SX.sym('omega', 3);
% rho     = SX.sym('rho', 3);
% rho_dot = SX.sym('rho_dot', 3);  % Optional, unused
% xi      = SX.sym('xi', 3);
% x       = [reshape(Rct, 9, 1); omega; rho; rho_dot; xi];
% 
% u       = SX.sym('u', 8);  % <-- 8 scalar thrust inputs
% 
% % ------------------ COST WEIGHTS ------------------------
% QR   = 70;
% Qw   = 5 * eye(3);
% Qr   = 1 * eye(3);
% Qxi  = 1e-6 * eye(3);
% %Wu   = 1e-6 * eye(8);    % thrust penalty
% % ---- add just below Qr ----
% Qv   = 0.05 * eye(3);     % NEW  velocity-error weight (ρ̇)
% 
% % ---- change thrust weight ----
% Wu   = 1e-1  * eye(8);    % was 1e-3 or 1e-6
% 
% 
% 
% D = eye(3);
% 
% % Compute net force from thrusters
% r_tk = getThrusterPositions();     % 3x8
% fc = SX.zeros(3,1);
% for k = 1:8
%     fk = u(3*(k-1)+1:3*k);
%     fc = fc + fk;
% end
% Wfc = 0.1 * eye(3);                % new weight on net force
% 
% 
% % ---- cost expression ----
% cost = QR * trace(D - D*Rct) ...
%      + omega' * Qw * omega ...
%      + rho'   * Qr * rho ...
%      + rho_dot' * Qv * rho_dot ...   % NEW
%      + xi'    * Qxi * xi ...
%      + u'     * Wu * u...
%      + fc' * Wfc * fc;  
%               % desired orientation (identity)
% 
% % % ------------------ COST FUNCTION -----------------------
% % cost = QR * trace(D - D * Rct) + ...
% %        omega' * Qw * omega + ...
% %        rho'   * Qr * rho   + ...
% %        xi'    * Qxi * xi   + ...
% %        u'     * Wu * u;      % thrust cost
% 
% % ------------------ CASADI FUNCTION ---------------------
% f_cost = Function('f_cost', {x, u}, {cost});
% end

function f_cost = f_cost()
% Cost function for NMPC with 8 scalar thrusters (u : 8×1)

import casadi.*

% ------------------------------------------------------------------
%                         SYMBOLIC STATE
% ------------------------------------------------------------------
Rct     = SX.sym('Rct',3,3);
omega   = SX.sym('omega',3);
rho     = SX.sym('rho',3);
rho_dot = SX.sym('rho_dot',3);
xi      = SX.sym('xi',3);
x       = [reshape(Rct,9,1); omega; rho; rho_dot; xi];

u  = SX.sym('u',8);          % 8 scalar thrust magnitudes

% ------------------------------------------------------------------
%                         WEIGHTS (tune here)
% ------------------------------------------------------------------
QR  = 70;                    % attitude trace weight
Qw  = 5    * eye(3);         % angular velocity
Qr  = 1    * eye(3);         % position  (↑ from 0.01)
Qv  = 0.05 * eye(3);         % velocity  (new, damps oscill.)
Qxi = 1e-6 * eye(3);         % integral term

Wu  = 1e-1  * eye(8);        % thrust effort (↑ makes jets expensive)
Wfc = 1     * eye(3);        % net-force penalty  (new term)

D = eye(3);                  % desired attitude (identity)

% ------------------------------------------------------------------
%                 NET FORCE   F_c  =  B * u
% ------------------------------------------------------------------
[B, ~] = getThrusterConfig();    % B : 3×8 direction matrix
Fc     = B * u;                  % 3×1 net body force

% ------------------------------------------------------------------
%                        COST EXPRESSION
% ------------------------------------------------------------------
cost =   QR * trace(D - D*Rct)          ...
       + omega'   * Qw  * omega         ...
       + rho'     * Qr  * rho           ...
       + rho_dot' * Qv  * rho_dot       ...
       + xi'      * Qxi * xi            ...
       + u'       * Wu  * u             ...  % thruster effort
       + Fc'      * Wfc * Fc;               % net-force penalty

f_cost = Function('f_cost',{x,u},{cost});
end

