function f_cost = f_cost()

import casadi.*

% ---------- symbolic state ----------
Rct     = SX.sym('Rct',3,3);
omega_e = SX.sym('omega_e',3);
rho     = SX.sym('rho',3);
rho_dot = SX.sym('rho_dot',3);
xi      = SX.sym('xi',3);
x       = [reshape(Rct,9,1); omega_e; rho; rho_dot; xi];

u       = SX.sym('u',24);                 % 24‑vector of forces

% ---------- weights (paper values) ----------
QR  = 70;         % attitude trace
Qw  = 5      * eye(3);
Qr  = 0.01   * eye(3);
Qxi = 1e-6   * eye(3);

Wf  = 1      * eye(24);     % thruster effort
Wtau = 9e13  * eye(3);      % torque penalty (huge)

% ---------- reshape thruster forces ----------
f_tk = reshape(u,3,8);       % 3×8, column k is f_tk

% ---------- torque τ_c = Σ r×f ----------
[r_tk, ~] = getThrusterConfig();
tau_c = SX.zeros(3,1);
for k = 1:8
    tau_c = tau_c + cross(r_tk(:,k), f_tk(:,k));
end

% ---------- cost expression ----------
D = eye(3);                  % desired Rc/t is identity
cost =  QR   * trace(D - D*Rct) ...
      + omega_e' * Qw  * omega_e ...
      + rho'     * Qr  * rho ...
      + xi'      * Qxi * xi ...
      + u'       * Wf  * u ...
      + tau_c'   * Wtau * tau_c;

f_cost = Function('f_cost',{x,u},{cost});

end
