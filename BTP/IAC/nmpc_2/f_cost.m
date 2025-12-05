function f_cost = f_cost()
import casadi.*

% Symbolic variables
Rct = SX.sym('Rct', 3, 3);
omega = SX.sym('omega', 3);
rho = SX.sym('rho', 3);
rho_dot = SX.sym('rho_dot', 3);  % Not used directly in cost
xi = SX.sym('xi', 3);
x = [reshape(Rct, 9, 1); omega; rho; rho_dot; xi];

fps = SX.sym('fps', 3*8);
u = fps;

% Weights
QR = 70; 
Qw = 5 * eye(3);
Qr = 0.01 * eye(3);
Qxi = 1e-6 * eye(3);
Wf = eye(24);
D = eye(3);

% Torque penalty (optional if torque computed inside f_dyn)
tau_c = SX.zeros(3,1);
r_tk = getThrusterPositions();
for k = 1:8
    fk = fps(3*(k-1)+1:3*k);
    tau_c = tau_c + cross(r_tk(:,k), fk);
end
Wt = 9e13 * eye(3);

% Cost terms
cost = QR * trace(D - D * Rct) + ...
       omega' * Qw * omega + ...
       rho' * Qr * rho + ...
       xi' * Qxi * xi + ...
       tau_c' * Wt * tau_c + ...
       fps' * Wf * fps;

f_cost = Function('f_cost', {x, u}, {cost});
end

function r_tk = getThrusterPositions()
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
