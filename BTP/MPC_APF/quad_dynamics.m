function state_dot = quad_dynamics(state, taux, tauy, tauz, T)
param;
%disp('into quad_dynamics\n');
phi = state(1);
theta = state(2);
psi = state(3);
phi_d = state(7);
theta_d = state(8);
psi_d = state(9);
x_d = state(10);
y_d = state(11);
z_d = state(12);
% phi_dd = state_dot_input(1);
% theta_dd = state_dot_input(2);
% psi_dd = state_dot_input(3);

% del1_dot = phi_dd + theta_dd*sin(phi)*tan(theta) + psi_dd*cos(phi)*tan(theta);
% del2_dot = theta_dd*cos(phi) - psi_dd*sin(phi);
% del3_dot = (sin(phi)/cos(theta))*theta_dd + (cos(phi)/cos(theta))*psi_dd;
state_dot(1) = phi_d;
state_dot(2) = theta_d;
state_dot(3) = psi_d;
state_dot(4) = x_d;
state_dot(5) = y_d;
state_dot(6) = z_d;

state_dot(7) = (taux/ix) - ((iz - iy)/ix)*theta_d*psi_d;
state_dot(8) = (tauy/iy) - ((ix - iz)/iy)*phi_d*psi_d;
state_dot(9) = (tauz/iz) - ((iy - ix)/iz)*phi_d*theta_d;
state_dot(10) = - (T/m)*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi));
state_dot(11) = - (T/m)*(cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi));
state_dot(12) = G - (T/m)*cos(phi)*cos(theta);

end