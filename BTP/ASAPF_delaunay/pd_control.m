function [taux, tauy, tauz, T] = pd_control(vx_apf, vy_apf, psi_ref, z_ref, del, del_dot, k)
param;

taux = vx_apf - k(2)*del_dot(4) + k(3)*del(2) - k(4)*del_dot(2)
tauy = vy_apf - k(6)*del_dot(5) + k(7)*del(1) - k(8)*del_dot(1)
% phi_des = -vy_apf / g;
% theta_des =  vx_apf / g;
% taux = k(1)*(phi_des - del(1)) - k(2)*del_dot(4);
% tauy = k(3)*(theta_des - del(2)) - k(4)*del_dot(5);

tauz = -k(9)*(del(3) - psi_ref) - k(10)*del_dot(3);
%T = k(11)*(z_ref - del(6)) - k(12)*del_dot(6);
T = m*g;
end
