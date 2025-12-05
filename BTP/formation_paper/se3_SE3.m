function exp_delta= se3_SE3(v,omega)
theta= sqrt(omega'*omega);
A=sin(theta)/theta;
B= (1-cos(theta))/(theta^2);
C= (1-A)/(theta^2);
R= 1+ A*cross_pdt(omega) + B*(cross_pdt(omega))^2;
V= 1+ B*cross_pdt(omega) + C*(cross_pdt(omega))^2;
exp_delta= [R, V*v; 0,0,0,1];
end