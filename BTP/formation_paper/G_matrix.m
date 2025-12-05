function G=G_matrix(eta)
theta=eta(1:3);
beta=eta(4:6);
theta_mod=norm(theta);
A= eye(3) + 0.5*cross_pdt(theta) + ((1/theta_mod^2) - ((1+cos(theta_mod))/(2*theta_mod*sin(theta_mod)))*cross_pdt(theta)*cross_pdt(theta));
S= eye(3) +((1-cos(theta_mod))/(theta_mod^2))*cross_pdt(theta) + ((theta_mod - sin(theta_mod))/(theta_mod^3))*cross_pdt(theta)*cross_pdt(theta);
L=cross_pdt(S*beta)*A;
M=theta*beta' +theta'*beta*A;
N=S*beta*theta';
O=theta'*beta*(theta*theta');
p=(1+cos(theta_mod))*(theta_mod +sin(theta_mod))/(2*theta_mod^3 * sin(theta_mod)* sin(theta_mod));
T=0.5*L + ((1/(theta_mod^2))-((1+cos(theta_mod))/(2*theta_mod*sin(theta_mod))))*M -((1+cos(theta_mod))*(theta_mod -sin(theta_mod))/(2*theta_mod*sin(theta_mod)*sin(theta_mod)))*N + (p-(2/theta_mod^4))*O;
G=[A ,zeros(3); T, A];
end