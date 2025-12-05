function I_epsilon_dot= dyn_follower(m,v,j,omega,r,b, tau_c, phi_c, tau_d, phi_d)
epsilon=[omega;v];
p=r'*b;
J=0.5*trace(j)*eye(3)+j;
f_g=-(m*mu/norm(b)^3)*p - 3*(mu/norm(b)^5)*J*p +7.5*(mu*(p^jp/norm(b)^7)*p);
m_g=3*(mu/norm(b)^5)*(cross_pdt(p)*j*p);
bi=b(1);
bj=b(2);
bk=b(3);
J2=1.08263 * 10^(-3);
R=6378 * 10^(6); %in m
x=3*mu*J2*R^2/(2*norm(b)^5);
y=5*bk^2/norm(b)^2;
aj2=[-x*bi*(1-y);
-x*bj*(1-y);
-x*bk*(3-y)];
phi_g=[m_g; f_g + m*r'*aj2];
phi_c_2=[tau_c; phi_c];
phi_d_2=[tau_d; phi_d];
I=[J 0; 0 m*eye(3)];
adj=[cross_pdt(omega) 0; cross_pdt(v) cross_pdt(omega)];
adj_t=adj';
I_epsilon_dot=adj_t*I*epsilon + phi_g + phi_c_2 + phi_d_2;
end