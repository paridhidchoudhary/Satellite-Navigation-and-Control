function I_not_epsilon_dot= dyn_leader(m,v,j,omega,mu,r,b)
epsilon=[omega;v];
p=r'*b;
J=0.5*trace(j)*eye(3)+j;
f_g=-(m*mu/norm(b)^3)*p - 3*(mu/norm(b)^5)*J*p +7.5*(mu*(p^jp/norm(b)^7)*p);
m_g=3*(mu/norm(b)^5)*(cross_pdt(p)*j*p);
phi_g=[m_g;f_g];
I_not=[j 0;0 m*eye(3)];
adj=[cross_pdt(omega) 0; cross_pdt(v) cross_pdt(omega)];
adj_t=adj';
I_not_epsilon_dot= adj_t*I_not*epsilon + phi_g; %equation 15
end