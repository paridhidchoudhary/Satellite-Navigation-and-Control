function G = G_matrix(eta)
Q=eta(1:3);
B=eta(4:6);

th = norm(Q);
A = eye(3,3) + (1/2)*skew(Q) + skew(Q) - ((1+cos(th))/(2*th*sin(th)))*(skew(Q)^2);
S = eye(3,3) + ((1-cos(th))/(th^2))*skew(Q) + ((th-sin(th))/(2*th*sin(th)))*(skew(Q)^2);
T = (1/2)*(skew(S*B)*A) + ((1/(th^2))-((1+cos(th))/(2*th*sin(th))))*((Q*B')+(Q'*B)*A) - (((1+cos(th))*(th-sin(th)))/(2*th*((sin(th))^2)))*S*(B*Q') + ((((1+cos(th))*(th+sin(th)))/(2*(th^3)*((sin(th))^2)))-(2/(th^4)))*(Q'*B)*(Q*Q');

G = [ A   zeros(3,3) ;...
      T       A      ];