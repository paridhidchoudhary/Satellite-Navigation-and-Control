m= 56.7; %mass in kg
J = diag([50, 100 , 50]);
J0 = diag([50 , 50 , 50]);

%control gains
c1=0.96;
c2=0.96;
p1=0.01;
p2=0.02;
k1=0.96;
k2=0.96;
C= diag([c1*ones(1,3), c2*ones(1,3)]);
P= diag([p1*ones(1,3), p2*ones(1,3)]);
K= diag([k1*ones(1,3), k2*ones(1,3)]);

%max control force
phi_f = 10;
phi_m = 20;



