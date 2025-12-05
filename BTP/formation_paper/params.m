m= [56.7, 56.7, 56.7]; %mass in kg
J= diag([4.85, 5.10, 4.76]); %moment of inertia in kg*m^2
%I=[J, zeros(3); zeros(3), m*eye(3)];

I = zeros(6, 6, 3);
I0= zeros(6,6);
for i = 1:3
    I(:,:,i) = [J, zeros(3); zeros(3), m(i)*eye(3)];
end
I0=I(:,:,1);
hf=[-1,0,0,0;0,-1,0,0;0,0,1,0;0,0,0,1];

%control gains
c1=1.3;
c2=0.6;
p1=0.0069;
p2=0.0275;
k1=0.96;
k2=0.96;
C= diag([c1*ones(1,3), c2*ones(1,3)]);
P= diag([p1*ones(1,3), p2*ones(1,3)]);
K= diag([k1*ones(1,3), k2*ones(1,3)]);

%max control force
phi_m= m*0.05;

%Simulation parmeters
dt=0.1;
%T=718*60;
T=0.5*718*60;
N=T/dt;

