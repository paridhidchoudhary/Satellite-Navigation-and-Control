function [X_dot,Tc,sms] = coupledEOM_smc_only(t,X)
 
%including parameters
param;
Er(1:6,1) = X(1:6,1);
N(1:6,1) = X(7:12,1);
E0(1:6,1) = X(13:18,1);
E(1:6,1) = X(19:24,1);
E0_dot(1:6,1) = X(25:30,1);
g01(1:4,1) = X(31:34,1);
g02(1:4,1) = X(35:38,1);
g03(1:4,1) = X(39:42,1);
g04(1:4,1) = X(43:46,1);
g1(1:4,1) = X(47:50,1);
g2(1:4,1) = X(51:54,1);
g3(1:4,1) = X(55:58,1);
g4(1:4,1) = X(59:62,1);

t_s = 1;

Q = N(1:3,1);             % relative attitude
B = N(4:6,1);             % relative distance
wr = Er(1:3,1);           % relative angular velocity
vr = Er(4:6,1);           % relative translational velocity
w0 = E0(1:3,1);           % leader angular velocity
v0 = E0(4:6,1);           % leader translation velocity
w = E(1:3,1);             % follower angular velocity
v = E(4:6,1);             % follower translation velocity


mu  = 398600;  


%% Leader spacecraft Dynamics-----------------------------------------------------------

g0 = [g01 , g02 , g03 , g04];   % leader pose matrix

R0 = g0(1:3,1:3);   % Leaders Rotation matrix
b0 = g0(1:3,4);     % Leaders position
  
E0_v = [  skew(w0)     v0;...            % se(3) form of velocity     4x4
         zeros(1,3)    0 ];

% adjoint operator of follower velocity
ad_E0 = adj_se3(E0_v);
      
% co-adjoint operator of follower velocity
ad_E0x = -ad_E0';

p0 = R0' * b0;
j0 = 0.5*trace(J0)*eye(3) + J0;

Mg0 = 3*(mu/(norm(b0)^5))*skew(p0)*J0*p0;   % Gravity Gradient Moment on the follower

Fg0 = -(m0*mu/(norm(b0)^3))*p0 - 3*(mu/(norm(b0)^5))*j0*p0 + (15/2)*((mu*p0'*J0*p0)/(norm(b0)^7))*p0;  % Gravity Force on follower

Tg0 = [ Mg0 ; Fg0];       % 6x1 vector of gravitational moment and force

% mass and inertia properties of follower
P0 = [     J0       zeros(3,3);...
       zeros(3,3)   m0*eye(3)];

g0_dot = g0*E0_v;
g01_dot = g0_dot(:,1);
g02_dot = g0_dot(:,2);
g03_dot = g0_dot(:,3);
g04_dot = g0_dot(:,4);

   
%% Follower spacecraft Dynamics-----------------------------------------------------------

g = [g1 , g2 , g3 , g4];   % leader pose matrix
R = g(1:3,1:3);   % Leaders Rotation matrix
b = g(1:3,4);     % Leaders position

E_v = [  skew(w)      v ;...           % SE(3) form of velocity
        zeros(1,3)    0 ];
   
p = R' * b;
j = 0.5*trace(J)*eye(3) + J;

Mg = 3*(mu/(norm(b)^5))*skew(p)*J*p;   % Gravity Gradient Moment on the follower

Fg = -(m*mu/(norm(b)^3))*p - 3*(mu/(norm(b)^5))*j *p + (15/2)*((mu*p'*J*p)/(norm(b)^7))*p;  % Gravity Force on follower

% perturbation due to earth's oblateness J2 in ECI coordinate
aj2 = [ -((3*mu*J2*b(1)*(Re^2))/(2*norm(b)^5))*(1-(5*b(3)^2)/(norm(b)^2));...
        -((3*mu*J2*b(2)*(Re^2))/(2*norm(b)^5))*(1-(5*b(3)^2)/(norm(b)^2));...
        -((3*mu*J2*b(3)*(Re^2))/(2*norm(b)^5))*(1-(5*b(3)^2)/(norm(b)^2))];


Tg = [ Mg ; Fg+(m*R'*aj2)];       % 6x1 vector of gravitational moment and force

% mass and inertia properties of follower
P = [      J       zeros(3,3);...
       zeros(3,3)   m*eye(3)];
   
% adjoint operator of follower velocity
ad_E = adj_se3(E_v);
      
% co-adjoint operator of follower velocity
ad_Ex = -ad_E';

% Disturbance force and moment on follower spacecraft
Md = 0.005*[sin(0.1*t) cos(0.1*t) -sin(0.1*t)]';                       
Fd = 0.005*[sin(0.1*t) cos(0.1*t) -sin(0.1*t)]';

% Disturbance moment and force on follower spacecraft
Td = [Md; Fd];

g_dot = g*E_v;
g1_dot = g_dot(:,1);
g2_dot = g_dot(:,2);
g3_dot = g_dot(:,3);
g4_dot = g_dot(:,4);


 
  
%% Relative coupled spacecraft Dynamics
Nv = R6_se3(N);     % se(3) form of N

% Desired relative configuration
hf = [-1  0  0  0;...
       0 -1  0  0;...
       0  0 1  0;...
       0  0  0  1];

% h = SE3_inv(g0)*g;        % Configuration error between follower and leader spacecraft
h = hf*se3_SE3(Nv);

h_inv = SE3_inv(h);

Ad_h_inv = Ad_SE3(h_inv);


% Kinematics of exponential coordinates
th = norm(Q);
A = eye(3,3) + (1/2)*skew(Q) + skew(Q) - ((1+cos(th))/(2*th*sin(th)))*(skew(Q)^2);
S = eye(3,3) + ((1-cos(th))/(th^2))*skew(Q) + ((th-sin(th))/(2*th*sin(th)))*(skew(Q)^2);
T = (1/2)*(skew(S*B)*A) + ((1/(th^2))-((1+cos(th))/(2*th*sin(th))))*((Q*B')+(Q'*B)*A) - (((1+cos(th))*(th-sin(th)))/(2*th*((sin(th))^2)))*S*(B*Q') + ((((1+cos(th))*(th+sin(th)))/(2*(th^3)*((sin(th))^2)))-(2/(th^4)))*(Q'*B)*(Q*Q');

G = [ A   zeros(3,3) ;...
      T       A      ];


%% Sliding Mode Controller-----------------------------------------------------------------

C = [ c1*eye(3,3)   zeros(3,3) ;...
       zeros(3,3)  c2*eye(3,3)];
   
F = [ p1*eye(3,3)   zeros(3,3) ;...
       zeros(3,3)  p2*eye(3,3)];
   
K = [ k1*eye(3,3)   zeros(3,3) ;...
       zeros(3,3)  k2*eye(3,3)];
   
sms = Er + C*N;    % Sliding Surface


% control moment and force vector Tc = [Mcx Mcy Mcz Fcx Fcy Fcz]
Tc = -Tg -Td - (ad_Ex*P*E) - P*(ad_E*Ad_h_inv*E0 - Ad_h_inv*E0_dot + C*G*Er) - F*sms - K*tanh(sms);

%% Calculation----------------------------------------------------------------------------

X_dot(1:6,1) = (inv(P)*(ad_Ex*P*E + Tg + Tc + Td) + ad_E*Ad_h_inv*E0 - Ad_h_inv*E0_dot)*t_s;     % Relative velocity error  (Er = Er_dot x T) 
% X_dot(1:6,1) = (-C*G*Er - inv(P)*(F*sms + K*sign(sms)))*t_s;

X_dot(7:12,1) = (-G*C*N)*t_s;                               % Relative configuration error  (N = N_dot x T)

X_dot(13:18,1) = (inv(P0)*(ad_E0x*P0*E0 + Tg0))*t_s;         % Leader unified velocity vector (E0 = E0_dot x T)

X_dot(19:24,1) = (Er + Ad_h_inv*E0 - E)*t_s;               % Follower unified velocity vector (E = E_dot x T)

X_dot(25:30,1) = ((inv(P0)*(ad_E0x*P0*E0 + Tg0))-E0_dot)*t_s;  % leader unified acceleration vector E0_dot

X_dot(31:34,1) = g01_dot * t_s;
X_dot(35:38,1) = g02_dot * t_s;
X_dot(39:42,1) = g03_dot * t_s;
X_dot(43:46,1) = g04_dot * t_s;

X_dot(47:50,1) = g1_dot * t_s;
X_dot(51:54,1) = g2_dot * t_s;
X_dot(55:58,1) = g3_dot * t_s;
X_dot(59:62,1) = g4_dot * t_s;

disp(t);

end