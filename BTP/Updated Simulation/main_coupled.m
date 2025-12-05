clearvars;
clear global;
close all;
clc;

% parmeters;
params;


%% Initial Conditions
xi_tilde = [-0.6; 0.26; 0.5; 3.17; 1.24; 1.82];

eta_tilde = [30*pi/180; 0; 0; -1000; 1500; -500];

xi_0 = [1.5 ; 1 ; 1.2 ; 100*[2 ; 1 ; 1.5]];

xi_0_dot = zeros(6,1);

% State Vector-----------------------------------------
X0 = [xi_tilde ; eta_tilde ; xi_0 ; xi_0_dot];                  

% tf inital guess
T = 1;



%% intial trajectory generation
disp('generating initial trajectory');
%assuming linear variation of sigma
dt = 0.001;
%code run time:
t = 10;

tspan = (0:dt:t)';


%Y = initial state trajectory
[Y,phi_c,s] = rk4(@Coupled_EOM,tspan,X0);

disp('done')


%% Plots

figure(1);
subplot(4,2,1); plot(tspan,Y(:,1),tspan,Y(:,2),tspan,Y(:,3),'linewidth',2); 
xlabel('Time [sec]');
ylabel('Relative Angular Velocity (rad/s)');
legend('Wx','Wy','Wz');
grid on


subplot(4,2,2); plot(tspan,Y(:,4),tspan,Y(:,5),tspan,Y(:,6),'linewidth',2);
xlabel('Time [sec]');
ylabel('Relative Translational Velocity (m/s)');
legend('Vx','Vy','Vz');
grid on

subplot(4,2,3); plot(tspan,Y(:,7),tspan,Y(:,8),tspan,Y(:,9),'linewidth',2);
xlabel('Time [sec]');
ylabel('Attitude Tracking Error (rad)');
legend('Th_x','Th_y','Th_z');
grid on

subplot(4,2,4);plot(tspan,Y(:,10),tspan,Y(:,11),tspan,Y(:,12),'linewidth',2);
xlabel('Time [sec]');
ylabel('Relative Position (m)');
legend('Bx','By','Bz');
grid on

subplot(4,2,5);plot(tspan,phi_c(:,1),tspan,phi_c(:,2),tspan,phi_c(:,3),'linewidth',2);
xlabel('Time [sec]');
ylabel('Control Torque (Nm)');
legend('Mx','My','Mz');
grid on

subplot(4,2,6);plot(tspan,phi_c(:,4),tspan,phi_c(:,5),tspan,phi_c(:,6),'linewidth',2);
xlabel('Time [sec]');
ylabel('Control Force (N)');
legend('Fx','Fy','Fz');
grid on

subplot(4,2,7); plot(tspan,s(:,1),tspan,s(:,2),tspan,s(:,3),'linewidth',2); 
xlabel('Time [sec]');
ylabel('Sliding surface (rotation)');
legend('Sx','Sy','Sz');
grid on

subplot(4,2,8); plot(tspan,s(:,4),tspan,s(:,5),tspan,s(:,6),'linewidth',2); 
xlabel('Time [sec]');
ylabel('Sliding surface (translation)');
legend('Sx','Sy','Sz');
grid on

figure(2);   
plot(Y(:,7),Y(:,1),'linewidth',2)
hold on
plot(Y(:,8),Y(:,2),'linewidth',2)
hold on
plot(Y(:,9),Y(:,3),'linewidth',2)
hold on
plot(Y(:,10),Y(:,4),'linewidth',2)
hold on
plot(Y(:,11),Y(:,5),'linewidth',2)
hold on
plot(Y(:,12),Y(:,6),'linewidth',2)
hold off
grid on