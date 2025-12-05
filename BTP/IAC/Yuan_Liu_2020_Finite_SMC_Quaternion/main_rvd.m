clearvars;
clear global;
close all;
clc;

% parmeters;
param;

import casadi.*

%% Initial Conditions
qv_t0 = [0;0;0];
omega_t0 = [0.01;0.01;0.01];
r_t0 = [1; 1; 1] * 7.078e8;
v_t0 = [2; 3; -2] * 1e4;
Xt0 = [qv_t0; omega_t0; r_t0; v_t0];  


qv0 = [sin(deg2rad(19.9984)/2);
       sin(deg2rad(-9.9987)/2);
       sin(deg2rad(15.0050)/2)];
omega0 = [0;0;0];
r0 = [10;10;-10];
v0 = [0;0;0];               

% Initial State Vector-----------------------------------------
X0 = [qv0; omega0; r0; v0];                 

% tf inital guess
T = 1;

%% intial trajectory generation
disp('generating initial trajectory');
%assuming linear variation of sigma
dt = 0.1;
%code run time:
t = 150;
tspan = (0:dt:t)';

[target_state,q_t,omega_t] = rk4t(@EOM_Target,tspan,Xt0);

[X,u,s] = rk4(@EOM_RvD,tspan,X0, target_state);
disp('done');


%% Plots

% Relative position
figure(1);
subplot(2,2,1);
plot(tspan, X(:,7), 'r', 'LineWidth',1.5); hold on;
plot(tspan, X(:,8), 'g', 'LineWidth',1.5);
plot(tspan, X(:,9), 'b', 'LineWidth',1.5); grid on;
title('Relative Position'); xlabel('Time [s]'); ylabel('Position Error [m]');
legend('x','y','z');


% Relative velocity
subplot(2,2,2);
plot(tspan, X(:,10), 'r', 'LineWidth',1.5); hold on;
plot(tspan, X(:,11), 'g', 'LineWidth',1.5);
plot(tspan, X(:,12), 'b', 'LineWidth',1.5); grid on;
title('Relative Velocity'); xlabel('Time [s]'); ylabel('Velocity Error [m/s]');
legend('x','y','z');


% Quaternion vector part
subplot(2,2,3);
plot(tspan, X(:,1), 'r', 'LineWidth',1.5); hold on;
plot(tspan, X(:,2), 'g', 'LineWidth',1.5);
plot(tspan, X(:,3), 'b', 'LineWidth',1.5); grid on;
title('Quaternion Vector Part Error'); xlabel('Time [s]'); ylabel('q_v');
legend('q_1','q_2','q_3');


% Relative angular velocity
subplot(2,2,4);
plot(tspan, X(:,4), 'r', 'LineWidth',1.5); hold on;
plot(tspan, X(:,5), 'g', 'LineWidth',1.5);
plot(tspan, X(:,6), 'b', 'LineWidth',1.5); grid on;
title('Relative Angular Velocity'); xlabel('Time [s]'); ylabel('Angular Velocity [rad/s]');
legend('x','y','z');


% Control torque tau_r
figure(2);
subplot(2,2,1);
plot(tspan, u(:,1), 'r', 'LineWidth',1.5); hold on;
plot(tspan, u(:,2), 'g', 'LineWidth',1.5);
plot(tspan, u(:,3), 'b', 'LineWidth',1.5); grid on;
title('Control Torque \tau_r'); xlabel('Time [s]'); ylabel('\tau_r [Nm]');
legend('\tau_{r1}','\tau_{r2}','\tau_{r3}');


% Control force tau_t
subplot(2,2,2);
plot(tspan, u(:,4), 'r', 'LineWidth',1.5); hold on;
plot(tspan, u(:,5), 'g', 'LineWidth',1.5);
plot(tspan, u(:,6), 'b', 'LineWidth',1.5); grid on;
title('Control Force \tau_t'); xlabel('Time [s]'); ylabel('\tau_t [N]');
legend('\tau_{t1}','\tau_{t2}','\tau_{t3}');

% Sliding variable for rotation
subplot(2,2,3);
plot(tspan, s(:,1), 'r', 'LineWidth',1.5); hold on;
plot(tspan, s(:,2), 'g', 'LineWidth',1.5);
plot(tspan, s(:,3), 'b', 'LineWidth',1.5); grid on;
title('Sliding variable for Rotation'); xlabel('Time [s]'); ylabel('Sliding Variable (S)');
legend('x','y','z');

% Sliding variable for translational
subplot(2,2,4);
plot(tspan, s(:,4), 'r', 'LineWidth',1.5); hold on;
plot(tspan, s(:,5), 'g', 'LineWidth',1.5);
plot(tspan, s(:,6), 'b', 'LineWidth',1.5); grid on;
title('Sliding variable for Translation'); xlabel('Time [s]'); ylabel('Sliding Variable (S)');
legend('x','y','z');


% figure(3);
% x_axis = linspace(-0.1,0.1);
% y_axis = linspace(-0.1,0.1);
% plot(-k1.*x_axis,y_axis, 'k', 'LineWidth', 1.5); hold on;
% plot(X(:,1),X(:,4), 'r', 'LineWidth',1.5); hold on;
% plot(X(:,2),X(:,5), 'g', 'LineWidth',1.5); hold on;
% plot(X(:,3),X(:,6), 'b', 'LineWidth',1.5); hold on;
% title('Phase portrait of rotational motion');
% grid on
% 
% figure(4);
% x_axis = linspace(-1,1);
% y_axis = linspace(-1,1);
% plot(-k1.*x_axis,y_axis, 'k', 'LineWidth', 1.5); hold on;
% plot(X(:,7),X(:,10), 'r', 'LineWidth',1.5); hold on;
% plot(X(:,8),X(:,11), 'g', 'LineWidth',1.5); hold on;
% plot(X(:,9),X(:,12), 'b', 'LineWidth',1.5); hold on;
% title('Phase portrait of translational motion');
% grid on