clearvars;
clear global;
close all;
clc;

% parmeters;
param;


%% Initial Conditions

X = [ 50 , 50 , 50]';                     % Initial relative position in meter

vr = [ 0.51 , 0.60 , -0.5]';         % Initial relative velocity

wr = (1)*[ 1.2 , 1.98 , -0.931]';     % Initial relative angular velocity

[s0,u0] = sv_from_coe([hs , e , RAAN , i_ang , omega , TA], mu);        % initial position and velocity of leader in eci frame
b0 = s0';                      
v0 = u0';

w0 = [0 0 0.0011]';                   % leader angular velocity

th0 = [60 ; 90 ; 60]*(pi/180);        % Principal rotation angle of leader
th =  [0 ; 60 ; 30]*(pi/180);         % Principal rotation angle of follower

R0 = so3_SO3(th0);                    % Leader rotation matrix
R = so3_SO3(th);                      % Follower rotation matrix

Rr = R0'*R;
Q = SO3_so3_back(Rr);

v = vr + Rr'*(v0 + skew(w0)*X);       % follower translation velocity

w = wr + Rr'*w0;                      % follower angular velocity

b = b0 + inv(R0')*X;                  % follower initial position           

w0_dot = [9 , 9 , 6]';

v0_dot = [1 1 1]';

theta = norm(Q);
if (theta==0)
    Von=eye(3);
else
    A = theta*cos(theta/2);
    C = sin(theta/2);
    Von = eye(3) - (1/2)*skew(Q) + (1/(theta^2))*(1-(A/(2*C)))*(skew(Q)*skew(Q));
end

B = Von*X;

g0 = [  R0  b0 ;...    % leader pose
      0 0 0 1 ];
  
g = [  R   b ;...      % follower pose
     0 0 0 1];
 
% Decomposition of g0 column wise
g01 = g0(:,1);
g02 = g0(:,2);
g03 = g0(:,3);
g04 = g0(:,4);

% Decomposition of g column wise
g1 = g(:,1);
g2 = g(:,2);
g3 = g(:,3);
g4 = g(:,4);

%% State Variables--------------------------------------
Er = [wr ; vr ];        
N  = [ Q ; B  ];         
E0 = [w0 ; v0 ];        
E  = [ w ; v  ];
E0_dot = [ w0_dot ; v0_dot ];

% State Vector-----------------------------------------
X0 = [Er ; N ; E0 ; E ; E0_dot ; g01 ; g02 ; g03 ; g04 ; g1 ; g2 ; g3 ; g4];                  

% tf inital guess
T = 1;

%% % Define time domain
t0 = 0;
tf = 10;
N = 100; % Number of collocation points
tTrain = linspace(t0, tf, N)';
dt = tTrain(2) - tTrain(1);

% Neural net setup
inputSize = 1;        % time
outputSize = 62;      % state vector
numHiddenUnits = 128;

% Define the network
layers = [
    featureInputLayer(inputSize,"Normalization","none","Name","input")
    fullyConnectedLayer(numHiddenUnits,"Name","fc1")
    tanhLayer("Name","tanh1")
    fullyConnectedLayer(numHiddenUnits,"Name","fc2")
    tanhLayer("Name","tanh2")
    fullyConnectedLayer(outputSize,"Name","output")];

dlnet = dlnetwork(layers);

% Training options
numEpochs = 2000;
learningRate = 1e-3;
gradDecay = 0.9;
sqGradDecay = 0.999;
trailingAvg = [];
trailingAvgSq = [];

% Training loop
lossHistory = zeros(numEpochs, 1);

for epoch = 1:numEpochs
    % Sample time points
    tBatch = tTrain(randperm(N, 64));  % mini-batch
    dlT = dlarray(tBatch','CB');       % format [1, batchSize]

    % Evaluate model gradients and loss
    [loss, gradients] = dlfeval(@modelGradients, dlnet, dlT);
    
    % Update network
    [dlnet, trailingAvg, trailingAvgSq] = adamupdate(dlnet, gradients, ...
        trailingAvg, trailingAvgSq, epoch, learningRate, gradDecay, sqGradDecay);
    
    % Store loss
    lossHistory(epoch) = double(gather(extractdata(loss)));
    
    if mod(epoch, 100) == 0
        disp("Epoch " + epoch + " | Loss = " + lossHistory(epoch));
    end
end

%% Predict over horizon for MPC
tHorizon = linspace(0, 2, 50)'; % 2-sec prediction
dlT_horizon = dlarray(tHorizon','CB');
dlX_pred = predict(dlnet, dlT_horizon);
X_pred = extractdata(dlX_pred)';  % Each row: state at a time point

% Now use X_pred in your MPC controller



%% intial trajectory generation
%pehle wala code ->modify post implementing mpc or just let it be idk
% disp('generating initial trajectory');
% %assuming linear variation of sigma
% dt = 0.01;
% %code run time:
% t = 10;
% 
% tspan = (0:dt:t)';
% 
% 
% %Y = initial state trajectory
% [Y,Tc,sms] = rk4(@coupledEOM_smc_only,tspan,X0);
% 
% disp('done');
% %% %% plot trajectory (3d)
% figure(1);
% hold on; grid on;
% title('Spacecraft Path with 3D Obstacles');
% xlabel('x');
% ylabel('y');
% zlabel('z');
% for i = 1:Obstacle_count
%     % Extract 3D obstacle data for plotting
%     x_obst = squeeze(obstacle(i, :, 1));
%     y_obst = squeeze(obstacle(i, :, 2));
%     z_obst = squeeze(obstacle(i, :, 3));
% 
%     x_obst=reshape(x_obst, size(phi));
%     y_obst=reshape(y_obst, size(phi));
%     z_obst=reshape(z_obst, size(phi));
% 
%     size(x_obst)
%     size(y_obst)
%     % Plot the obstacle
%     scatter3((obstacle(i,:,1)), (obstacle(i,:,2)), ...
%             (obstacle(i,:,3)), 'r.');
% 
% 
% end
% n(:,1:6) = Y(:,7:12);
% B=Y(:,10:12);
% x=B(:,1);
% y=B(:,2);
% z=B(:,3);
% size(x)
% size(y)
% size(z)
% 
% 
% plot3(x,y,z,'b');
% plot3(x(1),y(1),z(1),'ko'); % start point
% plot3(0,0,0,'rs'); % goal point
% view(3);
% hold off;
% 
% 
% %% 
% 
% 
% 
% %% ============================ plot ====================================
% figure(2);
% plot(tspan,Y(:,1),tspan,Y(:,2),tspan,Y(:,3),'linewidth',2); 
% xlabel('Time [sec]');
% ylabel('Relative Angular Velocity (rad/s)');
% legend('Wx','Wy','Wz');
% grid on
% 
% 
% figure(3);plot(tspan,Y(:,4),tspan,Y(:,5),tspan,Y(:,6),'linewidth',2);
% xlabel('Time [sec]');
% ylabel('Relative Translational Velocity (m/s)');
% legend('Vx','Vy','Vz');
% grid on
% 
% figure(4); plot(tspan,Y(:,7),tspan,Y(:,8),tspan,Y(:,9),'linewidth',2);
% xlabel('Time [sec]');
% ylabel('Attitude Tracking Error (rad)');
% legend('Th_x','Th_y','Th_z');
% grid on
% 
% figure(5); plot(tspan,Y(:,10),tspan,Y(:,11),tspan,Y(:,12),'linewidth',2);
% xlabel('Time [sec]');
% ylabel('Relative Position (m)');
% legend('Bx','By','Bz');
% grid on
% 
% figure(6); plot(tspan,Tc(:,1),tspan,Tc(:,2),tspan,Tc(:,3),'linewidth',2);
% xlabel('Time [sec]');
% ylabel('Control Torque (Nm)');
% legend('Mx','My','Mz');
% grid on
% 
% figure(7); plot(tspan,Tc(:,4),tspan,Tc(:,5),tspan,Tc(:,6),'linewidth',2);
% xlabel('Time [sec]');
% ylabel('Control Force (N)');
% legend('Fx','Fy','Fz');
% grid on
% 
% figure(8); plot(tspan,sms(:,1),tspan,sms(:,2),tspan,sms(:,3),'linewidth',2); 
% xlabel('Time [sec]');
% ylabel('Sliding surface (rotation)');
% legend('Sx','Sy','Sz');
% grid on
% 
% figure(9); plot(tspan,sms(:,4),tspan,sms(:,5),tspan,sms(:,6),'linewidth',2); 
% xlabel('Time [sec]');
% ylabel('Sliding surface (translation)');
% legend('Sx','Sy','Sz');
% grid on
% 
