clearvars; clear global; close all; clc;

% Load parameters
param;

%% Initial Conditions
X = [50, 50, 50]';                     % Initial relative position in meter
vr = [0.51, 0.60, -0.5]';              % Initial relative velocity
wr = [1.2, 1.98, -0.931]';             % Initial relative angular velocity

% Initial position and velocity of leader in ECI frame
[s0, u0] = sv_from_coe([hs, e, RAAN, i_ang, omega, TA], mu); 
b0 = s0';
v0 = u0';
w0 = [0, 0, 0.0011]';                  % Leader angular velocity

% Principal rotation angles (in radians)
th0 = [60; 90; 60] * (pi/180);         % Leader
th = [0; 60; 30] * (pi/180);           % Follower

% Rotation matrices
R0 = so3_SO3(th0);                     % Leader rotation matrix
R = so3_SO3(th);                       % Follower rotation matrix
Rr = R0' * R;
Q = SO3_so3_back(Rr);

% Follower states
v = vr + Rr' * (v0 + skew(w0) * X);    % Follower translation velocity
w = wr + Rr' * w0;                     % Follower angular velocity
b = b0 + R0' \ X;                      % Follower initial position

% Leader acceleration
w0_dot = [9, 9, 6]';
v0_dot = [1, 1, 1]';

% Von calculation for B
theta = norm(Q);
if (theta == 0)
    Von = eye(3);
else
    A = theta * cos(theta/2);
    C = sin(theta/2);
    Von = eye(3) - (1/2) * skew(Q) + (1/(theta^2)) * (1 - (A/(2*C))) * (skew(Q) * skew(Q));
end
B = Von * X;

% Pose matrices
g0 = [R0, b0; 0, 0, 0, 1];             % Leader pose
g = [R, b; 0, 0, 0, 1];                % Follower pose

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

%% State Variables
Er = [wr; vr];
N = [Q; B];
E0 = [w0; v0];
E = [w; v];
E0_dot = [w0_dot; v0_dot];

% Initial state vector
X0 = [Er; N; E0; E; E0_dot; g01; g02; g03; g04; g1; g2; g3; g4];

%% Define time domain for training
t0 = 0;
tf = 10;
num_train_points = 100;                % Number of collocation points
tTrain = linspace(t0, tf, num_train_points)';
dt = tTrain(2) - tTrain(1);

%% Neural Network for Prediction Horizon
% First network: State prediction
inputSize_state = 1;                   % Time input
outputSize_state = 62;                 % State vector
numHiddenUnits_state = 128;

% Define the state prediction network
layers_state = [
    featureInputLayer(inputSize_state, "Normalization", "none", "Name", "input")
    fullyConnectedLayer(numHiddenUnits_state, "Name", "fc1")
    tanhLayer("Name", "tanh1")
    fullyConnectedLayer(numHiddenUnits_state, "Name", "fc2")
    tanhLayer("Name", "tanh2")
    fullyConnectedLayer(outputSize_state, "Name", "output")
];
dlnet_state = dlnetwork(layers_state);

% Second network: MPC Horizon prediction
inputSize_horizon = 62;                % Current state vector
outputSize_horizon = 1;                % Prediction horizon N
numHiddenUnits_horizon = 64;

% Define the horizon prediction network
layers_horizon = [
    featureInputLayer(inputSize_horizon, "Normalization", "none", "Name", "input")
    fullyConnectedLayer(numHiddenUnits_horizon, "Name", "fc1")
    tanhLayer("Name", "tanh1")
    fullyConnectedLayer(numHiddenUnits_horizon, "Name", "fc2")
    tanhLayer("Name", "tanh2")
    fullyConnectedLayer(outputSize_horizon, "Name", "output")
    reluLayer("Name", "relu")          % Ensure positive horizon values
];
dlnet_horizon = dlnetwork(layers_horizon);

%% Training options
numEpochs = 200;
miniBatchSize = 64;
learningRate = 1e-3;
gradDecay = 0.9;
sqGradDecay = 0.999;

% Initialize optimizer state
trailingAvg_state = [];
trailingAvgSq_state = [];
trailingAvg_horizon = [];
trailingAvgSq_horizon = [];

% Training loop
lossHistory_state = zeros(numEpochs, 1);
lossHistory_horizon = zeros(numEpochs, 1);

% Get initial reference trajectory using ODE solver
[t_ref, X_ref] = ode45(@(t, X) coupledEOM_smc_only(t, X), tTrain, X0);
X_ref = X_ref';  % Convert to [state_dim, num_points]

%% Training Loop for State Prediction Network
fprintf('Training state prediction network...\n');
for epoch = 1:numEpochs
    % Sample mini-batch
    indices = randperm(num_train_points, miniBatchSize);
    t_batch = tTrain(indices);
    X_batch = X_ref(:, indices);
    
    % Format as dlarray
    dlT = dlarray(t_batch', 'CB');  % [1, batchSize]
    dlX_true = dlarray(X_batch, 'CB');  % [stateDim, batchSize]
    
    % Evaluate model gradients and loss
    [loss_state, gradients_state] = dlfeval(@modelGradients_state, dlnet_state, dlT, dlX_true);
    
    % Update network
    [dlnet_state, trailingAvg_state, trailingAvgSq_state] = adamupdate(dlnet_state, gradients_state, ...
        trailingAvg_state, trailingAvgSq_state, epoch, learningRate, gradDecay, sqGradDecay);
    
    % Store loss
    lossHistory_state(epoch) = double(gather(extractdata(loss_state)));
    
    if mod(epoch, 100) == 0
        fprintf('State Network - Epoch %d/%d | Loss = %.4e\n', epoch, numEpochs, lossHistory_state(epoch));
    end
end
%output of the above code -> dlnet_state (mainly) -> it is the learned
%model that maps t->X
%% Training Loop for Horizon Prediction Network
fprintf('\nTraining horizon prediction network...\n');

% Generate optimal horizon data (This would be your actual MPC optimization)
% For demonstration, we'll use a simple heuristic based on state error magnitude
N_optimal = zeros(num_train_points, 1);
for i = 1:num_train_points
    state_error = norm(X_ref(1:12, i));  % Error in relative states
    
    % Simple heuristic: larger error -> longer horizon (you would replace this)
    N_optimal(i) = 5 + 15 * min(state_error / 50, 1);  % N between 5 and 20
end

for epoch = 1:numEpochs
    % Sample mini-batch
    indices = randperm(num_train_points, miniBatchSize);
    X_batch = X_ref(:, indices);
    N_batch = N_optimal(indices);
    
    % Format as dlarray
    dlX = dlarray(X_batch, 'CB');  % [stateDim, batchSize]
    dlN_true = dlarray(N_batch', 'CB');  % [1, batchSize]
    
    % Evaluate model gradients and loss
    [loss_horizon, gradients_horizon] = dlfeval(@modelGradients_horizon, dlnet_horizon, dlX, dlN_true);
    
    % Update network
    [dlnet_horizon, trailingAvg_horizon, trailingAvgSq_horizon] = adamupdate(dlnet_horizon, gradients_horizon, ...
        trailingAvg_horizon, trailingAvgSq_horizon, epoch, learningRate, gradDecay, sqGradDecay);
    
    % Store loss
    lossHistory_horizon(epoch) = double(gather(extractdata(loss_horizon)));
    
    if mod(epoch, 100) == 0
        fprintf('Horizon Network - Epoch %d/%d | Loss = %.4e\n', epoch, numEpochs, lossHistory_horizon(epoch));
    end
end

%% Plot loss history
figure;
subplot(2,1,1);
semilogy(lossHistory_state);
title('State Prediction Network Loss');
xlabel('Epoch');
ylabel('Loss');
grid on;

subplot(2,1,2);
semilogy(lossHistory_horizon);
title('Horizon Prediction Network Loss');
xlabel('Epoch');
ylabel('Loss');
grid on;

%% MPC Implementation using Neural Network Horizon

% Current state (for demonstration)
current_state = X0;
current_time = 0;

% Simulation parameters
sim_time = 10;
dt_sim = 0.1;
num_steps = sim_time / dt_sim;

% Storage for results
state_history = zeros(62, num_steps+1);
horizon_history = zeros(1, num_steps+1);
control_history = zeros(6, num_steps);
time_history = (0:num_steps) * dt_sim;

state_history(:,1) = current_state;

% MPC loop
for step = 1:num_steps
    % Get current state
    current_state = state_history(:, step);
    current_time = time_history(step);
    
    % Predict optimal horizon using NN
    dlX_current = dlarray(current_state, 'CB');
    dlN_pred = predict(dlnet_horizon, dlX_current);
    N_mpc = max(5, round(extractdata(dlN_pred))); % Ensure minimum horizon of 5
    
    horizon_history(step) = N_mpc;
    
    % MPC optimization with dynamic horizon N_mpc
    [Tc, ~] = mpcController(current_state, current_time, N_mpc);
    
    % Apply control and simulate one step
    [t_step, X_step] = ode45(@(t, X) coupledEOM_smc_only(t, X), [current_time, current_time+dt_sim], current_state);
    
    % Store results
    control_history(:, step) = Tc;
    state_history(:, step+1) = X_step(end, :)';
end
fprintf('Final Prediction Horizon for mpc= %d \n',horizon_history(1,num_steps+1));
%% Plot Results
figure;
subplot(2,1,1);
plot(time_history, horizon_history);
title('Dynamic MPC Horizon');
xlabel('Time (s)');
ylabel('Prediction Horizon N');
grid on;

subplot(2,1,2);
plot(time_history(1:end-1), vecnorm(control_history(1:3,:)), 'b-', 'LineWidth', 1.5);
hold on;
plot(time_history(1:end-1), vecnorm(control_history(4:6,:)), 'r-', 'LineWidth', 1.5);
title('Control Effort');
xlabel('Time (s)');
ylabel('Magnitude');
legend('Torque', 'Force');
grid on;

%% Define gradient functions for PINN training

function [loss, gradients] = modelGradients_state(dlnet, dlT, dlX_true)
    % Forward pass to get predicted states

    dlX_pred = forward(dlnet, dlT);
    
    % Calculate MSE loss between predicted and true states
    prediction_loss = mean(sum((dlX_pred - dlX_true).^2, 1));
    
    % Compute dX/dt via automatic differentiation w.r.t. time
    dlDX = dlgradient(sum(dlX_pred, 'all'), dlT);
    
    % Evaluate true dynamics for each time point
    numSamples = size(dlT, 2);
    physics_residuals = zeros(size(dlX_pred), 'like', dlX_pred);
    
    for i = 1:numSamples
        t_val = extractdata(dlT(1,i));
        x_val = extractdata(dlX_pred(:,i));
        
        % Evaluate dynamics
        [x_dot_true, ~, ~] = coupledEOM_smc_only(t_val, x_val);
        physics_residuals(:,i) = dlDX(:,i) - x_dot_true;
    end
    
    % Compute physics loss (MSE over residuals)
    physics_loss = mean(sum(physics_residuals.^2, 1));
    
    % Total loss with weighting
    loss = prediction_loss + 0.1 * physics_loss;
    
    % Get gradients w.r.t. the learnable parameters
    gradients = dlgradient(loss, dlnet.Learnables);
end

function [loss, gradients] = modelGradients_horizon(dlnet, dlX, dlN_true)
    % Forward pass to get predicted horizon
    dlN_pred = forward(dlnet, dlX);
    
    % Mean squared error loss
    mse_loss = mean((dlN_pred - dlN_true).^2);
    
    % Smoothness regularization - penalize abrupt changes in horizon predictions
    smoothness_loss = 0;
    if size(dlX, 2) > 1
        dlX_sorted = dlX(:, [2:end, 1]);  % Shifted states (simulation of sequential states)
        dlN_next = forward(dlnet, dlX_sorted);
        smoothness_loss = mean((dlN_next - dlN_pred).^2);
    end
    
    % Total loss
    loss = mse_loss + 0.1 * smoothness_loss;
    
    % Get gradients w.r.t. the learnable parameters
    gradients = dlgradient(loss, dlnet.Learnables);
end

function [Tc, X_mpc] = mpcController(current_state, current_time, N)
    % This function would implement your MPC controller
    % Using the predicted horizon N
    % ...
    
    % For demonstration, we'll return a simple control based on current state
    Er = current_state(1:6);
    N_state = current_state(7:12);
    
    % Simple PD controller for demonstration
    Kp = diag([0.1, 0.1, 0.1, 0.05, 0.05, 0.05]);
    Kd = diag([0.5, 0.5, 0.5, 0.2, 0.2, 0.2]);
    
    Tc = -Kp * N_state - Kd * Er;
    
    % In a real implementation, you would solve the MPC optimization problem
    % using the dynamic horizon N
    
    % Placeholder for MPC predicted trajectory
    X_mpc = zeros(62, N);
end