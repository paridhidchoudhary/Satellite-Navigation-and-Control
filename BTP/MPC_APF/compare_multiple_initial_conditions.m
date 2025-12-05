function compare_multiple_initial_conditions()
    % Define different initial conditions to test
    % Each row represents a different initial state vector
    % Based on the error, I'm assuming the state vector has indices:
    % [roll, pitch, yaw, x, y, z, roll_rate, pitch_rate, yaw_rate, x_vel, y_vel, z_vel]
    initial_states = [
        % First initial state (original)
        0, 0, 0,          % roll, pitch, yaw
        0, 0, 0,          % x, y, z positions
        0, 0, 0,          % roll, pitch, yaw rates
        0, 0, 0;          % x, y, z velocities
        
        % Second initial state (different position)
        0, 0, 0,          % roll, pitch, yaw
        2, 1, 1,          % x, y, z positions
        0, 0, 0,          % roll, pitch, yaw rates
        0, 0, 0;          % x, y, z velocities
         
        % Third initial state (with initial velocities)
        0, 0, 0,          % roll, pitch, yaw
        -1, -1, 0,        % x, y, z positions
        0, 0, 0,          % roll, pitch, yaw rates
        0.5, 0.5, 0       % x, y, z velocities
    ];
    
    % Labels for the different initial conditions (for legends)
    condition_labels = {'Original', 'Offset Position', 'With Initial Velocity'};
    
    % Colors for plotting different conditions
    plot_colors = {'r', 'b', 'g'};
    
    % Run simulations for each initial condition
    num_conditions = size(initial_states, 1);
    all_results = cell(1, num_conditions);
    
    for i = 1:num_conditions
        fprintf('Simulating condition %d: %s\n', i, condition_labels{i});
        % Ensure the initial state is in the correct format
        current_initial_state = initial_states(i,:);
        
        [actual_trajectory, actual_trajectory_del, time_log, vapf_del] = ...
            simulate_quadrotor_with_mpc(current_initial_state);
        
        % Store results
        all_results{i}.actual = actual_trajectory;
        all_results{i}.actual_del = actual_trajectory_del;
        all_results{i}.time = time_log;
        all_results{i}.vapf_del = vapf_del;
        all_results{i}.label = condition_labels{i};
        all_results{i}.color = plot_colors{i};
    end
    
    % Plot comparative results
    plot_comparative_results(all_results);
end

function [actual_trajectory, actual_trajectory_del, time_log, vapf_del] = simulate_quadrotor_with_mpc(initial_state_override)
    % Load parameters
    param;
    
    % Override initial state if provided
    if nargin > 0
        % Make sure initial_state is properly sized/formatted
        % This assumes initial_state in param is a column vector with 12 elements
        if isrow(initial_state_override)
            initial_state = initial_state_override'; % Convert row to column if needed
        else
            initial_state = initial_state_override;
        end
    end

    % Generate reference trajectory using APF
    [xr, vapf_del, time_log] = generate_apf_trajectory(initial_state);

    % Initialize quadrotor state
    state = initial_state;
    state_del = initial_state;

    % Simulation parameters
    N_sim = length(time_log);
    N_horizon = 12;  % MPC prediction horizon

    % Preallocate arrays for results
    actual_trajectory = zeros(12, N_sim);
    control_inputs = zeros(4, N_sim);
    actual_trajectory_del = zeros(12, N_sim);
    control_inputs_del = zeros(4, N_sim);

    % Main simulation loop
    for i = 1:N_sim
        % Current reference
        ref_traj = xr(1:12, i:min(i+N_horizon, N_sim));
        ref_traj_del = xr(13:24, i:min(i+N_horizon, N_sim));

        % If reference is shorter than horizon, pad with last value
        if size(ref_traj, 2) < N_horizon+1
            ref_traj = [ref_traj, repmat(ref_traj(:,end), 1, N_horizon+1-size(ref_traj,2))];
            ref_traj_del = [ref_traj_del, repmat(ref_traj_del(:,end), 1, N_horizon+1-size(ref_traj_del,2))];
        end

        % Get MPC control inputs
        [u_opt, x_pred] = quad_mpc(state, ref_traj, N_horizon);
        [u_opt_del, x_pred] = quad_mpc(state_del, ref_traj_del, N_horizon);

        % Apply first control input
        u = u_opt(:,1);
        u_del = u_opt_del(:,1);

        % Simulate system dynamics
        state_dot = quad_dynamics(state, u(2), u(3), u(4), u(1));
        state_dot_del = quad_dynamics(state_del, u_del(2), u_del(3), u_del(4), u_del(1));

        % Euler integration
        state = state + dT * state_dot;
        state_del = state_del + dT * state_dot_del;

        % Store results
        actual_trajectory(:,i) = state;
        control_inputs(:,i) = u;
        actual_trajectory_del(:,i) = state_del;
        control_inputs_del(:,i) = u_del;

        % Display progress
        if mod(i, 10) == 0
            fprintf('Simulation progress: %.1f%%\n', 100*i/N_sim);
        end
    end
end

% Function to generate APF trajectory (modified to take initial state)
function [xr, vapf_del, time_log] = generate_apf_trajectory(initial_state_override)
    param;
    
    % Override initial state if provided
    if nargin > 0
        state = initial_state_override;
    else
        state = initial_state;
    end

    time_log = zeros(1, N);
    
    xr = zeros(24, N); % 1st 12, normal apf, next 12, D-APF
    
    % Ensure state is a column vector and properly formatted
    if isrow(state)
        state = state';
    end
    
    % Make sure state is at least 6 elements long
    if length(state) < 6
        error('State vector must have at least 6 elements');
    end
    
    % Access the position components of the state
    % Here we assume indices 4, 5, 6 correspond to x, y, z positions
    pos = [state(4); state(5); state(6); state(4); state(5); state(6)];
    vapf_del = zeros(3, N);
    
    for i=1:N
        artificial_pot = @(s) apf(s);
        l1 = artificial_pot(pos);
        l2 = artificial_pot(pos + 0.5*dT*l1);
        l3 = artificial_pot(pos + 0.5*dT*l2);
        l4 = artificial_pot(pos + dT*l3);
        pos = pos + (dT/6)*(l1 + 2*l2 + 2*l3 + l4);
        
        xr(:,i) = zeros(24,1);
        xr(4,i) = pos(1,1);
        xr(5,i) = pos(2,1);
        xr(6,i) = pos(3,1);
        xr(16,i) = pos(4,1);
        xr(17,i) = pos(5,1);
        xr(18,i) = pos(6,1);
        
        time_log(i) = i*dT;
        vapf_del(:,i) = l1(4:6,1);
    end
end

% Function to plot comparative results from multiple initial conditions
function plot_comparative_results(all_results)
    param;
    num_conditions = length(all_results);
    
    % Position error plots (APF)
    figure('Name', 'Position Error - APF', 'NumberTitle', 'off');
    
    % X position error
    subplot(3,1,1);
    hold on;
    for i = 1:num_conditions
        plot(all_results{i}.time, all_results{i}.actual(4,:) - goal(1), all_results{i}.color);
    end
    title('X Position Error (APF)');
    legend_entries = cellfun(@(x) x.label, all_results, 'UniformOutput', false);
    legend(legend_entries);
    xlabel('Time (s)');
    ylabel('Error (m)');
    grid on;
    
    % Y position error
    subplot(3,1,2);
    hold on;
    for i = 1:num_conditions
        plot(all_results{i}.time, all_results{i}.actual(5,:) - goal(2), all_results{i}.color);
    end
    title('Y Position Error (APF)');
    xlabel('Time (s)');
    ylabel('Error (m)');
    grid on;
    
    % Z position error
    subplot(3,1,3);
    hold on;
    for i = 1:num_conditions
        plot(all_results{i}.time, all_results{i}.actual(6,:) - goal(3), all_results{i}.color);
    end
    title('Z Position Error (APF)');
    xlabel('Time (s)');
    ylabel('Error (m)');
    grid on;
    
    % Position error plots (D-APF)
    figure('Name', 'Position Error - D-APF', 'NumberTitle', 'off');
    
    % X position error
    subplot(3,1,1);
    hold on;
    for i = 1:num_conditions
        plot(all_results{i}.time, all_results{i}.actual_del(4,:) - goal(1), all_results{i}.color);
    end
    title('X Position Error (D-APF)');
    legend(legend_entries);
    xlabel('Time (s)');
    ylabel('Error (m)');
    grid on;
    
    % Y position error
    subplot(3,1,2);
    hold on;
    for i = 1:num_conditions
        plot(all_results{i}.time, all_results{i}.actual_del(5,:) - goal(2), all_results{i}.color);
    end
    title('Y Position Error (D-APF)');
    xlabel('Time (s)');
    ylabel('Error (m)');
    grid on;
    
    % Z position error
    subplot(3,1,3);
    hold on;
    for i = 1:num_conditions
        plot(all_results{i}.time, all_results{i}.actual_del(6,:) - goal(3), all_results{i}.color);
    end
    title('Z Position Error (D-APF)');
    xlabel('Time (s)');
    ylabel('Error (m)');
    grid on;
    
    % 3D Trajectory plot (APF)
    figure('Name', '3D Trajectories - APF', 'NumberTitle', 'off');
    hold on;
    
    % Plot all trajectories
    for i = 1:num_conditions
        plot3(all_results{i}.actual(4,:), all_results{i}.actual(5,:), all_results{i}.actual(6,:), ...
              all_results{i}.color, 'LineWidth', 2);
    end
    
    % Mark the goal
    scatter3(goal(1), goal(2), goal(3), 100, 'k', 'filled');
    
    % Plot initial positions for each condition
    for i = 1:num_conditions
        scatter3(all_results{i}.actual(4,1), all_results{i}.actual(5,1), all_results{i}.actual(6,1), ...
                80, all_results{i}.color, 'filled', 'MarkerEdgeColor', 'k');
    end
    
    % Plot obstacles
    for i = 1:Obstacle_count
        obs_x = squeeze(obstacle(i,:,1));
        obs_y = squeeze(obstacle(i,:,2));
        obs_z = squeeze(obstacle(i,:,3));
        fill3(obs_x, obs_y, obs_z, [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'k');
    end
    
    % Aesthetics
    title('APF 3D Quadrotor Trajectories');
    legend([legend_entries, {'Goal'}]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on;
    axis equal;
    view(3);
    
    % 3D Trajectory plot (D-APF)
    figure('Name', '3D Trajectories - D-APF', 'NumberTitle', 'off');
    hold on;
    
    % Plot all trajectories
    for i = 1:num_conditions
        plot3(all_results{i}.actual_del(4,:), all_results{i}.actual_del(5,:), all_results{i}.actual_del(6,:), ...
              all_results{i}.color, 'LineWidth', 2);
    end
    
    % Mark the goal
    scatter3(goal(1), goal(2), goal(3), 100, 'k', 'filled');
    
    % Plot initial positions for each condition
    for i = 1:num_conditions
        scatter3(all_results{i}.actual_del(4,1), all_results{i}.actual_del(5,1), all_results{i}.actual_del(6,1), ...
                80, all_results{i}.color, 'filled', 'MarkerEdgeColor', 'k');
    end
    
    % Plot obstacles
    for i = 1:Obstacle_count
        obs_x = squeeze(obstacle(i,:,1));
        obs_y = squeeze(obstacle(i,:,2));
        obs_z = squeeze(obstacle(i,:,3));
        fill3(obs_x, obs_y, obs_z, [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'k');
    end
    
    % Aesthetics
    title('D-APF 3D Quadrotor Trajectories');
    legend([legend_entries, {'Goal'}]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on;
    axis equal;
    view(3);
    
    % Attitude comparison
    figure('Name', 'Attitude Comparison', 'NumberTitle', 'off');
    subplot(2,1,1);
    hold on;
    for i = 1:num_conditions
        plot(all_results{i}.time, all_results{i}.actual(1,:), all_results{i}.color);
    end
    title('Roll Angle Comparison (APF)');
    legend(legend_entries);
    xlabel('Time (s)');
    ylabel('Roll (rad)');
    grid on;
    
    subplot(2,1,2);
    hold on;
    for i = 1:num_conditions
        plot(all_results{i}.time, all_results{i}.actual_del(1,:), all_results{i}.color);
    end
    title('Roll Angle Comparison (D-APF)');
    xlabel('Time (s)');
    ylabel('Roll (rad)');
    grid on;
    
    % Translational velocity comparison
    figure('Name', 'Velocity Comparison', 'NumberTitle', 'off');
    
    % X velocity
    subplot(3,2,1);
    hold on;
    for i = 1:num_conditions
        plot(all_results{i}.time, all_results{i}.actual(10,:), all_results{i}.color);
    end
    title('X Velocity (APF)');
    legend(legend_entries);
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    grid on;
    
    subplot(3,2,2);
    hold on;
    for i = 1:num_conditions
        plot(all_results{i}.time, all_results{i}.actual_del(10,:), all_results{i}.color);
    end
    title('X Velocity (D-APF)');
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    grid on;
    
    % Y velocity
    subplot(3,2,3);
    hold on;
    for i = 1:num_conditions
        plot(all_results{i}.time, all_results{i}.actual(11,:), all_results{i}.color);
    end
    title('Y Velocity (APF)');
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    grid on;
    
    subplot(3,2,4);
    hold on;
    for i = 1:num_conditions
        plot(all_results{i}.time, all_results{i}.actual_del(11,:), all_results{i}.color);
    end
    title('Y Velocity (D-APF)');
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    grid on;
    
    % Z velocity
    subplot(3,2,5);
    hold on;
    for i = 1:num_conditions
        plot(all_results{i}.time, all_results{i}.actual(12,:), all_results{i}.color);
    end
    title('Z Velocity (APF)');
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    grid on;
    
    subplot(3,2,6);
    hold on;
    for i = 1:num_conditions
        plot(all_results{i}.time, all_results{i}.actual_del(12,:), all_results{i}.color);
    end
    title('Z Velocity (D-APF)');
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    grid on;
    
    % Distance to goal over time
    figure('Name', 'Distance to Goal', 'NumberTitle', 'off');
    
    % APF
    subplot(2,1,1);
    hold on;
    for i = 1:num_conditions
        distances = sqrt(sum((all_results{i}.actual(4:6,:) - repmat(goal, 1, size(all_results{i}.time, 2))).^2, 1));
        plot(all_results{i}.time, distances, all_results{i}.color);
    end
    title('Distance to Goal (APF)');
    legend(legend_entries);
    xlabel('Time (s)');
    ylabel('Distance (m)');
    grid on;
    
    % D-APF
    subplot(2,1,2);
    hold on;
    for i = 1:num_conditions
        distances = sqrt(sum((all_results{i}.actual_del(4:6,:) - repmat(goal, 1, size(all_results{i}.time, 2))).^2, 1));
        plot(all_results{i}.time, distances, all_results{i}.color);
    end
    title('Distance to Goal (D-APF)');
    xlabel('Time (s)');
    ylabel('Distance (m)');
    grid on;
end