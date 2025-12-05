function simulate_quadrotor_with_mpc()
    % Load parameters
    param;

    
    % Initialize quadrotor state
    state = initial_state1;
    state_del1 = initial_state1;
    state_del2 = initial_state2;
    state_del3 = initial_state3;

    % Generate reference trajectory using APF
    [xr, vapf, time_log] = generate_apf_trajectory(state);
    [xr1, vapf_del1, time_log] = generate_apf_trajectory(state_del1);
    [xr2, vapf_del2, time_log] = generate_apf_trajectory(state_del2);
    [xr3, vapf_del3, time_log] = generate_apf_trajectory(state_del3);

    % Simulation parameters
    N_sim = length(time_log);
    N_horizon = 12;  % MPC prediction horizon

    % Preallocate arrays for results
    actual_trajectory = zeros(12, N_sim);
    control_inputs = zeros(4, N_sim);
    actual_trajectory_del1 = zeros(12, N_sim);
    actual_trajectory_del2 = zeros(12, N_sim);
    actual_trajectory_del3 = zeros(12, N_sim);
    control_inputs_del1 = zeros(4, N_sim);
    control_inputs_del2 = zeros(4, N_sim);
    control_inputs_del3 = zeros(4, N_sim);

    % Main simulation loop
    for i = 1:N_sim
        if norm(state_del1(4:5) - goal(1:2))< position_accuracy && norm(state_del2(4:5) - goal(1:2))< position_accuracy && norm(state_del3(4:5) - goal(1:2))< position_accuracy
            disp('Goal reached');
            actual_trajectory = actual_trajectory(:,1:i-1);
            actual_trajectory_del1 = actual_trajectory_del1(:,1:i-1);
            actual_trajectory_del2 = actual_trajectory_del2(:,1:i-1);
            actual_trajectory_del3 = actual_trajectory_del3(:,1:i-1);
            time_log = time_log(1:i-1);
            break;
        end
        % Current reference
        ref_traj = xr(1:12, i:min(i+N_horizon, N_sim));
        ref_traj_del1 = xr1(13:24, i:min(i+N_horizon, N_sim));
        ref_traj_del2 = xr2(13:24, i:min(i+N_horizon, N_sim));
        ref_traj_del3 = xr3(13:24, i:min(i+N_horizon, N_sim));

        % If reference is shorter than horizon, pad with last value
        if size(ref_traj, 2) < N_horizon+1
            ref_traj = [ref_traj, repmat(ref_traj(:,end), 1, N_horizon+1-size(ref_traj,2))];
            ref_traj_del1 = [ref_traj_del1, repmat(ref_traj_del1(:,end), 1, N_horizon+1-size(ref_traj_del1,2))];
            ref_traj_del2 = [ref_traj_del2, repmat(ref_traj_del2(:,end), 1, N_horizon+1-size(ref_traj_del2,2))];
            ref_traj_del3 = [ref_traj_del3, repmat(ref_traj_del3(:,end), 1, N_horizon+1-size(ref_traj_del3,2))];
        end

        % Get MPC control inputs
        [u_opt, x_pred] = quad_mpc(state, ref_traj, N_horizon);
        [u_opt_del1, x_pred] = quad_mpc(state_del1, ref_traj_del1, N_horizon);
        [u_opt_del2, x_pred] = quad_mpc(state_del2, ref_traj_del2, N_horizon);
        [u_opt_del3, x_pred] = quad_mpc(state_del3, ref_traj_del3, N_horizon);

        % Apply first control input
        u = u_opt(:,1);
        u_del1 = u_opt_del1(:,1);
        u_del2 = u_opt_del2(:,1);
        u_del3 = u_opt_del3(:,1);

        % Simulate system dynamics
        state_dot = quad_dynamics(state, u(2), u(3), u(4), u(1));
        state_dot_del1 = quad_dynamics(state_del1, u_del1(2), u_del1(3), u_del1(4), u_del1(1));
        state_dot_del2 = quad_dynamics(state_del2, u_del2(2), u_del2(3), u_del2(4), u_del2(1));
        state_dot_del3 = quad_dynamics(state_del3, u_del3(2), u_del3(3), u_del3(4), u_del3(1));

        % Euler integration
        state = state + dT * state_dot;
        state_del1 = state_del1 + dT * state_dot_del1;
        state_del2 = state_del2 + dT * state_dot_del2;
        state_del3 = state_del3 + dT * state_dot_del3;

        % Store results
        actual_trajectory(:,i) = state';
        control_inputs(:,i) = u;
        actual_trajectory_del1(:,i) = state_del1';
        control_inputs_del1(:,i) = u_del1;
        actual_trajectory_del2(:,i) = state_del2';
        control_inputs_del2(:,i) = u_del2;
        actual_trajectory_del3(:,i) = state_del3';
        control_inputs_del3(:,i) = u_del3;

        % Display progress
        if mod(i, 10) == 0
            fprintf('Simulation progress: %.1f%%\n', 100*i/N_sim);
        end
    end


    % Plot results
    plot_results(time_log, actual_trajectory, actual_trajectory_del1, actual_trajectory_del2, actual_trajectory_del3, vapf_del1, vapf_del2, vapf_del3);
end



% Function to generate APF trajectory (your existing code)
function [xr, vapf_del, time_log] = generate_apf_trajectory(initial_state_arg)
    param;

    state = initial_state_arg;
    time_log =  zeros(1, N);
    
    xr = zeros(24, N); %1st 12, normal apf, next 12, D-APF
    pos = [state(4); state(5); state(6); state(4); state(5); state(6)];
    vapf_del = zeros(3,N);
    for i=1:N
        artificial_pot = @(s) apf(s);
        l1= artificial_pot(pos);
        l2= artificial_pot(pos + 0.5*dT*l1);
        l3= artificial_pot(pos + 0.5*dT*l2);
        l4= artificial_pot(pos + dT*l3);
        pos = pos + (dT/6)*(l1 + 2*l2 + 2*l3 + l4);
        xr(:,i) = zeros(24,1);
        xr(4,i) = pos(1,1);
        xr(5,i) = pos(2,1);
        xr(6,i) = pos(3,1);
        xr(16,i) = pos(4,1);
        xr(17,i) = pos(5,1);
        xr(18,i) = pos(6,1);
        time_log(i)= i*dT;
        vapf_del(:,i) = l1(4:6,1);
    end
end

% Function to plot results
function plot_results(time, actual, actual_del1, actual_del2, actual_del3, vapf_del1, vapf_del2, vapf_del3)
    figure(1);
    param;
    % Position plots
    pos_err = actual(4:6, :) - goal(1:3)';
    size(time)
    size(vecnorm(pos_err))
    plot(time, vecnorm(pos_err),'r', 'LineWidth',1.5);
    

    figure(2);
    hold on;    
    % Plot the actual 3D trajectory
    plot3(actual(4,:), actual(5,:), actual(6,:), 'r', 'LineWidth', 2);    
    % Mark the start and goal
    scatter3(initial_state(4), initial_state(5), initial_state(6), 100, 'g', 'filled');
    scatter3(goal(1), goal(2), goal(3), 100, 'r', 'filled');
    for i = 1:Obstacle_count
        obs_x = squeeze(obstacle(i,:,1));
        obs_y = squeeze(obstacle(i,:,2));
        obs_z = squeeze(obstacle(i,:,3));
        fill3(obs_x, obs_y, obs_z, [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'k');
    end    
    % Aesthetics
    title('APF 3D Quadrotor Trajectory');
    legend('Actual', 'Start', 'Goal');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on;
    axis equal;
    view(3);

    
    figure(3);
    subplot(3,1,1);
    plot(time, vecnorm(actual(1:3,:)), 'g','LineWidth', 1.5);
    
    title('attitude APF');
    subplot(3,1,2);
    plot(time, vecnorm(actual(7:9,:)),'g','LineWidth',1.5);
    
    title('Angular Velocity APF');
    subplot(3,1,3);
    plot(time, vecnorm(actual(10:12,:)), 'g','LineWidth',1.5);
    
    title('Translational Velocity APF');


    %DAPF plots
    figure(4);
    param;
    % Position plots
    subplot(3,1,1);
    plot(time, vecnorm(actual_del1(4:6,:) - goal(1:3)'), 'r','LineWidth', 1.5);
    title('Initial State 1 Position error norm(D-APF)');    
    subplot(3,1,2);
    plot(time, vecnorm(actual_del2(4:6,:) - goal(1:3)'), 'r','LineWidth', 1.5);
    title('Initial State 2 Position error norm(D-APF)');
    subplot(3,1,3);
    plot(time, vecnorm(actual_del3(4:6,:) - goal(1:3)'), 'r','LineWidth', 1.5);
    title('Initial State 3 Position error norm(D-APF)');
    

    figure(5);
    hold on;    
    % Plot the actual 3D trajectory
    plot3(actual_del1(4,:), actual_del1(5,:), actual_del1(6,:), 'r', 'LineWidth', 2); hold on;
    plot3(actual_del2(4,:), actual_del2(5,:), actual_del2(6,:), 'g', 'LineWidth', 2); hold on;
    plot3(actual_del3(4,:), actual_del3(5,:), actual_del3(6,:), 'b', 'LineWidth', 2); hold on;
    
    % Mark the start and goal
    scatter3(initial_state1(4), initial_state1(5), initial_state1(6), 100, 'yellow', 'filled');
    scatter3(initial_state2(4), initial_state2(5), initial_state2(6), 100, 'magenta', 'filled');
    scatter3(initial_state3(4), initial_state3(5), initial_state3(6), 100, 'cyan', 'filled');
    scatter3(goal(1), goal(2), goal(3), 100, 'r', 'filled');
    for i = 1:Obstacle_count
        obs_x = squeeze(obstacle(i,:,1));
        obs_y = squeeze(obstacle(i,:,2));
        obs_z = squeeze(obstacle(i,:,3));
        fill3(obs_x, obs_y, obs_z, [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'k');
    end    
    % Aesthetics
    title('D-APF 3D Quadrotor Trajectory');
    legend('Actual1','Actual2','Actual3', 'Start1','Start2','Start3', 'Goal');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on;
    axis equal;
    view(3);

    
    figure(6);
    subplot(3,1,1);
    plot(time, vecnorm(actual_del1(1:3,:)),'g', time, vecnorm(actual_del2(1:3,:)), 'b', time, vecnorm(actual_del3(1:3,:)), 'r');
    legend('Initial State1', 'Initial State2', 'Initial State3');
    title('attitude D-APF');
    subplot(3,1,2);
    plot(time, vecnorm(actual_del1(7:9,:)),'g', time, vecnorm(actual_del2(7:9,:)), 'b', time, vecnorm(actual_del3(7:9,:)), 'r');
    legend('Initial State1', 'Initial State2', 'Initial State3');
    title('Angular Velocity D-APF');
    subplot(3,1,3);
    plot(time, vecnorm(actual_del1(10:12,:)),'g', time, vecnorm(actual_del2(10:12,:)), 'b', time, vecnorm(actual_del3(10:12,:)), 'r');
    legend('Initial State1', 'Initial State2', 'Initial State3');
    title('Translational Velocity D-APF');

    figure(7);
    plot(time, vecnorm(vapf_del1(1:3,:)),'g', time, vecnorm(vapf_del2(1:3,:)), 'b', time, vecnorm(vapf_del3(1:3,:)), 'r');
    legend('Initial State1', 'Initial State2', 'Initial State3');    title('Vapf ');


end