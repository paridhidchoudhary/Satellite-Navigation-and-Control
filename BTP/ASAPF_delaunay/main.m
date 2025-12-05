param;

state = initial_state;
state_dot = initial_derivative;
trajectory = zeros(12, N);
time_log =  zeros(1, N);

% Initialize logging arrays
error_x = [];
error_y = [];
error_z = [];

vel_x_log = [];
vel_y_log = [];
vel_z_log = [];

omega_x_log = [];
omega_y_log = [];
omega_z_log = [];

vattx_log = [];
vatty_log = [];
vrepx_log = [];
vrepy_log = [];


for i = 1:1:N
    %fprintf('norm error: %0.4f \n', norm(state(4:5)- goal(1:2)));
    if norm(state(4:5) - goal(1:2))< position_accuracy
        disp('Goal reached');
        trajectory = trajectory(:,1:i-1);
        time_log = time_log(1:i-1);
        break;
    end
    x = state(4);
    y = state(5);
    x_goal = goal(1);
    y_goal = goal(2);
    % Position errors
    error_x(end+1) = x_goal - x;
    error_y(end+1) = y_goal - y;
    error_z(end+1) = goal(3) - state(6);  % if you're working in 3D
    
    % Velocities (assumed from control logic or dynamics)
    vel_x_log(end+1) = state(10);
    vel_y_log(end+1) = state(11);
    vel_z_log(end+1) = state(12);  % if simulated
    
    % Angular velocities (based on attitude control or from sim state)
    omega_x_log(end+1) = state(7);
    omega_y_log(end+1) = state(8);
    omega_z_log(end+1) = state(9);

    %disp('now calc simout \n');
    %[simout, nablaU, v_x, v_y, psi_ref] = ASAPF_LocalPathPlanning(dvort, 4e0 *dT, x, y, psi, x_goal, y_goal, position_accuracy, obstacle,r,c);
    [vapf, vrep, vatt, psi_ref] = apf(state);
    vapf_x(end+1) = vapf(1);
    vapf_y(end+1) = vapf(2);
    vattx_log(end+1) = vatt(1);
    vatty_log(end+1) = vatt(2);
    vrepx_log(end+1) = vrep(1);
    vrepy_log(end+1) = vrep(2);
    %disp(nablaU);
    % [taux, tauy, tauz, T] = pd_control(v_x, v_y, psi_ref, goal(3), state, state_dot, k);
    [taux, tauy, tauz, T] = pd_control(vapf(1), vapf(2), goal(3), psi_ref, state, state_dot, k);
    %fprintf("v_x: %.4f, v_y: %.4f, taux: %.4f, tauy: %.4f, tauz: %.4f, T: %.4f\n", v_x, v_y, taux, tauy, tauz, T);

    dynamics =@(s) quad_dynamics(s, taux, tauy, tauz, T);
    %disp('now calculating k terms\n');
    l1= dynamics(state);
    l2= dynamics(state + 0.5*dT*l1);
    l3= dynamics(state + 0.5*dT*l2);
    l4= dynamics(state + dT*l3);
    state = state + (dT/6)*(l1 + 2*l2 + 2*l3 + l4);
    %state_dot = l1;

    trajectory(:,i) = state;
    time_log(i)= i*dT;
end


% Plot it
figure(1); subplot(1,1,1);
cla; hold on; grid on; box on;
daspect([1 1 1]); 
xlim([-1,4]);  ylim([-1 3]); xlabel("x [m]"); ylabel("y [m]");
box on; hold on;
plot(trajectory(4,:), trajectory(5,:), 'Color',[0 0.4470 0.7410], 'LineWidth', 1); % Plot traveled path
plot(goal(1), goal(2), 'xg');
for i=1:Obstacle_count
    plot(obstacle(i,:,1), obstacle(i,:,2), '-r');
    plot(obstacle_02(i,:,1), obstacle_02(i,:,2), '--r');
    plot(obstacle_04(i,:,1), obstacle_04(i,:,2), '--b');
end
drawnow;
title('ASAPF')


figure(2);

subplot(3,1,1);
plot(time_log, error_x, 'r', 'LineWidth', 1.5); hold on;
plot(time_log, error_y, 'g', 'LineWidth', 1.5);
plot(time_log, error_z, 'b', 'LineWidth', 1.5);
title('Position Error'); xlabel('Time Step'); ylabel('Error (m)');
legend('X error','Y error','Z error'); grid on;

subplot(3,1,2);
plot(time_log, vel_x_log, 'r', 'LineWidth', 1.5); hold on;
plot(time_log, vel_y_log, 'g', 'LineWidth', 1.5);
plot(time_log, vel_z_log, 'b', 'LineWidth', 1.5);
title('Velocity Components'); xlabel('Time Step'); ylabel('Velocity (m/s)');
legend('Vx','Vy','Vz'); grid on;

subplot(3,1,3);
plot(time_log, rad2deg(omega_x_log), 'r', 'LineWidth', 1.5); hold on;
plot(time_log, rad2deg(omega_y_log), 'g', 'LineWidth', 1.5);
plot(time_log, rad2deg(omega_z_log), 'b', 'LineWidth', 1.5);
title('Angular Velocity'); xlabel('Time Step'); ylabel('\omega (deg/s)');
legend('\omega_x','\omega_y','\omega_z'); grid on;

figure(3);

subplot(2,1,1);
plot(time_log, vattx_log, 'r', 'LineWidth', 1.5); hold on;
plot(time_log, vatty_log, 'g', 'LineWidth', 1.5);
title('APF attractive Force'); xlabel('Time Step'); ylabel('force');
legend('x comp', 'y comp'); grid on;

subplot(2,1,2);
plot(time_log, vrepx_log, 'r', 'LineWidth', 1.5); hold on;
plot(time_log, vrepy_log, 'g', 'LineWidth', 1.5);
title('APF repulsive Force'); xlabel('Time Step'); ylabel('force');
legend('x comp', 'y comp'); grid on;


   

