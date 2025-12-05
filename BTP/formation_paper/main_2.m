params;

% Initialize arrays for 3 satellites
eta_tilde = zeros(6, N, 3);
xi_tilde = zeros(6, N, 3);
eta=zeros(6,N,3);
xi=zeros(6,N,3);

% Define desired relative configurations for each spacecraft
% eta_desired = zeros(6, 3);
% eta_desired(:, 1) = [0; 0; 0; -1000; 1500; -500];  % Spacecraft 1
% eta_desired(:, 2) = [0; 0; 0; -350; 1500; 0];      % Spacecraft 2
% eta_desired(:, 3) = [0; 0; 0; 350; -1500; 500];    % Spacecraft 3

% Initial conditions for each satellite (as given in the paper)
eta_tilde(:,1,1) = [30*pi/180; 0; 0; -1000; 1500; -500];
eta_tilde(:,1,2) = [60*pi/180; 0; 0; -350; 1500; 0];
eta_tilde(:,1,3) = [60*pi/180; 0; 0; 350; -1500; 500];

xi_tilde(:,1,1) = [0.001*[-0.06; 0.26; 0.05]; 3.17; 1.24; 1.82];
xi_tilde(:,1,2) = [0.001*[-0.06; 0.12; -0.099]; -0.16; -0.24; 3.65];
xi_tilde(:,1,3) = [0.001*[-0.06; 0.12; -0.099]; -1.77; -0.28; -1.51];

%make remaining 3 columns of eta_tilde and xi_tilde 0
% for i=1:3 %loop over all the satellites
%     for k=2:4 %make remaining column 0
%         eta_tilde(:,k,1,i)=[0;0;0;0;0;0];
%         xi_tilde(:,k,1,i)=[0;0;0;0;0;0];
%     end
% end
%xi_0= zeros(6,4);
xi_0 = [0.001; 0.001; 0.001; 0.1; 0.1; 0.1];
% for i=2:4
%     xi_0(:,i)=[0;0;0;0;0;0];
% end
R=eye(3);
b=zeros(3,1);
g0=[R, b; 0,0,0,1];
tspan=[0,dt];

s_values = cell(1, 3);

% %for k = 1:N-1
%     % Loop over 3 satellites
% for i = 1:3
%     initial_state = [xi_tilde(:,1,i); eta_tilde(:,1,i); xi_0];
%     size(initial_state)
%     [t, state] = ode45(@(t, state) dynamics(t, state, i, dt), tspan, initial_state);
% 
%     % Store state values
%     xi_tilde(:,1:length(t),i) = state(:,1:6)';
%     eta_tilde(:,1:length(t),i) = state(:,7:12)';
% 
%     % Calculate sliding surface and norms
%     s_temp = zeros(1, length(t));
%     norm_eta = zeros(1, length(t));
%     norm_xi = zeros(1, length(t));
%     for j = 1:length(t)
%         xi_tilde_current = state(j, 1:6)';
%         eta_tilde_current = state(j, 7:12)';
%         s_temp(j) = norm(sliding_surface(eta_tilde_current, xi_tilde_current));
%         norm_eta(j) = norm(eta_tilde_current);
%         norm_xi(j) = norm(xi_tilde_current);
%     end
%     s_values{i} = s_temp;
% 
%     % Plot sliding variable vs. time for this satellite
%     figure;
%     %subplot(2,1,1);
%     plot(t, s_temp, 'LineWidth', 1.5);
%     hold on;
%     % xlabel('Time (s)');
%     % ylabel('Sliding Variable (s)');
%     % title(['Sliding Variable vs. Time for Satellite ', num2str(i)]);
% 
%     % Plot norm of eta vs. norm of xi for this satellite
%     %subplot(2,1,2);
% 
%     plot(norm_eta, norm_xi, 'LineWidth', 1.5);
%     % xlabel('||\eta||');
%     % ylabel('||\xi||');
%     % title(['Norm of \eta vs. Norm of \xi for Satellite ', num2str(i)]);
%     hold off;
%     grid on;
% end
%     hold off;
% %end
% Loop over 3 satellites (only plotting for satellite 3)
for i = 1:3
    % Initial state setup
    initial_state = [xi_tilde(:, 1, i); eta_tilde(:, 1, i); xi_0];
    [t, state] = ode45(@(t, state) dynamics(t, state, i, dt), tspan, initial_state);
    
    % Store state values
    xi_tilde(:, 1:length(t), i) = state(:, 1:6)';
    eta_tilde(:, 1:length(t), i) = state(:, 7:12)';
    
    % Initialize variables to store sliding surface and norms
    s_temp = zeros(1, length(t));
    norm_eta = zeros(1, length(t));
    norm_xi = zeros(1, length(t));
    
    % Compute sliding surface and norms for each time step
    for j = 1:length(t)
        xi_tilde_current = state(j, 1:6)';
        eta_tilde_current = state(j, 7:12)';
        s_temp(j) = norm(sliding_surface(eta_tilde_current, xi_tilde_current));
        norm_eta(j) = norm(eta_tilde_current);
        norm_xi(j) = norm(xi_tilde_current);
    end
    
    % Store s_temp for satellite 3 only
    if i == 3
        % Create a new figure for satellite 3 plots
        figure;
        
        % Plot sliding surface vs time
        subplot(2, 1, 1);
        plot(t, s_temp, 'LineWidth', 1.5, 'Color', 'b');
        xlabel('Time (s)');
        ylabel('Sliding Surface (s)');
        title('Sliding Surface vs. Time for Satellite 3');
        grid on;
        
        % Plot norm of xi vs norm of eta
        subplot(2, 1, 2);
        plot(norm_eta, norm_xi, 'LineWidth', 1.5, 'Color', 'r');
        xlabel('||\eta||');
        ylabel('||\xi||');
        title('Norm of \eta vs. Norm of \xi for Satellite 3');
        grid on;
        
        % Exit loop after plotting for satellite 3
        break;
    end
end


% Plotting results for 3 satellites
figure;
subplot(2,1,1);
for i = 1:3
    norm_att=zeros(1,N);
    for j=1:N
        norm_att(j)=norm(eta_tilde(1:3,j,i));
    end
    plot(1:N,norm_att);
    hold on;
end
title('Attitude Error');
xlabel('Time(s)');
ylabel('Error(rad)');
legend('\theta_1', '\theta_2', '\theta_3' );
hold off;

subplot(2,1,2);
for i = 1:3
    norm_pos=zeros(1,N);
    for j=1:N
        norm_pos(j)=norm(eta_tilde(4:6,j,i));
    end
    plot(1:N,norm_pos);
    hold on;
end
title('Position Error');
xlabel('Time(s)');
ylabel('Error(m)');
legend('x1', 'x2', 'x3');
hold off;

figure;
subplot(2,1,1);
for i = 1:3
    norm_ang=zeros(1,N);
    for j=1:N
        norm_ang(j)=norm(xi_tilde(1:3,j,i));
    end
    plot(1:N,norm_ang);
    hold on;
end
title('Angular velocity error');
xlabel('Time(s)');
ylabel('Error (rad/s)');
legend('\omega_x1', '\omega_y1', '\omega_z1', '\omega_x2', '\omega_y2', '\omega_z2', '\omega_x3', '\omega_y3', '\omega_z3');
hold off;

subplot(2,1,2);
for i = 1:3
    norm_vel=zeros(1,N);
    for j=1:N
        norm_vel(j)=norm(xi_tilde(4:6,j,i));
    end
    plot(1:N,norm_vel);
    hold on;
end
title('velocity error');
xlabel('Time(s)');
ylabel('Error (m/s)');
legend('\v_x1', '\v_y1', '\v_z1', '\v_x2', '\v_y2', '\v_z2', '\v_x3', '\v_y3', '\v_z3');
hold off;


% Add a plot for the formation shape
% figure;
% for k = 1:N
%     plot3(squeeze(eta_tilde(4,k,:)), squeeze(eta_tilde(5,k,:)), squeeze(eta_tilde(6,k,:)), 'o');
%     hold on;
%     plot3(squeeze(eta_tilde(4,k,:)), squeeze(eta_tilde(5,k,:)), squeeze(eta_tilde(6,k,:)));
%     title(['Formation at t = ', num2str((k-1)*dt), ' s']);
%     xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
%     %axis equal;
%     grid on;
%     if k == 1 || k == N
%         pause(1);
%     else
%         pause(0.01);
%     end
%     if k < N
%         clf;
%     end
% end
% Define Molniya orbit parameters
eccentricity = 0.75;          % Eccentricity for the ellipse
semi_major_axis = 26600e3;     % Semi-major axis in meters (approximate for Molniya)
inclination = 63.4;            % Inclination of Molniya orbit in degrees

% Calculate semi-minor axis from eccentricity and semi-major axis
semi_minor_axis = semi_major_axis * sqrt(1 - eccentricity^2);

% Parametric angle for the ellipse
theta = linspace(0, 2 * pi, 1000);

% 2D coordinates in the orbital plane
x = semi_major_axis * cos(theta) - semi_major_axis * eccentricity; % Center the ellipse
y = semi_minor_axis * sin(theta);

% Rotate the orbit by the inclination angle to project into 3D
inclination_rad = deg2rad(inclination);
z = y * sin(inclination_rad);
y = y * cos(inclination_rad);

% Plotting both the eta_tilde trajectory and the ellipse in the same figure
figure;

% Plot eta_tilde trajectory
plot3(squeeze(eta_tilde(4,k,:)), squeeze(eta_tilde(5,k,:)), squeeze(eta_tilde(6,k,:)), 'o');
hold on;
plot3(squeeze(eta_tilde(4,k,:)), squeeze(eta_tilde(5,k,:)), squeeze(eta_tilde(6,k,:)), 'b');
hold on;

% Plot the Molniya-like orbit ellipse
plot3(x, y, z, 'r', 'LineWidth', 1.5); % Using red color for the orbit line

% Plot the focus (representing Earth at one of the foci)
plot3(-semi_major_axis * eccentricity, 0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% Customize plot
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('3D Trajectory of eta\_tilde and Molniya-like Orbit');
grid on;
axis equal;
legend('eta\_tilde points', 'eta\_tilde trajectory', 'Molniya-like orbit', 'Earth (focus)');

hold off;

% Plot the trajectories of the three spacecraft
% for i = 1:3
%     plot3(squeeze(eta_tilde(4, :, i)), squeeze(eta_tilde(5, :, i)), squeeze(eta_tilde(6, :, i)), 'o-', 'DisplayName', ['Spacecraft ', num2str(i)]);
% end
% 
% title('Spacecraft Formation and Convergence to Elliptical Orbit');
% xlabel('X (m)');
% ylabel('Y (m)');
% zlabel('Z (m)');
% legend;
% grid on;
% axis equal;

