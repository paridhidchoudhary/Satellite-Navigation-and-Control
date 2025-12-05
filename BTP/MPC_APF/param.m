% Initial position and orientation 
x = 0;
y = 0;
z = 0;
phi = 0;
theta = 0;
psi = 0;
x_d = 10;
y_d = 10;
z_d = 0;
phi_d = 0.1;
theta_d = 0.1;
psi_d = 0.1;
initial_state1 = [phi, theta, psi, x, y, z, phi_d, theta_d, psi_d, x_d, y_d, z_d];
initial_state = initial_state1;
initial_state2 = [phi, theta, psi, x+20, y+20, z+20, phi_d, theta_d, psi_d, x_d, y_d, z_d];
initial_state3 = [phi, theta, psi, x+15, y+15, z+15, phi_d, theta_d, psi_d, x_d, y_d, z_d];
%initial_state = [phi, theta, psi, x, y, phi_d, theta_d, psi_d, x_d, y_d];
initial_derivative = zeros(12, 1);
dvort = 0.4;

G=10;
m=1.4;

%APF gains
ka = 20;
kr = 100;

% Goal position
x_goal = 10;
y_goal = 10;
z_ref = 10;
goal = [x_goal, y_goal, z_ref];
position_accuracy = 0.1;

% Sampling period
dT = 0.1;
t_total = 40;
N=t_total/dT;

ix= 0.08;
iy= 0.08;
iz= 0.08;

k=[0.001, 0.001, 0.001, 1, 1, 0 , 0.001 , 0.001, 0.001, 1, 1, 0];
    %% 
    

% Generate obstacles
% Obstacle_count = 8;
% angles = linspace(0, 2*pi, 360)';
% obstacle = zeros(Obstacle_count, length(angles), 2);
% obstacle_02 = zeros(Obstacle_count, length(angles), 2);
% obstacle_03 = zeros(Obstacle_count, length(angles), 2);
% obstacle_04 = zeros(Obstacle_count, length(angles), 2);
% obstacle_05 = zeros(Obstacle_count, length(angles), 2);
% c = zeros(Obstacle_count,2);
% r = zeros(Obstacle_count,1);
% for i=1:Obstacle_count
%     while 1
%         c(i,:) = 4*rand(1,2) - 1;
%         r(i) = 0.25*rand() + 0.15;
% 
%         if norm([x y] - c(i,:)) > (r(i) + 0.35) && norm([x_goal y_goal] - c(i,:)) > (r(i) + 0.35)
%             if i == 1, break; end
%             [idx, dist] = dsearchn([c(1:(i-1),1) c(1:(i-1),2)], c(i,:));
%             if dist > (r(idx)+r(i)+0.1)
%                 break;
%             end
%         end
%     end
%     obstacle(i,:,:) = [r(i) * cos(angles)+c(i,1) r(i)*sin(angles)+c(i,2) ];
%     r(i) = r(i) + 0.2;
%     obstacle_02(i,:,:) = [r(i) * cos(angles)+c(i,1) r(i)*sin(angles)+c(i,2) ];
%     r(i) = r(i) + 0.2;
%     obstacle_04(i,:,:) = [r(i) * cos(angles)+c(i,1) r(i)*sin(angles)+c(i,2) ];
% end
%% 
%% 


% Increase number of obstacles
Obstacle_count = 7;  % Increased number of obstacles

% Create meshgrid for sphere surface points
[phi, theta] = meshgrid(linspace(0, pi, 20), linspace(0, 2*pi, 36));
points_per_sphere = numel(phi);

% Initialize arrays for obstacles
obstacle = zeros(Obstacle_count, points_per_sphere, 3);
obstacle_02 = zeros(Obstacle_count, points_per_sphere, 3);
obstacle_04 = zeros(Obstacle_count, points_per_sphere, 3);

% Centers and radii
c = zeros(Obstacle_count, 3);
r = zeros(Obstacle_count, 1);

% Define workspace bounds based on start/goal points with some margin
x_bound = [-10, max(initial_state1(4), goal(1)) + 10];
y_bound = [-10, max(initial_state1(5), goal(2)) + 10];
z_bound = [-10, max(initial_state1(6), goal(3)) + 10];

for i = 1:Obstacle_count
    attempts = 0;
    max_attempts = 100;

    while attempts < max_attempts
        % Generate random position
        c(i,:) = [
            x_bound(1) + (x_bound(2)-x_bound(1))*rand()
            y_bound(1) + (y_bound(2)-y_bound(1))*rand()
            z_bound(1) + (z_bound(2)-z_bound(1))*rand()
        ];

        % Larger radius range appropriate for the workspace size
        r(i) = 3*rand() + 2;  % Radius between 2 and 5 units

        % Check collision with start and goal positions
        start_clear = norm([initial_state(4) initial_state(5) initial_state(6)] - c(i,:)) > (r(i) + 5);  % Larger clearance
        goal_clear = norm([goal(1) goal(2) goal(3)] - c(i,:)) > (r(i) + 5);

        if start_clear && goal_clear
            if i == 1
                break;
            end

            % Check collision with other obstacles
            distances = sqrt(sum((c(1:i-1,:) - c(i,:)).^2, 2));
            min_dist = min(distances);

            % Increased minimum spacing between obstacles
            if min_dist > (max(r(1:i-1)) + r(i) + 3)
                break;
            end
        end

        attempts = attempts + 1;
    end

    if attempts >= max_attempts
        warning('Could not place obstacle %d after %d attempts', i, max_attempts);
        Obstacle_count = i-1;
        break;
    end

    % Generate sphere surface points
    x_sphere = r(i) * sin(phi) .* cos(theta);
    y_sphere = r(i) * sin(phi) .* sin(theta);
    z_sphere = r(i) * cos(phi);

    % Store points for base obstacle
    obstacle(i,:,:) = [reshape(x_sphere + c(i,1), 1, []);
                      reshape(y_sphere + c(i,2), 1, []); 
                      reshape(z_sphere + c(i,3), 1, [])]';

    % Generate larger shells with increased shell thickness
    r_outer = r(i) + 2;  % Thicker shells
    x_sphere = r_outer * sin(phi) .* cos(theta);
    y_sphere = r_outer * sin(phi) .* sin(theta);
    z_sphere = r_outer * cos(phi);

    obstacle_02(i,:,:) = [reshape(x_sphere + c(i,1), 1, []);
                         reshape(y_sphere + c(i,2), 1, []); 
                         reshape(z_sphere + c(i,3), 1, [])]';

    r_outer = r_outer + 2;
    x_sphere = r_outer * sin(phi) .* cos(theta);
    y_sphere = r_outer * sin(phi) .* sin(theta);
    z_sphere = r_outer * cos(phi);

    obstacle_04(i,:,:) = [reshape(x_sphere + c(i,1), 1, []);
                         reshape(y_sphere + c(i,2), 1, []); 
                         reshape(z_sphere + c(i,3), 1, [])]';
end

