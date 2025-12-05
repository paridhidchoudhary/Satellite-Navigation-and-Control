% Follower Spacecraft mass
m = 50;  % kg

% Leader Spacecraft mass
m0 = 50;  % kg

% Follower Inertia matrix
J = diag([48.5 , 51.0 , 47.6]);

% Leader Inertia matrix
J0 = diag([48.5 , 51.0 , 47.6]);

%% Maximum magnitude of control thrust
F_max = 500000;  % Newton

%% Maximum magnitude of control moment
M_max = 80;  % Newton.meter

%% Control Parameters

c1 = 0.96;       % c1 and c2 are constants of C matrix
c2 = 0.96;

p1 = 0.055;       % p1 and p2 are constants of F matrix
p2 = 0.0275;

k1 = 0.5;          % k1 and k2 are constants of K matrix
k2 = 0.5;


%% Standard Parameters
J2 = 1.08263e-3;
mu = 3.986e14;
Re = 6378e3;

%% Leaders orbital elements parameters
a = 26628*1e3;               % semi major axis
e = 0.7417;                   % eccentricity
i_ang = 63.4*(pi/180);        % inclination
omega = 270*(pi/180);         % Argument of perigee
RAAN = 0*(pi/180);            % Right Ascesion of Ascending Node
TA = 120*(pi/180);            % True Anomaly
hs = sqrt(mu*a*(1-(e^2)));    % Angular Momentum
Tp = 2*pi*sqrt((a^3)/mu);     % Time Period of spacecraft

%% % Generate obstacles
% Number of obstacles
Obstacle_count = 8;

% Define the angular steps (for full 360° rotation)
theta = linspace(0, 2*pi, 360); % Azimuthal angle (0 to 2*pi)
phi = linspace(0, pi, 180); % Polar angle (0 to pi)

% Initialize obstacle arrays
obstacle = zeros(Obstacle_count, length(theta), length(phi), 3); % Store x, y, z for each obstacle

% Random centers and radii for obstacles
c = zeros(Obstacle_count, 3); % 3D centers (x, y, z)
r = zeros(Obstacle_count, 1); % Radii


%% 
% Define start and goal points
%x = 33.4582; y = 37.7976; z = 72.6427;  % Start point
x = 50; y = 50; z = 50;  % Start point
x_goal = 0; y_goal = 0; z_goal = 0;     % Goal point

% Increase number of obstacles
Obstacle_count = 30;  % Increased number of obstacles

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
x_bound = [-10, max(x, x_goal) + 10];
y_bound = [-10, max(y, y_goal) + 10];
z_bound = [-10, max(z, z_goal) + 10];

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
        start_clear = norm([x y z] - c(i,:)) > (r(i) + 5);  % Larger clearance
        goal_clear = norm([x_goal y_goal z_goal] - c(i,:)) > (r(i) + 5);
        
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


