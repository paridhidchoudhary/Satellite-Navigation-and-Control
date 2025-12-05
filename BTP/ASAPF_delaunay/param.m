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
initial_state = [phi, theta, psi, x, y, z, phi_d, theta_d, psi_d, x_d, y_d, z_d];
%initial_state = [phi, theta, psi, x, y, phi_d, theta_d, psi_d, x_d, y_d];
initial_derivative = zeros(12, 1);
dvort = 0.4;

g=10;
m=1.4;

%APF gains
ka = 20;
kr = 100;

% Goal position
x_goal = 2.5;
y_goal = 2.5;
z_ref = 10;
goal = [x_goal, y_goal, z_ref];
position_accuracy = 0.1;

% Sampling period
dT = 0.1;
t_total = 20;
N=t_total/dT;

ix= 0.08;
iy= 0.08;
iz= 0.08;

k=[0.001, 0.001, 0.001, 1, 1, 0 , 0.001 , 0.001, 0.001, 1, 1, 0];
    

% Generate obstacles
Obstacle_count = 8;
angles = linspace(0, 2*pi, 360)';
obstacle = zeros(Obstacle_count, length(angles), 2);
obstacle_02 = zeros(Obstacle_count, length(angles), 2);
obstacle_03 = zeros(Obstacle_count, length(angles), 2);
obstacle_04 = zeros(Obstacle_count, length(angles), 2);
obstacle_05 = zeros(Obstacle_count, length(angles), 2);
c = zeros(Obstacle_count,2);
r = zeros(Obstacle_count,1);
for i=1:Obstacle_count
    while 1
        c(i,:) = 4*rand(1,2) - 1;
        r(i) = 0.25*rand() + 0.15;

        if norm([x y] - c(i,:)) > (r(i) + 0.35) && norm([x_goal y_goal] - c(i,:)) > (r(i) + 0.35)
            if i == 1, break; end
            [idx, dist] = dsearchn([c(1:(i-1),1) c(1:(i-1),2)], c(i,:));
            if dist > (r(idx)+r(i)+0.1)
                break;
            end
        end
    end
    obstacle(i,:,:) = [r(i) * cos(angles)+c(i,1) r(i)*sin(angles)+c(i,2) ];
    r(i) = r(i) + 0.2;
    obstacle_02(i,:,:) = [r(i) * cos(angles)+c(i,1) r(i)*sin(angles)+c(i,2) ];
    r(i) = r(i) + 0.2;
    obstacle_04(i,:,:) = [r(i) * cos(angles)+c(i,1) r(i)*sin(angles)+c(i,2) ];
end
