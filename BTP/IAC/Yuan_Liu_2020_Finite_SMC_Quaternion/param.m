%% Chaser parameters
J_c = diag([20, 20, 15]);
m_c = 90;


%% Target parameters
J_t = diag([20,20,15]);
m_t = 90;
sigma_t = [5 ; 0 ; 0];   % docking port position


%% Control Gains
k1 = 2; 
k2 = 2;
k3 = 0.5;
k4 = 0.1;
k5 = 0.05;
k6 = 4;
k7 = 10;
k8 = 10;



J2 = 1.08263e-3; % Earth's J2
Re = 6378.137e3; % Earth's radius [m]
mu_earth = 3.986004418e14; % Earth's mu [m^3/s^2]