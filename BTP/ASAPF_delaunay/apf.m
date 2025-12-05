function [vapf,vrep, vatt, psi_ref] = apf(state)
param;
vatt_x = -ka*(state(4) - goal(1));
vatt_y = -ka*(state(5) - goal(2));
vatt(1) = vatt_x;
vatt(2) = vatt_y;

Obstacle_count = size(obstacle,1);
obst_idx = zeros(1,Obstacle_count);
obst_dist = zeros(1,Obstacle_count);
r0= 5;
vrep_x = 0;
vrep_y = 0;
x = state(4);
y = state(5);
for i=1:Obstacle_count
    P = [obstacle(i,:,1)' obstacle(i,:,2)'];
    [obst_idx(i), obst_dist(i)] = dsearchn(P, [x y]);
    xo = P(obst_idx(i), 1);
    yo = P(obst_idx(i), 2);
    xor = state(4) - xo;
    yor = state(5) - yo;
    dist_vec = [xor, yor];
    dist_mag = norm(dist_vec);
    if dist_mag<=r0
        vrep_x = vrep_x - kr*(1 - dist_mag/r0)*xor/(dist_mag^3);
        vrep_y = vrep_y - kr*(1 - dist_mag/r0)*yor/(dist_mag^3);
    else
        vrep_x = vrep_x;
        vrep_y = vrep_y;
    end
end
vrep(1) = vrep_x;
vrep(2) = vrep_y;
vapf(1) = vatt_x + vrep_x;
vapf(2) = vatt_y + vrep_y;
vatt
vrep
psi_ref = atan2(vapf(2), vapf(1));
end