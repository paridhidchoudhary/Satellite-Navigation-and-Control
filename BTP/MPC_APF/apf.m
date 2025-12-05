function vapf = apf(pos)
param;
x = pos(1);
y = pos(2);
z = pos(3);
vatt_x = -ka*(x - goal(1));
vatt_y = -ka*(y - goal(2));
vatt_z = -ka*(z - goal(3));

Obstacle_count = size(obstacle,1);
obst_idx = zeros(1,Obstacle_count);
obst_dist = zeros(1,Obstacle_count);
r0= 5;
vrep_x = 0;
vrep_y = 0;
vrep_z = 0;

vrep_x_del = 0;
vrep_y_del = 0;
vrep_z_del = 0;

%delaunay triangulation
delaunayTri = delaunayTriangulation(c(:,1), c(:,2), c(:,3));  % Create triangulation
neighbors = delaunayTri.vertexAttachments;  % Get neighbors for each obstacle

for i=1:Obstacle_count
    if norm([x y z] - [goal(1) goal(2) goal(3)]) <= position_accuracy
        vapf(1,1) = vatt_x + vrep_x;
        vapf(2,1) = vatt_y + vrep_y;
        vapf(3,1) = vatt_z + vrep_z;
        
        vapf(4,1) = vatt_x + vrep_x_del;
        vapf(5,1) = vatt_y + vrep_y_del;
        vapf(6,1) = vatt_z + vrep_z_del;
        break;
    end
    P = [obstacle(i,:,1)' obstacle(i,:,2)' obstacle(i,:,3)'];
    [obst_idx(i), obst_dist(i)] = dsearchn(P, [x y z]);
    xo = P(obst_idx(i), 1);
    yo = P(obst_idx(i), 2);
    zo = P(obst_idx(i), 3);
    xor = x - xo;
    yor = y - yo;
    zor = z - zo;
    dist_vec = [xor, yor, zor];
    dist_mag = norm(dist_vec);
    if dist_mag<=r0
        vrep_x = vrep_x - kr*(1 - dist_mag/r0)*xor/(dist_mag^3);
        vrep_y = vrep_y - kr*(1 - dist_mag/r0)*yor/(dist_mag^3);
        vrep_z = vrep_z - kr*(1 - dist_mag/r0)*zor/(dist_mag^3);
    end

    %D-APF Computation
    neighbor_indices = neighbors{i};
    valid_neighbors = neighbor_indices(neighbor_indices <= Obstacle_count);
    if isempty(valid_neighbors)
        continue; % Skip if no valid neighbors
    end
    distances = vecnorm(c(valid_neighbors,:) - [x, y, z], 2, 2);
    % Find nearest two neighbors
    [~, sorted_idx] = sort(distances);
    nearest_neighbors = valid_neighbors(sorted_idx(1:min(2, numel(sorted_idx))));
    for neighbor_idx = nearest_neighbors(:)'
        neighbor_pos = c(neighbor_idx, :);
        xor_del = x - neighbor_pos(1);
        yor_del = y - neighbor_pos(2);
        zor_del = z - neighbor_pos(3);
        radius_obst = r(neighbor_idx);
        dist_to_neighbor = norm([x, y, z] - neighbor_pos)-radius_obst;

        if dist_to_neighbor<=r0
            vrep_x_del = vrep_x_del - kr*(1 - dist_to_neighbor/r0)*xor_del/(dist_to_neighbor^3);
            vrep_y_del = vrep_y_del - kr*(1 - dist_to_neighbor/r0)*yor_del/(dist_to_neighbor^3);
            vrep_z_del = vrep_z_del - kr*(1 - dist_to_neighbor/r0)*zor_del/(dist_to_neighbor^3);
        end
    end


vapf(1,1) = vatt_x + vrep_x;
vapf(2,1) = vatt_y + vrep_y;
vapf(3,1) = vatt_z + vrep_z;

vapf(4,1) = vatt_x + vrep_x_del;
vapf(5,1) = vatt_y + vrep_y_del;
vapf(6,1) = vatt_z + vrep_z_del;
end