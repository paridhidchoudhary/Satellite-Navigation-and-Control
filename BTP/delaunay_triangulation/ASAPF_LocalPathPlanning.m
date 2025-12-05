function [nablaU] = ASAPF_LocalPathPlanning(x,y,z, ASAPF, obstacle, c, r)
    
    % Goal position
    x_goal = 0;
    y_goal = 0;
    z_goal= 0;
    position_accuracy = 0.1;
    
    % Sampling period
    dT = 0.1;
    N=10;
    
    % Generate obstacles 

    [nablaU] = AdaptiveSafeArtificialPotentialField(0.4, 4e0 *dT, x, y,z, x_goal, y_goal, z_goal, position_accuracy, obstacle,c,r, dT, N, ASAPF);
    %simout2 = ASAPF_using_MPC(0.4, 4e0 *dT, x, y, theta, x_goal, y_goal, position_accuracy, obstacle, dT, simTimeMax, N);

    
end

function [nablaU] = AdaptiveSafeArtificialPotentialField(dvort, adaptationGain, x, y, z, x_goal, y_goal, z_goal, position_accuracy, obstacle,c,r, dT, N, ASAPF)

    ASAPF.X = x;
    ASAPF.Y = y;
    ASAPF.Z = z;
    ASAPF.dvort= dvort;
    
    Obstacle_count = size(obstacle, 1);

    % Flatten obstacle array to 3D points and do Delaunay triangulation
    %obst_points = reshape(obstacle, [], 3);  % Reshape to 2D array of points (each row is an obstacle's [x, y, z])
    delaunayTri = delaunayTriangulation(c(:,1), c(:,2), c(:,3));  % Create triangulation
    neighbors = delaunayTri.vertexAttachments;  % Get neighbors for each obstacle

    % Calculate Attractive Potential
    if norm([x y z] - [x_goal y_goal z_goal]) <= ASAPF.dstar
        nablaU_att = ASAPF.zeta * ([x y z] - [x_goal y_goal z_goal]);
    else
        nablaU_att = (ASAPF.dstar / norm([x y z] - [x_goal y_goal z_goal])) * ASAPF.zeta * ([x y z] - [x_goal y_goal z_goal]);
    end

    % Initialize repulsive potential
    nablaU_obst = [0 0 0];

    % Loop through each obstacle
    for i = 1:Obstacle_count
        % Get neighbors for the current obstacle point
        neighbor_indices = neighbors{i};

        % Validate indices (ensure they're within bounds)
        %valid_neighbors = neighbor_indices(neighbor_indices <= size(obst_points, 1));
        valid_neighbors = neighbor_indices(neighbor_indices <= Obstacle_count);
        if isempty(valid_neighbors)
            continue; % Skip if no valid neighbors
        end

        % Compute distances to the valid neighbors
        %distances = vecnorm(obst_points(valid_neighbors, :) - [x, y, z], 2, 2);
        distances = vecnorm(c(valid_neighbors,:) - [x, y, z], 2, 2);

        % Find nearest two neighbors
        [~, sorted_idx] = sort(distances);
        nearest_neighbors = valid_neighbors(sorted_idx(1:min(2, numel(sorted_idx))));
        

        % Calculate repulsive potential only from the nearest two obstacles
        for neighbor_idx = nearest_neighbors(:)'
        
            neighbor_pos = c(neighbor_idx, :);
            radius_obst = r(neighbor_idx);
            dist_to_neighbor = norm([x, y, z] - neighbor_pos);

            
            if dist_to_neighbor <= radius_obst
                nablaU_rep_Oi = -ASAPF.eta*(1/radius_obst - 1/ASAPF.Qstar)*(1/radius_obst^2)*(neighbor_pos - [x y z])/dist_to_neighbor;
            elseif radius_obst< dist_to_neighbor < ASAPF.Qstar
                nablaU_rep_Oi = -ASAPF.eta*(1/dist_to_neighbor - 1/ASAPF.Qstar) *(neighbor_pos - [x y z])/(dist_to_neighbor^3);
            else
                nablaU_rep_Oi = [0 0 0];
            end
            nablaU_obst= nablaU_obst + nablaU_rep_Oi;
            
        end
    end

    
    % Calculate final potential
    nablaU = (nablaU_att + nablaU_obst);
    nablaU_diff =(nablaU_att - nablaU_obst);
    
    end




