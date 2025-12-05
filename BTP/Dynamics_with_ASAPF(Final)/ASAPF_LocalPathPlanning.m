function [simout, nablaU, theta_dot] = ASAPF_LocalPathPlanning(x,y,theta,z, ASAPF, obstacle)
    
    % Goal position
    x_goal = 0;
    y_goal = 0;
    z_goal= 0;
    position_accuracy = 0.1;
    
    % Sampling period
    dT = 0.1;
    N=10;
    
    % Generate obstacles 

    [simout, nablaU, theta_dot] = AdaptiveSafeArtificialPotentialField(0.4, 4e0 *dT, x, y,z, theta, x_goal, y_goal, z_goal, position_accuracy, obstacle, dT, N, ASAPF);
    %simout2 = ASAPF_using_MPC(0.4, 4e0 *dT, x, y, theta, x_goal, y_goal, position_accuracy, obstacle, dT, simTimeMax, N);

    
end

function [t, v, w] = CalculateVelocities(simout)

    t = simout.t;
    x = simout.X(1:simout.t); 
    y = simout.Y(1:simout.t);
    theta = simout.Theta(1:simout.t);

    dT = simout.dT;
    v = sqrt( (x(2:end)-x(1:(end-1))).^2 + (y(2:end)-y(1:(end-1))).^2) / dT;
    w = (theta(2:end)-theta(1:(end-1))) / dT;
    
    t = dT*(1:(t+1));
    v = [0 v 0];
    w = [0 w 0];
end
function MPC = ASAPF_using_MPC(dvort, adaptationGain, x, y, theta, x_goal, y_goal, position_accuracy, obstacle, dT, simTimeMax, N)
   
    % SAPF parameters
    MPC.dstar = 0.3;
    MPC.Qstar = 1.5;
    MPC.dsafe = 0.2;
    MPC.dvort = dvort;
    MPC.alpha_th = deg2rad(5);

    % Parameters related to kinematic model
    MPC.error_theta_max = deg2rad(45);
    MPC.v_max = 0.2;
    MPC.a_max = 0.2;
    MPC.omega_max = 0.5*pi; 

    %MPC Objective function constants:
    MPC.Q_pos = diag([10, 10]);    % Position error weight
    MPC.R_input = diag([1, 0.1]);  % Control input weight
    %MPC.Q_obs = 100;
    MPC.Q_theta=40;

    [MPC.zeta, MPC.eta_max, MPC.eta_min] = ASAPF_CalculateParameters(MPC.Qstar, MPC.dstar, MPC.dsafe, MPC.v_max, MPC.a_max);

    t = 1;
    MPC.dT = dT;
    MPC.X = zeros(1,simTimeMax);
    MPC.Y = zeros(1,simTimeMax);
    MPC.Theta = zeros(1,simTimeMax);
    MPC.eta = zeros(1,simTimeMax);
    MPC.X(1) = x;
    MPC.Y(1) = y;
    MPC.Theta(1) = theta;
    MPC.eta(1) = MPC.eta_min;

    t_max = -inf;
    t_min = inf;
    t_count = 0;
    t_sum = 0;

    Obstacle_count = size(obstacle,1);

    while norm([x_goal y_goal] - [x y]) > position_accuracy && t < simTimeMax  
        tic;
    
        % Calculate Attractive Potential
        if norm([x y]-[x_goal y_goal]) <= MPC.dstar
            nablaU_att =  MPC.zeta*([x y]-[x_goal y_goal]);
        else 
            nablaU_att = MPC.dstar/norm([x y]-[x_goal y_goal]) * MPC.zeta*([x y]-[x_goal y_goal]);
        end
    
        % Find the minimum distance from the obstacle & Calculate Repulsive Potential
        obst_idx = zeros(1,Obstacle_count);
        obst_dist = zeros(1,Obstacle_count);
        nablaU_obst = [0 0];
        for i=1:Obstacle_count
            [obst_idx(i), obst_dist(i)] = dsearchn([obstacle(i,:,1)' obstacle(i,:,2)'], [x y]);
            alpha = theta - atan2(obstacle(i,obst_idx(i),2)-y, obstacle(i,obst_idx(i),1)-x);
            alpha = atan2(sin(alpha), cos(alpha));
            if obst_dist(i) <= MPC.Qstar &&  abs(alpha) < deg2rad(150)
                nablaU_rep_Oi = (MPC.eta(t)*(1/MPC.Qstar - 1/obst_dist(i)) * 1/obst_dist(i)^2)*([x y] - [obstacle(i,obst_idx(i),1) obstacle(i,obst_idx(i),2)]);
                
                if obst_dist(i) <= MPC.dsafe 
                    drel_Oi = 0;
                elseif obst_dist(i) >= 2*MPC.dvort-MPC.dsafe 
                    drel_Oi = 1; 
                else
                    drel_Oi = (obst_dist(i)-MPC.dsafe) / (2*(MPC.dvort-MPC.dsafe));    
                end
                
                if rad2deg(alpha) <= MPC.alpha_th, D_alpha = +1; else, D_alpha = -1; end 

                if drel_Oi <= 0.5
                    gamma = pi * D_alpha * drel_Oi;
                else
                    gamma = pi * D_alpha * (1-drel_Oi);
                end

                R = [cos(gamma) -sin(gamma)
                     sin(gamma)  cos(gamma) ];

                nablaU_obst = nablaU_obst + (R*nablaU_rep_Oi')';
            end
        end

        % Widrow-Hoff rule
        minObstacleDistance = obst_dist(obst_dist < MPC.Qstar);  
        minObstacleDistance = minObstacleDistance(minObstacleDistance > 1e-6);  
        minObstacleDistance = min(minObstacleDistance);   
        try
            e = MPC.dvort - minObstacleDistance;
            MPC.eta(t+1) = MPC.eta(t) + adaptationGain * e * minObstacleDistance;
            MPC.eta(t+1) = max(min(MPC.eta(t+1), MPC.eta_max), MPC.eta_min);
        catch
            MPC.eta(t+1) = MPC.eta(t);
        end

        % Calculate final potential
        nablaU = nablaU_att+nablaU_obst;

        %using mpc
        % Calculate reference value of linear velocity (v_ref) and orientation (theta_ref)
        theta_ref = atan2(-nablaU(2), -nablaU(1));

        error_theta = theta_ref - theta;
        t_i = toc;
        t_max = max(t_max, t_i);
        t_min = min(t_min, t_i);
        t_sum = t_sum + t_i;
        t_count = t_count + 1;
    
        % Simple kinematic mobile robot model
        % Omitted dynamics.
        [v_ref,omega_ref]= solveMPC(x,y,theta, x_goal, y_goal, MPC, N, dT,theta_ref);
        theta = theta + omega_ref * dT;
        theta = atan2(sin(theta), cos(theta));
        x = x + v_ref*cos(theta) * dT;
        y = y + v_ref*sin(theta) * dT;
    
        % Archive and plot it
        t = t + 1;
        MPC.X(t) = x;
        MPC.Y(t) = y;
        MPC.Theta(t) = theta;
    end
    MPC.t = t;
    MPC.travelTime = (t-1)*dT;
    MPC.MeanCalculationTime = t_sum/t_count;
    MPC.MaxCalculationTime = t_max;
    MPC.MinCalculationTime = t_min;
end

function [asapf,nablaU, omega_ref] = AdaptiveSafeArtificialPotentialField(dvort, adaptationGain, x, y,z, theta, x_goal, y_goal, z_goal, position_accuracy, obstacle, dT, N, ASAPF)

    ASAPF.X = x;
    ASAPF.Y = y;
    ASAPF.Z = z;
    ASAPF.Theta = theta;
    Obstacle_count = size(obstacle,1);

    %while norm([x_goal y_goal z_goal] - [x y z]) > position_accuracy && t < simTimeMax  
        tic;

        % Calculate Attractive Potential
        if norm([x y z]-[x_goal y_goal z_goal]) <= ASAPF.dstar
            nablaU_att =  ASAPF.zeta*([x y z]-[x_goal y_goal z_goal]);
        else 
            nablaU_att = ASAPF.dstar/norm([x y z]-[x_goal y_goal z_goal]) * ASAPF.zeta*([x y z]-[x_goal y_goal z_goal]);
        end

        % Find the minimum distance from the obstacle & Calculate Repulsive Potential
        obst_idx = zeros(1,Obstacle_count);
        obst_dist = zeros(1,Obstacle_count);
        nablaU_obst = [0 0 0];
        for i=1:Obstacle_count
            [obst_idx(i), obst_dist(i)] = dsearchn([obstacle(i,:,1)' obstacle(i,:,2)' obstacle(i,:,3)'], [x y z]);
            % Current position in 3D
            pos = [x, y, z]; % Replace x, y, z with actual 3D current position coordinates

            % Obstacle position in 3D
            obst_pos = [obstacle(i, obst_idx(i), 1), obstacle(i, obst_idx(i), 2), obstacle(i, obst_idx(i), 3)];

            % Direction vector from current position to obstacle
            dir_to_obst = obst_pos - pos;

            % Desired direction vector (e.g., towards goal or along current velocity vector)
            dir_motion = [x_goal - x, y_goal - y, z_goal - z]; % Replace with actual direction vector

            % Normalize the vectors
            dir_to_obst = dir_to_obst / norm(dir_to_obst);
            dir_motion = dir_motion / norm(dir_motion);

            % Compute angle between direction vectors in 3D
            alpha = acos(dot(dir_to_obst, dir_motion));

            if obst_dist(i) <= ASAPF.Qstar &&  abs(alpha) < deg2rad(150)
                nablaU_rep_Oi = (ASAPF.eta*(1/ASAPF.Qstar - 1/obst_dist(i)) * 1/obst_dist(i)^2)*([x y z] - [obstacle(i,obst_idx(i),1) obstacle(i,obst_idx(i),2) obstacle(i,obst_idx(i),3)]);

                if obst_dist(i) <= ASAPF.dsafe 
                    drel_Oi = 0;
                elseif obst_dist(i) >= 2*ASAPF.dvort-ASAPF.dsafe 
                    drel_Oi = 1; 
                else
                    drel_Oi = (obst_dist(i)-ASAPF.dsafe) / (2*(ASAPF.dvort-ASAPF.dsafe));    
                end

                if rad2deg(alpha) <= ASAPF.alpha_th, D_alpha = +1; else, D_alpha = -1; end 

                dir_to_obst = obst_pos - pos; % Vector to the obstacle
                dir_motion = [x_goal - x, y_goal - y, z_goal - z]; % Desired direction vector
                rotation_axis = cross(dir_motion, dir_to_obst);
                rotation_axis = rotation_axis / norm(rotation_axis); % Normalize

                % Construct the 3D rotation matrix using Rodrigues' rotation formula
                K = [0, -rotation_axis(3), rotation_axis(2);
                     rotation_axis(3), 0, -rotation_axis(1);
                     -rotation_axis(2), rotation_axis(1), 0];
                % Replace 'gamma' with 'rotation_angle' to avoid MATLAB function conflict
                rotation_angle = pi * D_alpha * drel_Oi;

                if drel_Oi <= 0.5
                    rotation_angle = pi * D_alpha * drel_Oi;
                else
                    rotation_angle = pi * D_alpha * (1 - drel_Oi);
                end


                R = eye(3) + sin(rotation_angle) * K + (1 - cos(rotation_angle)) * (K * K);

                % Apply the rotation matrix to the repulsive force vector
                nablaU_obst = nablaU_obst + (R * nablaU_rep_Oi')';           
            end
        end

        % Widrow-Hoff rule
        minObstacleDistance = obst_dist(obst_dist < ASAPF.Qstar);  
        minObstacleDistance = minObstacleDistance(minObstacleDistance > 1e-6);  
        minObstacleDistance = min(minObstacleDistance);   
        try
            e = ASAPF.dvort - minObstacleDistance;
            asapf.eta = ASAPF.eta + adaptationGain * e * minObstacleDistance;
            asapf.eta = max(min(asapf.eta, ASAPF.eta_max), ASAPF.eta_min);
        catch
            asapf.eta = ASAPF.eta;
        end

        % Calculate final potential
        nablaU = nablaU_att+nablaU_obst;
        theta_ref = atan2(-nablaU(2), -nablaU(1));
        error_theta = theta_ref - theta;
        if abs(error_theta) <= ASAPF.error_theta_max
            alpha = (ASAPF.error_theta_max - abs(error_theta)) / ASAPF.error_theta_max;
            v_ref = min( alpha*norm(-nablaU), ASAPF.v_max );
        else
            v_ref = 0;
        end

        % Simple kinematic mobile robot model
        % Omitted dynamics.
        omega_ref = ASAPF.Kp_omega * error_theta;
        omega_ref = min( max(omega_ref, -ASAPF.omega_max), ASAPF.omega_max);

end


% function [asapf, nablaU, omega_ref] = AdaptiveSafeArtificialPotentialField(dvort, adaptationGain, x, y, z, theta, x_goal, y_goal, z_goal, position_accuracy, obstacle, dT, N, ASAPF)
% 
%     ASAPF.X = x;
%     ASAPF.Y = y;
%     ASAPF.Theta = theta;
% 
%     % Delaunay triangulation for obstacle neighbors
%     obst_points = unique(reshape(obstacle, [], 3), 'rows'); % Flatten obstacle array to points
%     fprintf("No of obstacles= %d\n",size(obst_points,1));
%     delaunayTri = delaunayTriangulation(obst_points); % Create triangulation
%     neighbors = delaunayTri.vertexAttachments; % Neighbors for each point %Returns a cell array where each element contains the indices of tetrahedra (or triangles) attached to the corresponding point in obst_points.
% %From these attachments, you can infer the neighboring points.
% 
%     % Calculate Attractive Potential
%     if norm([x y z] - [x_goal y_goal z_goal]) <= ASAPF.dstar
%         nablaU_att = ASAPF.zeta * ([x y z] - [x_goal y_goal z_goal]);
%     else
%         nablaU_att = (ASAPF.dstar / norm([x y z] - [x_goal y_goal z_goal])) * ASAPF.zeta * ([x y z] - [x_goal y_goal z_goal]);
%     end
% 
%     % Calculate Repulsive Potential using nearest two neighbors
%     nablaU_obst = [0 0 0];
%     %fprintf('Initial nablaU_obst size: [%d, %d]\n', size(nablaU_obst));
%     for i = 1:size(obst_points, 1)
%         % Get neighbors for the current obstacle point
%         neighbor_indices = neighbors{i};
% 
%         % Validate indices (ensure they're within bounds)
%         valid_neighbors = neighbor_indices(neighbor_indices <= size(obst_points, 1));
%         if isempty(valid_neighbors)
%             continue; % Skip if no valid neighbors
%         end
% 
%         % Compute distances to the valid neighbors
%         distances = vecnorm(obst_points(valid_neighbors, :) - [x, y, z], 2, 2);
% 
%         % Find nearest two neighbors
%         [~, sorted_idx] = sort(distances);
%         nearest_neighbors = valid_neighbors(sorted_idx(1:min(2, numel(sorted_idx))));
% 
%         % Calculate repulsive potential
%         % for neighbor_idx = nearest_neighbors'
%         %     neighbor_pos = obst_points(neighbor_idx, :);
%         %     dist_to_neighbor = norm([x, y, z] - neighbor_pos);
%         % 
%         %     if dist_to_neighbor <= ASAPF.Qstar
%         %         scalar_term = ASAPF.eta * (1 / ASAPF.Qstar - 1 / dist_to_neighbor) * (1 / dist_to_neighbor^2)
%         %         fprintf('Size of position vector: [%d, %d]\n', size([x, y, z]));
%         %         fprintf('Size of neighbor_pos: [%d, %d]\n', size(neighbor_pos));
%         % 
%         %         nablaU_rep_Oi = scalar_term .* ([x, y, z] - neighbor_pos);
%         %         fprintf('Size of nablaU_rep_0i: [%d, %d]\n', size(nablaU_rep_Oi));
%         %         %nablaU_rep_Oi = (ASAPF.eta * (1 / ASAPF.Qstar - 1 / dist_to_neighbor) * (1 / dist_to_neighbor^2)) * ([x, y, z] - neighbor_pos);
%         %         nablaU_obst = nablaU_obst + nablaU_rep_Oi;
%         %         fprintf('Size of nablaU_obst(inside loop): [%d, %d]\n', size(nablaU_obst));
%         %     end
%         % end
% 
%         % Calculate repulsive potential
%         for neighbor_idx = nearest_neighbors'
%             neighbor_pos = obst_points(neighbor_idx, :);
% 
%             % Ensure neighbor_pos has consistent dimensions
%             if size(neighbor_pos, 1) == 2  % Case: Two neighbors
%                 for j = 1:size(neighbor_pos, 1)
%                     current_neighbor_pos = neighbor_pos(j, :);
%                     dist_to_neighbor = norm([x, y, z] - current_neighbor_pos);
% 
%                     if dist_to_neighbor <= ASAPF.Qstar
%                         nablaU_rep_Oi = (ASAPF.eta * (1 / ASAPF.Qstar - 1 / dist_to_neighbor) * (1 / dist_to_neighbor^2)) ...
%                                         * ([x, y, z] - current_neighbor_pos);
%                         nablaU_obst = nablaU_obst + nablaU_rep_Oi;
%                     end
%                 end
%             else  % Case: Single neighbor
%                 dist_to_neighbor = norm([x, y, z] - neighbor_pos);
% 
%                 if dist_to_neighbor <= ASAPF.Qstar
%                     nablaU_rep_Oi = (ASAPF.eta * (1 / ASAPF.Qstar - 1 / dist_to_neighbor) * (1 / dist_to_neighbor^2)) ...
%                                     * ([x, y, z] - neighbor_pos);
%                     nablaU_obst = nablaU_obst + nablaU_rep_Oi;
%                 end
%             end
%         end
% 
%     end
% 
%     % Update potential field parameters using Widrow-Hoff rule
%     minObstacleDistance = vecnorm(obst_points - [x y z], 2, 2);
%     minObstacleDistance = minObstacleDistance(minObstacleDistance < ASAPF.Qstar & minObstacleDistance > 1e-6);
%     if ~isempty(minObstacleDistance)
%         minDistance = min(minObstacleDistance);
%         e = ASAPF.dvort - minDistance;
%         asapf.eta = ASAPF.eta + adaptationGain * e * minDistance;
%         asapf.eta = max(min(asapf.eta, ASAPF.eta_max), ASAPF.eta_min);
%     else
%         asapf.eta = ASAPF.eta;
%     end
% 
%     % Calculate final potential
%     nablaU = nablaU_att + nablaU_obst;
%     %fprintf('Size of nablaU_obst: [%d, %d]\n', size(nablaU_obst));
% 
%     theta_ref = atan2(-nablaU(2), -nablaU(1));
%     error_theta = theta_ref - theta;
% 
%     if abs(error_theta) <= ASAPF.error_theta_max
%         alpha = (ASAPF.error_theta_max - abs(error_theta)) / ASAPF.error_theta_max;
%         v_ref = min(alpha * norm(-nablaU), ASAPF.v_max);
%     else
%         v_ref = 0;
%     end
% 
%     % Simple kinematic mobile robot model
%     omega_ref = ASAPF.Kp_omega * error_theta;
%     omega_ref = min(max(omega_ref, -ASAPF.omega_max), ASAPF.omega_max);
% end


% function [asapf, nablaU, omega_ref] = AdaptiveSafeArtificialPotentialField(dvort, adaptationGain, x, y, z, theta, x_goal, y_goal, z_goal, position_accuracy, obstacle, dT, N, ASAPF)
% 
%     ASAPF.X = x;
%     ASAPF.Y = y;
%     ASAPF.Theta = theta;
%     Obstacle_count = size(obstacle, 1);
% 
%     % Flatten obstacle array to 3D points and do Delaunay triangulation
%     obst_points = reshape(obstacle, [], 3);  % Reshape to 2D array of points (each row is an obstacle's [x, y, z])
%     delaunayTri = delaunayTriangulation(obst_points);  % Create triangulation
%     neighbors = delaunayTri.vertexAttachments;  % Get neighbors for each obstacle
% 
%     % Calculate Attractive Potential
%     if norm([x y z] - [x_goal y_goal z_goal]) <= ASAPF.dstar
%         nablaU_att = ASAPF.zeta * ([x y z] - [x_goal y_goal z_goal]);
%     else
%         nablaU_att = (ASAPF.dstar / norm([x y z] - [x_goal y_goal z_goal])) * ASAPF.zeta * ([x y z] - [x_goal y_goal z_goal]);
%     end
% 
%     % Initialize repulsive potential
%     nablaU_obst = [0 0 0];
% 
%     % Loop through each obstacle
%     for i = 1:Obstacle_count
%         % Get neighbors for the current obstacle point
%         neighbor_indices = neighbors{i};
% 
%         % Validate indices (ensure they're within bounds)
%         valid_neighbors = neighbor_indices(neighbor_indices <= size(obst_points, 1));
%         if isempty(valid_neighbors)
%             continue; % Skip if no valid neighbors
%         end
% 
%         % Compute distances to the valid neighbors
%         distances = vecnorm(obst_points(valid_neighbors, :) - [x, y, z], 2, 2);
% 
%         % Find nearest two neighbors
%         [~, sorted_idx] = sort(distances);
%         nearest_neighbors = valid_neighbors(sorted_idx(1:min(2, numel(sorted_idx))));
% 
%         % Calculate repulsive potential only from the nearest two obstacles
%         for neighbor_idx = nearest_neighbors'
%             neighbor_pos = obst_points(neighbor_idx, :);
%             dist_to_neighbor = norm([x, y, z] - neighbor_pos);
% 
%             if dist_to_neighbor <= ASAPF.Qstar
%                 nablaU_rep_Oi = (ASAPF.eta * (1 / ASAPF.Qstar - 1 / dist_to_neighbor) * (1 / dist_to_neighbor^2)) * ([x, y, z] - neighbor_pos);
%                 nablaU_obst = nablaU_obst + nablaU_rep_Oi;
%             end
%         end
%     end
% 
%     % Widrow-Hoff rule
%     minObstacleDistance = vecnorm(obst_points - [x, y, z], 2, 2);
%     minObstacleDistance = minObstacleDistance(minObstacleDistance < ASAPF.Qstar & minObstacleDistance > 1e-6);
%     if ~isempty(minObstacleDistance)
%         minDistance = min(minObstacleDistance);
%         e = ASAPF.dvort - minDistance;
%         asapf.eta = ASAPF.eta + adaptationGain * e * minDistance;
%         asapf.eta = max(min(asapf.eta, ASAPF.eta_max), ASAPF.eta_min);
%     else
%         asapf.eta = ASAPF.eta;
%     end
% 
%     % Calculate final potential
%     nablaU = nablaU_att + nablaU_obst;
%     fprintf('Size of nablaU: [%d, %d]\n', size(nablaU));
%     fprintf('Size of nablaU_att: [%d, %d]\n', size(nablaU_att));
%     fprintf('Size of nablaU_obst: [%d, %d]\n', size(nablaU_obst));
% 
%     theta_ref = atan2(-nablaU(2), -nablaU(1));
%     error_theta = theta_ref - theta;
% 
%     if abs(error_theta) <= ASAPF.error_theta_max
%         alpha = (ASAPF.error_theta_max - abs(error_theta)) / ASAPF.error_theta_max;
%         v_ref = min(alpha * norm(-nablaU), ASAPF.v_max);
%     else
%         v_ref = 0;
%     end
% 
%     % Simple kinematic mobile robot model
%     omega_ref = ASAPF.Kp_omega * error_theta;
%     omega_ref = min(max(omega_ref, -ASAPF.omega_max), ASAPF.omega_max);
% end
% 



function [v_ref, omega_ref]=solveMPC(x, y, theta, x_goal, y_goal,MPC, N, dT, theta_ref)
num_vars=2*N ; %[v, omega] for each timestep
x0= zeros(num_vars, 1);

%bounds for control inputs
lb=repmat([-MPC.v_max; -MPC.omega_max],N, 1);
ub=repmat([MPC.v_max; MPC.omega_max], N,1);

%optimization options
options= optimoptions('fmincon', 'Display','off', 'Algorithm','sqp');

[u_opt, ~]=fmincon(@(u)objectiveFunction(u,x,y,theta,x_goal,y_goal, MPC, N, dT, theta_ref), x0, [], [], [], [], lb, ub, [], options);

v_ref=u_opt(1);
omega_ref=u_opt(2);

end

function J= objectiveFunction(u,x0,y0,theta0, x_goal, y_goal,MPC, N, dT, theta_ref)
J=0;
x=x0;
y=y0;
theta=theta0;

disp(size(theta_ref));

for k=1:2:2*N
    v=u(k);
    omega=u(k+1);
    error_theta=theta_ref-theta;
    %error_theta= atan2(cos(error_theta), sin(error_theta));
    theta_next= theta+ omega*dT;
    x_next= x+ v*cos(theta)*dT;
    y_next= y+ v*sin(theta)*dT;

    %position error term in obj function
    pos_err= [x_next-x_goal; y_next-y_goal];
    J=J+ (pos_err' * MPC.Q_pos * pos_err);

    %we want to maximise potential, or minimise (-potentai)
    %J=J+ lambda_rep*((nablau_obst-nablau_att)*ASAPF.Q_obs*(nablau_obst-nablau_att)');
    J=J+ error_theta'*MPC.Q_theta*error_theta;


    %control input
    J= J+ ([v;omega]' * MPC.R_input * [v;omega]);

    x=x_next;
    y=y_next;
    theta=theta_next;
end

end

