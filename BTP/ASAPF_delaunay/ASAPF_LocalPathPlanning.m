
% 
% function [ASAPF, nablaU, v_x, v_y, psi_ref] = ASAPF_LocalPathPlanning(dvort, adaptationGain, x, y, psi, x_goal, y_goal, position_accuracy, obstacle,r,c)
%     %disp('into asapf fn \n');
%     % SAPF parameters
%     ASAPF.dstar = 0.3;
%     ASAPF.Qstar = 1.5;
%     ASAPF.dsafe = 0.2;
%     ASAPF.dvort = dvort;
%     ASAPF.alpha_th = deg2rad(5);
% 
%     % Parameters related to kinematic model
%     ASAPF.error_psi_max = deg2rad(45);
%     ASAPF.v_max = 0.2;
%     ASAPF.a_max = 0.2;
%     ASAPF.Kp_omega = 1.5;
%     ASAPF.omega_max = 0.5*pi; 
% 
%     [ASAPF.zeta, ASAPF.eta_max, ASAPF.eta_min] = ASAPF_CalculateParameters(ASAPF.Qstar, ASAPF.dstar, ASAPF.dsafe, ASAPF.v_max, ASAPF.a_max);
% 
%     %ASAPF.dT = dT;
%     % ASAPF.X = zeros(1,simTimeMax);
%     % ASAPF.Y = zeros(1,simTimeMax);
%     % ASAPF.Theta = zeros(1,simTimeMax);
%     %ASAPF.eta = zeros(1,simTimeMax);
%     ASAPF.X = x;
%     ASAPF.Y = y;
%     ASAPF.Theta = psi;
%     ASAPF.eta = ASAPF.eta_min;
%     %ASAPF.eta
% 
%     % t_max = -inf;
%     % t_min = inf;
%     % t_count = 0;
%     % t_sum = 0;
% 
%     Obstacle_count = size(obstacle,1);
%     delaunayTri = delaunayTriangulation(c(:,1), c(:,2));
%     neighbors = delaunayTri.vertexAttachments;
%     %tic;
%     %disp('into while loop \n');
%     % Calculate Attractive Potential
%     if norm([x y]-[x_goal y_goal]) <= ASAPF.dstar
%         nablaU_att =  ASAPF.zeta*([x y]-[x_goal y_goal]);
%     else 
%         nablaU_att = ASAPF.dstar/norm([x y]-[x_goal y_goal]) * ASAPF.zeta*([x y]-[x_goal y_goal]);
%     end
% 
%     nablaU_obst = [0 0];
%     for i=1:Obstacle_count
% 
%         neighbor_indices = neighbors{i};
% 
%         valid_neighbors = neighbor_indices(neighbor_indices<=Obstacle_count);
%         if isempty(valid_neighbors)
%             continue;
%         end
% 
%         distances = vecnorm(c(valid_neighbors, :) - [x y], 2 , 2);
%         [~, sorted_idx] = sort(distances);
%         nearest_neighbors = valid_neighbors(sorted_idx(1:min(2, numel(sorted_idx))));
% 
%         for neighbor_idx = nearest_neighbors(:)'
%             neighbor_pos = c(neighbor_idx, :);
%             radius_obst = r(neighbor_idx);
%             dist_to_neighbor = norm([x, y] - neighbor_pos);
% 
%             if dist_to_neighbor <= radius_obst
%                 nablaU_rep_Oi = -ASAPF.eta*(1/radius_obst - 1/ASAPF.Qstar)*(1/radius_obst^2)*(neighbor_pos - [x y])/dist_to_neighbor;
%             elseif radius_obst< dist_to_neighbor < ASAPF.Qstar
%                 nablaU_rep_Oi = -ASAPF.eta*(1/dist_to_neighbor - 1/ASAPF.Qstar) *(neighbor_pos - [x y])/(dist_to_neighbor^3);
%             else
%                 nablaU_rep_Oi = [0 0];
%             end
%             nablaU_obst= nablaU_obst + nablaU_rep_Oi;
%         end
%         % 
%         % [obst_idx(i), obst_dist(i)] = dsearchn([obstacle(i,:,1)' obstacle(i,:,2)'], [x y]);
%         % alpha = theta - atan2(obstacle(i,obst_idx(i),2)-y, obstacle(i,obst_idx(i),1)-x);
%         % alpha = atan2(sin(alpha), cos(alpha));
%         % if obst_dist(i) <= ASAPF.Qstar &&  abs(alpha) < deg2rad(150)
%         %     nablaU_rep_Oi = (ASAPF.eta(t)*(1/ASAPF.Qstar - 1/obst_dist(i)) * 1/obst_dist(i)^2)*([x y] - [obstacle(i,obst_idx(i),1) obstacle(i,obst_idx(i),2)]);
%         % 
%         %     if obst_dist(i) <= ASAPF.dsafe 
%         %         drel_Oi = 0;
%         %     elseif obst_dist(i) >= 2*ASAPF.dvort-ASAPF.dsafe 
%         %         drel_Oi = 1; 
%         %     else
%         %         drel_Oi = (obst_dist(i)-ASAPF.dsafe) / (2*(ASAPF.dvort-ASAPF.dsafe));    
%         %     end
%         % 
%         %     if rad2deg(alpha) <= ASAPF.alpha_th, D_alpha = +1; else, D_alpha = -1; end 
%         % 
%         %     if drel_Oi <= 0.5
%         %         gamma = pi * D_alpha * drel_Oi;
%         %     else
%         %         gamma = pi * D_alpha * (1-drel_Oi);
%         %     end
%         % 
%         %     R = [cos(gamma) -sin(gamma)
%         %          sin(gamma)  cos(gamma) ];
%         % 
%         %     nablaU_obst = nablaU_obst + (R*nablaU_rep_Oi')';
%         % end
%     end
% 
%     % Widrow-Hoff rule
%     minObstacleDistance = distances(distances < ASAPF.Qstar);  
%     minObstacleDistance = minObstacleDistance(minObstacleDistance > 1e-6);  
%     minObstacleDistance = min(minObstacleDistance);   
%     try
%         e = ASAPF.dvort - minObstacleDistance;
%         ASAPF.eta = ASAPF.eta + adaptationGain * e * minObstacleDistance;
%         ASAPF.eta = max(min(ASAPF.eta, ASAPF.eta_max), ASAPF.eta_min);
%     catch
%         ASAPF.eta = ASAPF.eta;
%     end
%     ASAPF.eta;
% 
%     % Calculate final potential
%     nablaU = nablaU_att+nablaU_obst;
%     psi_ref = atan2(-nablaU(2), -nablaU(1));
%     % v_ref_vec = -nablaU;  % already pointing in the descent direction
%     % if norm(v_ref_vec) > 1e-3  % to avoid division by zero
%     %     v_ref_vec = v_ref_vec / norm(v_ref_vec);  % normalize direction
%     % end
% 
%     error_psi = psi_ref - psi;
%     %if abs(error_psi) <= ASAPF.error_psi_max
%         alpha = (ASAPF.error_psi_max - abs(error_psi)) / ASAPF.error_psi_max;
%         v_ref = min( alpha*norm(-nablaU), ASAPF.v_max );
%     %else
%     %    v_ref = 0;
%     %end
%     v_x = v_ref * cos(psi_ref);
%     v_y = v_ref * sin(psi_ref);
% 
% end


function [ASAPF, nablaU, v_x, v_y, psi_ref] = ASAPF_LocalPathPlanning(dvort, adaptationGain, x, y, psi, x_goal, y_goal, position_accuracy, obstacle,r,c)
    %disp('into asapf fn \n');
    % SAPF parameters
    ASAPF.dstar = 0.3;
    ASAPF.Qstar = 1.5;
    ASAPF.dsafe = 0.2;
    ASAPF.dvort = dvort;
    ASAPF.alpha_th = deg2rad(5);

    % Parameters related to kinematic model
    ASAPF.error_psi_max = deg2rad(45);
    ASAPF.v_max = 0.2;
    ASAPF.a_max = 0.2;
    ASAPF.Kp_omega = 1.5;
    ASAPF.omega_max = 0.5*pi; 

    [ASAPF.zeta, ASAPF.eta_max, ASAPF.eta_min] = ASAPF_CalculateParameters(ASAPF.Qstar, ASAPF.dstar, ASAPF.dsafe, ASAPF.v_max, ASAPF.a_max);

    %ASAPF.dT = dT;
    % ASAPF.X = zeros(1,simTimeMax);
    % ASAPF.Y = zeros(1,simTimeMax);
    % ASAPF.Theta = zeros(1,simTimeMax);
    %ASAPF.eta = zeros(1,simTimeMax);
    ASAPF.X = x;
    ASAPF.Y = y;
    ASAPF.Theta = psi;
    ASAPF.eta = ASAPF.eta_min;
    %ASAPF.eta

    % t_max = -inf;
    % t_min = inf;
    % t_count = 0;
    % t_sum = 0;

    Obstacle_count = size(obstacle,1);
    delaunayTri = delaunayTriangulation(c(:,1), c(:,2));
    neighbors = delaunayTri.vertexAttachments;
    %tic;
    %disp('into while loop \n');
    % Calculate Attractive Potential
    if norm([x y]-[x_goal y_goal]) <= ASAPF.dstar
        nablaU_att =  ASAPF.zeta*([x y]-[x_goal y_goal]);
    else 
        nablaU_att = ASAPF.dstar/norm([x y]-[x_goal y_goal]) * ASAPF.zeta*([x y]-[x_goal y_goal]);
    end

    % Find the minimum distance from the obstacle & Calculate Repulsive Potential
    obst_idx = zeros(1,Obstacle_count);
    obst_dist = zeros(1,Obstacle_count);
    nablaU_obst = [0 0];
    for i=1:Obstacle_count
        [obst_idx(i), obst_dist(i)] = dsearchn([obstacle(i,:,1)' obstacle(i,:,2)'], [x y]);
        alpha = psi - atan2(obstacle(i,obst_idx(i),2)-y, obstacle(i,obst_idx(i),1)-x);
        alpha = atan2(sin(alpha), cos(alpha));
        if obst_dist(i) <= ASAPF.Qstar &&  abs(alpha) < deg2rad(150)
            nablaU_rep_Oi = (ASAPF.eta*(1/ASAPF.Qstar - 1/obst_dist(i)) * 1/obst_dist(i)^2)*([x y] - [obstacle(i,obst_idx(i),1) obstacle(i,obst_idx(i),2)]);
            
            if obst_dist(i) <= ASAPF.dsafe 
                drel_Oi = 0;
            elseif obst_dist(i) >= 2*ASAPF.dvort-ASAPF.dsafe 
                drel_Oi = 1; 
            else
                drel_Oi = (obst_dist(i)-ASAPF.dsafe) / (2*(ASAPF.dvort-ASAPF.dsafe));    
            end
            
            if rad2deg(alpha) <= ASAPF.alpha_th, D_alpha = +1; else, D_alpha = -1; end 

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
    minObstacleDistance = obst_dist(obst_dist < ASAPF.Qstar);  
    minObstacleDistance = minObstacleDistance(minObstacleDistance > 1e-6);  
    minObstacleDistance = min(minObstacleDistance);   
    try
        e = ASAPF.dvort - minObstacleDistance;
        ASAPF.eta = ASAPF.eta + adaptationGain * e * minObstacleDistance;
        ASAPF.eta = max(min(ASAPF.eta, ASAPF.eta_max), ASAPF.eta_min);
    catch
        ASAPF.eta = ASAPF.eta;
    end
    ASAPF.eta;

    % Calculate final potential
    nablaU = nablaU_att+nablaU_obst;
    psi_ref = atan2(-nablaU(2), -nablaU(1));
    % v_ref_vec = -nablaU;  % already pointing in the descent direction
    % if norm(v_ref_vec) > 1e-3  % to avoid division by zero
    %     v_ref_vec = v_ref_vec / norm(v_ref_vec);  % normalize direction
    % end

    error_psi = psi_ref - psi;
    %if abs(error_psi) <= ASAPF.error_psi_max
        alpha = (ASAPF.error_psi_max - abs(error_psi)) / ASAPF.error_psi_max;
        v_ref = min( alpha*norm(-nablaU), ASAPF.v_max );
    %else
    %    v_ref = 0;
    %end
    v_x = v_ref * cos(psi_ref);
    v_y = v_ref * sin(psi_ref);

end


