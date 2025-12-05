

function [simout, nablaU, Moment, omega_not_dot] = ASAPF_LocalPathPlanning(x,y,theta,z,ve, ve_dot, Rr, omegae_dot, theta_e, omega_not_prev, ASAPF, obstacle)
    
    % Goal position
    x_goal = 0;
    y_goal = 0;
    z_goal= 0;
    position_accuracy = 0.1;
    
    % Sampling period
    dT = 0.1;
    N=10;
    
    % Generate obstacles 

    [simout, nablaU, Moment, omega_not_dot] = AdaptiveSafeArtificialPotentialField(0.4, 4e0 *dT, x, y,z,ve, ve_dot, Rr, omegae_dot, theta_e, omega_not_prev, theta, x_goal, y_goal, z_goal, position_accuracy, obstacle, dT, N, ASAPF);
    %simout2 = ASAPF_using_MPC(0.4, 4e0 *dT, x, y, theta, x_goal, y_goal, position_accuracy, obstacle, dT, simTimeMax, N);

    
end



function [asapf,nablaU, Moment, omega_not_dot_values] = AdaptiveSafeArtificialPotentialField(dvort, adaptationGain, x, y,z,ve,ve_dot, Rr, omegae_dot, theta_e, omega_not_prev, theta, x_goal, y_goal, z_goal, position_accuracy, obstacle, dT, N, ASAPF)
    psi = 1e-6;
    ASAPF.X = x;
    ASAPF.Y = y;
    ASAPF.Z = z;
    ASAPF.Theta = theta;
    Obstacle_count = size(obstacle,1);

    %while norm([x_goal y_goal z_goal] - [x y z]) > position_accuracy && t < simTimeMax  
        tic;
        % disp(['size of omegae_dot = ', num2str(size(omegae_dot))]);
        % disp(['size of theta_e = ', num2str(size(theta_e))]);
        theta_e= theta_e';
        % Calculate Attractive Potential
        %disp(['size of ve_Dot= ', num2str(size(ve_dot))]);
        nablaU_att =  ASAPF.zeta*(ve_dot' + [x y z]);
        %disp(['size of nablau_att= ', num2str(size(nablaU_att))]);
        Ma= ASAPF.zeta*[omegae_dot(1) + theta_e(1), omegae_dot(2) + theta_e(2), omegae_dot(3) + theta_e(3)];
        %disp(['size of Ma = ', num2str(size(Ma))]);

        % Find the minimum distance from the obstacle & Calculate Repulsive Potential
        obst_idx = zeros(1,Obstacle_count);
        obst_dist = zeros(1,Obstacle_count);
        nablaU_obst = [0 0 0];
        Mr= [0 0 0];
        omega_not_dot_values= zeros(3*Obstacle_count, 1);
        k=1;
        for i=1:Obstacle_count
            [obst_idx(i), obst_dist(i)] = dsearchn([obstacle(i,:,1)' obstacle(i,:,2)' obstacle(i,:,3)'], [x y z]);
            % Current position in 3D
            pos = [x, y, z]; % Replace x, y, z with actual 3D current position coordinates

            % Obstacle position in 3D
            obst_pos = [obstacle(i, obst_idx(i), 1), obstacle(i, obst_idx(i), 2), obstacle(i, obst_idx(i), 3)];

            % Direction vector from current position to obstacle
            dir_to_obst = obst_pos - pos;
            xeo= dir_to_obst;
            % size(xeo)
            % size(Rr)
            veo = ve; % obstacle velocity =0
            veo_dot= ve_dot;

            xoe= (- Rr'*xeo')';
            xoe_cap= xoe/ norm(xoe);
            theta_0= [acos(xoe_cap(1)), acos(xoe_cap(2)), acos(xoe_cap(3))];

            voe= -ve;
            % disp(['size of voe= ',num2str(size(voe))]);
            % disp(['size of xoe_cap= ', num2str(size(xoe_cap))]);
            omega_0= [-voe(1)/(sqrt(1- xoe_cap(1)^2)), -voe(2)/sqrt(1- xoe_cap(2)^2), -voe(3)/sqrt(1-xoe_cap(3)^2)];
            % omega_not_values(i) = omega_0(1);
            % omega_not_values(i+1) = omega_0(2);
            % omega_not_values(i+2) = omega_0(3);

            omega_0_dot= [(omega_0(1)-omega_not_prev(k))/dT, (omega_0(2)-omega_not_prev(k+1))/dT, (omega_0(3)-omega_not_prev(k+2))/dT];
            omega_not_dot_values(k,1)= omega_0_dot(1);
            omega_not_dot_values(k+1,1)= omega_0_dot(2);
            omega_not_dot_values(k+2,1)= omega_0_dot(3);
            k= k+3;

            nablaU_rep = (-ASAPF.eta)*([veo_dot(1)/(veo(1)^2 + psi), veo_dot(2)/(veo(2)^2+psi), veo_dot(3)/(veo(3)^2+psi)] + [1/(xeo(1)^3+psi), 1/(xeo(2)^3+psi), 1/(xeo(3)^3+psi)]);
            % disp(['veo_Dot= ', num2str(veo_dot(:).')]);
            % disp(['veo= ', num2str(veo(:).')]);
            % disp(['xeo =', num2str(xeo(:).')]);
            %disp(['size of nablaU_rep = ', num2str(size(nablaU_rep))]);
            nablaU_obst = nablaU_obst + nablaU_rep;
            %disp(['size of nablaU_obst = ', num2str(size(nablaU_obst))]);
            
            Mrep= -ASAPF.eta*([omega_0_dot(1)/(omega_0(1)^2 + psi), omega_0_dot(2)/(omega_0(2)^2 + psi), omega_0_dot(3)/(omega_0(3)^2 + psi)]- [1/(theta_0(1)^3 + psi), 1/(theta_0(2)^3+ psi), 1/(theta_0(3)^3+ psi)]);
            %disp(['size of mrep =', num2str(size(Mrep))]);
            Mr= Mr + Mrep;
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
        % theta_ref = atan2(-nablaU(2), -nablaU(1));
        % error_theta = theta_ref - theta;
        % if abs(error_theta) <= ASAPF.error_theta_max
        %     alpha = (ASAPF.error_theta_max - abs(error_theta)) / ASAPF.error_theta_max;
        %     v_ref = min( alpha*norm(-nablaU), ASAPF.v_max );
        % else
        %     v_ref = 0;
        % end

        % Simple kinematic mobile robot model
        % Omitted dynamics.
        % omega_ref = ASAPF.Kp_omega * error_theta;
        % omega_ref = min( max(omega_ref, -ASAPF.omega_max), ASAPF.omega_max);

        %moment
        Moment = Ma + Mr;

end




