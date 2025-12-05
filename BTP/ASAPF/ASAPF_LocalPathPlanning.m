function [simout, simout2] = ASAPF_LocalPathPlanning()
   
    % Initial position and orientation 
    x = -0.5;
    y = 0.5;
    theta = 0;
    
    % Goal position
    x_goal = 3.5;
    y_goal = 2.75;
    position_accuracy = 0.1;
    
    % Sampling period
    dT = 0.1;
    N=10;
    
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

    % Simulation
    simTimeMax = 600; 

    simout = AdaptiveSafeArtificialPotentialField(0.4, 4e0 *dT, x, y, theta, x_goal, y_goal, position_accuracy, obstacle, dT, simTimeMax, N);
    simout2 = ASAPF_using_MPC(0.4, 4e0 *dT, x, y, theta, x_goal, y_goal, position_accuracy, obstacle, dT, simTimeMax, N);

    % Plot it
    figure(1); subplot(1,1,1);
    cla; hold on; grid on; box on;
    daspect([1 1 1]); 
    xlim([-1,4]);  ylim([-1 3]); xlabel("x [m]"); ylabel("y [m]");
    box on; hold on;
    plot(simout.X(1:simout.t), simout.Y(1:simout.t), 'Color',[0 0.4470 0.7410], 'LineWidth', 1); % Plot traveled path
    plot(x_goal, y_goal, 'xg');
    for i=1:Obstacle_count
        plot(obstacle(i,:,1), obstacle(i,:,2), '-r');
        plot(obstacle_02(i,:,1), obstacle_02(i,:,2), '--r');
        plot(obstacle_04(i,:,1), obstacle_04(i,:,2), '--b');
    end
    drawnow;
    title('ASAPF without MPC')
       
    figure(2);
    sgtitle('Plots without using MPC')
    subplot(311); cla; hold on; grid on; box on; ylabel("v [m/s]")
    subplot(312); cla; hold on; grid on; box on; ylabel("\omega [rad/s]")
    subplot(313); cla; hold on; grid on; box on; ylabel("\eta"); xlabel("t [s]")
    [t, v, w] = CalculateVelocities(simout);
    subplot(311);  xlim([0 max(t)]); ylim([0 0.22])
    plot(t,v, 'Color',[0 0.4470 0.7410]);
    subplot(312);  xlim([0 max(t)]); 
    plot(t,w, 'Color',[0 0.4470 0.7410]);
    subplot(313); xlim([0 max(t)]); 
    plot(dT*(1:simout.t),simout.eta(1:simout.t), 'Color', [0 0.4470 0.7410]);
    

    figure(3); subplot(1,1,1);
    cla; hold on; grid on; box on;
    daspect([1 1 1]); 
    xlim([-1,4]);  ylim([-1 3]); xlabel("x [m]"); ylabel("y [m]");
    box on; hold on;
    plot(simout2.X(1:simout2.t), simout2.Y(1:simout2.t), 'Color',[0 0.4470 0.7410], 'LineWidth', 1); % Plot traveled path
    plot(x_goal, y_goal, 'xg');
    for i=1:Obstacle_count
        plot(obstacle(i,:,1), obstacle(i,:,2), '-r');
        plot(obstacle_02(i,:,1), obstacle_02(i,:,2), '--r');
        plot(obstacle_04(i,:,1), obstacle_04(i,:,2), '--b');
    end
    drawnow;
    title('ASAPF with MPC')

    figure(4);
    subplot(311); cla; hold on; grid on; box on; ylabel("v [m/s]")
    subplot(312); cla; hold on; grid on; box on; ylabel("\omega [rad/s]")
    subplot(313); cla; hold on; grid on; box on; ylabel("\eta"); xlabel("t [s]")
    [t, v, w] = CalculateVelocities(simout2);
    subplot(311);  xlim([0 max(t)]); ylim([0 0.22])
    plot(t,v, 'Color',[0 0.4470 0.7410]);
    subplot(312);  xlim([0 max(t)]); 
    plot(t,w, 'Color',[0 0.4470 0.7410]);
    subplot(313); xlim([0 max(t)]); 
    plot(dT*(1:simout2.t),simout2.eta(1:simout2.t), 'Color', [0 0.4470 0.7410]);
    sgtitle('Plots with mpc');


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

function ASAPF = AdaptiveSafeArtificialPotentialField(dvort, adaptationGain, x, y, theta, x_goal, y_goal, position_accuracy, obstacle, dT, simTimeMax, N)
   
    % SAPF parameters
    ASAPF.dstar = 0.3;
    ASAPF.Qstar = 1.5;
    ASAPF.dsafe = 0.2;
    ASAPF.dvort = dvort;
    ASAPF.alpha_th = deg2rad(5);

    % Parameters related to kinematic model
    ASAPF.error_theta_max = deg2rad(45);
    ASAPF.v_max = 0.2;
    ASAPF.a_max = 0.2;
    ASAPF.Kp_omega = 1.5;
    ASAPF.omega_max = 0.5*pi; 

    [ASAPF.zeta, ASAPF.eta_max, ASAPF.eta_min] = ASAPF_CalculateParameters(ASAPF.Qstar, ASAPF.dstar, ASAPF.dsafe, ASAPF.v_max, ASAPF.a_max);

    t = 1;
    ASAPF.dT = dT;
    ASAPF.X = zeros(1,simTimeMax);
    ASAPF.Y = zeros(1,simTimeMax);
    ASAPF.Theta = zeros(1,simTimeMax);
    ASAPF.eta = zeros(1,simTimeMax);
    ASAPF.X(1) = x;
    ASAPF.Y(1) = y;
    ASAPF.Theta(1) = theta;
    ASAPF.eta(1) = ASAPF.eta_min;

    t_max = -inf;
    t_min = inf;
    t_count = 0;
    t_sum = 0;

    Obstacle_count = size(obstacle,1);

    while norm([x_goal y_goal] - [x y]) > position_accuracy && t < simTimeMax  
        tic;
    
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
            alpha = theta - atan2(obstacle(i,obst_idx(i),2)-y, obstacle(i,obst_idx(i),1)-x);
            alpha = atan2(sin(alpha), cos(alpha));
            if obst_dist(i) <= ASAPF.Qstar &&  abs(alpha) < deg2rad(150)
                nablaU_rep_Oi = (ASAPF.eta(t)*(1/ASAPF.Qstar - 1/obst_dist(i)) * 1/obst_dist(i)^2)*([x y] - [obstacle(i,obst_idx(i),1) obstacle(i,obst_idx(i),2)]);
                
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
            ASAPF.eta(t+1) = ASAPF.eta(t) + adaptationGain * e * minObstacleDistance;
            ASAPF.eta(t+1) = max(min(ASAPF.eta(t+1), ASAPF.eta_max), ASAPF.eta_min);
        catch
            ASAPF.eta(t+1) = ASAPF.eta(t);
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

        t_i = toc;
        t_max = max(t_max, t_i);
        t_min = min(t_min, t_i);
        t_sum = t_sum + t_i;
        t_count = t_count + 1;
    
        % Simple kinematic mobile robot model
        % Omitted dynamics.
        omega_ref = ASAPF.Kp_omega * error_theta;
        omega_ref = min( max(omega_ref, -ASAPF.omega_max), ASAPF.omega_max);
        theta = theta + omega_ref * dT;
        theta = atan2(sin(theta), cos(theta));
        x = x + v_ref*cos(theta) * dT;
        y = y + v_ref*sin(theta) * dT;
    
        % Archive and plot it
        t = t + 1;
        ASAPF.X(t) = x;
        ASAPF.Y(t) = y;
        ASAPF.Theta(t) = theta;
    end

    ASAPF.t = t;
    ASAPF.travelTime = (t-1)*dT;
    ASAPF.MeanCalculationTime = t_sum/t_count;
    ASAPF.MaxCalculationTime = t_max;
    ASAPF.MinCalculationTime = t_min;
end

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

