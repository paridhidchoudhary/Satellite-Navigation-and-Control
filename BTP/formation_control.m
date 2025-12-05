function [phi_f] = formation_control(eta, xi, desired_formation)
    % eta: 6xN matrix of relative configurations for 3 satellites
    % xi: 6xN matrix of relative velocities for 3 satellites
    % desired_formation: 3x3 matrix of desired relative positions

    current_positions = eta(4:6, :);
    relative_positions = current_positions - current_positions(:, 1);
    
    errors = relative_positions - desired_formation;
    
    % Simple P controller for formation
    K_f = 0.1;  % Formation control gain
    phi_f = -K_f * errors;
end