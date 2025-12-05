function [u_opt, x_pred] = quad_mpc(current_state, reference_trajectory, N_horizon)
    import casadi.*
    
    param; 
    
    dt = 0.1;  % Sampling time
    
    % Define state and control dimensions
    nx = 12;  % State dimension
    nu = 4;   % Control inputs: [thrust, tau_x, tau_y, tau_z]
    
    % Define optimization variables
    x = SX.sym('x', nx, N_horizon+1);  % States over horizon
    u = SX.sym('u', nu, N_horizon);    % Controls over horizon
    
    % Define the MPC objective function weights (as in the paper)
    Q = diag([200, 200, 2000, 10, 10, 10, 1, 1, 1, 0.2, 0.2, 0.2]);  % State weights
    R = diag([0.5, 500, 500, 500]);  % Control weights
    
    % Define state bounds
    x_min = [-pi/4; -pi/4; -inf; 0; 0; -10; -inf; -inf; -inf; -1; -1; -1];
    x_max = [pi/4; pi/4; inf; 300; 300; 10; inf; inf; inf; 8; 8; 8];
    
    % Define control bounds
    u_min = [0; -2; -2; -2];     % [thrust, tau_x, tau_y, tau_z]
    u_max = [40; 2; 2; 2];
    
    % Initialize objective function and constraints
    obj = 0;
    g = [];  % Constraints
    
    % Initial state constraint
    x_current = current_state;
    g = [g; x(:,1) - x_current'];
    
    % Loop over control horizon
    for k = 1:N_horizon
        % Get reference for this step
        x_ref = reference_trajectory(:, min(k, size(reference_trajectory, 2)));
        u_ref = [m*G; 0; 0; 0]; % Hovering reference control
        
        % Add stage cost: weighted error from reference
        obj = obj + (x(:,k) - x_ref)'*Q*(x(:,k) - x_ref) + ...
                    (u(:,k) - u_ref)'*R*(u(:,k) - u_ref);
        
        % State equation constraint (system dynamics)
        f = quad_dynamics_casadi(x(:,k), u(2,k), u(3,k), u(4,k), u(1,k));
        x_next = x(:,k) + dt*f;
        g = [g; x(:,k+1) - x_next];
    end
    
    % Terminal cost
    x_ref_final = reference_trajectory(:, min(N_horizon+1, size(reference_trajectory, 2)));
    obj = obj + (x(:,N_horizon+1) - x_ref_final)'*Q*(x(:,N_horizon+1) - x_ref_final);
    
    % Define NLP
    opts = struct;
    opts.ipopt.max_iter = 100;
    opts.ipopt.print_level = 0;
    opts.print_time = 0;
    opts.ipopt.acceptable_tol = 1e-8;
    opts.ipopt.acceptable_obj_change_tol = 1e-6;
    
    nlp = struct('f', obj, 'x', [reshape(x, (N_horizon+1)*nx, 1); reshape(u, N_horizon*nu, 1)], ...
                 'g', g);
    solver = nlpsol('solver', 'ipopt', nlp, opts);
    
    % Set bounds for variables and constraints
    x_bounds_lower = repmat(x_min, N_horizon+1, 1);
    x_bounds_upper = repmat(x_max, N_horizon+1, 1);
    u_bounds_lower = repmat(u_min, N_horizon, 1);
    u_bounds_upper = repmat(u_max, N_horizon, 1);
    
    vars_lb = [x_bounds_lower; u_bounds_lower];
    vars_ub = [x_bounds_upper; u_bounds_upper];
    
    % Equality constraints (system dynamics and initial state)
    constraints_lb = zeros(size(g));
    constraints_ub = zeros(size(g));
    
    % Initial guess for the optimization problem
    vars_init = zeros((N_horizon+1)*nx + N_horizon*nu, 1);
    
    % Solve the NLP
    sol = solver('x0', vars_init, 'lbx', vars_lb, 'ubx', vars_ub, ...
                  'lbg', constraints_lb, 'ubg', constraints_ub);
    
    % Extract the solution
    solution = full(sol.x);
    
    % Reshape the solution to get states and controls
    x_opt = reshape(solution(1:(N_horizon+1)*nx), nx, N_horizon+1);
    u_opt = reshape(solution((N_horizon+1)*nx+1:end), nu, N_horizon);
    
    % Return the optimal control sequence and predicted states
    x_pred = x_opt;
end

