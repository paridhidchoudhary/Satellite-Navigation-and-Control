function u_opt = solve_nmpc_step(xk, omega_cce, omega_tte, thruster_positions, params)
    import casadi.*

    % Create new Opti stack each time
    opti = Opti();

    % Constants
    Jc = params.Jc; Jt = params.Jt; m_c = params.m_c; fmax = params.fmax;
    r_pc_c = params.r_pc_c;

    % Parameters
    omega_cce_param = opti.parameter(3,1);
    omega_tte_param = opti.parameter(3,1);
    thruster_param = opti.parameter(3,8);
    x0_param = opti.parameter(21,1);

    % Variables
    N = params.N; Ts = params.Ts;
    X = opti.variable(21, N+1);
    U = opti.variable(24, N);
    opti.set_value(x0_param, xk);
    opti.set_value(omega_cce_param, omega_cce);
    opti.set_value(omega_tte_param, omega_tte);
    opti.set_value(thruster_param, thruster_positions);
    opti.subject_to(X(:,1) == x0_param);

    total_cost = 0;

    Jc_inv = MX(inv(Jc));
    Jt_inv = MX(inv(Jt));

    for k = 1:N
        xk_var = X(:,k); uk = U(:,k);
        R_ct = reshape(xk_var(1:9), 3, 3);
        omega_ct = xk_var(10:12);
        rho_ct = xk_var(13:15);
        rhodot_ct = xk_var(16:18);
        xi = xk_var(19:21);
        fps = reshape(uk, 3, 8);

        % Dynamics
        F_total = sum(fps, 2);
        Tau_total = MX.zeros(3,1);
        for i = 1:8
            Tau_total = Tau_total + cross(thruster_param(:,i), fps(:,i));
        end

        omega_skew = [0 -omega_ct(3) omega_ct(2);
                      omega_ct(3) 0 -omega_ct(1);
                      -omega_ct(2) omega_ct(1) 0];
        R_dot = R_ct * omega_skew;

        omega_dot = R_ct*Jc_inv*(cross_casadi(Jc*omega_cce_param, omega_cce_param) + Tau_total) ...
                    - Jt_inv*(cross_casadi(Jt*omega_tte_param, omega_tte_param)) ...
                    + cross_casadi(omega_ct, R_ct * omega_cce_param);
        rot_pc = R_ct * r_pc_c;
        rho_ddot = F_total/m_c + cross_casadi(omega_dot, rot_pc) + cross_casadi(omega_ct, cross_casadi(omega_ct, rot_pc));
        xi_dot = rho_ct;
        x_dot = [reshape(R_dot,9,1); omega_dot; rhodot_ct; rho_ddot; xi_dot];

        % Cost
        F_cost = uk.' * uk; % Simplified cost
        total_cost = total_cost + F_cost;

        x_next = xk_var + Ts * x_dot;
        opti.subject_to(X(:,k+1) == x_next);
    end

    opti.minimize(total_cost);

    % Constraints
    for k = 1:N
        for j = 1:8
            force_j = U((j-1)*3+1:j*3, k);
            % opti.subject_to(norm_2(U((j-1)*3+1:j*3,k)) <= fmax);
            opti.subject_to(sumsqr(force_j)<=fmax^2);
        end
    end

    for k = 1:N+1
        % Bound relative position (within 200m)
        opti.subject_to(-200 <= X(13:15,k) <= 200);
        % Bound relative velocity (within 10 m/s)
        opti.subject_to(-10 <= X(16:18,k) <= 10);
        % Bound angular velocity (within 0.1 rad/s)
        opti.subject_to(-0.1 <= X(10:12,k) <= 0.1);
    end

    opts = struct();
    opts.ipopt.print_level = 0;
    opts.print_time = false;
    opts.ipopt.max_iter = 1000;
    opts.ipopt.tol = 1e-6;
    opti.solver('ipopt', opts);

    %initial guess
    opti.set_initial(X, repmat(xk, 1, N+1));
    opti.set_initial(U, 0.1*ones(24, N));

    try
        sol = opti.solve();
        u_opt = sol.value(U(:,1));
        fprintf('NMPC solved successfully. Max force: %.3f N\n', max(abs(u_opt)));

    catch ME
        disp(['❌ NMPC solve failed: ' ME.message]);
        try
            valX = opti.debug.value(X); valU = opti.debug.value(U);
            disp('⚠️  Debug X(:,1):'); disp(valX(:,1));
            disp('⚠️  Debug U(:,1):'); disp(valU(:,1));
        catch inner
            disp('⚠️  Could not extract debug values.');
        end
        u_opt = zeros(24,1); % fallback
    end
end
