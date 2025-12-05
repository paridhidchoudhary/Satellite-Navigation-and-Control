% CasADi compatible version of dynamics
function dxdt = quad_dynamics_casadi(state, taux, tauy, tauz, T)
    import casadi.*
    
    param;
    
    % Extract state variables
    phi = state(1);
    theta = state(2);
    psi = state(3);
    phi_d = state(7);
    theta_d = state(8);
    psi_d = state(9);
    
    % Initialize state derivative vector
    dxdt = SX.zeros(12,1);
    
    % State derivatives (same as in your quad_dynamics function)
    dxdt(1) = phi_d;
    dxdt(2) = theta_d;
    dxdt(3) = psi_d;
    dxdt(4) = state(10); % x_dot
    dxdt(5) = state(11); % y_dot
    dxdt(6) = state(12); % z_dot (was set to 0 in your function)
    
    % Angular accelerations
    dxdt(7) = (taux/ix) - ((iz - iy)/ix)*theta_d*psi_d;
    dxdt(8) = (tauy/iy) - ((ix - iz)/iy)*phi_d*psi_d;
    dxdt(9) = (tauz/iz) - ((iy - ix)/iz)*phi_d*theta_d;
    
    % Linear accelerations
    dxdt(10) = - (T/m)*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi));
    dxdt(11) = - (T/m)*(cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi));
    dxdt(12) = G - (T/m)*cos(phi)*cos(theta);
end