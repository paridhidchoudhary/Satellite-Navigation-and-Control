function I_xi_tilde_dot = relative_dynamics(I, xi, xi_0, h, phi_g, phi_c, phi_d)
    % I: Inertia matrix
    % xi: Velocity of follower
    % xi_0: Velocity of leader
    % h: Relative configuration
    % phi_g: Gravity forces and moments
    % phi_c: Control forces and moments
    % phi_d: Disturbance forces and moments

    % Calculate Ad(h)^-1
    R = h(1:3, 1:3);
    p = h(1:3, 4);
    Ad_h_inv = [R', -R'*cross_pdt(p); zeros(3), R'];

    % Calculate adξ
    omega = xi(1:3);
    v = xi(4:6);
    ad_xi = [cross_pdt(omega), zeros(3); cross_pdt(v), cross_pdt(omega)];

    % Calculate ad*ξ
    ad_xi_star = ad_xi';

    % Calculate ˙ξ0 (you need to provide this or approximate it)
    xi_0_dot = zeros(6,1); % This is a placeholder. You should calculate or approximate this.

    % Calculate the relative dynamics
    I_xi_tilde_dot = ad_xi_star * I * xi + phi_g + phi_c + phi_d + ...
                     I * (ad_xi * Ad_h_inv * xi_0 - Ad_h_inv * xi_0_dot);
end
