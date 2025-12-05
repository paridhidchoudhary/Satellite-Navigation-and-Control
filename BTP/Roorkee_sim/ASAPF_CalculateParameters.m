function [zeta, eta_max, eta_min] = ASAPF_CalculateParameters(Q_star, dg_star, d_safe, v_max, a_max)   
    zeta = sqrt(2 * a_max * dg_star) / dg_star;

    eta_min = CalculateEta(d_safe, zeta, Q_star, dg_star, (-1.0)*v_max);
    eta_max = CalculateEta(Q_star-d_safe, zeta, Q_star, dg_star, v_max);
end

function eta = CalculateEta(d_obstacle, zeta, Q_star, dg_star, v_max)   
    eta = (d_obstacle^2 * Q_star * (v_max - dg_star * zeta)) / (d_obstacle - Q_star);
end