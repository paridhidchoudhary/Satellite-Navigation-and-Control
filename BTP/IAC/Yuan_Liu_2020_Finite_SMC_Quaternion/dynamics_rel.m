function xdot = dynamics_rel(x, u, target)
% x      = [qv; ω; r; v]      relative state (12x1)
% u      = [τ;  F]            control input (6x1)
% target = 12x1 vector, target state at this horizon step

param;  % brings J_c, m_c, mu_earth, etc. into scope

% Unpack states and inputs
qv    = x(1:3);     omega = x(4:6);
r     = x(7:9);     v     = x(10:12);
tau   = u(1:3);     F     = u(4:6);

% Extract target position (inertial)
if numel(target) >= 9
    r_t = target(7:9);
else
    error('Target state vector must have at least 9 elements (position).');
end

% Attitude dynamics
q0      = sqrt(max(1 - dot(qv, qv), 0));
qv_dot  = 0.5 * (q0 * eye(3) + skew(qv)) * omega;
omega_dot = J_c \ (tau - cross(omega, J_c*omega));

% Differential gravity
mu      = mu_earth;
a_c     = -mu / norm(r_t + r)^3 * (r_t + r);   % chaser gravity
a_t     = -mu / norm(r_t)^3 * r_t;             % target gravity
a_rel   = a_c - a_t;

% Rotational coupling terms
coriolis      = -2 * cross(omega, v);
%centrifugal   = -cross(omega, cross(omega, r));

% Translational dynamics (with rotational coupling)
r_dot   = v - cross(omega, r);  % you may omit -cross(omega, r) if not needed
v_dot   = F / m_c + a_rel + coriolis ;

xdot    = [qv_dot; omega_dot; r_dot; v_dot];
end
