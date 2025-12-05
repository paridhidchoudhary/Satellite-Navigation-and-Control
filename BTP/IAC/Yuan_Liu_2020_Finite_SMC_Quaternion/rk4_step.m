function x_next = rk4_step(f, x, u, h)
% One classic 4th‑order Runge–Kutta step.
%   f   – function handle  f(x,u)  →  ẋ
%   x   – current state
%   u   – control input held constant over [0,h]
%   h   – step size
k1 = f(           x,            u);
k2 = f( x + 0.5*h*k1,           u);
k3 = f( x + 0.5*h*k2,           u);
k4 = f( x +     h*k3,           u);
x_next = x + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end
