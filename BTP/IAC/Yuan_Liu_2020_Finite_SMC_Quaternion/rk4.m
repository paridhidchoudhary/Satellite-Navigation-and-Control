% function [Y,u,S] = rk4(X,tspan,y0, )
% % rk4 function for solving equations using 4TH ORDER RK METHOD.
% % INPUTS
% % X function handle
% % tspan simulation running time
% % y0 intial state
% %
% % OUTPUT
% % y solution matrix
% %
% h = diff(tspan);
% y0 = y0(:);   % Make a column vector.
% neq = length(y0);
% N = length(tspan);
% Y = zeros(neq,N);%output vector
% F = zeros(neq,4);
% u = zeros(N,6);
% S = zeros(N,6);
% Y(:,1) = y0;
% 
% 
% for i = 2:N
%  ti = tspan(i-1);
%  hi = h(i-1);
%  yi = Y(:,i-1);
%  [F(:,1),u(i-1,:),S(i-1,:)] = feval(X,ti,yi);
%  [F(:,2),~,~] = feval(X,ti+0.5*hi,yi+0.5*hi*F(:,1));
%  [F(:,3),~,~] = feval(X,ti+0.5*hi,yi+0.5*hi*F(:,2)); 
%  [F(:,4),~,~] = feval(X,tspan(i),yi+hi*F(:,3));
%   Y(:,i) = yi + (hi/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4));
% end
%  [~,u(N,:),~] = (feval(X,tspan(N),Y(:,N)));
%  [~,~,S(N,:)] = (feval(X,tspan(N),Y(:,N)));
%  Y = Y.';
% end

function [Y, u, S] = rk4(X, tspan, y0, extra_input)
% rk4: 4th-order Runge-Kutta solver for dynamic systems
% INPUTS:
%   X            - function handle for dynamics: X(t, y, extra_input)
%   tspan        - vector of time steps
%   y0           - initial state (column vector)
%   extra_input  - static input (e.g., 3x1 vector or struct) used in dynamics
%
% OUTPUTS:
%   Y            - solution trajectory (each row = state at time i)
%   u            - control inputs at each time step
%   S            - sliding surfaces at each time step

    h = diff(tspan);
    y0 = y0(:);   % Ensure y0 is a column vector
    neq = length(y0);
    N = length(tspan);
    Y = zeros(neq, N); % State trajectory
    F = zeros(neq, 4); % RK4 slopes
    u = zeros(N, 6);   % Control input storage
    S = zeros(N, 6);   % Sliding surfaces
    Y(:,1) = y0;

    for i = 2:N
        ti = tspan(i-1);
        hi = h(i-1);
        yi = Y(:,i-1);

        [F(:,1), u(i-1,:), S(i-1,:)] = X(ti, yi, extra_input);
        [F(:,2), ~, ~] = X(ti + 0.5*hi, yi + 0.5*hi*F(:,1), extra_input);
        [F(:,3), ~, ~] = X(ti + 0.5*hi, yi + 0.5*hi*F(:,2), extra_input);
        [F(:,4), ~, ~] = X(tspan(i),   yi + hi*F(:,3), extra_input);

        Y(:,i) = yi + (hi/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4));
    end

    [~, u(N,:), ~] = X(tspan(N), Y(:,N), extra_input);
    [~, ~, S(N,:)] = X(tspan(N), Y(:,N), extra_input);

    Y = Y.';
end
% function x_next = rk4(f, x, u, h, param)
% % Fourth‑order Runge‑Kutta step
% k1 = f(x,           u, param);
% k2 = f(x+0.5*h*k1,  u, param);
% k3 = f(x+0.5*h*k2,  u, param);
% k4 = f(x+    h*k3,  u, param);
% x_next = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
% end
