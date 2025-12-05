import casadi.*
% ---------- CONSTANTS ----------
T   = 960;             % prediction horizon
N   = 320;             % grid points  (dt = T/N)
dt  = T / N;

nx  = 21;
nu  = 8;

% ---------- BUILD THE NLP ONCE ----------
fdyn  = f_dyn();
fcost = f_cost();
[X,U,J,g,lbg,ubg] = build_nlp(fdyn,fcost,nx,nu,N,dt);   % your existing code, wrapped in a function

w     = [X(:); U(:)];
nlp   = struct('x', w, 'f', J, 'g', g);
opts  = struct('ipopt', struct('print_level',0, 'max_iter', 500));
solver = nlpsol('solver','ipopt', nlp, opts);

% ------- CONTROL BOUNDS -------
u_max   = 15;
lbx     = -inf(size(w));
ubx     =  inf(size(w));
u_start = nx*(N+1) + 1;
lbx(u_start:end) = 0;
ubx(u_start:end) = u_max;

% ------- TERMINAL COST WEIGHT -------
Pterm   = 1e4 * eye(6);   % (ρ,ρdot) only

% ------- SIMULATION ARRAYS -------
nSteps  = round(T / dt);          % Tsim = total *real‑time* you want to fly
Xhist   = zeros(nx, nSteps+1);
Uhist   = zeros(nu, nSteps);

% ------- INITIAL STATE -------
x_now   = [1 0 0 0 1 0 0 0 1  ...   % identity C_bn in quaternion/9‑DOF form
          zeros(1,3) ...            % ω
          -1.1394 113.0847 -7.6204  % ρ
           0.0081  0.0039  1.3490]; % ρdot
Xhist(:,1) = x_now.';

% ---------- GUESS FOR FIRST OPTIMISATION ----------
Xguess = repmat(x_now.', 1, N+1);
Uguess = 0.5 + 0.5*rand(nu,N);      % small bias → helps IPOPT

for kSim = 1:nSteps
    
    % --- PACK THE WARM‑START VECTOR ---
    w0 = [Xguess(:); Uguess(:)];
    
    % --- SOLVE NLP WITH *CURRENT* STATE CONSTRAINT ---
    % Impose x(0) = x_now  (overwrite first nx variables in lbx/ubx)
    lbx(1:nx) = x_now.';           % equality: lower = upper
    ubx(1:nx) = x_now.';
    
    sol = solver('x0', w0, ...
                 'lbx', lbx, 'ubx', ubx, ...
                 'lbg', lbg, 'ubg', ubg);
    
    w_opt = full(sol.x);
    X_opt = reshape(w_opt(1:nx*(N+1)), nx, N+1);
    U_opt = reshape(w_opt(nx*(N+1)+1:end), nu, N);
    
    % --- APPLY FIRST CONTROL MOVE ---
    u_now         = U_opt(:,1);
    Uhist(:,kSim) = u_now;
    
    % --- PLANT INTEGRATION (1 step of “true” dynamics) ---
    k1 = fdyn(x_now,           u_now);
    k2 = fdyn(x_now+0.5*dt*k1, u_now);
    k3 = fdyn(x_now+0.5*dt*k2, u_now);
    k4 = fdyn(x_now+    dt*k3, u_now);
    x_next  = x_now + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    Xhist(:,kSim+1) = x_next.';
    
    % --- SHIFT INITIAL GUESS FOR NEXT ITERATION ----------
    %  • drop first column, append last column
    Xguess = [X_opt(:,2:end) , X_opt(:,end)];
    %  • same for controls, keep last value as tail (or zeros)
    Uguess = [U_opt(:,2:end) , U_opt(:,end)];
    
    % --- UPDATE FOR NEXT LOOP ---
    x_now = x_next;
end
%% ==========================  P L O T S  ============================ %
%  Assumes you have                                                  
%       • Xhist  (nx × [nSteps+1])   – closed‑loop state history      
%       • Uhist  (nu ×  nSteps)      – control history                
%       • dt                          – sampling period               
%       • nSteps                      – # of control intervals        
%  and that state indices are:                                       
%       10:12  → ω (angular‑vel error)                               
%       13:15  → ρ  (position error)                                 
%       16:18  → ρ̇ (velocity error)                                 
% =================================================================== %

t_state   = 0 : dt :  nSteps*dt;          % length = nSteps+1
t_control = 0 : dt : (nSteps-1)*dt;       % length = nSteps

%% 1. Error plots (ρ, ρ̇, ω)
figure('Name','Error histories','NumberTitle','off');

subplot(3,1,1);
plot(t_state, Xhist(13:15,:));
title('\rho – Relative Position Error'); ylabel('[m]');
legend('\rho_x','\rho_y','\rho_z','Location','best'); grid on;

subplot(3,1,2);
plot(t_state, Xhist(16:18,:));
title('\dot{\rho} – Relative Velocity Error'); ylabel('[m/s]');
legend('\dot{\rho}_x','\dot{\rho}_y','\dot{\rho}_z','Location','best'); grid on;

subplot(3,1,3);
plot(t_state, Xhist(10:12,:));
title('\omega – Angular Velocity Error'); xlabel('Time [s]'); ylabel('[rad/s]');
legend('\omega_x','\omega_y','\omega_z','Location','best'); grid on;

%% 2. 3‑D trajectory in the target frame
figure('Name','Chaser trajectory','NumberTitle','off');
plot3(Xhist(13,:), Xhist(14,:), Xhist(15,:), 'LineWidth',1.2); hold on;
plot3(0,0,0,'r*','MarkerSize',10,'DisplayName','Target');
grid on; axis equal;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Chaser Trajectory in Target Frame');
legend('Path','Target','Location','best'); view(135,30);

%% 3. Thrust history per jet
figure('Name','Thrust per thruster','NumberTitle','off'); hold on;
for jj = 1:nu
    plot(t_control, Uhist(jj,:), 'DisplayName', sprintf('f_{%d}',jj));
end
hold off; grid on;
xlabel('Time [s]'); ylabel('Force [N]');
title('Thrust per Thruster');
legend('Location','bestoutside');


%% 
function [X,U,J,g,lbg,ubg] = build_nlp(fdyn, fcost, nx, nu, N, dt)
% BUILD_NLP  Constructs the symbolic variables, cost and constraints
%            for the coupled‑spacecraft NMPC problem.
%
%   INPUTS
%     fdyn   – CasADi SX‑function  ẋ = fdyn(x,u)
%     fcost  – CasADi SX‑function  ℓ = fcost(x,u)
%     nx     – # state variables
%     nu     – # control inputs
%     N      – # control intervals in the prediction horizon
%     dt     – time step (seconds)
%
%   OUTPUTS
%     X,U    – symbolic matrices (nx×N+1) and (nu×N)
%     J      – total cost (SX scalar)
%     g      – stacked dynamics equality constraints
%     lbg    – lower bounds for g   (vector)
%     ubg    – upper bounds for g   (vector)
%
%   The function reproduces exactly what you had inline in your
%   open‑loop script, just wrapped up so the main file reads cleaner.

    import casadi.*

    %------------- decision variables -----------------------------------%
    X = SX.sym('X', nx, N+1);     % states  (columns 0…N)
    U = SX.sym('U', nu, N);       % inputs  (columns 0…N‑1)

    %------------- stage & running cost ---------------------------------%
    J   = 0;
    Wr  = 1e-2 * eye(nu);         % Δu regularisation weight
    gcs = {};   lbg = [];   ubg = [];

    for k = 1:N
        % stage cost ℓ(x_k,u_k)
        J = J + fcost(X(:,k), U(:,k));

        % Δu smoothing
        if k > 1
            du = U(:,k) - U(:,k-1);
            J  = J + du.' * Wr * du;
        end

        % RK‑4 integration
        xk = X(:,k);   uk = U(:,k);
        k1 = fdyn(xk,               uk);
        k2 = fdyn(xk + 0.5*dt*k1,   uk);
        k3 = fdyn(xk + 0.5*dt*k2,   uk);
        k4 = fdyn(xk +     dt*k3,   uk);
        x_next = xk + dt/6*(k1 + 2*k2 + 2*k3 + k4);

        % dynamics equality constraint  x_{k+1} − ϕ(x_k,u_k) = 0
        gcs{end+1} = X(:,k+1) - x_next;
        lbg = [lbg; zeros(nx,1)];
        ubg = [ubg; zeros(nx,1)];
    end

    %------------- terminal cost on (ρ, ρ̇) -----------------------------%
    Pterm  = 1e4 * eye(6);
    rhoT   = X(13:15,end);
    rhodT  = X(16:18,end);
    J      = J + [rhoT; rhodT].' * Pterm * [rhoT; rhodT];

    %------------- stack constraints ------------------------------------%
    g = vertcat(gcs{:});
end


