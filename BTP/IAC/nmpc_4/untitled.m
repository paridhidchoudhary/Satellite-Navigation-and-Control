
clear; clc; import casadi.*

% ------------- USER PARAMETERS --------------------------------
N     = 20;                   % shooting intervals   (dt = 8 s)
dt    = 8;                    % [s] control period   (T_h = 480 s)
T_sim = 240;                  % simulate 10 minutes
fmax  = 20;                   % N   per thruster (component box)
omega_t = [0;0;0];            % target angular rate

% ------------- DYNAMICS / COST PLACEHOLDERS -------------------
fdyn  = @fdyn_user;           
fcost = @fcost_user;          

nx = 21;  nu = 24;            % state / input sizes
Wr = 1e-2*ones(nu,1);         % slew vector (diagonal)

% ------------- THRUSTER CONE MATRICES -------------------------
[r_tk, j_nom] = getThrusterConfig();         
[Acell, bcell] = buildCones(j_nom, pi/4);    % helper below

% ------------- NLP VARIABLES (MX) -----------------------------
X = SX.sym('X', nx, N+1);
U = SX.sym('U', nu, N);

J  = 0;   G = {};  lbg = [];  ubg = [];

for k = 1:N
    disp(k);
    xk = X(:,k);  uk = U(:,k);

    J = J + fcost(xk,uk);
    if k>1
        Gslew = uk - U(:,k-1);
        J = J + sum(Wr.*(Gslew.^2));
    end

    % thrust limits & cone
    fMat = reshape(uk,3,8);
    for j = 1:8
        fj = fMat(:,j);
        G{end+1} = fj.'*fj;                        lbg=[lbg;0]; ubg=[ubg;fmax^2];
        G{end+1} = Acell{j}*fj;                   lbg=[lbg;-inf(4,1)]; ubg=[ubg;bcell{j}];
    end

    % Euler forward model
    x_next = xk + dt * fdyn(xk,uk,omega_t);
    G{end+1} = X(:,k+1) - x_next;                lbg=[lbg;zeros(nx,1)]; ubg=[ubg;zeros(nx,1)];
end

% terminal quadratic on ρ and ρ̇
P = 10*eye(6);  J = J + X(13:18,end).'*P*X(13:18,end);

% ------------- SOLVER ----------------------------------------
w  = [X(:);U(:)];   G = vertcat(G{:});
lbx = -inf(size(w)); ubx = inf(size(w));
lbx(nx*(N+1)+1:end) = -fmax;   ubx(nx*(N+1)+1:end) = fmax;

solver = nlpsol('s','ipopt',struct('x',w,'f',J,'g',G),...
        struct('ipopt',struct('print_level',0,'hessian_approximation','limited-memory')));

% ------------- CLOSED‑LOOP SIM --------------------------------
Xsim = zeros(nx,T_sim/dt+1); Ulog = zeros(nu,T_sim/dt);
x0 = zeros(nx,1); x0([1 5 9]) = 1;    % identity attitude
x0(13:15) = [-1.1394  113.0847  -7.6204];   

Xsim(:,1)=x0;  Xguess=repmat(x0,1,N+1); Uguess=zeros(nu,N);

fd_num = @(x,u) full(fdyn(x,u,omega_t));   

fprintf('Running NMPC …\n');
for k=1:size(Ulog,2)
    lbx(1:nx)=Xsim(:,k); ubx(1:nx)=Xsim(:,k);

    sol = solver('x0',[Xguess(:);Uguess(:)],'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg);
    wopt = full(sol.x);
    Xopt = reshape(wopt(1:nx*(N+1)),nx,N+1);
    Uopt = reshape(wopt(nx*(N+1)+1:end),nu,N);

    u = max(min(Uopt(:,1),fmax),-fmax);          
    Ulog(:,k)=u;

    % 4 sub‑steps Euler
    x = Xsim(:,k);
    for s=1:4, x = x + (dt/4)*fd_num(x,u); end
    x(1:9)=gs_project(x(1:9)); Xsim(:,k+1)=x;

    Xguess=[Xopt(:,2:end) Xopt(:,end)];
    Uguess=[Uopt(:,2:end) Uopt(:,end)];
end
fprintf('Done.\n');

% figure; plot(0:dt:T_sim,Xsim(13:15,:).'); title('\rho');

t = 0:dt:T_sim;                  
t_ctrl = 0:dt:T_sim-dt;       

%% Translation error ρ and velocity error ρ̇
figure;
subplot(2,1,1);
plot(t, Xsim(13:15,:)); grid on;
title('\rho – position error'); ylabel('m');
legend('\rho_x','\rho_y','\rho_z');

subplot(2,1,2);
plot(t, Xsim(16:18,:)); grid on;
title('\rhȯ – velocity error'); ylabel('m/s');
legend('\rhȯ_x','\rhȯ_y','\rhȯ_z');
xlabel('time (s)');

%% Attitude error (angle) and angular‑rate error ω̃
theta = zeros(1,length(t));                
for k = 1:length(t)
    Rct = reshape(Xsim(1:9,k),3,3);        % 3×3 rotation
    theta(k) = acos( 0.5*(trace(Rct) - 1) ); % radians
end
figure;
subplot(2,1,1);
plot(t, rad2deg(theta)); grid on;
title('\theta – attitude error'); ylabel('deg');

subplot(2,1,2);
plot(t, Xsim(10:12,:)); grid on;
title('\omegã – angular‑rate error'); ylabel('rad/s');
legend('\omega_x','\omega_y','\omega_z');
xlabel('time (s)');

%% Thrust magnitude per thruster
figure; hold on;
for j = 1:8
    idx = 3*(j-1)+(1:3);
    plot(t_ctrl, vecnorm(Ulog(idx,:),2,1));
end
grid on; title('‖f_{tk}‖ per thruster'); ylabel('N');
legend(arrayfun(@(j)sprintf('T%d',j),1:8,'uni',0));

% %% 3‑D trajectory (relative frame)
% target_pt  = [0 0 0];          
% chaser_pts = Xsim(13:15,:).' + target_pt;   
% 
% figure;
% plot3(chaser_pts(:,1), chaser_pts(:,2), chaser_pts(:,3),'b-','LineWidth',1.5);
% hold on;
% plot3(chaser_pts(1,1), chaser_pts(1,2), chaser_pts(1,3),'go','MarkerSize',8,'LineWidth',2);
% plot3(target_pt(1),target_pt(2),target_pt(3),'rx','MarkerSize',10,'LineWidth',2);
% grid on; axis equal;
% xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
% title('3‑D rendezvous trajectory');
% legend('path','start','target rock');


% ===== helper functions =====
function [Acell,bcell]=buildCones(J,ang)
for k=1:8
    a=null(J(:,k).'); dirs=[a,-a]; Ak=zeros(4,3);
    for p=1:4, Ak(p,:)=J(:,k).'*cos(ang)+dirs(:,p).'*sin(ang); end
    Acell{k}=Ak; bcell{k}=zeros(4,1);
end, end

function r9=gs_project(r9)
R=reshape(r9,3,3); c1=R(:,1)/norm(R(:,1));
c2=R(:,2)-c1*(c1.'*R(:,2)); c2=c2/norm(c2); c3=cross(c1,c2);
r9=reshape([c1 c2 c3],9,1); end


function xdot = fdyn_user(x,u,omega_t)
    persistent F
    if isempty(F), F = f_dyn(); end
    xdot = F(x,u,omega_t);
end

function L = fcost_user(x,u)
    persistent F           
    if isempty(F), F = f_cost(); end
    L = F(x,u);            
end

