function [X_dot,phi_c,s] = Coupled_EOM(t,X)

    params;
    
    xi_tilde(1:6,1) = X(1:6,1);
    eta_tilde(1:6,1) = X(7:12,1);
    xi_0(1:6,1) = X(13:18,1);
    xi_0_dot(1:6,1) = X(19:24,1);


    
    %% Disturbance and gravity components
    phi_g = [1.2; 1.8; 2.1; zeros(3,1)];
    phi_d = [1.92; -1.986; -1.517; 2; -2; -1.5];

    
    I = [   J          zeros(3,3);...
         zeros(3,3)    m*eye(3,3)];
     
    I0 = [   J0         zeros(3,3);...
         zeros(3,3)    m*eye(3,3)];
    
    
    %% Equation of Motions
    eta_tilde_v = R6_se3(eta_tilde);  
    
    hf = [-1 ,  0 ,  0 ,  0;...
           0 , -1 ,  0 ,  0;...
           0 ,  0 , -1 ,  0;...
           0 ,  0 ,  0 ,  1];
       
    h = hf*se3_SE3(eta_tilde_v);
    
    h_inv = SE3_inv(h);
    
    Ad_h_inv = Ad_SE3(h_inv); 
    
    G_eta = G_matrix(eta_tilde);
    
    adj_xi_0 = adj_se3(xi_0);
    adj_star_xi_0 = adj_xi_0';

    xi = xi_tilde + Ad_h_inv*xi_0;
    
    adj_xi = adj_se3(xi);
    adj_star_xi = adj_xi';
    
    
    %% Sliding Mode Control Input
    s = xi_tilde + C*eta_tilde;
   
    
    phi_c = -phi_g - adj_star_xi*I*xi - I*(adj_xi*Ad_h_inv*xi_0 - Ad_h_inv*xi_0_dot + C*G_eta*xi_tilde) - P*s - K*tanh(s);

%     norm_phi_c = norm(phi_c(4:6));
%     size(phi_c)
% 
%     if norm_phi_c > phi_m(i)
%         phi_c(4:6) = phi_m (i)* phi_c(4:6)/norm_phi_c;
%     end

    xi_0_dot_0 = inv(I0)*(phi_g + adj_star_xi_0*I0*xi_0);

    
    %% Differential Equations
     X_dot(1:6,1) = inv(I)*(phi_g + phi_c + phi_d + adj_star_xi*I*xi) + adj_xi*Ad_h_inv*xi_0 - Ad_h_inv*xi_0_dot;    % Relative velocity error xi_tilde
    %X_dot(1:6,1) = -C*G_eta*xi_tilde - inv(I)*(P*s + K*tanh(s));

    X_dot(7:12,1) = -G_eta*C*eta_tilde;     % Relative configuration error eta_tilde

    X_dot(13:18,1) = inv(I0)*(phi_g + adj_star_xi_0*I0*xi_0);     % Leader unified velocity vector xi_0
    
    X_dot(19:24,1) = xi_0_dot_0-xi_0_dot;
    
    disp(t);

end