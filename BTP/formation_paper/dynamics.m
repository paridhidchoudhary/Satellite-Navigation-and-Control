function state_dot= dynamics(t,state, i,dt)
    params;
    
    %state = reshape(state, 22, 4);
    xi_tilde=state(1:6);
    eta_tilde=state(7:12);
    xi_0= state(13:18);
    %g0=state(19:22,1:4);

    phi_g = [0; 0; -9.81*m(i); zeros(3,1)];
    phi_d = 1e-5 * sin(2*pi*round(t/dt)*dt/5400)*[1.92; -1.986; -1.517; 2; -2; -1.5];

    s = sliding_surface(eta_tilde, xi_tilde);
    plot(t,s);
    xlabel('time')
    ylabel('sliding variable')
    % ad_eta = adjoint_op(eta_tilde);
    % 
    % exp_ad_eta = expm(ad_eta);
    eta_tilde_v= v_op(eta_tilde);
    ad_h_inv=Ad_h_inv(eta_tilde, hf);
    size(ad_h_inv)
    G_eta = G_matrix(eta_tilde);

    % xi_0_dot = zeros(6,1);

    %xi_0_dot:
    xi_0_dot=I0\ (phi_g + adj_xi_star(xi_0)*I0*xi_0);


    %g0_dot= g0*v_op(xi_0);

    xi_k= xi_tilde + ad_h_inv*xi_0;

    phi_c = -phi_g - adj_xi_star(xi_k)*I(:,:,i)*xi_k - I(:,:,i)*(adjoint_op(xi_k)*ad_h_inv ...
        *xi_0 - ad_h_inv*xi_0_dot + C*G_eta*(xi_tilde)) - P*s - K*tanh(s);

    norm_phi_c = norm(phi_c(4:6));
    size(phi_c)

    if norm_phi_c > phi_m(i)
        phi_c(4:6) = phi_m (i)* phi_c(4:6)/norm_phi_c;
    end

    xi_tilde_dot = I(:,:,i)\(phi_g + phi_c + phi_d - adj_xi_star(xi_k)*I(:,:,i)*xi_k);
    %xi_tilde_dot= -C*G_eta*xi_tilde(:,k,i) - inv(I(:,:,i))*(P*s + K*sat(s));
    %eta_tilde_dot = G_eta * xi_tilde(:,k,i);
    eta_tilde_dot= G_eta * (-C*eta_tilde);

    state_dot=[xi_tilde_dot; eta_tilde_dot; xi_0_dot]

end