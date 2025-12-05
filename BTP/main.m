params;
eta=zeros(6,N);
xi=zeros(6,N);
eta(:,1)=[30*pi/180;0;0;-1000;1500;-500];
xi(:,1)= [0.001*[-0.06;0.26;0.05]; 3.17;1.24;1.82];

for k=1:N-1
    s=sliding_surface(eta(:,k),xi(:,k));
    ad_eta=adjoint_op(eta(:,k));
    exp_ad_eta=expm(ad_eta);
    G_eta=G_matrix(eta(:,k));

    xi_0=[0.001; 0.001; 0.001; 0.1; 0.1; 0.1];
    xi_0_dot=zeros(6,1);

    phi_g=[0; 0; -9.81*m; zeros(3,1)];
    phi_d= 1e-5 * sin(2*pi*k*dt/5400)*[1.92; -1.986; -1.517; 2; -2; -1.5];

    phi_c= -phi_g - adj_xi_star(xi(:,k))*I*xi(:,k) - I*(adjoint_op(xi(:,k))*exp_ad_eta ...
        *xi_0 - exp_ad_eta*xi_0_dot + C*G_eta*(xi(:,k)- exp_ad_eta*xi_0)) - P*s - K*sign(s);

    norm_phi_c=norm(phi_c(4:6));
    if norm_phi_c> phi_m
        phi_c(4:6)= phi_m * phi_c(4:6)/norm_phi_c;
    end

    xi_dot=I\(phi_g + phi_c + phi_d -adj_xi_star(xi(:,k))*I*xi(:,k));
    eta_dot= G_eta *xi(:,k);

    xi(:,k+1)= xi(:,k)+xi_dot*dt;
    eta(:,k+1)=eta(:,k) +eta_dot*dt;
end

figure;
subplot(2,1,1);
plot(0:dt:T-dt, eta(1:3,:));
title('Attitude Error');
xlabel('Time(s)');
ylabel('Error(rad)');
legend('\theta_x', '\theta_y', '\theta_z' );

subplot(2,1,2);
plot(0:dt:T-dt, eta(4:6,:));
title('Position Error');
xlabel('Time(s');
ylabel('Error(m)');
legend('x','y','z');
    








