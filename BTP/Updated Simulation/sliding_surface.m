function s = sliding_surface(eta_tilde, xi_tilde)
params

%C=[c1*eye(3), zeros(3); zeros(3),c2*eye(3)];

s= xi_tilde + C*eta_tilde;
end