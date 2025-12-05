function tau_GG = computeTauGG(r_body, r_inertial, J)

param;

r_norm = norm(r_inertial);
rb_norm = norm(r_body);
r_hat = r_body / rb_norm;
    
tau_GG = (3 * mu_earth / (r_norm^3)) * cross(r_hat, J * r_hat);

end
