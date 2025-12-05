function a_g = compute_gravity_acceleration(r_body, r_inertial)

param;

r_norm = norm(r_inertial);
a_g = -(mu_earth/(r_norm^3))*r_body;

end