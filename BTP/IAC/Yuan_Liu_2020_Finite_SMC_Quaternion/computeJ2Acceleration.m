function a_J2 = computeJ2Acceleration(r)
    
param;

    x = r(1); y = r(2); z = r(3);
    r_norm = norm(r);
    
    factor = -(1.5 * J2 * mu_earth * (Re^2)) / (r_norm^5);
    
    a_x = factor * x * ((1 - 5 * (z^2)) / (r_norm^2));
    a_y = factor * y * ((1 - 5 * (z^2)) / (r_norm^2));
    a_z = factor * z * ((3 - 5 * (z^2)) / (r_norm^2));
    
    a_J2 = [a_x; a_y; a_z];
end
