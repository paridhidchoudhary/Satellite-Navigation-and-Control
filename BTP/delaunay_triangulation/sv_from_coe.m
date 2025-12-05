function [r, v] = sv_from_coe(coe,mu)
h    = coe(1);
e    = coe(2);
RA   = coe(3);
incl = coe(4);
w    = coe(5);
TA   = coe(6);

% position vector in the perifocal frame in km
rp = (h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
% velocity vector in the perifocal frame in km/s
vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);

% Rotation matrix about the z-axis through Right Accession of Ascending Node angle
R3_W = [ cos(RA)  sin(RA)  0
        -sin(RA)  cos(RA)  0
            0        0     1];

% Rotation matrix about the x-axis through inclination angle
R1_i = [1       0          0
        0   cos(incl)  sin(incl)
        0  -sin(incl)  cos(incl)];

% Rotation matrix about the z-axis through argument of perigee angle
R3_w = [ cos(w)  sin(w)  0 
        -sin(w)  cos(w)  0
           0       0     1];

% Transformation Rotation Matrix from perifocal frame to ECI frame 
Q_pX = (R3_w*R1_i*R3_W)';

% position vector in the ECI frame in km
r = Q_pX*rp;
% VELOCITY vector in the ECI frame in km/s
v = Q_pX*vp;

r = r';
v = v';

end