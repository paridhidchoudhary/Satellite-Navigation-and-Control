function [ se3 ] = R6_se3( R6 )
% R6 is a 6x1 Column Vector of the form = [wx ; wy ; wz ; vx ; vy ; vz]

w = R6(1:3);
v = R6(4:6);

R = [   skew(w)     v  ;...
      zeros(1,3)    0 ];

se3 = R;     % Rotation Matrix
end