function [ R6 ] = se3_R6_back( se3 )
% R6 is a 6x1 Column Vector of the form = [wx ; wy ; wz ; vx ; vy ; vz]

R = se3(1:3,1:3);
v = se3(1:3,4);

theta = acos((trace(R)-1)/2);

if theta == 0
    w = zeros(3,1);
else
    w = (theta/(2*sin(theta)))*[ R(3,2)-R(2,3) ; R(1,3)-R(3,1) ; R(2,1)-R(1,2)];
end

R6 = [ w ; v];

end