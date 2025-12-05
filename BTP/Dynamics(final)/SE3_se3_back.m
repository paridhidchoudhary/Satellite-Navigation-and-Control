function [se3] = SE3_se3_back(SE3)
% SE3_se3_back Logarithm Mapping from Lie Group to Lie Algebra
% se3 is a 4x4 square matrix of the form = [wx v ; 0 1]
% SE3 is a 4x4 matrix of the form=[R u ; 0 1]

R = SE3(1:3,1:3);
theta = acos((trace(R)-1)/2);

if theta == 0
    wx = zeros(3,3);
else
    wx = (theta/(2*sin(theta)))*(R-R');
end

if (theta==0)
    Vin=eye(3);
else
    A = theta*cos(theta/2);
    B = sin(theta/2);
    Vin = eye(3) - (1/2)*wx + (1/(theta^2))*(1-(A/(2*B)))*(wx*wx);
end

u = Vin*SE3(1:3,4);
se3 = [     wx      u  ;...
        zeros(1,3)  1 ];
end

