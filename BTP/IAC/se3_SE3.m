function [ SE3 ] = se3_SE3( se3 )
% se3_SE3 Exponential Mapping from Lie Algebra to Lie Group
% se3 is a 4x4 square matrix of the form = [ wx v'; 0 1]

wx = se3(1:3,1:3);   % wx is skew symmetricmatrix form of w
u = se3(1:3,4);      % u is translational velocity

w = [ -wx(2,3) ; wx(1,3) ; -wx(1,2)];      

theta = sqrt(w'*w);

if(theta~=0)
    A=sin(theta)/theta;
    B=(1-cos(theta))/(theta^2);
    C=(theta-sin(theta))/(theta^3);
else
    A=0;
    B=0;
    C=0;  
end
R = eye(3) + (A*skew(w)) + (B*(skew(w)*skew(w)));
V = eye(3) + B*skew(w) + C*(skew(w)*skew(w));

Vp = V*u;

SE3 = zeros(4);
SE3(1:3,1:3) = R;
SE3(1:3,4) = Vp;
SE3(4,4) = 1;
end