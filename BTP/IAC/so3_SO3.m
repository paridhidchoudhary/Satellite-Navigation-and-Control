function [ SO3 ] = so3_SO3( so3 )
% so3_SO3 Exponential Mapping
% so3 is a 1x3 Column Vector of the form = [thx ; thy ; thz]

th = so3'*so3;

A = sin(th)/th;
B = (1-cos(th))/(th^2);
thx = [  0    -so3(3)    so3(2);...
        so3(3)    0      -so3(1);...
       -so3(2)   so3(1)    0];
   
R = eye(3) + A*thx + B*thx*thx;
SO3 = R;     % Rotation Matrix
end