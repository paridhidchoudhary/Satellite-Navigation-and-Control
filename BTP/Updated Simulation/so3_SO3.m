function [ SO3 ] = so3_SO3( so3 )
% so3_SO3 Exponential Mapping
% so3 is a 1x6 Column Vector of the form = [thx ; thy ; thz]

% n = so3/norm(so3);
th = so3'*so3;

% R = [   cos(th)+(n(1)^2)*(1-cos(th))          n(1)*n(2)*(1-cos(th))-n(3)*sin(th)    n(1)*n(3)*(1-cos(th))+n(2)*sin(th)  ;...
%       n(1)*n(2)*(1-cos(th))+n(3)*sin(th)        cos(th)+(n(2)^2)*(1-cos(th))        n(1)*n(3)*(1-cos(th))-n(1)*sin(th)  ;...
%       n(1)*n(3)*(1-cos(th))-n(2)*sin(th)      n(2)*n(3)*(1-cos(th))+n(1)*sin(th)       cos(th)+(n(3)^2)*(1-cos(th))    ];

A = sin(th)/th;
B = (1-cos(th))/(th^2);
thx = [  0    -so3(3)    so3(2);...
        so3(3)    0      -so3(1);...
       -so3(2)   so3(1)    0];
   
R = eye(3) + A*thx + B*thx*thx;
SO3 = R;     % Rotation Matrix
end