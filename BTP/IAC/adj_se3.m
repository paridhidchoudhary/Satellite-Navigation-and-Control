function [ adj ] = adj_se3( se3 )

wx = se3(1:3,1:3);
v = se3(1:3,4);

% Adjoint operation of se(3) 4x4 matrix to 6x6 matrix
adj = [    wx  ,  zeros(3,3)  ;...
        skew(v)  ,   wx      ];
end