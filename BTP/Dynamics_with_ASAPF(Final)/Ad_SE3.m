function [ Ad ] = Ad_SE3( SE3 )

R = SE3(1:3,1:3);
b = SE3(1:3,4);

% Adjoint operation of SE(3) 4x4 matrix to 6x6 matrix
Ad = [     R      ,  zeros(3,3)  ;...
       skew(b)*R  ,      R      ];
   
   
end