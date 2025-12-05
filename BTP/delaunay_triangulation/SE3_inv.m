function [inv_SE3] = SE3_inv(SE3)
% SE3_se3_back Logarithm Mapping from Lie Group to Lie Algebra
% se3 is a 1x6 Column Vector of the form=[v1 v2 v3 w1 w2 w3]
% SE3 is a 4x4 matrix of the form=[R u ; 0 1]

R = SE3(1:3,1:3);
t = SE3(1:3,4);

inv_SE3= [    R'       -R'*t  ;...
          zeros(1,3)      1  ];
      
end
