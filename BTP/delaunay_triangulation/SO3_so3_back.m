function [so3] = SO3_so3_back(SO3)
% SO3_so3_back Logarithm Mapping from Lie Group to Lie Algebra
% so3 is a Column Vector of the form = [thx ; thy ; thz]
% SO3 is a 3x3 rotation matrix

theta = acos((trace(SO3)-1)*(1/2));

if theta == 0
    so3 = zeros(3,1);
else
    so3 = (theta/(2*sin(theta)))*[ SO3(3,2)-SO3(2,3) ; SO3(1,3)-SO3(3,1) ; SO3(2,1)-SO3(1,2)];
end
end