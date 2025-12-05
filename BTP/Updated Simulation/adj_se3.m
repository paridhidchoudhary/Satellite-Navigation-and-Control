function adj = adj_se3(se3)

w = se3(1:3,1);
v = se3(4:6,1);

% Adjoint operation
adj = [ skew(w)  ,  zeros(3,3)  ;...
        skew(v)  ,   skew(w)   ];
end