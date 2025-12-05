function v=v_op(m)
theta=m(1:3);
beta=m(4:6);
v=[cross_pdt(theta), beta;0,0,0,1];
end