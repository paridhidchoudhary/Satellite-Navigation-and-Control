function inv=lie_inv(m)
R=m(1:3,1:3);
t=m(1:3,4);
inv=[R.', -R'*t;0,0,0,1];
end