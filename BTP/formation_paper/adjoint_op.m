function ad_xi=adjoint_op(xi)
omega = xi(1:3);
 v = xi(4:6);
 ad_xi = [cross_pdt(omega), zeros(3); cross_pdt(v), cross_pdt(omega)];
end