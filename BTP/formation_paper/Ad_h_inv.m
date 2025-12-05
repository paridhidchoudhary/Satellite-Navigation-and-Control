function ad_h_inv= Ad_h_inv(eta_tilde, hf)
h=hf*se3_SE3(eta_tilde(4:6), eta_tilde(1:3));
h_inv=lie_inv(h);
ad_h_inv=adjoint_op(h_inv);
end
