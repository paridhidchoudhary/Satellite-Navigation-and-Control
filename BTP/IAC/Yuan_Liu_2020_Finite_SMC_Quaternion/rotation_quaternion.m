function R = rotation_quaternion(q0,qv)

R = (q0*q0 - transpose(qv)*qv)*eye(3) + 2*qv*transpose(qv) - 2*q0*skew(qv);

end