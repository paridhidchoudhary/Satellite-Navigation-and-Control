function Rproj = local_projectR(R)
    c1 = R(:,1); c1 = c1/norm(c1);
    c2 = R(:,2) - c1*(c1.'*R(:,2)); c2 = c2/norm(c2);
    c3 = cross(c1, c2);
    Rproj = [c1, c2, c3];
end