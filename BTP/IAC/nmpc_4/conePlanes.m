function [A,b] = conePlanes(j_nom)
% j_nom : 3×1 nominal inward unit vector
% Build four planes delimiting a π/4 rad circular cone around j_nom

cone_ang = pi/4;
ax = null(j_nom.');        % two orthonormal vectors orthogonal to j_nom
dirs = [ ax(:,1)   ax(:,2)   -ax(:,1)  -ax(:,2) ];   % 4 radial dirs

A = zeros(4,3);  b = zeros(4,1);
for p = 1:4
    n = (j_nom*cos(cone_ang) + dirs(:,p)*sin(cone_ang));  % plane normal
    A(p,:) = n.';                                         % row
    b(p)   = 0;                                           % origin side
end
end
