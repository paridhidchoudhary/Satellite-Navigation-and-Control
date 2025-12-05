% ======================================================================
%            AUXILIARY ROUTINE  (replaces getThrusterPositions)
% ======================================================================
function [B,r_tk] = getThrusterConfig()
% 3×8 positions & 3×8 unit directions (pointing inward toward CoM)
x = [0.75, -0.75];
y = [2.0, -2.0];
z = [0.75, -0.75];

r_tk = [];
for i = 1:2, for j = 1:2, for k = 1:2
        r_tk = [r_tk, [x(i); y(j); z(k)]];
end, end, end

B = -r_tk ./ vecnorm(r_tk);   % unit vectors toward vehicle centre


end
