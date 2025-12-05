% ======================================================================
% 8 gimballed thrusters on the two 1.5×4 m faces (Malladi et al.)
% Returns:
%   r_tk   3×8   position of thruster k in Fc  (m)
%   j_nom  3×8   nominal inward direction  (unit)   ˆj = –r_tk/‖r_tk‖
% ======================================================================
function [r_tk, j_nom] = getThrusterConfig()

hx = 0.75;                 % half‑lengths of the cuboid [m]
hy = 2.0;
hz = 0.75;

% eight corner points
r_tk = [  hx  hx  hx  hx  -hx -hx -hx -hx ;   % x
          hy  hy -hy -hy   hy  hy -hy -hy ;   % y
          hz -hz  hz -hz   hz -hz  hz -hz ];  % z   (3×8)

j_nom = -r_tk ./ vecnorm(r_tk);              % nominal inward ˆj

end
