function [g_dot, epsilon_v]= kin_leader_follower(r,v,omega,b)
g=[r b;
0 1];
epsilon=[omega;v];
epsilon_v=[cross_pdt(omega) v; %equation10b
0 0];
g_dot=g*epsilon_v;    %equation10a
end