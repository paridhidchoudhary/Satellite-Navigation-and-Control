function Tcs = sat_Tc(Tc,M_max,F_max)
Mc = Tc(1:3,1);
Fc = Tc(4:6,1);

Fcs = zeros(3,1);
Mcs = zeros(3,1);

if norm(Fc)<F_max
    Fcs = Fc;
else
    Fcs = F_max.*sign(Fc);
end

if norm(Mc)<M_max
    Mcs = Mc;
else
    Mcs = M_max.*sign(Mc);
end

Tcs = [Mcs ; Fcs];
end


% if norm(Fc)<F_max
%     Fcs = Fc;
% else
%     Fcs = Fc.*sign(Fc)*(1/(F_max/norm(Fc)));
% end
% 
% if norm(Mc)<M_max
%     Mcs = Mc;
% else
%     Mcs = Mc.*sign(Mc)*(1/(M_max/norm(Mc)));
% end
% 
% Tcs = [Mcs ; Fcs];
% end