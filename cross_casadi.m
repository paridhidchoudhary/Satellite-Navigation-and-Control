function c = cross_casadi(a,b)
    c = [a(2)*b(3)-a(3)*b(2);
         a(3)*b(1)-a(1)*b(3);
         a(1)*b(2)-a(2)*b(1)];
end



% function out = cross_casadi(a, b)
%     % if isa(a, 'casadi.SX') && isa(b, 'casadi.MX')
%     %     a = casadi.MX(a);  % promote to MX
%     % elseif isa(a, 'casadi.MX') && isa(b, 'casadi.SX')
%     %     b = casadi.MX(b);  % promote to MX
%     % end
%     % out = cross(a, b);
% end
