function [dxdt, u] = singular_interval_model(t, x, data)
dxdt = zeros(size(x));
Pm = data.Pm;
Pe = data.Pe;
C1 = data.C1;
C2 = data.C2;
D = data.D;
dep = data.dep;
H = data.H;
x1 = x(1);
x2 = x(2);
%u
% u = -(C1^2.*D.*x2 - C1^2*Pm + 2.*C1^2.*H.*x2 + 2.*C2^2.*H.*x1)/(C1^2.*Pe.*sin(dep + x1));
u = -(C1^2.*D.*x2 - C1^2*Pm + 2.*C1^2.*H.*x2 + 2.*C2^2.*H.*x1)/(C1^2);

% dxdt(1) = x2;
% dxdt(2) = ((Pm-D*x2) - Pe *sin(x1+dep).*u)/(2*H);
dxdt(1) = x2;
dxdt(2) = ((Pm-D*x2) - u)/(2*H);


end