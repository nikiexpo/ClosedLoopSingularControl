function dxdt = indirect_model(t, x, data)
dxdt = zeros([1 4]);
Pm = data.Pm;
Pe = data.Pe;
C1 = data.C1;
C2 = data.C2;
D = data.D;
dep = data.dep;
H = data.H;
%u
threshold = 10^-9;
if x(4) < -threshold
    u = +1;
elseif x(4) > +threshold
    u = -1;
else
    u = -(C1^2.*D.*x(2) - C1^2*Pm + 2.*C1^2.*H.*x(2) + 2.*C2^2.*H.*x(1))./(C1^2);
end

dxdt(1) = x(2);
dxdt(2) = (Pm - D*x(2))/(2*H) - (u/(2*H));
dxdt(3) = ((1/tan(dep + x(1))) *x(4)*u)/(2*H) - (2*x(1))/C1^2;
dxdt(4) = -x(3) - (2*x(2))/C2^2 + (D*x(4))/(2*H);
end