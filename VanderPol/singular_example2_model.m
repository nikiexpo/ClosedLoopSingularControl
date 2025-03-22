function [dxdt, u] = singular_example2_model(t, x, data)
dxdt = zeros(size(x));
alpha = data.alpha;
beta = data.beta;
c1 = data.c1;
c2 = data.c2;
c3 = data.c3;
x1 = x(1);
x2 = x(2);
%u
u = x1 - alpha*cos(beta*t) + c3*x1 + alpha*beta*cos(beta*t) - c1*x2*(c2-x1);

dxdt(1) = x2;
dxdt(2) = c1*x2*(c2 - x1) - c3*x1 + u;


end