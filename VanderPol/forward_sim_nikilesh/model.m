function dx = model(t, x, a, b, c)
dx = zeros(2,1);
u = x(1)*(1+c(1)*x(2)+c(3))-c(1)*c(2)*x(2)-a*cos(b*t)+x(1)/(x(2)+0.01) *a*sin(b*t)
dx(1) = x(2);
dx(2) = c(1)*x(2)*(c(2)-x(1))-c(3)*x(1)+u;
end