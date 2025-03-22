figure(1)
t = 0:0.001:16;
a = 0.2;
b = -1;
x1 = linspace(-0.4,1,length(t));
x2 = sqrt(x1.^2 - 2.*x1.*a.*cos(b.*t) - a^2);
plot(x1,x2)

c=[14.39999,0.14583333,0.8];
tspan = [0 4];
x0 = [0.2; -0.198];
[t, x] = ode45(@(t, x) model(t,x,a,b,c), tspan, x0);

figure(2)
plot(t, x(:,1))

figure(3)
plot(t, x(:,2))