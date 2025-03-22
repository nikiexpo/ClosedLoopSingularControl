
% parameters
data.H=0.0106;
data.Xt=0.28;
data.Pm=1;
data.Es=1.21;
data.V=1;
data.Pe=data.Es*data.V/data.Pm/data.Xt;
data.D=0.03;
data.dep=asin(1/data.Pe);
data.C1 = 3;
data.C2 = 30;

%boundary value
solinit = bvpinit(linspace(0,4.5,2500), [1.5 15 -0.84 -0.06]);
options = bvpset('Stats','on','RelTol',1e-1);
sol = bvp4c(@(t,x)indirect_model(t,x,data), @boundary_condition, solinit);
t = sol.x';
x = sol.y';
figure
plot(t,[x(:,1), x(:,2)])
xlabel("Time")
ylabel("States")
legend(["x1", "x2"])
figure
plot(t, [x(:,3), x(:,4)])
xlabel("Time")
ylabel("States")
legend(["p1", "p2"])
ylim([-2, 6])
xlim([0 4.5])