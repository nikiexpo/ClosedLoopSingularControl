
solinit = bvpinit(linspace(0,pi/2,2500), [0 1 0 0 0.1]);
options = bvpset('Stats','on','RelTol',1e-1);
sol = bvp4c(@(t,x)AlyChan_Dynamics_Singular(t,x), @AlyChan_bounds, solinit);
t = sol.x';
x = sol.y';
figure
plot(t,[x(:,1), x(:,2), x(:,3)])
grid on
legend(["x1", "x2", "x3"])
xlabel("Time")
ylabel("States")
title('bvp4c solution')
% figure
% plot(t, [x(:,3), x(:,4)])
optimality_gap = x(end,3)

function [dx] = AlyChan_Dynamics_Singular(t,x)

x1 = x(1);x2 = x(2);x4=x(4);

dx(:,1) = x2;
dx(:,2) = -x1;
dx(:,3) = 0.5*(x2.^2-x1.^2);
dx(:,4) = -x1;
dx(:,5) = -x4 + x2;

end

function res = AlyChan_bounds(x0, xf)
res = [x0(1) - 0;
        x0(2) - 1;
        x0(3) - 0;
        xf(4) - 0;
        xf(5) - 0];
end