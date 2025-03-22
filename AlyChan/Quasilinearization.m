

%% ALY CHAN
syms x_1 x_2 x_3 x_4 x_5 
dx = [x_2, -x_1, (1/2)*(x_2)^2 - x_1^2, -x_1, -x_4+x_2];
A_linear = matlabFunction(jacobian(dx, [x_1 x_2 x_3 x_4 x_5]))
%% QUASILINEARIZATION 
n=2;
t_0 = 0; t_f = pi/2;
tspan = [t_0 t_f];
tpoints = linspace(t_0, t_f, 1000)';
xh_t = zeros(1000,5,n);
x0_hi = [0 1 0 0 0.1];
stop_iter = 0;
max_iter = 100;
currentIter = 1;
ode_h_first = @(t,x)A_linear(x(1),x(2))*x;
[th_i,xh_i] = ode45(@(th_i, xh_i)ode_h_first(th_i,xh_i),tspan,x0_hi);
x_iter = @(t, tr, xr) interp1(tr, xr, t);
% ode_h = @(t,x,xg)A(xg(1),xg(2))*x;
error_vec = [];
while (stop_iter==0 && max_iter>currentIter)
    pf = [0 0]';
    ph = zeros(n);
    for i = 1:n
        x0_h = [0 0 0 0 0];
        x0_h(i+3) = 1;
        [th,xh] = ode45(@(th, xh)ode_h(th,xh, th_i, xh_i, A_linear),tspan,x0_h); %% HERE: MAKE A SEP FUNCITON 
        
        ph(:,i) = xh(end,4:5);
        xh_t(:,:,i) = interp1(th, xh(:,:), tpoints);
    end
    % ode_p = @(t,x)A_linear(x(1), x(2))*x + AlyChan(t,x);
    x0_p = [0 1 0 0 0];
    [tp,xp] = ode45(@(tp, xp)ode_p(tp,xp, th_i, xh_i, A_linear),tspan,x0_p);
    xp_t = interp1(tp, xp, tpoints);
    pp_tf = xp(end, 4:5)';
    c =ph\(pf - pp_tf)
    x_tp1 = zeros(1000,5);
    for i=1:n
        x_tp1 = x_tp1 + c(i).*xh_t(:,:,i);
    end
    x_tp1 = x_tp1 + xp_t;
    x_prev = interp1(th_i, xh_i, tpoints);
    error = sum(max(x_tp1 - x_prev));
    if (error < 0.1)
        stop_iter = 1;
    else
        stop_iter = 0;
        currentIter = currentIter + 1;
        xh_i = x_tp1;
        th_i = tpoints;
        error_vec = [error_vec error];
    end
end

t = tpoints;
x = x_tp1; 

optimality_gap = x(end,3) - 0 % zero is analytic min
%%
figure
plot(t,[x(:,1), x(:,2), x(:,3)])
xlabel("Time")
ylabel("States")
xlim(tspan)
ylim([0 1])
legend(["x1", "x2", "x3"])
title('Quasilinearization solution')
grid on
figure
plot(t, [x(:,3), x(:,4)])
xlabel("Time")
ylabel("Co-states")
legend(["p1", "p2"])
ylim([-2, 6])
xlim(tspan)
grid on
%%
function dx = ode_p(t, x, tp, xp, A_linear)
xp_t = interp1(tp, xp, t)';
A = A_linear(xp_t(1), xp_t(2));
dx = A*x - A*xp_t + AlyChan(t,xp_t);
end

function dx = ode_h(t, x, tp, xp, A_linear)
xp_t = interp1(tp, xp, t);
A = A_linear(xp_t(1), xp_t(2));
dx = A*x;
end

function dx = AlyChan(t, x)

    if x(5) < -1e-09
        u = +1;
    elseif x(5) > -1e-09
            u = -1;
    else
        u = -x(1);
    end
    u = -x(1); % for now
    dx = [x(2), u, (1/2)*(x(2)^2 - x(1)^2), -x(1), -x(4)+x(2)]';
end

function res = AlyChanRes(x0, xf)
    res = [x0(1) - 0, x0(2) - 1, x0(3) - 0, xf(1) - 0, xf(2) - 0]';
end