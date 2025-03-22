% REALLY ACCURATE ANALYTICAL SOLUTION TO TPBVP
%% BVP4C MIGHT NOT BE POSSIBLE FOR SINGULAR INTERVALS
% solinit = bvpinit(linspace(0,pi/2,2500), [0 1 0 0 0.1]);
% options = bvpset('Stats','on','RelTol',1e-1);
% sol = bvp4c(@(t,x)AlyChan(t,x), @AlyChanRes, solinit);
% t = sol.x';
% x = sol.y';

%% For SMIB sys
% clear; 
% syms t
% syms x_1(t) x_2(t)  u(t)
% syms Pe H Pm D c_1 c_2 dep
% syms p_1(t) p_2(t)  
% x_dot_1 = x_2;
% x_dot_2 = (Pm - D*x_2)/(2*H) - Pe/(2*H)*sin(x_1 + dep)*u;
% % L = exp(-t)* ((x_1/c_1)^2 + (x_2/c_2)^2)
% L = ((x_1/c_1)^2 + (x_2/c_2)^2);
% Hamiltonian = p_1*x_dot_1 + p_2*x_dot_2 + L;
% dx1_dt = diff(x_1, t);
% dx2_dt = diff(x_2,t);
% dp1_dt = diff(p_1,t);
% ddx1_dt2 = diff(dx1_dt, t);
% p_dot_1 = -diff(Hamiltonian, sym(x_1)) ;
% p_dot_1 = subs(p_dot_1, p_2, 0);
% 
% p_dot_2 = -diff(Hamiltonian,sym(x_2));
% %subs p2 = 0
% p_dot_2 = subs(p_dot_2, p_2, 0) ;
% p_ddot_2 = diff(p_dot_2, t) ;
% p_ddot_2 = subs(p_ddot_2, {dp1_dt, dx2_dt}, {p_dot_1, x_dot_2});
% 
% u_solve = isolate(p_ddot_2 == 0, u);
% u_sol = simplify(rhs(u_solve));
% 
% % parameters
% data.H=0.0106;
% data.Xt=0.28;
% data.Pm=1;
% data.Es=1.21;
% data.V=1;
% data.Pe=data.Es*data.V/data.Pm/data.Xt;
% data.D=0.03;
% data.dep=asin(1/data.Pe);
% data.C1 = 3;
% data.C2 = 30;
% 
% % getting the linearization reduced ODE
% x_dot_2 = subs(x_dot_2, u, u_sol);
% p_dot_1 = -diff(Hamiltonian, sym(x_1)) ;
% p_dot_1 = simplify(subs(p_dot_1, u, u_sol));
% p_dot_2 = -diff(Hamiltonian,sym(x_2));
% syms y_1 y_2 y_3 y_4
% f_rODE = [x_dot_1, simplify(x_dot_2), p_dot_1, p_dot_2 ];
% f_rODE = subs(f_rODE, [c_1, c_2, H, Pm, Pe, dep, D], [data.C1, data.C2, data.H, data.Pm, data.Pe, data.dep, data.D]);
% f_rODE = subs(f_rODE, [x_1, x_2, p_1, p_2], [y_1, y_2, y_3, y_4]);
% A_linear = matlabFunction(jacobian(f_rODE, [y_1, y_2, y_3, y_4]))

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
legend(["x1", "x2", "x3"])
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