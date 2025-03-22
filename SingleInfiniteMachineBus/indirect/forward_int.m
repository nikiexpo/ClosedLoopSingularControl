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

%ode integration [Open Loop]
tspan = [0.32 0.32+0.5];
x0 = [1.3049  -12.4312]'; 

[t,x] = ode78(@(t, x) singular_interval_model(t, x, data), tspan, x0);


for k=1:numel(t)
    [~,control(k,:)] = singular_interval_model(t(k), x(k,:), data);
end
figure
plot(t,control)
figure
plot(t, x(:,1))
% closed 
% t_step=0.25;
% t_end=4;
% t_update=transpose(0:t_step:t_end);
% sim_step=0.01;
% x_sol=[];
% u_sol=[];
% t_sol=0;
% comp_time=[];
% tspanCL = [0 0.5];
% for i=1:t_end/t_step
%     [ tv, xv, uv ] = simulateSolutionSegment( problem, solution, 'ode113', 0:sim_step:t_step );
%     [ tv, xv, uv ] = ode113(@(tv, xv) singular_interval_model(tv, xv,
%     data), tspanCL, ) %% CONTINUE HERE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%     for k=1:numel(tv)
%         [~,control(k,:)] = singular_interval_model(tv(k), xv(k,:), data);
%     end
%     t_sol=[t_sol;t_sol(end)+tv(2:end)];
%     u_sol=[u_sol;uv(1:end-1,:)];
%     x_sol=[x_sol;xv(1:end-1,:)];
% end