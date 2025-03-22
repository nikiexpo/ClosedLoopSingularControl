% Main script to solve the Optimal Control Problem 
%
% Single Infinite Machine Bus Problem in Closed-loop
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk


clear all
close all
clc

%%

[problem,guess]=AlyChan;
options=problem.settings(25);

%% Start simulation
% t_step=0.1;
% t_end=pi/2;
% t_update=transpose(0:t_step:t_end);
% sim_step=0.01;
% x_sol=[];
% u_sol=[];
% p2_sol = [];
% p3_sol = [];
% t_sol=0;
% comp_time=[];
% tpoints=linspace(0,t_end,1000)';
% for i=1:t_end/t_step
% 
% 
%     if i==1
%         [solution,MRHistory,OCP]=solveMyProblem(problem,guess,options);
%         [ tv, xv, uv ] = simulateSolutionSegment( problem, solution, 'ode45', 0:sim_step:t_end );
% 
%         t_sol=tpoints;
%         u_sol=interp1(tv,uv,tpoints);
%         x_sol=interp1(tv,xv,tpoints);
%     else
%         vdat.ini_t=t_update(i);
%         problem.time.t0_min = t_update(i);
%         problem.time.t0_max = t_update(i);
%         % solution.z(1)=t_update(i);
%         curr_time = find(tv >= (i-1)*t_step,1);
%         [OCP] = updateMyProblem(OCP,'z0',solution.z,'x0',xv(curr_time,:),'userdata',vdat); % look into the OCP structure
%         [solution,MRHistory,OCP]=solveMyProblem(problem,guess,options,OCP);
%         [ tv, xv, uv ] = simulateSolutionSegment( problem, solution, 'ode45', (i-1)*t_step:sim_step:t_end );
%         tpoints_cl = tpoints(curr_time:end);
%         u_sol(curr_time:end,:)=interp1(tv,uv,tpoints_cl);
%         x_sol(curr_time:end,:)=interp1(tv,xv,tpoints_cl);
%     end
    
    % [ tv, xv, uv ] = simulateSolutionSegment( problem, solution, 'ode113', 0:sim_step:t_end );
    % t_sol=[t_sol;t_sol(end)+tv(2:end)];
    % u_sol=[u_sol;uv(1:end-1,:)];
    % x_sol=[x_sol;xv(1:end-1,:)];
    % comp_time=[comp_time;MRHistory.cpu];
    % p2_sol = [p2_sol; solution.multipliers.lambda_1toN(:,2)];
    % p3_sol = [p3_sol; solution.multipliers.lambda_1toN(:,3)];
% end
% x_sol(end+1,:)=xv(end,:);

%%
t_end = pi/2;
t_step = 0.05;
t_update=transpose(0:t_step:t_end);
sim_step = 0.01;
ini_t = 0;
[solution,MRHistory,OCP]=solveMyProblem(problem,guess,options);
[ tv, xv, uv ] = simulateSolutionSegment( problem, solution, 'ode45', ini_t:sim_step:t_end );
%%
% t_sol = solution.T;
% x_sol = solution.X;
% u_sol = solution.U;\
t_sol = [0:t_step:t_end]';
x_sol = interp1(tv,xv,t_sol);
u_sol = interp1(tv,uv,t_sol);
figure
subplot(2,1,1)
hold on
plot(t_sol,x_sol(:,1))
plot(t_sol,x_sol(:,2))
plot(t_sol,x_sol(:,3))
% plot(tv,xv(:,1), "--b")
% plot(tv,xv(:,2), "--b")
% plot(tv,xv(:,3), "--b")
legend('x_1','x_2', 'x_3')
xlabel('Time')
ylabel('States')
grid on
hold off

subplot(2,1,2)
hold on
plot(t_sol,u_sol )
% plot(tv,uv, "--b")
xlabel('Time')
grid on
ylabel('Control Input')

ol_sol = solution;
optimality_gap_ol = x_sol(end,3)
opt_traj_gap_ol = -x_sol(:,1) - u_sol;

%% 
while ini_t < t_end - t_step
    ini_t = ini_t + t_step;
    problem.time.t0_max = ini_t;
    problem.time.t0_min = ini_t;
    indx_curT = find(solution.T >= ini_t,1);
    interp1(solution.T,solution.X,ini_t);
    problem.states.x0 = interp1(solution.T,solution.X,ini_t);
    problem.states.x0l = interp1(solution.T,solution.X,ini_t);
    problem.states.x0u = interp1(solution.T,solution.X,ini_t);
    problem.inputs.u0l = interp1(solution.T,solution.U,ini_t);
    problem.inputs.u0u = interp1(solution.T,solution.U,ini_t);
    
    guess.inputs = [problem.inputs.u0l; 0];
    guess.states = [problem.states.x0 ; solution.X(end,:)];
    guess.t0 = ini_t;
    
    [solution,MRHistory,OCP]=solveMyProblem(problem,guess,options);
    [ tv, xv, uv ] = simulateSolutionSegment( problem, solution, 'ode45', ini_t:sim_step:t_end );
    
    idx = find(t_sol == ini_t,1);
    x_sol(idx:end,:) = interp1(solution.T, solution.X, t_sol(idx:end));
    u_sol(idx:end,:) = interp1(solution.T, solution.U, t_sol(idx:end));
end


%% 
figure
subplot(2,1,1)
hold on
plot(t_sol,x_sol(:,1))
plot(t_sol,x_sol(:,2))
plot(t_sol,x_sol(:,3))
% plot(tv,xv(:,1), "--b")
% plot(tv,xv(:,2), "--b")
% plot(tv,xv(:,3), "--b")
legend('x_1','x_2', 'x_3')
xlabel('Time [s]')
ylabel('States')
grid on
xlim([0 t_end])

subplot(2,1,2)
hold on
plot(t_sol,u_sol )
% plot(tv,uv, "--b")
xlabel('Time [s]')
grid on
ylabel('Control Input')
xlim([0 t_end])

cost = x_sol(end,3)
acost = 0
optimality_gap_cl = cost

figure
opt_traj_gap_cl = -x_sol(:,1) - u_sol;
hold on
plot(t_sol, opt_traj_gap_cl)
plot(t_sol, opt_traj_gap_ol)
yline(0,'--r')
legend(["Opt gap CL", "Opt gap OL"])
grid on

