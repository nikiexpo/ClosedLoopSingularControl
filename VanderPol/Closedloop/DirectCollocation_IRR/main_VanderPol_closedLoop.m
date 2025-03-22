% Main script to solve the Optimal Control Problem 
%
% Van der Pol Problem in Closed-loop
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


[problem,guess]=VanderPol;
options=problem.settings(51);

%% Start simulation
t_step=0.5;
t_end=16;
t_update=transpose(0:t_step:t_end);
sim_step=0.01;
x_sol=[];
u_sol=[];
t_sol=0;
comp_time=[];
for i=1:t_end/t_step

    
    if i==1
        [solution,MRHistory,OCP]=solveMyProblem(problem,guess,options);
    else
        vdat.ini_t=t_update(i);
        [OCP] = updateMyProblem(OCP,'z0',solution.z,'x0',xv(end,:),'userdata',vdat); % look into the OCP structure
        [solution,MRHistory,OCP]=solveMyProblem(problem,guess,options,OCP);
    end
    
    [ tv, xv, uv ] = simulateSolutionSegment( problem, solution, 'ode113', 0:sim_step:t_step );
    t_sol=[t_sol;t_sol(end)+tv(2:end)];
    u_sol=[u_sol;uv(1:end-1,:)];
    x_sol=[x_sol;xv(1:end-1,:)];
    comp_time=[comp_time;MRHistory.cpu];

end
x_sol(end+1,:)=xv(end,:);
%%
solutionCL_IRR.t = t_sol;
solutionCL_IRR.x = x_sol;
solutionCL_IRR.u = u_sol;

figure
subplot(4,1,1)
hold on
plot(t_sol,x_sol(:,1) )
plot(t_sol,0.2*cos(-t_sol),'k--')
legend('x_1','x_{1_{ref}}')
xlabel('Time [s]')
ylabel('States')
ylim([-1 1])
grid on

subplot(4,1,2)
hold on
plot(t_sol,x_sol(:,2) )
plot(t_sol,0.2*sin(-t_sol),'k-.')
legend('x_2','x_{2_{ref}}')
xlabel('Time [s]')
ylabel('States')
ylim([-1 1])
grid on


subplot(4,1,3)
hold on
plot(t_sol(1:end-1),u_sol(:,1) )
xlabel('Time [s]')
grid on
ylabel('Control Input')


subplot(4,1,4)
hold on
plot(t_sol,abs(x_sol(:,1)-0.2*cos(-t_sol)) )
plot(t_sol,abs(x_sol(:,2)-0.2*sin(-t_sol)) )
xlabel('Time [s]')
ylabel('Tracking Error')
legend('y_1','y_2')
grid on

figure
hold on
plot(x_sol(:,1),x_sol(:,2) )
plot(0.2*cos(-t_sol),0.2*sin(-t_sol),'k--')
legend('Controlled','Reference')
xlabel('Time [s]')
ylabel('States')
grid on