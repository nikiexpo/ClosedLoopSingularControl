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


[problem,guess]=SIMBS;
options=problem.settings(51);

%% Start simulation
t_step=0.25;
t_end=4;
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
figure
subplot(3,1,1)
hold on
plot(t_sol,x_sol(:,1)/3 )
plot(t_sol,x_sol(:,2)/30 )
legend('x_1','x_2')
xlabel('Time [s]')
ylabel('States')
grid on


subplot(3,1,2)
hold on
plot(t_sol(1:end-1),u_sol(:,1) )
xlabel('Time [s]')
grid on
ylabel('Control Input')



subplot(3,1,3)
hold on
plot(t_sol,abs(x_sol(:,1)/3) )
plot(t_sol,abs(x_sol(:,2)/30) )
xlabel('Time [s]')
ylabel('Tracking Error')
legend('y_1','y_2')
grid on