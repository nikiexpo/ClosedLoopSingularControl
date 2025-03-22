% Main script to solve the Optimal Control Problem 
%
% Single Infinite Machine Bus Problem Problem in Open-loop
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk


%--------------------------------------------------------

clear all;close all;format compact;

[problem,guess]=SIMBS_OL;          % Fetch the problem definition
problem.settings=@settings_SIMBS_OL_meshrefinement;
options= problem.settings(50);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode45', 0.001 );


%% figure
xx=linspace(solution.T(1,1),solution.tf,10000);


figure
subplot(3,1,1)
hold on
plot(xx,speval(solution,'X',1,xx)/3 )
plot(xx,speval(solution,'X',2,xx)/30 )
% plot(tv,xv(:,1)/3,'k-.' )
% plot(tv,xv(:,2)/30,'k-.' )
xlabel('Time [s]')
ylabel('States')
grid on


subplot(3,1,2)
hold on
plot(xx,speval(solution,'U',1,xx))
% plot(tv,uv(:,1),'k-.' )
plot([solution.T(1,1); solution.tf],[problem.inputs.ul, problem.inputs.ul],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu, problem.inputs.uu],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
grid on
ylabel('Control Input')



subplot(3,1,3)
hold on
plot(xx,abs(speval(solution,'X',1,xx))/3)
plot(xx,abs(speval(solution,'X',2,xx))/30)
% plot(tv,abs(xv(:,1)/3),'k-.' )
% plot(tv,abs(xv(:,2)/30),'k-.' )
xlabel('Time [s]')
ylabel('Tracking Error')
legend('y_1','y_2')
grid on


figure
hold on


[sw_idx] = find(solution.multipliers.lambda_1toN(:,2) > -10e-09,1);
xregion(solution.T(sw_idx), solution.T(end), 'FaceColor','r', 'FaceAlpha',0.25)
xlim([0 4.5])
ylim([-0.06 0.01])
plot(solution.T, solution.multipliers.lambda_1toN(:,2), "LineWidth", 2)
grid on