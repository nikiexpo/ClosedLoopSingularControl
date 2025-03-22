% Main script to solve the Optimal Control Problem 
%
% Van der Pol Problem in Open-loop
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

[problem,guess]=VanderPol_OL;          % Fetch the problem definition
options= problem.settings(51);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);


%% figure
t_sol=linspace(solution.T(1,1),solution.tf,5000)';
x_sol(:,1)=speval(solution,'X',1,t_sol);
x_sol(:,2)=speval(solution,'X',2,t_sol);
u_sol(:,1)=speval(solution,'U',1,t_sol);

% figure
% subplot(3,1,1)
% hold on
% plot(xx,speval(solution,'X',1,xx) )
% plot(xx,speval(solution,'X',2,xx) )
% xlabel('Time [s]')
% ylabel('States')
% grid on
% 
% 
% subplot(3,1,2)
% hold on
% plot(xx,speval(solution,'U',1,xx))
% xlim([0 solution.tf])
% xlabel('Time [s]')
% grid on
% ylabel('Control Input')
% 
% 
% 
% subplot(3,1,3)
% hold on
% plot(xx,abs(speval(solution,'X',1,xx)-0.2*cos(-xx')))
% plot(xx,abs(speval(solution,'X',2,xx)-0.2*sin(-xx')))
% xlabel('Time [s]')
% ylabel('Tracking Error')
% legend('y_1','y_2')
% grid on


figure
subplot(4,1,1)
hold on
plot(t_sol,x_sol(:,1) )
plot(t_sol,0.2*cos(-t_sol),'k--')
legend('x_1','x_{1_{ref}}')
xlabel('Time [s]')
ylabel('States')
ylim([-1 1])
xlim([0 solution.tf])
grid on

subplot(4,1,2)
hold on
plot(t_sol,x_sol(:,2) )
plot(t_sol,0.2*sin(-t_sol),'k-.')
legend('x_2','x_{2_{ref}}')
xlabel('Time [s]')
ylabel('States')
ylim([-1 1])
xlim([0 solution.tf])
grid on


subplot(4,1,3)
hold on
plot(t_sol,u_sol(:,1) )
xlabel('Time [s]')
grid on
ylabel('Control Input')
xlim([0 solution.tf])

subplot(4,1,4)
hold on
plot(t_sol,abs(x_sol(:,1)-0.2*cos(-t_sol)) )
plot(t_sol,abs(x_sol(:,2)-0.2*sin(-t_sol)) )
xlabel('Time [s]')
ylabel('Tracking Error')
legend('y_1','y_2')
xlim([0 solution.tf])
grid on

figure
hold on
plot(x_sol(:,1),x_sol(:,2) )
plot(0.2*cos(-t_sol),0.2*sin(-t_sol),'k--')
legend('Controlled','Reference')
xlabel('Time [s]')
ylabel('States')
grid on