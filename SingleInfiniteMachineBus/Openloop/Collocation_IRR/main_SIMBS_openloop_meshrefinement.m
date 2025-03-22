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
plot(tv,xv(:,1)/3,'k-.' )
plot(tv,xv(:,2)/30,'k-.' )
xlabel('Time [s]')
ylabel('States')
grid on


subplot(3,1,2)
hold on
plot(xx,speval(solution,'U',1,xx))
plot(tv,uv(:,1),'k-.' )
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
plot(tv,abs(xv(:,1)/3),'k-.' )
plot(tv,abs(xv(:,2)/30),'k-.' )
xlabel('Time [s]')
ylabel('Tracking Error')
legend('y_1','y_2')
grid on

figure
H=0.0106;
Xt=0.28;
Pm=1;
Es=1.21;
V=1;
Pe=Es*V/Pm/Xt;
D=0.03;
dep=asin(1/Pe);
c1 = 3;
c2=30;
x1 = speval(solution,'X',1,xx);
x2 = speval(solution,'X',2,xx);
u_an = (-H*c2^2)./(Pe.*sin(x1+dep).*x2) .*((x2./c1).^2 .*(1 - (Pm - D)/H) - (x1./c1).^2 + 2/c1^2 .*x1.*x2);
u_an_2 = -(c1^2.*D.*x2 - c1^2*Pm + 2.*c1^2.*H.*x2 + 2.*c2^2.*H.*x1)./(c1^2.*Pe.*sin(dep + x1));
plot(xx, u_an_2)
hold on
plot(xx,speval(solution,'U',1,xx))
hold off
grid on
legend(["analytic", "numerical"])


figure
plot(xx, speval(solution,'U',1,xx) - u_an_2, 'LineWidth',2)
hold on
[sw_idx] = find(xx >= 0.32,1);
xregion(xx(sw_idx), xx(end), 'FaceColor','r', 'FaceAlpha',0.25)
grid on
xlabel('time [s]')
ylabel('error')
hold off

figure
GLC = -Pe.*sin(x1(sw_idx:end)' + dep).*exp(-xx(sw_idx:end)); %% SIGN CHANGED
plot(xx(sw_idx:end), GLC, 'LineWidth',2)
grid on
xlabel('time [s]')
ylabel('LHS of Generalized Legendre-Clebch')

