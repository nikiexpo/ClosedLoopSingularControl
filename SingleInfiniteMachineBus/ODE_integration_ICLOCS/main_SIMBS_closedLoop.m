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

[problem,guess]=SIMBS;
options=problem.settings(51);

%% Start simulation
t_step=0.25;
t_end=4;
t_update=transpose(0:t_step:t_end);
sim_step=0.01;
x_sol=[];
u_sol=[];
p2_sol = [];
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
    p2_sol = [p2_sol; solution.multipliers.lambda_1toN(:,2)];
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
x1 = x_sol(:,1);
x2 = x_sol(:,2);
u_an = (-H*c2^2)./(Pe.*sin(x1+dep).*x2) .*((x2./c1).^2 .*(1 - (Pm - D)/H) - (x1./c1).^2 + 2/c1^2 .*x1.*x2);
u_an_2 = -(c1^2.*D.*x2 - c1^2*Pm + 2.*c1^2.*H.*x2 + 2.*c2^2.*H.*x1)./(c1^2.*Pe.*sin(dep + x1));
plot(t_sol, u_an_2)
hold on
plot(t_sol(1:end-1),u_sol(:,1) )
hold off
grid on
legend(["analytic", "numerical"])
xlabel('time [s]')
ylabel('control [u]')

figure
plot(t_sol(1:end-1), u_sol - u_an_2(1:end-1), 'LineWidth',2)
hold on
[sw_idx] = find(t_sol >= 0.32,1);
xregion(t_sol(sw_idx), t_sol(end), 'FaceColor','r', 'FaceAlpha',0.25)
grid on
xlabel('time [s]')
ylabel('error')

figure
GLC = -Pe.*sin(x1(sw_idx:end) + dep).*exp(-t_sol(sw_idx:end)); %% SIGN CHANGED
plot(t_sol(sw_idx:end), GLC, 'LineWidth',2)
grid on
xlabel('time [s]')
ylabel('LHS of Generalized Legendre-Clebch')

figure
plot(p2_sol)
% xlabel("time")
ylabel("p_2 multiplier")
grid on
