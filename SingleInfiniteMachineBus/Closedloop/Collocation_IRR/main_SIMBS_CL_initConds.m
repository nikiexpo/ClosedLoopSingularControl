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




%%
conds = [1.5 15; 2 20; -1 5; -1.5 -15; 3 -3; -0.25 1; 0 5; 0 -5; -3 13; -1 17];
[iter, ~] = size(conds);
sol.t_sol = {};
sol.u_sol = {};
sol.x_sol = {};

for j = 1:iter
[problem,guess]=SIMBS_initConds(conds(j,:));
options=problem.settings(51);

%% Start simulation
t_step=0.2;
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
    tp = t_sol(end);
    t_sol=[t_sol;t_sol(end)+tv(2:end)];
    u_sol=[u_sol;uv(1:end-1,:)];
    x_sol=[x_sol;xv(1:end-1,:)];
    comp_time=[comp_time;MRHistory.cpu];
    % p2_sol = [p2_sol; interp1(solution.T,solution.multipliers.lambda_1toN(:,2),tp+tv(2:end))];
end
x_sol(end+1,:)=xv(end,:);

sol.t_sol(j) = {t_sol};
sol.u_sol(j) = {u_sol};
sol.x_sol(j) = {x_sol};
end