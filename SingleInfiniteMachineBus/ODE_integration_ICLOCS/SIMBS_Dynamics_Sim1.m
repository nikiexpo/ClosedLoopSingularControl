function dx = SIMBS_Dynamics_Sim1(x,t,vdat)
%Single Infinite Machine Bus Dynamics for Simulation
%
% Syntax:  
%          [dx] = Dynamics(x,u,p,t,vdat)
% 
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%    vdat - structured variable containing the values of additional data used inside
%          the function%      
% Output:
%    dx - time derivative of x
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk
%
%------------- BEGIN CODE --------------


x1_sol=speval(vdat.solution,'X',1,t);
x2_sol=speval(vdat.solution,'X',2,t);

u = -(vdat.C1^2.*vdat.D.*x2_sol - vdat.C1^2*vdat.Pm + 2.*vdat.C1^2.*vdat.H.*x2_sol + 2.*vdat.C2^2.*vdat.H.*x1_sol)/(vdat.C1^2.*vdat.Pe.*sin(vdat.dep + x1_sol));

x1 = x(1);x2 = x(2);
dx = [x2,(vdat.Pm-vdat.D*x2)/2/vdat.H-vdat.Pe/2/vdat.H*sin(x1+vdat.dep).*u];


%------------- END OF CODE --------------