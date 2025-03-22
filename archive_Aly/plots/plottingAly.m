load("aly_OL_sol.mat")
load("aly_RH_sol.mat")
load("sol_Aly_CLDC_16_03_25.mat")
load("sol_Aly_OLDC_16_03_25.mat")
load("sol_Aly_OLDCMR_16_03_25.mat")

%OL
xx=linspace(solutionOL.T(1,1),solutionOL.T(end,1),5000);
figure
hold on
plot(xx,speval(solutionOL,'X',1,xx),'b-' )
plot(xx,speval(solutionOL,'X',2,xx),'r-' )

xlim([0 solutionOL.tf])
xlabel('Time [s]')
ylabel('States')
legend('x1 [-]','x2 [-]')
grid on

% figure
% subplot(2,1,1)
% plot(xx,speval(solutionOL,'U',1,xx),'b-' )
% hold on
% plot([solutionOL.T(1,1); solutionOL.tf],[problem.inputs.ul, problem.inputs.ul],'r-' )
% plot([solutionOL.T(1,1); solutionOL.tf],[problem.inputs.uu, problem.inputs.uu],'r-' )
% xlim([0 solutionOL.tf])
% ylim([-1 1])
% xlabel('Time [s]')
% grid on
% ylabel('Control Input')
% legend({"$u_{Open Loop}$"}, "Interpreter","latex", "FontSize", 12)
% hold off
% 
% subplot(2,1,2)
% plot(solutionCL.t_sol(1:end-1,:), solutionCL.u_sol(1:end, :), 'b-' )
% hold on
% plot([solutionOL.T(1,1); solutionOL.tf],[problem.inputs.ul, problem.inputs.ul],'r-' )
% plot([solutionOL.T(1,1); solutionOL.tf],[problem.inputs.uu, problem.inputs.uu],'r-' )
% xlim([0 solutionOL.tf])
% ylim([-1 1])
% xlabel('Time [s]')
% grid on
% hold off
% ylabel('Control Input')
% legend({"$u_{Receding Horizon}$"}, "Interpreter","latex", "FontSize", 12)

%% figure for paper
f= figure(4);
f.Position = [100, 100, 900, 500];
a1 = axes('FontSize', 14, 'NextPlot', 'add');
a1.GridLineStyle = '--';
a1.GridColor = 'k';
a1.GridAlpha = 0.9 ;
hold on
plot(xx,speval(solutionOLDC,'U',1,xx),'b-',LineWidth=1.5 )
plot(xx,speval(solutionOLDCMR,'U',1,xx),'g-',LineWidth=1.5 )
plot(xx,speval(solutionOL,'U',1,xx),'k-', LineWidth=2)
plot([solutionOL.T(1,1); solutionOL.tf],[-1, -1],'r-' )
plot([solutionOL.T(1,1); solutionOL.tf],[1, 1],'r-' )
xlim([0 solutionOL.tf])
ylim([-1 1])
xlabel('Time [s]',"Interpreter","latex")
grid on
ylabel('Control Input', "Interpreter","latex")
legend({"DC", "DC with MR", "IRR-DC with MR", "Bounds"}, "Interpreter","latex", "FontSize", 14)
hold off
exportgraphics(f, "ModAly_OL.eps", "Resolution",600);

f = figure(3);
f.Position = [100, 100, 900, 500];
a1 = axes('FontSize', 14, 'NextPlot', 'add');
a1.GridLineStyle = '--';
a1.GridColor = 'k';
a1.GridAlpha = 0.9 ;
hold on
plot(solutionCLDC.t_sol(1:end-1,:), solutionCLDC.u_sol(1:end, :), 'b-',LineWidth=1.5 )
plot(solutionCL.t_sol(1:end-1,:), solutionCL.u_sol(1:end, :), 'k-',LineWidth=2 )
plot([solutionOL.T(1,1); solutionOL.tf],[-1, -1],'r-' )
plot([solutionOL.T(1,1); solutionOL.tf],[1, 1],'r-' )
xlim([0 solutionOL.tf])
ylim([-1 1])
xlabel('Time [s]',"Interpreter","latex")
grid on
ylabel('Control Input',"Interpreter","latex")
legend({"DC", "IRR-DC", "Bounds"}, "Interpreter","latex", "FontSize", 14)
hold off
exportgraphics(f, "ModAly_CL.eps", "Resolution",600);



