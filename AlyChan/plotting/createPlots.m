load("CL_sol_11_03_25.mat")
load("OL_sol_15_03_25.mat")
load("QL_sol_15_03_25.mat")
load("OL_DC_sol_15_03_25.mat")
load("sol_AlyChanCLDC_21_03_25.mat")
figure(1)
hold on
p1 = plot(solutionCL(:,1), solutionCL(:,2:end-1),"^", Color=[0.25 0.5 0.25], LineWidth=2);
p2 = plot(solution.T, solution.X(:,1:3), "-r");
p3 = plot(solutionQuasi(:,1), solutionQuasi(:,2:end-1), "--k", LineWidth=2);
p4 = plot(solutionDC(:,1), solutionDC(:,2:4), "-", Color = [0.1 0 0.9 0.4], LineWidth=0.5);
hold off
legend([p3(1) p4(1) p2(1) p1(1)  ], { "Analytical", "Open loop DC", "Open loop IRR-DC", "Receding horizon IRR-DC" })
grid on
xlim([0 pi/2])
ylim([0 1])
figure(2)
hold on
p4 = plot(solutionDC(:,1), solutionDC(:,end), "-", Color = [0.1 0 0.9 0.4], LineWidth=0.5);
p1 = plot(solutionCL(:,1), solutionCL(:,end), "^", Color=[0.25 0.5 0.25], LineWidth=2);
p2 = plot(solution.T, solution.U, "-r");
p3 = plot(solutionQuasi(:,1), -solutionQuasi(:,2), "--k", LineWidth=2);

hold off
% legend([ p2(1) p3(1) p4(1)], { "Open loop IRR-DC", "Analytical", "Open loop DC"})
legend([p3(1) p4(1) p2(1) p1(1)  ], { "Analytical", "Open loop DC", "Open loop IRR-DC", "Receding horizon IRR-DC" })
grid on
xlim([0 pi/2])
ylim([-1 0])


%% trying tiled layout
f = figure(3);
f.Position = [100, 100, 900, 800];

t = tiledlayout(3,1);
t.TileSpacing = 'tight';
a1 = nexttile;
a1.XAxis.FontSize = 16;
a1.YAxis.FontSize = 16;
a1.GridLineStyle = '--';
a1.GridColor = 'k';
a1.GridAlpha = 0.9 ;
hold on
p1 = plot(solutionCL(:,1), solutionCL(:,2:end-1),"m^", LineWidth=2, MarkerSize=10);
p2 = plot(solution.T, solution.X(:,1:3), "-r");
p3 = plot(solutionQuasi(:,1), solutionQuasi(:,2:end-1), "--k", LineWidth=2);
p4 = plot(solutionDC(:,1), solutionDC(:,2:4), "-", Color = [0.1 0 0.9 0.4], LineWidth=0.5);
p5 = plot(solutionCLDC.t_sol, solutionCLDC.x_sol, "*", Color=[0.25 0.5 0.25], LineWidth=1.5, MarkerSize=10);
hold off
xlabel("Time[s]", "Interpreter", "latex")
ylabel("States", "Interpreter", "latex")
% legend([p3(1) p4(1) p2(1) p1(1)  ], { "Analytical", "Open loop DC", "Open loop IRR-DC", "Receding horizon IRR-DC" })
grid on
xlim([0 pi/2])
ylim([0 1])

a2 = nexttile([2 1]);
a2.XAxisLocation = "top";
a2.XAxis.FontSize = 16;
a2.YAxis.FontSize = 16;
a2.GridLineStyle = '--';
a2.GridColor = 'k';
a2.GridAlpha = 0.9 ;
hold on
p4 = plot(solutionDC(:,1), solutionDC(:,end), "-", Color = [0.1 0 0.9 0.8], LineWidth=0.5);
p1 = plot(solutionCL(:,1), solutionCL(:,end), "m^", LineWidth=2, MarkerSize=10);
p2 = plot(solution.T, solution.U, "-r");
p3 = plot(solutionQuasi(:,1), -solutionQuasi(:,2), "--k", LineWidth=2);
p5 = plot(solutionCLDC.t_sol, solutionCLDC.u_sol, "*", Color=[0.25 0.5 0.25], LineWidth=1.5, MarkerSize=10);
hold off
ylabel("Control", "Interpreter", "latex")
% legend([ p2(1) p3(1) p4(1)], { "Open loop IRR-DC", "Analytical", "Open loop DC"})
lgd = legend([p3(1) p4(1) p2(1) p5(1) p1(1) ], { "Analytical", "Open loop DC", "Open loop IRR-DC", "Shrinking-horizon DC","Shrinking-horizon IRR-DC" }, "Interpreter", "latex");
grid on
xlim([0 pi/2])
ylim([-1 0])

lgd.Layout.Tile = 'south';
lgd.FontSize = 18;
lgd.NumColumns = 2;

exportgraphics(f,"AlyChan.eps", "Resolution",600)