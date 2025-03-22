load("solutionOL_DC.mat")
load("solutionOL_IRR.mat")
load("solutionOL_IRR_MR.mat")
load("solutionCL_IRR.mat")

%% OL DC
f1 = figure;
f1.Position = [100, 100, 900, 400];
a1 = axes('FontSize', 14, 'NextPlot', 'add');

hold on
plot(solutionOL_DC.t,solutionOL_DC.u(:,1), "LineWidth",2)
xlabel('Time [s]', Interpreter='latex', FontSize=16)
grid on
a1.GridLineStyle = '--';
a1.GridColor = 'k';
a1.GridAlpha = 1 ;% maximum line opacity;
ylabel('Control Input', Interpreter='latex', FontSize=16)
xlim([0 16])

exportgraphics(f1,"vanderOL_CD.eps", "Resolution",600)

%%
f1 = figure;
f1.Position = [100, 100, 900, 400];
a1 = axes('FontSize', 14, 'NextPlot', 'add');

hold on
plot(solutionOL_IRR.t,solutionOL_IRR.u(:,1), "LineWidth",2)
xlabel('Time [s]', Interpreter='latex', FontSize=16)
grid on
a1.GridLineStyle = '--';
a1.GridColor = 'k';
a1.GridAlpha = 1 ;% maximum line opacity;
ylabel('Control Input', Interpreter='latex', FontSize=16)
xlim([0 16])

exportgraphics(f1,"vanderOL_IRR.eps", "Resolution",600)

%%
f1 = figure;
f1.Position = [100, 100, 900, 400];
a1 = axes('FontSize', 14, 'NextPlot', 'add');

hold on
plot(solutionOL_IRR_MR.t,solutionOL_IRR_MR.u(:,1), "LineWidth",2)
xlabel('Time [s]', Interpreter='latex', FontSize=16)
grid on
a1.GridLineStyle = '--';
a1.GridColor = 'k';
a1.GridAlpha = 1 ;% maximum line opacity;
ylabel('Control Input', Interpreter='latex', FontSize=16)
xlim([0 16])

exportgraphics(f1,"vanderOL_IRR_MR.eps", "Resolution",600)

%%

f1 = figure;
f1.Position = [100, 100, 900, 400];
a1 = axes('FontSize', 14, 'NextPlot', 'add');

hold on
plot(solutionCL_IRR.t(1:end-1,:),solutionCL_IRR.u(:,1), "LineWidth",2)
xlabel('Time [s]', Interpreter='latex', FontSize=16)
grid on
a1.GridLineStyle = '--';
a1.GridColor = 'k';
a1.GridAlpha = 1 ;% maximum line opacity;
ylabel('Control Input', Interpreter='latex', FontSize=16)
xlim([0 16])

exportgraphics(f1,"vanderCL.eps", "Resolution",600)