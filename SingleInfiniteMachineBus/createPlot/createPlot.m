load("sol_SMIB_CL_16_03_25.mat")
load("sol_SMIB_OL_17_03_25.mat")
load("sol_SMIB_CL_multiConds2.mat")

xx=linspace(solutionOL.T(1,1),solutionOL.tf,10000);
%%
data.H=0.0106;
data.Xt=0.28;
data.Pm=1;
data.Es=1.21;
data.V=1;
data.Pe=data.Es*data.V/data.Pm/data.Xt;
data.D=0.03;
data.dep=asin(1/data.Pe);
data.C1 = 3;
data.C2 = 30;
%% two plot
f = figure(1);
f.Position = [100, 100, 900, 500];
a1 = axes('FontSize', 14, 'NextPlot', 'add');
a1.GridLineStyle = '--';
a1.GridColor = 'k';
a1.GridAlpha = 0.9 ;
hold on

p3 = plot(solutionCL.t_sol, solutionCL.x_sol(:,1)/3, 'mo');
p4 = plot(solutionCL.t_sol, solutionCL.x_sol(:,2)/30, 'mo');
p1 = plot(xx,speval(solutionOL,'X',1,xx)/3 , '-k',LineWidth=2);
p2 = plot(xx,speval(solutionOL,'X',2,xx)/30 , '-k',LineWidth=2);
xlabel('Time [s]',"Interpreter","latex")
ylabel('States',"Interpreter","latex")
xlim([0 4])
legend([p1, p3], {"Open Loop", "Closed Loop"}, "Interpreter","latex", "FontSize", 14)
grid on
exportgraphics(f,"SMIBStateSubfig.eps", "Resolution",600)

%analytic curve
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
x_1 = solutionCL.x_sol(:,1);
x_2 = solutionCL.x_sol(:,2);
tsw = find(solutionCL.t_sol >= 0.32, 1);
u_analytic = -(H*c2^2.*((2.*x_1)./c1^2 - (Pm - D.*x_2)./(H*c2^2)))./(Pe*sin(dep + x_1));
GLC = -Pe.*sin(x_1(tsw:end) + dep)./(c2^2 * H);


f = figure(2);
f.Position = [100, 100, 900, 500];
a1 = axes('FontSize', 14, 'NextPlot', 'add');
a1.GridLineStyle = '--';
a1.GridColor = 'k';
a1.GridAlpha = 0.9 ;
hold on
b1 = plot([solutionOL.T(1,1); solutionOL.tf],[-1, -1],'r-' );
b2 = plot([solutionOL.T(1,1); solutionOL.tf],[1, 1],'r-' );
p2 = plot(solutionCL.t_sol(1:end-1,:), solutionCL.u_sol(1:end, :), 'mo',LineWidth=1 );
p1 = plot(xx,speval(solutionOL,'U',1,xx), 'k-',LineWidth=2);
p3 = plot(solutionCL.t_sol(tsw:end), u_analytic(tsw:end), 'g-', LineWidth=3);

xlim([0 4])
ylim([-1.2 1.2])
xlabel('Time [s]',"Interpreter","latex")
legend([p1,p2,p3, b1], {"Open Loop", "Closed Loop", "Analytic curve", "Control bounds"}, "Interpreter","latex", "FontSize", 14)
grid on
ylabel('Control',"Interpreter","latex")
exportgraphics(f,"SMIBControlSubfig.eps", "Resolution",600)
%% tiled single plot
f = figure(3);
f.Position = [100, 100, 900, 900];
t = tiledlayout(2,1);
t.TileSpacing = 'tight';
a1 = nexttile;
a1.XAxis.FontSize = 16;
a1.YAxis.FontSize = 16;
a1.GridLineStyle = '--';
a1.GridColor = 'k';
a1.GridAlpha = 0.9 ;
hold on
p3 = plot(solutionCL.t_sol, solutionCL.x_sol(:,1)/3, 'mo');
p4 = plot(solutionCL.t_sol, solutionCL.x_sol(:,2)/30, 'mo');
p1 = plot(xx,speval(solutionOL,'X',1,xx)/3 , '-k',LineWidth=2);
p2 = plot(xx,speval(solutionOL,'X',2,xx)/30 , '-k',LineWidth=2);
% xlabel('Time [s]',"Interpreter","latex")
ylabel('States',"Interpreter","latex")
xlim([0 4])
% legend([p1, p3], {"Open Loop", "Closed Loop"}, "Interpreter","latex", "FontSize", 14)
grid on

a2 = nexttile;
a2.XAxisLocation = "top";
a2.XAxis.FontSize = 16;
a2.YAxis.FontSize = 16;
a2.GridLineStyle = '--';
a2.GridColor = 'k';
a2.GridAlpha = 0.9 ;
hold on
b1 = plot([solutionOL.T(1,1); solutionOL.tf],[-1, -1],'r-' );
b2 = plot([solutionOL.T(1,1); solutionOL.tf],[1, 1],'r-' );
p2 = plot(solutionCL.t_sol(1:end-1,:), solutionCL.u_sol(1:end, :), 'mo',LineWidth=1 );
p1 = plot(xx,speval(solutionOL,'U',1,xx), 'k-',LineWidth=2);
p3 = plot(solutionCL.t_sol(tsw:end), u_analytic(tsw:end), 'g-', LineWidth=3);
xlim([0 4])
ylim([-1.2 1.2])
xlabel('Time [s]',"Interpreter","latex")

lgd = legend([p1,p2,p3, b1], {"Open Loop", "Closed Loop", "Analytic curve", "Control bounds"}, "Interpreter","latex");
grid on
ylabel('Control',"Interpreter","latex")

lgd.Layout.Tile = 'south';
lgd.FontSize = 18;
lgd.NumColumns = 3;

exportgraphics(f,"SMIB_Tiled.eps", "Resolution",600)

%% 


figure(5)
plot(solutionCL.t_sol(tsw:end), GLC)
grid on 
xlabel('Time [s]',"Interpreter","latex")
ylabel('GLC',"Interpreter","latex")
hold off


%phase potrait
x1_range = linspace(-2.5,1.25,20);
x2_range = linspace(-15,15, 20);
[x,y] = meshgrid(x1_range, x2_range);
u = zeros(size(x));
v = zeros(size(x));
t=0;
for i = 1:numel(x)
    derivs = singular_interval_model(t, [x(i), y(i)], data);
    u(i) = derivs(1);
    v(i) = derivs(2);
end

figure(7)
quiver(x,y,u,v,'r')
xlabel('x_1')
ylabel('x_2')

figure(6)

hold on
% quiver(x,y,u,v,'r')
xlabel('x_1')
ylabel('x_2')
% plot(solutionCL.x_sol(:,1), solutionCL.x_sol(:,2), "m-", LineWidth=2);
% plot(solutionCL.x_sol(1,1), solutionCL.x_sol(1,2), "m*", LineWidth=2);
for j = 1:5
    xtemp = cell2mat(sol.x_sol(j));
    plot(xtemp(:,1), xtemp(:,2), LineWidth=2);
    plot(xtemp(1,1), xtemp(1,2), "m*", LineWidth=2);
end
% psw = plot(solutionCL.x_sol(tsw,1), solutionCL.x_sol(tsw,2), "k^", MarkerSize=12, LineWidth=2);
hold off
xlim([-3.5, 5])
ylim([-20,20])
% legend(p1, {"Closed loop solution"})
grid on
%%

function [dxdt, u] = singular_interval_model(t, x, data)
dxdt = zeros(size(x));
Pm = data.Pm;
Pe = data.Pe;
C1 = data.C1;
C2 = data.C2;
D = data.D;
dep = data.dep;
H = data.H;
x1 = x(1);
x2 = x(2);
%u
% u = -(C1^2.*D.*x2 - C1^2*Pm + 2.*C1^2.*H.*x2 + 2.*C2^2.*H.*x1)/(C1^2.*Pe.*sin(dep + x1));
u = -(C1^2.*D.*x2 - C1^2*Pm + 2.*C1^2.*H.*x2 + 2.*C2^2.*H.*x1)/(C1^2);

% dxdt(1) = x2;
% dxdt(2) = ((Pm-D*x2) - Pe *sin(x1+dep).*u)/(2*H);
dxdt(1) = x2;
dxdt(2) = ((Pm-D*x2) - u)/(2*H);


end