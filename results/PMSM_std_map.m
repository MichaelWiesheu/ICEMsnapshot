

load('InitValues.mat')
TmeanMapInit = mean(TorquesMapInit, 3);
TstdMapInit = std(TorquesMapInit, 1, 3);

%%
[xx,yy] = meshgrid(currentsMapInit, PhaseAnglesInit);
figure
s = surf(xx,yy,TstdMapInit');
% shading interp
xlim([min(currentsMapInit), max(currentsMapInit)])
ylim([min(PhaseAnglesInit), max(PhaseAnglesInit)])
xlabel("Phase Current")
ylabel("Phase Angle")
title("Initial Motor Torque std")
colorbar
view(2)

%%
figure
surf(xx,yy,TmeanMapInit')
shading interp
xlim([min(currentsMapInit), max(currentsMapInit)])
ylim([min(PhaseAnglesInit), max(PhaseAnglesInit)])
title("Initial Motor Torque mean")
colorbar
view(2)

%%

load('OptValues.mat')
TmeanMapOpt = mean(TorquesMapOpt, 3);
TstdMapOpt = std(TorquesMapOpt, 1, 3);

%%
[xx,yy] = meshgrid(currentsMapOpt, PhaseAnglesOpt);
figure
s = surf(xx,yy,TstdMapOpt');
% shading interp
xlim([min(currentsMapOpt), max(currentsMapOpt)])
ylim([min(PhaseAnglesOpt), max(PhaseAnglesOpt)])
xlabel("Phase Current")
ylabel("Phase Angle")
title("Optimized Motor Torque std")
colorbar
view(2)

%%
figure
surf(xx,yy,TmeanMapOpt')
shading interp
xlim([min(currentsMapOpt), max(currentsMapOpt)])
ylim([min(PhaseAnglesOpt), max(PhaseAnglesOpt)])
title("Optimized Motor Torque mean")
colorbar
view(2)

%%
set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex');
f = figure("Position",[0,0,500,300]);
currentsRel = currentsMapOpt/max(currentsMapOpt);
[xx,yy] = meshgrid(currentsRel, PhaseAnglesOpt);

Improvement = (TstdMapInit-TstdMapOpt)./TstdMapInit*100;
surf(xx,yy,Improvement')

% colormap winter
npts = 20;

c1 = [185, 15, 34];
c2 = [153, 192, 0];
c2 = [127, 171, 22];
% colors = [linspace(c1(1),c2(1),npts)', linspace(c1(2),c2(2),npts)', linspace(c1(3),c2(3),npts)']/255;
colors1 = [linspace(c1(1),255, npts)', linspace(c1(2),255,npts)', linspace(c1(3),255,npts)']/255;
colors2 = [linspace(255,c2(1), npts)', linspace(255,c2(2),npts)', linspace(255,c2(3),npts)']/255;
colors = [colors1; colors2];

colormap(colors)
view(2)
xticks(linspace(0,1,11))

shading interp
c = colorbar;
c.TickLabelInterpreter = 'latex';
c.Label.Interpreter = 'latex';
c.Label.String = 'Torque std improvement (\%)';
c.FontSize = 16;
c.Label.FontSize = 18;
clim([-100,100])

ylim([-90,90])
yticks([-90:15:90])

set(gca, 'FontSize', 16)
xlabel('$J/J_{\mathrm{0}}$', 'Interpreter','latex')
ylabel("Electric phase angle ($^\circ$)", 'Interpreter','latex')

exportgraphics(f, 'TorqueRippleReduction.png')
exportgraphics(f, 'TorqueRippleReduction.pdf')

