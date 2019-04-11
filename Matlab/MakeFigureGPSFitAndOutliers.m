% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

shouldSaveFigures = 0;

% Drifter to highlight in the final plots
choiceDrifter = 6;

drifters = open('sample_data/raw_rho1_drifters.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We need a consistent lon0 for all drifters
lon = drifters.lon{1};
lon0 = min(lon)+(max(lon)-min(lon))/2;

t_drifter = (drifters.date{choiceDrifter}-drifters.lastDeployment)*24*60*60;
lat = drifters.lat{choiceDrifter};
lon = drifters.lon{choiceDrifter};
spline = GPSTensionSpline(t_drifter,lat,lon,'shouldUseRobustFit',1,'lon0',lon0);

spline7 = GPSTensionSpline((drifters.date{7}-drifters.lastDeployment)*24*60*60,drifters.lat{7},drifters.lon{7},'shouldUseRobustFit',1,'lon0',lon0);

% Pull out the data of interest
x_data = spline.x;
y_data = spline.y;
t_data = spline.t;

t=linspace(min(t_data),max(t_data),length(t_data)*10).';
[xq,yq] = spline.xyAtTime(t);
[xq7,yq7] = spline7.xyAtTime(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position fit figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
FigureSize = [50 50 figure_width_1col+8 250*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

s = 1/1000; % scale for x and y axis
x0 = min(xq7);
y0 = yq7(1);

plot(s*(xq7-x0),s*(yq7-y0), 'LineWidth', 1.0*scaleFactor, 'Color',0.5*[1.0 1.0 1.0]), hold on
plot(s*(xq-x0),s*(yq-y0), 'LineWidth', 1.0*scaleFactor, 'Color',0.0*[1.0 1.0 1.0])
scatter(s*(spline.x(spline.outlierIndices)-x0),s*(spline.y(spline.outlierIndices)-y0),(6.5*scaleFactor)^2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(s*(spline.x-x0),s*(spline.y-y0),(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('y (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
axis equal
xlim([-1 13])
ylim([22 36])

xlim([-0.5 3.0])
ylim([-0.5 4.0])

fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';

if shouldSaveFigures == 1
    print('-depsc2', '../figures/gpsfit_xy.eps')
end

FigureSize = [50 50 figure_width_2col+8 250*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

sp1 = subplot(2,1,1);
plot(t/3600,s*(xq7-x0), 'LineWidth', 1.0*scaleFactor, 'Color',0.5*[1.0 1.0 1.0]), hold on
plot(t/3600,s*(xq-x0), 'LineWidth', 1.0*scaleFactor, 'Color',0.0*[1.0 1.0 1.0])
scatter(spline.t(spline.outlierIndices)/3600,s*(spline.x(spline.outlierIndices)-x0),(6.5*scaleFactor)^2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(spline.t/3600,s*(spline.x-x0),(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
ylabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp1.YAxisLocation = 'right';
sp1.XTickLabel = [];
xlim([0 40])

sp2 = subplot(2,1,2);
p1 = plot(t/3600,s*(yq7-y0), 'LineWidth', 1.0*scaleFactor, 'Color',0.5*[1.0 1.0 1.0]); hold on
p2 = plot(t/3600,s*(yq-y0), 'LineWidth', 1.0*scaleFactor, 'Color',0.0*[1.0 1.0 1.0]);
p3 = scatter(spline.t(spline.outlierIndices)/3600,s*(spline.y(spline.outlierIndices)-y0),(6.5*scaleFactor)^2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
p4 = scatter(spline.t/3600,s*(spline.y-y0),(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
xlabel('t (hours)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('y (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp2.YAxisLocation = 'right';
xlim([0 40])

legend([p1 p2 p4 p3],{'drifter 7 fit', 'drifter 6 fit', 'drifter 6 data','outlier'},'Location','southeast')

packfig(2,1)

fig2 = tightfig;
fig2.Position = FigureSize;
fig2.PaperPosition = FigureSize;
fig2.PaperSize = [FigureSize(3) FigureSize(4)];
fig2.PaperPositionMode = 'auto';

if shouldSaveFigures == 1
    print('-depsc2', '../figures/gpsfit_xtyt.eps')
end

return


ylim([5.58 9.18])



return



figure


% xlim([124 149])

figure
