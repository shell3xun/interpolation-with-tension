% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

shouldSaveFigures = 1;

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
spline = GPSSmoothingSpline(t_drifter,lat,lon,'shouldUseRobustFit',1,'lon0',lon0);
% spline.lambda

spline_nonrobust = GPSSmoothingSpline(t_drifter,lat,lon,'shouldUseRobustFit',0,'lon0',lon0);
% spline_nonrobust.lambda
spline_nonrobust.minimize( @(spline) spline.expectedMeanSquareErrorFromCV );
% spline_nonrobust.lambda

spline7 = GPSSmoothingSpline((drifters.date{7}-drifters.lastDeployment)*24*60*60,drifters.lat{7},drifters.lon{7},'shouldUseRobustFit',1,'lon0',lon0);

% Pull out the data of interest
x_data = spline.x;
y_data = spline.y;
t_data = spline.t;

t=linspace(min(t_data),max(t_data),length(t_data)*10).';
[xq,yq] = spline.xyAtTime(t);
[xq_nr,yq_nr] = spline_nonrobust.xyAtTime(t);
[xq7,yq7] = spline7.xyAtTime(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position fit figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
datapointColor = 0.5*[1.0 1.0 1.0];
datapointSize = 2.0*scaleFactor;
outlierFaceColor = 1.0*[1.0 1.0 1.0];
outlierEdgeColor = 0.2*[1.0 1.0 1.0];
outlierSize = 5.5;
nonOutlierIndices = 1:length(spline.t);
nrColor = [0.8500    0.3250    0.0980];

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

scatter(s*(spline.x(spline.outlierIndices)-x0),s*(spline.y(spline.outlierIndices)-y0),(outlierSize)^2, 'MarkerEdgeColor', outlierEdgeColor, 'MarkerFaceColor', outlierFaceColor); hold on
scatter(s*(spline.x(nonOutlierIndices)-x0),s*(spline.y(nonOutlierIndices)-y0),(datapointSize)^2,'filled', 'MarkerEdgeColor', datapointColor, 'MarkerFaceColor', datapointColor)
ax1 = plot(s*(xq_nr-x0),s*(yq_nr-y0), 'LineWidth', 1.0*scaleFactor, 'Color',nrColor);
plot(s*(xq7-x0),s*(yq7-y0), 'LineWidth', 1.0*scaleFactor, 'Color',0.4*[1.0 1.0 1.0])
plot(s*(xq-x0),s*(yq-y0), 'LineWidth', 1.0*scaleFactor, 'Color',0.0*[1.0 1.0 1.0])

xlabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('y (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
axis equal
xlim([-1 13])
ylim([22 36])

xlim([-0.5 3.5])
ylim([-0.5 4.5])

fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';
set(gca,'Box','on')

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
scatter(spline.t(spline.outlierIndices)/3600,s*(spline.x(spline.outlierIndices)-x0),(outlierSize)^2, 'MarkerEdgeColor', outlierEdgeColor, 'MarkerFaceColor', outlierFaceColor);  hold on
scatter(spline.t(nonOutlierIndices)/3600,s*(spline.x(nonOutlierIndices)-x0),(datapointSize)^2,'filled', 'MarkerEdgeColor', datapointColor, 'MarkerFaceColor', datapointColor)
plot(t/3600,s*(xq_nr-x0), 'LineWidth', 1.0*scaleFactor, 'Color',nrColor), hold on
plot(t/3600,s*(xq7-x0), 'LineWidth', 1.0*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]), hold on
plot(t/3600,s*(xq-x0), 'LineWidth', 1.0*scaleFactor, 'Color',0.0*[1.0 1.0 1.0])

ylabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp1.YAxisLocation = 'right';
sp1.XTickLabel = [];
xlim([0 45])

sp2 = subplot(2,1,2);
p3 = scatter(spline.t(spline.outlierIndices)/3600,s*(spline.y(spline.outlierIndices)-y0),(outlierSize)^2, 'MarkerEdgeColor', outlierEdgeColor, 'MarkerFaceColor', outlierFaceColor); hold on
p4 = scatter(spline.t(nonOutlierIndices)/3600,s*(spline.y(nonOutlierIndices)-y0),(datapointSize)^2,'filled', 'MarkerEdgeColor', datapointColor, 'MarkerFaceColor', datapointColor);
p0 = plot(t/3600,s*(yq_nr-y0), 'LineWidth', 1.0*scaleFactor, 'Color',nrColor); hold on
p1 = plot(t/3600,s*(yq7-y0), 'LineWidth', 1.0*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]); hold on
p2 = plot(t/3600,s*(yq-y0), 'LineWidth', 1.0*scaleFactor, 'Color',0.0*[1.0 1.0 1.0]);

xlabel('t (hours)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('y (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp2.YAxisLocation = 'right';
xlim([0 45])

legend([p1 p0 p2 p4 p3],{'drifter 7 fit', 'drifter 6 cv fit', 'drifter 6 ranged fit', 'drifter 6 data','outlier'},'Location','southeast')

packfig(2,1)

fig2 = tightfig;
fig2.Position = FigureSize;
fig2.PaperPosition = FigureSize;
fig2.PaperSize = [FigureSize(3) FigureSize(4)];
fig2.PaperPositionMode = 'auto';

sp1.Box = 'on';
sp2.Box = 'on';

if shouldSaveFigures == 1
    print('-depsc2', '../figures/gpsfit_xtyt.eps')
end

return


ylim([5.58 9.18])



return



figure


% xlim([124 149])

figure
