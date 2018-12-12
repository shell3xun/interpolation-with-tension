% Site 1 or site 2?
Site = 1;


% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

shouldSaveFigures = 0;

% Drifter to highlight in the final plots
choiceDrifter = 6;

if Site == 1
    drifters = load('sample_data/rho1_drifters_projected_ungridded.mat');
else
    drifters = load('sample_data/rho2_drifters_projected_ungridded.mat');
end

% Pull out the data of interest
x_data = drifters.x{choiceDrifter};
y_data = drifters.y{choiceDrifter};
t_data = drifters.t{choiceDrifter};



splinefit = TensionSpline(t_data,x_data,10);

gpsfit = GPSTensionSpline(t_data,x_data,y_data, 'shouldIdentifyOutliers', 0);

t_tension = gpsfit.spline_x.t_knot(gpsfit.spline_x.K:1:end-gpsfit.spline_x.K+1);
t_tension = t_tension(1:end-1) + diff(t_tension)/2;
[x_T,I] = sort(gpsfit.spline_x(t_tension,3));
x_T_Q1 = median(x_T(1:floor(length(x_T)/2)));
x_T_Q3 = median(x_T(ceil(length(x_T)/2)+1:end));
sigma_T = (x_T_Q3-x_T_Q1)/1.349;
p_x_T = 1-erf(abs(x_T)/sigma_T/sqrt(2));
T_x_T = t_tension(I);

indicesOfOutliers = p_x_T < 0.001;

% for odd K the knots are between the points, even they're on the points

figure, histogram(x_T,100), vlines([x_T_Q1,x_T_Q3],'k')

t=linspace(min(t_data),max(t_data),length(t_data)*10).';
[x,y] = gpsfit(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position fit figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t2 = sort(T_x_T(indicesOfOutliers));

figure
s = 1/1000; % scale
plot(s*x,s*y, 'LineWidth', 0.5*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]), hold on
% scatter(s*x_data(gpsfit.indicesOfOutliers),s*y_data(gpsfit.indicesOfOutliers),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(s*gpsfit.spline_x(t2),s*gpsfit.spline_y(t2),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(s*x_data,s*y_data,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position fit figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigureSize = [50 50 figure_width_2col+8 150*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

plot(t/3600,s*x, 'LineWidth', 1.0*scaleFactor, 'Color',0.0*[1.0 1.0 1.0]), hold on
plot(t/3600,s*splinefit(t), 'LineWidth', 1.0*scaleFactor, 'Color',0.4*[1.0 1.0 1.0])
scatter(drifters.t{choiceDrifter}(gpsfit.indicesOfOutliers)/3600,s*drifters.x{choiceDrifter}(gpsfit.indicesOfOutliers),(6.5*scaleFactor)^2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(drifters.t{choiceDrifter}/3600,s*drifters.x{choiceDrifter},(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('t (hours)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)

xlim([124 149])
ylim([5.58 9.18])

fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';

if shouldSaveFigures == 1
    print('-depsc2', '../figures/gpsfit.eps')
end

figure
plot(t/3600,s*y, 'LineWidth', 1.0*scaleFactor, 'Color',0.0*[1.0 1.0 1.0]), hold on
scatter(drifters.t{choiceDrifter}(gpsfit.indicesOfOutliers)/3600,s*drifters.y{choiceDrifter}(gpsfit.indicesOfOutliers),(6.5*scaleFactor)^2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(drifters.t{choiceDrifter}/3600,s*drifters.y{choiceDrifter},(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('t (hours)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)

xlim([124 149])