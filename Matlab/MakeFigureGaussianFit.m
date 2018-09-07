% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('./support');

drifters = load('sample_data/rho1_drifters_projected_ungridded.mat');

% Drifter to highlight in the final plots
choiceDrifter = 6;

S = 3; % order of the spline
K = S+1;
T = 2; % order of the tension

% Pull out the data of interest
x = drifters.x{choiceDrifter};
y = drifters.y{choiceDrifter};
t = drifters.t{choiceDrifter};
N = length(t);

% In my original version of this figure, I scaled the lambda for
% spline_fit_big_error by a factor of 1e4 manually. This gives quite a
% different answer than this approach, because here we (probably) aren't
% estimating u_rms correctly, because the error is incorrect.
spline_fit_small_error = TensionSpline(t,[x y],10);
spline_fit_big_error = TensionSpline(t,[x y],50);

tq = linspace(t(1),t(end),10*length(t))';

[x_fit_small, y_fit_small] = spline_fit_big_error(tq);
[x_fit_big, y_fit_big] = spline_fit_small_error(tq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position fit figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(2,2,[1 3])
s = 1/1000;
plot(s*x_fit_small,s*y_fit_small, 'LineWidth', 0.5*scaleFactor, 'Color', 'b'), hold on
plot(s*x_fit_big,s*y_fit_big, 'LineWidth', 1.0*scaleFactor, 'Color','k')
scatter(s*drifters.x{choiceDrifter},s*drifters.y{choiceDrifter},(5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('y (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
% xlim([-10 -7.4])
% ylim([13.5 19])

sp1 = subplot(2,2,2);
plot(tq/3600,s*x_fit_small, 'LineWidth', 0.5*scaleFactor, 'Color','b'), hold on
plot(tq/3600,s*x_fit_big, 'LineWidth', 1.0*scaleFactor, 'Color','k')
scatter(drifters.t{choiceDrifter}/3600,s*drifters.x{choiceDrifter},(5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
ylabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
set(sp1,'yaxislocation','right')
% xlim([min(t) 40])
% ylim([-10 -7.4])

sp1 = subplot(2,2,4);
plot(tq/3600,s*y_fit_small, 'LineWidth', 0.5*scaleFactor, 'Color','b'), hold on
plot(tq/3600,s*y_fit_big, 'LineWidth', 1.0*scaleFactor, 'Color','k')
scatter(drifters.t{choiceDrifter}/3600,s*drifters.y{choiceDrifter},(5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('t (hours)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('y (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
set(sp1,'yaxislocation','right')
% xlim([min(t) 40])
% ylim([13.5 19])


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

plot(tq/3600,s*x_fit_small, 'LineWidth', 0.5*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]), hold on
plot(tq/3600,s*x_fit_big, 'LineWidth', 1.0*scaleFactor, 'Color',0.0*[1.0 1.0 1.0])
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

print('-depsc2', '../figures/gaussianfit.eps')



