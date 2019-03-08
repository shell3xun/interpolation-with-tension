% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

shouldSaveFigures = 0;

% Drifter to highlight in the final plots. Drifter 7 has no outliers
choiceDrifter = 6;

drifters = open('sample_data/raw_rho1_drifters.mat');

t_drifter = (drifters.date{choiceDrifter}-drifters.lastDeployment)*24*60*60;
spline = GPSTensionSpline(t_drifter,drifters.lat{choiceDrifter},drifters.lon{choiceDrifter},'shouldUseRobustFit',1);

t = t_drifter;
x = spline.x;
y = spline.y;

spline_x_constrained = TensionSpline(t,x,spline.distribution,'S',spline.spline_x.S+1);
spline_x_constrained.lambda = 1e35;
spline_y_constrained = TensionSpline(t,y,spline.distribution,'S',spline.spline_y.S+1);
spline_y_constrained.lambda = 1e32;

x_rms = x-spline_x_constrained(t);
y_rms = y-spline_y_constrained(t);

spline_x_rms = RobustTensionSpline(t,x_rms,spline.distribution,'S',spline.spline_x.S);
spline_y_rms = RobustTensionSpline(t,y_rms,spline.distribution,'S',spline.spline_y.S);

tq = linspace(min(spline.t),max(spline.t),10*length(spline.t));

figure
scatter(x,y,(2.5)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
plot(spline_x_constrained(tq) + spline_x_rms(tq),spline_y_constrained(tq) + spline_y_rms(tq)),axis equal