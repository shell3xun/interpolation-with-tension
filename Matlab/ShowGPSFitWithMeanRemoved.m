% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

shouldSaveFigures = 0;

% Drifter to highlight in the final plots. Drifter 7 has no outliers
choiceDrifter = 6;

drifters = open('sample_data/raw_rho1_drifters.mat');

t_drifter = (drifters.date{choiceDrifter}-drifters.lastDeployment)*24*60*60;
lat = drifters.lat{choiceDrifter};
lon = drifters.lon{choiceDrifter};
lon0 = min(lon)+(max(lon)-min(lon))/2;
spline = GPSTensionSpline(t_drifter,lat,lon,'shouldUseRobustFit',1);

spline.estimateOutlierDistribution();
mse1 = spline.minimizeExpectedMeanSquareErrorInNoiseRange();

tq = linspace(min(spline.t),max(spline.t),10*length(spline.t));
[xq,yq] = spline.xyAtTime(tq);

figure
scatter(spline.x,spline.y,(2.5)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
plot(xq,yq,'k','LineWidth',1.5),axis equal

spline.setSigmaFromFullTensionSolution();
mse2 = spline.minimizeExpectedMeanSquareErrorInNoiseRange();

[xq,yq] = spline.xyAtTime(tq);
plot(xq,yq)

t_drifter = (drifters.date{7}-drifters.lastDeployment)*24*60*60;
spline7 = GPSTensionSpline(t_drifter,drifters.lat{7},drifters.lon{7},'lon0',spline.lon0);

tq = linspace(min(spline7.t),max(spline7.t),10*length(spline7.t));
[xq,yq] = spline7.xyAtTime(tq);
plot(xq,yq)

% figure
% scatter(t,y,(2.5)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
% plot(tq,yq)