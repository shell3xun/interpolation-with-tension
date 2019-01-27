% Site 1 or site 2?
Site = 1;


% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

shouldSaveFigures = 0;

% Drifter to highlight in the final plots. Drifter 7 has no outliers
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

tq = linspace(min(t_data),max(t_data),10*length(t_data));

% These are our working definitions for the noise
noiseDistribution = StudentTDistribution(8.5,4.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make a figure with the data
figure
sp1=subplot(2,1,1);
scatter(t_data,x_data,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
% errorbar(t_data,x_data,zmax*ones(size(t_data))), hold on
sp2=subplot(2,1,2);
scatter(t_data,y_data,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
% errorbar(t_data,y_data,zmax*ones(size(t_data))), hold on

spline_x = RobustTensionSpline(t_data,x_data,noiseDistribution, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
spline_y = RobustTensionSpline(t_data,y_data,noiseDistribution, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);

beta = 1/50;
spline_x.firstIteration(beta);
spline_y.firstIteration(beta);

% spline_x.removeOutlierKnotsAndRetension(1/100);
% spline_y.removeOutlierKnotsAndRetension(1/100);

% spline_x.secondIteration();
% spline_y.secondIteration();

% spline_x.secondIteration();
% spline_y.secondIteration();

fprintf('robust lambda: (%.3g, %.3g)\n', spline_x.lambda,spline_y.lambda );
fprintf('robust e-mse: (%.3g, %.3g)\n', spline_x.expectedMeanSquareErrorInPercentileRange(beta/2,1-beta/2),spline_y.expectedMeanSquareErrorInPercentileRange(beta/2,1-beta/2));
fprintf('robust a_rms: (%.3g, %.3g)\n', std(spline_x.uniqueValuesAtHighestDerivative), std(spline_y.uniqueValuesAtHighestDerivative) );
fprintf('robust a_rms: (%.3g, %.3g)\n', TensionSpline.StandardDeviationAndMeanOfInterquartileRange(spline_x.uniqueValuesAtHighestDerivative), TensionSpline.StandardDeviationAndMeanOfInterquartileRange(spline_y.uniqueValuesAtHighestDerivative) );

%%%%%%%%%%%%%%%%%%%%%
% Grab drifter 7 and plot that
spline_7x = TensionSpline(drifters.t{7},drifters.x{7},noiseDistribution);

subplot(sp1)
scatter(t_data(spline_x.indicesOfOutliers),x_data(spline_x.indicesOfOutliers),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(t_data(spline_x.indicesOfOutliers),x_data(spline_x.indicesOfOutliers),(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline_x(tq))
plot(tq,spline_7x(tq))
subplot(sp2)
scatter(t_data(spline_y.indicesOfOutliers),y_data(spline_y.indicesOfOutliers),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(t_data(spline_y.indicesOfOutliers),y_data(spline_y.indicesOfOutliers),(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline_y(tq))