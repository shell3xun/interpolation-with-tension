% Site 1 or site 2?
Site = 1;


% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

shouldSaveFigures = 0;

% Drifter to highlight in the final plots. Drifter 7 has no outliers
choiceDrifter = 7;

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

spline_x = RobustTensionSpline(t_data,x_data,noiseDistribution);
spline_y = RobustTensionSpline(t_data,y_data,noiseDistribution);

%%%%%%%%%%%%%%%%%%%%%
% Grab drifter 7 and plot that
spline_7x = TensionSpline(drifters.t{7},drifters.x{7},noiseDistribution);

subplot(sp1)
scatter(t_data(spline_x.outlierIndices),x_data(spline_x.outlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(t_data(spline_x.outlierIndices),x_data(spline_x.outlierIndices),(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline_x(tq))
plot(tq,spline_7x(tq))


fprintf('robust a_rms: (%.3g, %.3g)\n', std(spline_x.uniqueValuesAtHighestDerivative), std(spline_y.uniqueValuesAtHighestDerivative) );
fprintf('robust a_rms: (%.3g, %.3g)\n', TensionSpline.StandardDeviationAndMeanOfInterquartileRange(spline_x.uniqueValuesAtHighestDerivative), TensionSpline.StandardDeviationAndMeanOfInterquartileRange(spline_y.uniqueValuesAtHighestDerivative) );


subplot(sp2)
scatter(t_data(spline_y.outlierIndices),y_data(spline_y.outlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(t_data(spline_y.outlierIndices),y_data(spline_y.outlierIndices),(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline_y(tq))


return

[fn,sn] = powspec(tq(2)-tq(1),diff(spline_x2(tq)')/(tq(2)-tq(1)));
figure
plot(fn,sn), xlog, ylog
