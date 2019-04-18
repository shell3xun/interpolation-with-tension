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
outlierDistribution = StudentTDistribution(200,3.0);
distribution = AddedDistribution(0.01,outlierDistribution,noiseDistribution);
zmin = noiseDistribution.locationOfCDFPercentile(1/100/2);
zmax = noiseDistribution.locationOfCDFPercentile(1-1/100/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make a figure with the data
figure
sp1=subplot(2,1,1);
% scatter(t_data,x_data,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
errorbar(t_data,x_data,zmax*ones(size(t_data))), hold on
sp2=subplot(2,1,2);
errorbar(t_data,y_data,zmax*ones(size(t_data))), hold on

% scatter(t_data,y_data,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Start with the usual optimal iterated---no outliers

spline_x = SmoothingSpline(t_data,x_data,noiseDistribution,'lambda',Lambda.optimalIterated);
spline_y = SmoothingSpline(t_data,y_data,noiseDistribution,'lambda',Lambda.optimalIterated);

fprintf('optimal iterated no-outliers lambda: (%.3g, %.3g)\n', spline_x.lambda,spline_y.lambda );
fprintf('optimal iterated no-outliers a_rms: (%.3g, %.3g)\n', std(spline_x.uniqueValuesAtHighestDerivative), std(spline_y.uniqueValuesAtHighestDerivative) );
fprintf('optimal iterated no-outliers a_rms: (%.3g, %.3g)\n', SmoothingSpline.StandardDeviationAndMeanOfInterquartileRange(spline_x.uniqueValuesAtHighestDerivative), SmoothingSpline.StandardDeviationAndMeanOfInterquartileRange(spline_y.uniqueValuesAtHighestDerivative) );

subplot(sp1)
plot(tq,spline_x(tq))
subplot(sp2)
plot(tq,spline_y(tq))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Start with the usual optimal iterated---with outliers (this does badly)

tensionDistribution = NormalDistribution(SmoothingSpline.StandardDeviationOfInterquartileRange(spline_x.uniqueValuesAtHighestDerivative));
logLikelihood = @(spline) -sum(noiseDistribution.logPDF( spline.epsilon ) ) - sum(tensionDistribution.logPDF(spline.uniqueValuesAtHighestDerivative));
spline_x.minimize( logLikelihood );

tensionDistribution = NormalDistribution(SmoothingSpline.StandardDeviationOfInterquartileRange(spline_y.uniqueValuesAtHighestDerivative));
logLikelihood = @(spline) -sum(noiseDistribution.logPDF( spline.epsilon ) ) - sum(tensionDistribution.logPDF(spline.uniqueValuesAtHighestDerivative));
spline_y.minimize( logLikelihood );

fprintf('log-likelihood no-outliers lambda: (%.3g, %.3g)\n', spline_x.lambda,spline_y.lambda );
fprintf('log-likelihood no-outliers a_rms: (%.3g, %.3g)\n', std(spline_x.uniqueValuesAtHighestDerivative), std(spline_y.uniqueValuesAtHighestDerivative) );
fprintf('log-likelihood no-outliers a_rms: (%.3g, %.3g)\n', SmoothingSpline.StandardDeviationAndMeanOfInterquartileRange(spline_x.uniqueValuesAtHighestDerivative), SmoothingSpline.StandardDeviationAndMeanOfInterquartileRange(spline_y.uniqueValuesAtHighestDerivative) );

subplot(sp1)
plot(tq,spline_x(tq))
subplot(sp2)
plot(tq,spline_y(tq))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Then iterate used the expected MSE within some range (this does well)

spline_x = SmoothingSpline(t_data,x_data,distribution,'lambda',Lambda.optimalIterated);
spline_y = SmoothingSpline(t_data,y_data,distribution,'lambda',Lambda.optimalIterated);

spline_x.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
spline_y.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );

fprintf('ranged iterated lambda: (%.3g, %.3g)\n', spline_x.lambda,spline_y.lambda );
fprintf('ranged iterated a_rms: (%.3g, %.3g)\n', std(spline_x.uniqueValuesAtHighestDerivative), std(spline_y.uniqueValuesAtHighestDerivative) );
fprintf('ranged iterated a_rms: (%.3g, %.3g)\n', SmoothingSpline.StandardDeviationAndMeanOfInterquartileRange(spline_x.uniqueValuesAtHighestDerivative), SmoothingSpline.StandardDeviationAndMeanOfInterquartileRange(spline_y.uniqueValuesAtHighestDerivative) );

subplot(sp1)
plot(tq,spline_x(tq))
subplot(sp2)
plot(tq,spline_y(tq))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Identify the points that look like outliers
outlierThreshold = noiseDistribution.locationOfCDFPercentile(1-1/10000/2);

spline_x.outlierIndices = find(abs(spline_x.epsilon) > outlierThreshold);
spline_y.outlierIndices = find(abs(spline_y.epsilon) > outlierThreshold);

t_knot_x = InterpolatingSpline.KnotPointsForPoints(spline_x.t(spline_x.nonOutlierIndices),spline_x.K,1);
t_knot_y = InterpolatingSpline.KnotPointsForPoints(spline_y.t(spline_y.nonOutlierIndices),spline_y.K,1);

fprintf('identified outliers: (%d, %d)\n',length(spline_x.outlierIndices), length(spline_y.outlierIndices));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Remove spline support from outlier knots, then minimize in range again

distribution = AddedDistribution(length(spline_x.outlierIndices)/length(spline_x.t),outlierDistribution,noiseDistribution);
spline_x = SmoothingSpline(t_data,x_data,distribution,'lambda',Lambda.fullTensionExpected,'t_knot',t_knot_x);

distribution = AddedDistribution(length(spline_y.outlierIndices)/length(spline_y.t),outlierDistribution,noiseDistribution);
spline_y = SmoothingSpline(t_data,y_data,distribution,'lambda',Lambda.fullTensionExpected,'t_knot',t_knot_y);

spline_x.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
spline_y.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );

spline_x.outlierIndices = find(abs(spline_x.epsilon) > outlierThreshold);
spline_y.outlierIndices = find(abs(spline_y.epsilon) > outlierThreshold);

fprintf('ranged reduced iterated lambda: (%.3g, %.3g)\n', spline_x.lambda,spline_y.lambda );
fprintf('ranged reduced iterated a_rms: (%.3g, %.3g)\n', std(spline_x.uniqueValuesAtHighestDerivative), std(spline_y.uniqueValuesAtHighestDerivative) );
fprintf('ranged reduced iterated a_rms: (%.3g, %.3g)\n', SmoothingSpline.StandardDeviationAndMeanOfInterquartileRange(spline_x.uniqueValuesAtHighestDerivative), SmoothingSpline.StandardDeviationAndMeanOfInterquartileRange(spline_y.uniqueValuesAtHighestDerivative) );

subplot(sp1)
plot(tq,spline_x(tq))
subplot(sp2)
plot(tq,spline_y(tq))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now use maximum likelihood from the IQR of acceleration

% a = cat(1,spline_x.uniqueValuesAtHighestDerivative,spline_y.uniqueValuesAtHighestDerivative);
% tensionDistribution = NormalDistribution(SmoothingSpline.StandardDeviationOfInterquartileRange(a));

distribution = AddedDistribution(length(spline_x.outlierIndices)/length(spline_x.t),outlierDistribution,noiseDistribution);
tensionDistribution = NormalDistribution(SmoothingSpline.StandardDeviationOfInterquartileRange(spline_x.uniqueValuesAtHighestDerivative));
logLikelihood = @(spline) -sum(distribution.logPDF( spline.epsilon ) ) - sum(tensionDistribution.logPDF(spline.uniqueValuesAtHighestDerivative));
spline_x.minimize( logLikelihood );

distribution = AddedDistribution(length(spline_y.outlierIndices)/length(spline_y.t),outlierDistribution,noiseDistribution);
tensionDistribution = NormalDistribution(SmoothingSpline.StandardDeviationOfInterquartileRange(spline_y.uniqueValuesAtHighestDerivative));
logLikelihood = @(spline) -sum(distribution.logPDF( spline.epsilon ) ) - sum(tensionDistribution.logPDF(spline.uniqueValuesAtHighestDerivative));
spline_y.minimize( logLikelihood );

fprintf('log-likelihood lambda: (%.3g, %.3g)\n', spline_x.lambda,spline_y.lambda );
fprintf('log-likelihood a_rms: (%.3g, %.3g)\n', std(spline_x.uniqueValuesAtHighestDerivative), std(spline_y.uniqueValuesAtHighestDerivative) );
fprintf('log-likelihood a_rms: (%.3g, %.3g)\n', SmoothingSpline.StandardDeviationAndMeanOfInterquartileRange(spline_x.uniqueValuesAtHighestDerivative), SmoothingSpline.StandardDeviationAndMeanOfInterquartileRange(spline_y.uniqueValuesAtHighestDerivative) );

spline_x.outlierIndices = find(abs(spline_x.epsilon) > zmax);
spline_y.outlierIndices = find(abs(spline_y.epsilon) > zmax);

fprintf('identified outliers: (%d, %d)\n',length(spline_x.outlierIndices), length(spline_y.outlierIndices));

subplot(sp1)
scatter(t_data(spline_x.outlierIndices),x_data(spline_x.outlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(t_data(spline_x.outlierIndices),x_data(spline_x.outlierIndices),(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline_x(tq))
subplot(sp2)
scatter(t_data(spline_y.outlierIndices),y_data(spline_y.outlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(t_data(spline_y.outlierIndices),y_data(spline_y.outlierIndices),(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline_y(tq))

legend('data','expected mse n/o','expected mse','ranged mse', 'ranged reduced', 'outliers','data', 'log-likelihood')

return

% spline_x.minimize( @(spline) spline.expectedMeanSquareErrorFromGCV );
% spline_y.minimize( @(spline) spline.expectedMeanSquareErrorFromGCV );

a_rms = 6.1e-09;
tensionDistribution = NormalDistribution(a_rms*1.6);
logLikelihood = @(spline) -sum(noiseDistribution.logPDF( spline.epsilon ) ) - sum(tensionDistribution.logPDF(spline.uniqueValuesAtHighestDerivative));

spline_x.minimize( logLikelihood );
spline_y.minimize( logLikelihood );

% spline_x.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
% spline_y.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );

fprintf('(lambda_x, lambda_y)=(%g, %g)\n',spline_x.lambda, spline_y.lambda);
fprintf('(n_eff_x, n_eff_y)=(%.2f, %.2f)\n',spline_x.effectiveSampleSizeFromVarianceOfTheMean, spline_y.effectiveSampleSizeFromVarianceOfTheMean);
fprintf('(a_std_x, a_std_y)=(%g, %g)\n',std(spline_x.uniqueValuesAtHighestDerivative),std(spline_y.uniqueValuesAtHighestDerivative));



% spline.minimize( @(spline) noiseDistribution.kolmogorovSmirnovError(spline.epsilon,zmin,zmax) )
% spline.minimize( @(spline) noiseDistribution.kolmogorovSmirnovError(spline.epsilon,zmin,zmax) )

spline_x.outlierIndices = find(spline_x.epsilon < zmin | spline_x.epsilon > zmax);
spline_x.nonOutlierIndices = setdiff((1:length(spline_x.x))',spline_x.outlierIndices);

spline_y.outlierIndices = find(spline_y.epsilon < zmin | spline_y.epsilon > zmax);
spline_y.nonOutlierIndices = setdiff((1:length(spline_y.x))',spline_y.outlierIndices);

figure
subplot(2,1,1)
scatter(t_data(spline_y.outlierIndices),y_data(spline_y.outlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(t_data,y_data,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline_y(tq),'r')

subplot(2,1,2)
scatter(t_data(spline_x.outlierIndices),x_data(spline_x.outlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(t_data,x_data,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline_x(tq),'r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

t_knot = InterpolatingSpline.KnotPointsForPoints(t_data(spline.nonOutlierIndices),spline.K,1);
spline2 = SmoothingSpline(t_data,y_data,distribution,'lambda',Lambda.optimalIterated,'t_knot',t_knot);
spline2.minimizeExpectedMeanSquareErrorInPercentileRange(pct/2,1-pct/2);
plot(tq,spline2(tq),'b')

return



% splinefit = SmoothingSpline(t_data,x_data,10);
% 
% gpsfit = GPSSmoothingSpline(t_data,x_data,y_data, 'shouldIdentifyOutliers', 0);
% 
% t_tension = gpsfit.spline_x.t_knot(gpsfit.spline_x.K:1:end-gpsfit.spline_x.K+1);
% t_tension = t_tension(1:end-1) + diff(t_tension)/2;
% [x_T,I] = sort(gpsfit.spline_x(t_tension,3));
% x_T_Q1 = median(x_T(1:floor(length(x_T)/2)));
% x_T_Q3 = median(x_T(ceil(length(x_T)/2)+1:end));
% sigma_T = (x_T_Q3-x_T_Q1)/1.349;
% p_x_T = 1-erf(abs(x_T)/sigma_T/sqrt(2));
% T_x_T = t_tension(I);
% 
% outlierIndices = p_x_T < 0.001;
% 
% % for odd K the knots are between the points, even they're on the points
% 
% figure, histogram(x_T,100), vlines([x_T_Q1,x_T_Q3],'k')
% 
% 
% [x,y] = gpsfit(t);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Position fit figure
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% t2 = sort(T_x_T(outlierIndices));
% 
% figure
% s = 1/1000; % scale
% plot(s*x,s*y, 'LineWidth', 0.5*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]), hold on
% % scatter(s*x_data(gpsfit.outlierIndices),s*y_data(gpsfit.outlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
% scatter(s*gpsfit.spline_x(t2),s*gpsfit.spline_y(t2),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
% scatter(s*x_data,s*y_data,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Repeat, but now with student-t distribution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma = 8.5;
nu = 4.5;
variance_of_the_noise = sigma*sigma*nu/(nu-2);

nu_noise = 4.5;
sigma_noise = 40;

w_gps = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
w_outlier = @(z)((nu_noise/(nu_noise+1))*sigma_noise^2*(1+z.^2/(nu_noise*sigma_noise^2)));

% epsilon/sigma_g^2 = epsilon/sigma_t (1-x^2/c^2)^2;
sigma_g = 10; c = 2.1;
psi_tukey = @(z) (z/sigma_g^2)*((1-(z/(sigma_g*c)).^2).^2) .* (abs(z)<c*sigma_g);
w_tukey = @(z) (sigma_g^2*(1-(z/(sigma_g*c)).^2).^(-2)).*(abs(z)<c*sigma_g) + 4e5*(abs(z)>=c*sigma_g);

z = linspace(-50,50,1000)';
figure, plot(z,w_gps(z))
hold on, plot(z,w_tukey(z))
return
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));

% splinefit_t = SmoothingSpline(t_data, y_data, sqrt(variance_of_the_noise), 'K', 4, 'weightFunction', w, 'lambda', Lambda.fullTensionExpected);
% [x_T, t_T] = splinefit_t.UniqueValuesAtHighestDerivative();
% 
% [x_T,I] = sort(x_T);
% x_T_Q1 = median(x_T(1:floor(length(x_T)/2)));
% x_T_Q3 = median(x_T(ceil(length(x_T)/2)+1:end));
% sigma_T = (x_T_Q3-x_T_Q1)/1.349;
% 
% penaltyFunction = @(spline) KolmogorovSmirnovErrorForGaussian( spline.UniqueValuesAtHighestDerivative(), sigma_T );
% splinefit_t.Minimize(penaltyFunction);
% 
% figure
% plot(t/3600,s*splinefit_t(t)), hold on
% scatter(t_data/3600,s*y_data)
% 
% figure
% plot(t/3600,s*splinefit_t(t,splinefit_t.K-1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Repeat, but now with autocorrelation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho = @(dt) exp(-abs(dt)/550.);

spline_x = SmoothingSpline(t_data, x_data, sqrt(variance_of_the_noise), 'K', 4, 'weightFunction', w, 'lambda', Lambda.fullTensionExpected,'autocorrelationFunction',rho);
spline_y = SmoothingSpline(t_data, y_data, sqrt(variance_of_the_noise), 'K', 4, 'weightFunction', w, 'lambda', Lambda.fullTensionExpected,'autocorrelationFunction',rho);

x_T = spline_x.UniqueValuesAtHighestDerivative();
y_T = spline_y.UniqueValuesAtHighestDerivative();

sigma_Tx = SmoothingSpline.StandardDeviationOfInterquartileRange(x_T);
sigma_Ty = SmoothingSpline.StandardDeviationOfInterquartileRange(y_T);
sigma_T = SmoothingSpline.StandardDeviationOfInterquartileRange(cat(1,x_T,y_T));

sigma_T = 8e-9;

figure
histogram(x_T,50), hold on
histogram(y_T,50)

cdf = @(r) 1-exp(-r.*r/(2*sigma_T*sigma_T));
pdf = @(r) r.*exp(-r.*r/(2*sigma_T*sigma_T))/(sigma_T*sigma_T);

r_T = sqrt(x_T.^2 + y_T.^2);

figure, histogram(r_T,100,'Normalization','pdf')
r = linspace(0,max(r_T))';
hold on, plot(r,pdf(r))


% sigma_T = 2.9e-9;

% penaltyFunction = @(spline) LogLikelihoodOfTandG(spline,sigma,nu,sigma_T);
% splinefit_tac_x.Minimize(penaltyFunction);
% splinefit_tac_y.Minimize(penaltyFunction);

penaltyFunction = @(spline1,spline2) GPSSmoothingSpline.LogLikelihood2D(spline1,spline2,sigma,nu,sigma_T);
GPSSmoothingSpline.MinimizeFunctionOfTwoSplines(spline_x,spline_y,penaltyFunction);

figure
s = 1/1000; % scale
plot(s*spline_x(t),s*spline_y(t), 'LineWidth', 0.5*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]), hold on
scatter(s*x_data,s*y_data,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
axis equal

return

[x_T,t_T] = spline_x.UniqueValuesAtHighestDerivative();
y_T = spline_y.UniqueValuesAtHighestDerivative();

figure
histogram(x_T,50), hold on
histogram(y_T,50)

cdf = @(r) 1-exp(-r.*r/(2*sigma_T*sigma_T));
pdf = @(r) r.*exp(-r.*r/(2*sigma_T*sigma_T))/(sigma_T*sigma_T);

r_T = sqrt(x_T.^2 + y_T.^2);
totalOutliers = sum( (1-cdf(r_T)) < 0.01)

figure, histogram(r_T,100,'Normalization','pdf')
r = linspace(0,max(r_T))';
hold on, plot(r,pdf(r))

return

penaltyFunction = @(spline) KolmogorovSmirnovErrorForGaussian( spline.UniqueValuesAtHighestDerivative(), sigma_T );
splinefit_tac_x.Minimize(penaltyFunction);

penaltyFunction = @(spline) KolmogorovSmirnovErrorForGaussian( spline.UniqueValuesAtHighestDerivative(), sigma_T );
splinefit_tac_y.Minimize(penaltyFunction);



cdf = @(r) 1-exp(-r.*r/(2*sigma_T*sigma_T));
pdf = @(r) r.*exp(-r.*r/(2*sigma_T*sigma_T))/(sigma_T*sigma_T);

splinefit_tac_x.lambda = 22e14;
splinefit_tac_y.lambda = splinefit_tac_x.lambda;

[x_T,t_T] = splinefit_tac_x.UniqueValuesAtHighestDerivative();
y_T = splinefit_tac_y.UniqueValuesAtHighestDerivative();
r_T = sqrt(x_T.^2 + y_T.^2);
totalOutliers = sum( (1-cdf(r_T)) < 0.01)

figure, histogram(r_T,100,'Normalization','pdf')
r = linspace(0,max(r_T))';
hold on, plot(r,pdf(r))

[r, ~, cdf] = GPSSmoothingSpline.TwoDimStudentTProbabilityDistributionFunction(sigma, nu);

distanceError = sqrt(splinefit_tac_x.epsilon.^2 + splinefit_tac_y.epsilon.^2);
outlierCut = interp1(cdf,r,1-0.01,'spline');
outliers = distanceError > outlierCut;

totalOutliers = sum(outliers)

figure
histogram(x_T), hold on
histogram(y_T)

figure
s = 1/1000; % scale
plot(s*splinefit_tac_x(t),s*splinefit_tac_y(t), 'LineWidth', 0.5*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]), hold on
% scatter(s*x_data(gpsfit.outlierIndices),s*y_data(gpsfit.outlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(s*x_data(outliers),s*y_data(outliers),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(s*x_data,s*y_data,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')

% figure
% plot(t/3600,s*splinefit_tac_x(t)), hold on
% scatter(t_data/3600,s*y_data)
% 
% figure
% plot(t/3600,s*splinefit_tac_x(t,splinefit_tac_x.K-1))

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
scatter(drifters.t{choiceDrifter}(gpsfit.outlierIndices)/3600,s*drifters.x{choiceDrifter}(gpsfit.outlierIndices),(6.5*scaleFactor)^2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
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
scatter(drifters.t{choiceDrifter}(gpsfit.outlierIndices)/3600,s*drifters.y{choiceDrifter}(gpsfit.outlierIndices),(6.5*scaleFactor)^2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(drifters.t{choiceDrifter}/3600,s*drifters.y{choiceDrifter},(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('t (hours)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)

xlim([124 149])