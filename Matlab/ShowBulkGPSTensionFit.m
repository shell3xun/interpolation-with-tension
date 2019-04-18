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


pct = 0.10;
distribution = AddedDistribution(pct,NormalDistribution(400),StudentTDistribution(8.5,4.5));

pct = 0.100;
distribution = AddedDistribution(pct,StudentTDistribution(300,3.0),StudentTDistribution(8.5,4.5));

nDrifters = length(drifters.x);
splines_x = cell(nDrifters,1);
splines_y = cell(nDrifters,1);
for iDrifter = 1:nDrifters
    splines_x{iDrifter} = SmoothingSpline(drifters.t{iDrifter},drifters.x{iDrifter},distribution,'lambda',Lambda.fullTensionExpected);
    splines_y{iDrifter} = SmoothingSpline(drifters.t{iDrifter},drifters.y{iDrifter},distribution,'lambda',Lambda.fullTensionExpected);
end

zmin = distribution.locationOfCDFPercentile(pct/2);
zmax = distribution.locationOfCDFPercentile(1-pct/2);

% SmoothingSpline.minimizeFunctionOfSplines( splines_x, @(splines) SmoothingSpline.expectedMeanSquareErrorOfSplines(splines,zmin,zmax) );
% SmoothingSpline.minimizeFunctionOfSplines( splines_y, @(splines) SmoothingSpline.expectedMeanSquareErrorOfSplines(splines,zmin,zmax) );

splines = cell(2*nDrifters,1);
splines(1:nDrifters) = splines_x;
splines((nDrifters+1):(2*nDrifters)) = splines_y;
SmoothingSpline.minimizeFunctionOfSplines( splines, @(splines) SmoothingSpline.expectedMeanSquareErrorOfSplines(splines,zmin,zmax) );

epsilon_d = [];
epsilon = [];
for iDrifter = 1:nDrifters
    epsilon = cat(1,epsilon,splines_x{iDrifter}.epsilon,splines_y{iDrifter}.epsilon);
    ep = sqrt(splines_x{iDrifter}.epsilon.^2 + splines_y{iDrifter}.epsilon.^2);
    epsilon_d = cat(1,epsilon_d,ep);
    
    splines_x{iDrifter}.outlierIndices = find(ep > 200)';
    splines_y{iDrifter}.outlierIndices = find(ep > 200)';
    
    splines_x{iDrifter}.nonOutlierIndices = setdiff((1:length(splines_x{iDrifter}.t))',splines_x{iDrifter}.outlierIndices);
    splines_y{iDrifter}.nonOutlierIndices = setdiff((1:length(splines_y{iDrifter}.t))',splines_y{iDrifter}.outlierIndices);
end

z = linspace(0,1000,500);
distDistribution = TwoDimDistanceDistribution(distribution);
figure, histogram(epsilon_d,z,'Normalization','cdf')
hold on, plot(z,distDistribution.cdf(z))

splines_x2 = cell(nDrifters,1);
splines_y2 = cell(nDrifters,1);
for iDrifter = 1:nDrifters
    nonOutlierIndices = splines_x{iDrifter}.nonOutlierIndices;
    if nonOutlierIndices(1) ~= 1
        nonOutlierIndices = cat(1,1,nonOutlierIndices);
    end
    if nonOutlierIndices(end) ~= length(splines_x{iDrifter}.t)
        nonOutlierIndices = cat(1,nonOutlierIndices,length(splines_x{iDrifter}.t));
    end
    
    t_knot = InterpolatingSpline.KnotPointsForPoints(splines_x{iDrifter}.t(nonOutlierIndices),splines_x{iDrifter}.K,1);
    splines_x2{iDrifter} = SmoothingSpline(drifters.t{iDrifter},drifters.x{iDrifter},distribution,'lambda',Lambda.fullTensionExpected,'t_knot',t_knot);
    splines_y2{iDrifter} = SmoothingSpline(drifters.t{iDrifter},drifters.y{iDrifter},distribution,'lambda',Lambda.fullTensionExpected,'t_knot',t_knot);
end

splines = cell(2*nDrifters,1);
splines(1:nDrifters) = splines_x2;
splines((nDrifters+1):(2*nDrifters)) = splines_y2;
SmoothingSpline.minimizeFunctionOfSplines( splines, @(splines) SmoothingSpline.expectedMeanSquareErrorOfSplines(splines,zmin,zmax) );

% splines_no = SmoothingSpline(drifters.t{7},drifters.x{7},StudentTDistribution(8.5,4.5),'lambda',Lambda.optimalIterated);


spline = splines_x2{choiceDrifter};

% spline.outlierIndices = find(spline.epsilon < zmin | spline.epsilon > zmax);
spline.nonOutlierIndices = setdiff((1:length(spline.x))',spline.outlierIndices);

tq=linspace(min(spline.t),max(spline.t),length(spline.t)*10).';

figure
scatter(spline.t(spline.outlierIndices),spline.x(spline.outlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(spline.t,spline.x,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline(tq),'r')

figure
scatter(splines_x2{choiceDrifter}(spline.outlierIndices),splines_y2{choiceDrifter}(spline.outlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(splines_x2{choiceDrifter}.x,splines_y2{choiceDrifter}.x,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(splines_x2{choiceDrifter}(tq),splines_y2{choiceDrifter}(tq),'r')

return

xt = [];
for i=1:length(splines)
    xt = cat(1,xt,splines{i}.uniqueValuesAtHighestDerivative);
end
figure
histogram(xt)

return

% Pull out the data of interest
x_data = drifters.x{choiceDrifter};
y_data = drifters.y{choiceDrifter};
t_data = drifters.t{choiceDrifter};

t=linspace(min(t_data),max(t_data),length(t_data)*10).';


pct = 0.05;
distribution = AddedDistribution(pct,NormalDistribution(800),StudentTDistribution(8.5,4.5));
distribution = AddedDistribution(pct,StudentTDistribution(300,3.0),StudentTDistribution(8.5,4.5));

spline = SmoothingSpline(t_data,y_data,distribution,'lambda',Lambda.optimalIterated);
spline.minimizeExpectedMeanSquareErrorInPercentileRange(pct/2,1-pct/2);

tq = linspace(min(t_data),max(t_data),10*length(t_data));

zmin = spline.distribution.locationOfCDFPercentile(pct/2);
zmax = spline.distribution.locationOfCDFPercentile(1-pct/2);

spline.outlierIndices = find(spline.epsilon < zmin | spline.epsilon > zmax);
spline.nonOutlierIndices = setdiff((1:length(spline.x))',spline.outlierIndices);

figure
% scatter(t_data(spline.outlierIndices),x_data(spline.outlierIndices),(8.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'), hold on
scatter(t_data(spline.outlierIndices),y_data(spline.outlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(t_data,y_data,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline(tq),'r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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