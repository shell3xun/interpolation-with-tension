scaleFactor = 1;
LoadFigureDefaults

shouldUseStudentTDistribution = 1;

percentOutliers = 0.0;
outlierFactor = 50;

slopes = [-2; -3; -4];
slopes = -4;
totalSlopes = length(slopes);

S = 2;
T = S;
K = S+1;

result_stride = 2*[1;4;16;64];
result_stride = 150;
totalStrides = length(result_stride);

% matern signal parameters
sigma_u = 0.20;
base_dt = 5; % for whatever reason, we chose this as the primary dt
t_damp = 30*60;
n = 250;

totalEnsembles = 1;

for iSlope = 1:length(slopes)
    slope = slopes(iSlope);
    fprintf('slope %d, ',slope);
    
    for iStride=1:length(result_stride)
        stride = result_stride(iStride);
        dt = stride*base_dt;
        
        for iEnsemble = 1:totalEnsembles
            fprintf('..%d',iEnsemble);
            
            cv=maternoise(dt,n,sigma_u*sqrt(2),abs(slope),1/t_damp);
            cx = cumtrapz(cv)*dt;
            data = struct('t',dt*(0:n-1)','x',real(cx));
            
            if shouldUseStudentTDistribution == 1
                nu = 4.5; sigma =  8.5;        
                noiseDistribution = StudentTDistribution(sigma,nu);
                outlierDistribution = StudentTDistribution(outlierFactor*sigma,3.0);
            else
                sigma = 10;
                noiseDistribution = NormalDistribution(sigma);
                outlierDistribution = NormalDistribution(outlierFactor*sigma);
            end
            outlierIndices = rand(n,1)<=percentOutliers;
            outlierIndices([1 n]) = 0;
            
            epsilon_x = zeros(n,1);
            epsilon_x_outlier = zeros(n,1);
            
            epsilon_x(~outlierIndices) = noiseDistribution.rand(sum(~outlierIndices));
            epsilon_x_outlier(outlierIndices) = outlierDistribution.rand(sum(outlierIndices));
            
            epsilon = epsilon_x + epsilon_x_outlier;
            
            % One measure of outlier
            outlierThreshold = noiseDistribution.locationOfCDFPercentile(1-1/10000/2);
            trueOutliers = find(abs(epsilon) > outlierThreshold);
            nonOutlierIndices = setdiff(1:n,trueOutliers);
             
            f = @(z) abs( (1-percentOutliers)*noiseDistribution.pdf(z) - percentOutliers*outlierDistribution.pdf(z) );
            z_crossover = abs(fminsearch(f,sqrt(noiseDistribution.variance)));
            crossoverOutliers = find(abs(epsilon) > outlierThreshold);
            crossovernonOutlierIndices = setdiff(1:n,crossoverOutliers);
            
            x_obs = data.x + epsilon;
            t_obs = data.t;  
        end
        fprintf('\n');
    end
    
end
trueAddedDistribution = AddedDistribution(percentOutliers,outlierDistribution,noiseDistribution);
fprintf('sqrt(variance) (noise, outlier)=(%.1f, %.1f)\n',sqrt(noiseDistribution.variance),sqrt(outlierDistribution.variance));
fprintf('sqrt(sample variance) (noise, outlier)=(%.1f, %.1f)\n',sqrt(mean(epsilon(~outlierIndices).^2)),sqrt(mean(epsilon(outlierIndices).^2)));
fprintf('true alpha: %.2f\n',sum(outlierIndices)/n);

% to be used for ranged minimization when blind
shouldMinimizeBlind = 1;
f = @(z) abs( (percentOutliers/(1-percentOutliers))*outlierDistribution.cdf(-abs(z))/noiseDistribution.cdf(-abs(z)) - 200);
z_crossover = abs(fminsearch(f,sqrt(noiseDistribution.variance)));

tq = linspace(min(t_obs),max(t_obs),10*length(t_obs));
% [D,tD] = SmoothingSpline.FiniteDifferenceMatrixNoBoundary(2,data.t(indices),2);
% spline_clean = SmoothingSpline(data.t(indices),data.x(indices),noiseDistribution, 'S', S, 'lambda',0);
% 
% figure
% plot(tq,spline_clean(tq,2)), hold on
% plot(tD,D*data.x(indices))

dt = t_obs(2)-t_obs(1);
a = diff(diff(data.x))/(dt^2);
fprintf('true strided acceleration: %.3g\n', sqrt( mean( a.*a ) ));

compute_ms_error = @(spline) (mean(mean(  (data.x - spline(data.t)).^2,2 ),1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Toss out the points that are outliers, and perform the best fit.

spline_optimal = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',0/10000);
if shouldMinimizeBlind == 1
    spline_optimal.minimizeExpectedMeanSquareError();
    %     spline_optimal.minimize( @(spline) spline.expectedMeanSquareErrorInterquartile() );
    %     spline_optimal.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-z_crossover,z_crossover) );
else
    spline_optimal.minimizeMeanSquareError(data.t,data.x);
end

baseline_optimal_mse = compute_ms_error(spline_optimal);
fprintf('optimal mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', baseline_optimal_mse, std(spline_optimal.uniqueValuesAtHighestDerivative),spline_optimal.lambda );
falseNegativeOutliers = setdiff(trueOutliers,spline_optimal.outlierIndices);
falsePositiveOutliers = setdiff(spline_optimal.outlierIndices,trueOutliers);
fprintf('False positives: %d, negatives: %d\n',length(falsePositiveOutliers),length(falseNegativeOutliers))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Try our old method

% spline_robust_old = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
% spline_robust_old.rebuildOutlierDistributionAndAdjustWeightingsOldMethod();
% spline_robust_old.minimizeMeanSquareError(data.t,data.x);
% 
% fprintf('robust mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', compute_ms_error(spline_robust_old), std(spline_robust_old.uniqueValuesAtHighestDerivative),spline_robust_old.lambda );
% falseNegativeOutliers = setdiff(trueOutliers,spline_robust_old.outlierIndices);
% falsePositiveOutliers = setdiff(spline_robust_old.outlierIndices,trueOutliers);
% fprintf('False positives: %d, negatives: %d\n',length(falsePositiveOutliers),length(falseNegativeOutliers))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Try our 1 pass optimal *unblinded* method

% spline_robust_old = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
% lambda_opt = spline_robust_old.minimize( @(spline) trueAddedDistribution.andersonDarlingError(spline.epsilon) );
% zmax = max(abs(spline_robust_old.epsilonAtIndices(crossovernonOutlierIndices)));
% zmin = min(abs(spline_robust_old.epsilonAtIndices(crossoverOutliers)));
% 
% z_noise = sort(abs(spline_robust_old.epsilonAtIndices(crossovernonOutlierIndices)));
% noise_cdf = 1-(1:length(z_noise))'/length(z_noise);
% z_outlier = sort(abs(spline_robust_old.epsilonAtIndices(crossoverOutliers)));
% outlier_cdf = (1:length(z_outlier))'/length(z_outlier);
% 
% figure
% subplot(1,2,1)
% plot(z_noise,noise_cdf), hold on
% plot(z_outlier,outlier_cdf), xlim([0 max(z_noise)])
% 
% spline_robust_old.rebuildOutlierDistributionAndAdjustWeightings();
% fprintf('zmax: %.1f, pdfs: noise likelihood %.2g\n',zmax,noiseDistribution.cdf(-zmax));
% fprintf('zmin: %.1f, pdfs: noise likelihood %.2g\n',zmin,noiseDistribution.cdf(-zmin));
% 
% spline_robust_old.minimizeMeanSquareError(data.t,data.x);
% epsilon_min = spline_robust_old.epsilon;
% 
% fprintf('robust mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', compute_ms_error(spline_robust_old), std(spline_robust_old.uniqueValuesAtHighestDerivative),spline_robust_old.lambda );
% falseNegativeOutliers = setdiff(trueOutliers,spline_robust_old.outlierIndices);
% falsePositiveOutliers = setdiff(spline_robust_old.outlierIndices,trueOutliers);
% fprintf('False positives: %d, negatives: %d\n',length(falsePositiveOutliers),length(falseNegativeOutliers))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Try our 1 pass optimal method

RobustSmoothingSpline.estimateOutlierDistributionFromKnownNoise(epsilon,noiseDistribution);

spline_robust = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',0/10000);
spline_robust.setToFullTensionWithIteratedIQAD();
epsilon_full = spline_robust.epsilon;
lambda_full = spline_robust.lambda;

[empiricalOutlierDistribution,empiricalAlpha] = spline_robust.estimateOutlierDistribution();
fprintf('(alpha,sqrt(var))=(%.2f,%.1f m)\n',empiricalAlpha,sqrt(empiricalOutlierDistribution.variance));
newAddedDistribution = AddedDistribution(empiricalAlpha,empiricalOutlierDistribution,noiseDistribution);

outlierOddsArray = 10.^linspace(log10(10),log10(1e8),15)';
% outlierOddsArray = cat(1,outlierOddsArray,1e14);
z_mse = zeros(size(outlierOddsArray));
mse = zeros(size(outlierOddsArray));
lambdas = zeros(size(outlierOddsArray));
for iOdds = 1:length(outlierOddsArray)
    outlierOdds = outlierOddsArray(iOdds);
    f = @(z) abs( (empiricalAlpha/(1-empiricalAlpha))*empiricalOutlierDistribution.cdf(-abs(z))/noiseDistribution.cdf(-abs(z)) - outlierOdds);
    z_mse(iOdds) = fminsearch(f,sqrt(noiseDistribution.variance));
    
    spline_robust.distribution = newAddedDistribution;
    spline_robust.outlierDistribution = empiricalOutlierDistribution;
    spline_robust.alpha = empiricalAlpha;
    
    noiseIndices = abs(epsilon_full) <= z_mse(iOdds);
    spline_robust.sigma = noiseIndices .* sqrt(noiseDistribution.variance) + (~noiseIndices) .* sqrt(empiricalOutlierDistribution.variance);
    spline_robust.distribution.w = @(z) noiseIndices .* noiseDistribution.w(z) + (~noiseIndices) .* empiricalOutlierDistribution.w(z);
    spline_robust.lambda = lambda_full;

    if shouldMinimizeBlind == 1
        spline_robust.minimizeExpectedMeanSquareError();
%         spline_optimal.minimize( @(spline) spline.expectedMeanSquareErrorInterquartile() );
%         spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-z_crossover,z_crossover) );
    else
        spline_robust.minimizeMeanSquareError(data.t,data.x);
    end
    
    lambdas(iOdds) = spline_robust.lambda;
    mse(iOdds) = compute_ms_error(spline_robust);
end

figure
subplot(1,2,1)
plot(z_mse,mse), hold on
plot(z_mse,baseline_optimal_mse*ones(size(z_mse)));
xlog
ylim([0.8*min(min(mse),baseline_optimal_mse) 2*baseline_optimal_mse])
subplot(1,2,2)
plot(outlierOddsArray,mse), hold on
plot(outlierOddsArray,baseline_optimal_mse*ones(size(z_mse)));
xlog
ylim([0.8*min(min(mse),baseline_optimal_mse) 2*baseline_optimal_mse])

[minMSE,index] = min(mse,[],'omitnan');

fprintf('optimal mse: %.2f m\n', minMSE );
fprintf('z_opt: %.1f, odds %.2g\n',z_mse(index),outlierOddsArray(index));

noiseIndices = abs(epsilon_full) <= z_mse(index);
spline_robust.distribution = newAddedDistribution;
spline_robust.outlierDistribution = empiricalOutlierDistribution;
spline_robust.alpha = empiricalAlpha;
spline_robust.sigma = noiseIndices .* sqrt(noiseDistribution.variance) + (~noiseIndices) .* sqrt(empiricalOutlierDistribution.variance);
spline_robust.distribution.w = @(z) noiseIndices .* noiseDistribution.w(z) + (~noiseIndices) .* empiricalOutlierDistribution.w(z);
spline_robust.lambda = lambda_full;

if shouldMinimizeBlind == 1
    spline_robust.minimizeExpectedMeanSquareError();
%     spline_optimal.minimize( @(spline) spline.expectedMeanSquareErrorInterquartile() );
%     spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-z_crossover,z_crossover) );
else
    spline_robust.minimizeMeanSquareError(data.t,data.x);
end

fprintf('robust mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', compute_ms_error(spline_robust), std(spline_robust.uniqueValuesAtHighestDerivative),spline_robust.lambda );
falseNegativeOutliers = setdiff(trueOutliers,spline_robust.outlierIndices);
falsePositiveOutliers = setdiff(spline_robust.outlierIndices,trueOutliers);
fprintf('False positives: %d, negatives: %d\n',length(falsePositiveOutliers),length(falseNegativeOutliers))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure

figure
scatter(t_obs(spline_robust.outlierIndices),x_obs(spline_robust.outlierIndices),(8.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'), hold on
scatter(t_obs(trueOutliers),x_obs(trueOutliers),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
scatter(t_obs,x_obs,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline_optimal(tq),'k')
plot(tq,spline_robust(tq),'b')

return



zmax = max(abs(spline_robust.epsilonAtIndices(crossovernonOutlierIndices)));
zmin = min(abs(spline_robust.epsilonAtIndices(crossoverOutliers)));
fprintf('zmax: %.1f, pdfs: noise likelihood %.2g\n',zmax,(empiricalAlpha/(1-empiricalAlpha))*empiricalOutlierDistribution.cdf(-zmax)/noiseDistribution.cdf(-zmax));
fprintf('zmin: %.1f, pdfs: noise likelihood %.2g\n',zmin,(empiricalAlpha/(1-empiricalAlpha))*empiricalOutlierDistribution.cdf(-zmin)/noiseDistribution.cdf(-zmin));

z_noise = sort(abs(spline_robust.epsilonAtIndices(crossovernonOutlierIndices)));
noise_cdf = 1-(1:length(z_noise))'/length(z_noise);
z_outlier = sort(abs(spline_robust.epsilonAtIndices(crossoverOutliers)));
outlier_cdf = (1:length(z_outlier))'/length(z_outlier);

z = linspace(min(z_noise),max(z_noise),100)';
cd2 = interp1(z_noise,noise_cdf,z)+interp1(z_outlier,outlier_cdf,z);


%%%%%%%%%%
% Let actually compute the mse at different crossover points
% z_mse = linspace(min(zmax,zmin),max(zmax,zmin),10)';
z_mse = z_outlier( z_outlier < zmax );
cat(1,z_mse,zmax);
if isempty(z_mse)
   z_mse = zmin + (zmax-zmin)/2; 
end
mse = zeros(size(z_mse));



figure
subplot(1,2,1)
% figure
plot(z_noise,noise_cdf), hold on
plot(z_outlier,outlier_cdf), xlim([0 max(z_mse)])
plot(z,cd2)
[~,index] = min(cd2,[],'omitnan');
fprintf('z_opt: %.1f, pdfs: noise likelihood %.2g\n',zmin,(empiricalAlpha/(1-empiricalAlpha))*empiricalOutlierDistribution.cdf(-z(index))/noiseDistribution.cdf(-z(index)));

vlines(z(index),'g--')
mse0 = max(mse);
plot(z_mse,mse/mse0,'k')
[minMSE,index] = min(mse,[],'omitnan');

fprintf('optimal mse: %.2f m\n', minMSE );
fprintf('z_opt: %.1f, pdfs: noise likelihood %.2g\n',zmin,(empiricalAlpha/(1-empiricalAlpha))*empiricalOutlierDistribution.cdf(-z_mse(index))/noiseDistribution.cdf(-z_mse(index)));

%%%%%%%%%%
% figure out how well the two distributions are split
spline_robust = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
% spline_robust.setInterquartileSampleVarianceToExpectedValue(alpha,2.0*noiseDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
spline_robust.setToFullTensionWithIteratedIQSV(1.0);
lambda_full = spline_robust.lambda
[empiricalOutlierDistribution,empiricalAlpha] = spline_robust.estimateOutlierDistribution();
fprintf('(alpha,sqrt(var))=(%.2f,%.1f m)\n',empiricalAlpha,sqrt(empiricalOutlierDistribution.variance));
newAddedDistribution = AddedDistribution(empiricalAlpha,empiricalOutlierDistribution,noiseDistribution);

%%%%%%%%%%
% figure out how well the two distributions are split
zmax = max(abs(spline_robust.epsilonAtIndices(crossovernonOutlierIndices)));
zmin = min(abs(spline_robust.epsilonAtIndices(crossoverOutliers)));
fprintf('zmax: %.1f, pdfs: noise likelihood %.2g\n',zmax,(empiricalAlpha/(1-empiricalAlpha))*empiricalOutlierDistribution.cdf(-zmax)/noiseDistribution.cdf(-zmax));
fprintf('zmin: %.1f, pdfs: noise likelihood %.2g\n',zmin,(empiricalAlpha/(1-empiricalAlpha))*empiricalOutlierDistribution.cdf(-zmin)/noiseDistribution.cdf(-zmin));

z_noise = sort(abs(spline_robust.epsilonAtIndices(crossovernonOutlierIndices)));
noise_cdf = 1-(1:length(z_noise))'/length(z_noise);
z_outlier = sort(abs(spline_robust.epsilonAtIndices(crossoverOutliers)));
outlier_cdf = (1:length(z_outlier))'/length(z_outlier);

z = linspace(min(z_noise),max(z_noise),100)';
cd2 = interp1(z_noise,noise_cdf,z)+interp1(z_outlier,outlier_cdf,z);

%%%%%%%%%%
% Let actually compute the mse at different crossover points
% z_mse = linspace(min(zmax,zmin),max(zmax,zmin),10)';
z_mse = z_outlier( z_outlier < zmax );
cat(1,z_mse,zmax);
if isempty(z_mse)
   z_mse = zmin + (zmax-zmin)/2; 
end
mse = zeros(size(z_mse));

for iZ = 1:length(z_mse)
    z_crossover = z_mse(iZ);
    spline_robust.distribution = newAddedDistribution;
    epsilonSpline = spline_robust.epsilon;
    noiseIndices = epsilonSpline > -z_crossover & epsilonSpline < z_crossover;
    spline_robust.distribution.w = @(z) noiseIndices .* noiseDistribution.w(z) + (~noiseIndices) .* outlierDistribution.w(z);
    spline_robust.lambda = lambda_full;
    spline_robust.minimizeMeanSquareError(data.t,data.x);
    mse(iZ) = compute_ms_error(spline_robust);
end



subplot(1,2,2)
% figure
plot(z_noise,noise_cdf), hold on
plot(z_outlier,outlier_cdf), xlim([0 max(z_mse)])
plot(z,cd2)
[~,index] = min(cd2,[],'omitnan');
vlines(z(index),'g--')
fprintf('z_opt: %.1f, pdfs: noise likelihood %.2g\n',zmin,(empiricalAlpha/(1-empiricalAlpha))*empiricalOutlierDistribution.cdf(-z(index))/noiseDistribution.cdf(-z(index)));
plot(z_mse,mse/mse0,'k')
[minMSE,index] = min(mse,[],'omitnan');
fprintf('optimal mse: %.2f m\n', minMSE );
fprintf('z_opt: %.1f, pdfs: noise likelihood %.2g\n',z_mse(index),(empiricalAlpha/(1-empiricalAlpha))*empiricalOutlierDistribution.cdf(-z_mse(index))/noiseDistribution.cdf(-z_mse(index)));

spline_robust = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
z_crossover = z_mse(index);
spline_robust.distribution = newAddedDistribution;
epsilonSpline = spline_robust.epsilon;
noiseIndices = epsilonSpline > -z_crossover & epsilonSpline < z_crossover;
spline_robust.distribution.w = @(z) noiseIndices .* noiseDistribution.w(z) + (~noiseIndices) .* outlierDistribution.w(z);
spline_robust.lambda = lambda_full;
spline_robust.minimizeMeanSquareError(data.t,data.x);

% plot(z,100*cd2_alt)

% fprintf('lambda_opt = %.1g, lambda actual: %.1g\n', lambda_opt, spline_robust.lambda);
% spline_robust.rebuildOutlierDistributionAndAdjustWeightings();
% spline_robust.minimizeMeanSquareError(data.t,data.x);

fprintf('robust mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', compute_ms_error(spline_robust), std(spline_robust.uniqueValuesAtHighestDerivative),spline_robust.lambda );
falseNegativeOutliers = setdiff(trueOutliers,spline_robust.outlierIndices);
falsePositiveOutliers = setdiff(spline_robust.outlierIndices,trueOutliers);
fprintf('False positives: %d, negatives: %d\n',length(falsePositiveOutliers),length(falseNegativeOutliers))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Try our 2 pass optimal method

% alpha = 0.5;
% spline_robust2x = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
% spline_robust2x.setInterquartileSampleVarianceToExpectedValue(alpha,noiseDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
% [outlierDistribution, alpha_outliers] = spline_robust2x.estimateOutlierDistribution();
% addedDistribution = AddedDistribution(alpha_outliers,outlierDistribution,noiseDistribution);
% spline_robust2x.setInterquartileSampleVarianceToExpectedValue(alpha,addedDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
% 
% spline_robust2x.rebuildOutlierDistributionAndAdjustWeightings();
% spline_robust2x.minimizeMeanSquareError(data.t,data.x);
% 
% fprintf('robust mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', compute_ms_error(spline_robust2x), std(spline_robust2x.uniqueValuesAtHighestDerivative),spline_robust2x.lambda );
% falseNegativeOutliers = setdiff(trueOutliers,spline_robust2x.outlierIndices);
% falsePositiveOutliers = setdiff(spline_robust2x.outlierIndices,trueOutliers);
% fprintf('False positives: %d, negatives: %d\n',length(falsePositiveOutliers),length(falseNegativeOutliers))




% plot(tq,spline_robust2x(tq),'r')
% plot(tq,spline_robust_old(tq),'g')


