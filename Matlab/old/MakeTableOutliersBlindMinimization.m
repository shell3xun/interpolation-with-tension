%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script generates a signal using the Matern with three different
% slopes. It then adds noise to the signal that include a known
% (quantified) portion and an outlier portion. The script tries a variety
% of added distributions to see which distribution performs the best.
%
% Note that this script uses *unblinded* data to minimize. So here we're
% looking to see which distribution performs best. We later test minimizers
% to see which performs best.
%
% Overall, in terms of the smallest 90th percentile mean-square-error, the
% added distribution with 1/10000 works best.


shouldUseStudentTDistribution = 1;

if shouldUseStudentTDistribution == 1
    filename = 'MSETableOutliersBlindMinimizationStudentT.mat';
else
    filename = 'MSETableOutliersBlindMinimizationNormal.mat';
end

if exist(filename,'file')
    load(filename);
    totalOutlierRatios = size(total_outliers,1);
else
    slopes = [-2; -3; -4];
    slopes = -3;
    totalSlopes = length(slopes);

    strides = [5;20;80;200];
%     strides = 200;
    totalStrides = length(strides);
    totalEnsembles = 21; % best to choose an odd number for median
    
    outlierRatios = [0.0 0.05 0.15 0.25];
    outlierRatios = 0.15;
    totalOutlierRatios = length(outlierRatios);
    
    % spline fit parameters
    S = 2;
    T = S;
    K = S+1;
    
    vars = {'S', 'T', 'slopes', 'strides', 'outlierRatios', 'totalEnsembles'};
    
    % matern signal parameters
    sigma_u = 0.20;
    base_dt = 5; % for whatever reason, we chose this as the primary dt
    t_damp = 30*60;
    n = 250;

    % outlier parameters
   
    outlierFactor = 40;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preallocate the variables we need to save
    %
    nothing = nan(totalOutlierRatios, totalStrides, totalSlopes, totalEnsembles);
    nothing_struct = struct('mse',nothing,'neff_se',nothing,'lambda',nothing,'nonOutlierEffectiveSampleSize',nothing,'nonOutlierSampleVariance',nothing,'false_negatives', nothing, 'false_positives', nothing);
    
    total_outliers = nothing; vars{end+1} = 'total_outliers';
    
    optimal = nothing_struct; vars{end+1} = 'optimal';
    robust_beta2_optimal = nothing_struct; vars{end+1} = 'robust_beta2_optimal';
    robust_beta50_optimal = nothing_struct; vars{end+1} = 'robust_beta50_optimal';
    robust_beta100_optimal = nothing_struct; vars{end+1} = 'robust_beta100_optimal';
    robust_beta200_optimal = nothing_struct; vars{end+1} = 'robust_beta200_optimal';
    robust_beta400_optimal = nothing_struct; vars{end+1} = 'robust_beta400_optimal';
    robust_beta800_optimal = nothing_struct; vars{end+1} = 'robust_beta800_optimal';
    robust_noiseOdds1over10 = nothing_struct; vars{end+1} = 'robust_noiseOdds1over10';
    robust_noiseOdds1over3 = nothing_struct; vars{end+1} = 'robust_noiseOdds1over3';
    robust_noiseOdds1 = nothing_struct; vars{end+1} = 'robust_noiseOdds1';
    robust_noiseOdds2 = nothing_struct; vars{end+1} = 'robust_noiseOdds2';
    robust_noiseOdds3 = nothing_struct; vars{end+1} = 'robust_noiseOdds3';
    robust_noiseOdds4 = nothing_struct; vars{end+1} = 'robust_noiseOdds4';
    robust_noiseOdds10 = nothing_struct; vars{end+1} = 'robust_noiseOdds10';
    robust_varianceCrossover = nothing_struct; vars{end+1} = 'robust_varianceCrossover';
    robust_weighted = nothing_struct; vars{end+1} = 'robust_weighted';
                    
    for iOutlierRatio = 1:totalOutlierRatios
        percentOutliers = outlierRatios(iOutlierRatio);
        fprintf('percentOutliers %.1f, ',100*percentOutliers);
        
        for iSlope = 1:length(slopes)
            slope = slopes(iSlope);
            fprintf('slope %d, ',slope);
            
            for iStride=1:length(strides)
                stride = strides(iStride);
                dt = stride*base_dt;
                fprintf('stride %d, ',stride);
                
                for iEnsemble = 1:totalEnsembles
                    fprintf('..%d',iEnsemble);
                    
                    % Generate the signal
                    cv=maternoise(dt,n,sigma_u*sqrt(2),abs(slope),1/t_damp);
                    cx = cumtrapz(cv)*dt;
                    data = struct('t',dt*(0:n-1)','x',real(cx));
                    
                    compute_ms_error = @(spline) (mean(mean(  (data.x - spline(data.t)).^2,2 ),1));
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Generate the noise
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
                    
                    epsilon_noise = zeros(n,1);
                    epsilon_outlier = zeros(n,1);
                    
                    epsilon_noise(~outlierIndices) = noiseDistribution.rand(sum(~outlierIndices));
                    epsilon_outlier(outlierIndices) = outlierDistribution.rand(sum(outlierIndices));
                    
                    epsilon = epsilon_noise + epsilon_outlier;
                    
                    outlierThreshold = noiseDistribution.locationOfCDFPercentile(1-1/10000/2);
                    trueOutlierIndices = find(abs(epsilon) > outlierThreshold);
                    truenonOutlierIndices = setdiff(1:n,trueOutlierIndices);
                    
                    total_outliers(iOutlierRatio,iStride,iSlope,iEnsemble) = length(trueOutlierIndices);
                    
                    x_obs = data.x + epsilon;
                    t_obs = data.t;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Unblinded best fit with standard tension spline
                    
                    linearIndex = sub2ind(size(nothing),iOutlierRatio,iStride,iSlope,iEnsemble);
                    
                    alpha = 1/10000;
                    spline_robust = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'alpha',alpha);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    optimal = LogStatisticsFromSplineForOutlierTable(optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Blind best fit by minimizing the expected
                    % variance of the noise within different ranges.
                    
                    beta = 1/2;
                    zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                    zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                    robust_beta2_optimal = LogStatisticsFromSplineForOutlierTable(robust_beta2_optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    beta = 1/50;
                    zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                    zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                    robust_beta50_optimal = LogStatisticsFromSplineForOutlierTable(robust_beta50_optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    beta = 1/100;
                    zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                    zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                    robust_beta100_optimal = LogStatisticsFromSplineForOutlierTable(robust_beta100_optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    beta = 1/200;
                    zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                    zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                    robust_beta200_optimal = LogStatisticsFromSplineForOutlierTable(robust_beta200_optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    beta = 1/400;
                    zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                    zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                    robust_beta400_optimal = LogStatisticsFromSplineForOutlierTable(robust_beta400_optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    beta = 1/800;
                    zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                    zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                    robust_beta800_optimal = LogStatisticsFromSplineForOutlierTable(robust_beta800_optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Blind best fit by minimizing the expected
                    % variance of the noise within different cutoffs based
                    % on the empirical fit
                    
                    spline_robust.setToFullTensionWithIteratedIQAD();
                    [empiricalOutlierDistribution,empiricalAlpha] = spline_robust.estimateOutlierDistribution();
                    epsilon_full = spline_robust.epsilon;
                    addedDistribution = AddedDistribution(empiricalAlpha,empiricalOutlierDistribution,noiseDistribution);
                    
                    expectedVarianceDistribution = addedDistribution;
                    
                    if empiricalAlpha == 0 || isempty(empiricalOutlierDistribution)
                        continue;
                    end
                    
                    noiseOdds = 1/10;
                    f = @(z) abs((1-empiricalAlpha)*noiseDistribution.pdf(z)./(empiricalAlpha*empiricalOutlierDistribution.pdf(z))-noiseOdds);
                    zoutlier = abs(fminsearch(f,sqrt(noiseDistribution.variance)));
                    expectedVariance = expectedVarianceDistribution.varianceInRange(-zoutlier,zoutlier);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-zoutlier,zoutlier,expectedVariance) );
                    robust_noiseOdds1over10 = LogStatisticsFromSplineForOutlierTable(robust_noiseOdds1over10,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    noiseOdds = 1/3;
                    f = @(z) abs((1-empiricalAlpha)*noiseDistribution.pdf(z)./(empiricalAlpha*empiricalOutlierDistribution.pdf(z))-noiseOdds);
                    zoutlier = abs(fminsearch(f,sqrt(noiseDistribution.variance)));
                    expectedVariance = expectedVarianceDistribution.varianceInRange(-zoutlier,zoutlier);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-zoutlier,zoutlier,expectedVariance) );
                    robust_noiseOdds1over3 = LogStatisticsFromSplineForOutlierTable(robust_noiseOdds1over3,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);

                    noiseOdds = 1;
                    f = @(z) abs((1-empiricalAlpha)*noiseDistribution.pdf(z)./(empiricalAlpha*empiricalOutlierDistribution.pdf(z))-noiseOdds);
                    zoutlier = abs(fminsearch(f,sqrt(noiseDistribution.variance)));
                    expectedVariance = expectedVarianceDistribution.varianceInRange(-zoutlier,zoutlier);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-zoutlier,zoutlier,expectedVariance) );
                    robust_noiseOdds1 = LogStatisticsFromSplineForOutlierTable(robust_noiseOdds1,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    noiseOdds = 2;
                    f = @(z) abs((1-empiricalAlpha)*noiseDistribution.pdf(z)./(empiricalAlpha*empiricalOutlierDistribution.pdf(z))-noiseOdds);
                    zoutlier = abs(fminsearch(f,sqrt(noiseDistribution.variance)));
                    expectedVariance = expectedVarianceDistribution.varianceInRange(-zoutlier,zoutlier);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-zoutlier,zoutlier,expectedVariance) );
                    robust_noiseOdds2 = LogStatisticsFromSplineForOutlierTable(robust_noiseOdds2,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    noiseOdds = 3;
                    f = @(z) abs((1-empiricalAlpha)*noiseDistribution.pdf(z)./(empiricalAlpha*empiricalOutlierDistribution.pdf(z))-noiseOdds);
                    zoutlier = abs(fminsearch(f,sqrt(noiseDistribution.variance)));
                    expectedVariance = expectedVarianceDistribution.varianceInRange(-zoutlier,zoutlier);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-zoutlier,zoutlier,expectedVariance) );
                    robust_noiseOdds3 = LogStatisticsFromSplineForOutlierTable(robust_noiseOdds3,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    noiseOdds = 4;
                    f = @(z) abs((1-empiricalAlpha)*noiseDistribution.pdf(z)./(empiricalAlpha*empiricalOutlierDistribution.pdf(z))-noiseOdds);
                    zoutlier = abs(fminsearch(f,sqrt(noiseDistribution.variance)));
                    expectedVariance = expectedVarianceDistribution.varianceInRange(-zoutlier,zoutlier);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-zoutlier,zoutlier,expectedVariance) );
                    robust_noiseOdds4 = LogStatisticsFromSplineForOutlierTable(robust_noiseOdds4,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    noiseOdds = 10;
                    f = @(z) abs((1-empiricalAlpha)*noiseDistribution.pdf(z)./(empiricalAlpha*empiricalOutlierDistribution.pdf(z))-noiseOdds);
                    zoutlier = abs(fminsearch(f,sqrt(noiseDistribution.variance)));
                    expectedVariance = expectedVarianceDistribution.varianceInRange(-zoutlier,zoutlier);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-zoutlier,zoutlier,expectedVariance) );
                    robust_noiseOdds10 = LogStatisticsFromSplineForOutlierTable(robust_noiseOdds10,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                                        
                    
                    weakAddedDistribution = AddedDistribution(empiricalAlpha/1.5,StudentTDistribution(sqrt((addedDistribution.variance-(1-empiricalAlpha/1.5)*noiseDistribution.variance)/(3*empiricalAlpha/2)),3.0),noiseDistribution);
                    strongAddedDistribution = AddedDistribution(1.5*empiricalAlpha,StudentTDistribution(sqrt((addedDistribution.variance-(1-empiricalAlpha*1.5)*noiseDistribution.variance)/(3*empiricalAlpha*1.5)),3.0),noiseDistribution);
                    f = @(z) abs(weakAddedDistribution.varianceInRange(0,z)-strongAddedDistribution.varianceInRange(0,z));
                    zvariance_crossover = fminbnd(f,sqrt(noiseDistribution.variance),sqrt(strongAddedDistribution.variance));
                    expectedVariance = expectedVarianceDistribution.varianceInRange(-zvariance_crossover,zvariance_crossover);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-zvariance_crossover,zvariance_crossover,expectedVariance) );
                    robust_varianceCrossover = LogStatisticsFromSplineForOutlierTable(robust_varianceCrossover,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Blind best fit by minimizing the expected mean square
                    % error weighted by probability.
                    
                    spline_robust.w_epsilon = RobustSmoothingSpline.generateEpsilonWeightingFromOutlierDistribution(epsilon_full,empiricalAlpha,empiricalOutlierDistribution,noiseDistribution);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorWithWeighting(noiseDistribution.variance) );
                    robust_weighted = LogStatisticsFromSplineForOutlierTable(robust_weighted,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                end
                fprintf('\n');
            end
        end
    end
    
    save(filename, vars{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure

% figure
% scatter(t_obs(spline_robust_optimal.outlierIndices),x_obs(spline_robust_optimal.outlierIndices),(8.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'), hold on
% scatter(t_obs(trueOutlierIndices),x_obs(trueOutlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
% scatter(t_obs,x_obs,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
% tq = linspace(min(t_obs),max(t_obs),10*length(t_obs));
% plot(tq,spline_optimal(tq),'k')
% plot(tq,spline_robust_optimal(tq),'b')
% % plot(tq,spline_robust_cv(tq),'m')
% return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analyze the mse

dmse = @(mse) mse./optimal.mse-1;
pct_range = 0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*sum(~isnan(values))));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*sum(~isnan(values))));

optimal.dmse = dmse(optimal.mse);
robust_beta2_optimal.dmse = dmse(robust_beta2_optimal.mse);
robust_beta50_optimal.dmse = dmse(robust_beta50_optimal.mse);
robust_beta100_optimal.dmse = dmse(robust_beta100_optimal.mse);
robust_beta200_optimal.dmse = dmse(robust_beta200_optimal.mse);
robust_beta400_optimal.dmse = dmse(robust_beta400_optimal.mse);
robust_beta800_optimal.dmse = dmse(robust_beta800_optimal.mse);
robust_noiseOdds1over10.dmse = dmse(robust_noiseOdds1over10.mse);
robust_noiseOdds1over3.dmse = dmse(robust_noiseOdds1over3.mse);
robust_noiseOdds1.dmse = dmse(robust_noiseOdds1.mse);
robust_noiseOdds2.dmse = dmse(robust_noiseOdds2.mse);
robust_noiseOdds3.dmse = dmse(robust_noiseOdds3.mse);
robust_noiseOdds4.dmse = dmse(robust_noiseOdds4.mse);
robust_noiseOdds10.dmse = dmse(robust_noiseOdds10.mse);
robust_varianceCrossover.dmse = dmse(robust_varianceCrossover.mse);
robust_weighted.dmse = dmse(robust_weighted.mse);

print_pct = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f (%.1f-%.1f) ',100*median(stats.dmse(iOutlierRatio,iStride,iSlope,:),'omitnan'),minpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))), maxpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))) );
printcol = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %#.3g/%#.3g m$^2$ (%#.3g) (%d/%d) ', median(stats.mse(iOutlierRatio,iStride,iSlope,:)),mean(stats.mse(iOutlierRatio,iStride,iSlope,:)), median(stats.neff_se(iOutlierRatio,iStride,iSlope,:)), median(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)) );
print_dmse_all = @(stats) fprintf('& %.2f (%.2f) (%.2f-%.2f)\n',100*mean(stats.dmse(:),'omitnan'),100*median(stats.dmse(:),'omitnan'),minpct(sort(stats.dmse(:))), maxpct(sort(stats.dmse(:))) );
print_dmse_outlier = @(stats,iOutlierRatio) fprintf('& %.2f (%.2f) (%.2f-%.2f)\n',mean(reshape(stats.dmse(iOutlierRatio,:,:,:),[],1),'omitnan'),median(reshape(stats.dmse(iOutlierRatio,:,:,:),[],1),'omitnan'),minpct(sort(stats.dlambda(:)))/100, maxpct(sort(stats.dlambda(:)))/100 );

print_dmse_all(robust_beta2_optimal);
print_dmse_all(robust_beta50_optimal);
print_dmse_all(robust_beta100_optimal);
print_dmse_all(robust_beta200_optimal);
print_dmse_all(robust_beta400_optimal);
print_dmse_all(robust_beta800_optimal);
print_dmse_all(robust_noiseOdds1over10);
print_dmse_all(robust_noiseOdds1over3);
print_dmse_all(robust_noiseOdds1);
print_dmse_all(robust_noiseOdds2);
print_dmse_all(robust_noiseOdds3);
print_dmse_all(robust_noiseOdds4);
print_dmse_all(robust_noiseOdds10);
print_dmse_all(robust_varianceCrossover);
print_dmse_all(robust_weighted);

% Analyze the lambda---negative values indicate undertensioning, positive
% values indicate over tensioning
dlambda = @(stats) log10(stats.lambda ./ optimal.lambda);
robust_beta2_optimal.dlambda = dlambda(robust_beta2_optimal);
robust_beta50_optimal.dlambda = dlambda(robust_beta50_optimal);
robust_beta100_optimal.dlambda = dlambda(robust_beta100_optimal);
robust_beta200_optimal.dlambda = dlambda(robust_beta200_optimal);
robust_beta400_optimal.dlambda = dlambda(robust_beta400_optimal);
robust_beta800_optimal.dlambda = dlambda(robust_beta800_optimal);
robust_noiseOdds1over10.dlambda = dlambda(robust_noiseOdds1over10);
robust_noiseOdds1over3.dlambda = dlambda(robust_noiseOdds1over3);
robust_noiseOdds1.dlambda = dlambda(robust_noiseOdds1);
robust_noiseOdds2.dlambda = dlambda(robust_noiseOdds2);
robust_noiseOdds3.dlambda = dlambda(robust_noiseOdds3);
robust_noiseOdds4.dlambda = dlambda(robust_noiseOdds4);
robust_noiseOdds10.dlambda = dlambda(robust_noiseOdds10);
robust_varianceCrossover.dlambda = dlambda(robust_varianceCrossover);
robust_weighted.dlambda = dlambda(robust_weighted);


print_dlambda = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.2f (%.2f) (%.2f-%.2f) ',mean(stats.dlambda(iOutlierRatio,iStride,iSlope,:),'omitnan'),median(stats.dlambda(iOutlierRatio,iStride,iSlope,:),'omitnan'),minpct(sort(stats.dlambda(iOutlierRatio,iStride,iSlope,:)))/100, maxpct(sort(stats.dlambda(iOutlierRatio,iStride,iSlope,:)))/100 );

print_dlambda_all = @(stats) fprintf('& %.2f (%.2f) (%.2f-%.2f)\n',mean(stats.dlambda(:),'omitnan'),median(stats.dlambda(:),'omitnan'),minpct(sort(stats.dlambda(:)))/100, maxpct(sort(stats.dlambda(:)))/100 );
print_dlambda_outlier = @(stats,iOutlierRatio) fprintf('& %.2f (%.2f) (%.2f-%.2f)\n',mean(reshape(stats.dlambda(iOutlierRatio,:,:,:),[],1),'omitnan'),median(reshape(stats.dlambda(iOutlierRatio,:,:,:),[],1),'omitnan'),minpct(sort(stats.dlambda(:)))/100, maxpct(sort(stats.dlambda(:)))/100 );

print_dlambda_all(robust_beta2_optimal);
print_dlambda_all(robust_beta50_optimal);
print_dlambda_all(robust_beta100_optimal);
print_dlambda_all(robust_beta200_optimal);
print_dlambda_all(robust_beta400_optimal);
print_dlambda_all(robust_beta800_optimal);
print_dlambda_all(robust_noiseOdds1over10);
print_dlambda_all(robust_noiseOdds1over3);
print_dlambda_all(robust_noiseOdds1);
print_dlambda_all(robust_noiseOdds2);
print_dlambda_all(robust_noiseOdds3);
print_dlambda_all(robust_noiseOdds4);
print_dlambda_all(robust_noiseOdds10);
print_dlambda_all(robust_varianceCrossover);
print_dlambda_all(robust_weighted);

% all_maxpct = @(a,b,c,d,e,f,iOutlierRatio,iStride,iSlope) [maxpct(sort(a.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(b.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(c.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(d.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(e.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(f.dmse(iOutlierRatio,iStride,iSlope,:)))];
% 
% ratioed_maxpct = @(a,b,c,d,e,f,iOutlierRatio) [maxpct(sort(reshape(a.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(b.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(c.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(d.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(e.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(f.dmse(iOutlierRatio,:,:,:),[],1)))];
% ratioed_maxpct(robust_beta50_optimal, robust_beta100_optimal, robust_beta200_optimal, robust_beta400_optimal, robust_beta800_optimal,robust_beta800_optimal,1)
% ratioed_maxpct(robust_beta50_optimal, robust_beta100_optimal, robust_beta200_optimal, robust_beta400_optimal, robust_beta800_optimal,robust_beta800_optimal,2)
% ratioed_maxpct(robust_beta50_optimal, robust_beta100_optimal, robust_beta200_optimal, robust_beta400_optimal, robust_beta800_optimal,robust_beta800_optimal,3)
% 
% 
% all_all_maxpct = @(a,b,c,d,e,f) [maxpct(sort(a.dmse(:))),maxpct(sort(b.dmse(:))),maxpct(sort(c.dmse(:))),maxpct(sort(d.dmse(:))),maxpct(sort(e.dmse(:))),maxpct(sort(f.dmse(:)))];
% all_all_maxpct(robust_beta50_optimal, robust_beta100_optimal, robust_beta200_optimal, robust_beta400_optimal, robust_beta800_optimal,robust_beta800_optimal)

print_something_outlier = @(stat,iOutlierRatio) print_dlambda_outlier(stat,iOutlierRatio);

fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & optimal mse ($n_{eff}$) & blind optimal & robust & false pos/neg & robust 2nd & false pos/neg \\\\ \\hline \\hline \n');
for iOutlierRatio = 1:totalOutlierRatios
    
    print_something_outlier(robust_beta50_optimal,iOutlierRatio);
    print_something_outlier(robust_beta200_optimal,iOutlierRatio);
    print_something_outlier(robust_beta800_optimal,iOutlierRatio);
    
    print_something_outlier(robust_noiseOdds1over10,iOutlierRatio);
    print_something_outlier(robust_noiseOdds1over3,iOutlierRatio);
    print_something_outlier(robust_noiseOdds1,iOutlierRatio);
    print_something_outlier(robust_noiseOdds2,iOutlierRatio);
    print_something_outlier(robust_noiseOdds3,iOutlierRatio);
    print_something_outlier(robust_noiseOdds4,iOutlierRatio);
    print_something_outlier(robust_noiseOdds10,iOutlierRatio);
    print_something_outlier(robust_varianceCrossover,iOutlierRatio);
    print_something_outlier(robust_weighted,iOutlierRatio);
    
    fprintf(' \\\\ \n');
    
    continue
    
    for iSlope = 1:length(slopes)
        fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
        for iStride=1:length(strides)
            fprintf('%d ', strides(iStride));
                        
%             print_pct(robust_beta50_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta100_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta200_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta400_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta800_optimal,iOutlierRatio,iStride,iSlope);
            print_pct(robust_noiseOdds1over10,iOutlierRatio,iStride,iSlope);
            print_pct(robust_noiseOdds1over3,iOutlierRatio,iStride,iSlope);
            print_pct(robust_noiseOdds1,iOutlierRatio,iStride,iSlope);
            print_pct(robust_noiseOdds2,iOutlierRatio,iStride,iSlope);
            print_pct(robust_noiseOdds3,iOutlierRatio,iStride,iSlope);
            print_pct(robust_noiseOdds4,iOutlierRatio,iStride,iSlope);
            print_pct(robust_noiseOdds10,iOutlierRatio,iStride,iSlope);
            print_pct(robust_weighted,iOutlierRatio,iStride,iSlope);
            
%             print_dlambda(robust_beta2_optimal,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_beta50_optimal,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_beta100_optimal,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_beta200_optimal,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_beta400_optimal,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_beta800_optimal,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_noiseOdds1over10,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_noiseOdds1over3,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_noiseOdds1,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_noiseOdds2,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_noiseOdds3,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_noiseOdds4,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_noiseOdds10,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_weighted,iOutlierRatio,iStride,iSlope);

            
            fprintf(' \\\\ \n');
            
        end
    end
    fprintf('\n\n');
end
fprintf('\\end{tabular} \n');

