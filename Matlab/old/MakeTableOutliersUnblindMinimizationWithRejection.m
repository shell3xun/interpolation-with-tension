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
    filename = 'MSETableOutliersUnblindMinimizationWithRejectionStudentT.mat';
else
    filename = 'MSETableOutliersUnblindMinimizationRejectionNormal.mat';
end

if exist(filename,'file')
    load(filename);
    totalOutlierRatios = size(total_outliers,1);
else
    slopes = [-2; -3; -4];
    slopes = -3;
    totalSlopes = length(slopes);

    strides = [5;20;80;200];
    strides = 200;
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
    nothing_struct = struct('mse',nothing,'neff_se',nothing,'lambda',nothing,'nonOutlierEffectiveSampleSize',nothing,'nonOutlierSampleVariance',nothing,'false_negatives', nothing, 'false_positives', nothing,'rejects',nothing);
    
    total_outliers = nothing; vars{end+1} = 'total_outliers';
    
    optimal = nothing_struct; vars{end+1} = 'optimal';
    robust_noiseOdds1 = nothing_struct; vars{end+1} = 'robust_noiseOdds1';
    robust_outlierOdds1e4 = nothing_struct; vars{end+1} = 'robust_outlierOdds1e4';
    robust_outlierOdds3e4 = nothing_struct; vars{end+1} = 'robust_outlierOdds3e4';
    robust_outlierOdds1e5 = nothing_struct; vars{end+1} = 'robust_outlierOdds1e5';
    robust_outlierOdds3e5 = nothing_struct; vars{end+1} = 'robust_outlierOdds3e5';
    robust_outlierOdds1e6 = nothing_struct; vars{end+1} = 'robust_outlierOdds1e6';
    robust_outlierOdds3e6 = nothing_struct; vars{end+1} = 'robust_outlierOdds3e6';
    robust_outlierOdds1e7 = nothing_struct; vars{end+1} = 'robust_outlierOdds1e7';
    robust_outlierDynamic = nothing_struct; vars{end+1} = 'robust_outlierDynamic';
                    
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
                    % variance of the noise within different cutoffs based
                    % on the empirical fit
                    
                    spline_robust.setToFullTensionWithIteratedIQAD();
                    [empiricalOutlierDistribution,empiricalAlpha] = spline_robust.estimateOutlierDistribution();
                    epsilon_full = spline_robust.epsilon;
                    
                    if empiricalAlpha == 0 || isempty(empiricalOutlierDistribution)
                        continue;
                    end
                    
                    % First we perform our best blind minimization
                    noiseOdds = 1;
                    f = @(z) abs( (1-empiricalAlpha)*noiseDistribution.pdf(z) - noiseOdds*empiricalAlpha*empiricalOutlierDistribution.pdf(z) );
                    z_minimization_cutoff = abs(fminsearch(f,sqrt(noiseDistribution.variance)));
                    expectedVariance = spline_robust.noiseDistribution.varianceInRange(-z_minimization_cutoff,z_minimization_cutoff);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-z_minimization_cutoff,z_minimization_cutoff,expectedVariance) );
                    robust_noiseOdds1 = LogStatisticsFromSplineForOutlierTable(robust_noiseOdds1,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    if 1
                        minimize = @(spline) spline.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-z_minimization_cutoff,z_minimization_cutoff,expectedVariance) );
                    else
                        
                    end
                    
                    spline_robust.outlierDistribution = empiricalOutlierDistribution;
%                     spline_robust.alpha = empiricalAlpha;
%                     spline_robust.distribution = AddedDistribution(empiricalAlpha,empiricalOutlierDistribution,noiseDistribution);
                                        
                    % Now, let's run through several different rejection
                    % possibilities
                    
                    outlierOdds = 1e4;
                    f = @(z) abs( (spline_robust.alpha/(1-spline_robust.alpha))*spline_robust.outlierDistribution.pdf(-abs(z))/spline_robust.noiseDistribution.pdf(-abs(z)) - outlierOdds);
                    z_outlier = fminsearch(f,sqrt(spline_robust.noiseDistribution.variance));
                    noiseIndices = epsilon_full >= -z_outlier & epsilon_full <= z_outlier;
%                     spline_robust.distribution.w = @(z) noiseIndices .* spline_robust.noiseDistribution.w(z) + (~noiseIndices) .* spline_robust.outlierDistribution.w(z);
%                     spline_robust.distribution = noiseDistribution;
                    spline_robust.sigma = noiseIndices .* sqrt(spline_robust.noiseDistribution.variance) + (~noiseIndices) .* sqrt(spline_robust.outlierDistribution.variance);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    robust_outlierOdds1e4 = LogStatisticsFromSplineForOutlierTable(robust_outlierOdds1e4,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices,sum(~noiseIndices));
                    
                    outlierOdds = 3e4;
                    f = @(z) abs( (spline_robust.alpha/(1-spline_robust.alpha))*spline_robust.outlierDistribution.pdf(-abs(z))/spline_robust.noiseDistribution.pdf(-abs(z)) - outlierOdds);
                    z_outlier = fminsearch(f,sqrt(spline_robust.noiseDistribution.variance));
                    noiseIndices = epsilon_full >= -z_outlier & epsilon_full <= z_outlier;
%                     spline_robust.distribution.w = @(z) noiseIndices .* spline_robust.noiseDistribution.w(z) + (~noiseIndices) .* spline_robust.outlierDistribution.w(z);
%                     spline_robust.distribution = noiseDistribution;
                    spline_robust.sigma = noiseIndices .* sqrt(spline_robust.noiseDistribution.variance) + (~noiseIndices) .* sqrt(spline_robust.outlierDistribution.variance);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    robust_outlierOdds3e4 = LogStatisticsFromSplineForOutlierTable(robust_outlierOdds3e4,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices,sum(~noiseIndices));
                    
                    outlierOdds = 1e5;
                    f = @(z) abs( (spline_robust.alpha/(1-spline_robust.alpha))*spline_robust.outlierDistribution.pdf(-abs(z))/spline_robust.noiseDistribution.pdf(-abs(z)) - outlierOdds);
                    z_outlier = fminsearch(f,sqrt(spline_robust.noiseDistribution.variance));
                    noiseIndices = epsilon_full >= -z_outlier & epsilon_full <= z_outlier;
%                     spline_robust.distribution.w = @(z) noiseIndices .* spline_robust.noiseDistribution.w(z) + (~noiseIndices) .* spline_robust.outlierDistribution.w(z);
%                     spline_robust.distribution = noiseDistribution;
                    spline_robust.sigma = noiseIndices .* sqrt(spline_robust.noiseDistribution.variance) + (~noiseIndices) .* sqrt(spline_robust.outlierDistribution.variance);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    robust_outlierOdds1e5 = LogStatisticsFromSplineForOutlierTable(robust_outlierOdds1e5,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices,sum(~noiseIndices));
                    
                    outlierOdds = 3e5;
                    f = @(z) abs( (spline_robust.alpha/(1-spline_robust.alpha))*spline_robust.outlierDistribution.pdf(-abs(z))/spline_robust.noiseDistribution.pdf(-abs(z)) - outlierOdds);
                    z_outlier = fminsearch(f,sqrt(spline_robust.noiseDistribution.variance));
                    noiseIndices = epsilon_full >= -z_outlier & epsilon_full <= z_outlier;
%                     spline_robust.distribution.w = @(z) noiseIndices .* spline_robust.noiseDistribution.w(z) + (~noiseIndices) .* spline_robust.outlierDistribution.w(z);
%                     spline_robust.distribution = noiseDistribution;
                    spline_robust.sigma = noiseIndices .* sqrt(spline_robust.noiseDistribution.variance) + (~noiseIndices) .* sqrt(spline_robust.outlierDistribution.variance);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    robust_outlierOdds3e5 = LogStatisticsFromSplineForOutlierTable(robust_outlierOdds3e5,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices,sum(~noiseIndices));   
                    
                    outlierOdds = 1e6;
                    f = @(z) abs( (spline_robust.alpha/(1-spline_robust.alpha))*spline_robust.outlierDistribution.pdf(-abs(z))/spline_robust.noiseDistribution.pdf(-abs(z)) - outlierOdds);
                    z_outlier = fminsearch(f,sqrt(spline_robust.noiseDistribution.variance));
                    noiseIndices = epsilon_full >= -z_outlier & epsilon_full <= z_outlier;
%                     spline_robust.distribution.w = @(z) noiseIndices .* spline_robust.noiseDistribution.w(z) + (~noiseIndices) .* spline_robust.outlierDistribution.w(z);
%                     spline_robust.distribution = noiseDistribution;
                    spline_robust.sigma = noiseIndices .* sqrt(spline_robust.noiseDistribution.variance) + (~noiseIndices) .* sqrt(spline_robust.outlierDistribution.variance);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    robust_outlierOdds1e6 = LogStatisticsFromSplineForOutlierTable(robust_outlierOdds1e6,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices,sum(~noiseIndices));
                    
                    outlierOdds = 3e6;
                    f = @(z) abs( (spline_robust.alpha/(1-spline_robust.alpha))*spline_robust.outlierDistribution.pdf(-abs(z))/spline_robust.noiseDistribution.pdf(-abs(z)) - outlierOdds);
                    z_outlier = fminsearch(f,sqrt(spline_robust.noiseDistribution.variance));
                    noiseIndices = epsilon_full >= -z_outlier & epsilon_full <= z_outlier;
%                     spline_robust.distribution.w = @(z) noiseIndices .* spline_robust.noiseDistribution.w(z) + (~noiseIndices) .* spline_robust.outlierDistribution.w(z);
%                     spline_robust.distribution = noiseDistribution;
                    spline_robust.sigma = noiseIndices .* sqrt(spline_robust.noiseDistribution.variance) + (~noiseIndices) .* sqrt(spline_robust.outlierDistribution.variance);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    robust_outlierOdds3e6 = LogStatisticsFromSplineForOutlierTable(robust_outlierOdds3e6,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices,sum(~noiseIndices));
                    
                    outlierOdds = 1e7;
                    f = @(z) abs( (spline_robust.alpha/(1-spline_robust.alpha))*spline_robust.outlierDistribution.pdf(-abs(z))/spline_robust.noiseDistribution.pdf(-abs(z)) - outlierOdds);
                    z_outlier = fminsearch(f,sqrt(spline_robust.noiseDistribution.variance));
                    noiseIndices = epsilon_full >= -z_outlier & epsilon_full <= z_outlier;
%                     spline_robust.distribution.w = @(z) noiseIndices .* spline_robust.noiseDistribution.w(z) + (~noiseIndices) .* spline_robust.outlierDistribution.w(z);
%                     spline_robust.distribution = noiseDistribution;
                    spline_robust.sigma = noiseIndices .* sqrt(spline_robust.noiseDistribution.variance) + (~noiseIndices) .* sqrt(spline_robust.outlierDistribution.variance);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    robust_outlierOdds1e7 = LogStatisticsFromSplineForOutlierTable(robust_outlierOdds1e7,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices,sum(~noiseIndices));

                    sigma = abs(epsilon_full);
                    sigma(sigma<sqrt(spline_robust.noiseDistribution.variance))=sqrt(spline_robust.noiseDistribution.variance);
                    sigma(sigma>sqrt(empiricalOutlierDistribution.variance))=sqrt(empiricalOutlierDistribution.variance);
                    spline_robust.sigma = sigma;
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    robust_outlierDynamic = LogStatisticsFromSplineForOutlierTable(robust_outlierDynamic,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices,sum(~noiseIndices));

                end
                fprintf('\n');
            end
        end
    end
    
    save(filename, vars{:});
end

% Analyze the mse

dmse = @(mse) mse./optimal.mse-1;
pct_range = 0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

optimal.dmse = dmse(optimal.mse);
robust_noiseOdds1.dmse = dmse(robust_noiseOdds1.mse);
robust_outlierOdds1e4.dmse = dmse(robust_outlierOdds1e4.mse);
robust_outlierOdds3e4.dmse = dmse(robust_outlierOdds3e4.mse);
robust_outlierOdds1e5.dmse = dmse(robust_outlierOdds1e5.mse);
robust_outlierOdds3e5.dmse = dmse(robust_outlierOdds3e5.mse);
robust_outlierOdds1e6.dmse = dmse(robust_outlierOdds1e6.mse);
robust_outlierOdds3e6.dmse = dmse(robust_outlierOdds3e6.mse);
robust_outlierOdds1e7.dmse = dmse(robust_outlierOdds1e7.mse);
robust_outlierDynamic.dmse = dmse(robust_outlierDynamic.mse);

print_pct = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f (%.1f-%.1f) ',100*mean(stats.dmse(iOutlierRatio,iStride,iSlope,:)),minpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))), maxpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))) );
printcol = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %#.3g/%#.3g m$^2$ (%#.3g) (%d/%d) ', median(stats.mse(iOutlierRatio,iStride,iSlope,:)),mean(stats.mse(iOutlierRatio,iStride,iSlope,:)), median(stats.neff_se(iOutlierRatio,iStride,iSlope,:)), median(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)) );
print_dmse_all = @(stats) fprintf('& %.2f (%.2f) (%.2f-%.2f)\n',100*mean(stats.dmse(:)),100*median(stats.dmse(:)),minpct(sort(stats.dmse(:))), maxpct(sort(stats.dmse(:))) );

print_dmse_all(robust_noiseOdds1);
print_dmse_all(robust_outlierOdds1e4);
print_dmse_all(robust_outlierOdds3e4);
print_dmse_all(robust_outlierOdds1e5);
print_dmse_all(robust_outlierOdds3e5);
print_dmse_all(robust_outlierOdds1e6);
print_dmse_all(robust_outlierOdds3e6);
print_dmse_all(robust_outlierOdds1e7);
print_dmse_all(robust_outlierDynamic);

% Analyze the lambda---negative values indicate undertensioning, positive
% values indicate over tensioning
dlambda = @(stats) log10(stats.lambda ./ optimal.lambda);
robust_noiseOdds1.dlambda = dlambda(robust_noiseOdds1);
robust_outlierOdds1e4.dlambda = dlambda(robust_outlierOdds1e4);
robust_outlierOdds3e4.dlambda = dlambda(robust_outlierOdds3e4);
robust_outlierOdds1e5.dlambda = dlambda(robust_outlierOdds1e5);
robust_outlierOdds3e5.dlambda = dlambda(robust_outlierOdds3e5);
robust_outlierOdds1e6.dlambda = dlambda(robust_outlierOdds1e6);
robust_outlierOdds3e6.dlambda = dlambda(robust_outlierOdds3e6);
robust_outlierOdds1e7.dlambda = dlambda(robust_outlierOdds1e7);

print_dlambda = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.2f (%.2f) (%.2f-%.2f) ',mean(stats.dlambda(iOutlierRatio,iStride,iSlope,:)),median(stats.dlambda(iOutlierRatio,iStride,iSlope,:)),minpct(sort(stats.dlambda(iOutlierRatio,iStride,iSlope,:)))/100, maxpct(sort(stats.dlambda(iOutlierRatio,iStride,iSlope,:)))/100 );

print_dlambda_all = @(stats) fprintf('& %.2f (%.2f) (%.2f-%.2f)\n',mean(stats.dlambda(:)),median(stats.dlambda(:)),minpct(sort(stats.dlambda(:)))/100, maxpct(sort(stats.dlambda(:)))/100 );

print_dlambda_all(robust_noiseOdds1);
print_dlambda_all(robust_outlierOdds1e4);
print_dlambda_all(robust_outlierOdds3e4);
print_dlambda_all(robust_outlierOdds1e5);
print_dlambda_all(robust_outlierOdds3e5);
print_dlambda_all(robust_outlierOdds1e6);
print_dlambda_all(robust_outlierOdds3e6);
print_dlambda_all(robust_outlierOdds1e7);

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

fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & optimal mse ($n_{eff}$) & blind optimal & robust & false pos/neg & robust 2nd & false pos/neg \\\\ \\hline \\hline \n');
for iOutlierRatio = 1:totalOutlierRatios
    for iSlope = 1:length(slopes)
        fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
        for iStride=1:length(strides)
            fprintf('%d ', strides(iStride));
                        
%             print_pct(robust_beta50_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta100_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta200_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta400_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta800_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_noiseOdds1over10,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_noiseOdds1over3,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_noiseOdds1,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_noiseOdds2,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_noiseOdds3,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_noiseOdds4,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_noiseOdds10,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_weighted,iOutlierRatio,iStride,iSlope);
            
            print_dlambda(robust_noiseOdds1,iOutlierRatio,iStride,iSlope);
            print_dlambda(robust_outlierOdds1e4,iOutlierRatio,iStride,iSlope);
            print_dlambda(robust_outlierOdds3e4,iOutlierRatio,iStride,iSlope);
            print_dlambda(robust_outlierOdds1e5,iOutlierRatio,iStride,iSlope);
            print_dlambda(robust_outlierOdds3e5,iOutlierRatio,iStride,iSlope);
            print_dlambda(robust_outlierOdds1e6,iOutlierRatio,iStride,iSlope);
            print_dlambda(robust_outlierOdds3e6,iOutlierRatio,iStride,iSlope);
            print_dlambda(robust_outlierOdds1e7,iOutlierRatio,iStride,iSlope);

            
            fprintf(' \\\\ \n');
            
        end
    end
    fprintf('\n\n');
end
fprintf('\\end{tabular} \n');

