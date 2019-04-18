%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script generates a signal using the Matern with three different
% slopes. It then adds noise to the signal that include a known
% (quantified) portion and an outlier portion.


scaleFactor = 1;
LoadFigureDefaults

shouldUseStudentTDistribution = 1;

if shouldUseStudentTDistribution == 1
    filename = 'MSETableOutlierDetectionUnblindStudentT.mat';
else
    filename = 'MSETableOutlierDetectionUnblindNormal.mat';
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
    totalEnsembles = 11; % best to choose an odd number for median
    
    outlierRatios = [0; 0.025; 0.15; 0.25];
    outlierRatios = 0.25;
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
    nothing = zeros(totalOutlierRatios, totalStrides, totalSlopes, totalEnsembles);
    nothing_struct = struct('mse',nothing,'neff_se',nothing,'false_negatives', nothing, 'false_positives', nothing);
    
    total_outliers = nothing; vars{end+1} = 'total_outliers';
    
    default_dist = nothing_struct; vars{end+1} = 'default_dist';
    optimal_dist = nothing_struct; vars{end+1} = 'optimal_dist';
    optimal_dist_optimal_weight = nothing_struct; vars{end+1} = 'optimal_dist_optimal_weight';
    optimal_dist_outlier_weight = nothing_struct; vars{end+1} = 'optimal_dist_outlier_weight';
    optimal_dist_crossover_weight = nothing_struct; vars{end+1} = 'optimal_dist_crossover_weight';
    optimal_dist_no_outliers = nothing_struct; vars{end+1} = 'optimal_dist_no_outliers';
    optimal_dist_no_true_outliers = nothing_struct; vars{end+1} = 'optimal_dist_no_true_outliers';
    optimal_dist_magic_cutoff = nothing_struct; vars{end+1} = 'optimal_dist_magic_cutoff';
        
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

                    %%%%%%%%%%%%%%%%%%
                    % Default distribution
                    %%%%%%%%%%%%%%%%%%
                    spline_robust = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    default_dist = LogStatisticsFromSplineForOutlierDetectionTable(default_dist,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);

                    %%%%%%%%%%%%%%%%%%
                    % Known distribution
                    %%%%%%%%%%%%%%%%%%
                    spline_robust = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected);
                    spline_robust.outlierDistribution = outlierDistribution;
                    spline_robust.distribution = AddedDistribution(percentOutliers,outlierDistribution,noiseDistribution);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    optimal_dist = LogStatisticsFromSplineForOutlierDetectionTable(optimal_dist,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);

                    %%%%%%%%%%%%%%%%%%
                    % Known distribution, with known weightings
                    %%%%%%%%%%%%%%%%%%
                    spline_robust.distribution.w = @(z) (~outlierIndices) .* spline_robust.noiseDistribution.w(z) + outlierIndices .* spline_robust.outlierDistribution.w(z);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    optimal_dist_optimal_weight = LogStatisticsFromSplineForOutlierDetectionTable(optimal_dist_optimal_weight,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    % Known distribution, with naive weightings
                    %%%%%%%%%%%%%%%%%%
                    naiveOutliers = (abs(epsilon) > outlierThreshold);
                    spline_robust.distribution.w = @(z) (~naiveOutliers) .* spline_robust.noiseDistribution.w(z) + naiveOutliers .* spline_robust.outlierDistribution.w(z);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    optimal_dist_outlier_weight = LogStatisticsFromSplineForOutlierDetectionTable(optimal_dist_outlier_weight,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    % Known distribution, with crossover weightings
                    %%%%%%%%%%%%%%%%%%
                    f = @(z) abs( (1-percentOutliers)*noiseDistribution.pdf(z) - percentOutliers*outlierDistribution.pdf(z) );
                    z_crossover = abs(fminsearch(f,sqrt(noiseDistribution.variance)));
                    naiveOutliers = (abs(epsilon) > z_crossover);
                    spline_robust.distribution.w = @(z) (~naiveOutliers) .* spline_robust.noiseDistribution.w(z) + naiveOutliers .* spline_robust.outlierDistribution.w(z);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    optimal_dist_crossover_weight = LogStatisticsFromSplineForOutlierDetectionTable(optimal_dist_crossover_weight,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    % Known distribution, discarded outliers
                    %%%%%%%%%%%%%%%%%%
                    spline_robust = RobustSmoothingSpline(t_obs(~outlierIndices),x_obs(~outlierIndices),noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    optimal_dist_no_outliers = LogStatisticsFromSplineForOutlierDetectionTable(optimal_dist_no_outliers,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    % Known distribution, discarded crossover points
                    %%%%%%%%%%%%%%%%%%
                    spline_robust = RobustSmoothingSpline(t_obs(truenonOutlierIndices),x_obs(truenonOutlierIndices),noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    optimal_dist_no_true_outliers = LogStatisticsFromSplineForOutlierDetectionTable(optimal_dist_no_true_outliers,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);

                    %%%%%%%%%%%%%%%%%%
                    % Known distribution, with estimated optimal cutoff and weightings
                    %%%%%%%%%%%%%%%%%%
                    spline_robust = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    addedDistribution = AddedDistribution(percentOutliers,outlierDistribution,noiseDistribution);
                    m=2/3;
                    f = @(alpha) ((addedDistribution.variance-2*addedDistribution.varianceInRange(-Inf,addedDistribution.locationOfCDFPercentile(alpha/2)))^(1-m/2))/(1-alpha)^(m-1) + ((2*addedDistribution.varianceInRange(-Inf,addedDistribution.locationOfCDFPercentile(alpha/2)))^(1-m/2))/(alpha)^(m-1);
                    alpha_crossover = fminsearch(f,0.01);
                    z_crossover = abs(addedDistribution.locationOfCDFPercentile(alpha_crossover/2))/2;
                    naiveOutliers = (abs(epsilon) > z_crossover);
                    spline_robust.distribution.w = @(z) (~naiveOutliers) .* spline_robust.noiseDistribution.w(z) + naiveOutliers .* spline_robust.outlierDistribution.w(z);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    optimal_dist_magic_cutoff = LogStatisticsFromSplineForOutlierDetectionTable(optimal_dist_magic_cutoff,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);

                end
                fprintf('\n');
            end
        end
    end
    
    save(filename, vars{:});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the change in means square error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dmse = @(mse) mse./default_dist.mse-1;
pct_range = 0.90;%0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

default_dist.dmse = dmse(default_dist.mse);
optimal_dist.dmse = dmse(optimal_dist.mse);
optimal_dist_optimal_weight.dmse = dmse(optimal_dist_optimal_weight.mse);
optimal_dist_outlier_weight.dmse = dmse(optimal_dist_outlier_weight.mse);
optimal_dist_crossover_weight.dmse = dmse(optimal_dist_crossover_weight.mse);
optimal_dist_no_outliers.dmse = dmse(optimal_dist_no_outliers.mse);
optimal_dist_no_true_outliers.dmse = dmse(optimal_dist_no_true_outliers.mse);

print_pct = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f (%.1f-%.1f) ',100*mean(stats.dmse(iOutlierRatio,iStride,iSlope,:)),minpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))), maxpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))) );

% printcol = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f/%.1f (%.1f/%.1f) ', mean(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), mean(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)) );
printcol = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %#.3g/%#.3g m$^2$ (%#.3g) (%d/%d) ', median(stats.mse(iOutlierRatio,iStride,iSlope,:)),mean(stats.mse(iOutlierRatio,iStride,iSlope,:)), median(stats.neff_se(iOutlierRatio,iStride,iSlope,:)), median(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)) );


fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & optimal mse ($n_{eff}$) & blind optimal & robust & false pos/neg & robust 2nd & false pos/neg \\\\ \\hline \\hline \n');
for iOutlierRatio = 1:totalOutlierRatios
    for iSlope = 1:length(slopes)
        fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
        for iStride=1:length(strides)
            fprintf('%d ', strides(iStride));
 
            print_pct(default_dist,iOutlierRatio,iStride,iSlope);
            print_pct(optimal_dist,iOutlierRatio,iStride,iSlope);
            print_pct(optimal_dist_optimal_weight,iOutlierRatio,iStride,iSlope);
            print_pct(optimal_dist_outlier_weight,iOutlierRatio,iStride,iSlope);
            print_pct(optimal_dist_crossover_weight,iOutlierRatio,iStride,iSlope);
            print_pct(optimal_dist_no_outliers,iOutlierRatio,iStride,iSlope);
            print_pct(optimal_dist_no_true_outliers,iOutlierRatio,iStride,iSlope);

            fprintf(' \\\\ \n');
            
        end
    end
    fprintf('\n\n');
end
fprintf('\\end{tabular} \n');

