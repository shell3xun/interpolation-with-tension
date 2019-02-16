%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script generates a signal using the Matern with three different
% slopes. It then adds noise to the signal that include a known
% (quantified) portion and an outlier portion.


scaleFactor = 1;
LoadFigureDefaults

shouldUseStudentTDistribution = 1;

if shouldUseStudentTDistribution == 1
    filename = 'MSETableOutlierDetectionStudentT.mat';
else
    filename = 'MSETableOutlierDetectionNormal.mat';
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
    
%     full_tension_sv_p6 = nothing_struct; vars{end+1} = 'full_tension_sv_p6';
%     full_tension_sv_p5 = nothing_struct; vars{end+1} = 'full_tension_sv_p5';
%     full_tension_sv_p4 = nothing_struct; vars{end+1} = 'full_tension_sv_p4';
%     full_tension_sv_p25 = nothing_struct; vars{end+1} = 'full_tension_sv_p25';
%     full_tension_sv_p125 = nothing_struct; vars{end+1} = 'full_tension_sv_p125';
    
    optimal = nothing_struct; vars{end+1} = 'optimal';
    
    full_tension_innerSV_p6 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p6';
    full_tension_innerSV_p5 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p5';
    full_tension_innerSV_p4 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p4';
    full_tension_innerSV_p25 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p25';
    full_tension_innerSV_p125 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p125';
    
    full_tension_innerSV_p6x2 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p6x2';
    full_tension_innerSV_p5x2 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p5x2';
    full_tension_innerSV_p4x2 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p4x2';
    full_tension_innerSV_p25x2 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p25x2';
    full_tension_innerSV_p125x2 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p125x2';
        
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
                    trueGoodIndices = setdiff(1:n,trueOutlierIndices);
                    
                    total_outliers(iOutlierRatio,iStride,iSlope,iEnsemble) = length(trueOutlierIndices);
                    
                    f = @(z) abs( (1-percentOutliers)*noiseDistribution.pdf(z) - percentOutliers*outlierDistribution.pdf(z) );
                    z_crossover = abs(fminsearch(f,sqrt(noiseDistribution.variance)));
                    crossoverOutliers = (abs(epsilon) > z_crossover);
                    
                    x_obs = data.x + epsilon;
                    t_obs = data.t;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Unblinded best fit with standard tension spline
                    
                    linearIndex = sub2ind(size(nothing),iOutlierRatio,iStride,iSlope,iEnsemble);

                    % Optimal MSE solution, given the default settings
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    optimal = LogStatisticsFromSplineForOutlierDetectionTable(optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    alpha = 3/5; % .6
                    %%%%%%%%%%%%%%%%%%
                    
                    % 1st strategy
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,noiseDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings();
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    full_tension_innerSV_p6 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p6,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    % 2nd strategy
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,noiseDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    [outlierDistribution, alpha_outliers] = spline_robust.estimateOutlierDistribution();
                    addedDistribution = AddedDistribution(alpha_outliers,outlierDistribution,noiseDistribution);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,addedDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings();
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    full_tension_innerSV_p6x2 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p6x2,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    alpha = 1/2; % .5
                    %%%%%%%%%%%%%%%%%%
                    
                    % 1st strategy
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,noiseDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings();
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    full_tension_innerSV_p5 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p5,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    % 2nd strategy
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,noiseDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    [outlierDistribution, alpha_outliers] = spline_robust.estimateOutlierDistribution();
                    addedDistribution = AddedDistribution(alpha_outliers,outlierDistribution,noiseDistribution);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,addedDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings();
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    full_tension_innerSV_p5x2 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p5x2,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);

                    
                    %%%%%%%%%%%%%%%%%%
                    alpha = 2/5; % .4
                    %%%%%%%%%%%%%%%%%%
                    
                    % 1st strategy
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,noiseDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings();
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    full_tension_innerSV_p4 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p4,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    % 2nd strategy
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,noiseDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    [outlierDistribution, alpha_outliers] = spline_robust.estimateOutlierDistribution();
                    addedDistribution = AddedDistribution(alpha_outliers,outlierDistribution,noiseDistribution);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,addedDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings();
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    full_tension_innerSV_p4x2 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p4x2,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);

                    
                    %%%%%%%%%%%%%%%%%%
                    alpha = 1/4; % .25
                    %%%%%%%%%%%%%%%%%%
                    
                    % 1st strategy
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,noiseDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings();
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    full_tension_innerSV_p25 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p25,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    % 2nd strategy
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,noiseDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    [outlierDistribution, alpha_outliers] = spline_robust.estimateOutlierDistribution();
                    addedDistribution = AddedDistribution(alpha_outliers,outlierDistribution,noiseDistribution);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,addedDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings();
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    full_tension_innerSV_p25x2 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p25x2,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);

                    
                    %%%%%%%%%%%%%%%%%%
                    alpha = 1/8; % .125
                    %%%%%%%%%%%%%%%%%%
                    
                    % 1st strategy
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,noiseDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings();
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    full_tension_innerSV_p125 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p125,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    % 2nd strategy
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,noiseDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    [outlierDistribution, alpha_outliers] = spline_robust.estimateOutlierDistribution();
                    addedDistribution = AddedDistribution(alpha_outliers,outlierDistribution,noiseDistribution);
                    spline_robust.setInnerVarianceToExpectedValue(alpha,addedDistribution.varianceInPercentileRange(alpha/2,1-alpha/2));
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings();
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    full_tension_innerSV_p125x2 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p125x2,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);

                end
                fprintf('\n');
            end
        end
    end
    
%     save(filename, vars{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure

% figure
% scatter(t_obs(spline_robust_optimal.indicesOfOutliers),x_obs(spline_robust_optimal.indicesOfOutliers),(8.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'), hold on
% scatter(t_obs(trueOutlierIndices),x_obs(trueOutlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
% scatter(t_obs,x_obs,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
% tq = linspace(min(t_obs),max(t_obs),10*length(t_obs));
% plot(tq,spline_optimal(tq),'k')
% plot(tq,spline_robust_optimal(tq),'b')
% % plot(tq,spline_robust_cv(tq),'m')
% return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the change in means square error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dmse = @(mse) mse./optimal.mse-1;
pct_range = 0.90;%0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

optimal.dmse = dmse(optimal.mse);
full_tension_innerSV_p6.dmse = dmse(full_tension_innerSV_p6.mse);
full_tension_innerSV_p5.dmse = dmse(full_tension_innerSV_p5.mse);
full_tension_innerSV_p4.dmse = dmse(full_tension_innerSV_p4.mse);
full_tension_innerSV_p25.dmse = dmse(full_tension_innerSV_p25.mse);
full_tension_innerSV_p125.dmse = dmse(full_tension_innerSV_p125.mse);

full_tension_innerSV_p6x2.dmse = dmse(full_tension_innerSV_p6x2.mse);
full_tension_innerSV_p5x2.dmse = dmse(full_tension_innerSV_p5x2.mse);
full_tension_innerSV_p4x2.dmse = dmse(full_tension_innerSV_p4x2.mse);
full_tension_innerSV_p25x2.dmse = dmse(full_tension_innerSV_p25x2.mse);
full_tension_innerSV_p125x2.dmse = dmse(full_tension_innerSV_p125x2.mse);

print_pct = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f (%.1f-%.1f) ',100*mean(stats.dmse(iOutlierRatio,iStride,iSlope,:)),minpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))), maxpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute false positive/negative rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printall = @(stats) fprintf('%.1f/%.1f (%.1f/%.1f)\n', mean(stats.false_positives(:)), mean(stats.false_negatives(:)), median(stats.false_positives(:)), median(stats.false_negatives(:)));

printall(full_tension_innerSV_p6);
printall(full_tension_innerSV_p5);
printall(full_tension_innerSV_p4);
printall(full_tension_innerSV_p25);
printall(full_tension_innerSV_p125);

printall(full_tension_innerSV_p6x2);
printall(full_tension_innerSV_p5x2);
printall(full_tension_innerSV_p4x2);
printall(full_tension_innerSV_p25x2);
printall(full_tension_innerSV_p125x2);

% printcol = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f/%.1f (%.1f/%.1f) ', mean(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), mean(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)) );
printcol = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %#.3g/%#.3g m$^2$ (%#.3g) (%d/%d) ', median(stats.mse(iOutlierRatio,iStride,iSlope,:)),mean(stats.mse(iOutlierRatio,iStride,iSlope,:)), median(stats.neff_se(iOutlierRatio,iStride,iSlope,:)), median(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)) );


fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & optimal mse ($n_{eff}$) & blind optimal & robust & false pos/neg & robust 2nd & false pos/neg \\\\ \\hline \\hline \n');
for iOutlierRatio = 1:totalOutlierRatios
    for iSlope = 1:length(slopes)
        fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
        for iStride=1:length(strides)
            fprintf('%d ', strides(iStride));
            
%             printcol(full_tension_innerSV_p6,iOutlierRatio,iStride,iSlope);
%             printcol(full_tension_innerSV_p5,iOutlierRatio,iStride,iSlope);
%             printcol(full_tension_innerSV_p4,iOutlierRatio,iStride,iSlope);
%             printcol(full_tension_innerSV_p25,iOutlierRatio,iStride,iSlope);
%             printcol(full_tension_innerSV_p125,iOutlierRatio,iStride,iSlope);
%             
%             printcol(full_tension_innerSV_p6x2,iOutlierRatio,iStride,iSlope);
%             printcol(full_tension_innerSV_p5x2,iOutlierRatio,iStride,iSlope);
%             printcol(full_tension_innerSV_p4x2,iOutlierRatio,iStride,iSlope);
%             printcol(full_tension_innerSV_p25x2,iOutlierRatio,iStride,iSlope);
%             printcol(full_tension_innerSV_p125x2,iOutlierRatio,iStride,iSlope);
            
            print_pct(full_tension_innerSV_p6,iOutlierRatio,iStride,iSlope);
            print_pct(full_tension_innerSV_p5,iOutlierRatio,iStride,iSlope);
            print_pct(full_tension_innerSV_p4,iOutlierRatio,iStride,iSlope);
            print_pct(full_tension_innerSV_p25,iOutlierRatio,iStride,iSlope);
            print_pct(full_tension_innerSV_p125,iOutlierRatio,iStride,iSlope);
            
            print_pct(full_tension_innerSV_p6x2,iOutlierRatio,iStride,iSlope);
            print_pct(full_tension_innerSV_p5x2,iOutlierRatio,iStride,iSlope);
            print_pct(full_tension_innerSV_p4x2,iOutlierRatio,iStride,iSlope);
            print_pct(full_tension_innerSV_p25x2,iOutlierRatio,iStride,iSlope);
            print_pct(full_tension_innerSV_p125x2,iOutlierRatio,iStride,iSlope);
            
            fprintf(' \\\\ \n');
            
        end
    end
    fprintf('\n\n');
end
fprintf('\\end{tabular} \n');
