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
%     strides = 200;
    totalStrides = length(strides);
    totalEnsembles = 51; % best to choose an odd number for median
    
    outlierRatios = [0.05; 0.15; 0.25];
%     outlierRatios = 0.0;
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
      
    optimal = nothing_struct; vars{end+1} = 'optimal';
    outlier_cutoff_odds100 = nothing_struct; vars{end+1} = 'outlier_cutoff_odds100';
    outlier_cutoff_odds300 = nothing_struct; vars{end+1} = 'outlier_cutoff_odds300';
    outlier_cutoff_odds1000 = nothing_struct; vars{end+1} = 'outlier_cutoff_odds1000';
    outlier_cutoff_odds3000 = nothing_struct; vars{end+1} = 'outlier_cutoff_odds3000';
    outlier_cutoff_odds10000 = nothing_struct; vars{end+1} = 'outlier_cutoff_odds10000';    
    outlier_cutoff_odds30000 = nothing_struct; vars{end+1} = 'outlier_cutoff_odds30000';
    outlier_cutoff_odds100000 = nothing_struct; vars{end+1} = 'outlier_cutoff_odds100000';
    outlier_cutoff_odds300000 = nothing_struct; vars{end+1} = 'outlier_cutoff_odds300000';
    outlier_cutoff_odds1000000 = nothing_struct; vars{end+1} = 'outlier_cutoff_odds1000000';
    outlier_cutoff_odds3000000 = nothing_struct; vars{end+1} = 'outlier_cutoff_odds3000000';
        
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
                    
                    beta = 1/400;
                    zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                    zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                    expectedVariance = noiseDistribution.varianceInRange(zmin,zmax);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Unblinded best fit with standard tension spline
                    
                    linearIndex = sub2ind(size(nothing),iOutlierRatio,iStride,iSlope,iEnsemble);

                    % Optimal MSE solution, given the default settings
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.minimizeExpectedMeanSquareError();
%                     spline_robust.minimize(@(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax,expectedVariance));
                    optimal = LogStatisticsFromSplineForOutlierDetectionTable(optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    odds = 1e2;
                    %%%%%%%%%%%%%%%%%%

                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.setToFullTensionWithIteratedIQAD();
                    lambda_full = spline_robust.lambda;
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings(odds);
                    spline_robust.minimizeExpectedMeanSquareError();
                    outlier_cutoff_odds100 = LogStatisticsFromSplineForOutlierDetectionTable(outlier_cutoff_odds100,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    odds = 3e2;
                    %%%%%%%%%%%%%%%%%%

                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.lambda = lambda_full;
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings(odds);
                    spline_robust.minimizeExpectedMeanSquareError();
                    outlier_cutoff_odds300 = LogStatisticsFromSplineForOutlierDetectionTable(outlier_cutoff_odds300,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    odds = 1e3;
                    %%%%%%%%%%%%%%%%%%

                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.lambda = lambda_full;
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings(odds);
                    spline_robust.minimizeExpectedMeanSquareError();
                    outlier_cutoff_odds1000 = LogStatisticsFromSplineForOutlierDetectionTable(outlier_cutoff_odds1000,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    odds = 3e3;
                    %%%%%%%%%%%%%%%%%%

                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.lambda = lambda_full;
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings(odds);
                    spline_robust.minimizeExpectedMeanSquareError();
                    outlier_cutoff_odds3000 = LogStatisticsFromSplineForOutlierDetectionTable(outlier_cutoff_odds3000,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    odds = 1e4;
                    %%%%%%%%%%%%%%%%%%

                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.lambda = lambda_full;
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings(odds);
                    spline_robust.minimizeExpectedMeanSquareError();
                    outlier_cutoff_odds10000 = LogStatisticsFromSplineForOutlierDetectionTable(outlier_cutoff_odds10000,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    odds = 3e4;
                    %%%%%%%%%%%%%%%%%%

                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.lambda = lambda_full;
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings(odds);
                    spline_robust.minimizeExpectedMeanSquareError();
                    outlier_cutoff_odds30000 = LogStatisticsFromSplineForOutlierDetectionTable(outlier_cutoff_odds30000,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    odds = 1e5;
                    %%%%%%%%%%%%%%%%%%

                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.lambda = lambda_full;
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings(odds);
                    spline_robust.minimizeExpectedMeanSquareError();
                    outlier_cutoff_odds100000 = LogStatisticsFromSplineForOutlierDetectionTable(outlier_cutoff_odds100000,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    odds = 3e5;
                    %%%%%%%%%%%%%%%%%%
                    
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.lambda = lambda_full;
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings(odds);
                    spline_robust.minimizeExpectedMeanSquareError();
                    outlier_cutoff_odds300000 = LogStatisticsFromSplineForOutlierDetectionTable(outlier_cutoff_odds300000,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    odds = 1e6;
                    %%%%%%%%%%%%%%%%%%
                    
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.lambda = lambda_full;
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings(odds);
                    spline_robust.minimizeExpectedMeanSquareError();
                    outlier_cutoff_odds1000000 = LogStatisticsFromSplineForOutlierDetectionTable(outlier_cutoff_odds1000000,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%
                    odds = 3e6;
                    %%%%%%%%%%%%%%%%%%
                    
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    spline_robust.lambda = lambda_full;
                    spline_robust.rebuildOutlierDistributionAndAdjustWeightings(odds);
                    spline_robust.minimizeExpectedMeanSquareError();
                    outlier_cutoff_odds3000000 = LogStatisticsFromSplineForOutlierDetectionTable(outlier_cutoff_odds3000000,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
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
dmse = @(mse) log10(mse./optimal.mse);
pct_range = 0.90;%0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

optimal.dmse = dmse(optimal.mse);
outlier_cutoff_odds100.dmse = dmse(outlier_cutoff_odds100.mse);
outlier_cutoff_odds300.dmse = dmse(outlier_cutoff_odds300.mse);
outlier_cutoff_odds1000.dmse = dmse(outlier_cutoff_odds1000.mse);
outlier_cutoff_odds3000.dmse = dmse(outlier_cutoff_odds3000.mse);
outlier_cutoff_odds10000.dmse = dmse(outlier_cutoff_odds10000.mse);
outlier_cutoff_odds30000.dmse = dmse(outlier_cutoff_odds30000.mse);
outlier_cutoff_odds100000.dmse = dmse(outlier_cutoff_odds100000.mse);
outlier_cutoff_odds300000.dmse = dmse(outlier_cutoff_odds300000.mse);
outlier_cutoff_odds1000000.dmse = dmse(outlier_cutoff_odds1000000.mse);
outlier_cutoff_odds3000000.dmse = dmse(outlier_cutoff_odds3000000.mse);

print_pct = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f (%.1f-%.1f) ',100*median(stats.dmse(iOutlierRatio,iStride,iSlope,:)),minpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))), maxpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute false positive/negative rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% printall = @(stats) fprintf('%.1f/%.1f (%.1f/%.1f)\n', mean(stats.false_positives(:)), mean(stats.false_negatives(:)), median(stats.false_positives(:)), median(stats.false_negatives(:)));
% 
% printall(outlier_cutoff_odds100);
% printall(outlier_cutoff_odds300);
% printall(outlier_cutoff_odds1000);
% printall(outlier_cutoff_odds3000);
% printall(outlier_cutoff_odds10000);
% 
% printall(outlier_cutoff_odds30000);
% printall(outlier_cutoff_odds100000);
% printall(full_tension_innerSV_p4x2);
% printall(full_tension_innerSV_p25x2);
% printall(full_tension_innerSV_p125x2);
% 
% % printcol = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f/%.1f (%.1f/%.1f) ', mean(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), mean(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)) );
% printcol = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %#.3g/%#.3g m$^2$ (%#.3g) (%d/%d) ', median(stats.mse(iOutlierRatio,iStride,iSlope,:)),mean(stats.mse(iOutlierRatio,iStride,iSlope,:)), median(stats.neff_se(iOutlierRatio,iStride,iSlope,:)), median(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)) );


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
            
            print_pct(outlier_cutoff_odds100,iOutlierRatio,iStride,iSlope);
            print_pct(outlier_cutoff_odds300,iOutlierRatio,iStride,iSlope);
            print_pct(outlier_cutoff_odds1000,iOutlierRatio,iStride,iSlope);
            print_pct(outlier_cutoff_odds3000,iOutlierRatio,iStride,iSlope);
            print_pct(outlier_cutoff_odds10000,iOutlierRatio,iStride,iSlope);           
            print_pct(outlier_cutoff_odds30000,iOutlierRatio,iStride,iSlope);
            print_pct(outlier_cutoff_odds100000,iOutlierRatio,iStride,iSlope);
            print_pct(outlier_cutoff_odds300000,iOutlierRatio,iStride,iSlope);
            print_pct(outlier_cutoff_odds1000000,iOutlierRatio,iStride,iSlope);
           
            fprintf(' \\\\ \n');
            
        end
    end
    fprintf('\n\n');
end
fprintf('\\end{tabular} \n');

