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
    
    outlierRatios = [0; 0.025; 0.15; 0.25];
%     outlierRatios = 0.25;
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
   
    outlierFactor = 50;
    
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
    
    full_tension_innerSV_p6 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p6';
    full_tension_innerSV_p5 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p5';
    full_tension_innerSV_p4 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p4';
    full_tension_innerSV_p25 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p25';
    full_tension_innerSV_p125 = nothing_struct; vars{end+1} = 'full_tension_innerSV_p125';
        
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
                    
                    x_obs = data.x + epsilon;
                    t_obs = data.t;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Unblinded best fit with standard tension spline
                    
                    linearIndex = sub2ind(size(nothing),iOutlierRatio,iStride,iSlope,iEnsemble);
                    
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    lambda0 = spline_robust.lambda;
                    
                    alpha = 3/5; % .6
                    spline_robust.setToFullTensionWithInnerSV(alpha);
                    full_tension_innerSV_p6 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p6,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
%                     spline_robust.lambda = lambda0;
%                     spline_robust.setToFullTensionWithSV(alpha);
%                     full_tension_sv_p6 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_sv_p6,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    alpha = 1/2; % .5
                    spline_robust.lambda = lambda0;
                    spline_robust.setToFullTensionWithInnerSV(alpha);
                    full_tension_innerSV_p5 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p5,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
%                     spline_robust.lambda = lambda0;
%                     spline_robust.setToFullTensionWithSV(alpha);
%                     full_tension_sv_p5 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_sv_p5,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    alpha = 2/5; % .4
                    spline_robust.lambda = lambda0;
                    spline_robust.setToFullTensionWithInnerSV(alpha);
                    full_tension_innerSV_p4 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p4,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
%                     spline_robust.lambda = lambda0;
%                     spline_robust.setToFullTensionWithSV(alpha);
%                     full_tension_sv_p4 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_sv_p4,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    alpha = 1/4; % .25
                    spline_robust.lambda = lambda0;
                    spline_robust.setToFullTensionWithInnerSV(alpha);
                    full_tension_innerSV_p25 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p25,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
%                     spline_robust.lambda = lambda0;
%                     spline_robust.setToFullTensionWithSV(alpha);
%                     full_tension_sv_p25 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_sv_p25,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);

                    alpha = 1/8; % .125
                    spline_robust.lambda = lambda0;
                    spline_robust.setToFullTensionWithInnerSV(alpha);
                    full_tension_innerSV_p125 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_innerSV_p125,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);
%                     spline_robust.lambda = lambda0;
%                     spline_robust.setToFullTensionWithSV(alpha);
%                     full_tension_sv_p125 = LogStatisticsFromSplineForOutlierDetectionTable(full_tension_sv_p125,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices,outlierIndices);

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

printall = @(stats) fprintf('%.1f/%.1f (%.1f/%.1f)\n', mean(stats.false_positives(:)), mean(stats.false_negatives(:)), median(stats.false_positives(:)), median(stats.false_negatives(:)));
% printall(full_tension_sv_p6);
% printall(full_tension_sv_p5);
% printall(full_tension_sv_p4);
% printall(full_tension_sv_p25);
% printall(full_tension_sv_p125);

printall(full_tension_innerSV_p6);
printall(full_tension_innerSV_p5);
printall(full_tension_innerSV_p4);
printall(full_tension_innerSV_p25);
printall(full_tension_innerSV_p125);

printcol = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f/%.1f (%.1f/%.1f) ', mean(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), mean(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)) );

fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & optimal mse ($n_{eff}$) & blind optimal & robust & false pos/neg & robust 2nd & false pos/neg \\\\ \\hline \\hline \n');
for iOutlierRatio = 1:totalOutlierRatios
    for iSlope = 1:length(slopes)
        fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
        for iStride=1:length(strides)
            fprintf('%d ', strides(iStride));
            
%             printcol(full_tension_sv_p6,iOutlierRatio,iStride,iSlope);
%             printcol(full_tension_sv_p5,iOutlierRatio,iStride,iSlope);
%             printcol(full_tension_sv_p4,iOutlierRatio,iStride,iSlope);
%             printcol(full_tension_sv_p25,iOutlierRatio,iStride,iSlope);
%             printcol(full_tension_sv_p125,iOutlierRatio,iStride,iSlope);

            printcol(full_tension_innerSV_p6,iOutlierRatio,iStride,iSlope);
            printcol(full_tension_innerSV_p5,iOutlierRatio,iStride,iSlope);
            printcol(full_tension_innerSV_p4,iOutlierRatio,iStride,iSlope);
            printcol(full_tension_innerSV_p25,iOutlierRatio,iStride,iSlope);
            printcol(full_tension_innerSV_p125,iOutlierRatio,iStride,iSlope);
            
            fprintf(' \\\\ \n');
            
        end
    end
    fprintf('\n\n');
end
fprintf('\\end{tabular} \n');

