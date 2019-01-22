scaleFactor = 1;
LoadFigureDefaults

shouldUseStudentTDistribution = 1;

if shouldUseStudentTDistribution == 1
    filename = 'MSEComparisonTableForOutliersStudentT.mat';
else
    filename = 'MSEComparisonTableForOutliersNormal.mat';
end

if exist(filename,'file')
    load(filename);
else
    slopes = [-2; -3; -4];
    slopes = -3;
    totalSlopes = length(slopes);

    result_stride = [5;20;200];
    result_stride = 200;
    totalStrides = length(result_stride);
    totalEnsembles = 11; % best to choose an odd number for median

    % spline fit parameters
    S = 2;
    T = S;
    K = S+1;
    
    vars = {'S', 'T', 'slopes', 'result_stride'};
    
    % matern signal parameters
    sigma_u = 0.20;
    base_dt = 5; % for whatever reason, we chose this as the primary dt
    t_damp = 30*60;
    n = 250;

    % outlier parameters
    percentOutliers = 0.15;
    outlierFactor = 25;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preallocate the variables we need to save
    %
    
    total_outliers = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'total_outliers';
    
    mse_true_optimal= zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'mse_true_optimal';
    dof_se_true_optimal = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'dof_se_true_optimal';
    
    mse_blind_optimal= zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'mse_blind_optimal';
    dof_se_blind_optimal = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'dof_se_blind_optimal';
    
    mse_robust_blind_optimal_cv = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'mse_robust_blind_optimal_cv';
    dof_se_robust_blind_optimal_cv = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'dof_se_robust_blind_optimal_cv';
    false_positive_robust_blind_optimal_cv = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'false_positive_robust_blind_optimal_cv';
    false_negative_robust_blind_optimal_cv = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'false_negative_robust_blind_optimal_cv';
    
    mse_robust_blind_alpha100 = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'mse_robust_blind_alpha100';
    dof_se_robust_blind_alpha100 = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'dof_se_robust_blind_alpha100';
    false_positive_robust_blind_alpha100 = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'false_positive_robust_blind_alpha100';
    false_negative_robust_blind_alpha100 = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'false_negative_robust_blind_alpha100';
    
    mse_robust_blind_alpha10000 = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'mse_robust_blind_alpha10000';
    dof_se_robust_blind_alpha10000 = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'dof_se_robust_blind_alpha10000';
    false_positive_robust_blind_alpha10000 = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'false_positive_robust_blind_alpha10000';
    false_negative_robust_blind_alpha10000 = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'false_negative_robust_blind_alpha10000';    
    
    mse_robust_blind_alpha10000_rescale = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'mse_robust_blind_alpha10000_rescale';
    dof_se_robust_blind_alpha10000_rescale = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'dof_se_robust_blind_alpha10000_rescale';
    false_positive_robust_blind_alpha10000_rescale = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'false_positive_robust_blind_alpha10000_rescale';
    false_negative_robust_blind_alpha10000_rescale = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'false_negative_robust_blind_alpha10000_rescale';
    
    mse_robust_blind_alpha10000_remove_knot = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'mse_robust_blind_alpha10000_remove_knot';
    dof_se_robust_blind_alpha10000_remove_knot = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'dof_se_robust_blind_alpha10000_remove_knot';
    false_positive_robust_blind_alpha10000_remove_knot = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'false_positive_robust_blind_alpha10000_remove_knot';
    false_negative_robust_blind_alpha10000_remove_knot = zeros(totalStrides, totalSlopes, totalEnsembles); vars{end+1} = 'false_negative_robust_blind_alpha10000_remove_knot';
    
    for iSlope = 1:length(slopes)
        slope = slopes(iSlope);
        fprintf('slope %d, ',slope);
                
        for iStride=1:length(result_stride)
            stride = result_stride(iStride);
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
                else
                    sigma = 10;
                    noiseDistribution = NormalDistribution(sigma);
                end
                outlierIndices = rand(n,1)<=percentOutliers;
                outlierIndices([1 n]) = 0;
                
                epsilon_noise = zeros(n,1);
                epsilon_outlier = zeros(n,1);
                
                epsilon_noise(~outlierIndices) = noiseDistribution.rand(sum(~outlierIndices));
                epsilon_outlier(outlierIndices) = outlierFactor*noiseDistribution.rand(sum(outlierIndices));
                
                epsilon = epsilon_noise + epsilon_outlier;
                
                outlierThreshold = noiseDistribution.locationOfCDFPercentile(1-1/10000/2);
                trueOutlierIndices = find(abs(epsilon) > outlierThreshold);
                trueGoodIndices = setdiff(1:n,trueOutlierIndices);
                
                total_outliers(iStride,iSlope,iEnsemble) = length(trueOutlierIndices);
                
                x_obs = data.x + epsilon;
                t_obs = data.t;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Toss out the points that are outliers, and perform the best fit.
                spline_optimal = TensionSpline(t_obs(trueGoodIndices),x_obs(trueGoodIndices),noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected);
                spline_optimal.minimizeMeanSquareError(data.t,data.x);
                mse_true_optimal(iStride,iSlope,iEnsemble) = compute_ms_error(spline_optimal);
                dof_se_true_optimal(iStride,iSlope,iEnsemble) = spline_optimal.effectiveSampleSizeFromVarianceOfTheMean;
                
                spline_optimal.minimizeExpectedMeanSquareError();
                mse_blind_optimal(iStride,iSlope,iEnsemble) = compute_ms_error(spline_optimal);
                dof_se_blind_optimal(iStride,iSlope,iEnsemble) = spline_optimal.effectiveSampleSizeFromVarianceOfTheMean;
                
                
%                 spline_robust_cv = RobustTensionSpline(t_obs,x_obs,noiseDistribution, 'S', S);
%                 spline_robust_cv.firstIterationCV();
%                 mse_robust_blind_optimal_cv(iStride,iSlope,iEnsemble) = compute_ms_error(spline_robust_cv);
%                 dof_se_robust_blind_optimal_cv(iStride,iSlope,iEnsemble) = spline_robust_cv.effectiveSampleSizeFromVarianceOfTheMean;
%                 false_negative_robust_blind_optimal_cv(iStride,iSlope,iEnsemble) = length(setdiff(trueOutlierIndices,spline_robust_cv.indicesOfOutliers));
%                 false_positive_robust_blind_optimal_cv(iStride,iSlope,iEnsemble) = length(setdiff(spline_robust_cv.indicesOfOutliers,trueOutlierIndices));
                
                % Now repeat with the Robust Tension spline algorithm
                spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution, 'S', S);
                spline_robust.firstIteration(1/100);
                mse_robust_blind_alpha100(iStride,iSlope,iEnsemble) = compute_ms_error(spline_robust);
                dof_se_robust_blind_alpha100(iStride,iSlope,iEnsemble) = spline_robust.effectiveSampleSizeFromVarianceOfTheMean;
                false_negative_robust_blind_alpha100(iStride,iSlope,iEnsemble) = length(setdiff(trueOutlierIndices,spline_robust.indicesOfOutliers));
                false_positive_robust_blind_alpha100(iStride,iSlope,iEnsemble) = length(setdiff(spline_robust.indicesOfOutliers,trueOutlierIndices));
                
                spline_robust.firstIteration(1/10000);
                mse_robust_blind_alpha10000(iStride,iSlope,iEnsemble) = compute_ms_error(spline_robust);
                dof_se_robust_blind_alpha10000(iStride,iSlope,iEnsemble) = spline_robust.effectiveSampleSizeFromVarianceOfTheMean;
                false_negative_robust_blind_alpha10000(iStride,iSlope,iEnsemble) = length(setdiff(trueOutlierIndices,spline_robust.indicesOfOutliers));
                false_positive_robust_blind_alpha10000(iStride,iSlope,iEnsemble) = length(setdiff(spline_robust.indicesOfOutliers,trueOutlierIndices));
                   
%                 Repeat again after the second iteration is applied.
%                 spline_robust.rescaleDistributionAndRetension(1/10000);
%                 mse_robust_blind_alpha10000_rescale(iStride,iSlope,iEnsemble) = compute_ms_error(spline_robust);
%                 dof_se_robust_blind_alpha10000_rescale(iStride,iSlope,iEnsemble) = spline_robust.effectiveSampleSizeFromVarianceOfTheMean;
%                 false_negative_robust_blind_alpha10000_rescale(iStride,iSlope,iEnsemble) = length(setdiff(trueOutlierIndices,spline_robust.indicesOfOutliers));
%                 false_positive_robust_blind_alpha10000_rescale(iStride,iSlope,iEnsemble) = length(setdiff(spline_robust.indicesOfOutliers,trueOutlierIndices));
                
                spline_robust.firstIteration(1/10000);
                spline_robust.removeOutlierKnotsAndRetension(1/10000);
                mse_robust_blind_alpha10000_remove_knot(iStride,iSlope,iEnsemble) = compute_ms_error(spline_robust);
                dof_se_robust_blind_alpha10000_remove_knot(iStride,iSlope,iEnsemble) = spline_robust.effectiveSampleSizeFromVarianceOfTheMean;
                false_negative_robust_blind_alpha10000_remove_knot(iStride,iSlope,iEnsemble) = length(setdiff(trueOutlierIndices,spline_robust.indicesOfOutliers));
                false_positive_robust_blind_alpha10000_remove_knot(iStride,iSlope,iEnsemble) = length(setdiff(spline_robust.indicesOfOutliers,trueOutlierIndices));
            end
            fprintf('\n');
        end
    end
    
    save(filename, vars{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure

% figure
% scatter(t_obs(spline_robust.indicesOfOutliers),x_obs(spline_robust.indicesOfOutliers),(8.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'), hold on
% scatter(t_obs(trueOutlierIndices),x_obs(trueOutlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
% scatter(t_obs,x_obs,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
% tq = linspace(min(t_obs),max(t_obs),10*length(t_obs));
% plot(tq,spline_optimal(tq),'k')
% plot(tq,spline_robust(tq),'b')
% % plot(tq,spline_robust_cv(tq),'m')
% return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dmse = @(mse) mse./mse_true_optimal-1;
dmse_blind_optimal = dmse(mse_blind_optimal);
dmse_robust_blind_optimal_cv = dmse(mse_robust_blind_optimal_cv);
dmse_robust_blind_alpha100 = dmse(mse_robust_blind_alpha100);
dmse_robust_blind_alpha10000 = dmse(mse_robust_blind_alpha10000);
dmse_robust_blind_alpha10000_rescale = dmse(mse_robust_blind_alpha10000_rescale);
dmse_robust_blind_alpha10000_remove_knot = dmse(mse_robust_blind_alpha10000_remove_knot);

pct_range = 0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).

minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & optimal mse ($n_{eff}$) & blind optimal & robust & false pos/neg & robust 2nd & false pos/neg \\\\ \\hline \\hline \n');
for iSlope = 1:length(slopes)
    fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
    for iStride=1:length(result_stride)
        fprintf('%d ', result_stride(iStride));
        fprintf('& %#.3g m$^2$ (%#.3g) ', median(mse_true_optimal(iStride,iSlope,:)), median(dof_se_true_optimal(iStride,iSlope,:)) );
        fprintf('&  %.1f-%.1f ',minpct(sort(dmse_blind_optimal(iStride,iSlope,:))), maxpct(sort(dmse_blind_optimal(iStride,iSlope,:))) );
        
        fprintf('&  %.1f-%.1f ',minpct(sort(dmse_robust_blind_alpha100(iStride,iSlope,:))), maxpct(sort(dmse_robust_blind_alpha100(iStride,iSlope,:))) );
        fprintf('& %d (%d/%d)',median(total_outliers(iStride,iSlope,:)),median(false_positive_robust_blind_alpha100(iStride,iSlope,:)),median(false_negative_robust_blind_alpha100(iStride,iSlope,:)) );
        
        fprintf('&  %.1f-%.1f ',minpct(sort(dmse_robust_blind_alpha10000(iStride,iSlope,:))), maxpct(sort(dmse_robust_blind_alpha10000(iStride,iSlope,:))) );
        fprintf('& %d (%d/%d)',median(total_outliers(iStride,iSlope,:)),median(false_positive_robust_blind_alpha10000(iStride,iSlope,:)),median(false_negative_robust_blind_alpha10000(iStride,iSlope,:)) );
        
%         fprintf('&  %.1f-%.1f ',minpct(sort(dmse_robust_blind_optimal_cv(iStride,iSlope,:))), maxpct(sort(dmse_robust_blind_optimal_cv(iStride,iSlope,:))) );
%         fprintf('& %d (%d/%d)',median(total_outliers(iStride,iSlope,:)),median(false_positive_robust_blind_optimal_cv(iStride,iSlope,:)),median(false_negative_robust_blind_optimal_cv(iStride,iSlope,:)) );
        
%         fprintf('&  %.1f-%.1f ',minpct(sort(dmse_robust_blind_alpha10000_rescale(iStride,iSlope,:))), maxpct(sort(dmse_robust_blind_alpha10000_rescale(iStride,iSlope,:))) );
%         fprintf('& %d (%d/%d) ',median(total_outliers(iStride,iSlope,:)),median(false_positive_robust_blind_alpha10000_rescale(iStride,iSlope,:)),median(false_negative_robust_blind_alpha10000_rescale(iStride,iSlope,:)) );
        
        fprintf('&  %.1f-%.1f ',minpct(sort(dmse_robust_blind_alpha10000_remove_knot(iStride,iSlope,:))), maxpct(sort(dmse_robust_blind_alpha10000_remove_knot(iStride,iSlope,:))) );
        fprintf('& %d (%d/%d) \\\\ \n',median(total_outliers(iStride,iSlope,:)),median(false_positive_robust_blind_alpha10000_remove_knot(iStride,iSlope,:)),median(false_negative_robust_blind_alpha10000_remove_knot(iStride,iSlope,:)) );

    end
    
end
fprintf('\\end{tabular} \n');

