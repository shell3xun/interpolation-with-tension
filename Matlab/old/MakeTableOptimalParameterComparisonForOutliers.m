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
%     result_stride = 200;
    totalStrides = length(result_stride);
    totalEnsembles = 51; % best to choose an odd number for median

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
    outlierFactor = 50;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preallocate the variables we need to save
    %
    nothing = zeros(totalStrides, totalSlopes, totalEnsembles);
    nothing_struct = struct('mse',nothing,'neff_se',nothing,'false_negatives', nothing, 'false_positives', nothing);
    
    total_outliers = nothing; vars{end+1} = 'total_outliers';
    
%     optimal = nothing_struct; vars{end+1} = 'optimal';
%     robust_alpha100_optimal = nothing_struct; vars{end+1} = 'robust_alpha100_optimal';
%     robust_alpha100_knots_removed_optimal = nothing_struct; vars{end+1} = 'robust_alpha100_knots_removed_optimal';
    robust_alpha10k_optimal = nothing_struct; vars{end+1} = 'robust_alpha10k_optimal';
    robust_alpha10k_knots_removed_optimal = nothing_struct; vars{end+1} = 'robust_alpha10k_knots_removed_optimal';
    
    robust_alpha10k_blind_betaInf = nothing_struct; vars{end+1} = 'robust_alpha10k_blind_betaInf';
    robust_alpha10k_blind_beta10k = nothing_struct; vars{end+1} = 'robust_alpha10k_blind_beta10k';
    robust_alpha10k_blind_beta1k = nothing_struct; vars{end+1} = 'robust_alpha10k_blind_beta1k';
    robust_alpha10k_blind_beta100 = nothing_struct; vars{end+1} = 'robust_alpha10k_blind_beta100';
    robust_alpha10k_blind_beta10 = nothing_struct; vars{end+1} = 'robust_alpha10k_blind_beta10';
        
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
                truenonOutlierIndices = setdiff(1:n,trueOutlierIndices);
                
                total_outliers(iStride,iSlope,iEnsemble) = length(trueOutlierIndices);
                
                x_obs = data.x + epsilon;
                t_obs = data.t;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Unblinded best fit with standard tension spline
                
                linearIndex = sub2ind(size(nothing),iStride,iSlope,iEnsemble);
                
%                 spline_optimal = SmoothingSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected);
%                 spline_optimal.minimizeMeanSquareError(data.t,data.x);
%                 optimal = LogStatisticsFromSplineForOutlierTable(optimal,linearIndex,spline_optimal,compute_ms_error,trueOutlierIndices);
% 
%                 spline_robust_optimal = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'alpha',1/100);
%                 spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
%                 robust_alpha100_optimal = LogStatisticsFromSplineForOutlierTable(robust_alpha100_optimal,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
%                 
%                 spline_robust_optimal.removeOutlierKnotsAndRetension(1/10000);
%                 spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
%                 robust_alpha100_knots_removed_optimal = LogStatisticsFromSplineForOutlierTable(robust_alpha100_knots_removed_optimal,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
%                                
                spline_robust_optimal = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
                spline_robust_optimal.removeOutlierKnotsAndRetension(1/10000);
                spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
                robust_alpha10k_optimal = LogStatisticsFromSplineForOutlierTable(robust_alpha10k_optimal,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                
                beta = 1/10000;
                zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                spline_robust_optimal.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                robust_alpha10k_blind_beta10k = LogStatisticsFromSplineForOutlierTable(robust_alpha10k_blind_beta10k,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                
                beta = 1/1000;
                zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                spline_robust_optimal.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                robust_alpha10k_blind_beta1k = LogStatisticsFromSplineForOutlierTable(robust_alpha10k_blind_beta1k,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                
                beta = 1/100;
                zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                spline_robust_optimal.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                robust_alpha10k_blind_beta100 = LogStatisticsFromSplineForOutlierTable(robust_alpha10k_blind_beta100,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                
                beta = 1/10;
                zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                spline_robust_optimal.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                robust_alpha10k_blind_beta10 = LogStatisticsFromSplineForOutlierTable(robust_alpha10k_blind_beta10,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                
                spline_robust_optimal.minimizeExpectedMeanSquareError;
                robust_alpha10k_blind_betaInf = LogStatisticsFromSplineForOutlierTable(robust_alpha10k_blind_betaInf,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                
%                 spline_robust_optimal.removeOutlierKnotsAndRetension(1/10000);
%                 spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
%                 robust_alpha10k_knots_removed_optimal = LogStatisticsFromSplineForOutlierTable(robust_alpha10k_knots_removed_optimal,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                
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
% scatter(t_obs(spline_robust_optimal.outlierIndices),x_obs(spline_robust_optimal.outlierIndices),(8.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'), hold on
% scatter(t_obs(trueOutlierIndices),x_obs(trueOutlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
% scatter(t_obs,x_obs,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
% tq = linspace(min(t_obs),max(t_obs),10*length(t_obs));
% plot(tq,spline_optimal(tq),'k')
% plot(tq,spline_robust_optimal(tq),'b')
% % plot(tq,spline_robust_cv(tq),'m')
% return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% dmse = @(mse) mse./mse_robust_optimal-1;
% dmse_blind_optimal = dmse(mse_optimal);
% dmse_blind_ranged = dmse(mse_blind_optimal_ranged);
% dmse_robust_blind_alpha100 = dmse(mse_robust_blind_alpha100);
% dmse_robust_blind_alpha10000 = dmse(mse_robust_blind_alpha10000);
% dmse_robust_blind_alpha10000_rescale = dmse(mse_robust_blind_alpha10000_rescale);
% dmse_robust_blind_alpha10000_remove_knot = dmse(mse_robust_blind_alpha10000_remove_knot);
% 
% pct_range = 0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
% 
% minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
% maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

printcol = @(stats,iStride,iSlope) fprintf('& %#.3g/%#.3g m$^2$ (%#.3g) (%d/%d) ', median(stats.mse(iStride,iSlope,:)),mean(stats.mse(iStride,iSlope,:)), median(stats.neff_se(iStride,iSlope,:)), median(stats.false_positives(iStride,iSlope,:)), median(stats.false_negatives(iStride,iSlope,:)) );

dmse = @(mse) mse./robust_alpha10k_optimal.mse-1;
pct_range = 0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

robust_alpha10k_blind_beta10k.dmse = dmse(robust_alpha10k_blind_beta10k.mse);
robust_alpha10k_blind_beta1k.dmse = dmse(robust_alpha10k_blind_beta1k.mse);
robust_alpha10k_blind_beta100.dmse = dmse(robust_alpha10k_blind_beta100.mse);
robust_alpha10k_blind_beta10.dmse = dmse(robust_alpha10k_blind_beta10.mse);
robust_alpha10k_blind_betaInf.dmse = dmse(robust_alpha10k_blind_betaInf.mse);

print_pct = @(stats,iStride,iSlope) fprintf('&  %.1f-%.1f ',minpct(sort(stats.dmse(iStride,iSlope,:))), maxpct(sort(stats.dmse(iStride,iSlope,:))) );

fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & optimal mse ($n_{eff}$) & blind optimal & robust & false pos/neg & robust 2nd & false pos/neg \\\\ \\hline \\hline \n');
for iSlope = 1:length(slopes)
    fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
    for iStride=1:length(result_stride)
        fprintf('%d ', result_stride(iStride));
%         printcol(optimal,iStride,iSlope);
        printcol(robust_alpha10k_optimal,iStride,iSlope);
        print_pct(robust_alpha10k_blind_beta10k,iStride,iSlope);
        print_pct(robust_alpha10k_blind_beta1k,iStride,iSlope);
        print_pct(robust_alpha10k_blind_beta100,iStride,iSlope);
        print_pct(robust_alpha10k_blind_beta10,iStride,iSlope);
        print_pct(robust_alpha10k_blind_betaInf,iStride,iSlope);
        
%         printcol(robust_alpha10k_knots_removed_optimal,iStride,iSlope);
%         printcol(robust_alpha100_optimal,iStride,iSlope);
%         printcol(robust_alpha100_knots_removed_optimal,iStride,iSlope);
                
        fprintf(' \\\\ \n');
%         continue
%         
%         fprintf('&  %.1f-%.1f ',minpct(sort(dmse_blind_optimal(iStride,iSlope,:))), maxpct(sort(dmse_blind_optimal(iStride,iSlope,:))) );
%         
%         fprintf('&  %.1f-%.1f ',minpct(sort(dmse_blind_ranged(iStride,iSlope,:))), maxpct(sort(dmse_blind_ranged(iStride,iSlope,:))) );
%         fprintf('& %d (%d/%d)',median(total_outliers(iStride,iSlope,:)),median(false_positive_blind_optimal_ranged(iStride,iSlope,:)),median(false_negative_blind_optimal_ranged(iStride,iSlope,:)) );
%         
%         fprintf('&  %.1f-%.1f ',minpct(sort(dmse_robust_blind_alpha100(iStride,iSlope,:))), maxpct(sort(dmse_robust_blind_alpha100(iStride,iSlope,:))) );
%         fprintf('& %d (%d/%d)',median(total_outliers(iStride,iSlope,:)),median(false_positive_robust_blind_alpha100(iStride,iSlope,:)),median(false_negative_robust_blind_alpha100(iStride,iSlope,:)) );
%         
%         fprintf('&  %.1f-%.1f ',minpct(sort(dmse_robust_blind_alpha10000(iStride,iSlope,:))), maxpct(sort(dmse_robust_blind_alpha10000(iStride,iSlope,:))) );
%         fprintf('& %d (%d/%d)',median(total_outliers(iStride,iSlope,:)),median(false_positive_robust_blind_alpha10000(iStride,iSlope,:)),median(false_negative_robust_blind_alpha10000(iStride,iSlope,:)) );
%         
%                 
%         fprintf('&  %.1f-%.1f ',minpct(sort(dmse_robust_blind_alpha10000_rescale(iStride,iSlope,:))), maxpct(sort(dmse_robust_blind_alpha10000_rescale(iStride,iSlope,:))) );
%         fprintf('& %d (%d/%d) ',median(total_outliers(iStride,iSlope,:)),median(false_positive_robust_blind_alpha10000_rescale(iStride,iSlope,:)),median(false_negative_robust_blind_alpha10000_rescale(iStride,iSlope,:)) );
%         
%         fprintf('&  %.1f-%.1f ',minpct(sort(dmse_robust_blind_alpha10000_remove_knot(iStride,iSlope,:))), maxpct(sort(dmse_robust_blind_alpha10000_remove_knot(iStride,iSlope,:))) );
%         fprintf('& %d (%d/%d) \\\\ \n',median(total_outliers(iStride,iSlope,:)),median(false_positive_robust_blind_alpha10000_remove_knot(iStride,iSlope,:)),median(false_negative_robust_blind_alpha10000_remove_knot(iStride,iSlope,:)) );

    end
    
end
fprintf('\\end{tabular} \n');

