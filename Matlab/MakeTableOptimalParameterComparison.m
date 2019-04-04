scaleFactor = 1;
LoadFigureDefaults

shouldUseStudentTDistribution = 0;

if shouldUseStudentTDistribution == 1
    filename = 'MSEComparisonTableStudentT.mat';
else
    filename = 'MSEComparisonTableNormal.mat';
end

if exist(filename,'file')
    load(filename);
else
    slopes = [-2; -3; -4];
    slopes = -3;
    totalSlopes = length(slopes);

    strides = [5;20;80;200];
%     strides = [80;200];
%      strides = 200;
    totalStrides = length(strides);
    totalEnsembles = 101; % best to choose an odd number for median
    
    % spline fit parameters
    S = 2;
    T = S;
    K = S+1;
     
    varnames = {'S', 'T', 'slopes', 'strides', 'outlierRatios', 'totalEnsembles'};
    
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
    
    u_rms_true_strided = nothing; varnames{end+1} = 'u_rms_true_strided';
    a_rms_true_strided = nothing; varnames{end+1} = 'a_rms_true_strided';
    u_estimate_spectral = nothing; varnames{end+1} = 'u_estimate_spectral';
    a_estimate_spectral = nothing; varnames{end+1} = 'a_estimate_spectral';
    expected_neff = nothing; varnames{end+1} = 'expected_neff';
 
    stat_structs = cell(1,1);
    stat_structs{1} = nothing_struct; stat_structs{end}.name = 'minimization_noise_range';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'minimization_noise_range_sigma';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'minimization_ratio_1';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'minimization_ratio_1_sigma';
    varnames{end+1} = 'stat_structs';
           
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
                
                % record the 'true' (but strided) values of the signal
                rms = @(x) sqrt( mean( x.*x ) );
                u_rms_true_strided(iStride,iSlope,iEnsemble) = rms( diff(data.x)/dt );
                a_rms_true_strided(iStride,iSlope,iEnsemble) = rms( diff(diff(data.x))/(dt^2) );
                
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
                epsilon = noiseDistribution.rand(sum(~outlierIndices),1);
                
                x_obs = data.x + epsilon;
                t_obs = data.t;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Unblinded best fit with standard tension spline
                
                linearIndex = sub2ind(size(nothing),iOutlierRatio,iStride,iSlope,iEnsemble);
                
                alpha = 0;
                spline = TensionSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.optimalExpected);

                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Range-restricted expected mean square error.
                beta = 1/100;
                [lambda1,mse1] = spline.minimizeExpectedMeanSquareErrorInPercentileRange(beta/2,1-beta/2);
                %                     mse1_true = compute_ms_error(spline);
                
                iStruct = 1;
                stat_structs{iStruct} = LogStatisticsFromSplineForTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Range-restricted expected mean square error, with
                % different initial sigma.
                spline.sigma = spline.sigmaAtFullTension;
                spline.lambda = spline.lambdaAtFullTension;
                [~,mse2] = spline.minimizeExpectedMeanSquareErrorInPercentileRange(beta/2,1-beta/2);
                %                     mse2_true = compute_ms_error(spline);
                if mse1 < mse2
                    spline.sigma = sqrt(spline.distribution.variance);
                    spline.lambda = lambda1;
                    %                         if (mse1_true/mse2_true > 1.1)
                    %                             fprintf('Oops, bad choice. mse1 was worse by 10 percent.')
                    %                         end
                    %                     elseif (mse2_true/mse1_true > 1.1)
                    %                         fprintf('Oops, bad choice. mse2 was worse by 10 percent.')
                end
                
                iStruct = iStruct+1;
                stat_structs{iStruct} = LogStatisticsFromSplineForOutlierTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error,trueOutlierIndices,outlierIndices);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Range-restricted expected mean square error using the
                % pdf crossover point.
                spline.sigma = sqrt(spline.distribution.variance);
                spline.lambda = spline.lambdaAtFullTension;
                
                minimizationPDFRatio = 1;
                [lambda1,mse1] = spline.minimizeExpectedMeanSquareErrorInNoiseRange(minimizationPDFRatio);
                
                iStruct = iStruct+1;
                stat_structs{iStruct} = LogStatisticsFromSplineForOutlierTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error,trueOutlierIndices,outlierIndices);
                
                %                     fprintf('expected (%.1f, %.1f) actual (%.1f, %.1f)\n',mse0,mse1,stat_structs{iStruct-1}.mse(linearIndex),stat_structs{iStruct}.mse(linearIndex));
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Using the full tension values, search for an
                % alternative global minimum
                
                spline.sigma = spline.sigmaAtFullTension;
                spline.lambda = spline.lambdaAtFullTension;
                
                [~,mse2] = spline.minimizeExpectedMeanSquareErrorInNoiseRange(minimizationPDFRatio);
                if mse1 < mse2
                    spline.sigma = sqrt(spline.distribution.variance);
                    spline.lambda = lambda1;
                end
                
                iStruct = iStruct+1;
                stat_structs{iStruct} = LogStatisticsFromSplineForOutlierTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error,trueOutlierIndices,outlierIndices);
                
            end
            fprintf('\n');
        end
    end
    
    save(filename, varnames{:});

  
    
    reduced_dof_knot_dof = zeros(totalStrides, totalSlopes, totalEnsembles);
    
    mse_full_dof_true_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_true_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_blind_initial = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_blind_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_blind_optimal_ranged = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_log_likelihood = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_log_likelihood_blind = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_cv= zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_gcv= zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_robust_blind_optimal= zeros(totalStrides, totalSlopes, totalEnsembles);
    
    dof_se_full_dof_true_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_true_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_blind_initial = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_blind_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_blind_optimal_ranged = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_log_likelihood = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_log_likelihood_blind = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_cv = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_gcv = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_robust_blind_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    
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
                
                % record the 'true' (but strided) values of the signal
                rms = @(x) sqrt( mean( x.*x ) );
                u_rms_true_strided(iStride,iSlope,iEnsemble) = rms( diff(data.x)/dt );
                a_rms_true_strided(iStride,iSlope,iEnsemble) = rms( diff(diff(data.x))/(dt^2) );
                
                tensionDistribution = NormalDistribution(a_rms_true_strided(iStride,iSlope,iEnsemble));
                
                % Generate the noise
                if shouldUseStudentTDistribution == 1
                    nu = 4.5; sigma =  8.5;
                    noiseDistribution = StudentTDistribution(sigma,nu);
                else
                    sigma = 10;
                    noiseDistribution = NormalDistribution(sigma);
                end 
                epsilon_x = noiseDistribution.rand(n);
                x_obs = data.x + epsilon_x;
                t_obs = data.t;
                
                u = diff(x_obs)/dt;
                a = diff(diff(x_obs))/(dt^2);
            
                % Compute the expected values from the contaminated data
                u_rms = TensionSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_obs,sqrt(noiseDistribution.variance),1);
                a_rms = TensionSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_obs,sqrt(noiseDistribution.variance),T);
                n_eff = TensionSpline.EffectiveSampleSizeFromUrms(u_rms,t_obs, sqrt(noiseDistribution.variance));
                
                % record
                u_estimate_spectral(iStride,iSlope,iEnsemble) = u_rms;
                a_estimate_spectral(iStride,iSlope,iEnsemble) = a_rms;
                expectedDOF(iStride,iSlope,iEnsemble) = n_eff;
                
                tensionDistributionEstimate = NormalDistribution(a_rms);
                lambda = (n_eff-1)/(n_eff*a_rms.^2);
                lambda_full = 1/a_rms^2;
                
                % Fit using 1 dof for each data point
                spline_full = TensionSpline(t_obs,x_obs,noiseDistribution, 'lambda', lambda_full, 'S', S, 'T', T, 'knot_dof', 1);
                spline_full.minimizeMeanSquareError(data.t,data.x); 
                mse_full_dof_true_optimal(iStride,iSlope,iEnsemble) = compute_ms_error(spline_full);
                dof_se_full_dof_true_optimal(iStride,iSlope,iEnsemble) = spline_full.effectiveSampleSizeFromVarianceOfTheMean;
                
                % Fit using the automatic knot dof algorithm
                spline = TensionSpline(t_obs,x_obs,noiseDistribution, 'lambda', lambda, 'S', S, 'T', T, 'knot_dof', 'auto');
                reduced_dof_knot_dof(iStride,iSlope,iEnsemble) = spline.knot_dof;
                
                % record the initial fit
                mse_reduced_dof_blind_initial(iStride,iSlope,iEnsemble) = compute_ms_error(spline);
                dof_se_reduced_dof_blind_initial(iStride,iSlope,iEnsemble) = spline.effectiveSampleSizeFromVarianceOfTheMean;
                
                % record the blind optimal fit
                spline.lambda = lambda_full;
                spline.minimizeExpectedMeanSquareError;
                mse_reduced_dof_blind_optimal(iStride,iSlope,iEnsemble) = compute_ms_error(spline);
                dof_se_reduced_dof_blind_optimal(iStride,iSlope,iEnsemble) = spline.effectiveSampleSizeFromVarianceOfTheMean;
                
%                 dmse = mse_reduced_dof_blind_optimal(iStride,iSlope,iEnsemble)/mse_full_dof_true_optimal(iStride,iSlope,iEnsemble);
%                 if dmse > 1.2
%                    fprintf('dmse %f\n',dmse);
%                    figure
%                    tq = linspace(min(spline.t),max(spline.t),10*length(spline.t));
%                    scatter(t_obs,x_obs,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
%                    plot(tq,spline_full(tq))
%                    plot(tq,spline(tq))
%                 end
                
                % record the blind optimal fit
                pctmin = 1/100/2;
                pctmax = 1-1/100/2;
                zmin = noiseDistribution.locationOfCDFPercentile(pctmin);
                zmax = noiseDistribution.locationOfCDFPercentile(pctmax);
                spline.lambda = lambda_full;
                spline.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                mse_reduced_dof_blind_optimal_ranged(iStride,iSlope,iEnsemble) = compute_ms_error(spline);
                dof_se_reduced_dof_blind_optimal_ranged(iStride,iSlope,iEnsemble) = spline.effectiveSampleSizeFromVarianceOfTheMean;
                
                % record the log-likehood optimal fit
                logLikelihood = @(spline) -sum(noiseDistribution.logPDF( spline.epsilon ) ) - sum(tensionDistribution.logPDF(spline.uniqueValuesAtHighestDerivative));
                spline.lambda = lambda_full;
                spline.minimize( logLikelihood );
                mse_reduced_dof_log_likelihood(iStride,iSlope,iEnsemble) = compute_ms_error(spline);
                dof_se_reduced_dof_log_likelihood(iStride,iSlope,iEnsemble) = spline.effectiveSampleSizeFromVarianceOfTheMean;
                
                % record the log-likehood blind fit
                logLikelihood = @(spline) -sum(noiseDistribution.logPDF( spline.epsilon ) ) - sum(tensionDistributionEstimate.logPDF(spline.uniqueValuesAtHighestDerivative));
                spline.lambda = lambda_full;
                spline.minimize( logLikelihood );
                mse_reduced_dof_log_likelihood_blind(iStride,iSlope,iEnsemble) = compute_ms_error(spline);
                dof_se_reduced_dof_log_likelihood_blind(iStride,iSlope,iEnsemble) = spline.effectiveSampleSizeFromVarianceOfTheMean;
                
                % record the cross-validation fit
                spline.lambda = lambda_full;
                spline.minimize( @(spline) spline.expectedMeanSquareErrorFromCV );
                mse_reduced_dof_cv(iStride,iSlope,iEnsemble) = compute_ms_error(spline);
                dof_se_reduced_dof_cv(iStride,iSlope,iEnsemble) = spline.effectiveSampleSizeFromVarianceOfTheMean;
                
                % record the generalized cross-validation fit
                spline.lambda = lambda_full;
                spline.minimize( @(spline) spline.expectedMeanSquareErrorFromGCV );
                mse_reduced_dof_gcv(iStride,iSlope,iEnsemble) = compute_ms_error(spline);
                dof_se_reduced_dof_gcv(iStride,iSlope,iEnsemble) = spline.effectiveSampleSizeFromVarianceOfTheMean;
                
                % record the true optimal fit
                spline.lambda = lambda_full;
                spline.minimizeMeanSquareError(data.t,data.x);
                mse_reduced_dof_true_optimal(iStride,iSlope,iEnsemble) = compute_ms_error(spline);
                dof_se_reduced_dof_true_optimal(iStride,iSlope,iEnsemble) = spline.effectiveSampleSizeFromVarianceOfTheMean;
                
                % Now repeat with the Robust Tension spline algorithm
                spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'T', T);
                mse_robust_blind_optimal(iStride,iSlope,iEnsemble) = compute_ms_error(spline_robust);
                dof_se_robust_blind_optimal(iStride,iSlope,iEnsemble) = spline_robust.effectiveSampleSizeFromVarianceOfTheMean;
            end
            fprintf('\n');
        end
    end
    
    save(filename, 'S', 'T', 'slopes', 'result_stride', 'u_rms_true', 'a_rms_true', 'u_rms_true_strided','a_rms_true_strided', 'u_estimate_spectral', 'a_estimate_spectral', 'expectedDOF','reduced_dof_knot_dof','mse_full_dof_true_optimal','mse_reduced_dof_true_optimal','mse_reduced_dof_blind_initial','mse_reduced_dof_blind_optimal','mse_reduced_dof_blind_optimal_ranged','mse_reduced_dof_log_likelihood','mse_reduced_dof_log_likelihood_blind','mse_reduced_dof_cv','mse_reduced_dof_gcv','mse_robust_blind_optimal','dof_se_full_dof_true_optimal','dof_se_reduced_dof_true_optimal','dof_se_reduced_dof_blind_initial','dof_se_reduced_dof_blind_optimal','dof_se_reduced_dof_blind_optimal_ranged', 'dof_se_reduced_dof_log_likelihood', 'dof_se_reduced_dof_log_likelihood_blind','dof_se_reduced_dof_cv','dof_se_reduced_dof_gcv','dof_se_robust_blind_optimal');
end

reduced_dof_knot_dof_mean = mean(reduced_dof_knot_dof,3);

mse_absolute_optimal = min(mse_full_dof_true_optimal,mse_reduced_dof_true_optimal);

dmse_full_dof_true_optimal = mse_full_dof_true_optimal./mse_absolute_optimal-1;
mse_full_dof_true_optimal_mean = mean(mse_full_dof_true_optimal,3);

dmse_reduced_dof_true_optimal = mse_reduced_dof_true_optimal./mse_absolute_optimal-1;
dmse_reduced_dof_true_optimal_mean = median(dmse_reduced_dof_true_optimal,3);

dmse_reduced_dof_blind_initial = mse_reduced_dof_blind_initial./mse_absolute_optimal - 1;
dmse_reduced_dof_blind_initial_mean = median(dmse_reduced_dof_blind_initial,3);

dmse_reduced_dof_blind_optimal = mse_reduced_dof_blind_optimal./mse_absolute_optimal - 1;
dmse_reduced_dof_blind_optimal_mean = median(dmse_reduced_dof_blind_optimal,3);

dmse_reduced_dof_blind_optimal_ranged = mse_reduced_dof_blind_optimal_ranged./mse_absolute_optimal - 1;
dmse_reduced_dof_blind_optimal_ranged_mean = median(dmse_reduced_dof_blind_optimal_ranged,3);

dmse_reduced_dof_log_likelihood = mse_reduced_dof_log_likelihood./mse_absolute_optimal-1;
dmse_reduced_dof_log_likelihood_mean = median(dmse_reduced_dof_log_likelihood,3);

dmse_reduced_dof_log_likelihood_blind = mse_reduced_dof_log_likelihood_blind./mse_absolute_optimal-1;
dmse_reduced_dof_log_likelihood_blind_mean = median(dmse_reduced_dof_log_likelihood_blind,3);

dmse_reduced_dof_cv = mse_reduced_dof_cv./mse_absolute_optimal-1;
dmse_reduced_dof_cv_mean = median(dmse_reduced_dof_cv,3);

dmse_reduced_dof_gcv = mse_reduced_dof_gcv./mse_absolute_optimal-1;
dmse_reduced_dof_gcv_mean = median(dmse_reduced_dof_gcv,3);

dmse_robust_blind_optimal = mse_robust_blind_optimal./mse_absolute_optimal-1;
dmse_robust_blind_optimal_mean = median(dmse_robust_blind_optimal,3);

pct_range = 0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).

minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & optimal mse ($n_{eff}$) & reduced dof & blind initial & blind optimal & blind optimal ranged & log-likelihood & blind log-likelihood & GCV \\\\ \\hline \\hline \n');
for iSlope = 1:length(slopes)
    fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
    for iStride=1:length(result_stride)
        fprintf('%d ', result_stride(iStride));
        fprintf('& %#.3g m$^2$ (%#.3g) ', mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_full_dof_true_optimal(iStride,iSlope));
        fprintf('&  %.1f-%.1f ',minpct(sort(dmse_reduced_dof_true_optimal(iStride,iSlope,:))),maxpct(sort(dmse_reduced_dof_true_optimal(iStride,iSlope,:)) ));
        fprintf('&  %.1f-%.1f ',minpct(sort(dmse_reduced_dof_blind_initial(iStride,iSlope,:))),maxpct(sort(dmse_reduced_dof_blind_initial(iStride,iSlope,:)) ));
        fprintf('&  %.1f-%.1f ',minpct(sort(dmse_reduced_dof_blind_optimal(iStride,iSlope,:))),maxpct(sort(dmse_reduced_dof_blind_optimal(iStride,iSlope,:)) ));
        fprintf('&  %.1f-%.1f ',minpct(sort(dmse_reduced_dof_blind_optimal_ranged(iStride,iSlope,:))),maxpct(sort(dmse_reduced_dof_blind_optimal_ranged(iStride,iSlope,:)) ));
        fprintf('&  %.1f-%.1f ',minpct(sort(dmse_reduced_dof_log_likelihood(iStride,iSlope,:))),maxpct(sort(dmse_reduced_dof_log_likelihood(iStride,iSlope,:)) ));
        fprintf('&  %.1f-%.1f ',minpct(sort(dmse_reduced_dof_log_likelihood_blind(iStride,iSlope,:))),maxpct(sort(dmse_reduced_dof_log_likelihood_blind(iStride,iSlope,:)) ));
        fprintf('&  %.1f-%.1f ',minpct(sort(dmse_reduced_dof_cv(iStride,iSlope,:))),maxpct(sort(dmse_reduced_dof_cv(iStride,iSlope,:)) ));
        fprintf('&  %.1f-%.1f ',minpct(sort(dmse_reduced_dof_gcv(iStride,iSlope,:))),maxpct(sort(dmse_reduced_dof_gcv(iStride,iSlope,:)) ));
        fprintf('&  %.1f-%.1f \\\\ \n',minpct(sort(dmse_robust_blind_optimal(iStride,iSlope,:))),maxpct(sort(dmse_robust_blind_optimal(iStride,iSlope,:)) ));

    end
    
end
fprintf('\\end{tabular} \n');


minpct = @(values) values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) values(floor( ((1+pct_range)/2)*length(values)));
makeErrorPlot = @(x,data) errorbar(x,mean(data),mean(data)-minpct(sort(data)),maxpct(sort(data))-mean(data),'o','MarkerSize',5,'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
makeErrorArray = @(x,data) [x,mean(data),mean(data)-minpct(sort(data)),maxpct(sort(data))-mean(data)];
figure
for iSlope = 1:length(slopes)
    subplot(1,length(slopes),iSlope)
    for iStride=1:length(result_stride)      
        makeErrorPlot(1,mse_full_dof_true_optimal(iStride,iSlope,:)), hold on,
        makeErrorPlot(2,mse_reduced_dof_true_optimal(iStride,iSlope,:))
        makeErrorPlot(3,mse_reduced_dof_cv(iStride,iSlope,:))
        makeErrorPlot(4,mse_reduced_dof_gcv(iStride,iSlope,:))
        makeErrorPlot(5,mse_reduced_dof_blind_optimal(iStride,iSlope,:))
        makeErrorPlot(6,mse_robust_blind_optimal(iStride,iSlope,:))
        makeErrorPlot(7,mse_reduced_dof_blind_initial(iStride,iSlope,:))
        makeErrorPlot(8,mse_reduced_dof_log_likelihood(iStride,iSlope,:))
    end
    title(sprintf('slope %d',slopes(iSlope)))
    xlim([0 9])
end
packfig(1,length(slopes))

minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));
makeErrorPlot = @(x,data) errorbar(x,100*mean(data),100*mean(data)-minpct(sort(data)),maxpct(sort(data))-100*mean(data),'o','MarkerSize',5,'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
makeErrorArray = @(x,data) [x,mean(data),mean(data)-minpct(sort(data)),maxpct(sort(data))-mean(data)];

figure
for iSlope = 1:length(slopes)
    for iStride=1:length(result_stride)
        subplot(length(result_stride),length(slopes),(length(result_stride)-iStride)*length(slopes) + iSlope)
        makeErrorPlot(1,dmse_full_dof_true_optimal(iStride,iSlope,:)), hold on
        makeErrorPlot(2,dmse_reduced_dof_true_optimal(iStride,iSlope,:)), hold on
        makeErrorPlot(3,dmse_reduced_dof_blind_optimal(iStride,iSlope,:))
        makeErrorPlot(4,dmse_reduced_dof_blind_optimal_ranged(iStride,iSlope,:))
        makeErrorPlot(5,dmse_robust_blind_optimal(iStride,iSlope,:))
        
        makeErrorPlot(6,dmse_reduced_dof_cv(iStride,iSlope,:))
        makeErrorPlot(7,dmse_reduced_dof_gcv(iStride,iSlope,:))
        
        makeErrorPlot(8,dmse_reduced_dof_blind_initial(iStride,iSlope,:))
        makeErrorPlot(9,dmse_reduced_dof_log_likelihood(iStride,iSlope,:))
        makeErrorPlot(10,dmse_reduced_dof_log_likelihood_blind(iStride,iSlope,:))
        xlim([0 11])
        ylim([0 100])
        if iStride == length(result_stride)
           title(sprintf('slope %d',slopes(iSlope))) 
        end
        if iSlope == 1
           ylabel(sprintf('stride %d',result_stride(iStride))) 
        end
    end
    
end
packfig(length(result_stride),length(slopes))
return

% analysis of dof
% a = dof_se_reduced_dof_log_likelihood./dof_se_reduced_dof_true_optimal - 1;
% iStride=1; iSlope = 1; figure, histogram(a(iStride,iSlope,:))
% suggests that we underestimate dof at small strides, over estimate at
% large strides. big problem at w^-2, not a problem at omega^-4

fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & full dof & reduced dof & blind initial & blind optimal & log-likelihood \\\\ \\hline \\hline \n');
for iSlope = 1:length(slopes)
    
    fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
    for iStride=1:length(result_stride)
%         fprintf('%d & %#.3g m^2 (%#.3g) &  %#.3g m^2 (%#.3g) &  %#.3g m^2 (%#.3g) &  %#.3g m^2 (%#.3g) \\\\ \n', result_stride(iStride), mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_full_dof_true_optimal(iStride,iSlope), dmse_reduced_dof_true_optimal_mean(iStride,iSlope), dof_se_reduced_dof_true_optimal(iStride,iSlope), dmse_reduced_dof_blind_initial_mean(iStride,iSlope), dof_se_reduced_dof_blind_initial(iStride,iSlope), dmse_reduced_dof_blind_optimal_mean(iStride,iSlope), dof_se_reduced_dof_blind_optimal(iStride,iSlope) )  ;
        fprintf('%d & %#.3g m$^2$ (%#.3g) &  %+.1f\\%% (%#.3g) &  %+.1f\\%% (%#.3g) &  %+.1f\\%% (%#.3g) &  %+.1f\\%% (%#.3g) \\\\ \n', result_stride(iStride), mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_full_dof_true_optimal(iStride,iSlope), 100*dmse_reduced_dof_true_optimal_mean(iStride,iSlope), dof_se_reduced_dof_true_optimal(iStride,iSlope), 100*dmse_reduced_dof_blind_initial_mean(iStride,iSlope), dof_se_reduced_dof_blind_initial(iStride,iSlope), 100*dmse_reduced_dof_blind_optimal_mean(iStride,iSlope), dof_se_reduced_dof_blind_optimal(iStride,iSlope), 100*dmse_reduced_dof_gcv(iStride,iSlope), dof_se_reduced_dof_gcv(iStride,iSlope) )  ;
    end
    
end
fprintf('\\end{tabular} \n');

fprintf('\n\n');
fprintf('\\begin{tabular}{r | lll} stride & full dof & optimal log-likelihood & blind log-likelihood \\\\ \\hline \\hline \n');
for iSlope = 1:length(slopes)
    
    fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
    for iStride=1:length(result_stride)
        fprintf('%d & %#.3g m$^2$ (%#.3g) &  %+.1f\\%% (%#.3g) &  %+.1f\\%% (%#.3g) \\\\ \n', result_stride(iStride), mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_full_dof_true_optimal(iStride,iSlope), 100*dmse_reduced_dof_log_likelihood_mean(iStride,iSlope), dof_se_reduced_dof_log_likelihood(iStride,iSlope), 100*dmse_reduced_dof_log_likelihood_blind_mean(iStride,iSlope), dof_se_reduced_dof_log_likelihood_blind(iStride,iSlope) )  ;
    end
    
end
fprintf('\\end{tabular} \n');


% Note: my estimate of u_rms is pretty good,
% mean(u_estimate_spectral,3)./u_rms_true'

% but my estimate of a_rms is pretty bad,
%mean(a_estimate_spectral,3)./a_rms_true'