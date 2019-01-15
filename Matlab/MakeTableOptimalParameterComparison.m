scaleFactor = 1;
LoadFigureDefaults

shouldUseStudentTDistribution = 0;
shouldLoadExistingTable = 0;

if shouldUseStudentTDistribution == 1
    filename = 'MSEComparisonTableStudentT.mat';
else
    filename = 'MSEComparisonTable.mat';
end

if shouldLoadExistingTable == 1
    load(filename);
else
    slopes = [-2; -3; -4];
    totalSlopes = length(slopes);
    
    S = 2;
    T = S;
    K = S+1;
    
    result_stride = 2*[1;4;16;64];
    totalStrides = length(result_stride);
    
    totalEnsembles = 50;
    
%         result_stride = 2*[1;10;64];
%     totalStrides = length(result_stride);
%     
%     totalEnsembles = 5;
    
    % Do you want to assess the error using all the points from the signal
    % (which makes sense for an interpolation based metric) or just points from
    % the observed (decimated) signal only?
    shouldUseObservedSignalOnly = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preallocate the variables we need to save
    %
    % Table: rows are different strides, columns different optimal mse
    % columns: 1 knot dof,'auto' knot dof, initial guess, blind iteration
    u_rms_true = zeros(totalSlopes,1);
    a_rms_true = zeros(totalSlopes,1);
    
    u_rms_true_strided = zeros(totalStrides,totalSlopes);
    a_rms_true_strided = zeros(totalStrides,totalSlopes);
    
    u_estimate_spectral = zeros(totalStrides,totalSlopes, totalEnsembles);
    a_estimate_spectral = zeros(totalStrides,totalSlopes, totalEnsembles);
    expectedDOF = zeros(totalStrides,totalSlopes, totalEnsembles);
    
    reduced_dof_knot_dof = zeros(totalStrides, totalSlopes, totalEnsembles);
    
    mse_full_dof_true_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_true_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_blind_initial = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_blind_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_blind_optimal_ranged = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_log_likelihood = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_log_likelihood_blind = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_gcv= zeros(totalStrides, totalSlopes, totalEnsembles);
    
    dof_se_full_dof_true_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_true_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_blind_initial = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_blind_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_blind_optimal_ranged = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_log_likelihood = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_log_likelihood_blind = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_gcv = zeros(totalStrides, totalSlopes, totalEnsembles);
    
    for iSlope = 1:length(slopes)
        
        slope = slopes(iSlope);
        fprintf('slope %d, ',slope);
        
        if slope == -2
            data = load('sample_data/SyntheticTrajectories.mat');
            outputFile = 'OptimalParameters.mat';
        elseif slope == -3
            data = load('sample_data/SyntheticTrajectoriesSlope3.mat');
            outputFile = 'OptimalParametersSlope3.mat';
        elseif slope == -4
            data = load('sample_data/SyntheticTrajectoriesSlope4.mat');
            outputFile = 'OptimalParametersSlope4.mat';
        end
        sigma = data.position_error;
        dt = data.t(2)-data.t(1);
        u = diff(data.x)/dt;
        a = diff(diff(data.x))/(dt^2);
        u_rms_true(iSlope) = sqrt( mean( u.*u ) );
        a_rms_true(iSlope) = sqrt( mean( a.*a ) );
                
        for iStride=1:length(result_stride)
            stride = result_stride(iStride);
            
            % Reduce the total length in some cases
            if (stride < 10)
                shortenFactor = stride/10;
            else
                shortenFactor = 1;
            end
            
            indices = 1:stride:floor(shortenFactor*length(data.t));
            fprintf('dT = %d --- Using %d points with stride %d. Evaluating ensemble', stride*dt, length(indices), stride);
            if (shouldUseObservedSignalOnly == 1)
                indicesAll = indices;
            else
                indicesAll = 1:max(indices);
            end
            
            x_obs = data.x(indices);
            t_obs = data.t(indices);
            dt = t_obs(2)-t_obs(1);
            u = diff(x_obs)/dt;
            a = diff(diff(x_obs))/(dt^2);
            u_rms_true_strided(iStride,iSlope) = sqrt( mean( u.*u ) );
            a_rms_true_strided(iStride,iSlope) = sqrt( mean( a.*a ) );
            
            tensionDistribution = NormalDistribution(a_rms_true_strided(iStride,iSlope));
            
            for iEnsemble = 1:totalEnsembles
                fprintf('..%d',iEnsemble);
                
                if shouldUseStudentTDistribution == 1
                    nu = 4.5; sigma =  8.5;
                    noiseDistribution = StudentTDistribution(sigma,nu);
                    epsilon_x = randt(sigma,nu,length(indices));
                else
                    noiseDistribution = NormalDistribution(sigma);
                    epsilon_x = sigma*randn(length(indices),1);
                end 
                x_obs = data.x(indices) + epsilon_x;
                t_obs = data.t(indices);
                
                % Compute the expected tension before any fitting
                u_rms = TensionSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_obs,sqrt(noiseDistribution.variance),1);
                n_eff = TensionSpline.EffectiveSampleSizeFromUrms(u_rms, t_obs, sqrt(noiseDistribution.variance));
                
                % The idea here was to skip a bunch of unnecessary points
                % to trying to assess the noisy tension derivative... this
                % sort of works in that it no longer includes a bunch of
                % noise in the estimate, but it way under estimates the
                % power, especially in shallow slopes (where power is at
                % higher frequencies. The results is that we way over
                % tension.
                idx2 = 1:max(floor(n_eff/2),1):length(t_obs);
                idx2 = 1:length(t_obs);
                
                a_rms = TensionSpline.EstimateRMSDerivativeFromSpectrum(t_obs(idx2),x_obs(idx2),sqrt(noiseDistribution.variance),T);
                tensionDistributionEstimate = NormalDistribution(a_rms);
                lambda = (n_eff-1)/(n_eff*a_rms.^2);
                
                expectedDOF(iStride,iSlope,iEnsemble) = n_eff;
                u_estimate_spectral(iStride,iSlope,iEnsemble) = u_rms;
                a_estimate_spectral(iStride,iSlope,iEnsemble) = a_rms;

                % Fit using 1 dof for each data point
                spline_x = TensionSpline(t_obs,x_obs,noiseDistribution, 'lambda', lambda, 'S', S, 'T', T, 'knot_dof', 1);
                spline_x.minimizeMeanSquareError(data.t(indicesAll),data.x(indicesAll));
                compute_ms_error = @() (mean(mean(  (data.x(indicesAll) - spline_x(data.t(indicesAll))).^2,2 ),1));
                
                mse_full_dof_true_optimal(iStride,iSlope,iEnsemble) = compute_ms_error();
                dof_se_full_dof_true_optimal(iStride,iSlope,iEnsemble) = spline_x.effectiveSampleSizeFromVarianceOfTheMean;
                
                % Fit using the automatic knot dof algorithm
                spline_x = TensionSpline(t_obs,x_obs,noiseDistribution, 'lambda', lambda, 'S', S, 'T', T, 'knot_dof', 'auto');
                compute_ms_error = @() (mean(mean(  (data.x(indicesAll) - spline_x(data.t(indicesAll))).^2,2 ),1));
                reduced_dof_knot_dof(iStride,iSlope,iEnsemble) = spline_x.knot_dof;
                
                % record the initial fit
                mse_reduced_dof_blind_initial(iStride,iSlope,iEnsemble) = compute_ms_error();
                dof_se_reduced_dof_blind_initial(iStride,iSlope,iEnsemble) = spline_x.effectiveSampleSizeFromVarianceOfTheMean;
                
                % record the blind optimal fit
                spline_x.minimizeExpectedMeanSquareError;
                mse_reduced_dof_blind_optimal(iStride,iSlope,iEnsemble) = compute_ms_error();
                dof_se_reduced_dof_blind_optimal(iStride,iSlope,iEnsemble) = spline_x.effectiveSampleSizeFromVarianceOfTheMean;
                
                % record the blind optimal fit
                pctmin = 1/100/2;
                pctmax = 1-1/100/2;
                zmin = noiseDistribution.locationOfCDFPercentile(pctmin);
                zmax = noiseDistribution.locationOfCDFPercentile(pctmax);
                spline_x.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                mse_reduced_dof_blind_optimal_ranged(iStride,iSlope,iEnsemble) = compute_ms_error();
                dof_se_reduced_dof_blind_optimal_ranged(iStride,iSlope,iEnsemble) = spline_x.effectiveSampleSizeFromVarianceOfTheMean;
                
                % record the log-likehood optimal fit
                logLikelihood = @(spline) -sum(noiseDistribution.logPDF( spline.epsilon ) ) - sum(tensionDistribution.logPDF(spline.uniqueValuesAtHighestDerivative));
                spline_x.minimize( logLikelihood );
                mse_reduced_dof_log_likelihood(iStride,iSlope,iEnsemble) = compute_ms_error();
                dof_se_reduced_dof_log_likelihood(iStride,iSlope,iEnsemble) = spline_x.effectiveSampleSizeFromVarianceOfTheMean;
                
                % record the log-likehood blind fit
                logLikelihood = @(spline) -sum(noiseDistribution.logPDF( spline.epsilon ) ) - sum(tensionDistributionEstimate.logPDF(spline.uniqueValuesAtHighestDerivative));
                spline_x.minimize( logLikelihood );
                mse_reduced_dof_log_likelihood_blind(iStride,iSlope,iEnsemble) = compute_ms_error();
                dof_se_reduced_dof_log_likelihood_blind(iStride,iSlope,iEnsemble) = spline_x.effectiveSampleSizeFromVarianceOfTheMean;
                
                spline_x.minimize( @(spline) spline.expectedMeanSquareErrorFromGCV );
                mse_reduced_dof_gcv(iStride,iSlope,iEnsemble) = compute_ms_error();
                dof_se_reduced_dof_gcv(iStride,iSlope,iEnsemble) = spline_x.effectiveSampleSizeFromVarianceOfTheMean;
                
                % record the true optimal fit
                spline_x.minimizeMeanSquareError(data.t(indicesAll),data.x(indicesAll));
                mse_reduced_dof_true_optimal(iStride,iSlope,iEnsemble) = compute_ms_error();
                dof_se_reduced_dof_true_optimal(iStride,iSlope,iEnsemble) = spline_x.effectiveSampleSizeFromVarianceOfTheMean;
            end
            fprintf('\n');
        end
    end
    
    save(filename, 'S', 'T', 'slopes', 'result_stride', 'u_rms_true', 'a_rms_true', 'u_rms_true_strided','a_rms_true_strided', 'u_estimate_spectral', 'a_estimate_spectral', 'expectedDOF','reduced_dof_knot_dof','mse_full_dof_true_optimal','mse_reduced_dof_true_optimal','mse_reduced_dof_blind_initial','mse_reduced_dof_blind_optimal','mse_reduced_dof_blind_optimal_ranged','mse_reduced_dof_log_likelihood','mse_reduced_dof_log_likelihood_blind','mse_reduced_dof_gcv','dof_se_full_dof_true_optimal','dof_se_reduced_dof_true_optimal','dof_se_reduced_dof_blind_initial','dof_se_reduced_dof_blind_optimal','dof_se_reduced_dof_blind_optimal_ranged', 'dof_se_reduced_dof_log_likelihood', 'dof_se_reduced_dof_log_likelihood_blind','dof_se_reduced_dof_gcv');
end

reduced_dof_knot_dof_mean = mean(reduced_dof_knot_dof,3);

mse_full_dof_true_optimal_mean = mean(mse_full_dof_true_optimal,3);

dmse_reduced_dof_true_optimal = mse_reduced_dof_true_optimal./mse_full_dof_true_optimal-1;
dmse_reduced_dof_true_optimal_mean = median(dmse_reduced_dof_true_optimal,3);

dmse_reduced_dof_blind_initial = mse_reduced_dof_blind_initial./mse_full_dof_true_optimal - 1;
dmse_reduced_dof_blind_initial_mean = median(dmse_reduced_dof_blind_initial,3);

dmse_reduced_dof_blind_optimal = mse_reduced_dof_blind_optimal./mse_full_dof_true_optimal - 1;
dmse_reduced_dof_blind_optimal_mean = median(dmse_reduced_dof_blind_optimal,3);

dmse_reduced_dof_blind_optimal_ranged = mse_reduced_dof_blind_optimal_ranged./mse_full_dof_true_optimal - 1;
dmse_reduced_dof_blind_optimal_ranged_mean = median(dmse_reduced_dof_blind_optimal_ranged,3);

dmse_reduced_dof_log_likelihood = mse_reduced_dof_log_likelihood./mse_full_dof_true_optimal-1;
dmse_reduced_dof_log_likelihood_mean = median(dmse_reduced_dof_log_likelihood,3);

dmse_reduced_dof_log_likelihood_blind = mse_reduced_dof_log_likelihood_blind./mse_full_dof_true_optimal-1;
dmse_reduced_dof_log_likelihood_blind_mean = median(dmse_reduced_dof_log_likelihood_blind,3);

dmse_reduced_dof_gcv = mse_reduced_dof_gcv./mse_full_dof_true_optimal-1;
dmse_reduced_dof_gcv_mean = median(dmse_reduced_dof_gcv,3);

pct_range = 0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).

minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & full dof & reduced dof & blind initial & blind optimal & log-likelihood \\\\ \\hline \\hline \n');
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
        fprintf('&  %.1f-%.1f \\\\ \n',minpct(sort(dmse_reduced_dof_gcv(iStride,iSlope,:))),maxpct(sort(dmse_reduced_dof_gcv(iStride,iSlope,:)) ));
    end
    
end
fprintf('\\end{tabular} \n');

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