scaleFactor = 1;
LoadFigureDefaults

shouldUseStudentTDistribution = 0;

if shouldUseStudentTDistribution == 1
    filename = 'MSEComparisonTableStudentT.mat';
else
    filename = 'MSEComparisonTable.mat';
end

if exist(filename,'file')
    load(filename);
else
    slopes = [-2; -3; -4];
    totalSlopes = length(slopes);

%     result_stride = 2*[1;4;16;64];
%     totalStrides = length(result_stride);
%     totalEnsembles = 50;
    
    result_stride = 2*[1;10;64];
    totalStrides = length(result_stride);
    totalEnsembles = 5;
    
    
    % spline fit parameters
    S = 2;
    T = S;
    K = S+1;
    
    % matern signal parameters
    sigma_u = 0.20;
    base_dt = 5; % for whatever reason, we chose this as the primary dt
    t_damp = 30*60;
    n = 250;
    

    
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
    
    u_rms_true_strided = zeros(totalStrides,totalSlopes,totalEnsembles);
    a_rms_true_strided = zeros(totalStrides,totalSlopes,totalEnsembles);
    
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
                spline = TensionSpline(t_obs,x_obs,noiseDistribution, 'lambda', lambda, 'S', S, 'T', T, 'knot_dof', 1);
                spline.minimizeMeanSquareError(data.t,data.x);
                
                mse_full_dof_true_optimal(iStride,iSlope,iEnsemble) = compute_ms_error(spline);
                dof_se_full_dof_true_optimal(iStride,iSlope,iEnsemble) = spline.effectiveSampleSizeFromVarianceOfTheMean;
                
                % Fit using the automatic knot dof algorithm
                spline = TensionSpline(t_obs,x_obs,noiseDistribution, 'lambda', lambda, 'S', S, 'T', T, 'knot_dof', 'auto');
                reduced_dof_knot_dof(iStride,iSlope,iEnsemble) = spline.knot_dof;
                
                % record the initial fit
                mse_reduced_dof_blind_initial(iStride,iSlope,iEnsemble) = compute_ms_error(spline);
                dof_se_reduced_dof_blind_initial(iStride,iSlope,iEnsemble) = spline.effectiveSampleSizeFromVarianceOfTheMean;
                
                % record the blind optimal fit
                spline.minimizeExpectedMeanSquareError;
                mse_reduced_dof_blind_optimal(iStride,iSlope,iEnsemble) = compute_ms_error(spline);
                dof_se_reduced_dof_blind_optimal(iStride,iSlope,iEnsemble) = spline.effectiveSampleSizeFromVarianceOfTheMean;
                
                % record the blind optimal fit
                pctmin = 1/100/2;
                pctmax = 1-1/100/2;
                zmin = noiseDistribution.locationOfCDFPercentile(pctmin);
                zmax = noiseDistribution.locationOfCDFPercentile(pctmax);
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
                
                spline.lambda = lambda_full;
                spline.minimize( @(spline) spline.expectedMeanSquareErrorFromGCV );
                mse_reduced_dof_gcv(iStride,iSlope,iEnsemble) = compute_ms_error(spline);
                dof_se_reduced_dof_gcv(iStride,iSlope,iEnsemble) = spline.effectiveSampleSizeFromVarianceOfTheMean;
                
                % record the true optimal fit
                spline.lambda = lambda_full;
                spline.minimizeMeanSquareError(data.t,data.x);
                mse_reduced_dof_true_optimal(iStride,iSlope,iEnsemble) = compute_ms_error(spline);
                dof_se_reduced_dof_true_optimal(iStride,iSlope,iEnsemble) = spline.effectiveSampleSizeFromVarianceOfTheMean;
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