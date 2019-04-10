scaleFactor = 1;
LoadFigureDefaults

shouldUseStudentTDistribution = 1;

if shouldUseStudentTDistribution == 1
    filename = 'MSEComparisonTableStudentT.mat';
else
    filename = 'MSEComparisonTableNormal.mat';
end

if exist(filename,'file')
    load(filename);
else
    slopes = [-2; -3; -4];
    strides = (2.^(0:4)).';
    
    totalSlopes = length(slopes);
    totalStrides = length(strides);
    totalEnsembles = 200; % best to choose an odd number for median
    
    % spline fit parameters
    S = 3;
    T = S;
    K = S+1;
    
    varnames = {'S', 'T', 'slopes', 'strides', 'totalEnsembles'};
    
    % matern signal parameters
    sigma_u = 0.20;
    base_dt = 1.5*60; % chosen as the smallest interval considered, because anything shorter than this looks non-stationary... like a local polynomial fit is needed.
    t_damp = 30*60;
    n = 250;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preallocate the variables we need to save
    %
    nothing = nan(totalStrides, totalSlopes, totalEnsembles);
    nothing_struct = struct('mse',nothing,'neff_se',nothing,'lambda',nothing);
    
    u_rms_true_strided = nothing; varnames{end+1} = 'u_rms_true_strided';
    a_rms_true_strided = nothing; varnames{end+1} = 'a_rms_true_strided';
    u_estimate_spectral = nothing; varnames{end+1} = 'u_estimate_spectral';
    a_estimate_spectral = nothing; varnames{end+1} = 'a_estimate_spectral';
    
    stat_structs = cell(1,1);
    stat_structs{1} = nothing_struct; stat_structs{end}.name = 'true_optimal_knot_dof_1';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'optimal_expected';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'true_optimal_knot_dof_auto';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'optimal_iterated';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'optimal_ranged';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'optimal_cv';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'optimal_gcv';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'optimal_log_likelihood';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'optimal_iterated_mean_removed';
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
                epsilon = noiseDistribution.rand(size(data.x));
                
                x_obs = data.x + epsilon;
                t_obs = data.t;
                                
                % record the estimated values of the signal
                u_estimate_spectral(iStride,iSlope,iEnsemble) = TensionSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_obs,sqrt(noiseDistribution.variance),1);
                a_estimate_spectral(iStride,iSlope,iEnsemble) = TensionSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_obs,sqrt(noiseDistribution.variance),T);
                
                linearIndex = sub2ind(size(nothing),iStride,iSlope,iEnsemble);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Unblinded best fit with knot_dof==1
                spline = TensionSpline(t_obs,x_obs,noiseDistribution, 'lambda', Lambda.fullTensionExpected, 'S', S, 'knot_dof', 1);
                spline.minimizeMeanSquareError(data.t,data.x);
                
                iStruct = 1;
                stat_structs{iStruct} = LogStatisticsFromSplineForTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Automatic knot selection, optimal expected
                spline = TensionSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.optimalExpected, 'knot_dof', 'auto');
                
                iStruct = iStruct+1; % 2
                stat_structs{iStruct} = LogStatisticsFromSplineForTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Automatic knot selection, optimal actual
                spline = TensionSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected, 'knot_dof', 'auto');
                lambda_full = spline.lambda;
                spline.minimizeMeanSquareError(data.t,data.x);
                
                iStruct = iStruct+1; % 3
                stat_structs{iStruct} = LogStatisticsFromSplineForTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Automatic knot selection, optimal iterated
                spline.lambda = lambda_full;
                spline.minimizeExpectedMeanSquareError();
                
                iStruct = iStruct+1; % 4
                stat_structs{iStruct} = LogStatisticsFromSplineForTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Range-restricted expected mean square error.
                beta = 1/100;
                spline.lambda = lambda_full;
                spline.minimizeExpectedMeanSquareErrorInPercentileRange(beta/2,1-beta/2);
                
                iStruct = iStruct+1; % 5
                stat_structs{iStruct} = LogStatisticsFromSplineForTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % cross-validation
                spline.lambda = lambda_full;
                spline.minimize( @(spline) spline.expectedMeanSquareErrorFromCV );
                
                iStruct = iStruct+1; % 6
                stat_structs{iStruct} = LogStatisticsFromSplineForTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % generalized cross-validation
                spline.lambda = lambda_full;
                spline.minimize( @(spline) spline.expectedMeanSquareErrorFromGCV );
                
                iStruct = iStruct+1; % 7
                stat_structs{iStruct} = LogStatisticsFromSplineForTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % blind log-likelihood
                x_filtered = TensionSpline.RunningFilter(x_obs,11,'median');
                a_rms = TensionSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_filtered,sqrt(noiseDistribution.variance),T);
                tensionDistributionEstimate = NormalDistribution(a_rms);
                logLikelihood = @(spline) -sum(noiseDistribution.logPDF( spline.epsilon ) ) - sum(tensionDistributionEstimate.logPDF(spline.uniqueValuesAtHighestDerivative));
                spline.lambda = lambda_full;
                spline.minimize( logLikelihood );
                
                iStruct = iStruct+1; % 8
                stat_structs{iStruct} = LogStatisticsFromSplineForTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error);
                                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Fit a constrained tension spline (a polynomial fit), then
                % an optimal iterated spline. For stationary data, this
                % should make no appreciable difference.
                K = S+1;
                t_knot = cat(1,min(t_obs)*ones(K+1,1),max(t_obs)*ones(K+1,1));
                mean_spline = ConstrainedSpline(t_obs,x_obs,K+1,t_knot,noiseDistribution,[]);
                
                spline = TensionSpline(t_obs,x_obs-mean_spline(t_obs),noiseDistribution, 'S', S, 'lambda',Lambda.optimalIterated, 'knot_dof', 'auto');
                compute_ms_error = @(spline) (mean(mean(  (data.x - spline(data.t) - mean_spline(data.t)).^2,2 ),1));
                
                iStruct = iStruct+1; % 9
                stat_structs{iStruct} = LogStatisticsFromSplineForTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error);
                
            end
            fprintf('\n');
        end
    end
    
    save(filename, varnames{:});
end

% All d-mse values are given relative to the knot_dof=1 optimal solution
dmse = @(stats) log10(stats.mse./stat_structs{1}.mse);
for i=1:length(stat_structs)
    stat_structs{i}.dmse = dmse(stat_structs{i});
end

% This converts dmse into a percentile range
pct_range = 0.9; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
minpct = @(values) sign(values(ceil( ((1-pct_range)/2)*length(values))))*100*(10^abs(values(ceil( ((1-pct_range)/2)*length(values))))-1);
maxpct = @(values) sign(values(floor( ((1+pct_range)/2)*length(values))))*100*(10^abs(values(floor( ((1+pct_range)/2)*length(values))))-1);

% Show the aggregate total dmse for each solution
print_dmse_all = @(stats) fprintf('%s: & %.2f (%.2f) (%.2f-%.2f)\n',stats.name,100*(10^mean(stats.dmse(:))-1),100*(10^median(stats.dmse(:))-1),minpct(sort(stats.dmse(:))), maxpct(sort(stats.dmse(:))) );
fprintf('---dmse---\n');
for i=1:length(stat_structs)
    print_dmse_all( stat_structs{i} );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MSE Comparison table for manuscript
%
print_mse = @(stats,iStride,iSlope) sprintf('%#.3g m$^2$',mean(reshape(stats.mse(iStride,iSlope,:),[],1)));
print_dmse = @(stats,iStride,iSlope) sprintf('%.1f\\%%',100*(10^mean(reshape(stats.dmse(iStride,iSlope,:),[],1))-1));

fprintf('\n\n');
fprintf('\\begin{tabular}{r r p{1cm} | p{1cm}p{1cm}p{1cm}p{1cm}} stride & $n_\\textrm{eff}$ & optimal mse & reduced dof & blind initial & expected mse \\\\ \\hline \\hline \n');
for iSlope = 1:length(slopes)
    
    fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
    for iStride=1:length(strides)
                fprintf('%d & %#.1f & %s &  %s  &  %s  &  %s  \\\\ \n', strides(iStride),mean(reshape(stat_structs{1}.neff_se(iStride,iSlope,:),[],1)), print_mse(stat_structs{1},iStride,iSlope), print_dmse(stat_structs{3},iStride,iSlope), print_dmse(stat_structs{2},iStride,iSlope), print_dmse(stat_structs{4},iStride,iSlope) )  ;

%         fprintf('%#.1f (%d) & %s &  %s (%#.1f) &  %s (%#.1f) &  %s (%#.1f) \\\\ \n',mean(reshape(stat_structs{1}.neff_se(iStride,iSlope,:),[],1)), strides(iStride), print_mse(stat_structs{1},iStride,iSlope), print_dmse(stat_structs{3},iStride,iSlope),mean(reshape(stat_structs{3}.neff_se(iStride,iSlope,:),[],1)), print_dmse(stat_structs{2},iStride,iSlope),mean(reshape(stat_structs{2}.neff_se(iStride,iSlope,:),[],1)), print_dmse(stat_structs{4},iStride,iSlope),mean(reshape(stat_structs{4}.neff_se(iStride,iSlope,:),[],1)) )  ;
%         fprintf('%#.1f (%d) & %s &  %s (%#.1f) &  %s (%#.1f) &  %s (%#.1f) &  %s (%#.1f) \\\\ \n',mean(reshape(stat_structs{1}.neff_se(iStride,iSlope,:),[],1)), strides(iStride), print_mse(stat_structs{1},iStride,iSlope), print_dmse(stat_structs{3},iStride,iSlope),mean(reshape(stat_structs{3}.neff_se(iStride,iSlope,:),[],1)), print_dmse(stat_structs{2},iStride,iSlope),mean(reshape(stat_structs{2}.neff_se(iStride,iSlope,:),[],1)), print_dmse(stat_structs{4},iStride,iSlope),mean(reshape(stat_structs{4}.neff_se(iStride,iSlope,:),[],1)), print_dmse(stat_structs{9},iStride,iSlope),mean(reshape(stat_structs{9}.neff_se(iStride,iSlope,:),[],1)) )  ;

    end
    
end
fprintf('\\end{tabular} \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MSE Comparison table of alternative methods
%
% Note that we redefine the optimal mse to be relative to the blind optimal
dmse = @(stats) log10(stats.mse./stat_structs{4}.mse);
for i=1:length(stat_structs)
    stat_structs{i}.dmse = dmse(stat_structs{i});
end

fprintf('\n\n');
fprintf('\\begin{tabular}{r p{1cm} | p{1cm}p{1cm}p{1cm}p{1cm}} $n_\\textrm{eff}$ (stride) & expected mse & ranged & cv & gcv & log-likelihood \\\\ \\hline \\hline \n');
for iSlope = 1:length(slopes)
    
    fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
    for iStride=1:length(strides)
        fprintf('%#.1f (%d) & %s &  %s &  %s &  %s & %s \\\\ \n',mean(reshape(stat_structs{4}.neff_se(iStride,iSlope,:),[],1)), strides(iStride), print_mse(stat_structs{4},iStride,iSlope), print_dmse(stat_structs{5},iStride,iSlope), print_dmse(stat_structs{6},iStride,iSlope), print_dmse(stat_structs{7},iStride,iSlope), print_dmse(stat_structs{8},iStride,iSlope) )  ;
    end
    
end
fprintf('\\end{tabular} \n');


return

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