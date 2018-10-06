scaleFactor = 1;
LoadFigureDefaults

shouldLoadExistingTable = 1;
if shouldLoadExistingTable == 1
    load('MSEComparisonTable.mat');
else
    slopes = [-2; -3; -4];
    totalSlopes = length(slopes);
    
    S = 2;
    T = 2;
    K = S+1;
    
    result_stride = 2*[1;4;16;64];
    totalStrides = length(result_stride);
    
    totalEnsembles = 15;
    
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
    
    u_estimate_spectral = zeros(totalStrides,totalSlopes, totalEnsembles);
    a_estimate_spectral = zeros(totalStrides,totalSlopes, totalEnsembles);
    expectedDOF = zeros(totalStrides,totalSlopes, totalEnsembles);
    
    reduced_dof_knot_dof = zeros(totalStrides, totalSlopes, totalEnsembles);
    
    mse_full_dof_true_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_true_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_blind_initial = zeros(totalStrides, totalSlopes, totalEnsembles);
    mse_reduced_dof_blind_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    
    dof_se_full_dof_true_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_true_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_blind_initial = zeros(totalStrides, totalSlopes, totalEnsembles);
    dof_se_reduced_dof_blind_optimal = zeros(totalStrides, totalSlopes, totalEnsembles);
    
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
            sigma = data.position_error;
            
            for iEnsemble = 1:totalEnsembles
                fprintf('..%d',iEnsemble);
                
                epsilon_x = sigma*randn(length(indices),1);
                x_obs = data.x(indices) + epsilon_x;
                t_obs = data.t(indices);
                
                % Compute the expected tension before any fitting
                [lambda, expectedDOF(iStride,iSlope,iEnsemble), u_estimate_spectral(iStride,iSlope,iEnsemble), a_estimate_spectral(iStride,iSlope,iEnsemble)] = TensionSpline.ExpectedInitialTension(t_obs,x_obs,sigma,T);
                
                % Fit using 1 dof for each data point
                spline_x = TensionSpline(t_obs,x_obs,sigma, 'lambda', lambda, 'S', S, 'T', T, 'knot_dof', 1);
                TensionSpline.MinimizeMeanSquareError(spline_x,data.t(indicesAll),data.x(indicesAll));
                compute_ms_error = @() (mean(mean(  (data.x(indicesAll) - spline_x(data.t(indicesAll))).^2,2 ),1));
                
                mse_full_dof_true_optimal(iStride,iSlope,iEnsemble) = compute_ms_error();
                dof_se_full_dof_true_optimal(iStride,iSlope,iEnsemble) = spline_x.DOFFromVarianceOfTheMean;
                
                % Fit using the automatic knot dof algorithm
                spline_x = TensionSpline(t_obs,x_obs,sigma, 'lambda', lambda, 'S', S, 'T', T, 'knot_dof', 'auto');
                compute_ms_error = @() (mean(mean(  (data.x(indicesAll) - spline_x(data.t(indicesAll))).^2,2 ),1));
                reduced_dof_knot_dof(iStride,iSlope,iEnsemble) = spline_x.knot_dof;
                
                % record the initial fit
                mse_reduced_dof_blind_initial(iStride,iSlope,iEnsemble) = compute_ms_error();
                dof_se_reduced_dof_blind_initial(iStride,iSlope,iEnsemble) = spline_x.DOFFromVarianceOfTheMean;
                
                % record the blind optimal fit
                TensionSpline.MinimizeExpectedMeanSquareError(spline_x);
                mse_reduced_dof_blind_optimal(iStride,iSlope,iEnsemble) = compute_ms_error();
                dof_se_reduced_dof_blind_optimal(iStride,iSlope,iEnsemble) = spline_x.DOFFromVarianceOfTheMean;
                
                % record the true optimal fit
                TensionSpline.MinimizeMeanSquareError(spline_x,data.t(indicesAll),data.x(indicesAll));
                mse_reduced_dof_true_optimal(iStride,iSlope,iEnsemble) = compute_ms_error();
                dof_se_reduced_dof_true_optimal(iStride,iSlope,iEnsemble) = spline_x.DOFFromVarianceOfTheMean;
            end
            fprintf('\n');
        end
    end
    
    
    save('MSEComparisonTable.mat', 'S', 'T', 'slopes', 'result_stride', 'u_rms_true', 'a_rms_true', 'u_estimate_spectral', 'a_estimate_spectral', 'expectedDOF','reduced_dof_knot_dof','mse_full_dof_true_optimal','mse_reduced_dof_true_optimal','mse_reduced_dof_blind_initial','mse_reduced_dof_blind_optimal','dof_se_full_dof_true_optimal','dof_se_reduced_dof_true_optimal','dof_se_reduced_dof_blind_initial','dof_se_reduced_dof_blind_optimal');
end

reduced_dof_knot_dof_mean = mean(reduced_dof_knot_dof,3);

mse_full_dof_true_optimal_mean = mean(mse_full_dof_true_optimal,3);
mse_full_dof_true_optimal_std = std(mse_full_dof_true_optimal,0,3);

dmse_reduced_dof_true_optimal_mean = mean(mse_reduced_dof_true_optimal-mse_full_dof_true_optimal,3);
dmse_reduce_dof_true_optimal_std = std(mse_reduced_dof_true_optimal-mse_full_dof_true_optimal,0,3);

dmse_reduced_dof_blind_initial_mean = mean(mse_reduced_dof_blind_initial-mse_full_dof_true_optimal,3);
dmse_reduced_dof_blind_initial_std = std(mse_reduced_dof_blind_initial-mse_full_dof_true_optimal,0,3);

dmse_reduced_dof_blind_optimal_mean = mean(mse_reduced_dof_blind_optimal-mse_full_dof_true_optimal,3);
dmse_reduced_dof_blind_optimal_std = std(mse_reduced_dof_blind_optimal-mse_full_dof_true_optimal,0,3);

for iSlope = 1:length(slopes)
    fprintf('\n\n');
    fprintf('\\begin{tabular}{c | cccc} stride & optimal full dof & optimal reduced dof & initial estimate & blind optimal \\\\ \\hline \\hline \n');
    for iStride=1:length(result_stride)
%         fprintf('%d & %#.3g m^2 (%#.3g) &  %#.3g m^2 (%#.3g) &  %#.3g m^2 (%#.3g) &  %#.3g m^2 (%#.3g) \\\\ \n', result_stride(iStride), mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_full_dof_true_optimal(iStride,iSlope), dmse_reduced_dof_true_optimal_mean(iStride,iSlope), dof_se_reduced_dof_true_optimal(iStride,iSlope), dmse_reduced_dof_blind_initial_mean(iStride,iSlope), dof_se_reduced_dof_blind_initial(iStride,iSlope), dmse_reduced_dof_blind_optimal_mean(iStride,iSlope), dof_se_reduced_dof_blind_optimal(iStride,iSlope) )  ;
        fprintf('%d & %#.3g m$^2$ (%#.3g) &  %+.1f\\%% (%#.3g) &  %+.1f\\%% (%#.3g) &  %+.1f\\%% (%#.3g) \\\\ \n', result_stride(iStride), mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_full_dof_true_optimal(iStride,iSlope), 100*dmse_reduced_dof_true_optimal_mean(iStride,iSlope)./mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_reduced_dof_true_optimal(iStride,iSlope), 100*dmse_reduced_dof_blind_initial_mean(iStride,iSlope)./mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_reduced_dof_blind_initial(iStride,iSlope), 100*dmse_reduced_dof_blind_optimal_mean(iStride,iSlope)./mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_reduced_dof_blind_optimal(iStride,iSlope) )  ;
    end
    fprintf('\\end{tabular} \n');
end