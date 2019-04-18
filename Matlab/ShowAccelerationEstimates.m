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
    slopes = -2;
    totalSlopes = length(slopes);
    
    S = 2;
    T = S;
    K = S+1;
    
%     result_stride = 2*[1;4;16;64];
%     totalStrides = length(result_stride);
%     
%     totalEnsembles = 25;
    
        result_stride = 2*[1];
    totalStrides = length(result_stride);
    
    totalEnsembles = 5;
    
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
                u_rms = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_obs,sqrt(noiseDistribution.variance),1);
                n_eff = SmoothingSpline.EffectiveSampleSizeFromUrms(u_rms, t_obs, sqrt(noiseDistribution.variance));
                
                % The idea here was to skip a bunch of unnecessary points
                % to trying to assess the noisy tension derivative... this
                % sort of works in that it no longer includes a bunch of
                % noise in the estimate, but it way under estimates the
                % power, especially in shallow slopes (where power is at
                % higher frequencies. The results is that we way over
                % tension.
                idx2 = 1:max(floor(n_eff/2),1):length(t_obs);
                idx2 = 1:length(t_obs);
                
                a_rms = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(t_obs(idx2),x_obs(idx2),sqrt(noiseDistribution.variance),T,1);

                
                u_estimate_spectral(iStride,iSlope,iEnsemble) = u_rms;
                a_estimate_spectral(iStride,iSlope,iEnsemble) = a_rms;

 
            end
            fprintf('\n');
        end
    end
    
%     save(filename, 'S', 'T', 'slopes', 'result_stride', 'u_rms_true', 'a_rms_true', 'u_rms_true_strided','a_rms_true_strided', 'u_estimate_spectral', 'a_estimate_spectral', 'expectedDOF','reduced_dof_knot_dof','mse_full_dof_true_optimal','mse_reduced_dof_true_optimal','mse_reduced_dof_blind_initial','mse_reduced_dof_blind_optimal','mse_reduced_dof_log_likelihood','mse_reduced_dof_log_likelihood_blind','dof_se_full_dof_true_optimal','dof_se_reduced_dof_true_optimal','dof_se_reduced_dof_blind_initial','dof_se_reduced_dof_blind_optimal', 'dof_se_reduced_dof_log_likelihood', 'dof_se_reduced_dof_log_likelihood_blind');
end
