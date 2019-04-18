scaleFactor = 1;
LoadFigureDefaults

slopes = [-2; -3; -4];
totalSlopes = length(slopes);

result_stride = [1;2;4;8;16;32;64];

optimal_se_dof = zeros(length(result_stride),totalSlopes);
last_good_knot_dof = zeros(length(result_stride),totalSlopes);

shouldUseObservedSignalOnly = 1;

for iSlope = 1:length(slopes)
    
    slope = slopes(iSlope);
    
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
        
    for iStride=1:length(result_stride)
        stride = result_stride(iStride);
        
        % Reduce the total length in some cases
        if (stride < 10)
            shortenFactor = stride/10;
        else
            shortenFactor = 1;
        end
        
        indices = 1:stride:floor(shortenFactor*length(data.t));
        fprintf('dT = %d --- Using %d points with stride %d\n', stride*dt, length(indices), stride);
        if (shouldUseObservedSignalOnly == 1)
            indicesAll = indices;
        else
            indicesAll = 1:max(indices);
        end
        
        sigma = data.position_error;
        epsilon_x = data.epsilon_x(indices);
        epsilon_x = randn(length(indices),1);
        x_obs = data.x(indices) + epsilon_x;
        t_obs = data.t(indices);
        
        S = 2;
        T = 2;
        K = S+1;
        
        spline_x = SmoothingSpline(t_obs,x_obs,sigma, 'S', S, 'T', T, 'knot_dof', 1);
        optimal_lambda = SmoothingSpline.MinimizeMeanSquareError(spline_x,data.t(indicesAll),data.x(indicesAll));
        compute_rms_error = @() sqrt(mean(mean(  (data.x(indicesAll) - spline_x(data.t(indicesAll))).^2,2 ),1));
        
        optimal_rms_error = compute_rms_error();
        optimal_se_dof(iStride,iSlope) = spline_x.DOFFromVarianceOfTheMean;
        fprintf('DOF from SE: %.2f---', optimal_se_dof(iStride,iSlope));
        
        delta_relative_rms_error = 0;
        knot_dof = 1;
        fprintf('(knot dof, delta error): ');
        while delta_relative_rms_error < 0.01
            knot_dof = knot_dof+1;
            spline_x = SmoothingSpline(t_obs,x_obs,sigma, 'lambda', optimal_lambda, 'S', S, 'T', T, 'knot_dof', knot_dof);
            SmoothingSpline.MinimizeMeanSquareError(spline_x,data.t(indicesAll),data.x(indicesAll));
            compute_rms_error = @() sqrt(mean(mean(  (data.x(indicesAll) - spline_x(data.t(indicesAll))).^2,2 ),1));
            
            delta_relative_rms_error = (compute_rms_error()/optimal_rms_error) - 1;
            fprintf('(%d, %.4f) ', knot_dof, delta_relative_rms_error);
        end
        fprintf('\n');
        
        last_good_knot_dof(iStride,iSlope) = knot_dof-1;
    end
    
end

% The goal is to choose as high as knot_dof as possible based on the
% optimal_se_dof, without going past the last_god_knot_dof.
% This seems to be a pretty good estimate.
max(1,floor(0.66*optimal_se_dof))-last_good_knot_dof