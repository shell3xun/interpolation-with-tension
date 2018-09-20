scaleFactor = 1;
LoadFigureDefaults
addpath('./support');

slope = -2;

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

% Do you want to assess the error using all the points from the signal
% (which makes sense for an interpolation based metric) or just points from
% the observed (decimated) signal only?
shouldUseObservedSignalOnly = 1;

result_stride = 2.^(0:9)';
result_stride = [16];

expectedDOF = zeros(length(result_stride),1);

lambda_blind_initial = zeros(length(result_stride),1);
rms_error_blind_initial = zeros(length(result_stride),1);
dof_out_blind_initial = zeros(length(result_stride),1);

lambda_true_optimal = zeros(length(result_stride),1);
rms_error_true_optimal = zeros(length(result_stride),1);
dof_out_true_optimal = zeros(length(result_stride),1);

lambda_blind_expectedMSE = zeros(length(result_stride),1);
rms_error_blind_expectedMSE = zeros(length(result_stride),1);
dof_out_blind_expectedMSE = zeros(length(result_stride),1);

for i=1:length(result_stride)
    stride = result_stride(i);
    
    % Reduce the total length in some cases
    if (stride < 10)
        shortenFactor = stride/10;
    else
        shortenFactor = 1;
    end
       
    indices = 1:stride:floor(shortenFactor*length(data.t));
    fprintf('Using %d points with stride %d\n', length(indices), stride);
    if (shouldUseObservedSignalOnly == 1)
        indicesAll = indices;
    else
        indicesAll = 1:max(indices);
    end
    
    % Create some noise with a t-distribution
    nu = 4.5; sigma =  8.5;
    variance_of_the_noise = sigma*sigma*nu/(nu-2);
    w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
    epsilon_x = randt(sigma,nu,length(indices));
    epsilon_y = randt(sigma,nu,length(indices));
    
    x_obs = data.x(indices) + epsilon_x;
    y_obs = data.y(indices) + epsilon_y;
    t_obs = data.t(indices);
    
    S = 3;
    T = 2;
    K = S+1;
    isIsotropic = 1;
        
        % Compute the expected tension before any fitting
        [lambda_blind_initial(i), expectedDOF(i)] = TensionSpline.ExpectedInitialTension(t_obs,x_obs,sqrt(variance_of_the_noise),T,isIsotropic);
        
        % Fit
        spline_fit = TensionSpline(t_obs,x_obs,sqrt(variance_of_the_noise), 'lambda', lambda_blind_initial(i), 'S', S, 'T', T,'weightFunction',w);
        compute_rms_error = @() sqrt(mean(mean(  (data.x(indicesAll) - spline_fit(data.t(indicesAll))).^2,2 ),1));
        
        % record how well the initial fit did
        rms_error_blind_initial(i) = compute_rms_error();
        dof_out_blind_initial(i) = spline_fit.DOFFromExpectedMeanSquareError;
        
        fprintf('initial (rms,lambda,dof)=(%#.3g m, %#.3g, %#.3g)\n',rms_error_blind_initial(i),lambda_blind_initial(i),dof_out_blind_initial(i));
        
        % now optimize the tension using the true value
        lambda_true_optimal(i) = TensionSpline.MinimizeMeanSquareError(spline_fit,data.t(indices),data.x(indices));
        rms_error_true_optimal(i) = compute_rms_error();
        dof_out_true_optimal(i) = spline_fit.DOFFromExpectedMeanSquareError;
        
        fprintf('optimal (rms,lambda,dof)=(%#.3g m, %#.3g, %#.3g)\n',rms_error_true_optimal(i),lambda_true_optimal(i),dof_out_true_optimal(i));
        
        lambda_blind_expectedMSE(i) = TensionSpline.MinimizeExpectedMeanSquareError(spline_fit);
        rms_error_blind_expectedMSE(i) = compute_rms_error();
        dof_out_blind_expectedMSE(i) = spline_fit.DOFFromExpectedMeanSquareError;
        
        fprintf('expected mse (rms,lambda,dof)=(%#.3g m, %#.3g, %#.3g)\n',rms_error_blind_expectedMSE(i),lambda_blind_expectedMSE(i),dof_out_blind_expectedMSE(i));
end