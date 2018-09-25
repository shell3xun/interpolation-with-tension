scaleFactor = 1;
LoadFigureDefaults
addpath('./support');

slope = -4;

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

decorrelation_u = find(Autocorrelation(diff(data.x)/dt,5000)<0.2,1,'first')*dt;
decorrelation_a = find(Autocorrelation(diff(diff(data.x))/(dt^2),1000)<0.2,1,'first')*dt;
decorrelation_j = find(Autocorrelation(diff(diff(diff(data.x)))/(dt^3),1000)<0.2,1,'first')*dt;

fprintf('Decorrelation time for (u,a,j) = (%.0f, %.0f, %.0f)\n', decorrelation_u, decorrelation_a, decorrelation_j);

% Do you want to assess the error using all the points from the signal
% (which makes sense for an interpolation based metric) or just points from
% the observed (decimated) signal only?
shouldUseObservedSignalOnly = 0;

S_range = 1:4;
result_stride = 2.^(0:9)';
result_stride = [16];

expectedDOF = zeros(length(result_stride),length(S_range),length(S_range));

lambda_blind_initial = zeros(length(result_stride),length(S_range),length(S_range));
rms_error_blind_initial = zeros(length(result_stride),length(S_range),length(S_range));
dof_out_blind_initial = zeros(length(result_stride),length(S_range),length(S_range));

lambda_true_optimal = zeros(length(result_stride),length(S_range),length(S_range));
rms_error_true_optimal = zeros(length(result_stride),length(S_range),length(S_range));
dof_out_true_optimal = zeros(length(result_stride),length(S_range),length(S_range));

lambda_blind_expectedMSE = zeros(length(result_stride),length(S_range),length(S_range));
rms_error_blind_expectedMSE = zeros(length(result_stride),length(S_range),length(S_range));
dof_out_blind_expectedMSE = zeros(length(result_stride),length(S_range),length(S_range));

for i=1:length(result_stride)
    stride = result_stride(i);
    
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
    
    % Create some noise with a t-distribution
    nu = 4.5; sigma =  8.5;
    variance_of_the_noise = sigma*sigma*nu/(nu-2);
    w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
    epsilon_x = randt(sigma,nu,length(indices));
    epsilon_y = randt(sigma,nu,length(indices));
    
    x_obs = data.x(indices) + epsilon_x;
    y_obs = data.y(indices) + epsilon_y;
    t_obs = data.t(indices);
    
    for iS = 1:length(S_range)
        S = S_range(iS);
        for T = 1:S
            fprintf('\tS=%d, T=%d\n',S,T)
            K = S+1;
            isIsotropic = 1;
            
            % Compute the expected tension before any fitting
            [lambda_blind_initial(i,iS,T), expectedDOF(i,iS,T)] = TensionSpline.ExpectedInitialTension(t_obs,x_obs,sqrt(variance_of_the_noise),T,isIsotropic);
            
            % Fit
            spline_fit = TensionSpline(t_obs,x_obs,sqrt(variance_of_the_noise), 'lambda', lambda_blind_initial(i,iS,T), 'S', S, 'T', T,'weightFunction',w);
            compute_rms_error = @() sqrt(mean(mean(  (data.x(indicesAll) - spline_fit(data.t(indicesAll))).^2,2 ),1));
            
            % record how well the initial fit did
            rms_error_blind_initial(i,iS,T) = compute_rms_error();
            dof_out_blind_initial(i,iS,T) = spline_fit.DOFFromVarianceOfTheMean;
            
            fprintf('\t\tinitial (rms,lambda,dof)=(%#.3g m, %#.3g, %#.3g)\n',rms_error_blind_initial(i,iS,T),lambda_blind_initial(i,iS,T),dof_out_blind_initial(i,iS,T));
            
            % now optimize the tension using the true value
            lambda_true_optimal(i,iS,T) = TensionSpline.MinimizeMeanSquareError(spline_fit,data.t(indicesAll),data.x(indicesAll));
            rms_error_true_optimal(i,iS,T) = compute_rms_error();
            dof_out_true_optimal(i,iS,T) = spline_fit.DOFFromVarianceOfTheMean;
            
            fprintf('\t\toptimal (rms,lambda,dof)=(%#.3g m, %#.3g, %#.3g)\n',rms_error_true_optimal(i,iS,T),lambda_true_optimal(i,iS,T),dof_out_true_optimal(i,iS,T));
            
            lambda_blind_expectedMSE(i,iS,T) = TensionSpline.MinimizeExpectedMeanSquareError(spline_fit);
            rms_error_blind_expectedMSE(i,iS,T) = compute_rms_error();
            dof_out_blind_expectedMSE(i,iS,T) = spline_fit.DOFFromVarianceOfTheMean;
            
            fprintf('\t\texpected mse (rms,lambda,dof)=(%#.3g m, %#.3g, %#.3g)\n',rms_error_blind_expectedMSE(i,iS,T),lambda_blind_expectedMSE(i,iS,T),dof_out_blind_expectedMSE(i,iS,T));
        end
    end
end