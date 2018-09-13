scaleFactor = 1;
LoadFigureDefaults

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
shouldUseObservedSignalOnly = 0;

result_stride = 2.^(0:9)';
% result_stride = [2];

u_estimate_spectral = zeros(length(result_stride),1);
a_estimate_spectral = zeros(length(result_stride),1);
expectedDOF = zeros(length(result_stride),1);

lambda_blind_initial = zeros(length(result_stride),1);
rms_error_blind_initial = zeros(length(result_stride),1);
dof_out_blind_initial = zeros(length(result_stride),1);

lambda_blind_optimal = zeros(length(result_stride),1);
rms_error_blind_optimal = zeros(length(result_stride),1);
dof_out_blind_optimal = zeros(length(result_stride),1);

lambda_true_optimal = zeros(length(result_stride),1);
rms_error_true_optimal = zeros(length(result_stride),1);
dof_out_true_optimal = zeros(length(result_stride),1);
dof_var_out_true_optimal = zeros(length(result_stride),1);

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
    x_obs = data.x(indices) + data.epsilon_x(indices);
    y_obs = data.y(indices) + data.epsilon_y(indices);
    t_obs = data.t(indices);
    sigma = data.position_error;
    
    S = 3;
    T = 2;
    K = S+1;
    isIsotropic = 1;
        
    % Compute the expected tension before any fitting
    [lambda_blind_initial(i), expectedDOF(i)] = TensionSpline.ExpectedInitialTension(t_obs,[x_obs,y_obs],sigma,T,isIsotropic);
    
    % Fit
    spline_fit = TensionSpline(t_obs,[x_obs,y_obs],sigma, 'lambda', lambda_blind_initial(i), 'S', S, 'T', T);
    compute_rms_error = @() sqrt(mean(mean(  ([data.x(indicesAll),data.y(indicesAll)] - spline_fit(data.t(indicesAll))).^2,2 ),1));
    
    % record how well the initial fit did
    rms_error_blind_initial(i) = compute_rms_error();
    dof_out_blind_initial(i) = spline_fit.IsotropicDOF;
    
    % now optimize the tension by matching degrees-of-freedom
    lambda_blind_optimal(i) = TensionSpline.MatchedDOFSolution(spline_fit,expectedDOF(i));
    rms_error_blind_optimal(i) = compute_rms_error();
    dof_out_blind_optimal(i) = spline_fit.IsotropicDOF;
    
    % now optimize the tension using the true value
    lambda_true_optimal(i) = TensionSpline.OptimalTensionSolution(spline_fit,data.t(indicesAll),[data.x(indicesAll),data.y(indicesAll)]);
    rms_error_true_optimal(i) = compute_rms_error();
    dof_out_true_optimal(i) = spline_fit.IsotropicDOF;
    dof_var_out_true_optimal(i) = spline_fit.IsotropicVarianceDOF;
    
%     fprintf('S=%d, T=2, stride=%d, rms_error=%g, rms_error_blind_initial=%g, rms_error_blind_optimal=%g,\n', S, stride, rms_error_true_optimal(i), rms_error_blind_initial(i), rms_error_blind_optimal(i) );
    fprintf('%d & %#.3g m (%#.3g/%#.3g) &  %#.3g m (%#.3g) &  %#.3g m (%#.3g) \\\\ \n', result_stride(i), rms_error_true_optimal(i), dof_out_true_optimal(i), dof_var_out_true_optimal(i), rms_error_blind_optimal(i), dof_out_blind_optimal(i), rms_error_blind_initial(i), dof_out_blind_initial(i) )  ;
end

fprintf('\n\n');
fprintf('\\begin{tabular}{c | ccc} stride & optimal & iterated estimate & initial estimate \\\\ \\hline \\hline \n');
for i=1:length(result_stride)
    fprintf('%d & %#.3g m (%#.3g) &  %#.3g m (%#.3g) &  %#.3g m (%#.3g) \\\\ \n', result_stride(i), rms_error_true_optimal(i), dof_out_true_optimal(i), rms_error_blind_optimal(i), dof_out_blind_optimal(i), rms_error_blind_initial(i), dof_out_blind_initial(i) )  ;
end
fprintf('\\end{tabular} \n');