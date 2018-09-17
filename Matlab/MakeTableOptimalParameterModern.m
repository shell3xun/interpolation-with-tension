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
shouldUseObservedSignalOnly = 1;

result_stride = 2.^(0:9)';
result_stride = [8];

u_estimate_spectral = zeros(length(result_stride),1);
a_estimate_spectral = zeros(length(result_stride),1);
expectedDOF = zeros(length(result_stride),1);

lambda_blind_initial = zeros(length(result_stride),1);
rms_error_blind_initial = zeros(length(result_stride),1);
dof_out_blind_initial = zeros(length(result_stride),1);

lambda_blind_optimal = zeros(length(result_stride),1);
rms_error_blind_optimal = zeros(length(result_stride),1);
dof_out_blind_optimal = zeros(length(result_stride),1);

lambda_blind_expectedMSE = zeros(length(result_stride),1);
rms_error_blind_expectedMSE = zeros(length(result_stride),1);
dof_out_blind_expectedMSE = zeros(length(result_stride),1);

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
    epsilon_x = data.epsilon_x(indices);
    epsilon_y = data.epsilon_y(indices);
    epsilon_x = data.position_error*randn(length(indices),1);
    epsilon_y = data.position_error*randn(length(indices),1);
    x_obs = data.x(indices) + epsilon_x;
    y_obs = data.y(indices) + epsilon_y;
    t_obs = data.t(indices);
    sigma = data.position_error;
    
    S = 3;
    T = 2;
    K = S+1;
    isIsotropic = 1;
        
    if 1 == 0
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
    else
        % Compute the expected tension before any fitting
        [lambda_blind_initial(i), expectedDOF(i)] = TensionSpline.ExpectedInitialTension(t_obs,x_obs,sigma,T,isIsotropic);
        
        % Fit
        spline_fit = TensionSpline(t_obs,x_obs,sigma, 'lambda', lambda_blind_initial(i), 'S', S, 'T', T);
        compute_rms_error = @() sqrt(mean(mean(  (data.x(indicesAll) - spline_fit(data.t(indicesAll))).^2,2 ),1));
        % record how well the initial fit did
        rms_error_blind_initial(i) = compute_rms_error();
        dof_out_blind_initial(i) = spline_fit.IsotropicDOF;
        
        % now optimize the tension by matching degrees-of-freedom
        lambda_blind_optimal(i) = TensionSpline.MatchedDOFSolution(spline_fit,expectedDOF(i));
        rms_error_blind_optimal(i) = compute_rms_error();
        dof_out_blind_optimal(i) = spline_fit.IsotropicDOF;
        
        lambda_blind_expectedMSE(i) = TensionSpline.MinimizeExpectedMeanSquareError(spline_fit);
        rms_error_blind_expectedMSE(i) = compute_rms_error();
        dof_out_blind_expectedMSE(i) = spline_fit.IsotropicDOF;
        
        % now optimize the tension using the true value
        lambda_true_optimal(i) = TensionSpline.MinimizeMeanSquareError(spline_fit,data.t(indices),data.x(indices));
        rms_error_true_optimal(i) = compute_rms_error();
        dof_out_true_optimal(i) = spline_fit.IsotropicDOF;
        dof_var_out_true_optimal(i) = spline_fit.IsotropicVarianceDOF;
        
%         % Compute RMS using observed signal
%         S = spline_fit.SmoothingMatrix;
%         SI = (S-eye(size(S)));
%         a = mean((SI*x_obs).^2);
%         b = mean(2*(SI*x_obs).*epsilon_x);
%         c = mean(epsilon_x.^2);
%         b1 = mean(2*(S*x_obs).*epsilon_x);
%         b2 = -mean(2*(eye(size(S))*x_obs).*epsilon_x);
%         fprintf('obs. def. rms: %#.3g m, (a,b(b1,b2),c) = (%#.3g m^2, %#.3g m^2=%#.3g m^2 + %#.3g m^2, %#.3g m^2), sum=%#.3g m\n', rms_error_true_optimal(i),a,b,b1,b2,c,sqrt(a+b+c));
%         
%         % Compute RMS using true signal
%         a = mean((SI*data.x(indices)).^2);
%         b = mean(2*(SI*data.x(indices)).*(S*epsilon_x));
%         c = mean((S*epsilon_x).^2);
%         c_alt = sigma*sigma*trace(S*S)/length(S);
%         fprintf('true. def. rms: %#.3g m, (a,b,c) = (%#.3g m^2, %#.3g m^2, %#.3g m^2), sum=%#.3g m, sum_wo_b=%#.3g m, sum_expectation=%#.3g m\n', rms_error_true_optimal(i),a,b,c,sqrt(a+b+c),sqrt(a+c),sqrt(a+c_alt));
% 
%         % Compute RMS using true signal
%         a = mean((SI*(x_obs - epsilon_x )).^2);
%         b = mean(2*(SI*data.x(indices)).*(S*epsilon_x));
%         c = mean((S*epsilon_x).^2);
%         fprintf('true. def. rms: %#.3g m, (a,b,c) = (%#.3g m^2, %#.3g m^2, %#.3g m^2), sum=%#.3g m\n', rms_error_true_optimal(i),a,b,c,sqrt(a+b+c));
%         
%         % Compute RMS using true signal
% %         SIx = SI*x_obs;
% %         SIe = SI*epsilon_x;
%         a1 = mean((SI*x_obs).^2);
%         a2 = mean((SI*epsilon_x).^2);
%         a3 = -2*mean((SI*x_obs).*(SI*epsilon_x));
%         c = mean((S*epsilon_x).^2);
%         fprintf('true. def. expanded. rms: %#.3g m, (a1,a2,a3,c) = (%#.3g m^2, %#.3g m^2, %#.3g m^2, %#.3g m^2), sum=%#.3g m\n', rms_error_true_optimal(i),a1,a2,a3,c,sqrt(a1+a2+a3+c));
%         
%         % now under expectation
%         a1 = mean((SI*x_obs).^2);
%         a2 = sigma*sigma*trace(SI*SI)/length(S);
%         a3 = -2*mean((SI*x_obs).*(SI*epsilon_x));
%         c = sigma*sigma*trace(S*S)/length(S);
%         fprintf('true. def. expectation. rms: %#.3g m, (a1,a2,a3,c) = (%#.3g m^2, %#.3g m^2, %#.3g m^2, %#.3g m^2), sum=%#.3g m\n', rms_error_true_optimal(i),a1,a2,a3,c,sqrt(a1+a2+a3+c));
% 
% %         a1 = mean((SI*x_obs).^2);
% %         a2 = 0;
% %         a3 = -2*mean((SI*x_obs).*(SI*data.epsilon_x(indices)));
% %         c = sigma*sigma*(1+2*trace(S*S-S)/length(S));
% %         fprintf('true. def. expanded. rms: %#.3g m, (a1,a2,a3,c) = (%#.3g m^2, %#.3g m^2, %#.3g m^2, %#.3g m^2), sum=%#.3g m\n', rms_error_true_optimal(i),a1,a2,a3,c,sqrt(a1+a2+a3+c));
% %         
% %         % Craven and Wahba unbiassed estimator
%         a = mean((SI*x_obs).^2);
%         b = -sigma*sigma*trace(SI*SI)/length(S);
%         c = sigma*sigma*trace(S*S)/length(S);
%         fprintf('craven-wahba rms: %#.3g m, (a,b,c) = (%#.3g m^2, %#.3g m^2, %#.3g m^2), sum=%#.3g m\n', rms_error_true_optimal(i),a,b,c,sqrt(a+b+c));
% 
%         % Craven and Wahba unbiassed estimator
%         a = mean((SI*x_obs).^2);
%         b = 2*sigma*sigma*trace(S)/length(S);
%         c = -sigma*sigma;
%         fprintf('craven-wahba rms: %#.3g m, (a,b,c) = (%#.3g m^2, %#.3g m^2, %#.3g m^2), sum=%#.3g m\n', rms_error_true_optimal(i),a,b,c,sqrt(a+b+c));

        % 
%         % Craven and Wahba unbiassed estimator, all terms
%         a = mean((SI*x_obs).^2);
%         b = mean(2*(SI*x_obs).*data.epsilon_x(indices));
%         c = 2*sigma*sigma*trace(S)/length(S);
%         d = sigma*sigma;
% %         d = 2*mean( (SI*x_obs).*data.epsilon_x(indices) );
%         fprintf('rms: %#.3g m, (a,b,c,d) = (%#.3g m^2, %#.3g m^2, %#.3g m^2, %#.3g m^2), sum=%#.3g m\n', rms_error_true_optimal(i),a,b,c,d,sqrt(a+b+c+d));
        
    end
    

    
    
    
%     fprintf('S=%d, T=2, stride=%d, rms_error=%g, rms_error_blind_initial=%g, rms_error_blind_optimal=%g,\n', S, stride, rms_error_true_optimal(i), rms_error_blind_initial(i), rms_error_blind_optimal(i) );
    fprintf('%d & %#.3g m (%#.3g/%#.3g) &  %#.3g m (%#.3g) &  %#.3g m (%#.3g) &  %#.3g m (%#.3g) \\\\ \n', result_stride(i), rms_error_true_optimal(i), dof_out_true_optimal(i), dof_var_out_true_optimal(i), rms_error_blind_expectedMSE(i),dof_out_blind_expectedMSE(i), rms_error_blind_optimal(i), dof_out_blind_optimal(i), rms_error_blind_initial(i), dof_out_blind_initial(i) )  ;
end

fprintf('\n\n');
fprintf('\\begin{tabular}{c | ccc} stride & optimal & iterated estimate & initial estimate \\\\ \\hline \\hline \n');
for i=1:length(result_stride)
    fprintf('%d & %#.3g m (%#.3g) &  %#.3g m (%#.3g) &  %#.3g m (%#.3g) \\\\ \n', result_stride(i), rms_error_true_optimal(i), dof_out_true_optimal(i), rms_error_blind_optimal(i), dof_out_blind_optimal(i), rms_error_blind_initial(i), dof_out_blind_initial(i) )  ;
end
fprintf('\\end{tabular} \n');