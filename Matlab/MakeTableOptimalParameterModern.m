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

dt = data.t(2)-data.t(1);

decorrelation_u = find(Autocorrelation(diff(data.x)/dt,5000)<0.2,1,'first')*dt;
decorrelation_a = find(Autocorrelation(diff(diff(data.x))/(dt^2),1000)<0.2,1,'first')*dt;
decorrelation_j = find(Autocorrelation(diff(diff(diff(data.x)))/(dt^3),1000)<0.2,1,'first')*dt;

fprintf('Decorrelation time for (u,a,j) = (%.0f, %.0f, %.0f)\n', decorrelation_u, decorrelation_a, decorrelation_j);

% Do you want to assess the error using all the points from the signal
% (which makes sense for an interpolation based metric) or just points from
% the observed (decimated) signal only?
shouldUseObservedSignalOnly = 1;

S_range = 1:4;
result_stride = 2.^(3:9)';
result_stride = 8;

u_estimate_spectral = zeros(length(result_stride),1);
a_estimate_spectral = zeros(length(result_stride),1);
expectedDOF = zeros(length(result_stride),length(S_range),length(S_range));

lambda_blind_initial = zeros(length(result_stride),length(S_range),length(S_range));
rms_error_blind_initial = zeros(length(result_stride),length(S_range),length(S_range));
dof_out_blind_initial = zeros(length(result_stride),length(S_range),length(S_range));
Q_blind_initial = zeros(length(result_stride),length(S_range),length(S_range));

lambda_blind_expectedMSE = zeros(length(result_stride),length(S_range),length(S_range));
rms_error_blind_expectedMSE = zeros(length(result_stride),length(S_range),length(S_range));
dof_out_blind_expectedMSE = zeros(length(result_stride),length(S_range),length(S_range));
Q_blind_blind_expectedMSE = zeros(length(result_stride),length(S_range),length(S_range));

lambda_true_optimal = zeros(length(result_stride),length(S_range),length(S_range));
rms_error_true_optimal = zeros(length(result_stride),length(S_range),length(S_range));
dof_out_true_optimal = zeros(length(result_stride),length(S_range),length(S_range));
Q_blind_true_optimal = zeros(length(result_stride),length(S_range),length(S_range));

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
    epsilon_x = data.epsilon_x(indices);
    epsilon_y = data.epsilon_y(indices);
    epsilon_x = data.position_error*randn(length(indices),1);
    epsilon_y = data.position_error*randn(length(indices),1);
    x_obs = data.x(indices) + epsilon_x;
    y_obs = data.y(indices) + epsilon_y;
    t_obs = data.t(indices);
    sigma = data.position_error;
    
    for iS = 1:length(S_range)
        S = S_range(iS);
        for T = 1:S
            fprintf('\tS=%d, T=%d\n',S,T)
            K = S+1;
            isIsotropic = 1;
            
            
            % Compute the expected tension before any fitting
            [lambda_blind_initial(i,iS,T), expectedDOF(i,iS,T)] = TensionSpline.ExpectedInitialTension(t_obs,x_obs,sigma,T,isIsotropic);
            
            % Fit
            spline_fit = TensionSpline(t_obs,x_obs,sigma, 'lambda', lambda_blind_initial(i,iS,T), 'S', S, 'T', T);
            compute_rms_error = @() sqrt(mean(mean(  (data.x(indicesAll) - spline_fit(data.t(indicesAll))).^2,2 ),1));
            
            % record how well the initial fit did
            rms_error_blind_initial(i,iS,T) = compute_rms_error();
            dof_out_blind_initial(i,iS,T) = spline_fit.DOFFromVarianceOfTheMean;
            [rmse,norm] = TensionSpline.MeanSquareErrorAtAllOrders(spline_fit, data.t(indicesAll),data.x(indicesAll));
            Q_blind_initial(i,iS,T) = find(rmse./norm < 0.9,1,'last');
            
            fprintf('\t\tinitial (rms,lambda,dof,Q)=(%#.3g m, %#.3g, %#.3g, %#.3g)\n',rms_error_blind_initial(i,iS,T),lambda_blind_initial(i,iS,T),dof_out_blind_initial(i,iS,T),Q_blind_initial(i,iS,T));
            
            % now optimize the tension using the true value
            lambda_true_optimal(i,iS,T) = TensionSpline.MinimizeMeanSquareError(spline_fit,data.t(indicesAll),data.x(indicesAll));
            rms_error_true_optimal(i,iS,T) = compute_rms_error();
            dof_out_true_optimal(i,iS,T) = spline_fit.DOFFromVarianceOfTheMean;
            [rmse,norm] = TensionSpline.MeanSquareErrorAtAllOrders(spline_fit, data.t(indicesAll),data.x(indicesAll));
            Q_blind_true_optimal(i,iS,T) = find(rmse./norm < 0.9,1,'last');
            
            fprintf('\t\toptimal (rms,lambda,dof,Q)=(%#.3g m, %#.3g, %#.3g, %#.3g)\n',rms_error_true_optimal(i,iS,T),lambda_true_optimal(i,iS,T),dof_out_true_optimal(i,iS,T), Q_blind_true_optimal(i,iS,T));
            
%             lambda_blind_expectedMSE(i,iS,T) = TensionSpline.MinimizeExpectedMeanSquareError(spline_fit);
%             rms_error_blind_expectedMSE(i,iS,T) = compute_rms_error();
%             dof_out_blind_expectedMSE(i,iS,T) = spline_fit.DOFFromVarianceOfTheMean;
%             [rmse,norm] = TensionSpline.MeanSquareErrorAtAllOrders(spline_fit, data.t(indicesAll),data.x(indicesAll));
%             Q_blind_blind_expectedMSE(i,iS,T) = find(rmse./norm < 0.9,1,'last');
%             
%             fprintf('\t\texpected mse (rms,lambda,dof,Q)=(%#.3g m, %#.3g, %#.3g, %#.3g)\n',rms_error_blind_expectedMSE(i,iS,T),lambda_blind_expectedMSE(i,iS,T),dof_out_blind_expectedMSE(i,iS,T), Q_blind_blind_expectedMSE(i,iS,T));
            
        end
    end
    
    
    
%     fprintf('S=%d, T=2, stride=%d, rms_error=%g, rms_error_blind_initial=%g, rms_error_blind_optimal=%g,\n', S, stride, rms_error_true_optimal(i), rms_error_blind_initial(i), rms_error_blind_optimal(i) );
%     fprintf('%d & %#.3g m (%#.3g/%#.3g) &  %#.3g m (%#.3g) &  %#.3g m (%#.3g) &  %#.3g m (%#.3g) \\\\ \n', result_stride(i), rms_error_true_optimal(i), dof_out_true_optimal(i), dof_var_out_true_optimal(i), rms_error_blind_expectedMSE(i),dof_out_blind_expectedMSE(i), rms_error_blind_optimal(i), dof_out_blind_optimal(i), rms_error_blind_initial(i), dof_out_blind_initial(i) )  ;
end

SSER = spline_fit.SignalToStandardErrorRatio;

% [MSE,noise] = spline_fit.ExpectedMeanSquareErrorAtAllOrders();
% sqrt(MSE)./rmse
% sqrt(MSE./noise)

return

S = spline_fit.SmoothingMatrix;
sigma2 = sigma*sigma;
Diff = TensionSpline.FiniteDifferenceMatrixNoBoundary(K-1,t_obs,1);
DS = Diff*S;
A = (DS-Diff);
N = length(Diff);
MSE = ( sum( (A*data.x(indices)).^2) + sigma2*sum(sum(DS.^2)))/N
sqrt(MSE)/rmse(end)
MSE = (sum( (A*x_obs).^2) + sum((A*epsilon_x).^2) - 2*sum((A*x_obs).*(A*epsilon_x)) +  sigma2*sum(sum(DS.^2)))/N
sqrt(MSE)/rmse(end)

MSE = (sum( (A*x_obs).^2) + sigma2*sum(sum(A.^2)) - 2*sum((A*x_obs).*(A*epsilon_x)) +  sigma2*sum(sum(DS.^2)))/N
sqrt(MSE)/rmse(end)

MSE = (sum( (A*x_obs).^2) - sigma2*sum(sum(A.^2)) +  sigma2*sum(sum(DS.^2)))/N
sqrt(MSE)/rmse(end)

signal = spline_fit(t_obs);
dt = t_obs(2)-t_obs(1);
var_sig = var(signal);
var_noise = sigma*sigma*spline_fit.VarianceOfTheMean;
fprintf('var_sig: %.2g, var_noise: %.2g, ratio: %.2g\n',var_sig, var_noise, var_sig/var_noise);
a = [2;6;20;70;252];
for iDim = 1:S
    signal = diff(signal)/dt;
    var_sig = var(signal);
    var_noise = mean(diag(spline_fit.CovarianceMatrixForDerivative(iDim)));
%     var_noise = a(iDim)*sigma*sigma*spline_fit.VarianceOfTheMean/(dt^iDim);
    fprintf('var_sig: %.2g, var_noise: %.2g, ratio: %.2g\n',var_sig, var_noise, var_sig/var_noise);
end

fprintf('\n\n');
fprintf('\\begin{tabular}{c | ccc} stride & optimal & iterated estimate & initial estimate \\\\ \\hline \\hline \n');
for i=1:length(result_stride)
    fprintf('%d & %#.3g m (%#.3g) &  %#.3g m (%#.3g) &  %#.3g m (%#.3g) \\\\ \n', result_stride(i), rms_error_true_optimal(i), dof_out_true_optimal(i), rms_error_blind_optimal(i), dof_out_blind_optimal(i), rms_error_blind_initial(i), dof_out_blind_initial(i) )  ;
end
fprintf('\\end{tabular} \n');