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

result_stride = [4;8;16;32;64;128];
ensembles = 5*ones(size(result_stride));
totalEnsembles = sum(ensembles);

measured_dof = zeros(totalEnsembles,1);
d_gamma = zeros(totalEnsembles,1);

measured_dof_mean = zeros(length(result_stride),1);
measured_dof_median = zeros(length(result_stride),1);
measured_dof_std = zeros(length(result_stride),1);
d_gamma_mean = zeros(length(result_stride),1);

iteration = 0;
for i=1:length(result_stride)
    for j=1:ensembles(i)
        iteration = iteration + 1;
        
        stride = result_stride(i);
        
        % Reduce the total length in some cases
        if (stride < 10)
            shortenFactor = stride/10;
        else
            shortenFactor = 1;
        end
        
        indices = 1:stride:floor(shortenFactor*length(data.t));
        fprintf('Using %d points with stride %d\n', length(indices), stride);
        indicesAll = indices;
        
        
        epsilon_x = data.position_error*randn(length(indices),1);
        x_obs = data.x(indices) + epsilon_x;
        t_obs = data.t(indices);
        sigma = data.position_error;
        
        S = 3;
        T = 2;
        K = S+1;
        isIsotropic = 1;
        
        
        u_rms = sqrt(mean((diff( data.x(1:max(indices)) )/(data.t(2)-data.t(1))).^2));
        a_rms = sqrt(mean(diff(diff( data.x(indicesAll) ))/(data.t(2)-data.t(1))^2));
        spline_fit = TensionSpline(t_obs,x_obs,sigma, 'S', S, 'T', T);
        
        % Minimize to the observed points only
        TensionSpline.MinimizeMeanSquareError(spline_fit,data.t(indices),data.x(indices));
        
        measured_dof(iteration) = spline_fit.ExpectedMeanSquareErrorDOF;
        d_gamma(iteration) = 1 + sigma/(u_rms*(t_obs(2)-t_obs(1)));
    end
    
    cum_ensembles = cumsum(ensembles);
    startIndex = cum_ensembles(i) - ensembles(i) + 1;
    endIndex = cum_ensembles(i);
    measured_dof_mean(i) = mean( measured_dof(startIndex:endIndex) );
    measured_dof_median(i) = median( measured_dof(startIndex:endIndex) );
    measured_dof_std(i) = std( measured_dof(startIndex:endIndex) );
    d_gamma_mean(i) = mean( d_gamma(startIndex:endIndex) );
end


% Best fit to a straight line, using the estimate of the standard error for
% sigma. From Press et al page 662
sigma2 = (measured_dof_std.*measured_dof_std)./ensembles;
x = d_gamma_mean;
y = measured_dof_mean;

S = sum(1./sigma2);
Sx = sum(x./sigma2);
Sy = sum(y./sigma2);
Sxx = sum((x.*x)./sigma2);
Sxy = sum((x.*y)./sigma2);
Delta = S*Sxx - Sx*Sx;
a = (Sxx*Sy - Sx*Sxy)/Delta;
b = (S*Sxy - Sx*Sy)/Delta;

figure
scatter(d_gamma,measured_dof), hold on
plot(d_gamma_mean, b*d_gamma_mean + a)
xlabel('dof from \Gamma')
ylabel('measured dof')
title(sprintf('best fit slope: %f',b))


% D = 1;
% [p,~,mu]=polyfit(d_gamma,measured_dof,D);
% slope = factorial(D)*p(1)/mu(2)^D;
% 
% [p,~,mu]=polyfit(d_gamma_mean,measured_dof_median,D);
% slope = factorial(D)*p(1)/mu(2)^D
% intercept = p(2)-p(1)*mu(1)/mu(2)
% 
% figure
% scatter(d_gamma_mean,measured_dof_median), hold on
% plot(d_gamma_mean, slope*d_gamma_mean + intercept)
% plot(d_gamma_mean, b*d_gamma_mean + a)

