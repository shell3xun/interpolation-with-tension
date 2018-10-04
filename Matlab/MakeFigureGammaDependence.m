% \Gamma = \sigma/(u*dt)
% The hypothesis is that d = 1+scale*\Gamma
%
% We find that the scale is larger for steep slopes, and smaller for
% shallower slopes. This makes intuitive sense, because the underlying
% function isn't changing as much more smoother functions (steeper slopes).
% So, your degrees of freedom increase more quickly for steeper slopes.

scaleFactor = 1;
LoadFigureDefaults

shouldLoadExistingResults = 1;

if shouldLoadExistingResults == 1
    load('DegreesOfFreedomEstimates.mat');
else
    slopes = [-2; -3; -4];
    totalSlopes = length(slopes);
    
    result_stride = [2;4;8;16;32;64;128];
    ensembles = 10*ones(size(result_stride));
    totalEnsembles = sum(ensembles);
    
    measured_dof_mse = zeros(totalEnsembles,totalSlopes);
    measured_dof_se = zeros(totalEnsembles,totalSlopes);
    measured_dof_var = zeros(totalEnsembles,totalSlopes);
    gamma = zeros(totalEnsembles,totalSlopes);
    
    measured_dof_mse_mean = zeros(length(result_stride),totalSlopes);
    measured_dof_mse_median = zeros(length(result_stride),totalSlopes);
    measured_dof_mse_std = zeros(length(result_stride),totalSlopes);
    measured_dof_se_mean = zeros(length(result_stride),totalSlopes);
    measured_dof_se_median = zeros(length(result_stride),totalSlopes);
    measured_dof_se_std = zeros(length(result_stride),totalSlopes);
    measured_dof_var_mean = zeros(length(result_stride),totalSlopes);
    measured_dof_var_median = zeros(length(result_stride),totalSlopes);
    measured_dof_var_std = zeros(length(result_stride),totalSlopes);
    gamma_mean = zeros(length(result_stride),totalSlopes);
    
    
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
                
                S = 2;
                T = 2;
                K = S+1;
                isIsotropic = 1;
                
                
                u_rms = sqrt(mean((diff( data.x(1:max(indices)) )/(data.t(2)-data.t(1))).^2));
                a_rms = sqrt(mean((diff(diff( data.x(indicesAll) ))/(data.t(2)-data.t(1)))).^2);
                spline_fit = TensionSpline(t_obs,x_obs,sigma, 'S', S, 'T', T);
                
                % Minimize to the observed points only
                TensionSpline.MinimizeMeanSquareError(spline_fit,data.t(indices),data.x(indices));
                
                measured_dof_mse(iteration,iSlope) = spline_fit.DOFFromExpectedMeanSquareError;
                measured_dof_se(iteration,iSlope) = spline_fit.DOFFromVarianceOfTheMean;
                measured_dof_var(iteration,iSlope) = spline_fit.DOFFromSampleVariance;
                gamma(iteration,iSlope) = sigma/(u_rms*(t_obs(2)-t_obs(1)));
            end
            
            cum_ensembles = cumsum(ensembles);
            startIndex = cum_ensembles(i) - ensembles(i) + 1;
            endIndex = cum_ensembles(i);
            measured_dof_mse_mean(i,iSlope) = mean( measured_dof_mse(startIndex:endIndex,iSlope) );
            measured_dof_mse_median(i,iSlope) = median( measured_dof_mse(startIndex:endIndex,iSlope) );
            measured_dof_mse_std(i,iSlope) = std( measured_dof_mse(startIndex:endIndex,iSlope) );
            measured_dof_se_mean(i,iSlope) = mean( measured_dof_se(startIndex:endIndex,iSlope) );
            measured_dof_se_median(i,iSlope) = median( measured_dof_se(startIndex:endIndex,iSlope) );
            measured_dof_se_std(i,iSlope) = std( measured_dof_se(startIndex:endIndex,iSlope) );
            measured_dof_var_mean(i,iSlope) = mean( measured_dof_var(startIndex:endIndex,iSlope) );
            measured_dof_var_median(i,iSlope) = median( measured_dof_var(startIndex:endIndex,iSlope) );
            measured_dof_var_std(i,iSlope) = std( measured_dof_var(startIndex:endIndex,iSlope) );
            gamma_mean(i,iSlope) = mean( gamma(startIndex:endIndex,iSlope) );
        end
    end
    
    save('DegreesOfFreedomEstimates.mat','result_stride', 'slopes', 'ensembles', 'measured_dof_mse','measured_dof_se','measured_dof_var','gamma','measured_dof_mse_mean','measured_dof_mse_median','measured_dof_mse_std','measured_dof_se_mean','measured_dof_se_median', 'measured_dof_se_std','measured_dof_var_mean','measured_dof_var_median','measured_dof_var_std','gamma_mean');
end

figure

subplot(3,1,1)
for iSlope = 1:totalSlopes
    scatter(gamma(:,iSlope),measured_dof_mse(:,iSlope)), hold on
end
xlabel('dof from \Gamma')
ylabel('mse dof')


subplot(3,1,2)
for iSlope = 1:totalSlopes
    scatter(gamma(:,iSlope),measured_dof_se(:,iSlope)), hold on
end
xlabel('dof from \Gamma')
ylabel('se dof')


subplot(3,1,3)
for iSlope = 1:totalSlopes
    scatter(gamma(:,iSlope),measured_dof_var(:,iSlope)), hold on
end
xlabel('dof from \Gamma')
ylabel('var dof')


figure
m = zeros(length(totalSlopes),1);
b = zeros(length(totalSlopes),1);
for iSlope = 1:totalSlopes
    scatter(gamma(:,iSlope),measured_dof_se(:,iSlope)), hold on
    
    x = gamma_mean(:,iSlope);
    y = measured_dof_se_mean(:,iSlope);
    sigma2 = (measured_dof_se_std(:,iSlope).*measured_dof_se_std(:,iSlope))./ensembles;
    
%     [m(iSlope),b(iSlope)] = LinearBestFitWithVariableError(x,y,sigma2);
%     plot(x, m(iSlope)*x + b(iSlope))
    
    [m(iSlope),b(iSlope)] = LinearBestFitWithVariableError(log(x),log(y),ones(size(y)));
    plot(x, exp(b(iSlope))*(x).^m(iSlope))
end
xlabel('dof from \Gamma')
ylabel('se dof')


return
% Best fit to a straight line, using the estimate of the standard error for
% sigma. From Press et al page 662
sigma2 = (measured_dof_mse_std.*measured_dof_mse_std)./ensembles;
x = d_gamma_mean;
y = measured_dof_mse_mean;

S = sum(1./sigma2);
Sx = sum(x./sigma2);
Sy = sum(y./sigma2);
Sxx = sum((x.*x)./sigma2);
Sxy = sum((x.*y)./sigma2);
Delta = S*Sxx - Sx*Sx;
a = (Sxx*Sy - Sx*Sxy)/Delta;
b = (S*Sxy - Sx*Sy)/Delta;

figure

subplot(3,1,1)
scatter(d_gamma,measured_dof_mse), hold on
plot(d_gamma_mean, b*d_gamma_mean + a)
xlabel('dof from \Gamma')
ylabel('mse dof')
title(sprintf('best fit slope: %f',b))



sigma2 = (measured_dof_se_std.*measured_dof_se_std)./ensembles;
x = d_gamma_mean;
y = measured_dof_se_mean;

S = sum(1./sigma2);
Sx = sum(x./sigma2);
Sy = sum(y./sigma2);
Sxx = sum((x.*x)./sigma2);
Sxy = sum((x.*y)./sigma2);
Delta = S*Sxx - Sx*Sx;
a = (Sxx*Sy - Sx*Sxy)/Delta;
b = (S*Sxy - Sx*Sy)/Delta;

subplot(3,1,2)
scatter(d_gamma,measured_dof_se), hold on
plot(d_gamma_mean, b*d_gamma_mean + a)
xlabel('dof from \Gamma')
ylabel('se dof')
title(sprintf('best fit slope: %f',b))



sigma2 = (measured_dof_var_std.*measured_dof_var_std)./ensembles;
x = d_gamma_mean;
y = measured_dof_var_mean;

S = sum(1./sigma2);
Sx = sum(x./sigma2);
Sy = sum(y./sigma2);
Sxx = sum((x.*x)./sigma2);
Sxy = sum((x.*y)./sigma2);
Delta = S*Sxx - Sx*Sx;
a = (Sxx*Sy - Sx*Sxy)/Delta;
b = (S*Sxy - Sx*Sy)/Delta;

subplot(3,1,3)
scatter(d_gamma,measured_dof_var), hold on
plot(d_gamma_mean, b*d_gamma_mean + a)
xlabel('dof from \Gamma')
ylabel('var dof')
title(sprintf('best fit slope: %f',b))

%%%%%%%%%%%%%%%%%%%%%%%%%
figure

sigma2 = (measured_dof_se_std.*measured_dof_se_std)./ensembles;
x = d_gamma_mean.^(2/3);
y = measured_dof_se_mean;

S = sum(1./sigma2);
Sx = sum(x./sigma2);
Sy = sum(y./sigma2);
Sxx = sum((x.*x)./sigma2);
Sxy = sum((x.*y)./sigma2);
Delta = S*Sxx - Sx*Sx;
a = (Sxx*Sy - Sx*Sxy)/Delta;
b = (S*Sxy - Sx*Sy)/Delta;

scatter(d_gamma,measured_dof_se), hold on
plot(d_gamma_mean, b*(d_gamma_mean.^(2/3)) + a)
xlabel('\Gamma')
ylabel('se dof')
title(sprintf('best fit slope: %.2f + %.1f\Gamma^{(2/3)}',a,b))


% D = 2;
% [p,~,mu]=polyfit(x,y,D);
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

