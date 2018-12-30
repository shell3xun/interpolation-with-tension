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
    totalSlopes = length(slopes);
else
    slopes = [-2; -3; -4];
    totalSlopes = length(slopes);
    
    result_stride = [2;4;8;16;32;64;128];
    ensembles = 10*ones(size(result_stride));
    totalEnsembles = sum(ensembles);
    
    measurement_time = zeros(length(result_stride),totalSlopes);
    
    measured_u_rms = zeros(totalEnsembles,totalSlopes);
    measured_a_rms = zeros(totalEnsembles,totalSlopes);
    actual_u_rms = zeros(totalEnsembles,totalSlopes);
    actual_a_rms = zeros(totalEnsembles,totalSlopes);
    measured_u_rms_mean = zeros(length(result_stride),totalSlopes);
    measured_u_rms_median = zeros(length(result_stride),totalSlopes);
    measured_u_rms_std = zeros(length(result_stride),totalSlopes);
    measured_a_rms_mean = zeros(length(result_stride),totalSlopes);
    measured_a_rms_median = zeros(length(result_stride),totalSlopes);
    measured_a_rms_std = zeros(length(result_stride),totalSlopes);
    
    measured_lambda = zeros(totalEnsembles,totalSlopes);
    measured_lambda_mean = zeros(length(result_stride),totalSlopes);
    measured_lambda_median = zeros(length(result_stride),totalSlopes);
    measured_lambda_std = zeros(length(result_stride),totalSlopes);
    
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
        for iStride=1:length(result_stride)
            stride = result_stride(iStride);
            
            % Reduce the total length in some cases
            if (stride < 10)
                shortenFactor = stride/10;
            else
                shortenFactor = 1;
            end
            
            indices = 1:stride:floor(shortenFactor*length(data.t));
            fprintf('Using %d points with stride %d. Ensemble', length(indices), stride);
            indicesAll = indices;
            
            rms = @(z) sqrt( mean( z.^2 ) );
            
            u_rms = rms( diff( data.x(1:max(indices)) )/(data.t(2)-data.t(1)) );
            t_obs = data.t(indices);
            
            for j=1:ensembles(iStride)
                iteration = iteration + 1;
                fprintf('..%d',iteration);
                
                sigma = data.position_error;
                epsilon_x = sigma*randn(length(indices),1);
                x_obs = data.x(indices) + epsilon_x;
                

                S = 2;
                T = 2;
                distribution = NormalDistribution(sigma);
                spline_fit = TensionSpline(t_obs,x_obs,distribution, 'S', S, 'T', T);
                
                measured_u_rms(iteration,iSlope) = TensionSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_obs,sqrt(distribution.variance),1);
                measured_a_rms(iteration,iSlope) = TensionSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_obs,sqrt(distribution.variance),T);

                % Minimize to the true points---at the observation times only
                measured_lambda(iteration,iSlope) = spline_fit.minimizeMeanSquareError(data.t(indices),data.x(indices));
                
                measured_dof_mse(iteration,iSlope) = spline_fit.effectiveSampleSizeFromExpectedMeanSquareError;
                measured_dof_se(iteration,iSlope) = spline_fit.effectiveSampleSizeFromVarianceOfTheMean;
                measured_dof_var(iteration,iSlope) = spline_fit.effectiveSampleSizeFromSampleVariance;
                gamma(iteration,iSlope) = sigma/(u_rms*(t_obs(2)-t_obs(1)));
            end
            fprintf('..\n');
            
            dt = t_obs(2)-t_obs(1);
            actual_u_rms(iStride,iSlope) = rms( diff(data.x(indices))/dt );
            actual_a_rms(iStride,iSlope) = rms( diff(diff(data.x(indices)))/dt^2 );
            
            measurement_time(iStride,iSlope) = max(t_obs)-min(t_obs);
            
            cum_ensembles = cumsum(ensembles);
            startIndex = cum_ensembles(iStride) - ensembles(iStride) + 1;
            endIndex = cum_ensembles(iStride);
            measured_lambda_mean(iStride,iSlope) = mean( measured_lambda(startIndex:endIndex,iSlope) );
            measured_lambda_median(iStride,iSlope) = median( measured_lambda(startIndex:endIndex,iSlope) );
            measured_lambda_std(iStride,iSlope) = std( measured_lambda(startIndex:endIndex,iSlope) );
            
            measured_u_rms_mean(iStride,iSlope) = mean( measured_u_rms(startIndex:endIndex,iSlope) );
            measured_u_rms_median(iStride,iSlope) = median( measured_u_rms(startIndex:endIndex,iSlope) );
            measured_u_rms_std(iStride,iSlope) = std( measured_u_rms(startIndex:endIndex,iSlope) );
            
            measured_a_rms_mean(iStride,iSlope) = mean( measured_a_rms(startIndex:endIndex,iSlope) );
            measured_a_rms_median(iStride,iSlope) = median( measured_a_rms(startIndex:endIndex,iSlope) );
            measured_a_rms_std(iStride,iSlope) = std( measured_a_rms(startIndex:endIndex,iSlope) );
            
            measured_dof_mse_mean(iStride,iSlope) = mean( measured_dof_mse(startIndex:endIndex,iSlope) );
            measured_dof_mse_median(iStride,iSlope) = median( measured_dof_mse(startIndex:endIndex,iSlope) );
            measured_dof_mse_std(iStride,iSlope) = std( measured_dof_mse(startIndex:endIndex,iSlope) );
            
            measured_dof_se_mean(iStride,iSlope) = mean( measured_dof_se(startIndex:endIndex,iSlope) );
            measured_dof_se_median(iStride,iSlope) = median( measured_dof_se(startIndex:endIndex,iSlope) );
            measured_dof_se_std(iStride,iSlope) = std( measured_dof_se(startIndex:endIndex,iSlope) );
            
            measured_dof_var_mean(iStride,iSlope) = mean( measured_dof_var(startIndex:endIndex,iSlope) );
            measured_dof_var_median(iStride,iSlope) = median( measured_dof_var(startIndex:endIndex,iSlope) );
            measured_dof_var_std(iStride,iSlope) = std( measured_dof_var(startIndex:endIndex,iSlope) );
            
            gamma_mean(iStride,iSlope) = mean( gamma(startIndex:endIndex,iSlope) );
        end
    end
    
    save('DegreesOfFreedomEstimates.mat','result_stride', 'slopes', 'ensembles', 'measured_u_rms', 'measured_a_rms', 'measured_dof_mse','measured_dof_se','measured_dof_var','gamma','measured_dof_mse_mean','measured_dof_mse_median','measured_dof_mse_std','measured_dof_se_mean','measured_dof_se_median', 'measured_dof_se_std','measured_dof_var_mean','measured_dof_var_median','measured_dof_var_std','gamma_mean','measured_u_rms_mean','measured_u_rms_median','measured_u_rms_std','measured_a_rms_mean','measured_a_rms_median','measured_a_rms_std','measured_lambda_mean','measured_lambda_median','measured_lambda_std','actual_u_rms','actual_a_rms','measurement_time'); 
end

figure
for iSlope = 1:totalSlopes
    x = (gamma_mean(:,iSlope))./(measured_a_rms_mean(:,iSlope).^2 .* measurement_time(:,iSlope) );
    y = measured_lambda_mean(:,iSlope);
    scatter(x,y), hold on
end

return

% figure
% 
% subplot(3,1,1)
% for iSlope = 1:totalSlopes
%     scatter(gamma(:,iSlope),measured_dof_mse(:,iSlope)), hold on
% end
% xlabel('dof from \Gamma')
% ylabel('mse dof')
% 
% 
% subplot(3,1,2)
% for iSlope = 1:totalSlopes
%     scatter(gamma(:,iSlope),measured_dof_se(:,iSlope)), hold on
% end
% xlabel('dof from \Gamma')
% ylabel('se dof')
% 
% 
% subplot(3,1,3)
% for iSlope = 1:totalSlopes
%     scatter(gamma(:,iSlope),measured_dof_var(:,iSlope)), hold on
% end
% xlabel('dof from \Gamma')
% ylabel('var dof')


FigureSize = [50 50 figure_width_1col+7 150*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize,'Name','DOFvsGamma');
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

m = zeros(length(totalSlopes),1);
b = zeros(length(totalSlopes),1);
for iSlope = 1:totalSlopes
    plotHandles(iSlope) = scatter(gamma(:,iSlope),measured_dof_se(:,iSlope)); hold on
    
    x = gamma_mean(:,iSlope);
    y = measured_dof_se_mean(:,iSlope);
    sigma2 = (measured_dof_se_std(:,iSlope).*measured_dof_se_std(:,iSlope))./ensembles;
    
%     [m(iSlope),b(iSlope)] = LinearBestFitWithVariableError(x,y,sigma2);
%     plot(x, m(iSlope)*x + b(iSlope))
    
ax = gca;
ax.ColorOrderIndex = ax.ColorOrderIndex-1;
    [m(iSlope),b(iSlope)] = LinearBestFitWithVariableError(log(x),log(y),ones(size(y)));
    plot(x, exp(b(iSlope))*(x).^m(iSlope));
end
xlabel('\Gamma', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('n^{SE}_{eff}', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlog, ylog
leg = legend(plotHandles,'-2 slope','-3 slope','-4 slope');
leg.FontSize = figure_legend_size;
leg.FontName = figure_font;
leg.Location = 'northwest';
set( gca, 'FontSize', figure_axis_tick_size);

tightfig
print('-depsc2', '../figures/dofVsGamma.eps')

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

