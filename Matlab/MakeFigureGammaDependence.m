% \Gamma = \sigma/(u*dt)
% The hypothesis is that d = 1+scale*\Gamma
%
% We find that the scale is larger for steep slopes, and smaller for
% shallower slopes. This makes intuitive sense, because the underlying
% function isn't changing as much more smoother functions (steeper slopes).
% So, your degrees of freedom increase more quickly for steeper slopes.

scaleFactor = 1;
LoadFigureDefaults

shouldLoadExistingResults = 0;

rms = @(z) sqrt( mean( z.^2 ) );

if shouldLoadExistingResults == 1
    load('DegreesOfFreedomEstimates.mat');
    totalSlopes = length(slopes);
else
    slopes = [-2; -3; -4];
    totalSlopes = length(slopes);
    
    strides = [2;4;8;16;32;64;128];
    totalStrides = length(strides);
    
    totalEnsembles = 10; 
    
%     dt = zeros(totalStrides,1);
%     measurement_time = zeros(totalStrides,totalSlopes); 
%     actual_u_rms = zeros(totalStrides,totalSlopes);
%     actual_a_rms = zeros(totalStrides,totalSlopes);   
%     measured_u_rms = zeros(totalStrides,totalSlopes,totalEnsembles);
%     measured_a_rms = zeros(totalStrides,totalSlopes,totalEnsembles);
%     measured_lambda = zeros(totalStrides,totalSlopes,totalEnsembles);
%     measured_dof_mse = zeros(totalStrides,totalSlopes,totalEnsembles);
%     measured_dof_se = zeros(totalStrides,totalSlopes,totalEnsembles);
%     measured_dof_var = zeros(totalStrides,totalSlopes,totalEnsembles);
    
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
        
        for iStride=1:length(strides)
            stride = strides(iStride);
            
            % Reduce the total length in some cases
            if (stride < 10)
                shortenFactor = stride/10;
            else
                shortenFactor = 1;
            end
            
            indices = 1:stride:floor(shortenFactor*length(data.t));
            fprintf('Using %d points with stride %d. Ensemble', length(indices), stride);
            
            u_rms = rms( diff( data.x(1:max(indices)) )/(data.t(2)-data.t(1)) );
            sigma = data.position_error;
            t_obs = data.t(indices);
            
%             measurement_time(iStride,iSlope) = max(t_obs)-min(t_obs);
%             dt(iStride) = t_obs(2)-t_obs(1); 
%             actual_u_rms(iStride,iSlope) = rms( diff(data.x(indices))/dt(iStride) );
%             actual_a_rms(iStride,iSlope) = rms( diff(diff(data.x(indices)))/dt(iStride)^2 );
            
            for iEnsemble=1:totalEnsembles
                fprintf('..%d',iEnsemble);
                   
                epsilon_x = sigma*randn(length(indices),1);
                x_obs = data.x(indices) + epsilon_x;
                
                S = 2;
                T = 2;
                distribution = NormalDistribution(sigma);
%                 spline_fit = TensionSpline(t_obs,x_obs,distribution, 'S', S, 'T', T);
                
                measured_u_rms(iStride,iSlope,iEnsemble) = TensionSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_obs,sqrt(distribution.variance),1);
                measured_a_rms(iStride,iSlope,iEnsemble) = TensionSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_obs,sqrt(distribution.variance),T);
continue
                % Minimize to the true points---at the observation times only
                measured_lambda(iStride,iSlope,iEnsemble) = spline_fit.minimizeMeanSquareError(data.t(indices),data.x(indices));
                
                measured_dof_mse(iStride,iSlope,iEnsemble) = spline_fit.effectiveSampleSizeFromExpectedMeanSquareError;
                measured_dof_se(iStride,iSlope,iEnsemble) = spline_fit.effectiveSampleSizeFromVarianceOfTheMean;
                measured_dof_var(iStride,iSlope,iEnsemble) = spline_fit.effectiveSampleSizeFromSampleVariance;
            end
            fprintf('..\n');
        end
    end
    
    save('DegreesOfFreedomEstimates.mat','sigma','strides', 'slopes', 'totalEnsembles','dt','measurement_time','actual_u_rms','actual_a_rms','measured_u_rms','measured_a_rms','measured_lambda','measured_dof_mse','measured_dof_se','measured_dof_var'); 
end

measured_u_rms_mean = mean(measured_u_rms,3);
measured_u_rms_median = median(measured_u_rms,3);
measured_u_rms_std = std(measured_u_rms,0,3);

measured_a_rms_mean = mean(measured_a_rms,3);
measured_a_rms_median = median(measured_a_rms,3);
measured_a_rms_std = std(measured_a_rms,0,3);

measured_lambda_mean = mean(measured_lambda,3);
measured_lambda_median = median(measured_lambda,3);
measured_lambda_std = std(measured_lambda,0,3);

% measured_dof_mse_mean = zeros(length(strides),totalSlopes);
% measured_dof_mse_median = zeros(length(strides),totalSlopes);
% measured_dof_mse_std = zeros(length(strides),totalSlopes);
% 
% measured_dof_se_mean = zeros(length(strides),totalSlopes);
% measured_dof_se_median = zeros(length(strides),totalSlopes);
% measured_dof_se_std = zeros(length(strides),totalSlopes);
% 
% measured_dof_var_mean = zeros(length(strides),totalSlopes);
% measured_dof_var_median = zeros(length(strides),totalSlopes);
% measured_dof_var_std = zeros(length(strides),totalSlopes);

gamma_mean = mean(sigma./(measured_u_rms.*dt),3);

% lambda = (1-1/n)/(xt)
% xt*lambda = 1-1/n
% 1/n = 1-xt*lambda
figure
for iSlope = 1:totalSlopes
%     x = (gamma_mean(:,iSlope))./(measured_a_rms_mean(:,iSlope).^2 .* measurement_time(:,iSlope) );
%     y = measured_lambda_mean(:,iSlope);
    
    x = gamma_mean(:,iSlope);
    y = (measured_a_rms_mean(:,iSlope).^2 ).*measured_lambda_mean(:,iSlope);
    
    scatter(x,y), hold on
end

% biggest problem is estimating the rms value from a noisy derivative.
% using sleptap helps quite a bit, but the smallest strides still fail
figure
x = reshape(sigma./(measured_u_rms.*dt),[],1);
y = reshape((measured_a_rms.^2 ).*measured_lambda,[],1);
scatter(x,y)

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

