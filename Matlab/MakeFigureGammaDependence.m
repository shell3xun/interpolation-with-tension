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

addpath('support')

rms = @(z) sqrt( mean( z.^2 ) );

filename = 'DegreesOfFreedomEstimates.mat';

if exist(filename,'file')
    load(filename);   
else
    slopes = [-2; -3; -4];
    strides = (2.^(0:5)).';
    
%     slopes = [-4];
%     strides = [1;128];
    
    totalSlopes = length(slopes);
    totalStrides = length(strides);
    totalEnsembles = 50;
    
    S = 3;
    T = S;
    
    varnames = {'S', 'T', 'slopes', 'strides', 'totalEnsembles','sigma'};
    
    % matern signal parameters
    sigma_u = 0.20;
    base_dt = 3*60; % chosen as the smallest interval considered, because anything shorter than this looks non-stationary... like a local polynomial fit is needed.
    t_damp = 30*60;
    n = 250;

    nothing = nan(totalStrides,totalSlopes,totalEnsembles);
    actual_u_rms = nothing; varnames{end+1} = 'actual_u_rms';
    actual_a_rms = nothing; varnames{end+1} = 'actual_a_rms';
    actual_a_mean = nothing; varnames{end+1} = 'actual_a_mean';
    actual_a_std = nothing; varnames{end+1} = 'actual_a_std';
    actual_ad_error = nothing; varnames{end+1} = 'actual_ad_error';
    measurement_time = nan(totalStrides,totalSlopes); varnames{end+1} = 'measurement_time';
    measured_u_rms = nothing; varnames{end+1} = 'measured_u_rms';
    measured_a_rms = nothing; varnames{end+1} = 'measured_a_rms';
    measured_a_mean = nothing; varnames{end+1} = 'measured_a_mean';
    measured_a_std = nothing; varnames{end+1} = 'measured_a_std';
    measured_ad_error = nothing; varnames{end+1} = 'measured_ad_error';
    measured_u_rms_filtered = nothing; varnames{end+1} = 'measured_u_rms_filtered';
    measured_a_rms_filtered = nothing; varnames{end+1} = 'measured_a_rms_filtered';
    measured_a_mean_filtered = nothing; varnames{end+1} = 'measured_a_mean_filtered';
    measured_a_std_filtered = nothing; varnames{end+1} = 'measured_a_std_filtered';
    measured_lambda = nothing; varnames{end+1} = 'measured_lambda';
    measured_dof_mse = nothing; varnames{end+1} = 'measured_dof_mse';
    measured_dof_se = nothing; varnames{end+1} = 'measured_dof_se';
    measured_dof_var = nothing; varnames{end+1} = 'measured_dof_var';
    
    for iSlope = 1:length(slopes)  
        slope = slopes(iSlope);    
        fprintf('slope %d, ',slope); 
        
        for iStride=1:length(strides)
            stride = strides(iStride);
            dt = stride*base_dt;
            fprintf('stride %d, ',stride);
            
            measurement_time(iStride,iSlope) = dt;
            
            for iEnsemble=1:totalEnsembles
                fprintf('..%d',iEnsemble);
                
                linearIndex = sub2ind(size(nothing),iStride,iSlope,iEnsemble);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Generate the signal
                sampleFactor = 10; % 2500 points takes about 1/4 second to generate maternoise
                cv=maternoise(dt/sampleFactor,sampleFactor*n,sigma_u*sqrt(2),abs(slope),1/t_damp);
                cx = cumtrapz(cv)*dt/sampleFactor;
                t_all = (dt/sampleFactor)*(0:(sampleFactor*n-1))';
                x_all = real(cx);
                data = struct('t',dt*(0:n-1)','x',x_all(1:sampleFactor:end));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Generate the noise
                sigma = 10;
                distribution = NormalDistribution(sigma);
                epsilon = distribution.rand(size(data.x));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % t_obs and x_obs contain the times observed signal 
                x_obs = data.x + epsilon;
                t_obs = data.t;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Fit a constrained tension spline (a polynomial fit)
                K = S+1;
                t_knot = cat(1,min(t_obs)*ones(K+1,1),max(t_obs)*ones(K+1,1));
                mean_spline = ConstrainedSpline(t_obs,x_obs,K+1,t_knot,distribution,[]);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Remove the fit from the *true* signal, and compute stats.
                x = data.x - mean_spline(t_obs);
                actual_u_rms(linearIndex) = rms( diff(x)/dt );
                actual_a_rms(linearIndex) = rms( diff(x,T)/dt^T );
                actual_a_mean(linearIndex) = mean( diff(x,T)/dt^T );
                actual_a_std(linearIndex) = std( diff(x,T)/dt^T );
                
                actual_ad_error(linearIndex) = NormalDistribution(std(diff(x)/dt)).andersonDarlingError(diff(x)/dt);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Remove the fit from the observed signal and record stats
                x_obs = x_obs - mean_spline(t_obs);
                
                measured_u_rms(linearIndex) = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_obs,sqrt(distribution.variance),1);
                [measured_a_rms(linearIndex),measured_a_std(linearIndex),measured_a_mean(linearIndex)] = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_obs,sqrt(distribution.variance),T);
                
                x_filtered = SmoothingSpline.RunningFilter(x_obs,11,'median');
                measured_u_rms_filtered(linearIndex) = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_filtered,sqrt(distribution.variance),1);
                [measured_a_rms_filtered(linearIndex),measured_a_std_filtered(linearIndex),measured_a_mean_filtered(linearIndex)] = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_filtered,sqrt(distribution.variance),T);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Do a spline fit
                spline = SmoothingSpline(t_obs,x_obs,distribution, 'S', S, 'T', T);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Minimize to the true points---at the observation times only
                measured_lambda(linearIndex) = spline.minimizeMeanSquareError(data.t,x);
                
                measured_dof_mse(linearIndex) = spline.effectiveSampleSizeFromExpectedMeanSquareError;
                measured_dof_se(linearIndex) = spline.effectiveSampleSizeFromVarianceOfTheMean;
                measured_dof_var(linearIndex) = spline.effectiveSampleSizeFromSampleVariance;
            end
            fprintf('..\n');
        end
    end
    
    save('DegreesOfFreedomEstimates.mat',varnames{:}); 
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

gamma = sigma./(actual_u_rms.*measurement_time);
gamma_mean = mean(gamma,3);

% lambda = (1-1/n)/(xt)
% xt*lambda = 1-1/n
% 1/n = 1-xt*lambda
% figure
% for iSlope = 1:length(slopes)
% %     x = (gamma_mean(:,iSlope))./(measured_a_rms_mean(:,iSlope).^2 .* measurement_time(:,iSlope) );
% %     y = measured_lambda_mean(:,iSlope);
%     
%     x = gamma_mean(:,iSlope);
%     y = (measured_a_rms_mean(:,iSlope).^2 ).*measured_lambda_mean(:,iSlope);
%     
%     scatter(x,y), hold on
% end

% biggest problem is estimating the rms value from a noisy derivative.
% using sleptap helps quite a bit, but the smallest strides still fail
% figure
% x = reshape(sigma./(measured_u_rms.*dt),[],1);
% y = reshape((measured_a_rms.^2 ).*measured_lambda,[],1);
% scatter(x,y)

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

measured_dof_se_mean = mean(measured_dof_se,3);
measured_dof_se_std = std(measured_dof_se,1,3);

x = mean(log(gamma),3);
y = mean(log(measured_dof_se),3);
sigma2 = std(log(measured_dof_se),1,3).^2/totalEnsembles;

% figure
m = zeros(length(slopes),1);
b = zeros(length(slopes),1);

t = linspace(min(x(:)),max(x(:)),10)';
for iSlope = 1:length(slopes)
%     errorbar(x(:,iSlope),y(:,iSlope),sqrt(sigma2(:,iSlope))), hold on
    
    if slopes(iSlope) == -2
        strides_used = 1:(length(strides)-3);
    else
        strides_used = 1:(length(strides)-2);
    end
    
    [m(iSlope),b(iSlope)] = LinearBestFitWithVariableError(x(strides_used,iSlope),y(strides_used,iSlope),sigma2(strides_used,iSlope));
    
%     ax = gca;
%     ax.ColorOrderIndex = ax.ColorOrderIndex-1;
%     plot(t, max(0,m(iSlope)*t + b(iSlope)))
end

% The test statistic, z
z=(1./gamma).*(measured_dof_se.^(3/2))/sqrt(2);




FigureSize = [50 50 figure_width_1col+7 150*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize,'Name','DOFvsGamma');
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];



for iSlope = 1:length(slopes)
    scatter(reshape(gamma(:,iSlope,:),[],1),reshape(measured_dof_se(:,iSlope,:),[],1),2,'filled'); hold on
end

ax = gca;
ax.ColorOrderIndex = 1;
plotHandles = zeros(length(slopes),1);
x = 10.^linspace(log10(min(gamma(:))),log10(max(gamma(:))),50);
for iSlope = 1:length(slopes) 
    plotHandles(iSlope) = plot(x, max(1,exp(b(iSlope))*(x).^m(iSlope)), 'LineWidth',1.5);
end

xlabel('\Gamma', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('n^{SE}_{eff}', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlog, ylog
legend_labels = cell(1,length(slopes));
for iSlope = 1:length(slopes)
    legend_labels{iSlope} = sprintf('\\omega^{%d}, %.0f \\Gamma^{%.2f}',slopes(iSlope),exp(b(iSlope)),m(iSlope));
end
leg = legend(plotHandles,legend_labels);
leg.FontSize = figure_legend_size;
leg.FontName = figure_font;
leg.Location = 'southeast';
set( gca, 'FontSize', figure_axis_tick_size);
xlim([1e-2 1e1])
ylim([0.9 1e1])


z0 = exp(b./m)/sqrt(2);

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

