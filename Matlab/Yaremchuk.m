shouldUseStudentTDistribution = 1;

% slopes = [-2; -3; -4];
% strides = (2.^(0:4)).';

slopes = -3;
strides = 8;

totalSlopes = length(slopes);
totalStrides = length(strides);
totalEnsembles = 1; % best to choose an odd number for median

% spline fit parameters
S = 3;
T = S;
K = S+1;

varnames = {'S', 'T', 'slopes', 'strides', 'totalEnsembles'};

% matern signal parameters
sigma_u = 0.20;
base_dt = 1.5*60; % chosen as the smallest interval considered, because anything shorter than this looks non-stationary... like a local polynomial fit is needed.
t_damp = 30*60;
n = 250;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preallocate the variables we need to save
%
nothing = nan(totalStrides, totalSlopes, totalEnsembles);
nothing_struct = struct('mse',nothing,'neff_se',nothing,'lambda',nothing);

u_rms_true_strided = nothing; varnames{end+1} = 'u_rms_true_strided';
a_rms_true_strided = nothing; varnames{end+1} = 'a_rms_true_strided';
u_estimate_spectral = nothing; varnames{end+1} = 'u_estimate_spectral';
a_estimate_spectral = nothing; varnames{end+1} = 'a_estimate_spectral';

stat_structs = cell(1,1);
stat_structs{1} = nothing_struct; stat_structs{end}.name = 'true_optimal_knot_dof_1';
stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'optimal_expected';
stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'true_optimal_knot_dof_auto';
stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'optimal_iterated';
stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'optimal_ranged';
stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'optimal_cv';
stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'optimal_gcv';
stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'optimal_log_likelihood';
stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'optimal_iterated_mean_removed';
varnames{end+1} = 'stat_structs';

for iSlope = 1:length(slopes)
    slope = slopes(iSlope);
    fprintf('slope %d, ',slope);
    
    for iStride=1:length(strides)
        stride = strides(iStride);
        dt = stride*base_dt;
        fprintf('stride %d, ',stride);
        
        for iEnsemble = 1:totalEnsembles
            fprintf('..%d',iEnsemble);
            
            % Generate the signal
            cv=maternoise(dt,n,sigma_u*sqrt(2),abs(slope),1/t_damp);
            cx = cumtrapz(cv)*dt;
            data = struct('t',dt*(0:n-1)','x',real(cx));
            
            % record the 'true' (but strided) values of the signal
            rms = @(x) sqrt( mean( x.*x ) );
            u_rms_true_strided(iStride,iSlope,iEnsemble) = rms( diff(data.x)/dt );
            a_rms_true_strided(iStride,iSlope,iEnsemble) = rms( diff(diff(data.x))/(dt^2) );
            
            compute_ms_error = @(x) (mean(mean(  (data.x - x).^2,2 ),1));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Generate the noise
            if shouldUseStudentTDistribution == 1
                nu = 4.5; sigma =  8.5;
                noiseDistribution = StudentTDistribution(sigma,nu);
            else
                sigma = 10;
                noiseDistribution = NormalDistribution(sigma);
            end
            epsilon = noiseDistribution.rand(size(data.x));
            
            x_obs = data.x + epsilon;
            t_obs = data.t;
            
            % record the estimated values of the signal
            u_estimate_spectral(iStride,iSlope,iEnsemble) = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_obs,sqrt(noiseDistribution.variance),1);
            a_estimate_spectral(iStride,iSlope,iEnsemble) = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(t_obs,x_obs,sqrt(noiseDistribution.variance),T);
            
            linearIndex = sub2ind(size(nothing),iStride,iSlope,iEnsemble);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Unblinded best fit with knot_dof==1
            spline = SmoothingSpline(t_obs,x_obs,noiseDistribution, 'lambda', Lambda.fullTensionExpected, 'S', S, 'knot_dof', 1);
            spline.minimizeExpectedMeanSquareError;
            
            sigma_star = std(diff(diff(x_obs))/(dt^2));
            sigma_e = sqrt(noiseDistribution.variance); % std. position error
            sigma_p = sqrt(6)*sigma_e/(dt^2);
            mu = 0.91;
            C = gamma(1/mu)/gamma(2/mu);
            
            Wn = 1/(2*noiseDistribution.variance);
            Wa = (C*(dt^4)*(sigma_star^2 - sigma_p^2)).^mu;
            
            % I'm taking the abs here, which isn't in the paper, but must
            % be true
            J = @(x) ( Wa*sum(abs(diff(diff(x))).^(2*mu)) + Wn*sum( (x-x_obs).^2 ) )/length(x_obs);
            x_opt = fminsearch(J,x_obs,optimset('MaxFunEvals',1000*length(x_obs),'MaxIter',1000*length(x_obs)));
            
            fprintf('Yaremchuk rms error: %.2f\n',compute_ms_error(x_opt));
        end
        fprintf('\n');
    end
end

figure
scatter(t_obs,x_obs,(2.5)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
plot(t_obs,spline(t_obs),'k')
plot(t_obs,x_opt,'b')
legend('observed',sprintf('spline (%.2f m)\n',sqrt(compute_ms_error(spline(t_obs)))),sprintf('yaremchuk (%.2f m)\n',sqrt(compute_ms_error(x_opt))));

