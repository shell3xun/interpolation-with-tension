%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script generates a signal using the Matern with three different
% slopes. It then adds noise to the signal that include a known
% (quantified) portion and an outlier portion.

scaleFactor = 1;
LoadFigureDefaults

shouldUseStudentTDistribution = 1;

if shouldUseStudentTDistribution == 1
    filename = 'TableRobustAccelerationEstimateStudentT.mat';
else
    filename = 'TableRobustAccelerationEstimateNormal.mat';
end

if exist(filename,'file')
    load(filename);
else
    slopes = [-2; -3; -4];
%     slopes = -3;
    totalSlopes = length(slopes);

    strides = [5;20;200];
%     strides = 5;
%     strides = 200;
    totalStrides = length(strides);
    totalEnsembles = 101; % best to choose an odd number for median
    
    outlierRatios = [0; 0.015; 0.15; 0.25];
%     outlierRatios = 0.0;
    totalOutlierRatios = length(outlierRatios);
    
    % spline fit parameters
    S = 2;
    T = S;
    K = S+1;
    
    vars = {'S', 'T', 'slopes', 'strides', 'outlierRatios', 'totalEnsembles'};
    
    % matern signal parameters
    sigma_u = 0.20;
    base_dt = 5; % for whatever reason, we chose this as the primary dt
    t_damp = 30*60;
    n = 250;

    % outlier parameters
   
    outlierFactor = 50;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preallocate the variables we need to save
    %
    nothing = zeros(totalOutlierRatios, totalStrides, totalSlopes, totalEnsembles);
    
    total_outliers = nothing; vars{end+1} = 'total_outliers';
    u_rms_true = nothing; vars{end+1} = 'u_rms_true';
    a_rms_true = nothing; vars{end+1} = 'a_rms_true';
    
    u_rms_nofilter = nothing; vars{end+1} = 'u_rms_nofilter';
    a_rms_nofilter = nothing; vars{end+1} = 'a_rms_nofilter';
    
    u_rms_median5 = nothing; vars{end+1} = 'u_rms_median5';
    a_rms_median5 = nothing; vars{end+1} = 'a_rms_median5';
    
    u_rms_median11 = nothing; vars{end+1} = 'u_rms_median11';
    a_rms_median11 = nothing; vars{end+1} = 'a_rms_median11';
    
    u_rms_median15 = nothing; vars{end+1} = 'u_rms_median15';
    a_rms_median15 = nothing; vars{end+1} = 'a_rms_median15';
        
    for iOutlierRatio = 1:totalOutlierRatios
        percentOutliers = outlierRatios(iOutlierRatio);
        fprintf('percentOutliers %.1f, ',100*percentOutliers);
        
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
                    m=4; % number of extra points to use for integration
                    cv=maternoise(dt/m,m*n,sigma_u*sqrt(2),abs(slope),1/t_damp);
                    cx = cumtrapz(cv)*dt/m;
                    data = struct('t',dt*(0:n-1)','x',real(cx(1:m:end)));
                    
                    compute_ms_error = @(spline) (mean(mean(  (data.x - spline(data.t)).^2,2 ),1));
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Generate the noise
                    if shouldUseStudentTDistribution == 1
                        nu = 4.5; sigma =  8.5;
                        noiseDistribution = StudentTDistribution(sigma,nu);
                        outlierDistribution = StudentTDistribution(outlierFactor*sigma,3.0);
                    else
                        sigma = 10;
                        noiseDistribution = NormalDistribution(sigma);
                        outlierDistribution = NormalDistribution(outlierFactor*sigma);
                    end
                    outlierIndices = rand(n,1)<=percentOutliers;
                    outlierIndices([1 n]) = 0;
                    
                    epsilon_noise = zeros(n,1);
                    epsilon_outlier = zeros(n,1);
                    
                    epsilon_noise(~outlierIndices) = noiseDistribution.rand(sum(~outlierIndices));
                    epsilon_outlier(outlierIndices) = outlierDistribution.rand(sum(outlierIndices));
                    
                    epsilon = epsilon_noise + epsilon_outlier;
                    
                    outlierThreshold = noiseDistribution.locationOfCDFPercentile(1-1/10000/2);
                    trueOutlierIndices = find(abs(epsilon) > outlierThreshold);
                    truenonOutlierIndices = setdiff(1:n,trueOutlierIndices);
                    
                    total_outliers(iOutlierRatio,iStride,iSlope,iEnsemble) = length(trueOutlierIndices);
                    
                    x_obs = data.x + epsilon;
                    t_obs = data.t;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Unblinded best fit with standard tension spline
                    
                    D1 = SmoothingSpline.FiniteDifferenceMatrixNoBoundary(1,data.t,1);
                    DT = SmoothingSpline.FiniteDifferenceMatrixNoBoundary(T,data.t,1);
                    
                    rms = @(x) sqrt( mean( x.^2 ));
                    
                    u_rms_true(iOutlierRatio,iStride,iSlope,iEnsemble) = rms( D1*data.x );
                    a_rms_true(iOutlierRatio,iStride,iSlope,iEnsemble) = rms( DT*data.x );
                    
                    u_rms_nofilter(iOutlierRatio,iStride,iSlope,iEnsemble) = rms( D1*x_obs );
                    a_rms_nofilter(iOutlierRatio,iStride,iSlope,iEnsemble) = rms( DT*x_obs );
                    
                    x_filtered = RunningFilter(x_obs,5,'median');
                    u_rms_median5(iOutlierRatio,iStride,iSlope,iEnsemble) = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(data.t,x_filtered,sqrt(noiseDistribution.variance),1);
                    a_rms_median5(iOutlierRatio,iStride,iSlope,iEnsemble) = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(data.t,x_filtered,sqrt(noiseDistribution.variance),T);
                    
                    x_filtered = RunningFilter(x_obs,11,'median');                    
                    u_rms_median11(iOutlierRatio,iStride,iSlope,iEnsemble) = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(data.t,x_filtered,sqrt(noiseDistribution.variance),1);
                    a_rms_median11(iOutlierRatio,iStride,iSlope,iEnsemble) = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(data.t,x_filtered,sqrt(noiseDistribution.variance),T);
                    
                    x_filtered = RunningFilter(x_obs,15,'median');
                    u_rms_median15(iOutlierRatio,iStride,iSlope,iEnsemble) = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(data.t,x_filtered,sqrt(noiseDistribution.variance),1);
                    a_rms_median15(iOutlierRatio,iStride,iSlope,iEnsemble) = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(data.t,x_filtered,sqrt(noiseDistribution.variance),T);
                end
                fprintf('\n');
            end
        end
    end
    
     save(filename, vars{:});
end

du = @(u) u./u_rms_true-1;
da = @(a) a./a_rms_true-1;

du_rms_nofilter = du(u_rms_nofilter);
da_rms_nofilter = da(a_rms_nofilter);
du_rms_median5 = du(u_rms_median5);
da_rms_median5 = da(a_rms_median5);
du_rms_median11 = du(u_rms_median11);
da_rms_median11 = da(a_rms_median11);
du_rms_median15 = du(u_rms_median15);
da_rms_median15 = da(a_rms_median15);

pct_range = 0.6827;
minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

print_pct = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f (%.1f-%.1f) ',100*mean(stats(iOutlierRatio,iStride,iSlope,:)),minpct(sort(stats(iOutlierRatio,iStride,iSlope,:))), maxpct(sort(stats(iOutlierRatio,iStride,iSlope,:))) );

print_pct_all = @(stats) fprintf('& %.1f (%.1f-%.1f) ',100*median(stats),minpct(sort(stats)), maxpct(sort(stats)) );
print_pct_all(du_rms_nofilter(:));
print_pct_all(da_rms_nofilter(:));
print_pct_all(du_rms_median5(:));
print_pct_all(da_rms_median5(:));
print_pct_all(du_rms_median11(:));
print_pct_all(da_rms_median11(:));
print_pct_all(du_rms_median15(:));
print_pct_all(da_rms_median15(:));
fprintf('\n')

for iOutlierRatio = 1:totalOutlierRatios
    for iSlope = 1:length(slopes)
        fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
        for iStride=1:length(strides)
            fprintf('%d ', strides(iStride));
                        
            print_pct(du_rms_nofilter,iOutlierRatio,iStride,iSlope);
            print_pct(da_rms_nofilter,iOutlierRatio,iStride,iSlope);
            print_pct(du_rms_median5,iOutlierRatio,iStride,iSlope);
            print_pct(da_rms_median5,iOutlierRatio,iStride,iSlope);
            print_pct(du_rms_median11,iOutlierRatio,iStride,iSlope);
            print_pct(da_rms_median11,iOutlierRatio,iStride,iSlope);
            print_pct(du_rms_median15,iOutlierRatio,iStride,iSlope);
            print_pct(da_rms_median15,iOutlierRatio,iStride,iSlope);
     
            fprintf('\n');
        end
    end
end

