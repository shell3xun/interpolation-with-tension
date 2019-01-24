scaleFactor = 1;
LoadFigureDefaults

shouldUseStudentTDistribution = 1;

if shouldUseStudentTDistribution == 1
    filename = 'MSEComparisonTableForOutliersFullStackStudentT.mat';
else
    filename = 'MSEComparisonTableForOutliersFullStackNormal.mat';
end

if exist(filename,'file')
    load(filename);
else
    slopes = [-2; -3; -4];
%     slopes = -3;
    totalSlopes = length(slopes);

    result_stride = [5;20;200];
%     result_stride = 200;
    totalStrides = length(result_stride);
    totalEnsembles = 11; % best to choose an odd number for median

    % spline fit parameters
    S = 2;
    T = S;
    K = S+1;
    
    vars = {'S', 'T', 'slopes', 'result_stride'};
    
    % matern signal parameters
    sigma_u = 0.20;
    base_dt = 5; % for whatever reason, we chose this as the primary dt
    t_damp = 30*60;
    n = 250;

    % outlier parameters
    percentOutliers = 0.15;
    outlierFactor = 50;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preallocate the variables we need to save
    %
    nothing = zeros(totalStrides, totalSlopes, totalEnsembles);
    nothing_struct = struct('mse',nothing,'neff_se',nothing,'false_negatives', nothing, 'false_positives', nothing);
    
    total_outliers = nothing; vars{end+1} = 'total_outliers';
    
%     optimal = nothing_struct; vars{end+1} = 'optimal';
%     robust_alpha100_optimal = nothing_struct; vars{end+1} = 'robust_alpha100_optimal';
%     robust_alpha100_knots_removed_optimal = nothing_struct; vars{end+1} = 'robust_alpha100_knots_removed_optimal';
    robust_alpha10k_optimal = nothing_struct; vars{end+1} = 'robust_alpha10k_optimal';
    
    robust_alpha10k_blind_beta50 = nothing_struct; vars{end+1} = 'robust_alpha10k_blind_beta50';
    robust_alpha10k_blind_beta200 = nothing_struct; vars{end+1} = 'robust_alpha10k_blind_beta200';
    robust_alpha10k_blind_beta100 = nothing_struct; vars{end+1} = 'robust_alpha10k_blind_beta100';
    
    robust_alpha10k_blind_kr_beta100 = nothing_struct; vars{end+1} = 'robust_alpha10k_blind_kr_beta100';

        
    for iSlope = 1:length(slopes)
        slope = slopes(iSlope);
        fprintf('slope %d, ',slope);
                
        for iStride=1:length(result_stride)
            stride = result_stride(iStride);
            dt = stride*base_dt;
            fprintf('stride %d, ',stride);
            
            for iEnsemble = 1:totalEnsembles
                fprintf('..%d',iEnsemble);
                
                % Generate the signal
                cv=maternoise(dt,n,sigma_u*sqrt(2),abs(slope),1/t_damp);
                cx = cumtrapz(cv)*dt;
                data = struct('t',dt*(0:n-1)','x',real(cx));
                
                compute_ms_error = @(spline) (mean(mean(  (data.x - spline(data.t)).^2,2 ),1));
                                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Generate the noise
                if shouldUseStudentTDistribution == 1
                    nu = 4.5; sigma =  8.5;
                    noiseDistribution = StudentTDistribution(sigma,nu);
                else
                    sigma = 10;
                    noiseDistribution = NormalDistribution(sigma);
                end
                outlierIndices = rand(n,1)<=percentOutliers;
                outlierIndices([1 n]) = 0;
                
                epsilon_noise = zeros(n,1);
                epsilon_outlier = zeros(n,1);
                
                epsilon_noise(~outlierIndices) = noiseDistribution.rand(sum(~outlierIndices));
                epsilon_outlier(outlierIndices) = outlierFactor*noiseDistribution.rand(sum(outlierIndices));
                
                epsilon = epsilon_noise + epsilon_outlier;
                
                outlierThreshold = noiseDistribution.locationOfCDFPercentile(1-1/10000/2);
                trueOutlierIndices = find(abs(epsilon) > outlierThreshold);
                trueGoodIndices = setdiff(1:n,trueOutlierIndices);
                
                total_outliers(iStride,iSlope,iEnsemble) = length(trueOutlierIndices);
                
                x_obs = data.x + epsilon;
                t_obs = data.t;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Unblinded best fit with standard tension spline
                
                linearIndex = sub2ind(size(nothing),iStride,iSlope,iEnsemble);
                
%                 spline_optimal = TensionSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected);
%                 spline_optimal.minimizeMeanSquareError(data.t,data.x);
%                 optimal = LogStatisticsFromSplineForOutlierTable(optimal,linearIndex,spline_optimal,compute_ms_error,trueOutlierIndices);
% 
%                 spline_robust_optimal = RobustTensionSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'alpha',1/100);
%                 spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
%                 robust_alpha100_optimal = LogStatisticsFromSplineForOutlierTable(robust_alpha100_optimal,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
%                 
%                 spline_robust_optimal.removeOutlierKnotsAndRetension(1/10000);
%                 spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
%                 robust_alpha100_knots_removed_optimal = LogStatisticsFromSplineForOutlierTable(robust_alpha100_knots_removed_optimal,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
%                                
                spline_robust_optimal = RobustTensionSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
                spline_robust_optimal.removeOutlierKnotsAndRetension(1/10000);
                spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
                robust_alpha10k_optimal = LogStatisticsFromSplineForOutlierTable(robust_alpha10k_optimal,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                
                beta = 1/50;
                zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                spline_robust_optimal.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                robust_alpha10k_blind_beta50 = LogStatisticsFromSplineForOutlierTable(robust_alpha10k_blind_beta50,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                
                beta = 1/200;
                zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                spline_robust_optimal.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                robust_alpha10k_blind_beta200 = LogStatisticsFromSplineForOutlierTable(robust_alpha10k_blind_beta200,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                
                beta = 1/100;
                zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                spline_robust_optimal.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                robust_alpha10k_blind_beta100 = LogStatisticsFromSplineForOutlierTable(robust_alpha10k_blind_beta100,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                
                spline_robust_optimal.removeOutlierKnotsAndRetension(beta);
                robust_alpha10k_blind_kr_beta100 = LogStatisticsFromSplineForOutlierTable(robust_alpha10k_blind_kr_beta100,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);

%                 spline_robust_optimal.removeOutlierKnotsAndRetension(1/10000);
%                 spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
%                 robust_alpha10k_knots_removed_optimal = LogStatisticsFromSplineForOutlierTable(robust_alpha10k_knots_removed_optimal,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                
            end
            fprintf('\n');
        end
    end
    
    save(filename, vars{:});
end

dmse = @(mse) mse./robust_alpha10k_optimal.mse-1;
pct_range = 0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

robust_alpha10k_blind_beta50.dmse = dmse(robust_alpha10k_blind_beta50.mse);
robust_alpha10k_blind_beta200.dmse = dmse(robust_alpha10k_blind_beta200.mse);
robust_alpha10k_blind_beta100.dmse = dmse(robust_alpha10k_blind_beta100.mse);
robust_alpha10k_blind_kr_beta100.dmse = robust_alpha10k_blind_kr_beta100.mse./robust_alpha10k_blind_beta100.mse-1;

printcol = @(stats,iStride,iSlope) fprintf('& %#.3g/%#.3g m$^2$ (%#.3g) (%d/%d) ', median(stats.mse(iStride,iSlope,:)),mean(stats.mse(iStride,iSlope,:)), median(stats.neff_se(iStride,iSlope,:)), median(stats.false_positives(iStride,iSlope,:)), median(stats.false_negatives(iStride,iSlope,:)) );
print_pct = @(stats,iStride,iSlope) fprintf('&  %.1f-%.1f ',minpct(sort(stats.dmse(iStride,iSlope,:))), maxpct(sort(stats.dmse(iStride,iSlope,:))) );

fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & optimal mse ($n_{eff}$) & blind optimal & robust & false pos/neg & robust 2nd & false pos/neg \\\\ \\hline \\hline \n');
for iSlope = 1:length(slopes)
    fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
    for iStride=1:length(result_stride)
        fprintf('%d ', result_stride(iStride));

        printcol(robust_alpha10k_optimal,iStride,iSlope);
        print_pct(robust_alpha10k_blind_beta50,iStride,iSlope);
        print_pct(robust_alpha10k_blind_beta200,iStride,iSlope);
        print_pct(robust_alpha10k_blind_beta100,iStride,iSlope);
        print_pct(robust_alpha10k_blind_kr_beta100,iStride,iSlope);
                
        fprintf(' \\\\ \n');

    end
    
end
fprintf('\\end{tabular} \n');

