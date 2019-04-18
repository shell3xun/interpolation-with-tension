%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MakeTableOutliersWithAddedDistributionBlind
%
% This is the first (of two) scripts to test parameters used in the
% manuscript.
%
% This script generates a signal using the Matern with three different
% slopes. It then adds noise to the signal that include a known
% (quantified) portion and an outlier portion. The script tries a variety
% of added distributions to see which distribution performs the best in
% combination with a variety of range-restricted expected mean square error
% minimizers.
%
% For student-t errors:
% Overall, in terms of the smallest 90th percentile mean-square-error
% *not* using an added distribution and simply doing a range restricted
% minimization on the 99th percentile of expected noise is optimal. In
% fact, it appears to outperform the normal minimizer when there are no
% outliers.
%
% 2019/03/25 - JJE

scaleFactor = 1;
LoadFigureDefaults

shouldUseStudentTDistribution = 1;

if shouldUseStudentTDistribution == 1
    filename = 'MSETableOutliersWithAddedDistributionBlindStudentT.mat';
else
    filename = 'MSETableOutliersWithAddedDistributionBlindNormal.mat';
end

if exist(filename,'file')
    load(filename);
    totalOutlierRatios = size(total_outliers,1);
else
    slopes = [-2; -3; -4];
    slopes = -3;
    totalSlopes = length(slopes);

    strides = [5;20;80;200];
%     strides = 200;
    totalStrides = length(strides);
    totalEnsembles = 51; % best to choose an odd number for median
    
    outlierRatios = [0.0 0.05 0.15 0.25];
%     outlierRatios = 0.15;
    totalOutlierRatios = length(outlierRatios);
    
    % spline fit parameters
    S = 2;
    T = S;
    K = S+1;
    
    varnames = {'S', 'T', 'slopes', 'strides', 'outlierRatios', 'totalEnsembles'};
    
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
    nothing = nan(totalOutlierRatios, totalStrides, totalSlopes, totalEnsembles);
    nothing_struct = struct('mse',nothing,'neff_se',nothing,'lambda',nothing,'nonOutlierEffectiveSampleSize',nothing,'nonOutlierSampleVariance',nothing,'false_negatives', nothing, 'false_positives', nothing,'rejects',nothing);
    
    total_outliers = nothing; varnames{end+1} = 'total_outliers';
        
    stat_structs = cell(1,1);
    
    stat_structs{1} = nothing_struct; stat_structs{end}.name = 'nonrobust_expected_mse';
    
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha0_optimal';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha0_beta50';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha0_beta100';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha0_beta200';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha0_beta400';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha0_beta800';
    
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha100_optimal';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha100_beta50';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha100_beta100';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha100_beta200';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha100_beta400';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha100_beta800';
    
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha10000_optimal';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha10000_beta50';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha10000_beta100';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha10000_beta200';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha10000_beta400';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'robust_alpha10000_beta800';
    varnames{end+1} = 'stat_structs';
        
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
                    cv=maternoise(dt,n,sigma_u*sqrt(2),abs(slope),1/t_damp);
                    cx = cumtrapz(cv)*dt;
                    data = struct('t',dt*(0:n-1)','x',real(cx));
                    
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
                    
                    epsilon_noise(~outlierIndices) = noiseDistribution.rand(sum(~outlierIndices),1);
                    epsilon_outlier(outlierIndices) = outlierDistribution.rand(sum(outlierIndices),1);
                    
                    epsilon = epsilon_noise + epsilon_outlier;
                    
                    outlierThreshold = noiseDistribution.locationOfCDFPercentile(1-1/10000/2);
                    trueOutlierIndices = find(abs(epsilon) > outlierThreshold);
                    truenonOutlierIndices = setdiff(1:n,trueOutlierIndices);
                    
                    total_outliers(iOutlierRatio,iStride,iSlope,iEnsemble) = length(trueOutlierIndices);
                    
                    x_obs = data.x + epsilon;
                    t_obs = data.t;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Blinded best fit with standard tension spline
                    linearIndex = sub2ind(size(nothing),iOutlierRatio,iStride,iSlope,iEnsemble);
                    iStruct = 0;
                    
                    spline = SmoothingSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.optimalIterated);
                    iStruct = iStruct+1;
                    stat_structs{iStruct} = LogStatisticsFromSplineForOutlierTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Blinded best fit with robust tension splines
                    
                    alphaValues = [0, 1/100, 1/10000];
                    betaValues = 1./[50, 100, 200, 400, 800];
                    
                    for alpha=alphaValues                      
                        if alpha > 0
                            nu = 3.0;
                            sigma = sqrt(noiseDistribution.variance*1000*(nu-2)/nu);
                            distribution = AddedDistribution(alpha,StudentTDistribution(sigma,nu),noiseDistribution);
                        else
                            distribution = noiseDistribution;
                        end
                        
                        spline = SmoothingSpline(t_obs,x_obs,distribution, 'S', S, 'lambda',Lambda.fullTensionExpected);
                        spline.minimizeMeanSquareError(data.t,data.x);
                        iStruct = iStruct+1;
                        stat_structs{iStruct} = LogStatisticsFromSplineForOutlierTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error,trueOutlierIndices,outlierIndices);
                        
                        for beta=betaValues
                            zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                            zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                            expectedVariance = noiseDistribution.varianceInRange(zmin,zmax);
                            spline.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax,expectedVariance) );
                            
                            iStruct = iStruct+1;
                            stat_structs{iStruct} = LogStatisticsFromSplineForOutlierTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error,trueOutlierIndices,outlierIndices);
                        end
                    end

                end
                fprintf('\n');
            end
        end
    end
    
    save(filename, varnames{:});
end

% We find the lowest mean square error (mse) that was achieved with a blind
% minimization technique, and then use that as the basis for comparison.
optimal_blind_mse = min(cat(5,stat_structs{1}.mse,stat_structs{3}.mse,stat_structs{4}.mse,stat_structs{5}.mse,stat_structs{6}.mse,stat_structs{7}.mse,stat_structs{9}.mse,stat_structs{10}.mse,stat_structs{11}.mse,stat_structs{12}.mse,stat_structs{13}.mse,stat_structs{15}.mse,stat_structs{16}.mse,stat_structs{17}.mse,stat_structs{18}.mse,stat_structs{19}.mse),[],5);
dmse = @(stats) log10(stats.mse./optimal_blind_mse);
for i=1:length(stat_structs)
    stat_structs{i}.dmse = dmse(stat_structs{i});
end
pct_range = 0.9; %0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
minpct = @(values) 100*(10^values(ceil( ((1-pct_range)/2)*length(values)))-1);
maxpct = @(values) 100*(10^values(floor( ((1+pct_range)/2)*length(values)))-1);
print_dmse_outlier = @(stats,iOutlierRatio) fprintf('%s: %.1f (%.1f-%.1f)\n',stats.name ,100*(10^mean(reshape(stats.dmse(iOutlierRatio,:,:,:),[],1))-1),minpct(sort(reshape(stats.dmse(iOutlierRatio,:,:,:),[],1))), maxpct(sort(reshape(stats.dmse(iOutlierRatio,:,:,:),[],1))) );
for iOutlierRatio=1:length(outlierRatios)
    fprintf('---dmse for %.2f outliers---\n',outlierRatios(iOutlierRatio));
    for i=1:length(stat_structs)
        print_dmse_outlier( stat_structs{i},iOutlierRatio );
    end
    fprintf('\n\n');
end



return


pct_range = 0.90;%0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

minrange = @(values) values(ceil( ((1-pct_range)/2)*length(values)));
maxrange = @(values) values(floor( ((1+pct_range)/2)*length(values)));

print_mse = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f (%.1f-%.1f) ',100*mean(stats.mse(iOutlierRatio,iStride,iSlope,:)),minpct(sort(stats.mse(iOutlierRatio,iStride,iSlope,:))), maxpct(sort(stats.mse(iOutlierRatio,iStride,iSlope,:))) );

print_mse_all = @(stats) fprintf('%s: %.1f (%.1f-%.1f)\n',stats.name,sqrt(mean(stats.mse(:))),sqrt(minrange(sort(stats.mse(:)))), sqrt(maxrange(sort(stats.mse(:)))) );
fprintf('---rmse overall---\n');
for i=1:length(stat_structs)
    print_mse_all( stat_structs{i} );
end
fprintf('\n\n');

print_mse_outlier = @(stats,iOutlierRatio) fprintf('%s: %.1f (%.1f-%.1f)\n',stats.name ,sqrt(mean(reshape(stats.mse(iOutlierRatio,:,:,:),[],1))),sqrt(minrange(sort(reshape(stats.mse(iOutlierRatio,:,:,:),[],1)))), sqrt(maxrange(sort(reshape(stats.mse(iOutlierRatio,:,:,:),[],1))) ));

for iOutlierRatio=1:length(outlierRatios)
    fprintf('---rmse for %.2f outliers---\n',outlierRatios(iOutlierRatio));
    for i=1:length(stat_structs)
        print_mse_outlier( stat_structs{i},iOutlierRatio );
    end
    fprintf('\n\n');
end

fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & optimal mse ($n_{eff}$) & blind optimal & robust & false pos/neg & robust 2nd & false pos/neg \\\\ \\hline \\hline \n');
for iOutlierRatio = 1:totalOutlierRatios
    for iSlope = 1:length(slopes)
        fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
        for iStride=1:length(strides)
            fprintf('%d ', strides(iStride));
            
            for i=1:length(stat_structs)
                print_mse( stat_structs{i},iOutlierRatio,iStride,iSlope );
            end
            
            fprintf(' \\\\ \n');
            
        end
    end
    fprintf('\n\n');
end
fprintf('\\end{tabular} \n');

