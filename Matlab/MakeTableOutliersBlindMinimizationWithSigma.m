%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MakeTableOutliersBlindMinimizationWithSigma
%
% This is the second (of two) scripts to test parameters used in the
% manuscript.
%
% This script generates a signal using the Matern with three different
% slopes. It then adds noise to the signal that include a known
% (quantified) portion and an outlier portion.
%
% In the previous script (MakeTableOutliersWithAddedDistributionBlind) we
% showed that doing a ranged-expected mean square error minimization
% significantly improves the actual mse. 
%
% 2019/03/25 - JJE

shouldUseStudentTDistribution = 1;

if shouldUseStudentTDistribution == 1
    filename = 'MSETableOutliersBlindMinimizationWithSigmaStudentT.mat';
else
    filename = 'MSETableOutliersBlindMinimizationSigmaNormal.mat';
end

if exist(filename,'file')
    load(filename);
    totalOutlierRatios = size(total_outliers,1);
else
    slopes = [-2; -3; -4];
    slopes = -3;
    totalSlopes = length(slopes);

    strides = [5;20;80;200];
%      strides = 200;
    totalStrides = length(strides);
    totalEnsembles = 11; % best to choose an odd number for median
    
    outlierRatios = [0 0.05 0.15 0.25];
%      outlierRatios = 0.15;
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
   
    outlierFactor = 40;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preallocate the variables we need to save
    %
    nothing = nan(totalOutlierRatios, totalStrides, totalSlopes, totalEnsembles);
    nothing_struct = struct('mse',nothing,'neff_se',nothing,'lambda',nothing,'nonOutlierEffectiveSampleSize',nothing,'nonOutlierSampleVariance',nothing,'false_negatives', nothing, 'false_positives', nothing,'rejects',nothing);
    
    total_outliers = nothing; varnames{end+1} = 'total_outliers';
 
    stat_structs = cell(1,1);
    stat_structs{1} = nothing_struct; stat_structs{end}.name = 'minimization_noise_range';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'minimization_noise_range_sigma';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'minimization_ratio_1';
    stat_structs{end+1} = nothing_struct; stat_structs{end}.name = 'minimization_ratio_1_sigma';
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
                    % Unblinded best fit with standard tension spline
                    
                    linearIndex = sub2ind(size(nothing),iOutlierRatio,iStride,iSlope,iEnsemble);
                    
                    alpha = 0;
                    spline = TensionSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected);
                    spline.estimateOutlierDistribution();
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Range-restricted expected mean square error. This is
                    % our baseline.
                    beta = 1/100;
                    [lambda1,mse1] = spline.minimizeExpectedMeanSquareErrorInPercentileRange(beta/2,1-beta/2);
                    
                    iStruct = 1;
                    stat_structs{iStruct} = LogStatisticsFromSplineForOutlierTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Range-restricted expected mean square error, with
                    % different initial sigma.
                    spline.sigma = spline.sigmaAtFullTension;
                    spline.lambda = spline.lambdaAtFullTension;
                    [~,mse2] = spline.minimizeExpectedMeanSquareErrorInPercentileRange(beta/2,1-beta/2);
                    if mse1 < mse2
                        spline.sigma = sqrt(spline.distribution.variance);
                        spline.lambda = lambda1;
                    end
                    
                    iStruct = iStruct+1;
                    stat_structs{iStruct} = LogStatisticsFromSplineForOutlierTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Range-restricted expected mean square error using the
                    % pdf crossover point.
                    spline.sigma = sqrt(spline.distribution.variance);
                    spline.lambda = spline.lambdaAtFullTension;
                    
                    minimizationPDFRatio = 1;
                    [lambda1,mse1] = spline.minimizeExpectedMeanSquareErrorInNoiseRange(minimizationPDFRatio);
                    
                    iStruct = iStruct+1;
                    stat_structs{iStruct} = LogStatisticsFromSplineForOutlierTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
%                     fprintf('expected (%.1f, %.1f) actual (%.1f, %.1f)\n',mse0,mse1,stat_structs{iStruct-1}.mse(linearIndex),stat_structs{iStruct}.mse(linearIndex));

                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Using the full tension values, search for an
                    % alternative global minimum
                    
                    spline.sigma = spline.sigmaAtFullTension;
                    spline.lambda = spline.lambdaAtFullTension;
                    
                    [~,mse2] = spline.minimizeExpectedMeanSquareErrorInNoiseRange(minimizationPDFRatio);
                    if mse1 < mse2
                        spline.sigma = sqrt(spline.distribution.variance);
                        spline.lambda = lambda1;
                    end
                        
                    iStruct = iStruct+1;
                    stat_structs{iStruct} = LogStatisticsFromSplineForOutlierTable(stat_structs{iStruct},linearIndex,spline,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                end
                fprintf('\n');
            end
        end
    end
    
    save(filename, varnames{:});
end

% Analyze the mse

dmse = @(stats) log10(stats.mse./stat_structs{1}.mse);
pct_range = 0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
minpct = @(values) 100*(10^values(ceil( ((1-pct_range)/2)*length(values)))-1);
maxpct = @(values) 100*(10^values(floor( ((1+pct_range)/2)*length(values)))-1);

for i=1:length(stat_structs)
    stat_structs{i}.dmse = dmse(stat_structs{i});
end


print_dmse = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f (%.1f-%.1f) ',100*(10^mean(stats.dmse(iOutlierRatio,iStride,iSlope,:))-1),minpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))), maxpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))) );
printcol = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %#.3g/%#.3g m$^2$ (%#.3g) (%d/%d) ', median(stats.mse(iOutlierRatio,iStride,iSlope,:)),mean(stats.mse(iOutlierRatio,iStride,iSlope,:)), median(stats.neff_se(iOutlierRatio,iStride,iSlope,:)), median(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)) );
print_dmse_all = @(stats) fprintf('%s: & %.2f (%.2f) (%.2f-%.2f)\n',stats.name,100*(10^mean(stats.dmse(:))-1),100*(10^median(stats.dmse(:))-1),minpct(sort(stats.dmse(:))), maxpct(sort(stats.dmse(:))) );

fprintf('---dmse---\n');
for i=1:length(stat_structs)
    print_dmse_all( stat_structs{i} );
end

% Analyze the lambda---negative values indicate undertensioning, positive
% values indicate over tensioning
dlambda = @(stats) log10(stats.lambda ./ stat_structs{1}.lambda);
for i=1:length(stat_structs)
    stat_structs{i}.dlambda = dlambda(stat_structs{i});
end

print_dlambda = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.2f (%.2f) (%.2f-%.2f) ',mean(stats.dlambda(iOutlierRatio,iStride,iSlope,:)),median(stats.dlambda(iOutlierRatio,iStride,iSlope,:)),minpct(sort(stats.dlambda(iOutlierRatio,iStride,iSlope,:)))/100, maxpct(sort(stats.dlambda(iOutlierRatio,iStride,iSlope,:)))/100 );

print_dlambda_all = @(stats) fprintf('& %.2f (%.2f) (%.2f-%.2f)\n',mean(stats.dlambda(:)),median(stats.dlambda(:)),minpct(sort(stats.dlambda(:)))/100, maxpct(sort(stats.dlambda(:)))/100 );

fprintf('---dlambda---\n');
for i=1:length(stat_structs)
    print_dlambda_all( stat_structs{i} );
end

% all_maxpct = @(a,b,c,d,e,f,iOutlierRatio,iStride,iSlope) [maxpct(sort(a.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(b.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(c.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(d.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(e.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(f.dmse(iOutlierRatio,iStride,iSlope,:)))];
% 
% ratioed_maxpct = @(a,b,c,d,e,f,iOutlierRatio) [maxpct(sort(reshape(a.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(b.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(c.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(d.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(e.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(f.dmse(iOutlierRatio,:,:,:),[],1)))];
% ratioed_maxpct(robust_beta50_optimal, robust_beta100_optimal, robust_beta200_optimal, robust_beta400_optimal, robust_beta800_optimal,robust_beta800_optimal,1)
% ratioed_maxpct(robust_beta50_optimal, robust_beta100_optimal, robust_beta200_optimal, robust_beta400_optimal, robust_beta800_optimal,robust_beta800_optimal,2)
% ratioed_maxpct(robust_beta50_optimal, robust_beta100_optimal, robust_beta200_optimal, robust_beta400_optimal, robust_beta800_optimal,robust_beta800_optimal,3)
% 
% 
% all_all_maxpct = @(a,b,c,d,e,f) [maxpct(sort(a.dmse(:))),maxpct(sort(b.dmse(:))),maxpct(sort(c.dmse(:))),maxpct(sort(d.dmse(:))),maxpct(sort(e.dmse(:))),maxpct(sort(f.dmse(:)))];
% all_all_maxpct(robust_beta50_optimal, robust_beta100_optimal, robust_beta200_optimal, robust_beta400_optimal, robust_beta800_optimal,robust_beta800_optimal)

fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & optimal mse ($n_{eff}$) & blind optimal & robust & false pos/neg & robust 2nd & false pos/neg \\\\ \\hline \\hline \n');
for iOutlierRatio = 1:totalOutlierRatios
    for iSlope = 1:length(slopes)
        fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
        for iStride=1:length(strides)
            fprintf('%d ', strides(iStride));
                        
            for i=1:length(stat_structs)
                print_dmse( stat_structs{i},iOutlierRatio,iStride,iSlope );
            end
            
%             print_pct(robust_beta50_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta100_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta200_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta400_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta800_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_noiseOdds1over10,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_noiseOdds1over3,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_noiseOdds1,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_noiseOdds2,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_noiseOdds3,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_noiseOdds4,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_noiseOdds10,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_weighted,iOutlierRatio,iStride,iSlope);
            
%             print_dlambda(robust_noiseOdds1over3,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_noiseOdds1,iOutlierRatio,iStride,iSlope);
%             print_dlambda(robust_noiseOdds3,iOutlierRatio,iStride,iSlope);
            
            fprintf(' \\\\ \n');
            
        end
    end
    fprintf('\n\n');
end
fprintf('\\end{tabular} \n');

