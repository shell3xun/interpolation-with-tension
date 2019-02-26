%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script generates a signal using the Matern with three different
% slopes. It then adds noise to the signal that include a known
% (quantified) portion and an outlier portion. The script tries a variety
% of added distributions to see which distribution performs the best.
%
% Note that this script uses *unblinded* data to minimize. So here we're
% looking to see which distribution performs best. We later test minimizers
% to see which performs best.
%
% Overall, in terms of the smallest 90th percentile mean-square-error, the
% added distribution with 1/10000 works best.
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
%     slopes = -3;
    totalSlopes = length(slopes);

    strides = [5;20;200];
%     strides = 200;
    totalStrides = length(strides);
    totalEnsembles = 51; % best to choose an odd number for median
    
    outlierRatios = [0; 0.015; 0.15];
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
    nothing_struct = struct('mse',nothing,'neff_se',nothing,'false_negatives', nothing, 'false_positives', nothing);
    
    total_outliers = nothing; vars{end+1} = 'total_outliers';
    
    optimal = nothing_struct; vars{end+1} = 'optimal';
    robust_beta50_optimal = nothing_struct; vars{end+1} = 'robust_beta50_optimal';
    robust_beta100_optimal = nothing_struct; vars{end+1} = 'robust_beta100_optimal';
    robust_beta200_optimal = nothing_struct; vars{end+1} = 'robust_beta200_optimal';
    robust_beta400_optimal = nothing_struct; vars{end+1} = 'robust_beta400_optimal';
    robust_beta800_optimal = nothing_struct; vars{end+1} = 'robust_beta800_optimal';
        
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
                    truenonOutlierIndices = setdiff(1:n,trueOutlierIndices);
                    
                    total_outliers(iOutlierRatio,iStride,iSlope,iEnsemble) = length(trueOutlierIndices);
                    
                    x_obs = data.x + epsilon;
                    t_obs = data.t;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Unblinded best fit with standard tension spline
                    
                    linearIndex = sub2ind(size(nothing),iOutlierRatio,iStride,iSlope,iEnsemble);
                    
                    alpha = 1/10000;
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'alpha',alpha);
                    spline_robust.minimizeMeanSquareError(data.t,data.x);
                    optimal = LogStatisticsFromSplineForOutlierTable(optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices);
                    
                    beta = 1/50;
                    zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                    zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                    robust_beta50_optimal = LogStatisticsFromSplineForOutlierTable(robust_beta50_optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices);
                    
                    beta = 1/100;
                    zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                    zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                    robust_beta100_optimal = LogStatisticsFromSplineForOutlierTable(robust_beta100_optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices);
                    
                    beta = 1/200;
                    zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                    zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                    robust_beta200_optimal = LogStatisticsFromSplineForOutlierTable(robust_beta200_optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices);
                    
                    beta = 1/400;
                    zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                    zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                    robust_beta400_optimal = LogStatisticsFromSplineForOutlierTable(robust_beta400_optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices);
                    
                    beta = 1/800;
                    zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
                    zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
                    spline_robust.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
                    robust_beta800_optimal = LogStatisticsFromSplineForOutlierTable(robust_beta800_optimal,linearIndex,spline_robust,compute_ms_error,trueOutlierIndices);
                end
                fprintf('\n');
            end
        end
    end
    
    save(filename, vars{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure

% figure
% scatter(t_obs(spline_robust_optimal.outlierIndices),x_obs(spline_robust_optimal.outlierIndices),(8.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'), hold on
% scatter(t_obs(trueOutlierIndices),x_obs(trueOutlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
% scatter(t_obs,x_obs,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
% tq = linspace(min(t_obs),max(t_obs),10*length(t_obs));
% plot(tq,spline_optimal(tq),'k')
% plot(tq,spline_robust_optimal(tq),'b')
% % plot(tq,spline_robust_cv(tq),'m')
% return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dmse = @(mse) mse./optimal.mse-1;
pct_range = 0.90;%0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

optimal.dmse = dmse(optimal.mse);
robust_beta50_optimal.dmse = dmse(robust_beta50_optimal.mse);
robust_beta100_optimal.dmse = dmse(robust_beta100_optimal.mse);
robust_beta200_optimal.dmse = dmse(robust_beta200_optimal.mse);
robust_beta400_optimal.dmse = dmse(robust_beta400_optimal.mse);
robust_beta800_optimal.dmse = dmse(robust_beta800_optimal.mse);

print_pct = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f (%.1f-%.1f) ',100*mean(stats.dmse(iOutlierRatio,iStride,iSlope,:)),minpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))), maxpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))) );

printcol = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %#.3g/%#.3g m$^2$ (%#.3g) (%d/%d) ', median(stats.mse(iOutlierRatio,iStride,iSlope,:)),mean(stats.mse(iOutlierRatio,iStride,iSlope,:)), median(stats.neff_se(iOutlierRatio,iStride,iSlope,:)), median(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)) );

all_maxpct = @(a,b,c,d,e,f,iOutlierRatio,iStride,iSlope) [maxpct(sort(a.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(b.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(c.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(d.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(e.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(f.dmse(iOutlierRatio,iStride,iSlope,:)))];

ratioed_maxpct = @(a,b,c,d,e,f,iOutlierRatio) [maxpct(sort(reshape(a.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(b.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(c.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(d.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(e.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(f.dmse(iOutlierRatio,:,:,:),[],1)))];
ratioed_maxpct(robust_beta50_optimal, robust_beta100_optimal, robust_beta200_optimal, robust_beta400_optimal, robust_beta800_optimal,robust_beta800_optimal,1)
ratioed_maxpct(robust_beta50_optimal, robust_beta100_optimal, robust_beta200_optimal, robust_beta400_optimal, robust_beta800_optimal,robust_beta800_optimal,2)
ratioed_maxpct(robust_beta50_optimal, robust_beta100_optimal, robust_beta200_optimal, robust_beta400_optimal, robust_beta800_optimal,robust_beta800_optimal,3)


all_all_maxpct = @(a,b,c,d,e,f) [maxpct(sort(a.dmse(:))),maxpct(sort(b.dmse(:))),maxpct(sort(c.dmse(:))),maxpct(sort(d.dmse(:))),maxpct(sort(e.dmse(:))),maxpct(sort(f.dmse(:)))];
all_all_maxpct(robust_beta50_optimal, robust_beta100_optimal, robust_beta200_optimal, robust_beta400_optimal, robust_beta800_optimal,robust_beta800_optimal)

fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & optimal mse ($n_{eff}$) & blind optimal & robust & false pos/neg & robust 2nd & false pos/neg \\\\ \\hline \\hline \n');
for iOutlierRatio = 1:totalOutlierRatios
    for iSlope = 1:length(slopes)
        fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
        for iStride=1:length(strides)
            fprintf('%d ', strides(iStride));
            
%             print_pct(robust_beta50_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta100_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta200_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta400_optimal,iOutlierRatio,iStride,iSlope);
%             print_pct(robust_beta800_optimal,iOutlierRatio,iStride,iSlope);
            
            printcol(robust_beta50_optimal,iOutlierRatio,iStride,iSlope);
            printcol(robust_beta100_optimal,iOutlierRatio,iStride,iSlope);
            printcol(robust_beta200_optimal,iOutlierRatio,iStride,iSlope);
            printcol(robust_beta400_optimal,iOutlierRatio,iStride,iSlope);
            printcol(robust_beta800_optimal,iOutlierRatio,iStride,iSlope);
            
            
            the_maxpct = all_maxpct(robust_beta50_optimal, robust_beta100_optimal, robust_beta200_optimal, robust_beta400_optimal, robust_beta800_optimal,robust_beta800_optimal,iOutlierRatio,iStride,iSlope);
            [themin, indices] = sort(the_maxpct);
            
            fprintf('&\t %d, %d, %d ',indices(1),indices(2),indices(3));
            
            fprintf(' \\\\ \n');
            
        end
    end
    fprintf('\n\n');
end
fprintf('\\end{tabular} \n');

