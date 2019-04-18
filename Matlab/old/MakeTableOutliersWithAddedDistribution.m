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
    filename = 'MSETableOutliersWithAddedDistributionStudentT.mat';
else
    filename = 'MSETableOutliersWithAddedDistributionNormal.mat';
end

if exist(filename,'file')
    load(filename);
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
    robust_alpha10_optimal = nothing_struct; vars{end+1} = 'robust_alpha10_optimal';
    robust_alpha100_optimal = nothing_struct; vars{end+1} = 'robust_alpha100_optimal';
    robust_alpha1k_optimal = nothing_struct; vars{end+1} = 'robust_alpha1k_optimal';
    robust_alpha10k_optimal = nothing_struct; vars{end+1} = 'robust_alpha10k_optimal';
    robust_alpha100k_optimal = nothing_struct; vars{end+1} = 'robust_alpha100k_optimal';
        
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
                    
                    spline_optimal = SmoothingSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected);
                    spline_optimal.minimizeMeanSquareError(data.t,data.x);
                    optimal = LogStatisticsFromSplineForOutlierTable(optimal,linearIndex,spline_optimal,compute_ms_error,trueOutlierIndices);
                    
                    alpha = 1/10;
                    spline_robust_optimal = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'alpha',alpha);
                    spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
                    robust_alpha10_optimal = LogStatisticsFromSplineForOutlierTable(robust_alpha10_optimal,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                    
                    alpha = 1/100;
                    spline_robust_optimal = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'alpha',alpha);
                    spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
                    robust_alpha100_optimal = LogStatisticsFromSplineForOutlierTable(robust_alpha100_optimal,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                    
                    alpha = 1/1000;
                    spline_robust_optimal = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'alpha',alpha);
                    spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
                    robust_alpha1k_optimal = LogStatisticsFromSplineForOutlierTable(robust_alpha1k_optimal,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                    
                    alpha = 1/10000;
                    spline_robust_optimal = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'alpha',alpha);
                    spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
                    robust_alpha10k_optimal = LogStatisticsFromSplineForOutlierTable(robust_alpha10k_optimal,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
                    
                    alpha = 1/100000;
                    spline_robust_optimal = RobustSmoothingSpline(t_obs,x_obs,noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'alpha',alpha);
                    spline_robust_optimal.minimizeMeanSquareError(data.t,data.x);
                    robust_alpha100k_optimal = LogStatisticsFromSplineForOutlierTable(robust_alpha100k_optimal,linearIndex,spline_robust_optimal,compute_ms_error,trueOutlierIndices);
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

min_mse = min([optimal.mse(:), robust_alpha10_optimal.mse(:), robust_alpha100_optimal.mse(:), robust_alpha1k_optimal.mse(:), robust_alpha10k_optimal.mse(:), robust_alpha100k_optimal.mse(:)],[],2);
min_mse = reshape(min_mse,size(optimal.mse));

dmse = @(mse) mse./min_mse-1;
pct_range = 0.90;%0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

optimal.dmse = dmse(optimal.mse);
robust_alpha10_optimal.dmse = dmse(robust_alpha10_optimal.mse);
robust_alpha100_optimal.dmse = dmse(robust_alpha100_optimal.mse);
robust_alpha1k_optimal.dmse = dmse(robust_alpha1k_optimal.mse);
robust_alpha10k_optimal.dmse = dmse(robust_alpha10k_optimal.mse);
robust_alpha100k_optimal.dmse = dmse(robust_alpha100k_optimal.mse);

print_pct = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f (%.1f-%.1f) ',100*mean(stats.dmse(iOutlierRatio,iStride,iSlope,:)),minpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))), maxpct(sort(stats.dmse(iOutlierRatio,iStride,iSlope,:))) );

printcol = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %#.3g/%#.3g m$^2$ (%#.3g) (%d/%d) ', median(stats.mse(iOutlierRatio,iStride,iSlope,:)),mean(stats.mse(iOutlierRatio,iStride,iSlope,:)), median(stats.neff_se(iOutlierRatio,iStride,iSlope,:)), median(stats.false_positives(iOutlierRatio,iStride,iSlope,:)), median(stats.false_negatives(iOutlierRatio,iStride,iSlope,:)) );

all_maxpct = @(a,b,c,d,e,f,iOutlierRatio,iStride,iSlope) [maxpct(sort(a.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(b.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(c.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(d.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(e.dmse(iOutlierRatio,iStride,iSlope,:))),maxpct(sort(f.dmse(iOutlierRatio,iStride,iSlope,:)))];

ratioed_maxpct = @(a,b,c,d,e,f,iOutlierRatio) [maxpct(sort(reshape(a.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(b.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(c.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(d.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(e.dmse(iOutlierRatio,:,:,:),[],1))),maxpct(sort(reshape(f.dmse(iOutlierRatio,:,:,:),[],1)))];
ratioed_maxpct(robust_alpha10_optimal, robust_alpha100_optimal, robust_alpha1k_optimal, robust_alpha10k_optimal, robust_alpha100k_optimal,optimal,1)
ratioed_maxpct(robust_alpha10_optimal, robust_alpha100_optimal, robust_alpha1k_optimal, robust_alpha10k_optimal, robust_alpha100k_optimal,optimal,2)
ratioed_maxpct(robust_alpha10_optimal, robust_alpha100_optimal, robust_alpha1k_optimal, robust_alpha10k_optimal, robust_alpha100k_optimal,optimal,3)

all_all_maxpct = @(a,b,c,d,e,f) [maxpct(sort(a.dmse(:))),maxpct(sort(b.dmse(:))),maxpct(sort(c.dmse(:))),maxpct(sort(d.dmse(:))),maxpct(sort(e.dmse(:))),maxpct(sort(f.dmse(:)))];
all_all_maxpct(robust_alpha10_optimal, robust_alpha100_optimal, robust_alpha1k_optimal, robust_alpha10k_optimal, robust_alpha100k_optimal,optimal)

fprintf('\n\n');
fprintf('\\begin{tabular}{r | lllll} stride & optimal mse ($n_{eff}$) & blind optimal & robust & false pos/neg & robust 2nd & false pos/neg \\\\ \\hline \\hline \n');
for iOutlierRatio = 1:totalOutlierRatios
    for iSlope = 1:length(slopes)
        fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
        for iStride=1:length(strides)
            fprintf('%d ', strides(iStride));
            
            print_pct(robust_alpha10_optimal,iOutlierRatio,iStride,iSlope);
            print_pct(robust_alpha100_optimal,iOutlierRatio,iStride,iSlope);
            print_pct(robust_alpha1k_optimal,iOutlierRatio,iStride,iSlope);
            print_pct(robust_alpha10k_optimal,iOutlierRatio,iStride,iSlope);
            print_pct(robust_alpha100k_optimal,iOutlierRatio,iStride,iSlope);
            print_pct(optimal,iOutlierRatio,iStride,iSlope);
            
            the_maxpct = all_maxpct(robust_alpha10_optimal, robust_alpha100_optimal, robust_alpha1k_optimal, robust_alpha10k_optimal, robust_alpha100k_optimal,optimal,iOutlierRatio,iStride,iSlope);
            [themin, indices] = sort(the_maxpct);
            
            fprintf('&\t %d, %d, %d ',indices(1),indices(2),indices(3));
            
            fprintf(' \\\\ \n');
            
        end
    end
    fprintf('\n\n');
end
fprintf('\\end{tabular} \n');

