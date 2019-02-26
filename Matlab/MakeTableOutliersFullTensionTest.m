%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script generates a signal using the Matern with three different
% slopes. It then adds noise to the signal that include a known
% (quantified) portion and an outlier portion.

scaleFactor = 1;
LoadFigureDefaults

shouldUseStudentTDistribution = 1;

if shouldUseStudentTDistribution == 1
    filename = 'MSETableOutliersFullTensionTestStudentT.mat';
else
    filename = 'MSETableOutliersFullTensionTestNormal.mat';
end

if exist(filename,'file')
    load(filename);
else
    slopes = [-2; -3; -4];
    slopes = -3;
    totalSlopes = length(slopes);

    strides = [5;20;200];
    strides = 5;
     strides = 200;
    totalStrides = length(strides);
    totalEnsembles = 11; % best to choose an odd number for median
    
%     outlierRatios = [0; 0.015; 0.15; 0.25];
    outlierRatios = 0.25;
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
    nothing_struct = struct('mse',nothing,'neff_se',nothing,'false_negatives', nothing, 'false_positives', nothing, 'lambda', nothing,'nonOutlierSampleVariance',nothing,'nonOutlierEffectiveSampleSize',nothing);
    
    total_outliers = nothing; vars{end+1} = 'total_outliers';
    
    full_tension_ks = nothing_struct; vars{end+1} = 'full_tension_ks';
    full_tension_sv = nothing_struct; vars{end+1} = 'full_tension_sv';
    full_tension_blind_ks_alpha8 = nothing_struct; vars{end+1} = 'full_tension_blind_ks_alpha8';
    full_tension_blind_ks_alpha10 = nothing_struct; vars{end+1} = 'full_tension_blind_ks_alpha10';
    full_tension_blind_ks_alpha20 = nothing_struct; vars{end+1} = 'full_tension_blind_ks_alpha20';
    full_tension_blind_ks_alpha100 = nothing_struct; vars{end+1} = 'full_tension_blind_ks_alpha100';
    full_tension_blind_sv_alpha8 = nothing_struct; vars{end+1} = 'full_tension_blind_sv_alpha8';
    full_tension_blind_sv_alpha10 = nothing_struct; vars{end+1} = 'full_tension_blind_sv_alpha10';
    full_tension_blind_sv_alpha20 = nothing_struct; vars{end+1} = 'full_tension_blind_sv_alpha20';
    full_tension_blind_sv_alpha100 = nothing_struct; vars{end+1} = 'full_tension_blind_sv_alpha100';
        
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
                    
                    epsilon_noise(~outlierIndices) = noiseDistribution.rand(sum(~outlierIndices));
                    epsilon_outlier(outlierIndices) = outlierDistribution.rand(sum(outlierIndices));
                    
                    epsilon = epsilon_noise + epsilon_outlier;
                    
                    outlierThreshold = noiseDistribution.locationOfCDFPercentile(1-1/10000/2);
                    trueOutlierIndices = find(abs(epsilon) > outlierThreshold);
                    truenonOutlierIndices = setdiff(1:n,trueOutlierIndices);
                    
                    total_outliers(iOutlierRatio,iStride,iSlope,iEnsemble) = length(trueOutlierIndices);
                    
                    f = @(z) abs( (1-percentOutliers)*noiseDistribution.pdf(z) - percentOutliers*outlierDistribution.pdf(z) );
                    z_crossover = abs(fminsearch(f,sqrt(noiseDistribution.variance)));
                    outlierCrossoverIndices = find(abs(epsilon)>z_crossover);
                    
                    x_obs = data.x + epsilon;
                    t_obs = data.t;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Unblinded best fit with standard tension spline
                    
                    linearIndex = sub2ind(size(nothing),iOutlierRatio,iStride,iSlope,iEnsemble);
                    
                    spline_ks = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
%                     spline_ks.minimize( @(spline) noiseDistribution.kolmogorovSmirnovError(spline.epsilonAtIndices(~outlierIndices)) );
%                     full_tension_ks = LogStatisticsFromSplineForOutlierTable(full_tension_ks,linearIndex,spline_ks,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    spline_ks.minimize( @(spline) mean((spline.epsilonAtIndices(~outlierIndices)-epsilon(~outlierIndices)).^2) );
%                     spline_ks.minimize(@(spline) mean((epsilon-spline.epsilon).^2));
                    fprintf(' total sample variance (requested/actual): (%.1f/%.1f)\t',mean(epsilon.^2) , mean(spline_ks.epsilon.^2));
                    fprintf('noise sample variance (requested/actual): (%.1f/%.1f)\n',mean(epsilon(~outlierIndices).^2) , mean(spline_ks.epsilonAtIndices(~outlierIndices).^2));
                    full_tension_sv = LogStatisticsFromSplineForOutlierTable(full_tension_sv,linearIndex,spline_ks,compute_ms_error,trueOutlierIndices,outlierIndices);
                    
                    % Now we do the blind tests
                    spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
                    lambda0 = spline_robust.lambda;
                    
                    alpha = 0.4;
                    spline_robust.setToFullTensionWithIQSV(alpha);
                    full_tension_blind_sv_alpha8 = LogStatisticsFromSplineForFullTensionTable(full_tension_blind_sv_alpha8,linearIndex,spline_robust,compute_ms_error,z_crossover,outlierCrossoverIndices,outlierIndices);
                    spline_robust.lambda = lambda0;
                    spline_robust.setToFullTensionWithSV(alpha);
                    full_tension_blind_ks_alpha8 = LogStatisticsFromSplineForFullTensionTable(full_tension_blind_ks_alpha8,linearIndex,spline_robust,compute_ms_error,z_crossover,outlierCrossoverIndices,outlierIndices);
                    
                    alpha = 1/2;
                    spline_robust.lambda = lambda0;
                    spline_robust.setToFullTensionWithIQSV(alpha);
                    full_tension_blind_sv_alpha10 = LogStatisticsFromSplineForFullTensionTable(full_tension_blind_sv_alpha10,linearIndex,spline_robust,compute_ms_error,z_crossover,outlierCrossoverIndices,outlierIndices);
                    spline_robust.lambda = lambda0;
                    spline_robust.setToFullTensionWithSV(alpha);
                    full_tension_blind_ks_alpha10 = LogStatisticsFromSplineForFullTensionTable(full_tension_blind_ks_alpha10,linearIndex,spline_robust,compute_ms_error,z_crossover,outlierCrossoverIndices,outlierIndices);
                    
                    alpha = 1/4;
                    spline_robust.lambda = lambda0;
                    spline_robust.setToFullTensionWithIQSV(alpha);
                    full_tension_blind_sv_alpha20 = LogStatisticsFromSplineForFullTensionTable(full_tension_blind_sv_alpha20,linearIndex,spline_robust,compute_ms_error,z_crossover,outlierCrossoverIndices,outlierIndices);
                    spline_robust.lambda = lambda0;
                    spline_robust.setToFullTensionWithSV(alpha);
                    full_tension_blind_ks_alpha20 = LogStatisticsFromSplineForFullTensionTable(full_tension_blind_ks_alpha20,linearIndex,spline_robust,compute_ms_error,z_crossover,outlierCrossoverIndices,outlierIndices);

                    alpha = 1/8;
                    spline_robust.lambda = lambda0;
                    spline_robust.setToFullTensionWithIQSV(alpha);
                    full_tension_blind_sv_alpha100 = LogStatisticsFromSplineForFullTensionTable(full_tension_blind_sv_alpha100,linearIndex,spline_robust,compute_ms_error,z_crossover,outlierCrossoverIndices,outlierIndices);
                    spline_robust.lambda = lambda0;
                    spline_robust.setToFullTensionWithSV(alpha);
                    full_tension_blind_ks_alpha100 = LogStatisticsFromSplineForFullTensionTable(full_tension_blind_ks_alpha100,linearIndex,spline_robust,compute_ms_error,z_crossover,outlierCrossoverIndices,outlierIndices);
                end
                fprintf('\n');
            end
        end
    end
    
%     save(filename, vars{:});
end

% figure, histogram(spline_ks.epsilonAtIndices(~outlierIndices),100,'Normalization','pdf'); z=linspace(-200,200,1000)'; hold on, plot(z,noiseDistribution.pdf(z));
% figure, histogram(spline_robust.epsilonAtIndices(~outlierIndices),100,'Normalization','pdf'); z=linspace(-200,200,1000)'; hold on, plot(z,noiseDistribution.pdf(z));

printall = @(stats) fprintf('%.1f/%.1f (%.1f/%.1f)\n', mean(stats.false_positives(:)), mean(stats.false_negatives(:)), median(stats.false_positives(:)), median(stats.false_negatives(:)));

printall(full_tension_blind_ks_alpha8);
printall(full_tension_blind_ks_alpha10);
printall(full_tension_blind_ks_alpha10);
printall(full_tension_blind_ks_alpha100);
printall(full_tension_blind_sv_alpha8);
printall(full_tension_blind_sv_alpha10);
printall(full_tension_blind_sv_alpha20);
printall(full_tension_blind_sv_alpha100);


opt_sampleVariance = full_tension_sv.nonOutlierSampleVariance;
dsv = @(sv) sv./opt_sampleVariance-1;

full_tension_sv.dSampleVariance = dsv(full_tension_sv.nonOutlierSampleVariance);

full_tension_blind_ks_alpha8.dSampleVariance = dsv(full_tension_blind_ks_alpha8.nonOutlierSampleVariance);
full_tension_blind_ks_alpha10.dSampleVariance = dsv(full_tension_blind_ks_alpha10.nonOutlierSampleVariance);
full_tension_blind_ks_alpha20.dSampleVariance = dsv(full_tension_blind_ks_alpha20.nonOutlierSampleVariance);
full_tension_blind_ks_alpha100.dSampleVariance = dsv(full_tension_blind_ks_alpha100.nonOutlierSampleVariance);

full_tension_blind_sv_alpha8.dSampleVariance = dsv(full_tension_blind_sv_alpha8.nonOutlierSampleVariance);
full_tension_blind_sv_alpha10.dSampleVariance = dsv(full_tension_blind_sv_alpha10.nonOutlierSampleVariance);
full_tension_blind_sv_alpha20.dSampleVariance = dsv(full_tension_blind_sv_alpha20.nonOutlierSampleVariance);
full_tension_blind_sv_alpha100.dSampleVariance = dsv(full_tension_blind_sv_alpha100.nonOutlierSampleVariance);

pct_range = 0.6827;
minpct = @(values) 100*values(ceil( ((1-pct_range)/2)*length(values)));
maxpct = @(values) 100*values(floor( ((1+pct_range)/2)*length(values)));

print_pct = @(stats,iOutlierRatio,iStride,iSlope) fprintf('& %.1f (%.1f-%.1f) ',100*mean(stats.dSampleVariance(iOutlierRatio,iStride,iSlope,:)),minpct(sort(stats.dSampleVariance(iOutlierRatio,iStride,iSlope,:))), maxpct(sort(stats.dSampleVariance(iOutlierRatio,iStride,iSlope,:))) );

for iOutlierRatio = 1:totalOutlierRatios
    for iSlope = 1:length(slopes)
        fprintf('$\\omega^{%d}$ &&&&&  \\\\ \\hline \n',slopes(iSlope));
        for iStride=1:length(strides)
            fprintf('%d ', strides(iStride));
            
            print_pct(full_tension_sv,iOutlierRatio,iStride,iSlope);
            
            print_pct(full_tension_blind_ks_alpha8,iOutlierRatio,iStride,iSlope);
            print_pct(full_tension_blind_ks_alpha10,iOutlierRatio,iStride,iSlope);
            print_pct(full_tension_blind_ks_alpha20,iOutlierRatio,iStride,iSlope);
            print_pct(full_tension_blind_ks_alpha100,iOutlierRatio,iStride,iSlope);
            
            print_pct(full_tension_blind_sv_alpha8,iOutlierRatio,iStride,iSlope);
            print_pct(full_tension_blind_sv_alpha10,iOutlierRatio,iStride,iSlope);
            print_pct(full_tension_blind_sv_alpha20,iOutlierRatio,iStride,iSlope);
            print_pct(full_tension_blind_sv_alpha100,iOutlierRatio,iStride,iSlope);
            
            fprintf('\n');
        end
    end
end

% figure, histogram(spline_ks.epsilonAtIndices(~outlierIndices),100,'Normalization','pdf'); z=linspace(-200,200,1000)'; hold on, plot(z,noiseDistribution.pdf(z));


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure

figure
scatter(t_obs(trueOutlierIndices),x_obs(trueOutlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
scatter(t_obs,x_obs,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
tq = linspace(min(t_obs),max(t_obs),10*length(t_obs));
plot(tq,spline_ks(tq),'k')
plot(tq,spline_robust(tq),'b')
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

