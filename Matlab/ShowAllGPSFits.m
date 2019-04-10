% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

drifters = open('sample_data/raw_rho1_drifters.mat');

splines = cell(1,length(drifters.lon));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We need a consistent lon0 for all drifters
lon = drifters.lon{1};
lon0 = min(lon)+(max(lon)-min(lon))/2;

for iDrifter = 1:length(drifters.lon)
    t_drifter = (drifters.date{iDrifter}-drifters.lastDeployment)*24*60*60;
    lat = drifters.lat{iDrifter};
    lon = drifters.lon{iDrifter};
    spline = GPSTensionSpline(t_drifter,lat,lon,'shouldUseRobustFit',1,'lon0',lon0);
    
    fprintf('n_eff (x,y): (%.1f, %.1f)\n',spline.spline_x.effectiveSampleSizeFromVarianceOfTheMean,spline.spline_y.effectiveSampleSizeFromVarianceOfTheMean);
    fprintf('lambda: %.3g\n',spline.lambda);

    
    
    splines{iDrifter} = spline;
    tq = linspace(min(spline.t),max(spline.t),10*length(spline.t));
    [xq,yq] = spline.xyAtTime(tq);
    
    if iDrifter == 1
        figure
        plot(xq,yq,'LineWidth',1.5),axis equal, hold on
    else
        plot(xq,yq,'LineWidth',1.5)
    end
    
    spline.estimateOutlierDistribution();
    fprintf('alpha: %.2f',spline.alpha);
    if spline.alpha > 0
        fprintf('(nu,sigma): (%.2f,%.2f)\n',spline.outlierDistribution.nu,spline.outlierDistribution.sigma);
    else
        fprintf('\n');
    end
    
end

figure
subplot(3,1,1)
plot(tq,spline.spline_mean_x(tq)), hold on
scatter(spline.spline_mean_x.t,spline.spline_mean_x.x)
scatter(spline.spline_mean_x.t,spline.spline_mean_x.smoothingMatrix*spline.spline_mean_x.x)

subplot(3,1,2)
plot(tq,spline.spline_x(tq)), hold on
scatter(spline.spline_x.t,spline.spline_x.x)
scatter(spline.spline_x.t,spline.spline_x.smoothingMatrix*spline.spline_x.x)

subplot(3,1,3)
plot(tq,xq), hold on
scatter(spline.t,spline.x)
scatter(spline.t,spline.smoothingMatrixX*spline.x)