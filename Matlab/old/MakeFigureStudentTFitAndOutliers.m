%%%%%%%
%
% OLD --- use MakeFigureGPSFitAndOutliers.m

% Site 1 or site 2?
Site = 1;

% It appears to be the case that you do NOT want to look for outliers in X
% and Y independently. 
IndependentXYOutlierRejection = 0;

% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

% Drifter to highlight in the final plots
choiceDrifter = 6;

shouldSaveFigures = 0;

if Site == 1
    drifters = load('sample_data/rho1_drifters_projected_ungridded.mat');
else
    drifters = load('sample_data/rho2_drifters_projected_ungridded.mat');
end

S = 3; % order of the spline
K = S+1;
T = S; % order of the tension
nu = 4.5; sigma =  8.5;
variance_of_the_noise = sigma*sigma*nu/(nu-2);
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
t_pdf = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));

% Pull out the data of interest
x = drifters.x{choiceDrifter};
y = drifters.y{choiceDrifter};
t = drifters.t{choiceDrifter};

spline_x = SmoothingSpline(t,x,sqrt(variance_of_the_noise), 'S', S, 'T', T,'weightFunction',w,'lambda',Lambda.fullTensionExpected);
spline_y = SmoothingSpline(t,y,sqrt(variance_of_the_noise), 'S', S, 'T', T,'weightFunction',w,'lambda',Lambda.fullTensionExpected);

% Now we set a threshold for what constitutes an outlier. In this case we
% choose points that have 1 in 10000 odds of occurring.
outlierThreshold = 0.0036;
gps_cdf = @(z) abs(tcdf(z/sigma,nu) - outlierThreshold/2);
range(1) = fminsearch( gps_cdf, -50, optimset('TolX', 0.001, 'TolFun', 0.001) );
range(2) = -range(1); % it's symmetric

% apply "full" tension
if 1 == 0
    % This applies full tension by matching the sample variance to the
    % total expected variance, BUT this doesn't work well because the
    % outliers add significant variance... so ultimately it ends up
    % under-tensioned.
    spline_x.Minimize( @(spline) abs(spline.SampleVariance - variance_of_the_noise) );
    spline_y.Minimize( @(spline) abs(spline.SampleVariance - variance_of_the_noise) );
else
    % This works better because here we require the distribution of errors
    % match the expected distribution, but restrict ourselves to a limit
    % range (and ignore the outliers).
    spline_x.Minimize( @(spline) KolmogorovSmirnovErrorForTDistribution(spline.epsilon,sigma,nu,range))
    spline_y.Minimize( @(spline) KolmogorovSmirnovErrorForTDistribution(spline.epsilon,sigma,nu,range))
end

% pull out the observed errors
epsilon_x = spline_x.epsilon;
epsilon_y = spline_y.epsilon;

epsilon_high_tension = [];
epsilon_high_tension = [epsilon_high_tension; epsilon_x; epsilon_y];

Ndrifters = length(drifters.x);
rejectedPointIndices_x = cell(Ndrifters,1);
rejectedPointIndices_y = cell(Ndrifters,1);
iDrifter = choiceDrifter;

% % Now we set a threshold for what constitutes an outlier. In this case we
% % choose points that have 1 in 10000 odds of occurring.
if IndependentXYOutlierRejection == 1
    rejectedPointIndices_x{iDrifter} = find(epsilon_x(2:end-1) < range(1) | epsilon_x(2:end-1) > range(2) );
    rejectedPointIndices_y{iDrifter} = find(epsilon_y(2:end-1) < range(1) | epsilon_y(2:end-1) > range(2) );
else
    %     fprintf('Generating the 2D student t-distribution...\n')
%     [r, pdf1d] = TwoDimStudentTProbabilityDistributionFunction( sigma, nu, 150, 3001 );
%     cdf_2d = cumtrapz(r,pdf1d);
%     r(end+1)=20000; cdf_2d(end+1) = 1.0; % so the interpolation algorithm has something to hang its hat on.
%     outlierCut = interp1(cdf_2d,r, 0.9999);
%     fprintf('outlierCut=%g\n',outlierCut);
    
    outlierCut=131.966;
    
    dist_error = sqrt( epsilon_x.*epsilon_x + epsilon_y.*epsilon_y );
    rejectedPointIndices_x{iDrifter} = find(dist_error(2:end-1) >= outlierCut);
    rejectedPointIndices_y{iDrifter} = rejectedPointIndices_x{iDrifter};
end

if ~isempty(rejectedPointIndices_x{iDrifter})
    rejectedPointIndices_x{iDrifter} = rejectedPointIndices_x{iDrifter}+1;
end
if ~isempty(rejectedPointIndices_y{iDrifter})
    rejectedPointIndices_y{iDrifter} = rejectedPointIndices_y{iDrifter}+1;
end

t_x = t; t_x(rejectedPointIndices_x{iDrifter}) = [];
x_reduced = x; x_reduced(rejectedPointIndices_x{iDrifter}) = [];
t_y = t; t_y(rejectedPointIndices_y{iDrifter}) = [];
y_reduced = y; y_reduced(rejectedPointIndices_y{iDrifter}) = [];

spline_x = SmoothingSpline(t_x,x_reduced,sqrt(variance_of_the_noise), 'S', S, 'T', T,'weightFunction',w);
spline_y = SmoothingSpline(t_y,y_reduced,sqrt(variance_of_the_noise), 'S', S, 'T', T,'weightFunction',w);

spline_x.Minimize( @(spline) spline.ExpectedMeanSquareError );
spline_y.Minimize( @(spline) spline.ExpectedMeanSquareError );

% Create our high resolution temporal grid
res = 5*60;
tq = res*( ceil(min(t_x(1),t_y(2))/res):1:floor(max(t_x(end),t_y(end))/res) )';
if min(t) < min(tq)
    tq = [min(t); tq];
end
if max(t) > max(tq)
    tq(end+1) = max(t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
spline_x_gaussian = SmoothingSpline(t,x,10,'S',S,'lambda',spline_x.lambda);
SmoothingSpline.MinimizeExpectedMeanSquareError(spline_x_gaussian);
spline_x_studentt = SmoothingSpline(t,x,sqrt(variance_of_the_noise),'S',S,'weightFunction',w,'lambda',spline_x.lambda);
SmoothingSpline.MinimizeExpectedMeanSquareError(spline_x_studentt);

spline_y_gaussian = SmoothingSpline(t,y,10,'S',S,'lambda',spline_x.lambda);
SmoothingSpline.MinimizeExpectedMeanSquareError(spline_y_gaussian);
spline_y_studentt = SmoothingSpline(t,y,sqrt(variance_of_the_noise),'S',S,'weightFunction',w,'lambda',spline_x.lambda);
SmoothingSpline.MinimizeExpectedMeanSquareError(spline_y_studentt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position fit figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
s = 1/1000; % scale
plot(s*spline_x(tq),s*spline_y(tq), 'LineWidth', 0.5*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]), hold on
scatter(s*drifters.x{iDrifter},s*drifters.y{iDrifter},(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')

rejects_xy = intersect(rejectedPointIndices_x{iDrifter},rejectedPointIndices_y{iDrifter});
rejects_x = setdiff(rejectedPointIndices_x{iDrifter},rejects_xy);
rejects_y = setdiff(rejectedPointIndices_y{iDrifter},rejects_xy);
scatter(s*drifters.x{iDrifter}(rejects_xy),s*drifters.y{iDrifter}(rejects_xy),(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
scatter(s*drifters.x{iDrifter}(rejects_x),s*drifters.y{iDrifter}(rejects_x),(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g')
scatter(s*drifters.x{iDrifter}(rejects_y),s*drifters.y{iDrifter}(rejects_y),(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position fit figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigureSize = [50 50 figure_width_2col+8 150*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

s = 1/1000; % scale
plot(tq/3600,s*spline_x(tq), 'LineWidth', 1.0*scaleFactor, 'Color',0.0*[1.0 1.0 1.0]), hold on
plot(tq/3600,s*spline_x_gaussian(tq), 'LineWidth', 0.5*scaleFactor, 'Color',0.5*[1.0 1.0 1.0])
% plot(tq/3600,s*spline_x_studentt(tq), 'LineWidth', 2.5*scaleFactor, 'Color',0.4*[0.0 1.0 0.0]), hold on
scatter(drifters.t{choiceDrifter}(rejectedPointIndices_x{choiceDrifter})/3600,s*drifters.x{choiceDrifter}(rejectedPointIndices_x{choiceDrifter}),(6.5*scaleFactor)^2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(drifters.t{choiceDrifter}/3600,s*drifters.x{choiceDrifter},(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('t (hours)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)

xlim([124 149])
ylim([5.58 9.18])

fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';

if shouldSaveFigures == 1
print('-depsc2', 'figures/tdistributionfit.eps')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position fit figure - Y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigureSize = [50 50 figure_width_2col+8 150*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

s = 1/1000; % scale
plot(tq/3600,s*spline_y(tq), 'LineWidth', 1.0*scaleFactor, 'Color',0.0*[1.0 1.0 1.0]), hold on
plot(tq/3600,s*spline_y_gaussian(tq), 'LineWidth', 0.5*scaleFactor, 'Color',0.5*[1.0 1.0 1.0])
% plot(tq/3600,s*spline_x_studentt(tq), 'LineWidth', 2.5*scaleFactor, 'Color',0.4*[0.0 1.0 0.0]), hold on
scatter(drifters.t{choiceDrifter}(rejectedPointIndices_y{choiceDrifter})/3600,s*drifters.y{choiceDrifter}(rejectedPointIndices_y{choiceDrifter}),(6.5*scaleFactor)^2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(drifters.t{choiceDrifter}/3600,s*drifters.y{choiceDrifter},(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('t (hours)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)

xlim([124 149])
% ylim([5.58 9.18])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Error distribution figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigureSize = [50 50 figure_width_1col+4 150*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error PDF plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = epsilon_high_tension;
histwidth = 80;
nbins = 100;

% Create the bins for the data
binwidth = histwidth/nbins;
edges = [-histwidth*100;((-histwidth/2+binwidth):binwidth:(histwidth/2-binwidth))';histwidth*100];
binleft = linspace((-histwidth/2),(histwidth/2-binwidth),nbins)';

% plot the data

% this is the data that doesn't
count = histcounts(data,edges)';
g = bar(binleft, count/(length(data)*binwidth), 'histc'); hold on;
g.FaceColor = 0.8*[1.0 1.0 1.0];

% create bins for the analytical pdf
xi_left = linspace(-histwidth/2,-histwidth/2+binwidth,10)';
xi_mid = linspace(-histwidth/2+binwidth,histwidth/2-binwidth,100)';
xi_right = linspace(histwidth/2-binwidth,histwidth/2,10)';
xi = [xi_left;xi_mid;xi_right];

% plot the analytical pdf
pdf = t_pdf;
edgedensity = integral(pdf,(histwidth/2-binwidth),2*histwidth)/binwidth;
denfunc = edgedensity*ones(size(xi));
denfunc(11:110) = pdf(xi_mid);
plot(xi,denfunc,'LineWidth',2*scaleFactor,'Color',0.0*[1.0 1.0 1.0])

xlabel('meters', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
set( gca, 'FontSize', figure_axis_tick_size);
xlim([min(xi) max(xi)])
ylim([0 0.175])
