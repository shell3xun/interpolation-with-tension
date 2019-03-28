% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

% data = load('sample_data/motionless_garmin_epix.mat');
gpsdata = load('sample_data/motionless_garmin_edge_705_5_days.mat');

x = gpsdata.x;
y = gpsdata.y;
t = gpsdata.t;

% The GPS was motionless, so its position is the errors, but we should
% remove the mean.
x=x-mean(x);
y=y-mean(y);
errors = [x;y];

% Now we find the optimal parameters for a student t-distribution
sortedErrors = sort([x;y]);
sigma_g = std(sortedErrors);

% optParams = fminsearch( @(p) StudentTDistribution(p(1),p(2)).kolmogorovSmirnovError(sortedErrors) , [sigma_g 5], optimset('TolX',1e-1) );
% sigma_t = optParams(1);
% nu = optParams(2);
% fprintf('Optimal fit parameters for the t-distribution using KS, sigma_t=%f, nu=%f\n',sigma_t,nu)

optParams = fminsearch( @(p) StudentTDistribution(p(1),p(2)).andersonDarlingError(sortedErrors) , [sigma_g 5], optimset('TolX',1e-2) );
sigma_t = optParams(1);
nu = optParams(2);
fprintf('Optimal fit parameters for the t-distribution using AD, sigma_t=%f, nu=%f\n',sigma_t,nu)

% Create the 1D and 2D distributions of the errors
normalDistribution = NormalDistribution(sigma_g);
tDistribution = StudentTDistribution(sigma_t,nu);

normalDistanceDistribution =  RayleighDistribution(sigma_g);
tDistanceDistribution = TwoDimDistanceDistribution(tDistribution);

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

data = errors;
histwidth = 80;
nbins = 100;

% Create the bins for the data
binwidth = histwidth/nbins;
edges = [-histwidth*100;((-histwidth/2+binwidth):binwidth:(histwidth/2-binwidth))';histwidth*100];
binleft = linspace((-histwidth/2),(histwidth/2-binwidth),nbins)';

% plot the data
count = histcounts(data,edges)';
g = bar(binleft, count/(length(data)*binwidth), 'histc'); hold on;
g.FaceColor = 0.8*[1.0 1.0 1.0];

% create bins for the analytical pdf
xi_left = linspace(-histwidth/2,-histwidth/2+binwidth,10)';
xi_mid = linspace(-histwidth/2+binwidth,histwidth/2-binwidth,100)';
xi_right = linspace(histwidth/2-binwidth,histwidth/2,10)';
xi = [xi_left;xi_mid;xi_right];

% plot the analytical pdf
pdf = normalDistribution.pdf;
edgedensity = integral(pdf,(histwidth/2-binwidth),2*histwidth)/binwidth;
denfunc = edgedensity*ones(size(xi));
denfunc(11:110) = pdf(xi_mid);
plot(xi,denfunc,'LineWidth',2*scaleFactor,'Color',0.4*[1.0 1.0 1.0])

% plot the analytical pdf
pdf = tDistribution.pdf;
edgedensity = integral(pdf,(histwidth/2-binwidth),2*histwidth)/binwidth;
denfunc = edgedensity*ones(size(xi));
denfunc(11:110) = pdf(xi_mid);
plot(xi,denfunc,'LineWidth',2*scaleFactor,'Color',0.0*[1.0 1.0 1.0])

xlabel('meters', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
set( gca, 'FontSize', figure_axis_tick_size);
xlim([min(xi) max(xi)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tighten up the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%packfig(2,1)
fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';

print('-depsc2', '../figures/gps_error_distribution.eps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Error distance distribution figures
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

data =  sqrt(x.*x + y.*y);
histwidth = 50;
nbins = 100;

% Create the bins for the data
binwidth = histwidth/nbins;
edges = [(0:binwidth:(histwidth-binwidth))';histwidth*10];
binleft = edges(1:end-1);

% plot the data
count = histcounts(data,edges)';
g = bar(binleft, count/(length(data)*binwidth), 'histc'); hold on;
g.FaceColor = 0.8*[1.0 1.0 1.0];

% create bins for the analytical pdf
xi_mid = linspace(0,histwidth-binwidth,100)';
xi_right = linspace(histwidth-binwidth,histwidth,10)';
xi = [xi_mid;xi_right];

% plot the analytical pdf
pdf = normalDistanceDistribution.pdf;
edgedensity = integral(pdf,(histwidth-binwidth),20*histwidth)/binwidth;
denfunc = edgedensity*ones(size(xi));
denfunc(1:100) = pdf(xi_mid);
plot(xi,denfunc,'LineWidth',2*scaleFactor,'Color',0.4*[1.0 1.0 1.0])

% plot the analytical pdf
pdf = tDistanceDistribution.pdf;
edgedensity = integral(pdf,(histwidth-binwidth),20*histwidth)/binwidth;
denfunc = edgedensity*ones(size(xi));
denfunc(1:100) = pdf(xi_mid);
plot(xi,denfunc,'LineWidth',2*scaleFactor,'Color',0.0*[1.0 1.0 1.0])

v95 = tDistanceDistribution.locationOfCDFPercentile(0.95);
plot([v95 v95],get(gca,'ylim'), 'LineWidth', 0.5*scaleFactor, 'Color', 0.4*[1.0 1.0 1.0]);

xlabel('meters', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
set( gca, 'FontSize', figure_axis_tick_size);

fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';

print('-depsc2', '../figures/gps_distance_error_distribution.eps')



