function MakeFigureGPSErrorDistribution

% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

%  load('sample_data/motionless_garmin_epix.mat')
load('sample_data/motionless_garmin_edge_705_5_days.mat')

% The GPS was motionless, so its position is the errors, but we should
% remove the mean.
x=x-mean(x);
y=y-mean(y);
errors = [x;y];

% Now we find the optimal parameters for a student t-distribution
sortedErrors = sort([x;y]);
n = length(sortedErrors);
sigma_g = std(sortedErrors);
mu = mean(sortedErrors);

if 0 % hard code this in, because we already know the answer
    errorFunction = @(params) MotionlessDataMaximizationFunction( params, sortedErrors, mu);
    [optParams, lambda_t] = fminsearch( errorFunction, [sigma_g 10], optimset('TolX',1e-1, 'MaxFunEvals', 1e5, 'MaxFunEvals', 1e5) );
    sigma_t = optParams(1);
    nu = optParams(2);
    fprintf('Optimal fit parameters for the t-distribution, sigma_t=%f, nu=%f\n',sigma_t,nu)
else
    sigma_t=8.179396; nu=4.372602; % optimal for the Garmin Epix data set
    sigma_t=8.778964, nu=4.558391; % optimal for teh Garmin Edge data set
end


gaussian_pdf = @(z) exp(-(z.*z)/(2*sigma_g*sigma_g))/(sigma_g*sqrt(2*pi));
gaussian_cdf = @(z) 0.5*(1 + erf((z-mu)/(sigma_g*sqrt(2))));
rayleigh_pdf = @(z) (z/sigma_g^2) .* exp( -z.*z/(2*sigma_g^2));

t_pdf = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma_t*gamma(nu/2)*(1+(z.*z)/(nu*sigma_t*sigma_t)).^((nu+1)/2));
t_cdf = @(z) tcdf(z/sigma_t,nu);

if 0
    % used 20001 points for figure, but it took ~3 hours or so.
    [r, pdf1d] = TwoDimStudentTProbabilityDistributionFunction( sigma_t, nu, 100, 2001 );
    r(end+1)=1000; pdf1d(end+1) = 0.0;
    t2d_pdf = @(z) interp1(r,pdf1d,z);
else
    r=linspace(0,100,1000);
    pdf1d = tdistpdf(r/sigma_t,nu)/sigma_t;
    r(end+1)=1000; pdf1d(end+1) = 0.0;
    t2d_pdf = @(z) interp1(r,pdf1d,z);
end


cdf_2d = cumtrapz(r,pdf1d);
t2d_cdf = @(z) interp1(r,cdf_2d,z);

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
pdf = gaussian_pdf;
edgedensity = integral(pdf,(histwidth/2-binwidth),2*histwidth)/binwidth;
denfunc = edgedensity*ones(size(xi));
denfunc(11:110) = pdf(xi_mid);
plot(xi,denfunc,'LineWidth',2*scaleFactor,'Color',0.4*[1.0 1.0 1.0])

% plot the analytical pdf
pdf = t_pdf;
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
pdf = rayleigh_pdf;
edgedensity = integral(pdf,(histwidth-binwidth),20*histwidth)/binwidth;
denfunc = edgedensity*ones(size(xi));
denfunc(1:100) = pdf(xi_mid);
plot(xi,denfunc,'LineWidth',2*scaleFactor,'Color',0.4*[1.0 1.0 1.0])

% plot the analytical pdf
pdf = t2d_pdf;
edgedensity = integral(pdf,(histwidth-binwidth),20*histwidth)/binwidth;
denfunc = edgedensity*ones(size(xi));
denfunc(1:100) = pdf(xi_mid);
plot(xi,denfunc,'LineWidth',2*scaleFactor,'Color',0.0*[1.0 1.0 1.0])

v95 = interp1(t2d_cdf(r),r, 0.95);
plot([v95 v95],get(gca,'ylim'), 'LineWidth', 0.5*scaleFactor, 'Color', 0.4*[1.0 1.0 1.0]);

xlabel('meters', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
set( gca, 'FontSize', figure_axis_tick_size);

fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';

print('-depsc2', '../figures/gps_distance_error_distribution.eps')

return

stride = 10;
x_out = x(1:stride:end);
y_out = y(1:stride:end);

maxT = 1*3600;
dt = (0:stride:maxT)';
n = (0:(length(dt)-1))';
ACx = Autocorrelation(x_out, length(dt)-1);
ACy = Autocorrelation(y_out, length(dt)-1);

AC_95Confidence = (-1 + 1.645*sqrt(length(x_out)-n-1))./(length(x_out)-n);
SE =  sqrt((1 + 2*cumsum(ACx.^2))/length(x_out));
SE(1) = []; % first point is zero lag

figure
plot(dt,[ACx, ACy])
hold on
plot(dt, AC_95Confidence);
plot(dt(3:end), 2*SE(1:(end-1)), 'LineWidth', 1.5 )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error CDF plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2)
plot(sortedErrors,(1:n)'/n, 'k', 'LineWidth', 2)
hold on
z = linspace(-plotmax,plotmax,100)';
xlim([-plotmax plotmax])
plot(z, gaussian_cdf(z), 'LineWidth',1*scaleFactor,'Color',0.4*[1.0 1.0 1.0])
plot(z, tcdf((z-mu)/sigma_t,nu), 'LineWidth',1*scaleFactor,'Color',0.0*[1.0 1.0 1.0])

set(gca, 'XTick', []);
set( gca, 'FontSize', figure_axis_tick_size*scaleFactor);



end

% params should be [sigma nu], this is the KS test.
function lambda = MotionlessDataMaximizationFunction( params, sortedErrors, mu )

sigma = params(1);
nu = params(2);

x = sortedErrors;
n = length(x);
y_data = (1:n)'/n;
y = tcdf((x-mu)/sigma,nu);

D = max(abs(y-y_data));

lambda = (sqrt(n) + 0.12 + 0.11/sqrt(n))*D;
% j=1:25;
% p = sum(2*((-1).^(j-1)).*exp(-2*j.*j*lambda*lambda))

end



