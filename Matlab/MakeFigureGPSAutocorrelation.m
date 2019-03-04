% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

% Striding the data is a good way to check that our confidence intervals
% are robust.
stride = 1;
maxT = 1.0*3600;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start by using the epix data...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('sample_data/motionless_garmin_epix.mat')
x=x-mean(x);
y=y-mean(y);

x_out = x(1:stride:end);
y_out = y(1:stride:end);

dt = (0:stride:maxT)';
n = (0:(length(dt)-1))';
[ACx, DOFx] = Autocorrelation(x_out, length(dt)-1);
[ACy, DOFy] = Autocorrelation(y_out, length(dt)-1);

AC_epix = (ACx + ACy)/2;
DOF_epix = DOFx + DOFy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ...then using the edge data...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('sample_data/motionless_garmin_edge_705_5_days.mat')
% load('sample_data/motionless_garmin_edge_705.mat')

x=x-mean(x);
y=y-mean(y);

x_out = x(1:stride:end);
y_out = y(1:stride:end);

dt = (0:stride:maxT)';
n = (0:(length(dt)-1))';
[ACx, DOFx] = Autocorrelation(x_out, length(dt)-1);
[ACy, DOFy] = Autocorrelation(y_out, length(dt)-1);

AC_edge = (ACx + ACy)/2;
DOF_edge = DOFx + DOFy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ...then using the edge data...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('sample_data/motionless_garmin_edge_705.mat')

x=x-mean(x);
y=y-mean(y);

x_out = x(1:stride:end);
y_out = y(1:stride:end);

dt = (0:stride:maxT)';
n = (0:(length(dt)-1))';
[ACx, DOFx] = Autocorrelation(x_out, length(dt)-1);
[ACy, DOFy] = Autocorrelation(y_out, length(dt)-1);

AC_edge2 = (ACx + ACy)/2;
DOF_edge2 = DOFx + DOFy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ...and finally merge them together with a weighted average.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AC = (DOF_epix(1)*AC_epix + DOF_edge(1)*AC_edge + DOF_edge2(1)*AC_edge2)/(DOF_epix(1) + DOF_edge(1) +  DOF_edge2(1));
DOF = DOF_epix + DOF_edge + DOF_edge2;

SE_indep = dt(2:end);
SE =  sqrt((1 + 2*cumsum(AC.^2))./DOF);
SE(1) = sqrt(1/DOF(1)); % first point is lag 1
SE(end) = []; % there is no end point

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Autocorrelation figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigureSize = [50 50 figure_width_1col+4 150*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

s = 1/60;

plot(s*dt,AC, 'LineWidth',1*scaleFactor,'Color',0.0*[1.0 1.0 1.0])
hold on
plot(s*SE_indep, [3*SE,-3*SE], 'LineWidth', 1.5, 'Color',0.4*[1.0 1.0 1.0] )
xlabel('time lag (minutes)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('autocorrelation', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlim([0 s*maxT])
ylim([-0.2 1.0])

rho = @(dt) exp(max(-abs(dt)/100., - abs(dt)/760 -1.3415)) ;
hold on
plot(s*dt,rho(dt))

% 
% [p,S,mu] = polyfit(dt(1:1200),log(AC(1:1200)),5);
% hold on, plot(s*dt,exp(polyval(p,dt,[],mu)));

% print('-depsc2', '../figures/gps_autocorrelation.eps')

return

range = 1:find(dt < 22*60,1,'last');
[p,~,mu]=polyfit(dt(range),log(AC(range)),2);
slope = p(1)/mu(2)
figure, plot(dt(range),log(AC(range))), hold on, plot(dt(range),polyval(p,dt(range),[],mu))

figure, plot(dt(range),log(AC(range))), hold on,  plot( -(1/(5*60))*dt(range))

[TAU,R]=materncov(1,1500,1,0.6,1/(10*60));
R = R/R(1);
figure, plot(dt(range),log(AC(range))), hold on, plot(TAU,log(R))

figure, plot(dt(range),log(AC(range))), hold on, plot( (p(1)/mu(2))*dt(range) + p(2)-p(1)*mu(1)/mu(2))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cross correlation figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('sample_data/motionless_garmin_epix.mat')
x=x-mean(x);
y=y-mean(y);

x_out = x(1:stride:end);
y_out = y(1:stride:end);

dt = (0:stride:maxT)';
n = (0:(length(dt)-1))';

[C_epix, DOF_epix] = Crosscorrelation(x_out,y_out,length(dt)-1);

load('sample_data/motionless_garmin_edge_705.mat')
x=x-mean(x);
y=y-mean(y);

x_out = x(1:stride:end);
y_out = y(1:stride:end);

dt = (0:stride:maxT)';
n = (0:(length(dt)-1))';

[C_edge, DOF_edge] = Crosscorrelation(x_out,y_out,length(dt)-1);

C = (C_epix+C_edge)/2;
DOF = DOF_epix+DOF_edge;

SE_indep = dt(2:end);
SE =  sqrt((1 + 2*cumsum(C.^2))./DOF);
SE(1) = sqrt(1/DOF(1)); % first point is lag 1
SE(end) = []; % there is no end point

figure
s = 1/60;

plot(s*dt,C, 'LineWidth',1*scaleFactor,'Color',0.0*[1.0 1.0 1.0])
hold on
plot(s*SE_indep, [2*SE,-2*SE], 'LineWidth', 1.5, 'Color',0.4*[1.0 1.0 1.0] )
xlabel('lag (minutes)')
ylabel('cross correlation')