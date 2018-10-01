% Setting the order of the spline sets the 'interpolation fall off', by
% which I mean that it sets the slope that the recovered signal falls off
% at, once it becomes garbage.
%
% I think this falloff is with a slope *at least* as steep as the actual
% signal slope for the lower frequency garbage.

addpath('./support');
% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

slope = -2;

if slope == -2
    data = load('sample_data/SyntheticTrajectories.mat');
elseif slope == -3
    data = load('sample_data/SyntheticTrajectoriesSlope3.mat');
elseif slope == -4
    data = load('sample_data/SyntheticTrajectoriesSlope4.mat');
end

t = data.t;
x = data.x;
y = data.y;
epsilon_x = data.epsilon_x;
epsilon_y = data.epsilon_y;
position_error = data.position_error;

timescale = 60;
indices = find( t/timescale >= 0 & t/timescale <= 42);

markersize  = 2*scaleFactor;
stride = 10; % we start by decimating the signal
range = indices(1:stride:end);

FigureSize = [50 50 figure_width_large+7 300*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];



subplot(2,1,1)
errorbar(t(range)/timescale,y(range),position_error*ones(size(t(range))),'ko', 'LineWidth', 0.5*scaleFactor, 'MarkerSize', markersize^2, 'MarkerFaceColor', 'w')
hold on
xlim([min(t(range))/timescale max(t(range))/timescale]);
ylim([1.2*min(y(range)) 1.2*max(y(range))])

range = indices(1:100:end);
errorbar(t(range)/timescale,y(range),position_error*ones(size(t(range))),'ko', 'LineWidth', 1.5*scaleFactor, 'MarkerSize', markersize^2, 'MarkerFaceColor', 'k')

xlabel('minutes', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
ylabel('meters', 'FontSize', figure_axis_label_size, 'FontName', figure_font);



D = TensionSpline.FiniteDifferenceMatrixNoBoundary(1,t,1);

dt = t(2)-t(1);
cv = D*(x + sqrt(-1)*y);
cepsilon = D*(epsilon_x + sqrt(-1)*epsilon_y);
[psi,lambda]=sleptap(size(cv,1));
[f,spp,snn,spn]=mspec(dt,cv,psi, 'cyclic');
[f,spp_e,snn_e,spn_e]=mspec(dt,cepsilon,psi,'cyclic');

% sqrt(2)*position_error/dt
% std(imag(cepsilon))

ylimit = [1e-4 4e2];

subplot(2,1,2)
plot(f*timescale,vmean([snn, spp],2), 'k', 'LineWidth', 2)
hold on
plot(f*timescale,vmean([snn_e, spp_e],2), 'r', 'LineWidth', 2)
xlog, ylog
xlim([min(f*timescale) max(f*timescale)])
ylim(ylimit)

shouldUseObservedSignalOnly = 0;
result_stride = [1; 10; 100];
for i=1:length(result_stride)
    stride = result_stride(i);
    % Reduce the total length in some cases
    if (stride < 10)
        shortenFactor = stride/10;
        shortenFactor = 1;
    else
        shortenFactor = 1;
    end
    
    indices = 1:stride:floor(shortenFactor*length(data.t));
    fprintf('dT = %d --- Using %d points with stride %d\n', stride*dt, length(indices), stride);
    if (shouldUseObservedSignalOnly == 1)
        indicesAll = indices;
    else
        indicesAll = 1:max(indices);
    end
    epsilon_x = data.epsilon_x(indices);
    epsilon_y = data.epsilon_y(indices);
    x_obs = data.x(indices) + epsilon_x;
    y_obs = data.y(indices) + epsilon_y;
    t_obs = data.t(indices);
    sigma = data.position_error;
    
    S = 1;
    T = 1;
    K = S+1;
    
    spline_x = TensionSpline(t_obs,x_obs,sigma,'S', S, 'T', T);
    spline_y = TensionSpline(t_obs,y_obs,sigma,'S', S, 'T', T);
    TensionSpline.MinimizeExpectedMeanSquareError(spline_x);
    TensionSpline.MinimizeExpectedMeanSquareError(spline_y);
    
    cv = D*(spline_x(t)+sqrt(-1)*spline_y(t));
    [f,spp,snn,spn]=mspec(dt,cv,psi,'cyclic');
    plot(f*timescale,vmean([snn, spp],2), 'LineWidth', 0.5*scaleFactor)
    
    
end

legend('signal', 'noise', '1', '10', '100');

f_smooth = 0.5/(t_obs(2)-t_obs(1));
    plot([f_smooth f_smooth]*timescale,ylimit, 'LineWidth', 0.5*scaleFactor, 'Color', 0.4*[1.0 1.0 1.0]);


% dt_gamma1 = 3*position_error/sigma_u;
% f_gamma1 = timescale/dt_gamma1;
% % plot([f_gamma1 f_gamma1],ylimit, 'LineWidth', 0.5*scaleFactor, 'Color', 0.4*[1.0 1.0 1.0]);
% % 
% dt_gamma10 = 3*position_error/(sigma_u*10);
% f_gamma10 = timescale/dt_gamma10;
% % plot([f_gamma10 f_gamma10],ylimit, 'LineWidth', 0.5*scaleFactor, 'Color', 0.4*[1.0 1.0 1.0]);
% % 
% dt_gamma01 = 3*position_error/(sigma_u*0.1);
% f_gamma01 = timescale/dt_gamma01;
% % plot([f_gamma01 f_gamma01],ylimit, 'LineWidth', 0.5*scaleFactor, 'Color', 0.4*[1.0 1.0 1.0]);

xlabel('cycles per minute', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
ylabel('power (m^2/s)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

print('-depsc2', '../figures/synthetic_process_and_spectrum_with recovery.eps')

dof1 = 1 + 3*position_error/(sigma_u*dt);
dof10 = 1 + 3*position_error/(sigma_u*10*dt);
dof100 = 1 + 3*position_error/(sigma_u*100*dt);
