% Setting the order of the spline sets the 'interpolation fall off', by
% which I mean that it sets the slope that the recovered signal falls off
% at, once it becomes garbage.
%
% I think this falloff is with a slope *at least* as steep as the actual
% signal slope for the lower frequency garbage.
%
% The -2 slope, needs S=1
% The -4 slope, needs S=2
% Anything else is either too steep or too shallow in the spectrum
% The -3 slope is what you'd expect: neither works well
%
% Coherence would test signal to recovered signal

addpath('./support');
% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

slope = -2;
timescale = 60;

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

FigureSize = [50 50 figure_width_large+7 300*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the signal and noise separately

D = SmoothingSpline.FiniteDifferenceMatrixNoBoundary(1,t,1);

dt = t(2)-t(1);
cv = D*(x + sqrt(-1)*y);
[psi,lambda]=sleptap(size(cv,1));
[f,spp,snn,spn]=mspec(dt,cv,psi, 'cyclic');

ylimit = [1e-4 4e2];

sp1 = subplot(2,1,1);
sp2 = subplot(2,1,2);

set(fig1, 'currentaxes', sp1)
plot(f*timescale,vmean([snn, spp],2), 'k', 'LineWidth', 2);
hold on
xlog, ylog
xlim([min(f*timescale) max(f*timescale)])
ylim(ylimit)

subplot(sp2)
plotHandles(1) = plot(1, nan, 'k', 'LineWidth', 2); hold on;

rms_error_x = zeros(3,1);
rms_error_y = zeros(3,1);
for iSlope=1:4
    stride = 100;

    indices = 1:stride:floor(length(data.t));
    fprintf('dT = %d --- Using %d points with stride %d\n', stride*dt, length(indices), stride);
    if (shouldUseObservedSignalOnly == 1)
        indicesAll = indices;
    else
        indicesAll = 1:max(indices);
    end
    x_obs = data.x(indices);
    y_obs = data.y(indices);
    t_obs = data.t(indices);
        
    spline_x = InterpolatingSpline(t_obs,x_obs,'S', iSlope);
    spline_y = InterpolatingSpline(t_obs,y_obs,'S', iSlope);

    x_smooth = spline_x(t(indicesAll));
    y_smooth = spline_y(t(indicesAll));
    
    rms_error_x(iSlope) =  sqrt(mean(mean(  (data.x(indicesAll) - x_smooth).^2,2 ),1));
    rms_error_y(iSlope) =  sqrt(mean(mean(  (data.y(indicesAll) - y_smooth).^2,2 ),1));
    fprintf('rms error (x,y)=(%f,%f)\n',rms_error_x(iSlope), rms_error_y(iSlope));
    
    D = SmoothingSpline.FiniteDifferenceMatrixNoBoundary(1,t(indicesAll),1);
    u_smooth = D*x_smooth;
    v_smooth = D*y_smooth;
    
    f_limit = 1/(2*(t_obs(2)-t_obs(1)));
    
    subplot(sp1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot the spline fit
    cv = u_smooth+sqrt(-1)*v_smooth;
    [psi,lambda]=sleptap(size(cv,1));
    [f,spp,snn,spn]=mspec(dt,cv,psi,'cyclic');
    ax = gca;
    CO = ax.ColorOrderIndex;
    plot(f*timescale,vmean([snn, spp],2), 'LineWidth', 0.5*scaleFactor);
    ylimits = ylim;
    plot(timescale*f_limit*[1 1],ylimits,'k--')
    
    
    subplot(sp2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot the coherence of the spline fit with the true signal
    [f,sxx,syy,sxy]=mspec(dt,u_smooth,D*data.x(indicesAll),psi,'cyclic');
    gamma=frac(abs(sxy).^2,sxx.*syy);
    
    plotHandles(1+iSlope) = plot(f*timescale,RunningAverage(gamma,20)); hold on
    xlog
    xlim([min(f*timescale) max(f*timescale)])
    ylim([0 1])
    
    ylimits = ylim;
    plot(timescale*f_limit*[1 1],ylimits,'k--')
end

subplot(sp1)
ylabel('power (m^2/s)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

subplot(sp2)
ylabel('coherence', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
xlabel('cycles per minute', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
mylegend = legend(plotHandles,'signal', sprintf('S=1 (%.2f m rmse)',rms_error_x(1)), sprintf('S=2 (%.2f m rmse)',rms_error_x(2)), sprintf('S=3 (%.2f m rmse)',rms_error_x(3)), sprintf('S=4 (%.2f m rmse)',rms_error_x(4)),'Location','southwest');
xlabel('cycles per minute', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

packfig(2,1)

print('-depsc2', sprintf('../figures/interpolation_spectrum_slope%ddegreeVaried.eps',abs(slope)))
return

