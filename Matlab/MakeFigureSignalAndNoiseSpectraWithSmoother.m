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
numDerivs = 1;
S = 3;
T = 3;
K = S+1;

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

shouldUseObservedSignalOnly = 0;
result_stride = [1; 10; 100];
% result_stride = 100;
dof = zeros(length(result_stride),1);

FigureSize = [50 50 figure_width_large+7 300*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the signal and noise separately

D = SmoothingSpline.FiniteDifferenceMatrixNoBoundary(numDerivs,t,1);

dt = t(2)-t(1);
cv = D*(x + sqrt(-1)*y);
cepsilon = D*(epsilon_x + sqrt(-1)*epsilon_y);
[psi,lambda]=sleptap(size(cv,1));
[f,spp,snn,spn]=mspec(dt,cv,psi, 'cyclic');
[f,spp_e,snn_e,spn_e]=mspec(dt,cepsilon,psi,'cyclic');

ylimit = [1e-4 4e2];

sp1 = subplot(2,1,1);
sp2 = subplot(2,1,2);

set(fig1, 'currentaxes', sp1)
plot(f*timescale,vmean([snn, spp],2), 'k', 'LineWidth', 2);
hold on
plot(f*timescale,vmean([snn_e, spp_e],2), 'r', 'LineWidth', 2);
xlog, ylog
xlim([min(f*timescale) max(f*timescale)])
ylim(ylimit)

subplot(sp2)
plotHandles(1) = plot(1, nan, 'k', 'LineWidth', 2); hold on;
plotHandles(2) = plot(1, nan, 'r', 'LineWidth', 2);

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
        
    spline_x = SmoothingSpline(t_obs,x_obs,NormalDistribution(sigma),'S', S, 'T', T, 'knot_dof', 'auto');
    spline_y = SmoothingSpline(t_obs,y_obs,NormalDistribution(sigma),'S', S, 'T', T, 'knot_dof', 'auto');

    spline_x.minimizeMeanSquareError(data.t,data.x);
    spline_y.minimizeMeanSquareError(data.t,data.y);
 
    x_smooth = spline_x(t(indicesAll));
    y_smooth = spline_y(t(indicesAll));
    
    rms_error_x =  sqrt(mean(mean(  (data.x(indicesAll) - x_smooth).^2,2 ),1));
    rms_error_y =  sqrt(mean(mean(  (data.y(indicesAll) - y_smooth).^2,2 ),1));
    fprintf('rms error (x,y)=(%f,%f)\n',rms_error_x, rms_error_y);
    
    D = SmoothingSpline.FiniteDifferenceMatrixNoBoundary(numDerivs,t(indicesAll),1);
    u_smooth = D*x_smooth;
    v_smooth = D*y_smooth;
    
    dof(i) = (spline_x.effectiveSampleSizeFromVarianceOfTheMean + spline_y.effectiveSampleSizeFromVarianceOfTheMean)/2;
    f_limit = 1/(2*dof(i)*(t_obs(2)-t_obs(1)));
    
    
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
    
    plotHandles(2+i) = plot(f*timescale,RunningAverage(gamma,20)); hold on
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
mylegend = legend(plotHandles,'signal', 'noise', sprintf('stride 1 (n^{SE}_{eff}=%.2f)',dof(1)), sprintf('stride 10 (n^{SE}_{eff}=%.2f)',dof(2)), sprintf('stride 100 (n^{SE}_{eff}=%.2f)',dof(3)),'Location','southwest');
xlabel('cycles per minute', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

packfig(2,1)

print('-depsc2', sprintf('../figures/synthetic_process_and_spectrum_slope%ddegree%d.eps',abs(slope),S))
return

