% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('./support')

% This figure will create an interpolating spline (of increasing order)
% using a sampled sine function. The points chosen for the sine function
% here are artibrary, and simply meant to make it not obvious what the
% underlying function is.
t = 2*pi*2.5*[0;1;3;4;5;8;10]/10;
x = sin(t);

% tq is a high density data grid, simply used for plotting the
% interpolating function.
tq = linspace(t(1),t(end),100*length(t))';
ylimits = [-1.5 1.5];

maxK = 4;
maxD = 4;
iSpline = 2;
iSubplot = 1;
h = zeros(maxK*maxD,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position Figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigureSize = [50 50 figure_width_2col+8 225*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize, 'Name', 'BSpline interpolation example');
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

for K=1:maxK
    % All the interesting work takes place with this function call, which
    % returns an interpolating spline of order K.
    x_spline = InterpolatingSpline(t,x,K);
    
    t_knot = x_spline.t_knot;
    
    for D=1:maxK
        index=(K-1)*maxD+D;
        if D > min(K,maxD)
            % We create an empty subplot and set its visibility to zero,
            % otherwise packfig (below) creates a whitebox.
            h(index) = subplot(maxK,maxD,index);
            set(h(index),'visible','off')
        else
            h(index) = subplot(maxK,maxD,iSubplot);
            
            for i=1:length(t_knot)
                plot([t_knot(i) t_knot(i)],ylimits, 'LineWidth', 0.5*scaleFactor, 'Color', 0.4*[1.0 1.0 1.0]); hold on;
            end
            ylim(ylimits)
            xlim([min(t) max(t)])
            
            % here we plot the D-1th derivative of the interpolant.
            plot(tq,x_spline(tq,D-1), 'LineWidth', 1.0*scaleFactor, 'Color', 'k')
            
            if (D==1)
                scatter(t,x,(5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
            end
            
         
            set( gca, 'FontSize', figure_axis_tick_size);
            
            set(gca, 'XTick', []);
            if (K<maxK)
                set(gca, 'XTick', []);
            end
            
            set(gca, 'YTick', []);
            if (D>1)
                set(gca, 'YTick', []);
            else
                ylabel(sprintf('$K=%d$',K),'Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
            end
            
            if (D==K)
                if D==2
                    title(sprintf('$1^{st}$-deriv'),'Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
                elseif D==3
                    title(sprintf('$2^{nd}$-deriv'),'Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
                elseif D==4
                    title(sprintf('$3^{rd}$-deriv'),'Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
                end
                iSubplot = K*maxD + 1;
            else
                iSubplot = iSubplot + 1;
            end
        end
    end
end

packfig(maxK,maxD)
fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';

print('-depsc2', '../figures/interpolation.eps')

