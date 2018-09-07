% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

t_data = [0;1;3;4;5;8;10];

t = linspace(0,10,1001)';

maxK = 4;
maxD = 4;
iSpline = 2;
iSubplot = 1;
h = zeros(maxK*maxD,1);
ylimit = [-2 2];
FigureSize = [50 50 figure_width_1col+7 225*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize, 'Name', 'BSpline example');
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

for K=1:maxK
    
    % increase the multiplicity of the end knots for higher order splines
    t_knot = [repmat(t_knot(1),K-1,1); t_data; repmat(t_knot(end),K-1,1)];
    
    % And then let's plot one of the first splines that spans all the knots
    m = zeros(length(t_data)+K-2,1);
    m(iSpline+K-1) = 1;
    
    B = BSpline(t_data,zeros(size(t_data)),K,t_knot,m);
    
    for D=1:maxK
        index=(K-1)*maxD+D;
        if D > min(K,maxD)
            % We create an empty subplot and set its visibility to zero,
            % otherwise packfig (below) creates a whitebox.
            h(index) = subplot(maxK,maxD,index);
            set(h(index),'visible','off')
        else
            h(index) = subplot(maxK,maxD,iSubplot);
            
            for i=1:length(t_data)
                plot([t_data(i) t_data(i)],ylimit, 'LineWidth', 0.5*scaleFactor, 'Color', 0.4*[1.0 1.0 1.0]); hold on
            end
            ylim(ylimit)
            plot(t,B(t,D-1), 'LineWidth', 1.0*scaleFactor, 'Color', 'k');
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

print('-depsc2', 'bsplines.eps')