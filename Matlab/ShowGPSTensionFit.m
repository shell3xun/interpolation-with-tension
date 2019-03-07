% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

shouldSaveFigures = 0;

% Drifter to highlight in the final plots. Drifter 7 has no outliers
choiceDrifter = 6;

drifters = open('sample_data/raw_rho1_drifters.mat');

t_drifter = (drifters.date{choiceDrifter}-drifters.lastDeployment)*24*60*60;
spline = GPSTensionSpline(t_drifter,drifters.lat{choiceDrifter},drifters.lon{choiceDrifter},'shouldUseRobustFit',1);
% [outlierDistribution,alpha] = spline.setSigmaFromOutlierDistribution();

tq = linspace(min(t_drifter),max(t_drifter),10*length(t_drifter))';
[x,y] = spline.xyAtTime(tq);
figure
scatter(spline.x,spline.y,(2.5)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
plot(x,y),axis equal

figure
scatter(spline.q,spline.r,(2.5)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
plot(spline.spline_x(tq),spline.spline_y(tq))

figure
subplot(2,1,1)
scatter(spline.t,spline.spline_x.x,(2.5)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
plot(tq,spline.spline_x(tq))
subplot(2,1,2)
scatter(spline.t,spline.spline_y.x,(2.5)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
plot(tq,spline.spline_y(tq))

% q = x - mean(x);
% r = y - mean(y);
% M_qq = mean(q.*q);
% M_rr = mean(r.*r);
% M_qr = mean(q.*r);
% [A, B] = eig([M_qq, M_qr; M_qr, M_rr]);
% v = A(:,2);
% theta = atan(v(2)/v(1));
% minD = B(1,1);
% maxD = B(end,end);
% 
% theta = theta + pi/2;
% x_tilde = q*cos(theta) + r*sin(theta);
% y_tilde = -q*sin(theta) + r*cos(theta);



return
figure
subplot(2,2,[1 3])
plot(spline.q,spline.r), axis equal
subplot(2,2,2), plot(tq,x_tilde), subplot(2,2,4), plot(tq,y_tilde)
% figure
% plot(q,r), axis equal

return

% Pull out the data of interest
x_data = drifters.x{choiceDrifter};
y_data = drifters.y{choiceDrifter};
t_data = drifters.t{choiceDrifter};

tq = linspace(min(t_data),max(t_data),10*length(t_data));

% These are our working definitions for the noise
noiseDistribution = StudentTDistribution(8.5,4.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make a figure with the data
figure
sp1=subplot(2,1,1);
scatter(t_data,x_data,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
% errorbar(t_data,x_data,zmax*ones(size(t_data))), hold on
sp2=subplot(2,1,2);
scatter(t_data,y_data,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
% errorbar(t_data,y_data,zmax*ones(size(t_data))), hold on

spline_x = RobustTensionSpline(t_data,x_data,noiseDistribution,'lambda',Lambda.fullTensionIterated,'outlierMethod',OutlierMethod.sigmaMethod);
spline_y = RobustTensionSpline(t_data,y_data,noiseDistribution,'lambda',Lambda.fullTensionIterated,'outlierMethod',OutlierMethod.sigmaMethod);

return

%%%%%%%%%%%%%%%%%%%%%
% Grab drifter 7 and plot that
spline_7x = TensionSpline(drifters.t{7},drifters.x{7},noiseDistribution);
spline_7y = TensionSpline(drifters.t{7},drifters.y{7},noiseDistribution);

subplot(sp1)
scatter(t_data(spline_x.outlierIndices),x_data(spline_x.outlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(t_data(spline_x.outlierIndices),x_data(spline_x.outlierIndices),(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline_x(tq))
plot(tq,spline_7x(tq))


fprintf('robust a_rms: (%.3g, %.3g)\n', std(spline_x.uniqueValuesAtHighestDerivative), std(spline_y.uniqueValuesAtHighestDerivative) );
% fprintf('robust a_rms: (%.3g, %.3g)\n', TensionSpline.StandardDeviationAndMeanOfInterquartileRange(spline_x.uniqueValuesAtHighestDerivative), TensionSpline.StandardDeviationAndMeanOfInterquartileRange(spline_y.uniqueValuesAtHighestDerivative) );
fprintf('robust a_rms: (%.3g, %.3g)\n', std(spline_7x.uniqueValuesAtHighestDerivative), std(spline_7y.uniqueValuesAtHighestDerivative) );

subplot(sp2)
scatter(t_data(spline_y.outlierIndices),y_data(spline_y.outlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(t_data(spline_y.outlierIndices),y_data(spline_y.outlierIndices),(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline_y(tq))
plot(tq,spline_7y(tq))


return

[fn,sn] = powspec(tq(2)-tq(1),diff(spline_x2(tq)')/(tq(2)-tq(1)));
figure
plot(fn,sn), xlog, ylog
