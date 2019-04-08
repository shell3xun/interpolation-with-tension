% \Gamma = \sigma/(u*dt)
% The hypothesis is that d = 1+scale*\Gamma
%
% We find that the scale is larger for steep slopes, and smaller for
% shallower slopes. This makes intuitive sense, because the underlying
% function isn't changing as much more smoother functions (steeper slopes).
% So, your degrees of freedom increase more quickly for steeper slopes.

scaleFactor = 1;
LoadFigureDefaults

shouldLoadExistingResults = 1;

addpath('support')

rms = @(z) sqrt( mean( z.^2 ) );

filename = 'DegreesOfFreedomEstimates.mat';

if exist(filename,'file')
    load(filename);   
else
    error('You should create this file using MakeFigureGammaDependence');
end


gamma = sigma./(measured_u_rms.*measurement_time);
n_eff = 9*gamma.^0.69;
n_eff = 14*gamma.^0.73;
% n_eff = 18*gamma.^0.74;
n_eff(n_eff<1)=1;
lambda_est = (n_eff-1)./(n_eff.*measured_a_rms.^2);
% lambda_est = (n_eff-1)./(n_eff.*measured_a_rms_filtered.^2);

% figure
% subplot(2,1,1)
% scatter(reshape(measured_dof_se,[],1),reshape(n_eff,[],1))
% xlog, ylog
% title('n_{eff}')
% ylabel('actual')
% 
% subplot(2,1,2)
% scatter(reshape(measured_a_rms,[],1),reshape(actual_a_rms,[],1))
% xlog, ylog
% title('a_{rms}')
% ylabel('actual')

lambda_est_true = (measured_dof_se-1)./(measured_dof_se.*actual_a_rms.^2);
figure,
for iSlope=1:length(slopes)
%     scatter(reshape(measured_lambda(:,iSlope,:),[],1),reshape(lambda_est_true(:,iSlope,:),[],1)), hold on
measure = actual_ad_error(:,iSlope,:);
sz = max(1,round(abs(measure(:)))).^2;
max(sz)
scatter(reshape(measured_lambda(:,iSlope,:),[],1),reshape(lambda_est(:,iSlope,:),[],1),sz), hold on
end
hold on, plot(10.^linspace(12,15,100),10.^linspace(12,15,100),'k')
xlim([1e12 1e16])

% delta = log10(abs(lambda_est - measured_lambda));
delta = log10(lambda_est./measured_lambda);

sqrt(mean(delta(delta>-12).^2))
median(delta(delta>-12))

xlabel('true lambda')
ylabel('estimated lambda')
% hold on, plot(reshape(measured_lambda(:),[],1),reshape(measured_lambda(:),[],1))
xlog, ylog

return
figure,
scatter(reshape(measured_lambda,[],1),reshape(lambda_est,[],1))
xlog, ylog