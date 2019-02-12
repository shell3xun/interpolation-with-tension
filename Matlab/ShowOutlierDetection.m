scaleFactor = 1;
LoadFigureDefaults

shouldUseStudentTDistribution = 1;

percentOutliers = 0.25;
outlierFactor = 40;

slopes = [-2; -3; -4];
slope = -3;
totalSlopes = length(slopes);

S = 2;
T = S;
K = S+1;

result_stride = 2*[1;4;16;64];
result_stride = 200;
totalStrides = length(result_stride);

% matern signal parameters
sigma_u = 0.20;
base_dt = 5; % for whatever reason, we chose this as the primary dt
t_damp = 30*60;
n = 250;

totalEnsembles = 1;

for iSlope = 1:length(slopes)
    slope = slopes(iSlope);
    fprintf('slope %d, ',slope);
    
    for iStride=1:length(result_stride)
        stride = result_stride(iStride);
        dt = stride*base_dt;
        
        for iEnsemble = 1:totalEnsembles
            fprintf('..%d',iEnsemble);
            
            cv=maternoise(dt,n,sigma_u*sqrt(2),abs(slope),1/t_damp);
            cx = cumtrapz(cv)*dt;
            data = struct('t',dt*(0:n-1)','x',real(cx));
            
            if shouldUseStudentTDistribution == 1
                nu = 4.5; sigma =  8.5;        
                noiseDistribution = StudentTDistribution(sigma,nu);
                outlierDistribution = StudentTDistribution(outlierFactor*sigma,3.0);
            else
                sigma = 10;
                noiseDistribution = NormalDistribution(sigma);
                outlierDistribution = NormalDistribution(outlierFactor*sigma);
            end
            outlierIndices = rand(n,1)<=percentOutliers;
            outlierIndices([1 n]) = 0;
            
            epsilon_x = zeros(n,1);
            epsilon_x_outlier = zeros(n,1);
            
            epsilon_x(~outlierIndices) = noiseDistribution.rand(sum(~outlierIndices));
            epsilon_x_outlier(outlierIndices) = outlierDistribution.rand(sum(outlierIndices));
            
            epsilon = epsilon_x + epsilon_x_outlier;
            outlierThreshold = noiseDistribution.locationOfCDFPercentile(1-1/10000/2);
            trueOutliers = find(abs(epsilon) > outlierThreshold);
            goodIndices = setdiff(1:n,trueOutliers);
             
            x_obs = data.x + epsilon;
            t_obs = data.t;  
        end
        fprintf('\n');
    end
    
end
trueAddedDistribution = AddedDistribution(percentOutliers,outlierDistribution,noiseDistribution);
fprintf('sqrt(variance) (noise, outlier)=(%.1f, %.1f)\n',sqrt(noiseDistribution.variance),sqrt(outlierDistribution.variance));
tq = linspace(min(t_obs),max(t_obs),10*length(t_obs));
% [D,tD] = TensionSpline.FiniteDifferenceMatrixNoBoundary(2,data.t(indices),2);
% spline_clean = TensionSpline(data.t(indices),data.x(indices),noiseDistribution, 'S', S, 'lambda',0);
% 
% figure
% plot(tq,spline_clean(tq,2)), hold on
% plot(tD,D*data.x(indices))

dt = t_obs(2)-t_obs(1);
a = diff(diff(data.x))/(dt^2);
fprintf('true strided acceleration: %.3g\n', sqrt( mean( a.*a ) ));

compute_ms_error = @(spline) (mean(mean(  (data.x - spline(data.t)).^2,2 ),1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Toss out the points that are outliers, and perform the best fit.

spline_optimal = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
spline_optimal.minimizeMeanSquareError(data.t,data.x);

fprintf('optimal mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', compute_ms_error(spline_optimal), std(spline_optimal.uniqueValuesAtHighestDerivative),spline_optimal.lambda );
falseNegativeOutliers = setdiff(trueOutliers,spline_optimal.indicesOfOutliers);
falsePositiveOutliers = setdiff(spline_optimal.indicesOfOutliers,trueOutliers);
fprintf('False positives: %d, negatives: %d\n',length(falsePositiveOutliers),length(falseNegativeOutliers))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Try our old method

spline_robust_old = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
spline_robust_old.rebuildOutlierDistributionAndAdjustWeightingsOldMethod();
spline_robust_old.minimizeMeanSquareError(data.t,data.x);

fprintf('robust mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', compute_ms_error(spline_robust_old), std(spline_robust_old.uniqueValuesAtHighestDerivative),spline_robust_old.lambda );
falseNegativeOutliers = setdiff(trueOutliers,spline_robust_old.indicesOfOutliers);
falsePositiveOutliers = setdiff(spline_robust_old.indicesOfOutliers,trueOutliers);
fprintf('False positives: %d, negatives: %d\n',length(falsePositiveOutliers),length(falseNegativeOutliers))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Try our 1 pass optimal method

alpha = 0.5;
spline_robust = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
spline_robust.setToFullTensionWithInnerSVOnNoiseDistribution(alpha);
spline_robust.rebuildOutlierDistributionAndAdjustWeightings();
spline_robust.minimizeMeanSquareError(data.t,data.x);

fprintf('robust mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', compute_ms_error(spline_robust), std(spline_robust.uniqueValuesAtHighestDerivative),spline_robust.lambda );
falseNegativeOutliers = setdiff(trueOutliers,spline_robust.indicesOfOutliers);
falsePositiveOutliers = setdiff(spline_robust.indicesOfOutliers,trueOutliers);
fprintf('False positives: %d, negatives: %d\n',length(falsePositiveOutliers),length(falseNegativeOutliers))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Try our 2 pass optimal method

alpha = 0.5;
spline_robust2x = RobustTensionSpline(t_obs,x_obs,noiseDistribution,'S',S, 'lambda',Lambda.fullTensionExpected,'alpha',1/10000);
spline_robust2x.setToFullTensionWithInnerSVOnNoiseDistribution(alpha);
spline_robust2x.rebuildOutlierDistribution();
spline_robust2x.setToFullTensionWithInnerSV(alpha);
spline_robust2x.rebuildOutlierDistributionAndAdjustWeightings();
spline_robust2x.minimizeMeanSquareError(data.t,data.x);

fprintf('robust mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', compute_ms_error(spline_robust2x), std(spline_robust2x.uniqueValuesAtHighestDerivative),spline_robust2x.lambda );
falseNegativeOutliers = setdiff(trueOutliers,spline_robust2x.indicesOfOutliers);
falsePositiveOutliers = setdiff(spline_robust2x.indicesOfOutliers,trueOutliers);
fprintf('False positives: %d, negatives: %d\n',length(falsePositiveOutliers),length(falseNegativeOutliers))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure

figure
scatter(t_obs(spline_robust.indicesOfOutliers),x_obs(spline_robust.indicesOfOutliers),(8.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'), hold on
scatter(t_obs(trueOutliers),x_obs(trueOutliers),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
scatter(t_obs,x_obs,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline_optimal(tq),'k')
plot(tq,spline_robust(tq),'b')
plot(tq,spline_robust2x(tq),'r')
plot(tq,spline_robust_old(tq),'g')


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Ranged expected MSE

outlierDistribution = StudentTDistribution(200,3.0);
distribution = AddedDistribution(0.01,outlierDistribution,noiseDistribution);

spline = TensionSpline(t_obs,x_obs,distribution, 'S', S, 'lambda',Lambda.fullTensionExpected);
pctmin = 1/1000/2;
pctmax = 1-1/1000/2;
zmin = noiseDistribution.locationOfCDFPercentile(pctmin);
zmax = noiseDistribution.locationOfCDFPercentile(pctmax);

spline.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
compute_ms_error = @() (mean(mean(  (data.x(indices) - spline(data.t(indices))).^2,2 ),1));
fprintf('initial mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', compute_ms_error(), std(spline.uniqueValuesAtHighestDerivative), spline.lambda );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Toss outliers, do the same ranged MSE 

spline.indicesOfOutliers = find(abs(spline.epsilon) > outlierThreshold);
t_knot = InterpolatingSpline.KnotPointsForPoints(spline.t(spline.goodIndices),spline.K,1);
distribution = AddedDistribution(length(spline.indicesOfOutliers)/length(spline.t),outlierDistribution,noiseDistribution);

spline_reduced = TensionSpline(t_obs,x_obs,distribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'t_knot',t_knot);
compute_ms_error = @() (mean(mean(  (data.x(indices) - spline_reduced(data.t(indices))).^2,2 ),1));

spline_reduced.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
fprintf('reduced-range mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', compute_ms_error(), std(spline_reduced.uniqueValuesAtHighestDerivative), spline_reduced.lambda );

spline_reduced.minimize( @(spline) spline.expectedMeanSquareErrorFromGCV );
fprintf('reduced-gcv mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', compute_ms_error(), std(spline_reduced.uniqueValuesAtHighestDerivative), spline_reduced.lambda );

spline_reduced.minimize( @(spline) spline.expectedMeanSquareErrorFromCV );
fprintf('reduced-cv mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', compute_ms_error(), std(spline_reduced.uniqueValuesAtHighestDerivative), spline_reduced.lambda );

spline_reduced.indicesOfOutliers = find(abs(spline.epsilon) > outlierThreshold);

falseNegativeOutliers = setdiff(trueOutliers,spline_reduced.indicesOfOutliers);
falsePositiveOutliers = setdiff(spline_reduced.indicesOfOutliers,trueOutliers);
fprintf('False positives: %d, negatives: %d\n',length(falsePositiveOutliers),length(falseNegativeOutliers))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now re-adjust

distribution = AddedDistribution(length(spline_reduced.indicesOfOutliers)/length(spline_reduced.t),outlierDistribution,noiseDistribution);

spline_likelihood = TensionSpline(t_obs,x_obs,distribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'t_knot',t_knot);
tensionDistribution = NormalDistribution(TensionSpline.StandardDeviationOfInterquartileRange(spline.uniqueValuesAtHighestDerivative));
logLikelihood = @(spline) -sum(distribution.logPDF( spline.epsilon ) ) - sum(tensionDistribution.logPDF(spline.uniqueValuesAtHighestDerivative));
spline_likelihood.minimize( logLikelihood );
compute_ms_error = @() (mean(mean(  (data.x(indices) - spline_likelihood(data.t(indices))).^2,2 ),1));
fprintf('likelihood mse: %.2f m, acceleration: %.3g, lambda: %.3g\n', compute_ms_error(), std(spline_likelihood.uniqueValuesAtHighestDerivative), spline_likelihood.lambda );
spline_likelihood.indicesOfOutliers = find(abs(spline_likelihood.epsilon) > outlierThreshold);
falseNegativeOutliers = setdiff(trueOutliers,spline_likelihood.indicesOfOutliers);
falsePositiveOutliers = setdiff(spline_likelihood.indicesOfOutliers,trueOutliers);
fprintf('False positives: %d, negatives: %d\n',length(falsePositiveOutliers),length(falseNegativeOutliers))

figure
scatter(t_obs(spline_likelihood.indicesOfOutliers),x_obs(spline_likelihood.indicesOfOutliers),(8.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'), hold on
scatter(t_obs(trueOutliers),x_obs(trueOutliers),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(t_obs,x_obs,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline_optimal(tq),'k')
plot(tq,spline(tq),'b')
plot(tq,spline_reduced(tq),'r')
plot(tq,spline_likelihood(tq),'g')

[x_T, t_T] = spline_likelihood.uniqueValuesAtHighestDerivative;
hotspots = find( abs(x_T) > 3*tensionDistribution.sigma );
for i=1:length(hotspots)
    t_hot = linspace( t_obs(find(t_obs < t_T(hotspots(i)),1,'last')), t_obs(find(t_obs > t_T(hotspots(i)),1,'first')), 20 );
    plot(t_hot,spline_likelihood(t_hot),'r','LineWidth',2)
end

legend('flagged outliers', 'true outliers', 'data', 'optimal', 'initial', 'reduced', 'reduced optimal')

% lambda0 = spline.lambda;
% lambda = 10.^linspace(log10(lambda0/10),log10(10*lambda0),11)';
% expectedMSE = zeros(size(lambda));
% for iLambda = 1:length(lambda)
%    spline.lambda = lambda(iLambda);
%    expectedMSE(iLambda) = spline.expectedMeanSquareErrorInRange(zmin,zmax);
% end



return

spline.indicesOfOutliers = find( gaussian_cdf(spline.epsilon) < 0.05/2 | gaussian_cdf(spline.epsilon) > 1-0.05/2);
spline.goodIndices = setdiff(1:length(spline.x),spline.indicesOfOutliers);



t_knot = InterpolatingSpline.KnotPointsForPoints(t_obs(spline.goodIndices),spline.K,1);
spline_reduced = TensionSpline(t_obs,x_obs,sigma,'weightFunction', w_t, 'lambda',Lambda.fullTensionExpected,'t_knot',t_knot);
scatter(t_obs(spline_reduced.indicesOfOutliers),x_obs(spline_reduced.indicesOfOutliers),(3.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g'), hold on
plot(tq,spline_reduced(tq),'g')

return

w_t = @(z,sigma_t)((nu/(nu+1))*sigma_t^2*(1+z.^2/(nu*sigma_t^2)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nu = 4.5; sigma_t =  8.5;

a1 = gamma((nu+1)/2)./(sqrt(pi*nu)*sigma_t*gamma(nu/2));
c1 = nu*sigma_t*sigma_t;
m1 = (nu+1)/2;
variance_of_the_noise_1 = sigma_t*sigma_t*nu/(nu-2)

nu = 4.5; sigma_t =  200;
a2 = gamma((nu+1)/2)./(sqrt(pi*nu)*sigma_t*gamma(nu/2));
c2 = nu*sigma_t*sigma_t;
m2 = (nu+1)/2;
variance_of_the_noise_2 = sigma_t*sigma_t*nu/(nu-2)


alpha = 0.1;
a1 = (1-alpha)*a1;
a2 = alpha*a2;

variance_of_the_noise = (1-alpha)*variance_of_the_noise_1 + alpha*variance_of_the_noise_2;

w_tt = @(z) (a1*(1+z.*z/c1).^(-m1) + a2*(1+z.*z/c2).^(-m2))./(2*(a1*m1/c1)*(1+z.*z/c1).^(-m1-1) + 2*(a2*m2/c2)*(1+z.*z/c2).^(-m2-1));

nu = 4.5; sigma_t =  8.5;

a1 = gamma((nu+1)/2)./(sqrt(pi*nu)*sigma_t*gamma(nu/2));
c1 = nu*sigma_t*sigma_t;
m1 = (nu+1)/2;

sigma_t =  85;
a2 = 1/(sigma_t*sqrt(2*pi));
c2 = sigma_t*sigma_t;

alpha = 0.5;
a1 = (1-alpha)*a1;
a2 = alpha*a2;


w_tg = @(z) (a1*(1+z.*z/c1).^(-m1) + a2*exp(-z.*z/(2*c2)))./(2*(a1*m1/c1)*(1+z.*z/c1).^(-m1-1) + (a2/c2)*exp(-z.*z/(2*c2)));

sigma_t =  8.5;
a1 = 1/(sigma_t*sqrt(2*pi));
c1 = sigma_t*sigma_t;

sigma_t =  85;
a2 = 1/(sigma_t*sqrt(2*pi));
c2 = sigma_t*sigma_t;

alpha = 0.5;
a1 = (1-alpha)*a1;
a2 = alpha*a2;

w_gg = @(z) (a1*exp(-z.*z/(2*c1)) + a2*exp(-z.*z/(2*c2)))./((a1/c1)*exp(-z.*z/(2*c1)) + (a2/c2)*exp(-z.*z/(2*c2)));







