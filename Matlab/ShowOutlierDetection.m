scaleFactor = 1;
LoadFigureDefaults

shouldUseStudentTDistribution = 0;

percentOutliers = 0.1;
outlierFactor = 20;

slopes = [-2; -3; -4];
totalSlopes = length(slopes);

S = 2;
T = S;
K = S+1;

result_stride = 2*[1;4;16;64];
result_stride = 128;
totalStrides = length(result_stride);

totalEnsembles = 1;

for iSlope = 3:3;%length(slopes)
    slope = slopes(iSlope);
    fprintf('slope %d, ',slope);
    
    if slope == -2
        data = load('sample_data/SyntheticTrajectories.mat');
        outputFile = 'OptimalParameters.mat';
    elseif slope == -3
        data = load('sample_data/SyntheticTrajectoriesSlope3.mat');
        outputFile = 'OptimalParametersSlope3.mat';
    elseif slope == -4
        data = load('sample_data/SyntheticTrajectoriesSlope4.mat');
        outputFile = 'OptimalParametersSlope4.mat';
    end
    
    sigma = data.position_error;
    dt = data.t(2)-data.t(1);
    
    for iStride=1:length(result_stride)
        stride = result_stride(iStride);
        % Reduce the total length in some cases
        if (stride < 10)
            shortenFactor = stride/10;
        else
            shortenFactor = 1;
        end
        
        indices = 1:stride:floor(shortenFactor*length(data.t));
        fprintf('dT = %d --- Using %d points with stride %d. Evaluating ensemble', stride*dt, length(indices), stride);

        for iEnsemble = 1:totalEnsembles
            fprintf('..%d',iEnsemble);
            
            if shouldUseStudentTDistribution == 1
                nu = 4.5; sigma =  8.5;
                variance_of_the_noise = sigma*sigma*nu/(nu-2);
                w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
                epsilon_x = randt(sigma,nu,length(indices));
            else
                outlierIndices = rand(length(indices),1)<=percentOutliers;
                outlierIndices([1 length(indices)]) = 0;
                
                epsilon_x = zeros(length(indices),1);
                epsilon_x_outlier = zeros(length(indices),1);
                
                epsilon_x(~outlierIndices) = sigma*randn(sum(~outlierIndices),1);
                epsilon_x_outlier(outlierIndices) = outlierFactor*sigma*randn(sum(outlierIndices),1);
                                
                variance_of_the_noise = (1-percentOutliers)*sigma*sigma + percentOutliers*(outlierFactor*sigma)^2;
                w = [];
                
                epsilon = epsilon_x + epsilon_x_outlier;
                trueOutliers = find(abs(epsilon) > 4*sigma);
                goodIndices = setdiff(1:length(indices),trueOutliers);
            end
            
            x_obs = data.x(indices) + epsilon_x + epsilon_x_outlier;
            t_obs = data.t(indices);
            
            
        end
        fprintf('\n');
    end
    
end
tq = linspace(min(t_obs),max(t_obs),10*length(t_obs));
% [D,tD] = TensionSpline.FiniteDifferenceMatrixNoBoundary(2,data.t(indices),2);
% spline_clean = TensionSpline(data.t(indices),data.x(indices),noiseDistribution, 'S', S, 'lambda',0);
% 
% figure
% plot(tq,spline_clean(tq,2)), hold on
% plot(tD,D*data.x(indices))

dt = t_obs(2)-t_obs(1);
a = diff(diff(data.x(indices)))/(dt^2);
fprintf('true strided acceleration: %.3g\n', sqrt( mean( a.*a ) ));

% Toss out the points that are outliers, and perform the best fit.
noiseDistribution = NormalDistribution(sigma);
spline_optimal = TensionSpline(t_obs(goodIndices),x_obs(goodIndices),noiseDistribution, 'S', S, 'lambda',Lambda.fullTensionExpected);
spline_optimal.minimizeMeanSquareError(data.t(indices),data.x(indices));

compute_ms_error = @() (mean(mean(  (data.x(indices) - spline_optimal(data.t(indices))).^2,2 ),1));
fprintf('optimal mse: %.2f m, acceleration: %.3g\n', compute_ms_error(), std(spline_optimal.uniqueValuesAtHighestDerivative) );

% distribution = AddedDistribution(percentOutliers,NormalDistribution(outlierFactor*sigma),StudentTDistribution(8.5,4.5));

outlierDistribution = StudentTDistribution(200,3.0);
distribution = AddedDistribution(0.01,outlierDistribution,noiseDistribution);

spline = TensionSpline(t_obs,x_obs,distribution, 'S', S, 'lambda',Lambda.fullTensionExpected);
pctmin = 1/100/2;
pctmax = 1-1/100/2;
zmin = noiseDistribution.locationOfCDFPercentile(pctmin);
zmax = noiseDistribution.locationOfCDFPercentile(pctmax);

spline.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
compute_ms_error = @() (mean(mean(  (data.x(indices) - spline(data.t(indices))).^2,2 ),1));
fprintf('initial mse: %.2f m, acceleration: %.3g\n', compute_ms_error(), std(spline.uniqueValuesAtHighestDerivative) );


spline.indicesOfOutliers = find(abs(spline.epsilon) > 4*sigma);
spline.indicesOfOutliers = setdiff(spline.indicesOfOutliers,[1 length(spline.t)]);
spline.goodIndices = setdiff(1:length(indices),spline.indicesOfOutliers);
t_knot = InterpolatingSpline.KnotPointsForPoints(spline.t(spline.goodIndices),spline.K,1);

spline_reduced = TensionSpline(t_obs,x_obs,distribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'t_knot',t_knot);
spline_reduced.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax) );
compute_ms_error = @() (mean(mean(  (data.x(indices) - spline_reduced(data.t(indices))).^2,2 ),1));
fprintf('reduced mse: %.2f m, acceleration: %.3g\n', compute_ms_error(), std(spline_reduced.uniqueValuesAtHighestDerivative) );

falseNegativeOutliers = setdiff(trueOutliers,spline.indicesOfOutliers);
falsePositiveOutliers = setdiff(spline.indicesOfOutliers,trueOutliers);
fprintf('False positives: %d, negatives: %d\n',length(falsePositiveOutliers),length(falseNegativeOutliers))

spline_likelihood = TensionSpline(t_obs,x_obs,distribution, 'S', S, 'lambda',Lambda.fullTensionExpected,'t_knot',t_knot);
tensionDistribution = NormalDistribution(TensionSpline.StandardDeviationOfInterquartileRange(spline.uniqueValuesAtHighestDerivative));
logLikelihood = @(spline) -sum(distribution.logPDF( spline.epsilon ) ) - sum(tensionDistribution.logPDF(spline.uniqueValuesAtHighestDerivative));
spline_likelihood.minimize( logLikelihood );
compute_ms_error = @() (mean(mean(  (data.x(indices) - spline_likelihood(data.t(indices))).^2,2 ),1));
fprintf('initial mse: %.2f m, acceleration: %.3g\n', compute_ms_error(), std(spline_likelihood.uniqueValuesAtHighestDerivative) );
spline_likelihood.indicesOfOutliers = find(abs(spline_likelihood.epsilon) > 4*sigma);
falseNegativeOutliers = setdiff(trueOutliers,spline_likelihood.indicesOfOutliers);
falsePositiveOutliers = setdiff(spline_likelihood.indicesOfOutliers,trueOutliers);
fprintf('False positives: %d, negatives: %d\n',length(falsePositiveOutliers),length(falseNegativeOutliers))

figure
scatter(t_obs(spline.indicesOfOutliers),x_obs(spline.indicesOfOutliers),(8.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'), hold on
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







