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
result_stride = 64;
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
                
                epsilon_x = zeros(length(indices),1);
                epsilon_x_outlier = zeros(length(indices),1);
                
                epsilon_x(~outlierIndices) = sigma*randn(sum(~outlierIndices),1);
                epsilon_x_outlier(outlierIndices) = outlierFactor*sigma*randn(sum(outlierIndices),1);
                                
                variance_of_the_noise = (1-percentOutliers)*sigma*sigma + percentOutliers*(outlierFactor*sigma)^2;
                w = [];
            end
            
            x_obs = data.x(indices) + epsilon_x + epsilon_x_outlier;
            t_obs = data.t(indices);
            
            
        end
        fprintf('\n');
    end
    
end

% c = 2.1;
% w_tukey = @(epsilon) (sigma*sigma*(1-(epsilon/(c*sigma)).^2).^(-2)) .* (abs(epsilon) < c*sigma) + 1e8*sigma*sigma * (abs(epsilon) >= c*sigma);
% 
% nu = 4.5; sigma_t =  8.5;
% w_t = @(z)((nu/(nu+1))*sigma_t^2*(1+z.^2/(nu*sigma_t^2)));
% t_pdf = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma_t*gamma(nu/2)*(1+(z.*z)/(nu*sigma_t*sigma_t)).^((nu+1)/2));
% mu =0;
% sigma_g = 10;
% gaussian_cdf = @(z) 0.5*(1 + erf((z-mu)/(sigma_g*sqrt(2))));


sigma_t =  sigma;
a1 = 1/(sigma_t*sqrt(2*pi));
c1 = sigma_t*sigma_t;

sigma_t =  outlierFactor*sigma;
a2 = 1/(sigma_t*sqrt(2*pi));
c2 = sigma_t*sigma_t;

alpha = percentOutliers;
a1 = (1-alpha)*a1;
a2 = alpha*a2;

w_gg = @(z) (a1*exp(-z.*z/(2*c1)) + a2*exp(-z.*z/(2*c2)))./((a1/c1)*exp(-z.*z/(2*c1)) + (a2/c2)*exp(-z.*z/(2*c2)));

spline = TensionSpline(t_obs,x_obs,sqrt(variance_of_the_noise),'weightFunction', w_gg, 'lambda',Lambda.optimalIterated);
tq = linspace(min(t_obs),max(t_obs),10*length(t_obs));

figure
scatter(t_obs(spline.indicesOfOutliers),x_obs(spline.indicesOfOutliers),(8.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'), hold on
scatter(t_obs(outlierIndices),x_obs(outlierIndices),(6.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'), hold on
scatter(t_obs,x_obs,(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
plot(tq,spline(tq),'r')

lambda0 = spline.lambda;
lambda = 10.^linspace(log10(lambda0/10),log10(10*lambda0),11)';
expectedMSE = zeros(size(lambda));
for iLambda = 1:length(lambda)
   spline.lambda = lambda(iLambda);
   expectedMSE(iLambda) = spline.ExpectedMeanSquareError;
end

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

nu = 4.5; sigma_t =  85;
a2 = gamma((nu+1)/2)./(sqrt(pi*nu)*sigma_t*gamma(nu/2));
c2 = nu*sigma_t*sigma_t;
m2 = (nu+1)/2;

alpha = 0.5;
a1 = (1-alpha)*a1;
a2 = alpha*a2;

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







