% For fixed total variance, a larger alpha (percentOutliers) results in
% smaller variance near zero. Hence, over estimating alpha causes
percentOutliers = 0.15;
outlierFactor = 40;
nu = 4.5; sigma =  8.5;
noiseDistribution = StudentTDistribution(sigma,nu);
totalVariance = 1e5;
outlierDistribution = StudentTDistribution(sqrt((totalVariance-(1-percentOutliers)*noiseDistribution.variance)/(3*percentOutliers)),3.0);

addedDistribution = AddedDistribution(percentOutliers,outlierDistribution,noiseDistribution);
totalVariance = addedDistribution.variance

z = linspace(0,300,100)';
hold on
plot(z,cumtrapz(z,z.*z.*addedDistribution.pdf(z)))