% For fixed total variance, a larger alpha (percentOutliers) results in
% smaller variance near zero. Hence, over estimating alpha causes
percentOutliers = [0.10 0.2 0.3];
outlierFactor = 40;
nu = 4.5; sigma =  8.5;
noiseDistribution = StudentTDistribution(sigma,nu);
totalVariance = 1e5;

z = linspace(0,100,100)';
figure
plot(z,cumtrapz(z,z.*z.*noiseDistribution.pdf(z)))
for iPct = 1:length(percentOutliers)
    outlierDistribution = StudentTDistribution(sqrt((totalVariance-(1-percentOutliers(iPct))*noiseDistribution.variance)/(3*percentOutliers(iPct))),3.0);
    
    addedDistribution = AddedDistribution(percentOutliers(iPct),outlierDistribution,noiseDistribution);
    totalVariance = addedDistribution.variance
    
    hold on
    plot(z,cumtrapz(z,z.*z.*addedDistribution.pdf(z)))
end


