percentOutliers = 0.25;
outlierFactor = 40;
nu = 4.5; sigma =  8.5;
noiseDistribution = StudentTDistribution(sigma,nu);
outlierDistribution = StudentTDistribution(outlierFactor*sigma,3.0);

addedDistribution = AddedDistribution(percentOutliers,outlierDistribution,noiseDistribution);

alphaValues = 10.^(linspace(log10(1e-3),log10(.5),100)');
mse = zeros(size(alphaValues));
m = 2/3;
f = @(alpha) ((addedDistribution.variance-2*addedDistribution.varianceInRange(-Inf,addedDistribution.locationOfCDFPercentile(alpha/2)))^(1-m/2))/(1-alpha)^(m-1) + ((2*addedDistribution.varianceInRange(-Inf,addedDistribution.locationOfCDFPercentile(alpha/2)))^(1-m/2))/(alpha)^(m-1);

for iAlpha = 1:length(alphaValues)
    alpha = alphaValues(iAlpha);
%     var1 = 2*addedDistribution.varianceInRange(-Inf,addedDistribution.locationOfCDFPercentile(alpha/2));
%     var2 = addedDistribution.variance-var1;
%     mse(iAlpha) = (var2^(1-m/2))/(1-alpha)^(m-1) + (var1^(1-m/2))/(alpha)^(m-1);
    mse(iAlpha) = f(alpha);
end


figure
plot(alphaValues,mse),xlog
% hold on, plot(alphaValues,f(alphaValues))
% 
% [minMSE,index] = min(mse);
% minAlpha = alphaValues(index);
minAlpha = fminsearch(f,0.01);
zAlpha = addedDistribution.locationOfCDFPercentile(minAlpha/2);
cdfRatio = (percentOutliers/(1-percentOutliers))*outlierDistribution.cdf(zAlpha)/noiseDistribution.cdf(zAlpha);
varOutlier = 2*addedDistribution.varianceInRange(-Inf,zAlpha);
varNoise = addedDistribution.variance-varOutlier;
fprintf('Cutoff in z: %.1f m has a cdf ratio of: %.1g. There is %.0f m^2 of outlier variance and %.1f m^2 of noise variance.\n',zAlpha,cdfRatio,varOutlier,varNoise);