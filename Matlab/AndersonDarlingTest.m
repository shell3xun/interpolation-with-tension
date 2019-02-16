n = 250;
percentOutliers = 0.25;
outlierFactor = 40;
nu = 4.5; sigma =  8.5;
noiseDistribution = StudentTDistribution(sigma,nu);
outlierDistribution = StudentTDistribution(outlierFactor*sigma,3.0);

outlierIndices = rand(n,1)<=percentOutliers;
outlierIndices([1 n]) = 0;

epsilon_x = zeros(n,1);
epsilon_x_outlier = zeros(n,1);

epsilon_x(~outlierIndices) = noiseDistribution.rand(sum(~outlierIndices));
epsilon_x_outlier(outlierIndices) = outlierDistribution.rand(sum(outlierIndices));

epsilon = epsilon_x + epsilon_x_outlier;

addedDistribution = AddedDistribution(percentOutliers,outlierDistribution,noiseDistribution);
fprintf('sqrt(variance) (noise, outlier)=(%.1f, %.1f)\n',sqrt(noiseDistribution.variance),sqrt(outlierDistribution.variance));
fprintf('sqrt(sample variance) (noise, outlier)=(%.1f, %.1f)\n',sqrt(mean(epsilon(~outlierIndices).^2)),sqrt(mean(epsilon(outlierIndices).^2)));
fprintf('true alpha %.2f\n',sum(outlierIndices)/n);

s2_total = mean(epsilon.^2);
s2_noise = noiseDistribution.variance;

if s2_total/s2_noise < 1.1
    outlierDistribution = [];
    alpha = 0;
    return;
end

alpha_outlier = 10.^(linspace(log10(0.01),log10(0.5),100))';
sigma2_outlier = (s2_total-(1-alpha_outlier)*s2_noise)./(3*alpha_outlier);
var_total = zeros(size(alpha_outlier));
ks_error = zeros(size(alpha_outlier));

for iAlpha = 1:length(alpha_outlier)
    newAddedDistribution = AddedDistribution(alpha_outlier(iAlpha),StudentTDistribution(sqrt(sigma2_outlier(iAlpha)),3),noiseDistribution);
    var_total(iAlpha) = newAddedDistribution.variance;
    ks_error(iAlpha) = newAddedDistribution.andersonDarlingError(epsilon);
end

[minAlpha,minIndex] = min(ks_error);

fprintf('Best distribution match, alpha=%.2f and sqrt(var)=%.1f m\n',alpha_outlier(minIndex),sqrt(3*sigma2_outlier(minIndex)));