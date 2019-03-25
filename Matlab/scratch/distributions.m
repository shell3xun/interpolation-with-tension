nu = 4.5; sigma =  8.5;
noiseDistribution = StudentTDistribution(sigma,nu);

nu = 3.0;
sigma = sqrt(noiseDistribution.variance*1000*(nu-2)/nu);
outlierDistribution = StudentTDistribution(sigma,nu);


alphaValues = [0, 1/100, 1/10000];
betaValues = 1./[50, 100, 200, 400, 800];

for beta=betaValues
    for alpha=alphaValues
        
        distribution = AddedDistribution(alpha,outlierDistribution,noiseDistribution);
        
        zmin = noiseDistribution.locationOfCDFPercentile(beta/2);
        zmax = noiseDistribution.locationOfCDFPercentile(1-beta/2);
        var = distribution.varianceInRange(zmin,zmax);
        fprintf('alpha = %.2g, beta = %.2g, variance=%.2f\n',alpha,beta,var);
    end
end