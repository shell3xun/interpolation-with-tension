N = 100;
sigma = 8.5;
nu = 4.5;

numRepeats = 1000;
numRejections = zeros(numRepeats,1);

for iRepeat = 1:numRepeats
    
    epsilon = randt(sigma, nu, N);
    z = epsilon/sigma;
    cdf = 2*tcdf(abs(z),nu)-1;
    jpd = cdf.^N;
    p = 0.05;
    
    numRejections(iRepeat) = sum(jpd > 1-p);
end

fprintf('You set the p-value to 0.05. In practice, you discarded a point in %.1f%% of datasets',100*sum(numRejections)/numRepeats);