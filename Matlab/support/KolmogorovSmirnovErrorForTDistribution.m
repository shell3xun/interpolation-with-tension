function totalError = KolmogorovSmirnovErrorForTDistribution( epsilon, sigma, nu, range)
%%KolmogorovSmirnovErrorForTDistribution
%
% Return a measure of the Kolmogorov-Smirnov test metric for the
% distribution matched to the position errors
%
% the range allows you to restrict the range (in meters) over which the
% test is applied.

gps_cdf = @(z) tcdf(z/sigma,nu);

if nargin == 4
    x = sort(epsilon( epsilon > range(1) & epsilon < range(2) ));
else
    x = sort(epsilon);
end

n = length(x);
y_data = (1:n)'/n;
y = gps_cdf(x);

D = max(abs(y-y_data));

totalError = (sqrt(n) + 0.12 + 0.11/sqrt(n))*D;


end

