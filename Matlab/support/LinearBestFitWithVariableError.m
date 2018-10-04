function [m,b] = LinearBestFitWithVariableError(x,y,sigma2)
%LinearBestFitWithVariableError Fit y=m*x+b to the data with measurement
%error of sigma2 in y at each x.
S = sum(1./sigma2);
Sx = sum(x./sigma2);
Sy = sum(y./sigma2);
Sxx = sum((x.*x)./sigma2);
Sxy = sum((x.*y)./sigma2);
Delta = S*Sxx - Sx*Sx;
b = (Sxx*Sy - Sx*Sxy)/Delta;
m = (S*Sxy - Sx*Sy)/Delta;
end

