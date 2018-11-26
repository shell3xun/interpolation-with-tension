% Assuming you have two dimensions with independent t-distributions, this
% gives you the CDF in polar coordinate r, after integrating theta.
function [r, pdf1d] = TwoDimStudentTProbabilityDistributionFunction( sigma, nu, maxR, N )

if nargin < 4
    N = 1001;
end

pdf = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));

x = linspace(-maxR,maxR,N)';
y = linspace(-maxR,maxR,N)';
dR = x(2)-x(1);
r = (0:dR:maxR)';

% Radial vector, spaced unevenly
% r = (linspace(0,sqrt(maxR),1001).^2)';
% x = [-flip(r(2:end),1); r];
% y = x;

[X,Y] = meshgrid(x,y);
rho = sqrt( X.*X + Y.*Y);


pdf2d = pdf(X).*pdf(Y);
pdf2d_norm = pdf2d * dR;

% Now sum up the energy
pdf1d = zeros( size(r) );
for i = 1:length(r)
	pdf1d(i) = sum( pdf2d_norm(  rho >= r(i)-dR/2 & rho < r(i)+dR/2 ) );
end
