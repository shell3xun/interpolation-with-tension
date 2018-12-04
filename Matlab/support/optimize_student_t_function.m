nu = 4.5; sigma =  8.5;
% [r, pdf1d] = TwoDimStudentTProbabilityDistributionFunction( sigma, nu, 150, 1501 );

% trapz(r,pdf1d)

maxR = 150;
N = 100;

r = [0;10.^(linspace(log10(maxR/1e2),log10(maxR),N))'];
x = [-flip(r(2:end),1); r];
dx = [flip(diff(r)); r(2)-r(1); diff(r)];

trapz(x,pdf(x))

[X,Y] = ndgrid(x,x);
rho = sqrt( X.*X + Y.*Y);

[dX,dY] = ndgrid(dx,dx);

pdf = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));
pdf2d = pdf(X).*pdf(Y);
pdf2d_norm = pdf2d .* dX .* dY;

% Now sum up the energy
pdf1d = zeros( size(r) );
for i = 2:length(r)
	pdf1d(i) = sum( pdf2d_norm(  rho >= r(i-1) & rho < r(i) )/(r(i)-r(i-1)) );
end

tcdf(1.0/sigma,nu)