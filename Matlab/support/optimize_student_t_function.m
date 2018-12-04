nu = 4.5; sigma =  8.5;
% [r, pdf1d] = TwoDimStudentTProbabilityDistributionFunction( sigma, nu, 150, 1501 );

% scale our axis to be 20 times the standard deviation of the
% 1 dimensional pdf.
maxR = 20*sqrt(sigma*sigma*nu/(nu-2));

% This appears to be big enough, although we can go bigger if
% needed.
N = 150;

% create a logarithmic axis that includes zero
r = [0;10.^(linspace(log10(maxR/1e3),log10(maxR),N))'];
x = [-flip(r(2:end),1); r];

% evaluate the cdf on that axis
pdf_int = tcdf(x/sigma,nu);

% use that to define the pdf...
dx = diff(x);
pdf = diff(pdf_int)./dx;
%...so that sum(pdf.*dx)=1.00000

% now create a 2d version
pdf2d_norm = (pdf .* reshape(pdf,1,[])) .* (dx .* reshape(dx,1,[]));
% again, where sum(pdf2d_norm(:)) = 1.0000

% create a 2d distance metric, rho
y = [-flip(r(2:end)); r(2:end)];
rho = sqrt( y.*y + reshape(y.*y,1,[]));

% we are going to bin using the diagonal of rho, call it s.
s = diag(rho);
s = [0; s((N+1):end)];

% Now sum up the energy in each bin
pdf1d = zeros( size(s) );
midS = [0; s(1:(end-1))+diff(s)/2];
for i = 1:(length(s)-1)
    pdf1d(i) = sum( pdf2d_norm(  rho >= midS(i) & rho <= midS(i+1) ) );
end
% it must be true that sum(pdf1d) = 1.0000

% but we want that sum(pdf1d .* diff(s)) = 1.000
pdf1d(2:end) = pdf1d(2:end)./diff(s);
cdf1d = [0; cumsum(pdf1d(2:end).*diff(s))];

% now take care of the the fact that we used the diagnoal of a
% 2d rectangle, rather than the inscribed circle.

pdf = pdf1d(s<maxR);
cdf = cdf1d(s<maxR);
r = s(s<maxR);



return

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