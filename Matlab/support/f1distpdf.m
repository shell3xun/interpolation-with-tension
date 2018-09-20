function [y] =f1distpdf(x,nu)
% analytical PDF of a F(1,n) distribution
CC=gamma((1+nu)/2)/(sqrt(pi)*gamma(nu/2))*sqrt(1/nu);
y = CC./(sqrt(x).*(1+x/nu).^((1+nu)/2));