N=10000;
t=(0:N-1).';

D = 1;

if 1 == 1
    nu = 4.5; sigma =  8.5;
    variance_of_the_noise = sigma*sigma*nu/(nu-2);
    epsilon = randt(sigma,nu,N);
else
    sigma = 10;
    variance_of_the_noise = sigma*sigma;
    epsilon = sigma*randn(N,1);
end

x = epsilon;

% Now do some signal processing
dt = t(2) - t(1);
T = t(end)-t(1);
df = 1/T;
f = ([0:ceil(N/2)-1 -floor(N/2):-1]*df)';

ubar = fft(x);
s_signal = (ubar.*conj(ubar)) .* (2*pi*f).^(2*D) * (dt/N);
s_noise = variance_of_the_noise*dt*(2*pi*f).^(2*D);

alpha = 0.999;
dof = 2;
cutoff = TensionSpline.chi2inv(alpha,dof)/dof;

f = fftshift(f);
s_signal = fftshift(s_signal);
s_noise = fftshift(s_noise);

figure
plot(f,s_signal)
hold on
plot(f,cutoff*s_noise), ylog

actual_alpha = sum( (s_signal < cutoff*s_noise) )/N;

fprintf('Expected %.4f, actual %.4f\n',alpha,actual_alpha);