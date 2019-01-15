
position_error = 10; % let's assume 10 meter positioning errors
sigma_u = 0.20; % 20 cm/s rms speed
gamma = 10; % this is my interpolation condition parameter, over sample
dt = position_error/(sigma_u*gamma); % 8 minute 20 second intervals
n = round(1*86400/dt) % six days
t_damp = 30*60; % how long it takes to revert to the mean (of zero)
alpha = 2.0; % -2 slope

tic
cv=maternoise(dt,n,sigma_u*sqrt(2),alpha,1/t_damp);
cx = cumtrapz(cv)*dt;
toc

t = dt*(0:n-1)';
x = real(cx);
y = imag(cx);
u = real(cv);
v = imag(cv);

epsilon_x = position_error*randn(size(t));
epsilon_y = position_error*randn(size(t));

notes = 'The trajectories (t,x,y) are generated from a Matern using jlab with maternoise(dt,n,sigma_u*sqrt(2),alpha,1/t_damp). The position errors (epsilon_x, epsilon_y) are generated with position_error*randn(size(t)).';

% save('sample_data/SyntheticTrajectoriesSlope4.mat', 't', 'x', 'y', 'sigma_u', 't_damp', 'alpha', 'epsilon_x', 'epsilon_y', 'position_error', 'notes')

figure
subplot(3,1,1)
plot(x,y)
subplot(3,1,2)
plot(t,[u,v])


[psi,lambda]=sleptap(n);
[f,spp,snn,spn]=mspec(cv,psi);

subplot(3,1,3)
plot(f,[snn, spp])
xlog, ylog
xlim([min(f) max(f)])
ylim([0.9*min(spp) 1.1*max(spp)])