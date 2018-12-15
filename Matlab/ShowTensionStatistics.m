% Site 1 or site 2?
Site = 1;

% Drifter to highlight in the final plots. Drifter 7 has no outliers
choiceDrifter = 7;

if Site == 1
    drifters = load('sample_data/rho1_drifters_projected_ungridded.mat');
else
    drifters = load('sample_data/rho2_drifters_projected_ungridded.mat');
end

% Pull out the data of interest
x_data = drifters.x{choiceDrifter};
y_data = drifters.y{choiceDrifter};
t_data = drifters.t{choiceDrifter};

K = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Low order finite differencing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = TensionSpline.FiniteDifferenceMatrixNoBoundary(K-1,drifters.t{choiceDrifter},1);
DxDT = cat(1,D*drifters.x{choiceDrifter},D*drifters.y{choiceDrifter});

sigma = 8.5;
nu = 4.5;
variance_of_the_noise = sigma*sigma*nu/(nu-2);
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));

spline_x = TensionSpline(t_data, x_data, sqrt(variance_of_the_noise), 'K', K, 'weightFunction', w); %, 'lambda',0);
spline_y = TensionSpline(t_data, y_data, sqrt(variance_of_the_noise), 'K', K, 'weightFunction', w); %, 'lambda',0);

spline_x.Minimize( @(spline) LogLikelihoodOfTandG(spline,sigma,nu,sigma_FD))

x_T = spline_x.UniqueValuesAtHighestDerivative();
y_T = spline_y.UniqueValuesAtHighestDerivative();

sigma_Tx = TensionSpline.StandardDeviationOfInterquartileRange(x_T);
sigma_Ty = TensionSpline.StandardDeviationOfInterquartileRange(y_T);
sigma_T = TensionSpline.StandardDeviationOfInterquartileRange(cat(1,x_T,y_T));


figure
histogram(DxDT,100, 'Normalization', 'pdf'), hold on
histogram(cat(1,x_T,y_T),100, 'Normalization', 'pdf')

sigma_FD = std(DxDT)
sigma_spline = std(cat(1,x_T,y_T))

figure,
histogram(D*drifters.x{choiceDrifter},100, 'Normalization', 'pdf'), hold on
histogram(spline_x.UniqueValuesAtHighestDerivative(),100, 'Normalization', 'pdf')
% histogram(D*drifters.y{choiceDrifter},100, 'Normalization', 'pdf')

sigma_Dx = std(D*drifters.x{choiceDrifter});
sigma_Dy = std(D*drifters.y{choiceDrifter});