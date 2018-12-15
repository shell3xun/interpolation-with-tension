function likelihood = LogLikelihoodOfTandG(spline_x, sigma_t, nu, sigma_g)

% these are the logs of the pdfs
log_t_pdf = @(z) (-(nu+1)/2)*log(1+(z.*z)/(nu*sigma_t*sigma_t));
log_gaussian_pdf = @(z) -(z.*z)/(2*sigma_g*sigma_g);

error_likelihood = -sum(log_t_pdf(spline_x.epsilon));
acceleration_lielihood = -sum(log_gaussian_pdf(spline_x.UniqueValuesAtHighestDerivative()));

likelihood = error_likelihood + acceleration_lielihood;

end

