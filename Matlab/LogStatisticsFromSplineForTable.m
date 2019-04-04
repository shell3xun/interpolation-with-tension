function stats = LogStatisticsFromSplineForTable(stats, linearIndex, spline,compute_ms_error,varargin)

stats.mse(linearIndex) = compute_ms_error(spline);
stats.neff_se(linearIndex) = spline.effectiveSampleSizeFromVarianceOfTheMean;
stats.lambda(linearIndex) = spline.lambda;

end