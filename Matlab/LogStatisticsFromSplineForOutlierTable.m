function stats = LogStatisticsFromSplineForOutlierTable(stats, linearIndex, spline,compute_ms_error,trueOutlierIndices,outlierIndices,varargin)

stats.mse(linearIndex) = compute_ms_error(spline);
stats.neff_se(linearIndex) = spline.effectiveSampleSizeFromVarianceOfTheMean;
stats.lambda(linearIndex) = spline.lambda;
stats.nonOutlierEffectiveSampleSize(linearIndex) = spline.effectiveSampleSizeFromVarianceOfTheMeanForIndices(~outlierIndices);
stats.nonOutlierSampleVariance(linearIndex) = mean(spline.epsilonAtIndices(~outlierIndices).^2);
stats.false_negatives(linearIndex) = length(setdiff(trueOutlierIndices,spline.outlierIndices));
stats.false_positives(linearIndex) = length(setdiff(spline.outlierIndices,trueOutlierIndices));
if ~isempty(spline.alpha)
    stats.alpha(linearIndex) = spline.alpha;
end
if length(varargin) >= 1
    stats.rejects(linearIndex) = varargin{1};
end


end