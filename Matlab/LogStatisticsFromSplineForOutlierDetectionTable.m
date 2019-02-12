function stats = LogStatisticsFromSplineForOutlierDetectionTable(stats, linearIndex, spline,compute_ms_error,trueOutlierIndices,outlierIndices)

stats.mse(linearIndex) = compute_ms_error(spline);
stats.neff_se(linearIndex) = spline.effectiveSampleSizeFromVarianceOfTheMean;
stats.lambda(linearIndex) = spline.lambda;
stats.nonOutlierEffectiveSampleSize(linearIndex) = spline.effectiveSampleSizeFromVarianceOfTheMeanForIndices(~outlierIndices);
stats.nonOutlierSampleVariance(linearIndex) = mean(spline.epsilonAtIndices(~outlierIndices).^2);
stats.false_negatives(linearIndex) = length(setdiff(trueOutlierIndices,spline.indicesOfOutliers));
stats.false_positives(linearIndex) = length(setdiff(spline.indicesOfOutliers,trueOutlierIndices));

end