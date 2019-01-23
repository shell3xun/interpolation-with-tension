function stats = LogStatisticsFromSplineForOutlierTable(stats, linearIndex, spline,compute_ms_error,trueOutlierIndices)

stats.mse(linearIndex) = compute_ms_error(spline);
stats.neff_se(linearIndex) = spline.effectiveSampleSizeFromVarianceOfTheMean;
stats.false_negatives(linearIndex) = length(setdiff(trueOutlierIndices,spline.indicesOfOutliers));
stats.false_positives(linearIndex) = length(setdiff(spline.indicesOfOutliers,trueOutlierIndices));

end