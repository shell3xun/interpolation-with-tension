function stats = LogStatisticsFromSplineForFullTensionTable(stats, linearIndex, spline,compute_ms_error,zCrossover,outlierCrossoverIndices,outlierIndices)

stats.mse(linearIndex) = compute_ms_error(spline);
stats.neff_se(linearIndex) = spline.effectiveSampleSizeFromVarianceOfTheMean;
stats.lambda(linearIndex) = spline.lambda;
stats.nonOutlierEffectiveSampleSize(linearIndex) = spline.effectiveSampleSizeFromVarianceOfTheMeanForIndices(~outlierIndices);
stats.nonOutlierSampleVariance(linearIndex) = mean(spline.epsilonAtIndices(~outlierIndices).^2);

indices = find(abs(spline.epsilon)>zCrossover);
stats.false_negatives(linearIndex) = length(setdiff(outlierCrossoverIndices,indices));
stats.false_positives(linearIndex) = length(setdiff(indices,outlierCrossoverIndices));

end