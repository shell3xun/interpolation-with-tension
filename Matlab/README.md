Smoothing and interpolating with splines under tension, 2019-03-28
===========================

This directory contains all of the analysis files required to reproduce the figures and results in the manuscript.

This a minimized and cleaned-up subset of what was in the repository DrifterTrackSmoothing. That repository contains additional notes and many other experimental methods on this and related techniques.



# Manuscript figures and tables

## Figures

### MakeFigureInterpolation.m

Figure 1. Creates a figure to show how b-splines can be used to interpolate data.

### MakeFigureBSpline.m

Figure 2. Creates a figure to show an example b-spline and its derivatives at a few orders.

### MakeFigureSignalSpectraWithSmootherOfSplineDegrees.m

Figure 3. Creates a figure to how interpolating splines of different degrees are able to fit to an uncontaminated Matern signal.

### MakeFigureSignalAndNoiseSpectraWithSmoother.m

Figure 4 Shows the three fits to a noise Matern signal, with differing effective sample sizes (strides).

### MakeFigureGammaDependence.m

Figure 5. Shows the relationship between effective sample size and u_rms.

### MakeFigureGPSErrorDistribution.m

Figure 6. Plots the error distribution of the motionless GPS datasets, and fits the datasets to a Gaussian and Student-t distribution.

### MakeFigureGPSAutocorrelation.m

Figure 7. Plots the autocorrelation function of the motionless GPS datasets.

## Tables

### MakeTableOptimalTensionOrder.m

Reports RMS errors from the synthetic data. It finds the optimal tension parameter in an unblind search to minimize the RMS error. This is then compared to the blind estimate of the tension parameter, and the iterated blind estimate.

### MakeTableOptimalParameter.m

Reports RMS errors from the synthetic data. It finds the optimal tension parameter in an unblind search to minimize the RMS error. This is then compared to the blind estimate of the tension parameter, and the iterated blind estimate.


# Other figures and tables

### MakeFigureGaussianFit.m

Makes a figure showing the fit to the real drifter data, assuming a Gaussian distribution for the errors.

### MakeFigureSignalAndNoiseSpectra.m

Plots a sample of the synthetic data created by MakeAndSaveMaternNoise.m in the time domain, and also plots the power spectra of the signal and the noise.

### MakeTableRobustAccelerationEstimate.m and CheckNoiseThresholds.m

Used to make an estimated RMS for some derivative of the signal, given the known noise level. 

### MakeTableOutliersFullTensionTest.m

What method is best for getting to full tension?


# Support

### GenerateAndSaveSyntheticTrajectories.m

Generates synthetic trajectories using the Matern process for their velocity, and then adds Gaussian noise. The trajectories are saved to the sample_data folder.

### tcdf.m

Student-t cumulative distribution function


### LoadFigureDefaults.m

A file with the default figure size settings that is referenced by all MakeFigure scripts.
