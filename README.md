Smoothing and interpolating with splines under tension, 2018-09-07
===========================

This directory contains all of the analysis files required to reproduce the figures and results in the manuscript.

This a minimized and cleaned-up subset of what was in the repository DrifterTrackSmoothing. That repository contains additional notes and many other experimental methods on this and related techniques.



# Manuscript Figures and Tables

### GenerateAndSaveSyntheticTrajectories.m

Generates synthetic trajectories using the Matern process for their velocity, and then adds Gaussian noise. The trajectories are saved to the sample_data folder.

### MakeFigureBSpline.m

Creates a figure to show an example b-spline and its derivatives at a few orders.

### MakeFigureInterpolation.m

Creates a figure to show how b-splines can be used to interpolate data.

### MakeFigureSignalAndNoiseSpectra.m

Plots a sample of the synthetic data created by MakeAndSaveMaternNoise.m in the time domain, and also plots the power spectra of the signal and the noise.

### MakeTableOptimalParameter.m

Reports RMS errors from the synthetic data. It finds the optimal tension parameter in an unblind search to minimize the RMS error. This is then compared to the blind estimate of the tension parameter, and the iterated blind estimate.

### MakeFigureGaussianFit.m

Makes a figure showing the fit to the real drifter data, assuming a Gaussian distribution for the errors.

### MakeFigureGPSErrorDistribution.m

Plots the error distribution of the motionless GPS datasets, and fits the datasets to a Gaussian and Student-t distribution.

### MakeFigureGPSAutocorrelation.m

Plots the autocorrelation function of the motionless GPS datasets.

# Support

### tcdf.m

Student-t cumulative distribution function


### LoadFigureDefaults.m

A file with the default figure size settings that is referenced by all MakeFigure scripts.
