% Site 1 or site 2?
Site = 1;


% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

shouldSaveFigures = 0;

% Drifter to highlight in the final plots
choiceDrifter = 6;

if Site == 1
    drifters = load('sample_data/rho1_drifters_projected_ungridded.mat');
else
    drifters = load('sample_data/rho2_drifters_projected_ungridded.mat');
end

