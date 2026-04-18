function params = model_parameters_example()
% Model parameters used in the paper (example values)

params.mu0 = 4*pi*1e-7;

% Number of layers
params.N = 5;

% Radii of layers [m]
params.radius = [0.065 0.070 0.075 0.085 0.090];

% Radial and tangential conductivities [S/m]
params.sigma_r = [0.14 0.33 1.79 0.01 0.43];
params.sigma_t = [0.14 0.33 1.79 0.015 0.43];

% Eccentricity parameter (g in paper)
params.g = 0.0025;  % [m]

% Spherical harmonic truncation
params.n_max = 40;

% Sensor radius
params.sensor_radius = 0.09;

end
