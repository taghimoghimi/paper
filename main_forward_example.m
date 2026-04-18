% =========================================================
% Example script for MEG forward computation
% Analytical N-layer anisotropic eccentric spherical model
% =========================================================

clear; clc;

% Load model parameters
params = model_parameters_example();

% Define dipole (example used in the paper)
dipole.position = [0, 0, 0.068];   % [m]
dipole.moment   = [1, 0, 0];       % Qx dipole [nAm]

% Generate simple sensor configuration (spherical shell)
nSensors = 306;
theta = linspace(0, pi, nSensors);
phi   = linspace(0, 2*pi, nSensors);
[sTheta, sPhi] = meshgrid(theta, phi);

sensors.r = params.sensor_radius * ones(numel(sTheta),1);
sensors.theta = sTheta(:);
sensors.phi   = sPhi(:);

% Compute MEG forward field
B = compute_meg_forward_nlayer(dipole, sensors, params);

% Plot magnitude of magnetic field
figure;
plot(abs(B), 'LineWidth', 1.5);
xlabel('Sensor index');
ylabel('|B| (T)');
title('MEG Forward Field – N-layer Anisotropic Eccentric Model');
grid on;
