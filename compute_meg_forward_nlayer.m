function B = compute_meg_forward_nlayer(dipole, sensors, params)
% =========================================================
% Computes MEG forward field using analytical N-layer model
% =========================================================

nMax = params.n_max;
mu0  = params.mu0;

nSensors = length(sensors.r);
B = zeros(nSensors,1);

% Pre-compute coefficients
coeff = compute_coefficients_nlayer(params);

for i = 1:nSensors

    theta = sensors.theta(i);
    phi   = sensors.phi(i);
    r     = sensors.r(i);

    B_i = 0;

    for n = 1:nMax

        % Example term from spherical expansion
        Pn1 = legendre(n, cos(theta));
        Y = sqrt((2*n+1)/(4*pi)) * Pn1(2) * cos(phi);

 3️⃣ `compute_meg_forward_nlayer.m`
**Forward solver (سطح reference)**
```matlab
function B = compute_meg_forward_nlayer(dipole, sensors, params)
% =========================================================
% Computes MEG forward field using analytical N-layer model
% =========================================================

nMax = params.n_max;
mu0  = params.mu0;

nSensors = length(sensors.r);
B = zeros(nSensors,1);

% Pre-compute coefficients
coeff = compute_coefficients_nlayer(params);

for i = 1:nSensors

theta = sensors.theta(i);
phi   = sensors.phi(i);
r     = sensors.r(i);

B_i = 0;

for n = 1:nMax

% Example term from spherical expansion
Pn1 = legendre(n, cos(theta));
Y = sqrt((2*n+1)/(4*pi)) * Pn1(2) * cos(phi);

B_i = B_i + coeff.A(n) * Y / r^(n+1);

end

B(i) = mu0 * B_i;

end

end
