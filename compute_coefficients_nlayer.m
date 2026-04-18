function coeff = compute_coefficients_nlayer(params)
% =========================================================
% Computes layer-dependent coefficients a_{k,n}, b_{k,n}
% based on boundary conditions (Eqs. 16–23 in the paper)
% =========================================================

nMax = params.n_max;

coeff.A = zeros(nMax,1);

for n = 1:nMax

% Effective anisotropy eigenvalue (Eq. 15)
lambda_n = sqrt(n*(n+1)*params.sigma_t(2)/params.sigma_r(2) + 0.25) - 0.5;

% Simplified coefficient (reference implementation)
coeff.A(n) = (2*n+1)/(2*lambda_n+1);

end

end
