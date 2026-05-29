function c = damping(k, m, cfg)
%  Build a damping matrix using the selected formulation.
%
%   c = damping(k, m, cfg)
%
%   cfg.type    - 'modal' or 'proportional'
%   cfg.ratios  - 1 x Nc damping ratios (modal case)
%   cfg.alpha   - mass-proportional coefficient (Rayleigh case)
%   cfg.beta    - stiffness-proportional coefficient (Rayleigh case)
%
%   Augmented Modal Damping reference:
%     Craig & Kurdila (2006), Fundamentals of Structural Dynamics, §10.3.

switch lower(cfg.type)
    case 'modal'
        c = augmented_modal_damping(k, m, cfg.ratios);
    case 'proportional'
        c = cfg.alpha * m + cfg.beta * k;
end

end

function c = augmented_modal_damping(k, m, damp_ratios)

n_c = length(damp_ratios);
fprintf('\n%d modes damped (augmented modal damping).\n', n_c);

[phi, omega_sq] = eig(k, m);
omega = real(sqrt(diag(omega_sq)));

phi = phi ./ max(abs(phi), [], 1);
m_r = diag(phi.' * m * phi);

zeta        = zeros(size(omega));
zeta(1:n_c) = damp_ratios;

a1       = 2 * zeta(n_c) / omega(n_c);
zeta_hat = zeta - zeta(n_c) * (omega ./ omega(n_c));

c = a1 * k;
for i = 1:(n_c - 1)
    mphi_i = m * phi(:, i);
    c = c + (2 * zeta_hat(i) * omega(i) / m_r(i)) * (mphi_i * mphi_i.');
end

end
