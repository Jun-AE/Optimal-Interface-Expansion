function y_noisy = nnoise(status, y_clean, n1, n2, n3, n4, seed)
% NNOISE  Four-component pyFBS additive noise model (paper Eq. 17).
%
%   y_noisy = nnoise(status, y_clean, n1, n2, n3, n4)
%   y_noisy = nnoise(status, y_clean, n1, n2, n3, n4, seed)
%
%   Adds independent Gaussian noise to a clean FRF:
%
%       eta     = n1*|Y|.*G1 + 1i*n2*|Y|.*G2 + n3*G3 + 1i*n4*G4
%       y_noisy = y_clean + eta
%
%   Inputs:
%     status   - if false, returns y_clean unchanged
%     y_clean  - clean FRF (any shape)
%     n1, n2   - amplitude-dependent real/imag noise scales
%     n3, n4   - constant real/imag noise-floor magnitudes
%     seed     - optional RNG seed for reproducibility
%
%   References:
%     Bregar, T., Mahmoudi, A., Kodrič, M., Ocepek, D., Trainotti, F.,
%       Göldeli, M., Čepon, G., Boltežar, M., Rixen, D.J. (2022).
%       pyFBS: a python package for frequency based substructuring.
%       Journal of Open Source Software 7, 3399.
%       https://doi.org/10.21105/joss.03399
%
%     Equation reproduced in the form used here:
%       Junaid et al. (2026), Journal of Sound and Vibration, Eq. (17).
%       DOI: 10.1016/j.jsv.2026.001458

if ~status
    y_noisy = y_clean;
    return
end

if nargin >= 7 && ~isempty(seed)
    rng(seed, 'philox');
else
    rng('shuffle', 'philox');
end

absY = abs(y_clean);
g1   = randn(size(y_clean));
g2   = randn(size(y_clean));
g3   = randn(size(y_clean));
g4   = randn(size(y_clean));

eta     = n1 .* absY .* g1 + 1i * n2 .* absY .* g2 ...
        + n3 .* g3          + 1i * n4 .* g4;
y_noisy = y_clean + eta;

end
