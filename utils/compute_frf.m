function y = compute_frf(omega, k, m, c, frf_type)
% COMPUTE_FRF  Full FRF matrix via direct dynamic-stiffness inversion.
%
%   y = compute_frf(omega, k, m, c, frf_type)
%
%   omega    - 1 x nFreq angular frequency axis (rad/s)
%   k, m, c  - stiffness, mass, damping matrices (n x n)
%   frf_type - 'receptance' | 'mobility' | 'accelerance'
%
%   Returns y of size n x n x nFreq.
%  Only direct dynamic-stiffness inversion is implemented for now. 
%  Modal superposition method will be added to this function as an optional method (suitable for larger models).

n       = size(k, 1);
n_freq  = numel(omega);
y       = complex(zeros(n, n, n_freq));
I       = eye(n);

% parfor can be switched to the serial for loop if parallel toolbox is not available.
parfor j = 1:n_freq
    om_j     = omega(j);
    z_j      = k + (1i * om_j) * c - (om_j^2) * m;
    y(:,:,j) = z_j \ I;
end

switch lower(frf_type)
    case 'receptance'
        % y already holds receptance.
    case 'mobility'
        scale = reshape(1i * omega, 1, 1, []);
        y     = y .* scale;
    case 'accelerance'
        scale = reshape(-(omega.^2), 1, 1, []);
        y     = y .* scale;
end

end
