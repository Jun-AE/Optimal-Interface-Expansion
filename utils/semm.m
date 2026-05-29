function y_semm = semm(r, e, y_n, y_e, frequency_range, truncation_flag, reduction_num, trust_func_flag, trust_func_params)
% SEMM  Expand FRFs to unmeasured DoFs via System Equivalent Model Mixing.
%
%   y_semm = semm(r, e, y_n, y_e)
%   y_semm = semm(r, e, y_n, y_e, frequency_range, truncation_flag, ...
%                 reduction_num, trust_func_flag, trust_func_params)
%
%   Inputs:
%     r, e               - response and excitation DoF indices
%     y_n, y_e           - numerical and experimental FRFs (n x n x nFreq)
%     frequency_range    - rad/s axis (only used if trust_func_flag is true)
%     truncation_flag    - apply SVD truncation before pinv (default false)
%     reduction_num      - singular values to DROP per page (svd_truncation) when truncation_flag is true
%     trust_func_flag    - apply sigmoid blend along frequency (default false)
%     trust_func_params  - struct with .freq_Hz and .steepness
%
%   Reference:
%     Junaid et al. (2026), Journal of Sound and Vibration. DOI 10.1016/j.jsv.2026.001458

if nargin < 5 || isempty(frequency_range),   frequency_range   = [];    end
if nargin < 6 || isempty(truncation_flag),   truncation_flag   = false; end
if nargin < 7 || isempty(reduction_num),     reduction_num     = 0;     end
if nargin < 8 || isempty(trust_func_flag),   trust_func_flag   = false; end
if nargin < 9 || isempty(trust_func_params), trust_func_params = [];    end

r = round(r);
e = round(e);

y_e_re = y_e(r, e, :);
y_n_re = y_n(r, e, :);

switch truncation_flag
    case true
        y_n_rn = svd_truncation(y_n(r, :, :), reduction_num);
        y_n_ne = svd_truncation(y_n(:, e, :), reduction_num);
    case false
        y_n_rn = y_n(r, :, :);
        y_n_ne = y_n(:, e, :);
end

proj_l   = pagemtimes(y_n, pagepinv(y_n_rn));
delta_re = y_n_re - y_e_re;
proj_r   = pagemtimes(pagepinv(y_n_ne), y_n);
y_r      = pagemtimes(pagemtimes(proj_l, delta_re), proj_r);

switch trust_func_flag
    case true
        frequency_hz = frequency_range(:) / (2 * pi);
        f_c = trust_func_params.freq_Hz;
        k   = trust_func_params.steepness / f_c;
        w   = 1 ./ (1 + exp(-k .* (frequency_hz - f_c)));
        y_r = y_r .* reshape(w, 1, 1, []);
end

y_semm = y_n - y_r;

end
