function out = exhaustive_search(y_n, y_e, opts)
% EXHAUSTIVE_SEARCH  Generic optimal sensor/excitation search for SEMM expansion.
%
%   out = exhaustive_search(y_n, y_e, opts)
%
%   Exhaustively scores every sensor (and, for NDP/Both, excitation) DoF
%   combination by the SEMM-expanded FRF correlation on a validation DoF
%   subset, and returns the best configuration per correlation metric. One
%   geometry-agnostic function replaces the per-problem brute-force scripts
%   (beam, square plate, coupled plates, and any new geometry): the caller
%   passes a clean n x n x nFreq FRF and two DoF sets; no geometry lives here.
%
%   The SEMM expansion is delegated entirely to semm.m (the single source of
%   truth). For efficiency the expansion is computed once per combination and
%   scored on every requested metric, and the outer search runs under parfor.
%
%   Inputs:
%     y_n  - (complex, n x n x nFreq) numerical FRF [sensor x excitation x freq].
%     y_e  - (complex, n x n x nFreq) experimental (target) FRF, same size.
%     opts - struct of options (see below).
%
%   opts fields:
%     .num_sensors       (integer, required)  sensors picked per combination.
%     .search_type       'DP' | 'NDP' | 'Both'         (default 'NDP')
%                          DP   : excitations = sensors (drive-point).
%                          NDP  : excitations from candidate DoFs minus the
%                                 chosen sensors, size num_sensors+extra.
%                          Both : excitations from all candidate DoFs
%                                 (may overlap sensors), size num_sensors+extra.
%     .candidate_dofs    (1xC) pool to pick sensors/excitations from. (default 1:n)
%     .validation_dofs   (1xV) subset the metric is scored on (the interface /
%                          inaccessible DoFs).               (default 1:n)
%     .methods           char | cellstr of metrics {'COH','LAC'}.   (default 'COH')
%     .extra_excitations (integer) NDP/Both excitation surplus.      (default 1)
%     .extrema           'max' | 'min' best-pick rule.              (default 'max')
%     .truncation        (logical) SVD-truncate inside semm.        (default false)
%     .reduction         (integer) trailing singular values dropped (svd_truncation).
%                                                                    (default 0)
%     .trust_func        (logical) apply semm sigmoid trust weight. (default false)
%     .trust_func_val    (scalar) trust cut-off frequency [Hz].     (default [])
%     .frequency_range   (1xnFreq) rad/s axis (used by semm; required if trust_func).
%                                                                    (default [])
%     .use_parfor        (logical) parallelise the outer search.    (default true)
%     .verbose           (logical) print progress.                  (default true)
%
%   Output struct array `out` (one element per requested method):
%     .method          - (char) metric name.
%     .sensor_dofs     - (1 x num_sensors) best sensor DoFs.
%     .excitation_dofs - best excitation DoFs.
%     .cor_overall     - (scalar) metric value at the best configuration.
%     .cor_2d          - (V x V) frequency-averaged metric on the validation subset.
%     .cor_landscape   - (n_outer x n_inner) full metric over all combinations.
%     .ys              - (n x n x nFreq) full SEMM expansion at the best config.
%
%   INPUT CONTRACT (read before applying to a new structure):
%     * y_n and y_e must be the SAME size, n x n x nFreq, complex, with the
%       dimension order [response(sensor) x excitation x frequency], and must
%       share the SAME DoF indexing and the SAME frequency axis. y_n is the
%       numerical/clean model; y_e is the (noisy) experimental target.
%     * All DoF references are integer indices into 1:n. candidate_dofs and
%       validation_dofs are arbitrary subsets of 1:n and are INDEPENDENT  they
%       may overlap fully (whole-field search) or be disjoint (interface search):
%         - candidate_dofs  : where sensors/excitations may be placed.
%         - validation_dofs : the interface / inaccessible DoFs the expansion is
%                             scored on (the metric is evaluated on this block).
%     * GEOMETRY IS THE CALLER'S JOB. Any structure-specific DoF preparation
%       (extracting a Z-axis subset, node->DoF mapping, dropping CMS/fixed DoFs,
%       slicing to translational DoFs, etc.) must be done BEFORE calling. This
%       function never sees coordinates or connectivity  only the FRF tensors
%       and the two index sets. That is what makes it geometry-agnostic.
%     * frequency_range is the rad/s axis matching dim 3; it is passed to semm
%       and is REQUIRED only when trust_func is true (otherwise may be []).
%     * Feasibility: num_sensors <= numel(candidate_dofs); for NDP the pool must
%       leave >= num_sensors+extra_excitations DoFs after the sensors; for Both
%       numel(candidate_dofs) >= num_sensors+extra_excitations.
%
%   Examples (one per geometry  only the two index sets change):
%     % Beam (interface = last n_iface translational DoFs):
%     o.num_sensors=2; o.search_type='DP'; o.frequency_range=fr;
%     o.validation_dofs=(n-n_iface+1):n; o.candidate_dofs=setdiff(1:n,o.validation_dofs);
%     res = exhaustive_search(YA, YEA, o);
%
%     % Square plate (whole Z-axis field; caller already sliced to Z DoFs):
%     o.num_sensors=2; o.search_type='NDP'; o.extra_excitations=1; o.frequency_range=fr;
%     o.candidate_dofs=1:n; o.validation_dofs=1:n;
%     res = exhaustive_search(Yz, YEz, o);
%
%     % L-plate / coupled substructure P1 (interface = last 17 DoFs):
%     o.num_sensors=2; o.search_type='NDP'; o.frequency_range=fr;
%     o.candidate_dofs=1:(n-17); o.validation_dofs=(n-16):n;
%     res = exhaustive_search(YzA, YEzA, o);
%
%   Reference:
%     Junaid et al. (2026). Journal of Sound and Vibration.
%     DOI: 10.1016/j.jsv.2026.119782
%
%   See also: semm, func_coh, func_lac, svd_truncation, bruteforce_search_beam.

%% Option defaults

n = size(y_n, 1);

if ~isfield(opts, 'num_sensors') || isempty(opts.num_sensors)
    error('exhaustive_search:MissingNumSensors', 'opts.num_sensors is required.');
end
num_sensors = opts.num_sensors;

search_type      = get_opt(opts, 'search_type',      'NDP');
candidate_dofs   = get_opt(opts, 'candidate_dofs',   1:n);
validation_dofs  = get_opt(opts, 'validation_dofs',  1:n);
methods          = get_opt(opts, 'methods',          {'COH'});
extra_excitations= get_opt(opts, 'extra_excitations',1);
extrema          = get_opt(opts, 'extrema',          'max');
truncation       = get_opt(opts, 'truncation',       false);
reduction        = get_opt(opts, 'reduction',        0);
trust_func       = get_opt(opts, 'trust_func',       false);
trust_func_val   = get_opt(opts, 'trust_func_val',   []);
frequency_range  = get_opt(opts, 'frequency_range',  []);
use_parfor       = get_opt(opts, 'use_parfor',       true);
verbose          = get_opt(opts, 'verbose',          true);

if ischar(methods) || isstring(methods), methods = cellstr(methods); end
n_methods = numel(methods);

candidate_dofs  = candidate_dofs(:).';
validation_dofs = validation_dofs(:).';
I = validation_dofs;

if trust_func && isempty(frequency_range)
    error('exhaustive_search:NoFreq', 'trust_func requires opts.frequency_range.');
end
tf = build_tf(trust_func, trust_func_val);

%% Combination enumeration

n_cand = numel(candidate_dofs);
if num_sensors > n_cand
    error('exhaustive_search:TooManySensors', ...
          'num_sensors (%d) exceeds the candidate pool size (%d).', num_sensors, n_cand);
end

sensor_combs = nchoosek(candidate_dofs, num_sensors);
n_outer      = size(sensor_combs, 1);

switch upper(search_type)
    case 'DP'
        k_exc      = num_sensors;
        n_inner    = 1;
        exc_global = [];
    case 'NDP'
        k_exc      = num_sensors + extra_excitations;
        if (n_cand - num_sensors) < k_exc
            error('exhaustive_search:NDPInfeasible', ...
                  'NDP needs %d excitation DoFs but only %d remain after picking sensors.', ...
                  k_exc, n_cand - num_sensors);
        end
        n_inner    = nchoosek(n_cand - num_sensors, k_exc);
        exc_global = [];
    case 'BOTH'
        k_exc      = num_sensors + extra_excitations;
        if n_cand < k_exc
            error('exhaustive_search:BothInfeasible', ...
                  'Both needs %d excitation DoFs but the candidate pool has %d.', k_exc, n_cand);
        end
        exc_global = nchoosek(candidate_dofs, k_exc);
        n_inner    = size(exc_global, 1);
    otherwise
        error('exhaustive_search:BadType', ...
              'opts.search_type must be ''DP'', ''NDP'' or ''Both''.');
end

if verbose
    fprintf('exhaustive_search [%s]: %d sensor x %d excitation = %d combos, %d method(s)\n', ...
            upper(search_type), n_outer, n_inner, n_outer * n_inner, n_methods);
end

%% Search semm once per combination, scored on every metric

cor_all = zeros(n_outer, n_inner, n_methods);   % (outer x inner x method)

if use_parfor
    parfor i = 1:n_outer
        cor_all(i, :, :) = score_outer(i, sensor_combs, exc_global, y_n, y_e, I, ...
            k_exc, candidate_dofs, search_type, methods, frequency_range, ...
            truncation, reduction, trust_func, tf, n_inner, n_methods);
    end
else
    for i = 1:n_outer
        cor_all(i, :, :) = score_outer(i, sensor_combs, exc_global, y_n, y_e, I, ...
            k_exc, candidate_dofs, search_type, methods, frequency_range, ...
            truncation, reduction, trust_func, tf, n_inner, n_methods);
    end
end

%% Best configuration per method (recompute the winner's full SEMM)

y_e_II = y_e(I, I, :);

out = repmat(struct('method', '', 'sensor_dofs', [], 'excitation_dofs', [], ...
                    'cor_overall', [], 'cor_2d', [], 'cor_landscape', [], 'ys', []), ...
             1, n_methods);

for m = 1:n_methods
    landscape = cor_all(:, :, m);
    [best_row, best_col] = pick_best(landscape, extrema);

    r_best = sensor_combs(best_row, :);
    e_best = excitations_for(best_row, best_col, sensor_combs, exc_global, ...
                             candidate_dofs, k_exc, search_type);

    ys_best = semm(r_best, e_best, y_n, y_e, frequency_range, ...
                   truncation, reduction, trust_func, tf);

    [cor_2d, cor_overall] = metric_full(ys_best(I, I, :), y_e_II, methods{m});

    out(m).method          = methods{m};
    out(m).sensor_dofs     = r_best;
    out(m).excitation_dofs = e_best;
    out(m).cor_overall     = cor_overall;
    out(m).cor_2d          = cor_2d;
    out(m).cor_landscape   = landscape;
    out(m).ys              = ys_best;

    if verbose
        fprintf('  %-5s best %s = %.6f | sensors %s | excitations %s\n', ...
                methods{m}, extrema, cor_overall, ...
                mat2str(r_best), mat2str(e_best));
    end
end

end

% =========================================================================
% Outer-combination scorer: all excitations for one sensor set, all methods.
% =========================================================================
function row_scores = score_outer(i, sensor_combs, exc_global, y_n, y_e, I, ...
    k_exc, candidate_dofs, search_type, methods, frequency_range, ...
    truncation, reduction, trust_func, tf, n_inner, n_methods)

r = sensor_combs(i, :);

switch upper(search_type)
    case 'DP'
        exc_rows = r;                                    % drive-point
    case 'NDP'
        exc_rows = nchoosek(setdiff(candidate_dofs, r), k_exc);
    case 'BOTH'
        exc_rows = exc_global;
end

row_scores = zeros(1, n_inner, n_methods);

for j = 1:size(exc_rows, 1)
    e  = exc_rows(j, :);
    ys = semm(r, e, y_n, y_e, frequency_range, truncation, reduction, trust_func, tf);
    ys_sub = ys(I, I, :);
    for m = 1:n_methods
        row_scores(1, j, m) = metric_overall(ys_sub, y_e(I, I, :), methods{m});
    end
end

end

% =========================================================================
% Helpers
% =========================================================================
function v = get_opt(opts, name, default)
if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = default;
end
end


function e = excitations_for(row, col, sensor_combs, exc_global, candidate_dofs, k_exc, search_type)
r = sensor_combs(row, :);
switch upper(search_type)
    case 'DP'
        e = r;
    case 'NDP'
        exc_rows = nchoosek(setdiff(candidate_dofs, r), k_exc);
        e = exc_rows(col, :);
    case 'BOTH'
        e = exc_global(col, :);
end
end


function tf = build_tf(trust_func, trust_func_val)
if trust_func
    tf.freq_Hz   = trust_func_val;
    tf.steepness = 15;
else
    tf = [];
end
end


function [row, col] = pick_best(landscape, extrema)
switch lower(extrema)
    case 'max', [~, idx] = max(landscape(:));
    case 'min', [~, idx] = min(landscape(:));
    otherwise,  error('exhaustive_search:BadExtrema', 'opts.extrema must be ''max'' or ''min''.');
end
[row, col] = ind2sub(size(landscape), idx);
end


function v = metric_overall(ys_sub, ye_sub, method)
switch upper(method)
    case 'COH', [~, ~, v] = func_coh(ys_sub, ye_sub);
    case 'LAC', [~, ~, v] = func_lac(ys_sub, ye_sub);
    otherwise,  error('exhaustive_search:BadMethod', ...
                      'Unknown method ''%s''. Supported: COH, LAC.', method);
end
end


function [v2d, v] = metric_full(ys_sub, ye_sub, method)
switch upper(method)
    case 'COH', [~, v2d, v] = func_coh(ys_sub, ye_sub);
    case 'LAC', [~, v2d, v] = func_lac(ys_sub, ye_sub);
    otherwise,  error('exhaustive_search:BadMethod', ...
                      'Unknown method ''%s''. Supported: COH, LAC.', method);
end
end
