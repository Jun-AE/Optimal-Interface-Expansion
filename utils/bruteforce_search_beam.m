function out = bruteforce_search_beam(interface_dofs, num_sensors, y_n, y_e, opts)
% BRUTEFORCE_SEARCH_BEAM  Exhaustive sensor/excitation search for the best SEMM expansion.
%
%   out = bruteforce_search_beam(interface_dofs, num_sensors, y_n, y_e, opts)
%
%   Inputs:
%     interface_dofs - DoF indices excluded from the candidate set and used
%                      as the validation subset (the "inaccessible" DoFs).
%     num_sensors    - number of sensor DoFs to pick per combination.
%     y_n, y_e       - numerical and experimental FRF matrices (n x n x nFreq).
%     opts           - struct with fields:
%       .frequency_range  rad/s axis used by semm()
%       .method           'COH' or 'LAC'
%       .search_type      'DP' (sensors == excitations) |
%                         'NDP' (excitations from the DoFs minus sensors, with extra count opts.ec) |
%                         'Both' (excitations from all DoFs, may overlap sensors)
%       .truncation       SVD truncation flag for semm
%       .reduction        singular values to DROP (svd_truncation) when truncation is on
%       .ec               extra excitation count (NDP only)
%       .trust_func       sigmoid trust function flag for semm
%       .trust_func_val   trust function parameter
%       .extrema          'max' or 'min'
%
%   Output struct `out`:
%     .sensor_dofs        best sensor combination
%     .excitation_dofs    best excitation combination
%     .interface_cor_oval scalar correlation at the best configuration
%     .interface_cor_2d   2-D correlation matrix on the interface DoFs
%     .cor                full correlation array over all combinations searched
%     .ys                 SEMM-expanded FRF at the best configuration
%
%   Reference:
%     Junaid et al. (2026), Journal of Sound and Vibration.

cand_dofs    = setdiff(1:size(y_n, 1), interface_dofs);
y_e_sub      = y_e(interface_dofs, interface_dofs, :);
sensor_combs = nchoosek(cand_dofs, num_sensors);
n_sensor_combs = size(sensor_combs, 1);

correlation_metric = @(ys) corr_of(ys, y_e_sub, interface_dofs, opts.method);

switch upper(opts.search_type)
    case 'DP'
        cor = run_dp(sensor_combs, y_n, y_e, opts, correlation_metric);
        [best_row, best_col] = pick_best(cor, opts.extrema);
        best_sensors      = sensor_combs(best_row, :);
        best_excitations  = best_sensors;
    case 'NDP'
        exc_combs_per_row = ndp_excitations(sensor_combs, cand_dofs, opts.ec);
        cor = run_outer_inner(sensor_combs, exc_combs_per_row, y_n, y_e, opts, correlation_metric);
        [best_row, best_col] = pick_best(cor, opts.extrema);
        best_sensors      = sensor_combs(best_row, :);
        best_excitations  = exc_combs_per_row(:, :, best_row);
        best_excitations  = best_excitations(best_col, :);
    case 'BOTH'
        exc_combs_all = nchoosek(cand_dofs, num_sensors);
        cor = run_both(sensor_combs, exc_combs_all, y_n, y_e, opts, correlation_metric);
        [best_row, best_col] = pick_best(cor, opts.extrema);
        best_sensors      = sensor_combs(best_row, :);
        best_excitations  = exc_combs_all(best_col, :);
    otherwise
        error('bruteforce_search_beam:UnknownSearchType', ...
              'opts.search_type must be ''DP'', ''NDP'' or ''Both''.');
end

ys_best = semm(best_sensors, best_excitations, y_n, y_e, ...
               opts.frequency_range, opts.truncation, opts.reduction, ...
               opts.trust_func, opts.trust_func_val);

[~, cor_2d, cor_overall] = run_metric(ys_best, y_e, interface_dofs, opts.method);

out.sensor_dofs        = best_sensors;
out.excitation_dofs    = best_excitations;
out.interface_cor_oval = cor_overall;
out.interface_cor_2d   = cor_2d;
out.cor                = cor;
out.ys                 = ys_best;

end


function cor = run_dp(sensor_combs, y_n, y_e, opts, metric)
n = size(sensor_combs, 1);
cor = zeros(n, 1);
q = make_progress_queue(n);

parfor i = 1:n
    s   = sensor_combs(i, :);
    ys  = semm(s, s, y_n, y_e, opts.frequency_range, opts.truncation, ...
               opts.reduction, opts.trust_func, opts.trust_func_val);
    cor(i) = metric(ys);
    send(q, 1);
end
end


function cor = run_outer_inner(sensor_combs, exc_combs_per_row, y_n, y_e, opts, metric)
n_outer = size(sensor_combs, 1);
n_inner = size(exc_combs_per_row, 1);
total   = n_outer * n_inner;
cor_flat = zeros(total, 1);
q = make_progress_queue(total);

parfor k = 1:total
    [i, j] = ind2sub_local(k, n_outer);
    s   = sensor_combs(i, :);
    e   = exc_combs_per_row(j, :, i);
    ys  = semm(s, e, y_n, y_e, opts.frequency_range, opts.truncation, ...
               opts.reduction, opts.trust_func, opts.trust_func_val);
    cor_flat(k) = metric(ys);
    send(q, 1);
end

cor = reshape(cor_flat, n_outer, n_inner);
end


function cor = run_both(sensor_combs, exc_combs_all, y_n, y_e, opts, metric)
n_outer = size(sensor_combs, 1);
n_inner = size(exc_combs_all, 1);
total   = n_outer * n_inner;
cor_flat = zeros(total, 1);
q = make_progress_queue(total);

parfor k = 1:total
    [i, j] = ind2sub_local(k, n_outer);
    s   = sensor_combs(i, :);
    e   = exc_combs_all(j, :);
    ys  = semm(s, e, y_n, y_e, opts.frequency_range, opts.truncation, ...
               opts.reduction, opts.trust_func, opts.trust_func_val);
    cor_flat(k) = metric(ys);
    send(q, 1);
end

cor = reshape(cor_flat, n_outer, n_inner);
end


function exc_combs_per_row = ndp_excitations(sensor_combs, cand_dofs, ec)
n_outer = size(sensor_combs, 1);
k_exc   = size(sensor_combs, 2) + ec;
remaining = length(cand_dofs) - size(sensor_combs, 2);
n_inner = nchoosek(remaining, k_exc);
exc_combs_per_row = zeros(n_inner, k_exc, n_outer);
for i = 1:n_outer
    exc_combs_per_row(:, :, i) = nchoosek(setdiff(cand_dofs, sensor_combs(i, :)), k_exc);
end
end


function v = corr_of(ys, y_e_sub, interface_dofs, method)
ys_sub = ys(interface_dofs, interface_dofs, :);
switch upper(method)
    case 'COH', [~, ~, v] = func_coh(ys_sub, y_e_sub);
    case 'LAC', [~, ~, v] = func_lac(ys_sub, y_e_sub);
end
end


function [v_full, v_2d, v_overall] = run_metric(ys, y_e, interface_dofs, method)
ys_sub  = ys(interface_dofs, interface_dofs, :);
y_e_sub = y_e(interface_dofs, interface_dofs, :);
switch upper(method)
    case 'COH', [v_full, v_2d, v_overall] = func_coh(ys_sub, y_e_sub);
    case 'LAC', [v_full, v_2d, v_overall] = func_lac(ys_sub, y_e_sub);
end
end


function [row, col] = pick_best(cor, extrema)
switch lower(extrema)
    case 'max', [~, idx] = max(cor(:));
    case 'min', [~, idx] = min(cor(:));
end
[row, col] = ind2sub(size(cor), idx);
end


function [i, j] = ind2sub_local(k, n_outer)
i = mod(k - 1, n_outer) + 1;
j = floor((k - 1) / n_outer) + 1;
end


function q = make_progress_queue(total)
q = parallel.pool.DataQueue;
done   = 0;
last_p = -1;
afterEach(q, @(~) tick());
fprintf('bruteforce: %d combinations\n', total);

    function tick()
        done = done + 1;
        p = floor(100 * done / total);
        if p ~= last_p && mod(p, 5) == 0
            fprintf('  %3d%%\n', p);
            last_p = p;
        end
    end
end
