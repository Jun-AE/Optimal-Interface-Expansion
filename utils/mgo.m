function [best_f, best_x, cnvg, stats] = mgo(n, max_iter, lb_arg, ub_arg, dim, fobj, plot_flag, save_vid, solution_space, opts)
% MGO  Mountain Gazelle Optimizer (integer-bounded, vectorised batch form).
%
%   Baseline implementation with dedup + cache accelerations on by default.
%   For deterministic integer-bounded objectives (the project's GCCM search),
%   repeated candidates are evaluated once: duplicates within an iteration are
%   collapsed and results are memoised across iterations. The RNG sequence is
%   identical to the plain batched MGO, so optima are unchanged.
%
%   [best_f, best_x, cnvg, stats] = mgo(n, max_iter, lb, ub, dim, fobj)
%   [best_f, best_x, cnvg, stats] = mgo(..., plot_flag, save_vid, solution_space, opts)
%
%   Inputs:
%     n              - population size
%     max_iter       - iteration count
%     lb, ub         - bounds (scalar or 1xdim)
%     dim            - search-space dimensionality
%     fobj           - objective handle, fobj(x) returns scalar to minimise
%     plot_flag      - show live 2D progress plot (default false)
%     save_vid       - write progress animation to mgo_progress.mp4 (default false)
%     solution_space - precomputed fitness landscape for the contour overlay
%     opts (struct, all optional):
%       .use_dedup   - unique candidates per iteration (default true)
%       .use_cache   - cache fobj by rounded x (default true)
%       .use_parfor  - parfor for fobj batches (default false)
%
%   Outputs:
%     best_f - best objective value found
%     best_x - best position (1 x dim)
%     cnvg   - 1 x max_iter best-so-far convergence trace
%     stats  - struct: n_fobj_eval, n_cache_hit, n_dedup_saved
%
%   Reference:
%     Abdollahzadeh et al. (2022). Mountain gazelle optimizer.
%     Advances in Engineering Software 174, 103282.

if nargin < 7 || isempty(plot_flag), plot_flag = false; end
if nargin < 8 || isempty(save_vid),  save_vid  = false; end
if nargin < 9, solution_space = []; end
if nargin < 10 || isempty(opts), opts = struct(); end
if save_vid, plot_flag = true; end

use_dedup  = get_opt(opts, 'use_dedup',  true);
use_cache  = get_opt(opts, 'use_cache',  true);
use_parfor = get_opt(opts, 'use_parfor', false);

if use_parfor && isempty(gcp('nocreate'))
    use_parfor = false;
end

if isscalar(lb_arg), lb = repmat(lb_arg, 1, dim); else, lb = lb_arg(:).'; end
if isscalar(ub_arg), ub = repmat(ub_arg, 1, dim); else, ub = ub_arg(:).'; end
span = ub - lb;

stats = struct('n_fobj_eval', 0, 'n_cache_hit', 0, 'n_dedup_saved', 0);
cache = containers.Map('KeyType', 'char', 'ValueType', 'double');

%% Init

x    = round(rand(n, dim) .* span + lb);
[cost, stats] = evaluate_rows(x, fobj, cache, use_cache, use_parfor, stats);

[best_f, idx] = min(cost);
best_x = x(idx, :);
cnvg   = zeros(1, max_iter);

sub_n   = ceil(n / 3);
pool_sz = 4 * n;
cand    = zeros(pool_sz, dim);

%% Main loop

for iter = 1:max_iter

    m_all = zeros(n, dim);
    for i = 1:n
        sub = randperm(n, sub_n);
        m_all(i, :) = x(randi([sub_n, n]), :) * floor(rand * 2) ...
                    + mean(x(sub, :), 1)      * ceil(rand * 2);
    end

    decay = exp(2 - 2 * iter / max_iter);
    a_vec = randn(n, dim) .* decay;
    d_vec = (abs(x) + abs(best_x)) .* (2 * rand - 1);

    a2 = -1 - iter / max_iter;
    u  = randn(n, dim);
    v  = randn(n, dim);
    cofi = cat(3, ...
        rand(n, dim), ...
        ((a2 + 1) + rand(n, 1)) .* ones(1, dim), ...
        a2 .* randn(n, dim), ...
        u .* v.^2 .* cos((rand * 2) .* u));

    cofi_tsm  = pick_layer(cofi, randi(4, n, 1));
    cofi_mh_a = pick_layer(cofi, randi(4, n, 1));
    cofi_mh_b = pick_layer(cofi, randi(4, n, 1));
    cofi_bmh  = pick_layer(cofi, randi(4, n, 1));

    f1 = randi(2, n, dim);  f2 = randi(2, n, dim);
    f3 = randi(2, n, dim);  f4 = randi(2, n, dim);
    f5 = randi(2, n, dim);  f6 = randi(2, n, dim);
    partner_mh = randi(n, n, 1);

    c_msf = lb + span .* rand(n, dim);
    c_tsm = best_x - abs((f1 .* m_all - f2 .* x) .* a_vec) .* cofi_tsm;
    c_mh  = (m_all + cofi_mh_a) + (f3 .* best_x - f4 .* x(partner_mh, :)) .* cofi_mh_b;
    c_bmh = (x - d_vec) + (f5 .* best_x - f6 .* m_all) .* cofi_bmh;

    cand(1:n, :)       = c_msf;
    cand(n+1:2*n, :)   = c_tsm;
    cand(2*n+1:3*n, :) = c_mh;
    cand(3*n+1:4*n, :) = c_bmh;
    cand(:)            = round(max(min(cand, ub), lb));

    [new_cost, stats] = evaluate_rows(cand, fobj, cache, use_cache, use_parfor, stats, use_dedup);

    [all_cost, order] = sort([cost; new_cost]);
    all_x = [x; cand];
    x     = all_x(order(1:n), :);
    cost  = all_cost(1:n);

    if cost(1) < best_f
        best_f = cost(1);
        best_x = x(1, :);
    end
    cnvg(iter) = best_f;

    if mod(iter, 10) == 0
        fprintf('MGO  iter %4d/%d   best = %.6g   x* = %s\n', ...
            iter, max_iter, best_f, mat2str(best_x));
    end

    if plot_flag
        plot_progress(x, best_x, cnvg, lb, ub, iter, max_iter, save_vid, solution_space);
    end
end

end

%% -------------------------------------------------------------------------
function v = get_opt(opts, name, default)
if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = default;
end
end


function [costs, stats] = evaluate_rows(rows, fobj, cache, use_cache, use_parfor, stats, use_dedup)
if nargin < 7, use_dedup = false; end

n_rows = size(rows, 1);
if use_dedup
    [u_rows, ~, ic] = unique(rows, 'rows', 'stable');
    stats.n_dedup_saved = stats.n_dedup_saved + (n_rows - size(u_rows, 1));
else
    u_rows = rows;
    ic = (1:n_rows).';
end

[u_cost, stats] = evaluate_unique(u_rows, fobj, cache, use_cache, use_parfor, stats);
costs = u_cost(ic);
end


function [u_cost, stats] = evaluate_unique(u_rows, fobj, cache, use_cache, use_parfor, stats)
n_unique = size(u_rows, 1);
u_cost   = nan(n_unique, 1);
todo     = true(n_unique, 1);

if use_cache
    for k = 1:n_unique
        key = cache_key(u_rows(k, :));
        if isKey(cache, key)
            u_cost(k) = cache(key);
            todo(k)   = false;
            stats.n_cache_hit = stats.n_cache_hit + 1;
        end
    end
end

idx_eval = find(todo);
n_eval   = numel(idx_eval);

if n_eval == 0
    return;
end

if use_parfor
    chunk = nan(n_eval, 1);
    parfor j = 1:n_eval
        chunk(j) = fobj(u_rows(idx_eval(j), :));
    end
    for j = 1:n_eval
        k = idx_eval(j);
        u_cost(k) = chunk(j);
        if use_cache
            cache(cache_key(u_rows(k, :))) = u_cost(k);
        end
    end
else
    for j = 1:n_eval
        k = idx_eval(j);
        u_cost(k) = fobj(u_rows(k, :));
        if use_cache
            cache(cache_key(u_rows(k, :))) = u_cost(k);
        end
    end
end

stats.n_fobj_eval = stats.n_fobj_eval + n_eval;
end


function key = cache_key(x)
x = round(x(:)).';
key = sprintf('%d,', x);
end


function out = pick_layer(stack, picks)
n   = size(stack, 1);
dim = size(stack, 2);
out = zeros(n, dim);
for k = 1:n
    out(k, :) = stack(k, :, picks(k));
end
end


function plot_progress(x, best_x, cnvg, lb, ub, iter, max_iter, save_vid, solution_space)
persistent fig ax1 ax2 h_gaz h_best h_conv vid mcmap

if length(lb) ~= 2, return; end

if isempty(fig) || ~isvalid(fig)
    fig = figure('Name', 'mgo progress', 'Color', 'w', 'Position', [100 100 1200 500]);
    tlo = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    mcmap = local_magma(5);

    ax1 = nexttile(tlo, 1);
    hold(ax1, 'on');
    if ~isempty(solution_space)
        [r, c]   = size(solution_space);
        [xg, yg] = meshgrid(linspace(lb(1), ub(1), c), linspace(lb(2), ub(2), r));
        contourf(ax1, xg, yg, solution_space, 15, 'LineStyle', 'none', 'FaceAlpha', 0.85);
        colormap(ax1, local_cividis(256));
        cb = colorbar(ax1); cb.Label.String = 'Fitness';
    end
    h_gaz  = scatter(ax1, x(:, 1), x(:, 2), 70, mcmap(2, :), '^', 'filled', ...
        'MarkerEdgeColor', 'k', 'DisplayName', 'Herd');
    h_best = scatter(ax1, best_x(1), best_x(2), 280, mcmap(4, :), 'p', 'filled', ...
        'MarkerEdgeColor', 'k', 'DisplayName', 'Best');
    xlabel(ax1, 'x_1'); ylabel(ax1, 'x_2');
    title(ax1, 'Search-space exploration');
    legend(ax1, [h_gaz, h_best], 'Location', 'best');
    xlim(ax1, [lb(1), ub(1)]); ylim(ax1, [lb(2), ub(2)]);
    grid(ax1, 'on'); ax1.GridAlpha = 0.25;

    ax2 = nexttile(tlo, 2);
    hold(ax2, 'on'); grid(ax2, 'on'); ax2.GridAlpha = 0.25;
    h_conv = plot(ax2, 1:iter, cnvg(1:iter), 'LineWidth', 2, 'Color', mcmap(4, :));
    xlabel(ax2, 'Iteration'); ylabel(ax2, 'Best fitness');
    title(ax2, 'Convergence');
    xlim(ax2, [1, max_iter]);

    if save_vid
        vid           = VideoWriter('mgo_progress.mp4', 'MPEG-4');
        vid.FrameRate = 15; vid.Quality   = 100;
        open(vid);
    end
else
    set(h_gaz,  'XData', x(:, 1), 'YData', x(:, 2));
    set(h_best, 'XData', best_x(1), 'YData', best_x(2));
    set(h_conv, 'XData', 1:iter,    'YData', cnvg(1:iter));
    if iter > 1
        y_min = min(cnvg(1:iter));
        y_max = max(cnvg(1:iter));
        if y_min == y_max, y_max = y_max + eps(y_max) + 1e-12; end
        ylim(ax2, [y_min, y_max]);
    end
end

sgtitle(sprintf('mgo   iter %d/%d   best = %.6g', iter, max_iter, cnvg(iter)), ...
    'FontWeight', 'bold');
drawnow limitrate;

if ~isempty(vid)
    writeVideo(vid, getframe(fig));
    if iter == max_iter, close(vid); vid = []; end
end
end


function m = local_cividis(n)
ctrl = [0.000 0.135 0.305 ; 0.082 0.249 0.452 ; 0.346 0.408 0.450 ;
        0.555 0.560 0.451 ; 0.793 0.717 0.354 ; 1.000 0.948 0.038];
m = interp_cmap(ctrl, n);
end

function m = local_magma(n)
ctrl = [0.001 0.000 0.014 ; 0.207 0.072 0.388 ; 0.448 0.088 0.453 ;
        0.685 0.114 0.443 ; 0.881 0.207 0.367 ; 0.990 0.380 0.359 ;
        0.987 0.991 0.750];
m = interp_cmap(ctrl, n);
end

function m = interp_cmap(ctrl, n)
t_in  = linspace(0, 1, size(ctrl, 1));
t_out = linspace(0, 1, n);
m = zeros(n, 3);
for k = 1:3
    m(:, k) = interp1(t_in, ctrl(:, k), t_out, 'pchip');
end
m = max(0, min(1, m));
end
