function plot_frf(frequency_range_Hz, response_array, excitation_array, FRF_1, FRF_2, str_1, str_2, varargin)
%  Magnitude / phase / coherence comparison plot for 2-4 FRFs at one (r, e) pair.
%
%   PLOT_FRF(F_HZ, R, E, FRF_1, FRF_2, STR_1, STR_2)
%   PLOT_FRF(..., FRF_3, STR_3)
%   PLOT_FRF(..., FRF_3, STR_3, FRF_4, STR_4)
%
%   Three subplots stacked: |Y(r,e,f)|, angle(Y), coherence(FRF_1, FRF_2)
%   (plus pairwise coherences against FRF_3, FRF_4 if supplied).
%
%   Inputs:
%     F_HZ    - (double, 1×nFreq) [Hz] Frequency axis.
%     R       - (integer, scalar) Response DoF index into the FRF.
%     E       - (integer, scalar) Excitation DoF index into the FRF.
%     FRF_k   - (complex, n×n×nFreq) FRF matrix for the k-th curve.
%     STR_k   - (char) Legend label for the k-th curve.
%
%   Colour palette is sampled from cmap_magma (interior trim).
%
%   Reference:
%     Junaid et al. (2026). Journal of Sound and Vibration.
%     https://doi.org/10.1016/j.jsv.2026.119782
%
%   See also: func_coh, func_lac, cmap_magma, compute_frf.

%% Optional FRF_3 / FRF_4

FRF_extra = {};
str_extra = {};
if numel(varargin) >= 2
    FRF_extra{end+1} = varargin{1};
    str_extra{end+1} = varargin{2};
end
if numel(varargin) >= 4
    FRF_extra{end+1} = varargin{3};
    str_extra{end+1} = varargin{4};
end

nCurves = 2 + numel(FRF_extra);

%% Palette — magma sampled at nCurves interior points

mcmap = cmap_magma(nCurves + 2);
mcmap = mcmap(2:end-1, :);

%% Extract scalar FRF traces at (r, e)

p   = response_array;
q   = excitation_array;
om  = frequency_range_Hz;
nf  = length(om);

Y_all   = cell(1, nCurves);
str_all = [{str_1, str_2}, str_extra];
Y_all{1} = squeeze(FRF_1(p, q, 1:nf));
Y_all{2} = squeeze(FRF_2(p, q, 1:nf));
for k = 1:numel(FRF_extra)
    Y_all{2 + k} = squeeze(FRF_extra{k}(p, q, 1:nf));
end

%% Figure with three stacked tiles

fig = figure('Color', 'w', 'Renderer', 'opengl', 'Position', [716, 82, 864, 782]);
tlo = tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- Magnitude ---
ax1 = nexttile(tlo);
for k = 1:nCurves
    semilogy(ax1, om, abs(Y_all{k}), 'LineWidth', 2, 'Color', mcmap(k, :), ...
             'DisplayName', sprintf('$%s$', str_all{k}));
    hold(ax1, 'on');
end
ylabel(ax1, 'Magnitude (g/N)', 'Interpreter', 'latex');
legend(ax1, 'Location', 'southwest', 'Interpreter', 'latex');
grid(ax1, 'on'); ax1.GridAlpha = 0.25;
xlim(ax1, [om(1), om(end)]);
hold(ax1, 'off');

% --- Phase ---
ax2 = nexttile(tlo);
for k = 1:nCurves
    plot(ax2, om, rad2deg(angle(Y_all{k})), 'LineWidth', 1.8, 'Color', mcmap(k, :));
    hold(ax2, 'on');
end
plot(ax2, om, zeros(size(om)), '--', 'Color', [0 0 0 0.4], 'HandleVisibility', 'off');
ylabel(ax2, 'Phase ($^{\circ}$)', 'Interpreter', 'latex');
grid(ax2, 'on'); ax2.GridAlpha = 0.25;
xlim(ax2, [om(1), om(end)]);  ylim(ax2, [-200, 200]);
hold(ax2, 'off');

% --- Coherence (FRF_1 vs every other curve) ---
ax3 = nexttile(tlo);
hold(ax3, 'on');
for k = 2:nCurves
    [COH, ~, ~] = func_coh(FRF_1, getFRF(k, FRF_2, FRF_extra));
    coh         = squeeze(COH(p, q, 1:nf));

    % Translucent area under the curve, then the line on top.
    area(ax3, om, coh, ...
         'FaceColor', mcmap(k, :), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(ax3, om, coh, 'LineWidth', 1.6, 'Color', mcmap(k, :));
end
plot(ax3, om, 0.5 * ones(size(om)), '--', 'Color', [0 0 0 0.4]);
ylabel(ax3, 'Coherence',         'Interpreter', 'latex');
xlabel(ax3, 'Frequency (Hz)',    'Interpreter', 'latex');
grid(ax3, 'on'); ax3.GridAlpha = 0.25;
xlim(ax3, [om(1), om(end)]);  ylim(ax3, [-0.05, 1.05]);
hold(ax3, 'off');

end

% =========================================================================
% Local helper: pick the k-th FRF (k = 2..nCurves) from FRF_2 + extras
% =========================================================================
function F = getFRF(k, FRF_2, FRF_extra)
if k == 2
    F = FRF_2;
else
    F = FRF_extra{k - 2};
end
end
