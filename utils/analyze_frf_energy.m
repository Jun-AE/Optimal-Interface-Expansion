function E = analyze_frf_energy(Y, freq_Hz, topN, nModes, plotFlag, dB_floor)
% ANALYZE_FRF_ENERGY  Sensor / excitation / modal energy diagnostics on an FRF.
%
%   E = ANALYZE_FRF_ENERGY(Y, FREQ_HZ, TOPN, NMODES) computes per-sensor and
%   per-excitation FRF energy, identifies the TOPN highest-energy DoFs of
%   each, detects NMODES modal peaks from the frequency-averaged energy
%   profile, and lists the TOPN DoFs that participate most strongly in each
%   mode. A single diagnostic figure is produced (dB-scaled imagesc of the
%   row-wise energy density with overlaid top-DoF markers, plus a bar chart
%   of the per-DoF energy split).
%
%   E = ANALYZE_FRF_ENERGY(Y, FREQ_HZ, TOPN, NMODES, PLOTFLAG) suppresses the
%   diagnostic figure when PLOTFLAG is false (default true).
%
%   E = ANALYZE_FRF_ENERGY(Y, FREQ_HZ, TOPN, NMODES, PLOTFLAG, DB_FLOOR) sets
%   the lower clipping limit (in dB relative to the peak) for the heatmap.
%   Default: -60 dB. Lower values reveal more low-amplitude structure;
%   higher values (e.g. -30) emphasise only the strong modal ridges.
%
%   Inputs:
%     Y         - (complex, n×n×nFreq) FRF matrix [sensor × excitation × freq].
%     FREQ_HZ   - (double, 1×nFreq) [Hz] Frequency axis matching dim 3 of Y.
%     TOPN      - (integer, scalar) Number of top sensors / excitations / per-mode
%                 DoFs to report. Pass 0 to skip top-DoF selection.
%     NMODES    - (integer, scalar) Number of modal peaks to detect.
%     PLOTFLAG  - (logical) [optional] Show the diagnostic figure. Default: true.
%     DB_FLOOR  - (double, scalar) [optional] Lower clim of the dB heatmap.
%                 Default: -60 (dB rel. peak).
%
%   Output struct E with fields:
%     .sensor_energy     - (n×1) sum_{j,k} |Y(i,j,k)|^2          (per-sensor)
%     .excitation_energy - (n×1) sum_{i,k} |Y(i,j,k)|^2          (per-excitation)
%     .top_sensors       - (TOPN×1) indices of highest-energy sensors.
%     .top_excitations   - (TOPN×1) indices of highest-energy excitations.
%     .lambda_rows       - (n×nFreq) sqrt(sum_j |Y(i,j,k)|^2)    (row-wise energy density)
%     .lambda_cols       - (n×nFreq) sqrt(sum_i |Y(i,j,k)|^2)    (col-wise energy density)
%     .mode_freqs        - (NMODES×1) [Hz] frequencies of detected modal peaks.
%     .mode_dofs         - (TOPN×NMODES) top DoFs for each detected mode.
%     .common_dofs       - (variable) DoFs appearing in every mode's top-N list.
%
%   Notes:
%     This function consolidates three earlier helpers (compute_energy_importance,
%     SensorEnergy, getTopModalDoFs) that each recomputed the same row-wise
%     energy density and produced separate figures. All energy sums are
%     vectorised via sum(abs(Y).^2, [dim_list]) — no per-frequency loops.
%
%     Modal peak detection uses findpeaks (Signal Processing Toolbox) on the
%     mean of lambda_rows across sensors; peaks are returned in ascending
%     frequency order regardless of amplitude.
%
%   Example:
%     E = analyze_frf_energy(YA, beamA.frequency_range_Hz, 4, 3);
%     plot_beam(3, beamA.nnode-1, E.top_sensors, E.top_excitations, 'A');
%
%   Reference:
%     Junaid et al. (2026). Journal of Sound and Vibration.
%     DOI: 10.1016/j.jsv.2026.001458
%
%   See also: findpeaks, cmap_cividis, cmap_magma, plot_beam.

if nargin < 5 || isempty(plotFlag), plotFlag = true;  end
if nargin < 6 || isempty(dB_floor), dB_floor = -60;   end

%% Per-sensor / per-excitation total energy (squared)

absY2 = abs(Y).^2;                                  % [n × n × nFreq]
E.sensor_energy     = sum(absY2, [2 3]);            % [n × 1]  sum over (j, k)
E.excitation_energy = squeeze(sum(absY2, [1 3]));   % [n × 1]  sum over (i, k)

[~, sortS]        = sort(E.sensor_energy,     'descend');
[~, sortX]        = sort(E.excitation_energy, 'descend');
E.top_sensors     = sortS(1:min(topN, numel(sortS)));
E.top_excitations = sortX(1:min(topN, numel(sortX)));

%% Row-wise and column-wise energy density (kept for the heatmap)

E.lambda_rows = sqrt(squeeze(sum(absY2, 2)));       % [n × nFreq]
E.lambda_cols = sqrt(squeeze(sum(absY2, 1)));       % [n × nFreq]

%% Modal peak detection from the frequency-averaged row energy

avgResponse  = mean(E.lambda_rows, 1);              % [1 × nFreq]
[~, peakIdx] = findpeaks(avgResponse, ...
                         'NPeaks',  nModes, ...
                         'SortStr', 'descend');     % strongest peaks first
peakIdx      = sort(peakIdx);                       % then re-sort by frequency
E.mode_freqs = freq_Hz(peakIdx).';                  % [nModes × 1]

%% Top-N DoFs per mode and common-across-modes set

E.mode_dofs = zeros(topN, nModes);
for k = 1:nModes
    [~, sortedDofs]  = sort(E.lambda_rows(:, peakIdx(k)), 'descend');
    E.mode_dofs(:,k) = sortedDofs(1:topN);
end

E.common_dofs = E.mode_dofs(:, 1);
for k = 2:nModes
    E.common_dofs = intersect(E.common_dofs, E.mode_dofs(:, k));
end

%% Console summary

fprintf('\n--- FRF Energy Analysis ---\n');
for k = 1:nModes
    fprintf('Mode %d (%.2f Hz): Top DoFs = %s\n', ...
            k, E.mode_freqs(k), num2str(E.mode_dofs(:,k).'));
end
if isempty(E.common_dofs)
    fprintf('No common DoFs across all modes.\n');
else
    fprintf('Common DoFs: %s\n', strjoin(string(E.common_dofs), ', '));
end

%% Diagnostic figure

if plotFlag
    drawEnergyFigure(E, freq_Hz, topN, nModes, dB_floor);
end

end

% =========================================================================
% Local helper: combined diagnostic figure (heatmap + bar chart)
% =========================================================================
function drawEnergyFigure(E, freq_Hz, topN, nModes, dB_floor)

fig = figure('Position', [200 150 1100 750]);
tlo = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

%% Tile 1 — row-wise energy density heatmap (dB, normalised) with top-DoF markers

ax1   = nexttile(tlo, 1);
nRows = size(E.lambda_rows, 1);

% Convert to dB relative to the global peak — modal ridges and the
% inter-modal valleys then live in a comparable dynamic range, instead of
% the peaks saturating the colourmap and burying everything else in dark.
peak       = max(E.lambda_rows(:));
lambda_dB  = 20 * log10(E.lambda_rows / peak + eps);

imagesc(ax1, freq_Hz, 1:nRows, lambda_dB);
set(ax1, 'YDir', 'normal');
colormap(ax1, cmap_cividis(256));
clim(ax1, [dB_floor, 0]);
cb = colorbar(ax1);
cb.Label.String = '\Lambda_{row}(f)  (dB rel. peak)';
xlabel(ax1, 'Frequency (Hz)');
ylabel(ax1, 'DoF index');
title(ax1, 'Row-wise FRF energy density');
hold(ax1, 'on');

% Magma palette for the top-DoF markers (interior trim per skill convention).
mcmap = cmap_magma(max(topN, 3) + 2);
mcmap = mcmap(2:end-1, :);

for i = 1:numel(E.top_sensors)
    sIdx        = E.top_sensors(i);
    [~, fIdx]   = max(E.lambda_rows(sIdx, :));
    scatter(ax1, freq_Hz(fIdx), sIdx, 280, mcmap(i, :), 's', 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.2, 'MarkerFaceAlpha', 0.95, ...
            'DisplayName', sprintf('Top sensor %d', sIdx));
end

for i = 1:numel(E.top_excitations)
    eIdx        = E.top_excitations(i);
    [~, fIdx]   = max(E.lambda_cols(eIdx, :));
    scatter(ax1, freq_Hz(fIdx), eIdx, 220, mcmap(i, :), '^', 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.2, 'MarkerFaceAlpha', 0.95, ...
            'DisplayName', sprintf('Top excitation %d', eIdx));
end

% Vertical reference lines at each detected modal frequency.
for k = 1:nModes
    xline(ax1, E.mode_freqs(k), '--', sprintf('f_{%d} = %.1f Hz', k, E.mode_freqs(k)), ...
          'Color', [0.9 0.9 0.9], 'LineWidth', 1.0, 'Alpha', 0.6, ...
          'LabelHorizontalAlignment', 'left', ...
          'LabelVerticalAlignment',   'bottom');
end

grid(ax1, 'on'); ax1.GridAlpha = 0.25;
hold(ax1, 'off');

%% Tile 2 — per-DoF energy split (sensor vs excitation, normalised)

ax2   = nexttile(tlo, 2);
total = sum(E.sensor_energy) + sum(E.excitation_energy);
se_pct = 100 * E.sensor_energy     / total;
ee_pct = 100 * E.excitation_energy / total;

b = bar(ax2, [se_pct(:), ee_pct(:)], 'grouped', 'BarWidth', 1, 'EdgeColor', 'none');
b(1).FaceColor = mcmap(1, :);     b(1).FaceAlpha = 0.85;
b(2).FaceColor = mcmap(end, :);   b(2).FaceAlpha = 0.85;
xlabel(ax2, 'DoF index');
ylabel(ax2, 'Energy share (% of total)');
title(ax2, 'Per-DoF energy contribution');
legend(ax2, {'Sensor (response)', 'Excitation'}, 'Location', 'best');
grid(ax2, 'on'); ax2.GridAlpha = 0.25;

end
