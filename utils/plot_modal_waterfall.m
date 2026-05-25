function [fn, dr, mode_shapes_scaled] = plot_modal_waterfall(Y, excitation_dof, freq_Hz, fs, sensor_dofs, num_modes, mode_shapes_FE, fn_FE)
% PLOT_MODAL_WATERFALL  3D waterfall of imag(FRF) with sign-aligned mode shapes overlaid.
%
%   [FN, DR, MS] = PLOT_MODAL_WATERFALL(Y, EXCITATION_DOF, FREQ_HZ, FS,
%                                     SENSOR_DOFS, NUM_MODES,
%                                     MODE_SHAPES_FE, FN_FE)
%
%   Plots a 3D waterfall of the imaginary part of the FRF (highlighting
%   modal peaks) and overlays the mode shapes at each resonance as 3D
%   curves. The waterfall surface uses the cividis colormap; the mode
%   curves are sampled from magma for a coherent, perceptually-uniform
%   look that prints/scales correctly and is CVD-safe.
%
%   The function takes the FRF in its natural n×n×nFreq layout and slices
%   one excitation column internally — no caller-side permute/reshape
%   gymnastics required.
%
%   Inputs:
%     Y               - (complex, n×n×nFreq) Full FRF matrix [response × excitation × freq].
%     EXCITATION_DOF  - (integer, scalar) Index into the excitation dimension of Y.
%     FREQ_HZ         - (double, 1×nFreq) [Hz] Frequency axis.
%     FS              - (double, scalar) [Hz] Sampling frequency (used by MODALFIT
%                       when FE mode shapes are not supplied).
%     SENSOR_DOFS     - (integer vector) Response DoFs to display on the y-axis
%                       (typically the translational DoFs).
%     NUM_MODES       - (integer, scalar) Number of modes to extract and overlay.
%     MODE_SHAPES_FE  - (double, n×K) [optional] FE-derived eigenvector matrix.
%                       When supplied, MODALFIT is bypassed and these are used
%                       directly. First NUM_MODES columns are taken.
%     FN_FE           - (double, K×1) [optional] FE-derived natural frequencies [Hz].
%                       Required together with MODE_SHAPES_FE.
%
%   Outputs:
%     FN - (NUM_MODES×1 double) [Hz] Natural frequencies used.
%     DR - (NUM_MODES×1 double) Modal damping ratios. NaN(num_modes,1) when
%          FE mode shapes were supplied (no MODALFIT call performed).
%     MS - (n×NUM_MODES double) Sign-aligned, per-mode peak-scaled mode shapes.
%          A new array — the input MODE_SHAPES_FE is NOT mutated.
%
%   Scaling and sign convention:
%     - Sign alignment: each mode is flipped if needed so that real(MS(idx, k))
%       is non-negative at the DoF idx of largest |MS|. Replaces the previous
%       hardcoded "flip mode 1 and 3" special case (which was correct only for
%       the specific beam in the original demo).
%     - Amplitude: each mode shape is scaled so its peak |MS| equals the peak
%       |imag(FRF)| at that mode's resonance frequency. Uniform across all
%       modes — replaces the previous magic "× 40 boost on mode 1".
%
%   Example:
%     [fn, dr, ms] = plot_modal_waterfall( ...
%         YEA, 10, beamA.frequency_range_Hz, 1800, ...
%         1:2:beamA.ndof, 3, ...
%         vecA_exp(:, 1:3), freqA_exp(1:3));
%
%   Reference:
%     Junaid et al. (2026). Journal of Sound and Vibration.
%     DOI: 10.1016/j.jsv.2026.001458
%
%   See also: modalfit, cmap_magma, cmap_cividis, compute_frf.

%% Slice FRF for the chosen excitation DoF

% Y(:, excitation_dof, :) -> [sensor × 1 × freq] -> [freq × sensor]
FRF = squeeze(Y(:, excitation_dof, :)).';

% Defensive: drop sensor indices that fall outside the supplied FRF. This
% guards against the caller passing full-system DoF indices when Y has
% been reduced upstream (e.g. re-running a section after Y was sliced down
% to its translational DoFs elsewhere in the workflow).
nSensor     = size(FRF, 2);
sensor_dofs = sensor_dofs(sensor_dofs >= 1 & sensor_dofs <= nSensor);

%% Modal parameters: use FE-supplied or call MODALFIT

if nargin >= 8 && ~isempty(fn_FE) && ~isempty(mode_shapes_FE)
    fn              = fn_FE(1:num_modes);
    dr              = NaN(num_modes, 1);
    mode_shapes_raw = mode_shapes_FE(:, 1:num_modes);
else
    fmax = freq_Hz(end);
    fn   = modalfit(FRF, freq_Hz, fs, num_modes, ...
                    'FreqRange', [0 ceil(fmax)], 'FitMethod', 'lsrf');
    fn   = fn(1:num_modes);
    [~, dr, mode_shapes_raw] = modalfit(FRF, freq_Hz, fs, num_modes, ...
                    'PhysFreq', fn, 'FreqRange', [0 ceil(fmax)], 'FitMethod', 'lsrf');
end

%% Sign-align each mode to the excitation DoF

% At resonance omega_k, imag(Y_pq(omega_k)) is proportional to phi_k(p)*phi_k(q),
% so the waterfall ridge across sensors q carries the sign of phi_k(p). Aligning
% each mode shape to be non-negative at the excitation DoF p ensures the overlay
% z(q) = phi_k_aligned(q) tracks the waterfall ridge instead of its mirror.
if excitation_dof >= 1 && excitation_dof <= size(mode_shapes_raw, 1)
    for k = 1:num_modes
        if real(mode_shapes_raw(excitation_dof, k)) < 0
            mode_shapes_raw(:, k) = -mode_shapes_raw(:, k);
        end
    end
end

%% Uniform amplitude scaling across modes

% Global imag(FRF) peak across every resonance is used as a common scale for
% every mode. Per-mode peak matching would shrink low-frequency modes to near
% invisibility for accelerance/mobility (where peak heights vary by orders of
% magnitude between modes); using a single global peak keeps all modes equally
% readable without ad-hoc per-mode boost factors.
global_peak = 0;
for k = 1:num_modes
    [~, idx_freq] = min(abs(freq_Hz - fn(k)));
    global_peak   = max(global_peak, max(abs(imag(FRF(idx_freq, :)))));
end

mode_shapes_scaled = mode_shapes_raw;       % new array — do not mutate input
for k = 1:num_modes
    mode_amp = max(abs(mode_shapes_raw(:, k)));
    if mode_amp > 0
        mode_shapes_scaled(:, k) = mode_shapes_raw(:, k) / mode_amp * global_peak;
    end
end

%% Plot

figure;
ax = axes;
hold(ax, 'on');

drawImagWaterfall(ax, FRF, freq_Hz, sensor_dofs);
drawModeOverlay (ax, mode_shapes_scaled, sensor_dofs, fn, num_modes);
finaliseAxes    (ax);

end

% =========================================================================
% Local helpers
% =========================================================================

function drawImagWaterfall(ax, FRF, freq_Hz, sensor_dofs)
% Waterfall of imag(FRF) along (freq, sensor). Cividis colormap, alpha 0.45.

Z = imag(FRF(:, sensor_dofs)).';       % [sensor × freq] as waterfall expects
h = waterfall(ax, freq_Hz, sensor_dofs, Z);
alpha(h, 0.45);

colormap(ax, cmap_cividis(256));
end

function drawModeOverlay(ax, mode_shapes_scaled, sensor_dofs, fn, num_modes)
% Mode shape curves at resonance. Colours sampled from magma, avoiding the
% extreme dark and bright ends so each line stays distinguishable.

cmap = cmap_magma(num_modes + 2);
cmap = cmap(2:end-1, :);

for k = 1:num_modes
    msk_sub = mode_shapes_scaled(sensor_dofs, k);
    z       = abs(msk_sub) .* sign(real(msk_sub));
    plot3(ax, fn(k) * ones(size(sensor_dofs)), sensor_dofs, z, ...
          'LineWidth',   2.5, ...
          'Color',       cmap(k, :), ...
          'DisplayName', sprintf('Mode %d : %.2f Hz', k, fn(k)));
end
end

function finaliseAxes(ax)
xlabel(ax, 'Frequency (Hz)');
ylabel(ax, 'Sensor DoF');
zlabel(ax, 'imag(FRF) / Mode-shape amplitude');
set(ax, 'YDir', 'reverse');
grid(ax, 'on');
ax.GridAlpha = 0.25;
legend(ax, 'show', 'Location', 'best');
view(ax, 30, 30);
end
