function [naturalFrequencies, modeShapes, CMIF] = estimate_modal_parameters(FRF, f, varargin)
% ESTIMATE_MODAL_PARAMETERS  CMIF-based modal parameter extraction.
%
%   [FN, PHI, CMIF] = ESTIMATE_MODAL_PARAMETERS(FRF, F) returns natural
%   frequencies FN, mode shapes PHI, and the Complex Mode Indicator
%   Function CMIF computed from a 3D FRF array.
%
%   [FN, PHI, CMIF] = ESTIMATE_MODAL_PARAMETERS(FRF, F, NAME, VALUE, ...)
%   accepts the following Name-Value options:
%     'NumModes'            - (integer | []) If supplied, return exactly
%                             this many modes (the top-N peaks by
%                             prominence, re-sorted by frequency). Default
%                             [] means "all peaks above threshold". When
%                             comparing two FRFs via calculate_mac, set this
%                             to a fixed value to guarantee an N×N MAC.
%     'ProminenceThreshold' - (double) Min peak prominence. Default 0.1.
%     'MinPeakDistanceHz'   - (double) [Hz] Min peak separation. Default 5.
%     'Plot'                - (logical) Show a CMIF diagnostic figure.
%                             Default false.
%
%   Inputs:
%     FRF - (complex, n×n×nFreq) FRF matrix [response × excitation × freq].
%     F   - (double, 1×nFreq)    [Hz] Frequency axis.
%
%   Outputs:
%     FN   - (M×1 double) [Hz] Natural frequencies (M = number of modes
%            kept; equals 'NumModes' when supplied).
%     PHI  - (n×M double) Mode shapes — first left singular vector of FRF
%            at each peak frequency, sign-aligned to the dominant DoF.
%     CMIF - (nFreq×min(n,n) double) Singular values of FRF(:,:,k) at each
%            frequency, sorted descending.
%
%   Method:
%     For each frequency k, SVD of FRF(:,:,k) yields singular values
%     CMIF(k,:). Peaks of CMIF(:,1) (the leading singular value) identify
%     resonances. Mode shapes are the first left singular vector at each
%     peak; signs are flipped so the dominant entry has positive real part.
%
%   Example:
%     % Force exactly 3 modes so MACs are 3×3
%     [fn, Phi] = estimate_modal_parameters(YA, freq_Hz, 'NumModes', 3);
%
%   Reference:
%     Junaid et al. (2026). Journal of Sound and Vibration.
%     DOI: 10.1016/j.jsv.2026.001458
%
%   See also: calculate_mac, findpeaks, svd.

%% Parse options

p = inputParser;
addParameter(p, 'NumModes',            [],   @(x) isempty(x) || (isscalar(x) && x>0));
addParameter(p, 'ProminenceThreshold', 0.1,  @isnumeric);
addParameter(p, 'MinPeakDistanceHz',   5,    @isnumeric);
addParameter(p, 'Plot',                false, @islogical);
parse(p, varargin{:});

numModesTarget   = p.Results.NumModes;
prominenceThresh = p.Results.ProminenceThreshold;
minPeakDistHz    = p.Results.MinPeakDistanceHz;
plotFlag         = p.Results.Plot;

% Convert min peak distance from Hz to sample indices.
df                  = f(2) - f(1);
minPeakDistance_idx = max(1, round(minPeakDistHz / df));
maxAllowed          = floor(length(f) / 2);
if minPeakDistance_idx > maxAllowed
    warning('MinPeakDistance reduced from %d to %d samples', ...
            minPeakDistance_idx, maxAllowed);
    minPeakDistance_idx = maxAllowed;
end

%% CMIF — leading singular value of FRF at every frequency

[nSensors, nExcit, nFreq] = size(FRF);
CMIF = zeros(nFreq, min(nSensors, nExcit));
for k = 1:nFreq
    [~, S, ~]   = svd(FRF(:, :, k), 'econ');
    CMIF(k, :)  = diag(S);
end

%% Peak picking — optional NumModes truncation

if isempty(numModesTarget)
    [~, peakIdx] = findpeaks(CMIF(:, 1), ...
                             'MinPeakProminence', prominenceThresh, ...
                             'MinPeakDistance',   minPeakDistance_idx);
else
    % Top-N by prominence so two FRFs yield matching mode counts; then
    % re-sort by frequency for natural ascending order.
    [~, peakIdx] = findpeaks(CMIF(:, 1), ...
                             'MinPeakProminence', prominenceThresh, ...
                             'MinPeakDistance',   minPeakDistance_idx, ...
                             'NPeaks',  numModesTarget, ...
                             'SortStr', 'descend');
    peakIdx = sort(peakIdx);
end

naturalFrequencies = f(peakIdx);

%% Mode shapes — first left singular vector at each peak, sign-aligned

modeShapes = zeros(nSensors, length(peakIdx));
for m = 1:length(peakIdx)
    [U, ~, ~]          = svd(FRF(:, :, peakIdx(m)), 'econ');
    modeShapes(:, m)   = U(:, 1);

    [~, idxDom]        = max(abs(modeShapes(:, m)));
    modeShapes(:, m)   = modeShapes(:, m) * sign(real(modeShapes(idxDom, m)));
end

%% Optional CMIF diagnostic figure

if plotFlag
    drawCMIFDiagnostic(f, CMIF, peakIdx, naturalFrequencies);
end

end

% =========================================================================
% Local helper: CMIF magnitude + phase diagnostic with picked-peak markers
% =========================================================================
function drawCMIFDiagnostic(f, CMIF, peakIdx, fn)

fig  = figure('Name', 'CMIF — Modal Parameter Identification', 'Color', 'w');
tlo  = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

% Magnitude
ax1  = nexttile(tlo, 1);
mag  = 20 * log10(abs(CMIF(:, 1)));
plot(ax1, f, mag, 'LineWidth', 1.6, 'Color', [0.20 0.20 0.30]);
hold(ax1, 'on');
scatter(ax1, fn, mag(peakIdx), 80, 'v', 'filled', ...
        'MarkerFaceColor', [0.85 0.30 0.10], 'MarkerEdgeColor', 'k');
xlabel(ax1, 'Frequency (Hz)');  ylabel(ax1, '|CMIF_1| (dB)');
title(ax1, 'Complex Mode Indicator Function — magnitude');
grid(ax1, 'on'); ax1.GridAlpha = 0.25;
hold(ax1, 'off');

% Phase
ax2  = nexttile(tlo, 2);
phs  = unwrap(angle(CMIF(:, 1)));
plot(ax2, f, phs, 'LineWidth', 1.6, 'Color', [0.20 0.20 0.30]);
hold(ax2, 'on');
scatter(ax2, fn, phs(peakIdx), 80, 'v', 'filled', ...
        'MarkerFaceColor', [0.85 0.30 0.10], 'MarkerEdgeColor', 'k');
xlabel(ax2, 'Frequency (Hz)');  ylabel(ax2, 'arg(CMIF_1) (rad)');
title(ax2, 'Complex Mode Indicator Function — phase');
grid(ax2, 'on'); ax2.GridAlpha = 0.25;
hold(ax2, 'off');

end
