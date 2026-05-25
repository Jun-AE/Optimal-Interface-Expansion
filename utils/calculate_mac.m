function macMatrix = calculate_mac(modeShapes1, modeShapes2, plotFlag)
% CALCULATE_MAC  Modal Assurance Criterion (MAC) between two mode shape sets.
%
%   MAC = CALCULATE_MAC(MODESHAPES1, MODESHAPES2) returns an M1×M2 matrix
%   of MAC values where M1, M2 are the column counts of the two inputs.
%
%   MAC = CALCULATE_MAC(MODESHAPES1, MODESHAPES2, PLOTFLAG) optionally
%   shows an annotated heatmap (cividis colormap, MAC values overlaid).
%   Default PLOTFLAG = false.
%
%   Inputs:
%     MODESHAPES1 - (n×M1 double) Mode shape matrix of the first model.
%     MODESHAPES2 - (n×M2 double) Mode shape matrix of the second model.
%                   Both must have the same number of rows (DoFs).
%     PLOTFLAG    - (logical) [optional] Show the annotated heatmap.
%                   Default: false.
%
%   Output:
%     MAC - (M1×M2 double) MAC matrix. Element (i, j) is:
%
%             |phi1_i' * phi2_j|^2
%        MAC = -----------------------------------------
%             (phi1_i' * phi1_i) * (phi2_j' * phi2_j)
%
%        with values in [0, 1] (1 = perfect correlation, 0 = orthogonal).
%
%   Asymmetric MAC matrices:
%     If M1 ≠ M2 the matrix is rectangular — the inputs simply had
%     different mode counts. This is mathematically valid but usually
%     indicates that the upstream modal-extraction step returned different
%     numbers of peaks for the two models. Pass 'NumModes' to
%     ESTIMATEMODALPARAMETERS to force matching counts.
%
%   Example:
%     mac = calculate_mac(Phi_BeamA, Phi_BeamAExp, true);
%
%   Reference:
%     Junaid et al. (2026). Journal of Sound and Vibration.
%     DOI: 10.1016/j.jsv.2026.001458
%     Allemang, R.J. (2003). The Modal Assurance Criterion - Twenty Years
%     of Use and Abuse. Sound and Vibration 37(8): 14-23.
%
%   See also: estimate_modal_parameters, cmap_cividis.

if nargin < 3 || isempty(plotFlag)
    plotFlag = false;
end

if size(modeShapes1, 1) ~= size(modeShapes2, 1)
    error('calculate_mac:DoFMismatch', ...
        'Mode shape matrices must have the same number of rows (DoFs).');
end

%% MAC matrix

macMatrix = abs(modeShapes1' * modeShapes2).^2 ./ ...
            (diag(modeShapes1' * modeShapes1) * ...
             diag(modeShapes2' * modeShapes2)');

%% Optional annotated heatmap

if plotFlag
    drawMACHeatmap(macMatrix);
end

end

% =========================================================================
% Local helper: annotated MAC heatmap (cividis + value labels in each cell)
% =========================================================================
function drawMACHeatmap(M)

[M1, M2] = size(M);

figure('Name', 'Modal Assurance Criterion (MAC)', 'Color', 'w');
imagesc(M, [0 1]);
colormap(cmap_cividis(256));
cb = colorbar;
cb.Label.String = 'MAC';

% Keep the default imagesc YDir='reverse' so matrix row 1 sits at the top
% and the leading diagonal runs from the top-left to the bottom-right.
set(gca, 'XTick', 1:M2, 'YTick', 1:M1, 'TickLength', [0 0]);
xlabel('Mode index (set 2)');
ylabel('Mode index (set 1)');
title('Modal Assurance Criterion');
axis equal tight;

% Cell-wise text labels — white text on dark cells, black on bright.
% Threshold at the midpoint of the cividis luminance ramp.
for i = 1:M1
    for j = 1:M2
        val = M(i, j);
        if val > 0.5
            txtColor = [0 0 0];
        else
            txtColor = [1 1 1];
        end
        text(j, i, sprintf('%.2f', val), ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment',   'middle', ...
             'Color',               txtColor, ...
             'FontWeight',          'bold', ...
             'FontSize',            12);
    end
end

end
