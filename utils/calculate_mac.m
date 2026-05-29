function MAC = calculate_mac(mode_shape_1, mode_shape_2, plotFlag)
%   Modal Assurance Criterion (MAC) between two mode shape sets.
%
%   MAC = calculate_mac(mode_shape_1, mode_shape_2) returns an M1×M2 matrix
%   of MAC values where M1, M2 are the column counts of the two inputs.
%
%   MAC = calculate_mac(mode_shape_1, mode_shape_2, plotFlag) optionally
%   shows an annotated heatmap (cividis colormap, MAC values overlaid).
%   Default PLOTFLAG = false.
%
%   Inputs:
%     mode_shape_1 - (n×M1 double) Mode shape matrix of the first model.
%     mode_shape_2 - (n×M2 double) Mode shape matrix of the second model.
%                   Both must have the same number of rows (DoFs).
%     plotFlag    - (logical) [optional] Show the annotated heatmap.
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
%    Pass 'NumModes' to estimate_modal_parameters to force matching counts.
%
%   Example:
%     mac = calculate_mac(Phi_BeamA, Phi_BeamAExp, true);
%
%   Reference:
%     Allemang, R.J. (2003). The Modal Assurance Criterion - Twenty Years
%     of Use and Abuse. Sound and Vibration 37(8): 14-23.
%
%   See also: estimate_modal_parameters, cmap_cividis.

if nargin < 3 || isempty(plotFlag)
    plotFlag = false;
end

if size(mode_shape_1, 1) ~= size(mode_shape_2, 1)
    error('calculate_mac:DoFMismatch', ...
        'Mode shape matrices must have the same number of rows (DoFs).');
end

%% MAC matrix

MAC = abs(mode_shape_1' * mode_shape_2).^2 ./ (diag(mode_shape_1' * mode_shape_1) * diag(mode_shape_2' * mode_shape_2)');

%% Optional annotated heatmap

if plotFlag
    drawMACHeatmap(macMatrix);
end

end

% =========================================================================
% Local helper: annotated MAC heatmap
% =========================================================================

[M1, M2] = size(M);

figure('Name', 'Modal Assurance Criterion (MAC)', 'Color', 'w');
imagesc(M, [0 1]);
colormap(cmap_cividis(256));
cb = colorbar;
cb.Label.String = 'MAC';
set(gca, 'XTick', 1:M2, 'YTick', 1:M1, 'TickLength', [0 0]);
xlabel('Mode index (set 2)');
ylabel('Mode index (set 1)');
title('Modal Assurance Criterion');
axis equal tight;

% Cell-wise text labels — white text on dark cells, black on bright.
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
