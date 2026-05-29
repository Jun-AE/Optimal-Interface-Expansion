function correlation_plot(CorMatrix, caption)
% correlation_plot  Annotated heatmap of an FRF correlation matrix.
%
%   correlation_plot(CORMATRIX) shows a cividis heatmap of the input
%   correlation matrix with min, max and mean annotations, defaulting the
%   metric label to 'CCM'
%
%   CORRELATION_PLOT(CORMATRIX, CAPTION) overrides the metric label.
%
%   Inputs:
%     CORMATRIX - (N×N double) Correlation matrix in [0, 1] (COH, LAC, ...).
%     CAPTION   - (char) [optional] Metric name shown in the title.
%                 Default: 'CCM'
%
%   Reference:
%     Junaid et al. (2026). Journal of Sound and Vibration.
%     DOI: https://doi.org/10.1016/j.jsv.2026.119782
%
%   See also: func_coh, func_lac, cmap_cividis.

if nargin < 2 || isempty(caption)
    caption = 'CCM';
end

%% Summary statistics

avgVal             = mean(CorMatrix, 'all');
[maxVal, maxIdx]   = max(CorMatrix(:));
[minVal, minIdx]   = min(CorMatrix(:));
[maxRow, maxCol]   = ind2sub(size(CorMatrix), maxIdx);
[minRow, minCol]   = ind2sub(size(CorMatrix), minIdx);

%% Heatmap

figure('Name', caption, 'Color', 'w');
imagesc(CorMatrix, [0 1]);
colormap(cmap_cividis(256));
cb              = colorbar;
cb.Label.String = caption;

set(gca, 'YDir', 'normal', 'TickLength', [0 0]);
axis equal tight;
xlabel('Y_{i,j}');  ylabel('Y_{i,j}');
title(sprintf('%s  —  \\Gamma = %.5f', caption, avgVal));

%% Annotations

hold on;
text(maxCol, maxRow, sprintf('Max: %.3f', maxVal), ...
     'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
     'BackgroundColor', [1 1 1 0.7]);
text(minCol, minRow, sprintf('Min: %.3f', minVal), ...
     'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
     'BackgroundColor', [0 0 0 0.5]);
hold off;

end
