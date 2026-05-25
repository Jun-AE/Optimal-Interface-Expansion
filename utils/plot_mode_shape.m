function [V, freq] = plot_mode_shape(vec, fsol, n_start, n_end)
% PLOT_MODE_SHAPE  3D stacked plot of normalised translational mode shapes.
%
%   [V, FREQ] = PLOT_MODE_SHAPE(VEC, FSOL, N_START, N_END) extracts the
%   translational DoF rows from the eigenvector matrix VEC, normalises each
%   mode by its peak amplitude, prepends/appends boundary zeros, and plots
%   modes N_START through N_END as 3D curves with one mode per "row" along y.
%
%   Inputs:
%     VEC     - (ndof×ndof double) Eigenvectors from EIG(K, M).
%     FSOL    - (ndof×ndof double) Diagonal eigenvalue matrix from EIG.
%     N_START - (integer, scalar) First mode index to plot.
%     N_END   - (integer, scalar) Last mode index to plot.
%
%   Outputs:
%     V    - ((nTrans+2)×ndof double) Normalised mode shape matrix with
%            zero boundary rows prepended and appended.
%     FREQ - (ndof×1 double) [Hz] Sorted natural frequencies.
%
%   Mode line colours are sampled from cmap_magma (interior trim).
%
%   Reference:
%     Junaid et al. (2026). Journal of Sound and Vibration.
%     DOI: 10.1016/j.jsv.2026.001458
%
%   See also: cmap_magma, visualize_mode_shapes, plot_modal_waterfall.

%% Natural frequencies and translational mode shapes

fsol = diag(fsol);
fhz  = sqrt(fsol) / (2 * pi);
freq = sort(fhz);
fprintf('\nNatural frequencies (Hz) for modes %d..%d:\n', n_start, n_end);
disp(freq(n_start:n_end));

% Keep translational DoF rows only (every other entry in the 2-DoF/node beam).
v = vec(1:2:end, :);

% Pad with zero "boundary" rows top and bottom; normalise each column.
V                 = zeros(size(v, 1) + 2, size(v, 2));
modeAmp           = max(abs(v), [], 1);
modeAmp(modeAmp == 0) = 1;
V(2:end-1, :)     = v ./ modeAmp;

%% Magma palette for the mode curves

nModes = n_end - n_start + 1;
mcmap  = cmap_magma(nModes + 2);
mcmap  = mcmap(2:end-1, :);

%% Plot

figure('Name', 'Mode Shapes', 'Color', 'w');
hold on;

x = 0:(size(v, 1) + 1);             % node positions including boundaries
y = zeros(numel(x), 1);

for n = n_start:n_end
    c = mcmap(n - n_start + 1, :);
    plot3(x, y + 0.25 * (n - n_start), V(:, n), '-o', ...
          'MarkerFaceColor', c, 'MarkerEdgeColor', c, ...
          'Color',           c,  'LineWidth', 2, ...
          'DisplayName', sprintf('Mode %d : %.2f Hz', n, freq(n)));
end

% Zero-line in the y = 0 plane for reference.
plot3(x, y, zeros(size(x)), '--', 'Color', [0 0 0 0.4], 'HandleVisibility', 'off');

hold off;
axis tight;
xlabel('Node index'); ylabel('Mode'); zlabel('Normalised amplitude');
title('Mode shapes');
grid on; ax = gca; ax.GridAlpha = 0.25;
legend('Location', 'southeast');
view(0, 180);

end
