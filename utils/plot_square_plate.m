function plot_square_plate(node_coords, sensor_dofs, excitation_dofs, show_indices)
% Plot sensor / excitation placement on the square plate.
%
%   plot_square_plate(node_coords, sensor_dofs, excitation_dofs)
%   plot_square_plate(node_coords, sensor_dofs, excitation_dofs, show_indices)
%
%   Inputs:
%     node_coords     - (nNode x >=2) node coordinate table as returned by
%                       load_hcb_model (columns 1,2 are the in-plane x, z used
%                       for the plate face). Row i corresponds to measured DoF i.
%     sensor_dofs     - measured-DoF indices chosen as sensors.
%     excitation_dofs - measured-DoF indices chosen as excitations.
%     show_indices    - (logical) [optional] annotate node numbers. Default true.
%
%   See also: load_hcb_model, exhaustive_search, cmap_magma.

if nargin < 4 || isempty(show_indices), show_indices = true; end

x = node_coords(:, 1);
y = node_coords(:, 2);

mcmap    = cmap_magma(6);
c_sensor = mcmap(2, :);
c_excit  = mcmap(4, :);

figure('Color', 'w');
ax = axes;
hold(ax, 'on');

% All nodes.
scatter(ax, x, y, 60, [0.6 0.6 0.6], 'o', 'filled', 'DisplayName', 'Node');

% Plate boundary (1 m x 1 m, centred at the origin).
plot(ax, [-0.5 0.5 0.5 -0.5 -0.5], [-0.5 -0.5 0.5 0.5 -0.5], ...
     'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

sensor_dofs     = sensor_dofs(:).';
excitation_dofs = excitation_dofs(:).';

if ~isempty(sensor_dofs)
    scatter(ax, x(sensor_dofs), y(sensor_dofs), 200, c_sensor, 's', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'Sensor');
end
if ~isempty(excitation_dofs)
    scatter(ax, x(excitation_dofs), y(excitation_dofs), 200, c_excit, '^', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'Excitation');
end

if show_indices
    for i = 1:size(node_coords, 1)
        text(ax, x(i), y(i) - 0.035, num2str(i), 'FontSize', 8, ...
             'Color', [0.3 0.3 0.3], 'HorizontalAlignment', 'center', ...
             'HandleVisibility', 'off');
    end
end

hold(ax, 'off');
axis(ax, 'equal');
xlim(ax, [-0.75 0.75]);
ylim(ax, [-0.75 0.75]);
grid(ax, 'on'); ax.GridAlpha = 0.25;
xlabel(ax, 'x (m)');
ylabel(ax, 'z (m)');
title(ax, 'Square plate - sensor & excitation placement');
legend(ax, 'Location', 'bestoutside');

end
