function plot_beam(n_interface, n_node, sensor_dofs, excitation_dofs, side)
% PLOT_BEAM  Schematic of a cantilever beam with sensor / excitation markers.
%
%   plot_beam(n_interface, n_node, sensor_dofs, excitation_dofs, side)
%
%   n_interface     - number of inaccessible interface nodes at the free end
%   n_node          - total physical node count (including the fixed node)
%   sensor_dofs     - sensor node indices in the free-node numbering
%   excitation_dofs - excitation node indices in the free-node numbering
%   side            - 'A' (fixed at left) | 'B' (fixed at right)

side_char = upper(char(string(side)));
side_char = side_char(1);

n_free  = n_node - 1;
x_free  = 1:n_free;

switch side_char
    case 'A'
        x_boundary  = 0;
        inacc_nodes = (n_free - n_interface + 1):n_free;
    case 'B'
        x_boundary  = n_node;
        inacc_nodes = 1:n_interface;
    otherwise
        x_boundary  = 0;
        inacc_nodes = (n_free - n_interface + 1):n_free;
        side_char   = 'A';
end

sensor_dofs     = sensor_dofs(:).';
excitation_dofs = excitation_dofs(:).';
sensor_dofs     = sensor_dofs(    sensor_dofs    >= 1 & sensor_dofs    <= n_free);
excitation_dofs = excitation_dofs(excitation_dofs >= 1 & excitation_dofs <= n_free);

mcmap          = cmap_magma(6);
color_sensor   = mcmap(2, :);
color_excit    = mcmap(4, :);
color_inacc    = [0.85 0.15 0.15];
color_node     = [0.30 0.30 0.30];
color_boundary = [0.10 0.10 0.10];

figure('Color', 'w');
ax = axes;
hold(ax, 'on');

% Beam line spans the full extent including the boundary node.
beam_extent = [min(x_boundary, 1), max(x_boundary, n_free)];
plot(ax, beam_extent, [0 0], 'k', 'LineWidth', 2.5, 'HandleVisibility', 'off');

handles = gobjects(0);

h = scatter(ax, x_free, zeros(size(x_free)), 60, color_node, 'o', 'filled', ...
            'DisplayName', 'Beam node');
handles(end+1) = h;

if ~isempty(inacc_nodes)
    h = scatter(ax, inacc_nodes, zeros(size(inacc_nodes)), 220, color_inacc, 'x', ...
                'LineWidth', 3, 'DisplayName', 'Inaccessible DoF');
    handles(end+1) = h;
end

if ~isempty(sensor_dofs)
    h = scatter(ax, sensor_dofs, zeros(size(sensor_dofs)), 240, color_sensor, 's', 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'Sensor');
    handles(end+1) = h;
end

if ~isempty(excitation_dofs)
    h = scatter(ax, excitation_dofs, zeros(size(excitation_dofs)), 240, color_excit, '^', 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.8, 'DisplayName', 'Excitation');
    handles(end+1) = h;
end

h = scatter(ax, x_boundary, 0, 260, color_boundary, 's', 'filled', ...
            'DisplayName', 'Fixed boundary');
handles(end+1) = h;

hold(ax, 'off');

xlim(ax, [min(x_boundary, 1) - 1, max(x_boundary, n_free) + 1]);
ylim(ax, [-0.5, 0.5]);
yticks(ax, []);
ax.YColor = 'none';

xlabel(ax, 'Node index');
title(ax, sprintf('Beam %s  —  sensor & excitation layout', side_char));
grid(ax, 'on'); ax.GridAlpha = 0.25;
legend(ax, handles, 'Location', 'best', 'FontSize', 10);

end
