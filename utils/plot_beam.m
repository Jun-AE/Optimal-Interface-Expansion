function plot_beam(n_interface, n_node, sensor_dofs, excitation_dofs, side)
%  Schematic of a cantilever beam with sensor / excitation markers.
%
%   plot_beam(n_interface, n_node, sensor_dofs, excitation_dofs, side)
%
%   n_interface     - number of inaccessible interface nodes at the free end
%   n_node          - total physical node count (including the fixed node)
%   sensor_dofs     - sensor node indices in the free-node numbering
%   excitation_dofs - excitation node indices in the free-node numbering
%   side            - 'A' (fixed at left) | 'B' (fixed at right)
%


side_char = upper(char(string(side)));
side_char = side_char(1);
if ~ismember(side_char, {'A', 'B'}), side_char = 'A'; end

n_free = n_node - 1;
x_free = 1:n_free;

switch side_char
    case 'A'
        x_boundary  = 0;
        inacc_nodes = (n_free - n_interface + 1):n_free;
    case 'B'
        x_boundary  = n_node;
        inacc_nodes = 1:n_interface;
end

% Sanitise marker sets to valid free-node indices.
sensor_dofs     = sensor_dofs(:).';
excitation_dofs = excitation_dofs(:).';
sensor_dofs     = sensor_dofs(    sensor_dofs     >= 1 & sensor_dofs     <= n_free);
excitation_dofs = excitation_dofs(excitation_dofs >= 1 & excitation_dofs <= n_free);

% Plain nodes = free nodes that carry no special marker.
special    = unique([inacc_nodes, sensor_dofs, excitation_dofs]);
plain_nodes = setdiff(x_free, special);

mcmap          = cmap_magma(6);
color_sensor   = mcmap(2, :);
color_excit    = mcmap(4, :);
color_inacc    = [0.85 0.15 0.15];
color_node     = [0.45 0.45 0.45];
color_boundary = [0.10 0.10 0.10];

figure('Color', 'w');
ax = axes;
hold(ax, 'on');

% Continuous beam line through the boundary node and every free node.
x_line = sort([x_boundary, x_free]);
plot(ax, x_line, zeros(size(x_line)), '-', 'Color', color_boundary, ...
     'LineWidth', 2.5, 'HandleVisibility', 'off');

handles = gobjects(0);

if ~isempty(plain_nodes)
    h = scatter(ax, plain_nodes, zeros(size(plain_nodes)), 60, color_node, 'o', 'filled', ...
                'DisplayName', 'Beam node');
    handles(end+1) = h;
end

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

% Fixed support: distinct pentagram marker (not the sensor square).
h = scatter(ax, x_boundary, 0, 320, color_boundary, 'pentagram', 'filled', ...
            'DisplayName', 'Fixed boundary');
handles(end+1) = h;

hold(ax, 'off');

xlim(ax, [min(x_line) - 1, max(x_line) + 1]);
ylim(ax, [-0.5, 0.5]);
yticks(ax, []);
ax.YColor = 'none';

xlabel(ax, 'Node index');
title(ax, sprintf('Beam %s  —  sensor & excitation layout', side_char));
grid(ax, 'off');
legend(ax, handles, 'Location', 'best', 'FontSize', 10);

end
