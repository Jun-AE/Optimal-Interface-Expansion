% Define the number of nodes
num_nodes = 21;

% Equally spaced nodes along the x-axis for the first beam
x_nodes_beam1 = linspace(0, 1, num_nodes);

% Vertical positions of the nodes for the first beam (all at y=0)
y_nodes_beam1 = zeros(1, num_nodes);

% Equally spaced nodes along the x-axis for the second beam (shifted horizontally)
x_nodes_beam2 = x_nodes_beam1 + x_nodes_beam1(19); % Adjust the horizontal shift as needed

% Vertical positions of the nodes for the second beam (all at y=-0.2)
y_nodes_beam2 = -0.08 * ones(1, num_nodes);

% Define the new colors
accessible_color = [0.533, 0.847, 0.69]; % RGB: 136,216,176
inaccessible_color = [1, 0.435, 0.412]; % RGB: 255,111,105

% Create darker edges for nodes
darker_edge_accessible = accessible_color * 0.4;  % 70% of original color for darker edge
darker_edge_inaccessible = inaccessible_color * 0.4;  % 70% of original color for darker edge

% Calculate element length (distance between adjacent nodes)
element_length = x_nodes_beam1(2) - x_nodes_beam1(1);
line_length = element_length; % Line length is equal to element length
horizontal_length = line_length/2; % Horizontal lines half the length of vertical line

% Plot the beams - Set transparent figure background
figure('Position', [100, 100, 800, 400], 'Color', 'white');
hold on;

% Set transparent axes background
ax = gca;
ax.Color = 'none';
ax.XColor = 'none';
ax.YColor = 'none';
ax.Box = 'off';

% FIRST: Plot all beam lines to ensure they are behind markers
hBeam = plot([x_nodes_beam1(1), x_nodes_beam1(2)], [y_nodes_beam1(1), y_nodes_beam1(2)], 'k-', 'LineWidth', 2);
for i = 2:num_nodes-1
    plot([x_nodes_beam1(i), x_nodes_beam1(i+1)], [y_nodes_beam1(i), y_nodes_beam1(i+1)], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
end

for i = 1:num_nodes-1
    plot([x_nodes_beam2(i), x_nodes_beam2(i+1)], [y_nodes_beam2(i), y_nodes_beam2(i+1)], 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
end

% Draw vertical dotted lines between beams
for i = 1:3
    plot([x_nodes_beam1(end-i+1), x_nodes_beam1(end-i+1)], [y_nodes_beam1(end-i+1), y_nodes_beam2(i)], 'k:', 'LineWidth', 2, 'HandleVisibility', 'off');
end

% Define which nodes are inaccessible and fixed
inaccessible_beam1 = false(1, num_nodes);
inaccessible_beam1(end-2:end) = true;

inaccessible_beam2 = false(1, num_nodes);
inaccessible_beam2(1:3) = true;

% Define which nodes are fixed (have fixed supports)
fixed_beam1 = false(1, num_nodes);
fixed_beam1(1) = true;  % First node of beam 1 is fixed

fixed_beam2 = false(1, num_nodes);
fixed_beam2(end) = true;  % Last node of beam 2 is fixed

% FIRST: Plot cross markers for inaccessible nodes (so they appear behind circles)
% Plot crosses on inaccessible nodes for beam 1
inaccessible_indices = find(inaccessible_beam1 & ~fixed_beam1);
for i = 1:length(inaccessible_indices)
    plot(x_nodes_beam1(inaccessible_indices(i)), y_nodes_beam1(inaccessible_indices(i)), 'x', 'MarkerSize', 15, 'LineWidth', 2, 'Color', inaccessible_color, 'HandleVisibility', 'off');
end

% Plot crosses on inaccessible nodes for beam 2
inaccessible_indices = find(inaccessible_beam2 & ~fixed_beam2);
for i = 1:length(inaccessible_indices)
    plot(x_nodes_beam2(inaccessible_indices(i)), y_nodes_beam2(inaccessible_indices(i)), 'x', 'MarkerSize', 15, 'LineWidth', 2, 'Color', inaccessible_color, 'HandleVisibility', 'off');
end

% Create dummy objects for the combined legend entry
x_dummy = -100;
y_dummy = -100;
% Create cross first (to appear behind), then circle
hInaccessibleCross = plot(x_dummy, y_dummy, 'x', 'MarkerSize', 24, 'LineWidth', 1.5, 'Color', inaccessible_color);
hold on;
hInaccessibleCircle = plot(x_dummy, y_dummy, 'o', 'MarkerSize', 8, 'MarkerFaceColor', inaccessible_color, 'MarkerEdgeColor', darker_edge_inaccessible);
% Group them together
hInaccessible = [hInaccessibleCircle, hInaccessibleCross];

% Plot accessible DoFs for beam 1 (excluding fixed nodes)
accessible_indices = find(~inaccessible_beam1 & ~fixed_beam1);
if ~isempty(accessible_indices)
    hAccessible = plot(x_nodes_beam1(accessible_indices(1)), y_nodes_beam1(accessible_indices(1)), 'o', 'MarkerSize', 8, 'MarkerFaceColor', accessible_color, 'MarkerEdgeColor', darker_edge_accessible);
    for i = 2:length(accessible_indices)
        plot(x_nodes_beam1(accessible_indices(i)), y_nodes_beam1(accessible_indices(i)), 'o', 'MarkerSize', 8, 'MarkerFaceColor', accessible_color, 'MarkerEdgeColor', darker_edge_accessible, 'HandleVisibility', 'off');
    end
end

% Plot inaccessible DoFs for beam 1 (same color as crosses, excluding fixed nodes)
inaccessible_indices = find(inaccessible_beam1 & ~fixed_beam1);
if ~isempty(inaccessible_indices)
    for i = 1:length(inaccessible_indices)
        plot(x_nodes_beam1(inaccessible_indices(i)), y_nodes_beam1(inaccessible_indices(i)), 'o', 'MarkerSize', 8, 'MarkerFaceColor', inaccessible_color, 'MarkerEdgeColor', darker_edge_inaccessible, 'HandleVisibility', 'off');
    end
end

% Plot accessible DoFs for beam 2 (excluding fixed nodes)
accessible_indices = find(~inaccessible_beam2 & ~fixed_beam2);
for i = 1:length(accessible_indices)
    plot(x_nodes_beam2(accessible_indices(i)), y_nodes_beam2(accessible_indices(i)), 'o', 'MarkerSize', 8, 'MarkerFaceColor', accessible_color, 'MarkerEdgeColor', darker_edge_accessible, 'HandleVisibility', 'off');
end

% Plot inaccessible DoFs for beam 2 (same color as crosses, excluding fixed nodes)
inaccessible_indices = find(inaccessible_beam2 & ~fixed_beam2);
for i = 1:length(inaccessible_indices)
    plot(x_nodes_beam2(inaccessible_indices(i)), y_nodes_beam2(inaccessible_indices(i)), 'o', 'MarkerSize', 8, 'MarkerFaceColor', inaccessible_color, 'MarkerEdgeColor', darker_edge_inaccessible, 'HandleVisibility', 'off');
end

% Parameters for fixed supports
num_bars = 8;  % Number of slanted lines (excluding the main vertical line)
line_width = 1.2;  % Reduced line width
slant_angle = 18;  % Angle in degrees for the slant (from horizontal)

% First fixed node (beam A)
x_support1 = x_nodes_beam1(1);
y_support1 = y_nodes_beam1(1);

% Draw vertical line for beam A
plot([x_support1, x_support1], [y_support1+line_length/2, y_support1-line_length/2], 'k-', 'LineWidth', line_width, 'HandleVisibility', 'off');

% Calculate slant parameters
slant_dx = horizontal_length;
slant_dy = slant_dx * tan(slant_angle * pi/180);

% Draw slanted lines for beam A (extend to the left and slant downward)
y_positions = linspace(y_support1-line_length/2, y_support1+line_length/2, num_bars+2);
% Skip the first and last positions to ensure they don't overlap with the ends
for i = 2:num_bars+1
    % Start point on vertical line
    x_start = x_support1;
    y_start = y_positions(i);
    
    % End point: slanted
    x_end = x_support1 - slant_dx;
    y_end = y_start - slant_dy;
    
    plot([x_start, x_end], [y_start, y_end], 'k-', 'LineWidth', line_width, 'HandleVisibility', 'off');
end

% Second fixed node (beam B)
x_support2 = x_nodes_beam2(end);
y_support2 = y_nodes_beam2(end);

% Draw vertical line for beam B
plot([x_support2, x_support2], [y_support2+line_length/2, y_support2-line_length/2], 'k-', 'LineWidth', line_width, 'HandleVisibility', 'off');

% Draw slanted lines for beam B (extend to the right and slant downward)
y_positions = linspace(y_support2-line_length/2, y_support2+line_length/2, num_bars+2);
for i = 2:num_bars+1
    % Start point on vertical line
    x_start = x_support2;
    y_start = y_positions(i);
    
    % End point: slanted
    x_end = x_support2 + slant_dx;
    y_end = y_start - slant_dy;
    
    plot([x_start, x_end], [y_start, y_end], 'k-', 'LineWidth', line_width, 'HandleVisibility', 'off');
end

% Plot fixed nodes as very tiny circles (practically invisible)
fixed_indices_beam1 = find(fixed_beam1);
fixed_indices_beam2 = find(fixed_beam2);
plot(x_nodes_beam1(fixed_indices_beam1), y_nodes_beam1(fixed_indices_beam1), 'o', 'MarkerSize', 0.2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
plot(x_nodes_beam2(fixed_indices_beam2), y_nodes_beam2(fixed_indices_beam2), 'o', 'MarkerSize', 0.2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');

% Remove axis labels
set(gca, 'XTickLabel', [], 'YTickLabel', []);

% Update the legend to exclude the fixed support and merge inaccessible entries
legend([hAccessible, hInaccessible], { '- Accessible DoFs', '- Inaccessible DoFs'}, 'FontSize', 16, ...
    'Box','off',FontName='FixedWidth');

% Label Beam 1 and Beam 2 as substructures A and B
%text(0.5, 0.05, 'Substructure A (Cantilever beam)', 'HorizontalAlignment', 'center', 'FontSize', 15, 'FontWeight', 'bold');
%text(x_nodes_beam2(end)/2 + x_nodes_beam2(1)/2, y_nodes_beam2(1) - 0.05, 'Substructure B (Cantilever beam)', 'HorizontalAlignment', 'center', 'FontSize', 15, 'FontWeight', 'bold');

% Add DoF numbering for beam A (starting from 1, excluding fixed nodes)
dof_counter = 1;
for i = 1:num_nodes
    if ~fixed_beam1(i)
        text(x_nodes_beam1(i), y_nodes_beam1(i)+0.02, num2str(dof_counter), 'FontSize', 8, ...
             'Color', darker_edge_accessible, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        dof_counter = dof_counter + 1;
    end
end

% Add DoF numbering for beam B (also starting from 1, excluding fixed nodes)
dof_counter = 1;
for i = 1:num_nodes
    if ~fixed_beam2(i)
        text(x_nodes_beam2(i), y_nodes_beam2(i)-0.02, num2str(dof_counter), 'FontSize', 8, ...
             'Color', darker_edge_accessible, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        dof_counter = dof_counter + 1;
    end
end


% Add the title
title('Substructures A and B', 'FontSize', 18, 'FontWeight', 'bold');

% Adjust plot appearance
xlim([-0.1 2])
ylim([-0.4 0.2])

hold off;

% Expand axes to fill entire figure
ax = gca;
ax.Units = 'normalized';
ax.Position = [0 0 1 1]; % [left bottom width height]

% Hide figure margins
set(gcf, 'Units', 'normalized', 'Position', [0 0 0.42 0.3]);
% Corrected exportgraphics syntax
exportgraphics(gcf, 'beam_plot.pdf', ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'none', ...
    'Resolution', 800);



