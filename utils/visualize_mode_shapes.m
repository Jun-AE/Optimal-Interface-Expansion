function visualize_mode_shapes(NodeCoords, modeShapes, natFreq, varargin)
%  Visualise (and optionally animate) mode shapes on a 1D or 2D-regular grid.
%
%   VISUALIZE_MODE_SHAPES(..., NAME, VALUE, ...) accepts:
%     'Animate'      - (logical) Animate the sinusoidal motion. Default false.
%     'SaveVideo'    - (logical) Write the animation to disk. Implies 'Animate' true.
%                      Default false.
%     'VideoFile'    - (char)    Output file name. Default 'modeShapeAnimation.mp4'.
%     'NumFrames'    - (integer) Frames per single animation pass. Default 600.
%     'FrameRate'    - (double)  Playback / video frame rate [fps]. Default 30.
%     'Loop'         - (logical) When true, the live animation repeats
%                      continuously until the figure is closed. Default
%                      true when SaveVideo is false, false otherwise.
%                      (Single-pass live playback would otherwise finish
%                      in a few seconds and freeze on the last phase.)
%     'Exaggeration' - (double)  Amplitude scaling factor. Default 1.
%
%   Inputs:
%     NODECOORDS - (n×≥2 double) Spatial coordinates of each node. Columns
%                  1 and 2 are interpreted as (x, y). 1D layouts (one
%                  coordinate constant) and 2D rectangular grids are both
%                  supported.
%     MODESHAPES - (n×nModes double) Each column is a mode shape vector.
%     NATFREQ    - (nModes×1 double) [Hz] Natural frequency for each mode.
%
%   Coordinate layout detection:
%     - If all y-coordinates are equal (or all x), the layout is treated as 1D.
%     - If numel(unique x) * numel(unique y) == n, the layout is a regular
%       2D rectangular grid.
%
%   Example:
%     % Static figure
%     visualize_mode_shapes(NodeCoords, Phi, fn);
%
%     % Animated + saved to MP4
%     visualize_mode_shapes(NodeCoords, Phi, fn, ...
%         'Animate', true, 'SaveVideo', true, 'VideoFile', 'beam.mp4');
%
%   See also: cmap_cividis, calculate_mac, plot_mode_shape.

%% Parse options

opts = parseOptions(varargin{:});

%% Sort node coordinates and reorder mode shapes accordingly

[NodeCoords_sorted, sortIdx] = sortrows(NodeCoords, [2, 1]);
modeShapes_sorted            = modeShapes(sortIdx, :);
numModes                     = size(modeShapes_sorted, 2);

%% Detect grid layout: 1D, 2D-regular, or unsupported

tol      = 1e-4;
x_coords = uniquetol(NodeCoords_sorted(:, 1), tol, 'DataScale', 1);
y_coords = uniquetol(NodeCoords_sorted(:, 2), tol, 'DataScale', 1);

layout = detectLayout(NodeCoords_sorted, x_coords, y_coords);

%% Draw  dispatch on layout

[fig, plotHandles, modeData] = setUpFigure(numModes, opts);

switch layout.kind
    case '1D'
        plotHandles = draw1DModes(fig, plotHandles, modeData, ...
                                  layout, modeShapes_sorted, natFreq);
    case '2D'
        plotHandles = draw2DRegularModes(fig, plotHandles, modeData, ...
                                          layout, modeShapes_sorted, natFreq);
end

%% Optional animation / video

if opts.Animate || opts.SaveVideo
    runAnimation(fig, plotHandles, modeData, natFreq, layout, opts);
end

end

% =========================================================================
% Option parsing
% =========================================================================
function opts = parseOptions(varargin)
p = inputParser;
addParameter(p, 'Animate',      false,  @islogical);
addParameter(p, 'SaveVideo',    false,  @islogical);
addParameter(p, 'VideoFile',    'modeShapeAnimation.mp4', @(s) ischar(s) || isstring(s));
addParameter(p, 'NumFrames',    600,    @(x) isnumeric(x) && x > 0);
addParameter(p, 'FrameRate',    30,     @(x) isnumeric(x) && x > 0);
addParameter(p, 'Loop',         [],     @(x) isempty(x) || islogical(x));   % default below
addParameter(p, 'Exaggeration', 1.0,    @isnumeric);
parse(p, varargin{:});
opts = p.Results;
if opts.SaveVideo, opts.Animate = true; end             % video implies animation
if isempty(opts.Loop), opts.Loop = ~opts.SaveVideo; end % loop live, single-pass video
end

% =========================================================================
% Layout detection  1D, 2D-regular, or unsupported
% =========================================================================
function layout = detectLayout(NodeCoords_sorted, x_coords, y_coords)

if isscalar(y_coords)
    layout.kind        = '1D';
    layout.varyingAxis = 'x';
    layout.coords      = x_coords;
    layout.fixedVal    = y_coords(1);
elseif isscalar(x_coords)
    layout.kind        = '1D';
    layout.varyingAxis = 'y';
    layout.coords      = y_coords;
    layout.fixedVal    = x_coords(1);
elseif numel(x_coords) * numel(y_coords) == size(NodeCoords_sorted, 1)
    layout.kind     = '2D';
    layout.x_coords = x_coords;
    layout.y_coords = y_coords;
else
    error('visualize_mode_shapes:UnsupportedLayout', ...
        ['Node layout is neither 1D nor a regular 2D rectangular grid. ', ...
         'Irregular-grid support was removed in this version.']);
end
end

% =========================================================================
% Figure setup (shared)
% =========================================================================
function [fig, h, modeData] = setUpFigure(numModes, opts)

fig = figure('Name', 'Mode Shape Visualisation', 'Color', 'white');
if opts.SaveVideo
    screen = get(groot, 'ScreenSize');                 % [l b w h]
    figW   = min(1280, screen(3) - 100);
    figH   = min( 720, screen(4) - 150);
    fig.Position = [50, 50, figW, figH];
    fig.Resize   = 'off';                              % lock against window-manager nudges
end
colormap(fig, cmap_cividis(256));

h.tlo  = tiledlayout(fig, 1, numModes, ...
                     'Padding', 'compact', 'TileSpacing', 'compact');
h.ax   = gobjects(numModes, 1);
h.surf = gobjects(numModes, 1);
h.line = gobjects(numModes, 1);
h.dots = gobjects(numModes, 1);

modeData.smooth = cell(numModes, 1);
modeData.nodes  = cell(numModes, 1);
end

% =========================================================================
% 1D drawer
% =========================================================================
function h = draw1DModes(~, h, ~, layout, modeShapes_sorted, natFreq)

numModes            = size(modeShapes_sorted, 2);
interpolationFactor = 4;
coords_fine         = linspace(min(layout.coords), max(layout.coords), ...
                               numel(layout.coords) * interpolationFactor);

for m = 1:numModes
    h.ax(m) = nexttile(h.tlo);

    modeVec      = real(modeShapes_sorted(:, m));
    modeVec      = normaliseMode(modeVec);
    modeSmooth   = interp1(layout.coords, modeVec, coords_fine, 'spline');

    [Xline, Yline, Xnodes, Ynodes] = layoutCoords1D(layout, coords_fine);

    h.line(m) = surface([Xline; Xline], [Yline; Yline], ...
                        [modeSmooth; modeSmooth], [modeSmooth; modeSmooth], ...
                        'EdgeColor', 'interp', 'FaceColor', 'none', 'LineWidth', 6);
    hold(h.ax(m), 'on');
    h.dots(m) = scatter3(Xnodes, Ynodes, zeros(size(layout.coords)), 60, 'k', 'filled');

    titleAndAxes1D(h.ax(m), m, natFreq(m));
    hold(h.ax(m), 'off');

    h.smooth{m} = modeSmooth;
    h.modes{m}  = modeVec;
end

addSharedColorbar(h);
end

function [Xline, Yline, Xnodes, Ynodes] = layoutCoords1D(layout, coords_fine)
if layout.varyingAxis == 'x'
    Xline  = coords_fine;
    Yline  = layout.fixedVal * ones(size(coords_fine));
    Xnodes = layout.coords;
    Ynodes = layout.fixedVal * ones(size(layout.coords));
else
    Xline  = layout.fixedVal * ones(size(coords_fine));
    Yline  = coords_fine;
    Xnodes = layout.fixedVal * ones(size(layout.coords));
    Ynodes = layout.coords;
end
end

function titleAndAxes1D(ax, m, fn)
title(ax, sprintf('Mode %d: %.2f Hz', m, fn), 'FontSize', 14);
xlabel(ax, 'X'); ylabel(ax, 'Y'); zlabel(ax, 'Normalised displacement');
set(ax, 'ZDir', 'reverse', 'Color', 'white');
view(ax, 0, 180);
axis(ax, 'tight'); zlim(ax, [-1 1]);
grid(ax, 'on'); ax.GridAlpha = 0.25;
end

% =========================================================================
% 2D regular grid drawer
% =========================================================================
function h = draw2DRegularModes(~, h, ~, layout, modeShapes_sorted, natFreq)

numModes = size(modeShapes_sorted, 2);
nx       = numel(layout.x_coords);
ny       = numel(layout.y_coords);

x_fine = linspace(min(layout.x_coords), max(layout.x_coords), nx * 4);
y_fine = linspace(min(layout.y_coords), max(layout.y_coords), ny * 4);
[X_fine, Y_fine] = meshgrid(x_fine, y_fine);
[X, Y]           = meshgrid(layout.x_coords, layout.y_coords);

for m = 1:numModes
    h.ax(m) = nexttile(h.tlo);

    modeVec       = real(modeShapes_sorted(:, m));
    modeGrid      = reshape(modeVec, [ny, nx]);
    modeGrid      = modeGrid / max(abs(modeGrid(:)));
    modeSmooth    = interp2(X, Y, modeGrid, X_fine, Y_fine, 'spline');

    h.surf(m) = surf(X_fine, Y_fine, modeSmooth, 'EdgeColor', 'none');
    shading(h.ax(m), 'interp');
    hold(h.ax(m), 'on');
    mesh(X, Y, zeros(size(X)), 'EdgeColor', [0.5 0.5 0.5], 'FaceAlpha', 0.1);

    titleAndAxes2D(h.ax(m), m, natFreq(m));
    hold(h.ax(m), 'off');

    h.smooth{m} = modeSmooth;
end

addSharedColorbar(h);
end

function titleAndAxes2D(ax, m, fn)
title(ax, sprintf('Mode %d: %.2f Hz', m, fn), 'FontSize', 14);
xlabel(ax, 'X'); ylabel(ax, 'Y'); zlabel(ax, 'Normalised displacement');
set(ax, 'Color', 'white');
view(ax, -30, 30);
axis(ax, 'equal', 'tight');
zlim(ax, [-1 1]);
grid(ax, 'on'); ax.GridAlpha = 0.25;
end

% =========================================================================
% Shared cosmetic helpers
% =========================================================================
function addSharedColorbar(h)
for i = 1:numel(h.ax)
    if isgraphics(h.ax(i))
        clim(h.ax(i), [-1, 1]);
    end
end
cb              = colorbar(h.ax(end));
cb.Layout.Tile  = 'east';
cb.Label.String = 'Normalised displacement';
end

function v = normaliseMode(v)
m = max(abs(v));
if m > 1e-12, v = v / m; end
end

% =========================================================================
% Animation loop (+ optional video writing)
% =========================================================================
function runAnimation(fig, h, ~, natFreq, layout, opts)

timeScale = 1 / max(natFreq);
T_total   = opts.NumFrames / opts.FrameRate;
t_anim    = linspace(0, T_total, opts.NumFrames);
dt        = 1 / opts.FrameRate;

if opts.SaveVideo
    v           = VideoWriter(opts.VideoFile, 'MPEG-4');
    v.FrameRate = opts.FrameRate;
    v.Quality   = 100;
    open(v);
end

drawnow;  pause(0.1);

targetSize = [];                       % set when the first frame is captured
numModes   = numel(h.ax);

keepLooping = true;
while keepLooping
    for k = 1:opts.NumFrames
        for m = 1:numModes
            phase = 2 * pi * natFreq(m) * timeScale * t_anim(k);
            switch layout.kind
                case '1D'
                    Zline = opts.Exaggeration * h.smooth{m} * sin(phase);
                    set(h.line(m), 'ZData', [Zline; Zline], 'CData', [Zline; Zline]);
                    set(h.dots(m), 'ZData', h.modes{m} * sin(phase) * opts.Exaggeration);
                case '2D'
                    Zanim = opts.Exaggeration * h.smooth{m} * sin(phase);
                    set(h.surf(m), 'ZData', real(Zanim));
            end
        end

        drawnow;                                       % full repaint, no frame-skip
        if ~opts.SaveVideo
            pause(dt);                                 % live pacing → smooth real-time
        end

        if opts.SaveVideo
            frame = getframe(fig);
            if isempty(targetSize)
                targetSize = size(frame.cdata);
            else
                frame.cdata = forceFrameSize(frame.cdata, targetSize);
            end
            writeVideo(v, frame);
        end

        if ~isvalid(fig)                               % user closed the window
            keepLooping = false;
            break
        end
    end

    if ~opts.Loop, keepLooping = false; end
end

if opts.SaveVideo, close(v); end
end

% =========================================================================
% Defensive frame-size enforcement
% =========================================================================
function img = forceFrameSize(img, target)
[hImg, wImg, c] = size(img);
hT = target(1);  wT = target(2);

% Crop if too tall / wide
if hImg > hT, img = img(1:hT, :, :);     hImg = hT; end
if wImg > wT, img = img(:, 1:wT, :);     wImg = wT; end

% White-pad if too short / narrow
if hImg < hT
    pad = 255 * ones(hT - hImg, wImg, c, 'like', img);
    img = [img; pad];
end
if wImg < wT
    pad = 255 * ones(hT, wT - wImg, c, 'like', img);
    img = [img, pad];
end
end
