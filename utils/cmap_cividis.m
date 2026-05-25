function map = cmap_cividis(N)
% CMAP_CIVIDIS  Perceptually-uniform, CVD-safe cividis colormap.
%
%   MAP = CMAP_CIVIDIS(N) returns an N×3 RGB colormap closely approximating
%   the cividis colormap — a perceptually-uniform sequential map designed
%   to remain readable under colour-vision deficiency, running from dark
%   blue through neutral grey to bright yellow.
%
%   CMAP_CIVIDIS with no argument returns a 256×3 map.
%
%   Implementation: shape-preserving (pchip) interpolation through 9 RGB
%   control points sampled from the cividis table. Output is clamped to
%   [0, 1]. Visually indistinguishable from the original at typical
%   figure resolutions.
%
%   Example:
%     colormap(cmap_cividis);        % use as default 256-step colormap
%     C = cmap_cividis(7);           % seven discrete cividis colours
%
%   Reference:
%     Nuñez, Anderton & Renslow (2018). Optimizing colormaps with consideration
%     for color vision deficiency. PLoS ONE 13(7): e0199239.
%     Public-domain RGB table; reimplemented here as a compact approximation.
%
%   See also: cmap_magma, colormap, interp1.

if nargin < 1
    N = 256;
end

% 9 control points sampled at t = 0, 0.125, ..., 1 from the cividis table.
control = [ ...
    0.000  0.135  0.305 ;
    0.000  0.196  0.430 ;
    0.205  0.300  0.450 ;
    0.330  0.385  0.450 ;
    0.450  0.480  0.450 ;
    0.570  0.560  0.443 ;
    0.730  0.685  0.385 ;
    0.890  0.815  0.255 ;
    1.000  0.948  0.060];

t_ctrl = linspace(0, 1, size(control, 1));
t_out  = linspace(0, 1, N);

map = zeros(N, 3);
for k = 1:3
    map(:, k) = interp1(t_ctrl, control(:, k), t_out, 'pchip');
end
map = max(0, min(1, map));   % clamp to [0, 1]

end
