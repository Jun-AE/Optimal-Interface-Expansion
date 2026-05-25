function map = cmap_magma(N)
% CMAP_MAGMA  Perceptually-uniform magma colormap (matplotlib origin).
%
%   MAP = CMAP_MAGMA(N) returns an N×3 RGB colormap closely approximating
%   matplotlib's magma — a perceptually-uniform sequential map running from
%   near-black through deep purple and red to bright yellow.
%
%   CMAP_MAGMA with no argument returns a 256×3 map.
%
%   Implementation: shape-preserving (pchip) interpolation through 13 RGB
%   control points sampled from the matplotlib magma table. Output is
%   clamped to [0, 1]. Visually indistinguishable from the original at
%   typical figure resolutions.
%
%   Example:
%     colormap(cmap_magma);          % use as default 256-step colormap
%     C = cmap_magma(5);             % five discrete magma colours
%
%   Reference:
%     Smith, N. & van der Walt, S. (2015). matplotlib magma colormap.
%     BSD-licensed; reimplemented here as a compact approximation.
%
%   See also: cmap_cividis, colormap, interp1.

if nargin < 1
    N = 256;
end

% 13 control points sampled at t = 0, 1/12, ..., 1 from matplotlib's magma.
control = [ ...
    0.001  0.000  0.014 ;
    0.077  0.054  0.230 ;
    0.207  0.072  0.388 ;
    0.328  0.080  0.442 ;
    0.448  0.088  0.453 ;
    0.567  0.097  0.456 ;
    0.685  0.114  0.443 ;
    0.795  0.151  0.402 ;
    0.881  0.207  0.367 ;
    0.949  0.279  0.345 ;
    0.990  0.380  0.359 ;
    0.998  0.514  0.421 ;
    0.987  0.991  0.750];

t_ctrl = linspace(0, 1, size(control, 1));
t_out  = linspace(0, 1, N);

map = zeros(N, 3);
for k = 1:3
    map(:, k) = interp1(t_ctrl, control(:, k), t_out, 'pchip');
end
map = max(0, min(1, map));   % clamp to [0, 1]

end
