function map = cmap_cividis(N)
% CMAP_CIVIDIS  Perceptually-uniform, CVD-safe cividis colormap.
%
%   cmap_cividis(N) returns an N×3 RGB colormap
%
%   CMAP_CIVIDIS with no argument returns a 256×3 map.
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
