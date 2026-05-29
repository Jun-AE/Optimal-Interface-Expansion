function y_couple = couple_substructures(y_a, y_b, n_interface_nodes, coupling_mode)
% COUPLE_SUBSTRUCTURES  Dual LM-FBS coupling of two FRF substructures.
%
%   y_couple = couple_substructures(y_a, y_b, n_interface_nodes)
%   y_couple = couple_substructures(y_a, y_b, n_interface_nodes, coupling_mode)
%
%   y_a, y_b          - FRF arrays, each n x n x nFreq
%   n_interface_nodes - number of shared interface nodes
%   coupling_mode     - 'transverse' (1 DoF / node, z only) or
%                       'both'       (2 DoFs / node, z + rotation).
%                       Default 'transverse'.
%   Can be modified for 3 DoFs as well.
%
%   Algorithm (paper §2.2, Eq. 10–11):
%       Y_block = blkdiag(y_a, y_b)
%       Y_AB    = Y_block - Y_block B' (B Y_block B')^(-1) B Y_block
%       y_couple = Y_AB with B-side interface rows/cols removed
%
%   Assumes A's trailing interface DoFs match B's leading interface DoFs.
%

if nargin < 4 || isempty(coupling_mode), coupling_mode = 'transverse'; end

switch lower(coupling_mode)
    case {'transverse', 'z'}
        dof_per_node = 1;
    otherwise
        dof_per_node = 2;
end

n_int   = n_interface_nodes * dof_per_node;
[n_a, ~, n_freq] = size(y_a);
n_b     = size(y_b, 1);
n_total = n_a + n_b;

if dof_per_node == 2
    layoutTag = 'full';
else
    layoutTag = 'translational';
end
validate_interface_dofs(n_interface_nodes, n_a, n_b, ...
    'Layout', layoutTag, 'Caller', 'couple_substructures');

% Boolean compatibility: A's last n_int DoFs == B's first n_int DoFs.
b_mat = zeros(n_int, n_total);
b_mat(:, (n_a - n_int + 1):n_a)     = -eye(n_int);
b_mat(:, n_a + 1 : n_a + n_int)     =  eye(n_int);

% Block-diagonal Y stack.
y_block                              = zeros(n_total, n_total, n_freq);
y_block(1:n_a, 1:n_a, :)             = y_a;
y_block(n_a + 1:end, n_a + 1:end, :) = y_b;

% Dual LM-FBS assembly, batched across frequency.
yb       = pagemtimes(y_block, b_mat.');
byb      = pagemtimes(b_mat,  yb);
by       = pagemtimes(b_mat,  y_block);
y_couple = y_block - pagemtimes(yb, pagemldivide(byb, by));

% Drop B-side interface duplicates; A-side interface kept.
dup_idx                 = n_a + 1 : n_a + n_int;
y_couple(dup_idx, :, :) = [];
y_couple(:, dup_idx, :) = [];

end
