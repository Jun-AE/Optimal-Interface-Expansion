function model = load_hcb_model(paths, opts)
% LOAD_HCB_MODEL  Load and assemble an Ansys Craig-Bampton (HCB) reduced model.
%
%   model = load_hcb_model(paths, opts)
%
%   Self-contained reader for a Harwell-Boeing (HB) reduced K/M pair exported
%   from Ansys Mechanical (HCB / CMS super-element), plus its master-node
%   coordinate and DOF-mapping files. It ONLY loads, reads, and processes the
%   Ansys data: it assembles the sparse HB matrices to full, reorders them into
%   ascending master-node order with the CMS generalised modes last, and
%   reports the physical / measured DOF index sets. All dynamics (damping, FRF
%   generation, noise, frequency axes) are left to the caller, which uses the
%   existing utils (damping, compute_frf, add_noise, frequency_generation).
%
%   Inputs:
%     paths - struct of file addresses:
%       .nodes    node coordinate file: rows = "ID c1 c2 c3 ...".
%       .mapping  Ansys ".mapping" file (Matrix-Eqn / Node / DOF table).
%       .K, .M    Harwell-Boeing reduced stiffness / mass files (RSA format).
%     opts  - struct of options (all optional):
%       .dof_per_node       DOFs per master node.                (default 3)
%       .measured_dof_index which DOF within a node block is the measured /
%                           transverse one (1..dof_per_node).     (default 2)
%       .node_coord_cols    column reorder for NodeCoords (plot / aggregation
%                           only; not used in assembly).     (default [2 4 3 1])
%       .verbose            print progress.                       (default true)
%
%   Output struct `model`:
%     .K .M               reordered full reduced matrices (ndof x ndof).
%     .NodeID             sorted master node IDs (nnode x 1).
%     .NodeCoords         node coordinate table after .node_coord_cols reorder.
%     .NCMSMode           number of CMS generalised modes (ndof - dof_per_node*nnode).
%     .physical_DoF_array indices of the physical DOFs (1 : dof_per_node*nnode).
%     .measured_dof_array measured-DOF indices within the physical block.
%     .dof_per_node       echoed for the caller.
%     .n                  number of measured DOFs.
%
%   Example (caller builds the FRFs with the existing utils):
%     paths.nodes='Nodes_36.txt'; paths.mapping='KredHB.mapping';
%     paths.K='KredHB.txt'; paths.M='MredHB.txt';
%     m  = load_hcb_model(paths);
%     [w, w_hz] = frequency_generation(15, 40, 0.1);
%     C  = damping(m.K, m.M, struct('type','modal','ratios',0.0012*ones(1,12)));
%     Yn = compute_frf(w, m.K,            m.M, zeros(size(m.K)), 'accelerance');
%     Ye = compute_frf(w, 1.04 * m.K,     m.M, C,               'accelerance');
%     Yn = Yn(m.physical_DoF_array, m.physical_DoF_array, :);
%     Ye = Ye(m.physical_DoF_array, m.physical_DoF_array, :);
%     Ye = add_noise(true, Ye, 1e-3, 1e-3, 1e-3, 1e-3, 10);
%     Yz  = Yn(m.measured_dof_array, m.measured_dof_array, :);
%     YEz = Ye(m.measured_dof_array, m.measured_dof_array, :);
%
%   Reference:
%     Junaid et al. (2026). Journal of Sound and Vibration. DOI: 10.1016/j.jsv.2026.119782
%     HB reader after the assembled-symmetric (RSA) Harwell-Boeing convention.
%
%   See also: exhaustive_search, compute_frf, damping, add_noise, frequency_generation.

if nargin < 2, opts = struct(); end

dof_per_node       = get_opt(opts, 'dof_per_node',       3);
measured_dof_index = get_opt(opts, 'measured_dof_index', 2);
node_coord_cols    = get_opt(opts, 'node_coord_cols',    [2, 4, 3, 1]);
verbose            = get_opt(opts, 'verbose',            true);

for fld = {'nodes', 'mapping', 'K', 'M'}
    if ~isfield(paths, fld{1}) || isempty(paths.(fld{1}))
        error('load_hcb_model:MissingPath', 'paths.%s is required.', fld{1});
    end
end

%% Node coordinates + master-node count

coords_raw       = load(paths.nodes);
n_master         = size(coords_raw, 1);
model.NodeCoords = coords_raw(:, node_coord_cols);

%% DOF mapping -> sorted master node IDs

[NodeID, node_order] = read_mapping(paths.mapping, n_master, dof_per_node);
model.NodeID = NodeID;

%% Read + assemble K, M (Harwell-Boeing RSA), then reorder

if verbose, fprintf('load_hcb_model: reading HB matrices...\n'); end
K = read_hb_matrix(paths.K);
M = read_hb_matrix(paths.M);

ndof           = size(K, 1);
model.NCMSMode = ndof - dof_per_node * numel(NodeID);

id_dof = build_reorder(NodeID, node_order, dof_per_node, ndof, model.NCMSMode);
model.K = K(id_dof, id_dof);
model.M = M(id_dof, id_dof);

%% Index sets

model.dof_per_node       = dof_per_node;
model.physical_DoF_array = 1:(ndof - model.NCMSMode);
model.measured_dof_array = measured_dof_index:dof_per_node:numel(model.physical_DoF_array);
model.n                  = numel(model.measured_dof_array);

if verbose
    fprintf('  ndof=%d (physical=%d, CMS modes=%d), measured DOFs=%d\n', ...
            ndof, numel(model.physical_DoF_array), model.NCMSMode, model.n);
end

end

% =========================================================================
% Local helpers
% =========================================================================
function v = get_opt(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end


function [NodeID, node_order] = read_mapping(filename, n_master, dof_per_node)
% Parse the Ansys ".mapping" table: header line, then dof_per_node lines per
% master node. The node ID sits in columns 15:28 of each line.
fid = fopen(filename, 'rt');
if fid < 0, error('load_hcb_model:OpenFailed', 'Cannot open %s', filename); end
closer = onCleanup(@() fclose(fid));                                     %#ok<NASGU>

getline(fid);                                   % header row
node_order = zeros(n_master, 1);
for k = 1:n_master
    nid = NaN;
    for d = 1:dof_per_node                      %#ok<NASGU>
        ln  = padline(getline(fid), 28);
        nid = str2double(ln(15:28));
    end
    node_order(k) = nid;
end
NodeID = unique(node_order);                     % ascending master node IDs
end


function id_dof = build_reorder(NodeID, node_order, dof_per_node, ndof, n_cms)
% Permutation that places physical DOFs first, ordered by ascending NodeID,
% with the CMS generalised modes appended last (matches reorder_matrices.m).
nnode  = numel(NodeID);
id_dof = zeros(1, dof_per_node * nnode);
blk    = -(dof_per_node - 1):0;
for k = 1:nnode
    idx = find(NodeID(k) == node_order, 1);     % position in the matrix ordering
    id_dof(dof_per_node * k + blk) = dof_per_node * idx + blk;
end
id_dof = [id_dof, (ndof - n_cms + 1):ndof];
end


function A = read_hb_matrix(filename)
% Compact reader for an assembled real-symmetric (RSA) Harwell-Boeing matrix.
% Reads colptr / rowind / values by their header line counts (robust to the
% number of entries per line) and mirrors the stored triangle to full.
fid = fopen(filename, 'rt');
if fid < 0, error('load_hcb_model:OpenFailed', 'Cannot open %s', filename); end
closer = onCleanup(@() fclose(fid));                                     %#ok<NASGU>

getline(fid);                                   % line 1: title / key
l2 = padline(getline(fid), 70);                 % line 2: card counts
ptrcrd = str2double(l2(15:28));
indcrd = str2double(l2(29:42));
valcrd = str2double(l2(43:56));
rhscrd = str2double(l2(57:70));
l3 = padline(getline(fid), 56);                 % line 3: type + sizes
mxtype = strtrim(l3(1:3));
nrow   = str2double(l3(15:28));
ncol   = str2double(l3(29:42));
getline(fid);                                   % line 4: Fortran formats (unused)
if rhscrd > 0
    getline(fid);                               % line 5: RHS type (present, unused)
end

colptr = read_block(fid, ptrcrd, false);        % ncol+1 column pointers
rowind = read_block(fid, indcrd, false);        % nnzero row indices
values = read_block(fid, valcrd, true);         % nnzero values (D-exponent)

% Expand compressed-column pointers to explicit column indices.
colind = zeros(numel(rowind), 1);
for j = 1:ncol
    colind(colptr(j):colptr(j+1) - 1) = j;
end
A = sparse(rowind, colind, values, nrow, ncol);

% RSA: only one triangle (incl. diagonal) is stored — mirror to full.
if numel(mxtype) >= 2 && upper(mxtype(2)) == 'S'
    A = A + A.' - (A .* speye(nrow, ncol));
end
A = full(A);
end


function v = read_block(fid, n_lines, is_real)
% Read n_lines records and return all numeric tokens. Real values may use a
% Fortran 'D' exponent, which is converted to 'E' before parsing.
parts = cell(n_lines, 1);
for i = 1:n_lines
    ln = getline(fid);
    if is_real
        ln = strrep(strrep(ln, 'D', 'E'), 'd', 'e');
        parts{i} = sscanf(ln, '%f');
    else
        parts{i} = sscanf(ln, '%d');
    end
end
v = vertcat(parts{:});
end


function s = getline(fid)
s = fgetl(fid);
if ~ischar(s)
    error('load_hcb_model:UnexpectedEOF', 'Unexpected end of file while reading HB/mapping data.');
end
end


function s = padline(s, width)
if numel(s) < width
    s(end + 1:width) = ' ';
end
end
