function validate_interface_dofs(n_interface_nodes, n_a, n_b, varargin)
% validate_interface_dofs Assert LM-FBS interface DoF layout contract.
%
%   validate_interface_dofs(N_INTERFACE_NODES, N_A, N_B)
%   validate_interface_dofs(..., 'Layout', 'full' | 'translational', ...)
%
%   Contract (same as primal_coupling / couple_substructures):
%     - Substructure A: interface DoFs are the LAST  n_int rows/cols.
%     - Substructure B: interface DoFs are the FIRST n_int rows/cols.
%
%   'full' layout:          n_int = 2 * N_INTERFACE_NODES (transverse + rotation).
%   'translational' layout: n_int = N_INTERFACE_NODES (transverse only).
%
%   Optional name-value pairs:
%     'Layout'            - 'full' (default) or 'translational'
%     'InterfaceDofsA'    - user-defined interface index vector for beam A
%     'InterfaceDofsB'    - user-defined interface index vector for beam B
%     'Caller'            - string included in error messages (default '')
%

p = inputParser;
addParameter(p, 'Layout', 'full', @(s) any(strcmpi(s, {'full', 'translational'})));
addParameter(p, 'InterfaceDofsA', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'InterfaceDofsB', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'Caller', '', @(s) ischar(s) || isstring(s));
parse(p, varargin{:});

layout = lower(string(p.Results.Layout));
caller = char(p.Results.Caller);
if strlength(caller) > 0
    prefix = sprintf('%s: ', caller);
else
    prefix = '';
end

switch layout
    case "full"
        n_int = 2 * n_interface_nodes;
    case "translational"
        n_int = n_interface_nodes;
end

if n_int < 1 || n_int > n_a || n_int > n_b
    error('validate_interface_dofs:BadInterfaceSize', ...
        '%sinterface band size %d invalid for n_a=%d, n_b=%d.', prefix, n_int, n_a, n_b);
end

expected_a = (n_a - n_int + 1) : n_a;
expected_b = 1 : n_int;

if ~isempty(p.Results.InterfaceDofsA)
    ia = p.Results.InterfaceDofsA(:).';
    if ~isequal(ia, expected_a)
        error('validate_interface_dofs:InterfaceA', ...
            '%sSubstructure A interface DoFs [%s] do not match expected [%s].', ...
            prefix, mat2str(ia), mat2str(expected_a));
    end
end

if ~isempty(p.Results.InterfaceDofsB)
    ib = p.Results.InterfaceDofsB(:).';
    if ~isequal(ib, expected_b)
        error('validate_interface_dofs:InterfaceB', ...
            '%sSubstructure B interface DoFs [%s] do not match expected [%s].', ...
            prefix, mat2str(ib), mat2str(expected_b));
    end
end

end
