function [beam_num, beam_exp] = create_cantilever_beam(props, gamma, fixed_side)
% CREATE_CANTILEVER_BEAM  FE model of a cantilever beam (numerical + experimental).
%
%   [beam_num, beam_exp] = create_cantilever_beam(props, gamma, fixed_side)
%
%   props.Length, props.nelm, props.E, props.b, props.h, props.rho.
%   gamma       - stiffness scaling for the experimental beam.
%   fixed_side  - 'left' or 'right'.
%
%   Each returned struct carries: nnode, ndof, l, A, I, k, m, K, M, modes, freq.
%
%   Reference:
%     Junaid et al. (2026), Journal of Sound and Vibration.

nelm  = round(props.nelm);
nnode = nelm + 1;
ndof  = 2 * nnode;
l     = props.Length / nelm;
a     = props.b * props.h;
i_sec = (props.b * props.h^3) / 12;

k_elem = ((props.E * i_sec) / l^3) * ...
    [ 12      6*l    -12      6*l   ;
      6*l    4*l^2   -6*l    2*l^2  ;
     -12    -6*l     12     -6*l    ;
      6*l    2*l^2   -6*l    4*l^2  ];

m_elem = ((props.rho * a * l) / 420) * ...
    [ 156     22*l    54    -13*l   ;
      22*l   4*l^2   13*l   -3*l^2  ;
      54     13*l   156    -22*l    ;
     -13*l  -3*l^2  -22*l   4*l^2   ];

k_glob = zeros(ndof);
m_glob = zeros(ndof);
for ne = 1:nelm
    idx = (2*(ne-1) + 1) : (2*(ne-1) + 4);
    k_glob(idx, idx) = k_glob(idx, idx) + k_elem;
    m_glob(idx, idx) = m_glob(idx, idx) + m_elem;
end

beam_num = struct( ...
    'nnode', nnode, 'ndof', ndof, ...
    'l',     l,     'A',    a,    'I', i_sec, ...
    'k',     k_elem,'m',    m_elem, ...
    'K',     k_glob,'M',    m_glob);

beam_exp     = beam_num;
beam_exp.K   = gamma * beam_exp.K;

switch lower(fixed_side)
    case 'left'
        fixed_dofs = [1, 2];
    case 'right'
        fixed_dofs = [ndof - 1, ndof];
end

beam_num.K(fixed_dofs, :) = [];  beam_num.K(:, fixed_dofs) = [];
beam_num.M(fixed_dofs, :) = [];  beam_num.M(:, fixed_dofs) = [];
beam_exp.K(fixed_dofs, :) = [];  beam_exp.K(:, fixed_dofs) = [];
beam_exp.M(fixed_dofs, :) = [];  beam_exp.M(:, fixed_dofs) = [];

beam_num.ndof = size(beam_num.K, 1);
beam_exp.ndof = size(beam_exp.K, 1);

[beam_num.modes, omega_sq] = eig(beam_num.K, beam_num.M);
beam_num.freq = sqrt(diag(omega_sq)) / (2 * pi);

[beam_exp.modes, omega_sq] = eig(beam_exp.K, beam_exp.M);
beam_exp.freq = sqrt(diag(omega_sq)) / (2 * pi);

end
