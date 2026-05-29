% MAIN  Optimal interface expansion — worked examples.
%
%   Two self-contained, clearly separated examples for the methodology in:
%     Junaid et al. (2026). Journal of Sound and Vibration. DOI: 10.1016/j.jsv.2026.119782
%
%   EXAMPLE 1 — Beam (substructure coupling):
%     Two steel cantilevers (A fixed-left, B fixed-right) with Rayleigh damping.
%     Optimal sensor/excitation DoFs are found per beam by an NDP exhaustive
%     search (2 sensors + 3 excitations = 5 measured DoFs each); the resulting
%     SEMM expansions are coupled at the interface and compared against the full
%     coupled model. No plotting.
%
%   EXAMPLE 2 — Square plate (no coupling):
%     An Ansys Craig-Bampton reduced model is loaded, its FRFs built with the
%     same Rayleigh damping, and the optimal sensors/excitations are found by an
%     NDP exhaustive search over all valid DoFs. The placement is plotted.
%
%   Dependencies (utils/): create_cantilever_beam, damping, compute_frf,
%     add_noise, frequency_generation, exhaustive_search, semm, func_coh,
%     primal_coupling, couple_substructures, svd_truncation, load_hcb_model,
%     plot_square_plate.

clear; clc; close all;

if isempty(which('create_cantilever_beam'))
    addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));
end

%% Shared search / damping configuration

ray.type   = 'proportional';   % Rayleigh damping
ray.alpha  = 1.0;              % mass-proportional coefficient
ray.beta   = 1e-5;            % stiffness-proportional coefficient (≈0.4–2% over the bands used)

num_sensors       = 1;        % NDP: 1 sensor ...
extra_excitations = 0;        % ... + (1+0) = 1 excitation → 2 measured DoFs (small NDP, fast verification)


%% ========================================================================
%  EXAMPLE 1 — BEAM: optimal sensors + interface coupling
%  ========================================================================

%% Beam properties (steel)
props.Length = 3;       % [m]
props.nelm   = 20;      % elements
props.E      = 2.1e11;  % [Pa]      steel Young's modulus
props.b      = 0.05;    % [m]       breadth
props.h      = 0.2;     % [m]       height
props.rho    = 7850;    % [kg/m^3]  steel density

gamma_a = 1.05;         % numerical -> experimental stiffness scaling, beam A
gamma_b = 1.05;         % beam B

n_interface = 3;        % shared interface nodes
fb_start = 1;  fb_end = 700;  fb_step = 1;     % FRF band [Hz]
beam_noise_seed = 12;

%% Build cantilever pairs (numerical + experimental) and apply Rayleigh damping
[beam_a, beam_a_exp] = create_cantilever_beam(props, gamma_a, 'left');
[beam_b, beam_b_exp] = create_cantilever_beam(props, gamma_b, 'right');

beam_a.C     = damping(beam_a.K,     beam_a.M,     ray);
beam_b.C     = damping(beam_b.K,     beam_b.M,     ray);
beam_a_exp.C = damping(beam_a_exp.K, beam_a_exp.M, ray);
beam_b_exp.C = damping(beam_b_exp.K, beam_b_exp.M, ray);

%% FRFs (accelerance), reduced to translational DoFs
[wb, wb_hz] = frequency_generation(fb_start, fb_end, fb_step);
YA  = compute_frf(wb, beam_a.K,     beam_a.M,     beam_a.C,     'accelerance');
YB  = compute_frf(wb, beam_b.K,     beam_b.M,     beam_b.C,     'accelerance');
YEA = compute_frf(wb, beam_a_exp.K, beam_a_exp.M, beam_a_exp.C, 'accelerance');
YEB = compute_frf(wb, beam_b_exp.K, beam_b_exp.M, beam_b_exp.C, 'accelerance');

tdof = 1:2:beam_a.ndof;                          % transverse DoFs
YA  = YA(tdof, tdof, :);   YB  = YB(tdof, tdof, :);
YEA = YEA(tdof, tdof, :);  YEB = YEB(tdof, tdof, :);

%% Additive measurement noise on the experimental FRFs
YEA = add_noise(true, YEA, 0.005, 0.005, 1e-8, 1e-8, beam_noise_seed);
YEB = add_noise(true, YEB, 0.005, 0.005, 1e-8, 1e-8, beam_noise_seed);

%% Optimal sensor/excitation search (NDP, 2 sensors + 3 excitations)
% Interface = inaccessible DoFs scored on; candidates exclude the interface.
n = size(YA, 1);
interface_a = (n - n_interface + 1):n;          % beam A interface: free (right) end
interface_b = 1:n_interface;                    % beam B interface: free (left) end

optA = struct('num_sensors', num_sensors, 'search_type', 'NDP', ...
    'extra_excitations', extra_excitations, ...
    'validation_dofs', interface_a, 'candidate_dofs', setdiff(1:n, interface_a), ...
    'methods', 'COH', 'extrema', 'max', 'frequency_range', wb, 'verbose', true);
optB = optA;
optB.validation_dofs = interface_b;
optB.candidate_dofs  = setdiff(1:n, interface_b);

res_a = exhaustive_search(YA, YEA, optA);
res_b = exhaustive_search(YB, YEB, optB);

ys_a = res_a.ys;                                 % SEMM-expanded experimental FRF, beam A
ys_b = res_b.ys;

fprintf('\nBeam A: sensors %s | excitations %s | GCCM %.4f  (%d measured DoFs)\n', ...
    mat2str(res_a.sensor_dofs), mat2str(res_a.excitation_dofs), res_a.cor_overall, ...
    numel(res_a.sensor_dofs) + numel(res_a.excitation_dofs));
fprintf('Beam B: sensors %s | excitations %s | GCCM %.4f  (%d measured DoFs)\n', ...
    mat2str(res_b.sensor_dofs), mat2str(res_b.excitation_dofs), res_b.cor_overall, ...
    numel(res_b.sensor_dofs) + numel(res_b.excitation_dofs));

%% Full coupled (reference) model — primal coupling of the experimental beams
[~, k_full, m_full, c_full] = primal_coupling(n_interface, ...
    beam_a_exp.K, beam_b_exp.K, beam_a_exp.M, beam_b_exp.M, beam_a_exp.C, beam_b_exp.C);
y_full = compute_frf(wb, k_full, m_full, c_full, 'accelerance');
y_full = y_full(1:2:length(k_full), 1:2:length(k_full), :);

%% Optimally coupled SEMM model (transverse interface)
ys_a_t   = svd_truncation(ys_a, 15);             % rank-reduce noisy SEMM before coupling
ys_b_t   = svd_truncation(ys_b, 15);
y_couple = couple_substructures(ys_a_t, ys_b_t, n_interface, 'transverse');

%% Compare coupled SEMM vs full model
nf = min(size(y_couple, 3), size(y_full, 3));
[~, ~, gamma_couple] = func_coh(y_couple(:, :, 1:nf), y_full(:, :, 1:nf));
fprintf('\n[BEAM] coupled-SEMM vs full-model GCCM (Gamma) = %.6f\n', gamma_couple);


%% ========================================================================
%  EXAMPLE 2 — SQUARE PLATE: optimal sensors over all DoFs (no coupling)
%  ========================================================================

data_dir = fullfile(fileparts(mfilename('fullpath')), 'Data');

%% Load Ansys Craig-Bampton reduced model
paths = struct( ...
    'nodes',   fullfile(data_dir, 'Nodes_36.txt'), ...
    'mapping', fullfile(data_dir, 'KredHB.mapping'), ...
    'K',       fullfile(data_dir, 'KredHB.txt'), ...
    'M',       fullfile(data_dir, 'MredHB.txt'));
plate = load_hcb_model(paths, struct('dof_per_node', 3, 'measured_dof_index', 2));

%% FRFs — Rayleigh damping; experimental via stiffness shift (model error)
pf_start = 15;  pf_end = 40;  pf_step = 0.1;     % FRF band [Hz]
shift_value     = 1.04;                          % experimental stiffness shift
plate_noise_seed = 10;

[wp, wp_hz] = frequency_generation(pf_start, pf_end, pf_step);
Cp = damping(plate.K, plate.M, ray);             % Rayleigh, same config as the beam

Yp  = compute_frf(wp, plate.K,             plate.M, zeros(size(plate.K)), 'accelerance');
YEp = compute_frf(wp, shift_value*plate.K, plate.M, Cp,                   'accelerance');

ph  = plate.physical_DoF_array;
Yp  = Yp(ph, ph, :);
YEp = YEp(ph, ph, :);
YEp = add_noise(true, YEp, 1e-3, 1e-3, 1e-3, 1e-3, plate_noise_seed);

mz  = plate.measured_dof_array;                  % transverse (measured) DoF per node
Yz  = Yp(mz, mz, :);
YEz = YEp(mz, mz, :);

%% Optimal sensor/excitation search (NDP, 2 sensors + 3 excitations, all DoFs)
% NOTE: NDP(2,3) over all 36 DoFs is a large search (nchoosek(36,2)*nchoosek(34,3)
%       ~ 3.8M SEMM evaluations) — runs under parfor; expect a long run.
np = plate.n;
optP = struct('num_sensors', num_sensors, 'search_type', 'NDP', ...
    'extra_excitations', extra_excitations, ...
    'candidate_dofs', 1:np, 'validation_dofs', 1:np, ...
    'methods', 'COH', 'extrema', 'max', 'frequency_range', wp, 'verbose', true);
res_p = exhaustive_search(Yz, YEz, optP);

fprintf('\n[PLATE] sensors %s | excitations %s | GCCM %.4f\n', ...
    mat2str(res_p.sensor_dofs), mat2str(res_p.excitation_dofs), res_p.cor_overall);

%% Plot the optimal sensors / excitations on the plate
plot_square_plate(plate.NodeCoords, res_p.sensor_dofs, res_p.excitation_dofs);
