% MAIN  Optimal interface expansion — worked examples.
%
%   Two self-contained, end-to-end examples for the methodology in:
%     Junaid et al. (2026). Journal of Sound and Vibration. DOI: 10.1016/j.jsv.2026.119782
%
%   Each example: build the model -> FRFs (+ noise) -> find the optimal
%   sensor/excitation DoFs by BOTH exhaustive search and the Mountain Gazelle
%   Optimizer (mgo) on the same NDP space -> compare the two searches.
%
%   EXAMPLE 1 — Beam (substructure coupling):
%     Two steel cantilevers (A fixed-left, B fixed-right), Rayleigh damping.
%     After the optimal search, the SEMM expansions are coupled at the interface
%     and compared against the full coupled model (GCCM, Gamma).
%
%   EXAMPLE 2 — Square plate (standalone, NO coupling):
%     An Ansys Craig-Bampton reduced model; optimal sensor/excitation placement
%     over all valid DoFs (pure optimal expansion), with the placement plotted.
%
%   Dependencies (utils/): create_cantilever_beam, load_hcb_model, damping,
%     frequency_generation, compute_frf, add_noise, exhaustive_search,
%     objective_function, mgo, semm, func_coh, svd_truncation, primal_coupling,
%     couple_substructures, plot_square_plate.

clear; clc; close all;

if isempty(which('create_cantilever_beam'))
    addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));
end

%% Shared search / damping configuration

ray.type   = 'proportional';   % Rayleigh damping
ray.alpha  = 1.0;              % mass-proportional coefficient
ray.beta   = 1e-5;            % stiffness-proportional coefficient (≈0.4–2% over the bands used)

num_sensors       = 1;        % NDP: 1 sensor ...
extra_excitations = 0;        % ... + (1+0) = 1 excitation → 2 measured DoFs (small NDP, fast)

mgo_agents   = 5;             % mgo population size
mgo_max_iter = 50;            % mgo iterations
mgo_seed     = 42;            % RNG seed (philox) for reproducible mgo runs


%% ========================================================================
%  EXAMPLE 1 — BEAM: optimal search (exhaustive + mgo) + interface coupling
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

%% Exhaustive search (NDP) on each beam
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

%% MGO search on beam A — identical NDP space, solved with mgo.m directly
cand_a  = setdiff(1:n, interface_a);
k_exc   = num_sensors + extra_excitations;
r_combs = nchoosek(cand_a, num_sensors);
n_outer = size(r_combs, 1);
n_inner = nchoosek(numel(cand_a) - num_sensors, k_exc);
e_combs = zeros(n_inner, k_exc, n_outer);
for i = 1:n_outer
    e_combs(:, :, i) = nchoosek(setdiff(cand_a, r_combs(i, :)), k_exc);
end

% Maximise GCCM -> minimise its negative (objective_function: semm + func_coh).
fobj_a = @(x) -objective_function(x, wb, YA, YEA, r_combs, e_combs, ...
    false, 0, false, interface_a, 'COH', 3, []);

rng(mgo_seed, 'philox');
[mgo_score_a, mgo_pos_a] = mgo(mgo_agents, mgo_max_iter, [1 1], [n_outer n_inner], 2, fobj_a);
mgo_sensors_a     = r_combs(round(mgo_pos_a(1)), :);
mgo_excitations_a = e_combs(round(mgo_pos_a(2)), :, round(mgo_pos_a(1)));

%% Compare exhaustive vs mgo (beam A)
fprintf('\n=== BEAM A: exhaustive vs MGO (same NDP space) ===\n');
fprintf('  exhaustive: sensors %s exc %s GCCM %.6f\n', ...
    mat2str(res_a.sensor_dofs), mat2str(res_a.excitation_dofs), res_a.cor_overall);
fprintf('  MGO       : sensors %s exc %s GCCM %.6f\n', ...
    mat2str(mgo_sensors_a), mat2str(mgo_excitations_a), -mgo_score_a);
fprintf('  gap (exhaustive - MGO) = %.2e\n', res_a.cor_overall - (-mgo_score_a));

%% Full coupled (reference) model — primal coupling of the experimental beams
[~, k_full, m_full, c_full] = primal_coupling(n_interface, ...
    beam_a_exp.K, beam_b_exp.K, beam_a_exp.M, beam_b_exp.M, beam_a_exp.C, beam_b_exp.C);
y_full = compute_frf(wb, k_full, m_full, c_full, 'accelerance');
y_full = y_full(1:2:length(k_full), 1:2:length(k_full), :);

%% Optimally coupled SEMM model (transverse interface) + comparison vs full
ys_a_t   = svd_truncation(res_a.ys, 15);         % rank-reduce noisy SEMM before coupling
ys_b_t   = svd_truncation(res_b.ys, 15);
y_couple = couple_substructures(ys_a_t, ys_b_t, n_interface, 'transverse');

nf = min(size(y_couple, 3), size(y_full, 3));
[~, ~, gamma_couple] = func_coh(y_couple(:, :, 1:nf), y_full(:, :, 1:nf));
fprintf('\n[BEAM] coupled-SEMM vs full-model GCCM (Gamma) = %.6f\n', gamma_couple);


%% ========================================================================
%  EXAMPLE 2 — SQUARE PLATE: standalone optimal placement (NO coupling)
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
shift_value      = 1.04;                         % experimental stiffness shift
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

%% Exhaustive search (NDP) over all valid DoFs
np = plate.n;
optP = struct('num_sensors', num_sensors, 'search_type', 'NDP', ...
    'extra_excitations', extra_excitations, ...
    'candidate_dofs', 1:np, 'validation_dofs', 1:np, ...
    'methods', 'COH', 'extrema', 'max', 'frequency_range', wp, 'verbose', true);
res_p = exhaustive_search(Yz, YEz, optP);

%% MGO search on the plate — identical NDP space, solved with mgo.m directly
k_exc_p   = num_sensors + extra_excitations;
r_combs_p = nchoosek(1:np, num_sensors);
n_outer_p = size(r_combs_p, 1);
n_inner_p = nchoosek(np - num_sensors, k_exc_p);
e_combs_p = zeros(n_inner_p, k_exc_p, n_outer_p);
for i = 1:n_outer_p
    e_combs_p(:, :, i) = nchoosek(setdiff(1:np, r_combs_p(i, :)), k_exc_p);
end

fobj_p = @(x) -objective_function(x, wp, Yz, YEz, r_combs_p, e_combs_p, ...
    false, 0, false, 1:np, 'COH', 3, []);

rng(mgo_seed, 'philox');
[mgo_score_p, mgo_pos_p] = mgo(mgo_agents, mgo_max_iter, [1 1], [n_outer_p n_inner_p], 2, fobj_p);
mgo_sensors_p     = r_combs_p(round(mgo_pos_p(1)), :);
mgo_excitations_p = e_combs_p(round(mgo_pos_p(2)), :, round(mgo_pos_p(1)));

%% Compare exhaustive vs mgo (plate)
fprintf('\n=== SQUARE PLATE: exhaustive vs MGO (same NDP space) ===\n');
fprintf('  exhaustive: sensors %s exc %s GCCM %.6f\n', ...
    mat2str(res_p.sensor_dofs), mat2str(res_p.excitation_dofs), res_p.cor_overall);
fprintf('  MGO       : sensors %s exc %s GCCM %.6f\n', ...
    mat2str(mgo_sensors_p), mat2str(mgo_excitations_p), -mgo_score_p);
fprintf('  gap (exhaustive - MGO) = %.2e\n', res_p.cor_overall - (-mgo_score_p));

%% Plot the optimal sensors / excitations on the plate (exhaustive result)
plot_square_plate(plate.NodeCoords, res_p.sensor_dofs, res_p.excitation_dofs);
