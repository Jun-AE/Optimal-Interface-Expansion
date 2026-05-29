%% BEAM_MAIN_TEST — profiling, performance, and validation driver
%
% Mirrors the beam_coh_main workflow with PLAN sprint items applied here only:
%   A1  fixed RNG seeds (noise + MGO)
%   B4  do_plots gate (default false — fast runs)
%   E2  validate_interface_dofs at coupling checkpoints
%   E4  no figures inside search / MGO loops
%
% beam_coh_main.m is left unchanged. Utils updated for B2 (compute_frf) and E1 (svd_truncation).
%
% Usage:
%   addpath('utils')
%   BEAM_MAIN_TEST
%
% Optional: set profile_run = true, then profile viewer after completion.

clear; clc; close all;

if isempty(which('create_cantilever_beam'))
    addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));
end

%% Test / profiling controls (PLAN A1, B4)

do_plots         = false;   % false = skip all non-essential figures
rng_seed_noise   = 12;      % add_noise(..., seed)
rng_seed_mgo     = 12;      % rng before mgo(...)
profile_run      = false;  % true → profile on; results in profile viewer
test_exhaustive  = true;   % true → run exhaustive_search validation (beam)
test_square_plate= false;  % true → validate exhaustive_search + load_hcb_model (square plate, full NDP ~29 min)

if profile_run
    profile on;
end

t_total = tic;
fprintf('\n=== BEAM_MAIN_TEST (do_plots=%d) ===\n', do_plots);

%% Beam properties (shared by Beams A and B)

props.Length = 3;
props.nelm   = 20;
props.E      = 1.93e11;
props.b      = 0.05;
props.h      = 0.2;
props.rho    = 7872;

gamma_a = 1.05;
gamma_b = 1.05;

assembly.damping.type   = 'modal';
assembly.damping.alpha  = 0.1;
assembly.damping.beta   = 1e-6;
assembly.damping.ratios = [];

%% Build cantilever beam pairs

t_section = tic;
[beam_a, beam_a_exp] = create_cantilever_beam(props, gamma_a, 'left');
[beam_b, beam_b_exp] = create_cantilever_beam(props, gamma_b, 'right');

vec_a      = beam_a.modes;
vec_a_exp  = beam_a_exp.modes;
freq_a     = beam_a.freq;
freq_a_exp = beam_a_exp.freq;
freq_b     = beam_b.freq;
freq_b_exp = beam_b_exp.freq;

n_interface_nodes = 3;
validate_interface_dofs(n_interface_nodes, beam_a.ndof, beam_b.ndof, ...
    'Layout', 'full', 'Caller', 'BEAM_MAIN_TEST:beams-built');
fprintf('  build beams: %.2f s\n', toc(t_section));

if isequal(beam_a.ndof, beam_b.ndof)
    ndof_equal_flag = true;
    disp('Both beams have the same number of DoFs.');
else
    ndof_equal_flag = false;
end

%% Natural frequency comparison (optional plot)

if do_plots
    n_modes = 4;
    figure;
    subplot(2, 2, 1);
    bar([freq_a(1:n_modes), freq_a_exp(1:n_modes)]);
    title('Natural Frequencies of Beam A');
    xlabel('Mode Number'); ylabel('Frequency (Hz)');
    legend('Numerical', 'Experimental');
    subplot(2, 2, 2);
    bar([freq_b(1:n_modes), freq_b_exp(1:n_modes)]);
    title('Natural Frequencies of Beam B');
    xlabel('Mode Number'); ylabel('Frequency (Hz)');
    legend('Numerical', 'Experimental');
    subplot(2, 2, 3);
    bar(abs(freq_a(1:n_modes) - freq_a_exp(1:n_modes)));
    title('Difference in Natural Frequencies of Beam A');
    xlabel('Mode Number'); ylabel('Frequency Difference (Hz)');
    subplot(2, 2, 4);
    bar(abs(freq_b(1:n_modes) - freq_b_exp(1:n_modes)));
    title('Difference in Natural Frequencies of Beam B');
    xlabel('Mode Number'); ylabel('Frequency Difference (Hz)');
end

%% Coupled experimental model (primal coupling + eigen)

t_section = tic;
assembly.damping.ratios = 0.002 * ones(1, length(beam_a.K));
beam_a_exp.C = damping(beam_a_exp.K, beam_a_exp.M, assembly.damping);
beam_b_exp.C = damping(beam_b_exp.K, beam_b_exp.M, assembly.damping);

[L, k_full, m_full, c_full] = primal_coupling(n_interface_nodes, ...
    beam_a_exp.K, beam_b_exp.K, beam_a_exp.M, beam_b_exp.M, ...
    beam_a_exp.C, beam_b_exp.C);

[vec_exp, fsol_exp] = eig(k_full, m_full);
if do_plots
    [~, ~] = plot_mode_shape(vec_exp, fsol_exp, 1, 3);
else
  % E4: skip plot_mode_shape figure; eigenvalues only if needed later
end
fprintf('  primal coupling + eig: %.2f s\n', toc(t_section));

%% FRF generation

t_section = tic;
beam_a.freqstart = 1;
beam_a.freqend   = 700;
beam_a.stepsize  = 1;
[beam_a.frequency_range, beam_a.frequency_range_Hz] = ...
    frequency_generation(beam_a.freqstart, beam_a.freqend, beam_a.stepsize);

beam_b.freqstart = 1;
beam_b.freqend   = 700;
beam_b.stepsize  = 1;
[beam_b.frequency_range, beam_b.frequency_range_Hz] = ...
    frequency_generation(beam_b.freqstart, beam_b.freqend, beam_b.stepsize);

assembly.damping.ratios = 0 * ones(1, 10);
beam_a.C = damping(beam_a.K, beam_a.M, assembly.damping);
beam_b.C = damping(beam_b.K, beam_b.M, assembly.damping);
assembly.shift = false;

YA  = compute_frf(beam_a.frequency_range, beam_a.K, beam_a.M, beam_a.C, 'accelerance');
YB  = compute_frf(beam_b.frequency_range, beam_b.K, beam_b.M, beam_b.C, 'accelerance');
YEA = compute_frf(beam_a.frequency_range, beam_a_exp.K, beam_a_exp.M, beam_a_exp.C, 'accelerance');
YEB = compute_frf(beam_b.frequency_range, beam_b_exp.K, beam_b_exp.M, beam_b_exp.C, 'accelerance');
fprintf('  FRF generation (4 tensors): %.2f s\n', toc(t_section));

fs        = 1800;
num_modes = 3;

if ndof_equal_flag
    translational_dof_array = 1:2:beam_a.ndof;
end

if do_plots
    [modes, dampr, mode_shapes] = plot_modal_waterfall( ...
        YEA, 10, beam_a.frequency_range_Hz, fs, ...
        translational_dof_array, num_modes, ...
        vec_a_exp(:, 1:num_modes), freq_a_exp(1:num_modes));
end

YA  = YA(translational_dof_array, translational_dof_array, :);
YB  = YB(translational_dof_array, translational_dof_array, :);
YEA = YEA(translational_dof_array, translational_dof_array, :);
YEB = YEB(translational_dof_array, translational_dof_array, :);

noise_flag = true;
YEA = add_noise(noise_flag, YEA, 0.005, 0.005, 1e-8, 1e-8, rng_seed_noise);
YEB = add_noise(noise_flag, YEB, 0.005, 0.005, 1e-8, 1e-8, rng_seed_noise);

if do_plots
    plot_frf(beam_b.frequency_range_Hz, 4, 10, YB, YEB, '{Numerical}^B', '{Experimental}^B');
    plot_frf(beam_a.frequency_range_Hz, 11, 17, YA, YEA, '{Numerical}^A', '{Experimental}^A');
end

assembly.method = 'COH';
switch assembly.method
    case 'COH'
        [~, assembly.cor_matrix_b, ~] = func_coh(YB, YEB);
        [~, assembly.cor_matrix_a, ~] = func_coh(YA, YEA);
    case 'LAC'
        [~, assembly.cor_matrix_b, ~] = func_lac(YB, YEB);
        [~, assembly.cor_matrix_a, ~] = func_lac(YA, YEA);
end

if do_plots
    correlation_plot(assembly.cor_matrix_b);
    correlation_plot(assembly.cor_matrix_a);
end

%% Energy-based placement + modal ID (diagnostics optional)

numsens = 4;
e_a     = analyze_frf_energy(YA, beam_a.frequency_range_Hz, numsens, 3, do_plots);

best_sensor_dofs     = e_a.top_sensors;
best_excitation_dofs = e_a.top_excitations;
n_interface          = n_interface_nodes;

if do_plots
    plot_beam(n_interface, beam_a.nnode, ...
        round(best_sensor_dofs), round(best_excitation_dofs), 'A');
end

[nat_freq_beam_a, phi_beam_a] = estimate_modal_parameters(YA, beam_a.frequency_range_Hz, ...
    'NumModes', 3, 'ProminenceThreshold', 0, 'MinPeakDistanceHz', 1);
[nat_freq_beam_a_exp, phi_beam_a_exp] = estimate_modal_parameters(YEA, beam_a.frequency_range_Hz, ...
    'NumModes', 3, 'ProminenceThreshold', 0, 'MinPeakDistanceHz', 1);

mac_matrix = calculate_mac(phi_beam_a, phi_beam_a_exp, do_plots);

if do_plots
    node_coords = [1:size(YA, 1); zeros(1, size(YA, 1))]';
    visualize_mode_shapes(node_coords, phi_beam_a, nat_freq_beam_a, ...
        'Animate', false, 'SaveVideo', false, ...
        'VideoFile', 'beamModesAnimated.mp4');
    close all;
end

%% Exhaustive search (E4: no plots in bruteforce_search_beam)

t_section = tic;
num_sensors = 2;

opts.frequency_range = beam_a.frequency_range;
opts.method          = 'COH';
opts.search_type     = 'DP';
opts.truncation      = false;
opts.reduction       = 0;
opts.ec              = 1;
opts.trust_func      = false;
opts.trust_func_val  = [];
opts.extrema         = 'max';

interface_dofs_a = (size(YA, 1) - n_interface + 1) : size(YA, 1);
interface_dofs_b = 1 : n_interface;

validate_interface_dofs(n_interface_nodes, size(YA, 1), size(YB, 1), ...
    'Layout', 'translational', ...
    'InterfaceDofsA', interface_dofs_a, ...
    'InterfaceDofsB', interface_dofs_b, ...
    'Caller', 'BEAM_MAIN_TEST:translational-FRF');

out_a = bruteforce_search_beam(interface_dofs_a, num_sensors, YA, YEA, opts);
opts.frequency_range = beam_b.frequency_range;
out_b = bruteforce_search_beam(interface_dofs_b, num_sensors, YB, YEB, opts);
fprintf('  brute-force search (A+B): %.2f s\n', toc(t_section));

sensor_location_1    = out_a.sensor_dofs;
excitation_1         = out_a.excitation_dofs;
interface_cor_oval_1 = out_a.interface_cor_oval;
interface_cor_2d_1   = out_a.interface_cor_2d;
ys_a                 = out_a.ys;

sensor_location_2    = out_b.sensor_dofs;
excitation_2         = out_b.excitation_dofs;
interface_cor_2d_2   = out_b.interface_cor_2d;
ys_b                 = out_b.ys;

if do_plots
    correlation_plot(interface_cor_2d_1, 'Coherence: ys_a and YEA');
    correlation_plot(interface_cor_2d_2, 'Coherence: ys_b and YEB');
    plot_frf(beam_a.frequency_range_Hz, 18, 18, YEA, ys_a, 'semm', 'Experimental');
    plot_frf(beam_a.frequency_range_Hz, 2, 2, YB, ys_b, 'Numerical', 'semm', YEB, 'Experimental');
    close all;
    plot_beam(n_interface, beam_a.nnode, sensor_location_1, excitation_1, 'A');
    plot_beam(n_interface, beam_b.nnode, sensor_location_2, excitation_2, 'B');
end

%% MGO (A1: fixed seed; E4: mgo plot_flag stays default false)

t_section = tic;
numvars       = 2;
num_sensors_m = 1;
max_iter      = 50;
mgo_agents    = 5;

dofs           = 1:size(YA, 1);
interface_dofs = dofs;

r_combs = nchoosek(dofs, num_sensors_m);
n_outer = size(r_combs, 1);
n_inner = nchoosek(length(dofs) - num_sensors_m, num_sensors_m + 1);
e_combs_ndp = zeros(n_inner, num_sensors_m + 1, n_outer);
for i = 1:n_outer
    e_combs_ndp(:, :, i) = nchoosek(setdiff(dofs, r_combs(i, :)), num_sensors_m + 1);
end

fobj = @(x) -objective_function(x, beam_a.frequency_range, YA, YEA, ...
    r_combs, e_combs_ndp, ...
    false, 0, false, interface_dofs, ...
    'COH', 3, 130.125);

lb = [1, 1];
ub = [size(r_combs, 1), size(e_combs_ndp, 1)];

rng(rng_seed_mgo, 'philox');
[mgo_score, mgo_pos, ~, mgo_stats] = mgo(mgo_agents, max_iter, lb, ub, numvars, fobj);

mgo_sensor_dofs     = r_combs(round(mgo_pos(1)), :);
mgo_excitation_dofs = e_combs_ndp(round(mgo_pos(2)), :, round(mgo_pos(1)));

fprintf('  MGO + verification parfor: %.2f s\n', toc(t_section));
fprintf('\nMGO best  GCCM = %.6g\n', -mgo_score);
fprintf('MGO sensors    = %s\n', mat2str(mgo_sensor_dofs));
fprintf('MGO excitations= %s\n', mat2str(mgo_excitation_dofs));

%% MGO vs brute-force (console only — E4)

[ii, jj] = ndgrid(1:n_outer, 1:n_inner);
ii = ii(:);
jj = jj(:);
n_total = numel(ii);

fprintf('\nBrute-forcing the same NDP search space (%d combinations)...\n', n_total);
cor_flat = -inf(n_total, 1);
parfor k = 1:n_total
    cor_flat(k) = -fobj([ii(k), jj(k)]);
end

[bf_best_gccm, idx_lin] = max(cor_flat);
bf_sensor_dofs     = r_combs(ii(idx_lin), :);
bf_excitation_dofs = e_combs_ndp(jj(idx_lin), :, ii(idx_lin));

mgo_gccm = -mgo_score;
gap_abs  = bf_best_gccm - mgo_gccm;
gap_pct  = 100 * gap_abs / bf_best_gccm;

fprintf('\n=== MGO vs Brute-Force (identical NDP search space) ===\n');
fprintf('                 %-22s %-22s\n', 'Brute-force (global)', 'MGO');
fprintf('  sensors     :  %-22s %-22s\n', mat2str(bf_sensor_dofs), mat2str(mgo_sensor_dofs));
fprintf('  excitations :  %-22s %-22s\n', mat2str(bf_excitation_dofs), mat2str(mgo_excitation_dofs));
fprintf('  best GCCM   :  %-22.6f %-22.6f\n', bf_best_gccm, mgo_gccm);
fprintf('  gap (BF-MGO):  %.6f  (%.3f%%)\n', gap_abs, gap_pct);
fprintf('  evaluations :  %-22d %-22d\n', n_total, mgo_stats.n_fobj_eval);

%% Coupled system comparison

t_section = tic;
freqstart = 1;
freqend   = 185;
stepsize  = 1;

[frequency_range, frequency_range_Hz] = frequency_generation(freqstart, freqend, stepsize);
translational_dof_array = 1:2:length(k_full);

y_full = compute_frf(frequency_range, k_full, m_full, c_full, 'accelerance');
y_full = y_full(translational_dof_array, translational_dof_array, :);

y_full = add_noise(false, y_full, 0.005, 0.005, 1e-8, 1e-8, rng_seed_noise);

% E1: num_drop = 15 trailing singular values discarded per page
ys_a_c = svd_truncation(ys_a, 15);
ys_b_c = svd_truncation(ys_b, 15);

validate_interface_dofs(n_interface_nodes, size(ys_a_c, 1), size(ys_b_c, 1), ...
    'Layout', 'translational', 'Caller', 'BEAM_MAIN_TEST:pre-couple');

y_couple = couple_substructures(ys_a_c, ys_b_c, n_interface_nodes, 'transverse');
[~, cor_mat_couple, cor_couple] = func_coh(y_couple(:, :, 1:size(y_full, 3)), y_full);

fprintf('  coupled validation: %.2f s\n', toc(t_section));
fprintf('  final coupled GCCM (Gamma) = %.6f\n', cor_couple);

if do_plots
    correlation_plot(cor_mat_couple, 'COH');
    plot_frf(frequency_range_Hz, 23, 23, y_couple(:, :, 1:size(y_full, 3)), y_full, 'Coupled', 'Experimental');
    close all;
end

fprintf('\n=== BEAM_MAIN_TEST total: %.2f s ===\n', toc(t_total));

%% exhaustive_search validation (BEAM structure)
%
% Validates the generic utils/exhaustive_search.m against the proven beam
% references already computed above:
%   - bruteforce_search_beam (out_a)          — DP search
%   - the MGO-vs-BF brute force (bf_* / 3420)  — NDP global search
% Square plate and coupled-plate structures need experimental data + out-of-scope
% deps; their test procedure is documented in PLAN.md and run separately.

if test_exhaustive
    fprintf('\n=== exhaustive_search validation (BEAM) ===\n');
    verdict = {'FAIL', 'PASS'};
    tol_ref = 1e-6;     % vs independent reference implementation

    % --- A: DP (Beam A) vs bruteforce_search_beam (out_a) ---
    es = struct();
    es.num_sensors     = num_sensors;                       % 2 (from brute-force section)
    es.search_type     = 'DP';
    es.validation_dofs = interface_dofs_a;
    es.candidate_dofs  = setdiff(1:size(YA, 1), interface_dofs_a);
    es.frequency_range = beam_a.frequency_range;
    es.methods         = 'COH';
    es.extrema         = 'max';
    es.verbose         = false;
    a = exhaustive_search(YA, YEA, es);
    okA = isequal(sort(a.sensor_dofs), sort(out_a.sensor_dofs)) && ...
          abs(a.cor_overall - out_a.interface_cor_oval) < tol_ref;
    fprintf('  [%s] A  DP Beam A vs bruteforce_search_beam (sensors %s, GCCM %.6f)\n', ...
            verdict{1 + okA}, mat2str(a.sensor_dofs), a.cor_overall);

    % --- B: DP (Beam B) vs bruteforce_search_beam (out_b) ---
    es.validation_dofs = interface_dofs_b;
    es.candidate_dofs  = setdiff(1:size(YB, 1), interface_dofs_b);
    es.frequency_range = beam_b.frequency_range;
    b = exhaustive_search(YB, YEB, es);
    okB = isequal(sort(b.sensor_dofs), sort(out_b.sensor_dofs)) && ...
          abs(b.cor_overall - out_b.interface_cor_oval) < tol_ref;
    fprintf('  [%s] B  DP Beam B vs bruteforce_search_beam (sensors %s, GCCM %.6f)\n', ...
            verdict{1 + okB}, mat2str(b.sensor_dofs), b.cor_overall);

    % --- C: NDP global vs MGO-section brute force (1 sensor, all DoFs) ---
    % Mirrors the MGO-section search space exactly (3420 combos).
    es.num_sensors       = 1;
    es.search_type       = 'NDP';
    es.extra_excitations = 1;
    es.candidate_dofs    = 1:size(YA, 1);
    es.validation_dofs   = 1:size(YA, 1);
    es.frequency_range   = beam_a.frequency_range;
    tt = tic;
    nd = exhaustive_search(YA, YEA, es);
    t_ndp = toc(tt);
    okC = isequal(nd.sensor_dofs, bf_sensor_dofs) && ...
          isequal(sort(nd.excitation_dofs), sort(bf_excitation_dofs)) && ...
          abs(nd.cor_overall - bf_best_gccm) < tol_ref;
    fprintf('  [%s] C  NDP global vs MGO-section brute force (sensors %s, exc %s, GCCM %.6f, %.2f s)\n', ...
            verdict{1 + okC}, mat2str(nd.sensor_dofs), ...
            mat2str(nd.excitation_dofs), nd.cor_overall, t_ndp);

    % --- D: multi-method returns one struct per method ---
    es.num_sensors     = num_sensors;
    es.search_type     = 'DP';
    es.candidate_dofs  = setdiff(1:size(YA, 1), interface_dofs_a);
    es.validation_dofs = interface_dofs_a;
    es.methods         = {'COH', 'LAC'};
    rm = exhaustive_search(YA, YEA, es);
    okD = numel(rm) == 2 && strcmp(rm(1).method, 'COH') && strcmp(rm(2).method, 'LAC');
    fprintf('  [%s] D  multi-method returns %d entries\n', verdict{1 + okD}, numel(rm));

    all_ok = okA && okB && okC && okD;
    fprintf('=== exhaustive_search BEAM: %s ===\n', verdict{1 + all_ok});
end

%% SQUARE PLATE validation (exhaustive_search + load_hcb_model)
%
% Validates both new functions against the saved experimental run
% Square_Exhaustive_Result.mat (a MainSquare.m workspace dump):
%   A) exhaustive_search reproduces bruteforceSearchSquare's selected DoFs and
%      GCCM, run on the *exact* saved Yz/YEz (so noise realisation matches).
%   B) load_hcb_model reproduces the experimental sortKM-reordered K/M (and the
%      physical / measured DOF index sets) directly from the Ansys HCB files.
% Part B self-skips when the Ansys source files are not on this machine.
% Block self-skips (no error) if the .mat is absent.

if test_square_plate
    fprintf('\n=== validation (SQUARE PLATE) ===\n');
    verdict = {'FAIL', 'PASS'};
    sq_file = fullfile(fileparts(mfilename('fullpath')), 'Square_Exhaustive_Result.mat');

    if ~isfile(sq_file)
        fprintf('  (Square_Exhaustive_Result.mat not found — skipping)\n');
    else
        S    = load(sq_file);
        Yz   = S.PlateA.Yz;
        YEz  = S.PlateA.YEz;
        n_sq = size(Yz, 1);

        % --- A: exhaustive_search vs bruteforceSearchSquare (saved Yz/YEz) ---
        o = struct();
        o.num_sensors       = S.num_sensors;        % 2
        o.search_type       = S.BruteForce_Type;    % 'NDP'
        o.extra_excitations = S.c;                  % 0  (excitation size = num_sensors)
        o.candidate_dofs    = 1:n_sq;
        o.validation_dofs   = 1:n_sq;
        o.methods           = 'COH';
        o.extrema           = S.extrema_type;        % 'max'
        o.frequency_range   = S.Assembly.frequency_range;
        o.verbose           = false;

        res = exhaustive_search(Yz, YEz, o);

        ok_sens = isequal(sort(res.sensor_dofs(:).'),     sort(double(S.response(:).')));
        ok_exc  = isequal(sort(res.excitation_dofs(:).'), sort(double(S.excitation(:).')));
        val_gap = abs(res.cor_overall - S.InterfaceCorOVal);
        okA = ok_sens && ok_exc;
        fprintf('  [%s] A  exhaustive_search vs bruteforceSearchSquare\n', verdict{1 + okA});
        fprintf('         sensors     : got %-9s ref %s\n', mat2str(res.sensor_dofs),     mat2str(S.response));
        fprintf('         excitations : got %-9s ref %s\n', mat2str(res.excitation_dofs), mat2str(S.excitation));
        fprintf('         GCCM        : got %.6f  ref %.6f  (gap %.2e)\n', ...
                res.cor_overall, S.InterfaceCorOVal, val_gap);

        % --- B: load_hcb_model reproduces experimental sortKM K/M ---
        data_dir = fullfile(fileparts(mfilename('fullpath')), 'Data');
        paths  = struct( ...
            'nodes',   fullfile(data_dir, 'Nodes_36.txt'), ...
            'mapping', fullfile(data_dir, 'KredHB.mapping'), ...
            'K',       fullfile(data_dir, 'KredHB.txt'), ...
            'M',       fullfile(data_dir, 'MredHB.txt'));

        if all(cellfun(@isfile, struct2cell(paths)))
            m = load_hcb_model(paths, struct('dof_per_node', 3, 'measured_dof_index', 2, 'verbose', false));
            relK = max(abs(m.K(:) - S.PlateA.K(:))) / max(abs(S.PlateA.K(:)));
            relM = max(abs(m.M(:) - S.PlateA.M(:))) / max(abs(S.PlateA.M(:)));
            ok_meas = isequal(double(m.measured_dof_array(:).'), double(S.Assembly.Zaxis_DoF_Array(:).'));
            okB = relK < 1e-10 && relM < 1e-10 && ...
                  m.NCMSMode == S.Assembly.NCMSMode && ok_meas && m.n == n_sq;
            fprintf('  [%s] B  load_hcb_model vs experimental sortKM\n', verdict{1 + okB});
            fprintf('         relK %.1e  relM %.1e  NCMSMode %d (ref %d)  n %d  measured-DOFs match %d\n', ...
                    relK, relM, m.NCMSMode, S.Assembly.NCMSMode, m.n, ok_meas);
        else
            okB = true;   % not applicable on this machine
            fprintf('  [skip] B  Ansys source files not found — load_hcb_model check skipped\n');
        end

        fprintf('=== SQUARE PLATE: %s ===\n', verdict{1 + (okA && okB)});
    end
end

if profile_run
    profile off;
    profile viewer;
end
