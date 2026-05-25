clear; clc; close all;

%% Beam properties (shared by Beams A and B)

props.Length = 3;       % [m]      Beam length
props.nelm   = 20;      %          Number of elements
props.E      = 1.93e11; % [Pa]     Young's modulus
props.b      = 0.05;    % [m]      Breadth
props.h      = 0.2;     % [m]      Height
props.rho    = 7872;    % [kg/m^3] Density

%% Stiffness scaling (numerical -> experimental)

gamma_a = 1.05;
gamma_b = 1.05;

%% Damping configuration

% Single config struct consumed by damping(). Switch formulation by changing
% assembly.damping.type. Per-call-site overrides (e.g. .ratios) are applied
% just before each damping() call.
assembly.damping.type   = 'modal';   % 'modal' or 'proportional'
assembly.damping.alpha  = 0.1;       % Rayleigh mass-proportional coefficient
assembly.damping.beta   = 1e-6;      % Rayleigh stiffness-proportional coefficient
assembly.damping.ratios = [];        % set per call site for the 'modal' case

%% Build cantilever beam pairs (numerical + experimental)

% Beam A: fixed at left end.   Beam B: fixed at right end.
[beam_a, beam_a_exp] = create_cantilever_beam(props, gamma_a, 'left');
[beam_b, beam_b_exp] = create_cantilever_beam(props, gamma_b, 'right');

% Legacy variable aliases used by downstream sections of this script.
vec_a      = beam_a.modes;
vec_a_exp  = beam_a_exp.modes;
freq_a     = beam_a.freq;
freq_a_exp = beam_a_exp.freq;
freq_b     = beam_b.freq;
freq_b_exp = beam_b_exp.freq;

%% DoF parity check between the two beams

if isequal(beam_a.ndof, beam_b.ndof)
    ndof_equal_flag = true;
    disp('Both beams have the same number of DoFs.');
end

%% Natural frequency comparison plot (first 4 modes)

n_modes = 4;
figure;

subplot(2,2,1);
bar([freq_a(1:n_modes), freq_a_exp(1:n_modes)]);
title('Natural Frequencies of Beam A');
xlabel('Mode Number');  ylabel('Frequency (Hz)');
legend('Numerical', 'Experimental');

subplot(2,2,2);
bar([freq_b(1:n_modes), freq_b_exp(1:n_modes)]);
title('Natural Frequencies of Beam B');
xlabel('Mode Number');  ylabel('Frequency (Hz)');
legend('Numerical', 'Experimental');

subplot(2,2,3);
bar(abs(freq_a(1:n_modes) - freq_a_exp(1:n_modes)));
title('Difference in Natural Frequencies of Beam A');
xlabel('Mode Number');  ylabel('Frequency Difference (Hz)');

subplot(2,2,4);
bar(abs(freq_b(1:n_modes) - freq_b_exp(1:n_modes)));
title('Difference in Natural Frequencies of Beam B');
xlabel('Mode Number');  ylabel('Frequency Difference (Hz)');




%% assembly of both Beams for Experimental Model

% Damping Matrix C (Experimental beams) — zeta = 0.2% per mode (if 'modal')
assembly.damping.ratios = 0.002 * ones(1, length(beam_a.K));
beam_a_exp.C = damping(beam_a_exp.K, beam_a_exp.M, assembly.damping);
beam_b_exp.C = damping(beam_b_exp.K, beam_b_exp.M, assembly.damping);

% Primal Coupling for Experimental Beams
[L, k_full, m_full, c_full] = primal_coupling(3, beam_a_exp.K, beam_b_exp.K, beam_a_exp.M, beam_b_exp.M, beam_a_exp.C, beam_b_exp.C);

% Eigen Solution for Experimental Beams
[vec_exp, fsol_exp] = eig(k_full, m_full);
[v_exp, freq_exp] = plot_mode_shape(vec_exp, fsol_exp, 1, 3);


%% FRF Generation

% Frequency generation for beam_a
beam_a.freqstart = 1;
beam_a.freqend = 700;
beam_a.stepsize = 1;
[beam_a.frequency_range, beam_a.frequency_range_Hz] = frequency_generation(beam_a.freqstart, beam_a.freqend, beam_a.stepsize);

% Frequency generation for beam_b
beam_b.freqstart = 1;
beam_b.freqend = 700;
beam_b.stepsize = 1;
[beam_b.frequency_range, beam_b.frequency_range_Hz] = frequency_generation(beam_b.freqstart, beam_b.freqend, beam_b.stepsize);

% Damping Matrix C (Numerical beams) — Zero damping for the first 10 modes
assembly.damping.ratios = 0 * ones(1, 10);
beam_a.C = damping(beam_a.K, beam_a.M, assembly.damping);
beam_b.C = damping(beam_b.K, beam_b.M, assembly.damping);
assembly.shift = false;

% Full Harmonic FRF - Numerical (accelerance)
YA = compute_frf(beam_a.frequency_range, beam_a.K, beam_a.M, beam_a.C, 'accelerance');
YB = compute_frf(beam_b.frequency_range, beam_b.K, beam_b.M, beam_b.C, 'accelerance');

% Full Harmonic FRF - Experimental (accelerance)
YEA = compute_frf(beam_a.frequency_range, beam_a_exp.K, beam_a_exp.M, beam_a_exp.C, 'accelerance');
YEB = compute_frf(beam_b.frequency_range, beam_b_exp.K, beam_b_exp.M, beam_b_exp.C, 'accelerance');

% Plotting Mode Shapes from FRFs (excitation at DoF 10)
fs        = 1800;
num_modes = 3;

if ndof_equal_flag == true
    translational_dof_array = 1:2:beam_a.ndof;
end

[modes, dampr, mode_shapes] = plot_modal_waterfall( ...
    YEA, 10, beam_a.frequency_range_Hz, fs, ...
    translational_dof_array, num_modes, ...
    vec_a_exp(:, 1:num_modes), freq_a_exp(1:num_modes));

YA = YA(translational_dof_array, translational_dof_array, :);
YB = YB(translational_dof_array, translational_dof_array, :);
YEA = YEA(translational_dof_array, translational_dof_array, :);
YEB = YEB(translational_dof_array, translational_dof_array, :);

% Noise Addition to Experimental FRFs
noise_flag = true;
% pyFBS noise (paper Eq. 17): n1 / n2 are amplitude-dependent real / imag
% scales; n3 / n4 are constant real / imag noise floors.
YEA = nnoise(noise_flag, YEA, 0.005, 0.005, 1e-8, 1e-8, 12);
YEB = nnoise(noise_flag, YEB, 0.005, 0.005, 1e-8, 1e-8, 12);

plot_frf(beam_b.frequency_range_Hz, 4, 10, YB, YEB, '{Numerical}^B', '{Experimental}^B');
plot_frf(beam_a.frequency_range_Hz, 11, 17, YA, YEA, '{Numerical}^A', '{Experimental}^A');

assembly.method = 'COH'; % Correlation metric: 'COH' or 'LAC'
switch assembly.method
    case 'COH'
        [~, assembly.cor_matrix_b, ~] = func_coh(YB, YEB);
        [~, assembly.cor_matrix_a, ~] = func_coh(YA, YEA);
    case 'LAC'
        [~, assembly.cor_matrix_b, ~] = func_lac(YB, YEB);
        [~, assembly.cor_matrix_a, ~] = func_lac(YA, YEA);
end
correlation_plot(assembly.cor_matrix_b);
correlation_plot(assembly.cor_matrix_a);

%% Energy-based sensor / excitation placement (Beam A)

numsens = 4;
e_a     = analyze_frf_energy(YA, beam_a.frequency_range_Hz, numsens, 3);

best_sensor_dofs     = e_a.top_sensors;
best_excitation_dofs = e_a.top_excitations;
n_interface           = 3;
plot_beam(n_interface, beam_a.nnode, ...
         round(best_sensor_dofs), round(best_excitation_dofs), 'A');

[nat_freq_beam_a, phi_beam_a] = estimate_modal_parameters(YA, beam_a.frequency_range_Hz, ...
    'NumModes', 3, 'ProminenceThreshold', 0, 'MinPeakDistanceHz', 1);
[nat_freq_beam_a_exp, phi_beam_a_exp] = estimate_modal_parameters(YEA, beam_a.frequency_range_Hz, ...
    'NumModes', 3, 'ProminenceThreshold', 0, 'MinPeakDistanceHz', 1);

mac_matrix = calculate_mac(phi_beam_a, phi_beam_a_exp, true);

node_coords = [1:size(YA,1); zeros(1,size(YA,1))]';
visualize_mode_shapes(node_coords, phi_beam_a, nat_freq_beam_a, ...
                    'Animate', true, 'SaveVideo', false, ...
                    'VideoFile', 'beamModesAnimated.mp4');

%% Exhaustive Search

num_sensors = 2;
n_interface = 3;

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

out_a = bruteforce_search_beam(interface_dofs_a, num_sensors, YA, YEA, opts);
opts.frequency_range = beam_b.frequency_range;
out_b = bruteforce_search_beam(interface_dofs_b, num_sensors, YB, YEB, opts);

sensor_location_1     = out_a.sensor_dofs;
response_1            = out_a.sensor_dofs;
excitation_1          = out_a.excitation_dofs;
interface_cor_oval_1  = out_a.interface_cor_oval;
interface_cor_2d_1    = out_a.interface_cor_2d;
cor_1                 = out_a.cor;
ys_a                  = out_a.ys;

sensor_location_2     = out_b.sensor_dofs;
response_2            = out_b.sensor_dofs;
excitation_2          = out_b.excitation_dofs;
interface_cor_oval_2  = out_b.interface_cor_oval;
interface_cor_2d_2    = out_b.interface_cor_2d;
cor_2                 = out_b.cor;
ys_b                  = out_b.ys;


[~, cor_mat_semm_exp, cor_semm_exp] = func_coh(ys_a, YEA);
% correlation_plot(cor_mat_semm_exp, 'COH');
correlation_plot(interface_cor_2d_1, 'Coherence: ys_a and YEA');
correlation_plot(interface_cor_2d_2, 'Coherence: ys_b and YEB');

plot_frf(beam_a.frequency_range_Hz, 18, 18, YEA, ys_a, 'semm', 'Experimental');
plot_frf(beam_a.frequency_range_Hz, 2, 2, YB, ys_b, 'Numerical', 'semm', YEB, 'Experimental');

%% Plot Sensors & Excitations
close all;
plot_beam(n_interface,beam_a.nnode, sensor_location_1, excitation_1,'A')
plot_beam(n_interface,beam_b.nnode, sensor_location_2, excitation_2,'B')


%% Metaheuristic sensor / excitation search via MGO

numvars       = 2;
num_sensors_m = 1;
max_iter      = 100;
mgo_agents    = 5;

dofs            = 1:size(YA, 1);
interface_dofs  = dofs;                 % global GCCM (whole DoF range)

r_combs = nchoosek(dofs, num_sensors_m);

% NDP excitations: for each sensor combo, generate excitation combos from
% the remaining DoFs (size num_sensors+1 each).
n_outer     = size(r_combs, 1);
n_inner     = nchoosek(length(dofs) - num_sensors_m, num_sensors_m + 1);
e_combs_ndp = zeros(n_inner, num_sensors_m + 1, n_outer);
for i = 1:n_outer
    e_combs_ndp(:, :, i) = nchoosek(setdiff(dofs, r_combs(i, :)), num_sensors_m + 1);
end

% Maximise GCCM → minimise its negative.
fobj = @(x) -objective_function(x, beam_a.frequency_range, YA, YEA, ...
                                r_combs, e_combs_ndp, ...
                                false, 0, false, interface_dofs, ...
                                'COH', 3, 130.125);

lb = [1, 1];
ub = [size(r_combs, 1), size(e_combs_ndp, 1)];

[mgo_score, mgo_pos, mgo_cg] = mgo(mgo_agents, max_iter, lb, ub, numvars, fobj);

mgo_sensor_dofs     = r_combs(round(mgo_pos(1)), :);
mgo_excitation_dofs = e_combs_ndp(round(mgo_pos(2)), :, round(mgo_pos(1)));

fprintf('\nMGO best  GCCM = %.6g\n', -mgo_score);
fprintf('MGO sensors    = %s\n', mat2str(mgo_sensor_dofs));
fprintf('MGO excitations= %s\n', mat2str(mgo_excitation_dofs));

%% MGO vs brute-force one-to-one comparison

[ii, jj] = ndgrid(1:n_outer, 1:n_inner);
ii = ii(:);  jj = jj(:);
n_total  = numel(ii);

fprintf('\nBrute-forcing the same NDP search space (%d combinations)...\n', n_total);
cor_flat = -inf(n_total, 1);
parfor k = 1:n_total
    cor_flat(k) = -fobj([ii(k), jj(k)]);
end

[bf_best_gccm, idx_lin] = max(cor_flat);
bf_sensor_dofs          = r_combs(ii(idx_lin), :);
bf_excitation_dofs      = e_combs_ndp(jj(idx_lin), :, ii(idx_lin));

mgo_gccm = -mgo_score;
gap_abs  = bf_best_gccm - mgo_gccm;
gap_pct  = 100 * gap_abs / bf_best_gccm;

fprintf('\n=== MGO vs Brute-Force (identical NDP search space) ===\n');
fprintf('                 %-22s %-22s\n', 'Brute-force (global)', 'MGO');
fprintf('  sensors     :  %-22s %-22s\n', mat2str(bf_sensor_dofs),     mat2str(mgo_sensor_dofs));
fprintf('  excitations :  %-22s %-22s\n', mat2str(bf_excitation_dofs), mat2str(mgo_excitation_dofs));
fprintf('  best GCCM   :  %-22.6f %-22.6f\n', bf_best_gccm, mgo_gccm);
fprintf('  gap (BF-MGO):  %.6f  (%.3f%%)\n', gap_abs, gap_pct);
fprintf('  evaluations :  %-22d %-22d\n', n_total, mgo_agents + mgo_agents * 4 * max_iter);

%% Expanded & Full System Comparison
freqstart=1;
freqend=185;
stepsize=1;

[frequency_range,frequency_range_Hz] = frequency_generation(freqstart,freqend,stepsize);%,round(PlateA.fn),5);
translational_dof_array=1:2:length(k_full);


% 10.3.4 Damping Matrix C for Augmented Modal Damping (Fundamentals of Structural Dynamics, Kurdila & Craig, pg. 305)


% Full Harmonic FRF - Numerical
y_full = compute_frf(frequency_range, k_full, m_full, c_full, 'accelerance');
y_full = y_full(translational_dof_array,translational_dof_array,:);


% Noise Addition to Experimental FRFs

y_full = nnoise(false, y_full, 0.005, 0.005, 1e-8, 1e-8, 10);
ys_a_c = svd_truncation(ys_a, 15);
ys_b_c = svd_truncation(ys_b, 15);
y_couple = couple_substructures(ys_a_c, ys_b_c, n_interface, 'transverse');
[~, cor_mat_couple, cor_couple] = func_coh(y_couple, y_full);
correlation_plot(cor_mat_couple, 'COH');

plot_frf(frequency_range_Hz, 23, 23, y_couple, y_full, 'Coupled', 'Experimental');

