# PLAN — Beam SEMM / OSP codebase: status & roadmap

Companion code for:

> *Enhancing structural coupling: A frequency-based methodology for optimal interface expansion.*
> JSV 635 (2026) 119782.

**Guiding principle:** speed up and generalise the implementation while preserving the
physics and numerical results (no shortcut may change FRF quality, frequency resolution,
or the model definition). Every optimisation must be proven bit-equivalent (or
tolerance-equivalent) to the published workflow.

**Drivers:**

| Script | Role |
|--------|------|
| `MAIN.m`           | Two worked examples — beam coupling + square plate. Uses `exhaustive_search`, `load_hcb_model`, accelerated `mgo`, `add_noise`. |
| `BEAM_MAIN_TEST.m` | Profiling / validation harness (`do_plots`, fixed seeds, section timings). |

**Status key:** ✅ done · 🟡 in progress · 🔬 needs proof before implementation · ⏸ deferred · ❌ out of scope

---

## 1. Achieved

### Performance & reproducibility (first sprint)

| Item | Where | Notes |
|------|-------|-------|
| ✅ `do_plots` gate | `BEAM_MAIN_TEST.m` | All non-essential figures off by default; console + timings always run. |
| ✅ Reuse `eye(n)` in FRF solve | `utils/compute_frf.m` | `z_j \ I` with one preallocated identity inside `parfor`. |
| ✅ RNG seeds | `BEAM_MAIN_TEST.m` | `add_noise(...,12)`; `rng(42/12,'philox')` before `mgo`. Deterministic reruns. |
| ✅ `svd_truncation` semantics | `utils/svd_truncation.m` | 2nd arg `num_drop` = trailing singular values discarded (pyFBS-style). |
| ✅ Interface-DoF contract | `utils/validate_interface_dofs.m` | Shared checks called from `primal_coupling`, `couple_substructures`, test harness. |
| ✅ No plotting inside search loops | search + objective utils | Figures gated in the harness only. |

### MGO accelerated baseline (this sprint) — ✅ done

`mgo_v2` has been **promoted to the baseline `utils/mgo.m`** and retired:

- `utils/mgo.m` now does per-iteration **dedup** + cross-iteration **cache** (deterministic
  objectives), with optional `parfor` (`opts.use_parfor`, default false). RNG sequence is
  unchanged, so optima match the original MGO exactly.
- New 4th output `stats` (`n_fobj_eval`, `n_cache_hit`, `n_dedup_saved`). `MAIN.m`
  and `BEAM_MAIN_TEST.m` now report the **real** MGO evaluation count from `stats` instead
  of the theoretical `n + n*4*max_iter`.
- Retired: `utils/mgo_v2.m`, `utils/mgo_benchmark_compare.m`, `utils/expensive_rastrigin.m`,
  and the `benchmark_mgo` block in `BEAM_MAIN_TEST.m` (comparison complete — results below).

**Validated benchmark result (before retirement)** — dedup+cache reproduces the vanilla
optimum on both objectives and cuts work substantially on the expensive one:

| Objective | Variant | Wall | fobj evals | Speedup |
|-----------|---------|------|-----------|---------|
| Beam SEMM (GCCM) | vanilla | 7.53 s | 1005 | — |
| | dedup+cache (serial) | 4.22 s | 555 | 1.8× |
| | dedup+cache (parfor) | 3.06 s | 555 | 2.5× |
| Generic Rastrigin | vanilla | 0.136 s | 1005 | — |
| | dedup+cache (serial) | 0.081 s | 294 | 1.7× |

Same optimum everywhere (beam `[5 64]`, GCCM 0.990875; Rastrigin `[1 1]`, f=2).
`parfor` only helps when `fobj` is expensive (on cheap objectives, pool overhead dominates) →
keep `use_parfor=false` for cheap objectives.

### Noise function consolidation — ✅ done

`utils/nnoise.m` → **`utils/add_noise.m`**. Confirmed mathematically identical to the
experimental `addNoise.m` (same four-component model
`η = n1|Y|G1 + i·n2|Y|G2 + n3 G3 + i·n4 G4`). Our version keeps the `status` bypass and
`philox` seeding for reproducibility. All call sites updated (`MAIN.m`,
`BEAM_MAIN_TEST.m`, README pipeline).

### Baseline timing (live reference)

`BEAM_MAIN_TEST` (`do_plots=0`, pool on): total **~20 s**, dominated by the
**MGO + NDP verification `parfor` (~7–8 s)**. FRF generation is cheap (~0.12 s for four
700-point tensors); brute-force DP (136×2) ~0.5 s; coupled validation ~0.1 s. The dominant
cost is the optimisation/search phase — that is where the remaining roadmap is aimed.

Determinism check (two identical runs): sensors `5`, excitations `[4 18]`, GCCM
**0.990875**, gap 0%; final coupled Γ **0.916420**.

Re-profile: set `profile_run = true` at the top of `BEAM_MAIN_TEST.m`, run once, `profile viewer`.

---

## 2. Next goals

### 2.1 Generalised `exhaustive_search.m` — ✅ built, 🟡 validation pending  ← headline item

**Done.** `utils/exhaustive_search.m` replaces the three near-duplicate per-problem searches
(`bruteforce_search_beam.m`, experimental `bruteforceSearchSquare.m`, `bruteforceSearch.m`).
Studying all three showed they differ in only four parametrisable concepts; the third
(`bruteforceSearch.m`, coupled plates) added the key insight that **candidate and validation DoF
sets are fully independent** (its `Plate_name`/`node_config` switch just computes those two sets):

| Concept | Beam | Square | Coupling P1/P2 | Generalised opts field |
|---------|------|--------|----------------|------------------------|
| Candidate pool | `setdiff(all, interface)` | all | `1:n-17` / `18:n` | `candidate_dofs` |
| Scored subset | interface | all | last/first 17 | `validation_dofs` |
| Metrics | `'COH'` | cellstr | cellstr | `methods` (char/cellstr) |
| Extra excitations | `ec` | `c` | fixed +1 | `extra_excitations` |

Geometry (beam line vs plate grid, Z-axis extraction) **stays in the caller** — the function only
sees a clean `n×n×nFreq` tensor. `MainSquare.m` calls it with `candidate = validation = 1:n`;
the coupled plates with disjoint candidate/validation sets.

**I/O contract** — `out = exhaustive_search(y_n, y_e, opts)`; `out` is a struct array (one per
method) with `.method .sensor_dofs .excitation_dofs .cor_overall .cor_2d .cor_landscape .ys`.
See the function header for the full opts list.

**SEMM engine — `semm.m` only (single source of truth).** The expansion is delegated entirely to
`semm.m`; no SEMM math is duplicated in `exhaustive_search`. `semm.m` is left untouched. Kept the
performance wins that do **not** reimplement SEMM:

1. **SEMM computed once per combination, scored on all metrics** — the experimental versions
   recompute SEMM per method; this does not.
2. **`parfor` over the outer (sensor) loop**; degrades gracefully to serial (`use_parfor=false` /
   no pool).
3. Full `ys` is stored only for the per-method winner (recomputed once at the end), not for every
   combination — keeps memory flat.

(The interface-only factorisation / `proj_L` hoisting was dropped to avoid duplicating `semm.m`'s
math — a deliberate trust-over-speed choice.)

**Testing — one structure at a time, with real MATLAB runs.**

*Beam — automated in `BEAM_MAIN_TEST.m` (set `test_exhaustive = true`).* Validates against the
references already computed in the same run, so no extra brute force is needed:

| # | Check | Reference | Tol |
|---|-------|-----------|-----|
| A | DP Beam A | `bruteforce_search_beam` (`out_a`) | 1e-6 |
| B | DP Beam B | `bruteforce_search_beam` (`out_b`) | 1e-6 |
| C | NDP global (1 sensor, all DoFs) | MGO-section brute force (`bf_*`, 3420 combos) | 1e-6 |
| D | multi-method `{'COH','LAC'}` | returns 2 structs | — |

*Square plate — automated in `BEAM_MAIN_TEST.m` (`test_square_plate = true`).* Validates both new
functions against `Square_Exhaustive_Result.mat` (a `MainSquare.m` workspace dump, NDP `c=0`,
sensors `[20 23]`, excitations `[22 24]`, GCCM `0.990348`):
- **A** — `exhaustive_search` on the *saved* `PlateA.Yz/YEz` (so noise realisation matches the
  reference). **Verified PASS** (full NDP, 1749 s / 12 threads): sensors `[20 23]`, excitations
  `[22 24]`, GCCM `0.990348` — exact match to `bruteforceSearchSquare`. Flag `test_square_plate`
  defaults **off** (the NDP search is ~29 min); flip on for full validation.
- **B** — `load_hcb_model` rebuilt from the Ansys HCB files vs the experimental `sortKM` output.
  **Verified bit-exact: relK = relM = 0**, NCMSMode 200, n 36, measured-DOF set matches. Self-skips
  when the Ansys source files are absent.

*L-plate (coupled plates) — pending.* Same pattern from `Coupling_YA_YB.m`: `P1` →
`candidate=1:n-17`, `validation=(n-16):n`; `P2` → `candidate=18:n`, `validation=1:17`.

**Primary pass criterion = selected sensor/excitation DoFs** (the argmax), robust to small numeric
differences between our `semm.m`/`func_coh.m` and the experimental `SEMM`/`correlation`.

**Correlation metrics:** only `COH` and `LAC` are required/supported (`func_coh.m`, `func_lac.m`).
`FRAC/MFRAC/FAAC` from the experimental scripts are out of scope.

Then: retire `bruteforce_search_beam.m` (make it a shim or migrate `MAIN`).

### 2.2 Square plate case (Case 2) integration — 🟡 loader done

- ✅ **`utils/load_hcb_model.m`** — single self-contained loader for an Ansys Craig-Bampton (HCB)
  reduced model: inlines a compact assembled-symmetric (RSA) Harwell-Boeing reader, the `.mapping`
  parse, and the physical-DOF reorder (CMS modes last); reuses `frequency_generation`, `damping`
  (≡ experimental `modaldamping`), `compute_frf` (≡ `FEMgrecept`), and `add_noise`. Returns a
  `model` struct with reordered `K/M`, `NodeID`, `NodeCoords`, `physical_DoF_array`,
  `measured_dof_array`, `frequency_range(_Hz)`, `fn`, `C`, full `Y/YE/YE_noise`, and measured-DOF
  `Yz/YEz` (square plate: 36×36). Drop-in for `exhaustive_search` with `candidate=validation=1:n`.
- Remaining: optional `plot_square_plate` (geometry plotting / `AggregationIndex` equivalent).
  Dynamics path otherwise reuses existing utils — no `SEMM`/`correlation` re-port needed.

Caveats (documented in the function): `add_noise` uses philox vs the experimental Mersenne stream
(same model, different draw); `compute_frf` omits `FEMgrecept`'s negligible `eps*eye` term.

### 2.3 Modal-superposition FRF path — ⏸ deferred

Optional faster `compute_frf` via truncated modal superposition; only worth it if the frequency
band / mesh grows. Direct dynamic-stiffness inversion is already cheap at the current size.

### 2.4 MAPDL Craig-Bampton tutorial — ⏸ future

Supplementary script + tutorial for getting CB-reduced models out of Ansys Mechanical (per README).

---

## 3. Excluded / rejected

- Reduce `K,M,C` before `compute_frf` (slice `Y` after instead).
- `decomposition`/`linsolve` per frequency (current `z\I` sufficient).
- Coarser frequency grid during search (keep one grid for search + validation).
- Replacing `pagepinv` wholesale.

---

## 4. Reference

- Profiler (pre-sprint): `D:\NUST Data\MATLAB_CodeBase\Profiler_Analysis.pdf`
- Profiler (post-sprint): `D:\NUST Data\MATLAB_CodeBase\Profiler_Analysis_Test.pdf`
- Square plate workflow: `D:\NUST Data\MATLAB_CodeBase\Experimental Final Model Code\MainSquare.m`
- pyFBS: [doi:10.21105/joss.03399](https://doi.org/10.21105/joss.03399)
