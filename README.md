# Optimal Interface Expansion

Companion code for the paper:

> **Muhammad Junaid Ali, Muhammad Safdar, Zeeshan Saeed, Abdullah Jamil, Hafiz Mohammad Mutee Ur Rehman and Muhammad Umer (2026).**
> *Enhancing structural coupling: A frequency-based methodology for optimal interface expansion.*
> Journal of Sound and Vibration 635 (2026) 119782.
> [DOI](https://doi.org/10.1016/j.jsv.2026.119782)

The methodology combines **System Equivalent Model Mixing (SEMM)** for FRF expansion with an **Optimal Sensor Placement (OSP)** framework driven by the **Global Coherence Correlation Metric (GCCM, Γ)**. Optimisation is performed by exhaustive search where feasible and by the **Mountain Gazelle Optimizer (MGO)** for larger spaces.

The repository implements the **beam case study (Case 1)** end-to-end (for now).

---

## Quick start

1. Open MATLAB R2025b (tested older versions include 2024b only).
2. From this folder, add `utils/` to the path:
   ```matlab
   addpath('utils')
   ```
3. Run the main driver:
   ```matlab
   MAIN
   ```
`MAIN.m` holds two self-contained, clearly separated examples:
- **Example 1 — Beam (coupled):** builds two steel cantilevers with Rayleigh damping, generates clean/noisy FRFs, finds optimal sensors/excitations per beam by NDP exhaustive search (2 sensors + 3 excitations), couples the SEMM expansions at the interface, and compares against the full coupled model.
- **Example 2 — Square plate (no coupling):** loads an Ansys Craig-Bampton reduced model (`Data/`), builds FRFs, finds optimal sensors/excitations over all DoFs, and plots the placement.

For **profiling / fast validation** (plots off by default, fixed RNG seeds, section timings, optional beam / square-plate validation blocks), use:
   ```matlab
   BEAM_MAIN_TEST
   ```
See `PLAN.md` for the optimisation roadmap.

---

## Repository layout

```
MAIN.m            two worked examples (beam coupling + square plate)
BEAM_MAIN_TEST.m  profiling / validation harness
utils/            FEM, SEMM, search, FRF, plotting functions
Data/             Ansys HCB reduced-model files for the square plate
Scripts/          Ansys MAPDL model-reduction script
```

## Pipeline

```
MAIN.m
    │
    ├── create_cantilever_beam     beam FEM model, eigen-solve
    ├── load_hcb_model             Ansys Craig-Bampton reduced model loader
    ├── damping                    modal | proportional (Rayleigh) damping
    ├── frequency_generation       rad/s + Hz frequency axes
    ├── compute_frf                receptance | mobility | accelerance
    ├── add_noise                  pyFBS additive noise (paper Eq. 17)
    ├── func_coh / func_lac        correlation metrics (COH, LAC)
    ├── exhaustive_search          generic DP/NDP sensor/excitation search via semm
    ├── semm                       System Equivalent Model Mixing (paper Eq. 12)
    ├── svd_truncation             rank reduction for noisy SEMM
    ├── primal_coupling            substructure assembly (K, M, C)
    ├── couple_substructures       dual LM-FBS coupling (transverse or full)
    ├── plot_square_plate          plate sensor/excitation placement
    └── mgo                        Mountain Gazelle Optimizer (metaheuristic search)
```

Plotting uses the **magma** and **cividis** colormaps from `utils/cmap_magma.m` and `utils/cmap_cividis.m`. They are self-contained and require no native toolboxes.
Note: The colormaps are different from the published article. I changed them purely out of my love for them. :)

---

## Requirements

- MATLAB R2025b
- Parallel Computing Toolbox (`parfor`, `pagemtimes`, `pagemldivide`) — recommended for more efficient computation

---

## Status

The current state of the code is limited to the beam-case only. Future additions will include the square plate case.
A supplementary MAPDL script will also be provided with a tutorial (hopefully) on how to get Craig-Bampton reduced models from Ansys Mechanical (GUI method only) 

---

## Licence

Released under the [MIT License](LICENSE) — free to use, modify, and redistribute. Citation is appreciated but not legally required (see below).

---

## Citation

If you use this code, kindly cite the paper as follows.

```@article{ALI2026119782,
title = {Enhancing structural coupling: A frequency-based methodology for optimal interface expansion},
journal = {Journal of Sound and Vibration},
volume = {635},
pages = {119782},
year = {2026},
issn = {0022-460X},
doi = {https://doi.org/10.1016/j.jsv.2026.119782},
url = {https://www.sciencedirect.com/science/article/pii/S0022460X26001458},
author = {Muhammad Junaid Ali and Muhammad Safdar and Zeeshan Saeed and Abdullah Jamil and Hafiz Mohammad {Mutee Ur Rehman} and Muhammad Umer},
keywords = {Dynamic substructuring, Optimal sensor placement, Interface dynamics, Frequency response functions, System equivalent model mixing},
```
