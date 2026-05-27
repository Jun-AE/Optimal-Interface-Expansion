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
   beam_coh_main
   ```
The script builds two cantilever beams, generates clean and noisy FRFs, runs SEMM expansion, performs exhaustive and MGO-based sensor/excitation search, and compares the two.

---

## Pipeline

```
beam_coh_main.m
    │
    ├── create_cantilever_beam     FEM model, eigen-solve
    ├── damping                    modal | proportional damping
    ├── primal_coupling            substructure assembly (K, M, C)
    ├── frequency_generation       rad/s + Hz frequency axes
    ├── compute_frf                receptance | mobility | accelerance
    ├── nnoise                     pyFBS additive noise (paper Eq. 17)
    ├── plot_modal_waterfall       imag(FRF) waterfall + mode overlays
    ├── func_coh / func_lac        correlation metrics
    ├── correlation_plot           annotated heatmap
    ├── analyze_frf_energy         energy-based DoF selection
    ├── estimate_modal_parameters  CMIF + peak picking
    ├── calculate_mac              MAC matrix with annotated heatmap
    ├── visualize_mode_shapes      animated mode shapes (optional video)
    ├── bruteforce_search_beam     exhaustive sensor/excitation search via semm
    ├── semm                       System Equivalent Model Mixing (paper Eq. 12)
    ├── couple_substructures       dual LM-FBS coupling (transverse or full)
    ├── svd_truncation             rank reduction for noisy SEMM
    ├── mgo                        Mountain Gazelle Optimizer
    └── objective_function         GCCM fitness for MGO
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
