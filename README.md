# Lindbladian Simulation

Numerical experiments for Lindbladian product-formula simulation and Richardson extrapolation. The main benchmark is a dissipative transverse-field Ising model (TFIM); the current robustness checks add 1D Heisenberg chains with local amplitude damping and local dephasing.

## Overview

This repository stores reproducible data-generation scripts, raw density matrices, and plotting utilities for the PRX Quantum manuscript workflow. The main experiments compare Trotterized density-matrix evolution against exact Liouvillian evolution and then fit the resulting error scalings.

The numerical workflow intentionally separates data generation from plotting:

- scripts write `.npy` density matrices and `metadata.json` provenance under `data/`;
- plotting notebooks/scripts read those data and write PDFs under `plot/`;
- manuscript-ready figures are copied into `overleaf/figs/` only after inspection.

## Main Models

### Dissipative TFIM

The TFIM Hamiltonian is split into two coherent groups,

```text
H = H_X + H_Z,
H_X = -J sum_j X_j X_{j+1},
H_Z = -h sum_j Z_j.
```

The current manuscript convention uses local amplitude damping

```text
L_nu = sqrt(gamma) |0><1|_nu.
```

The main channel-error figures use `J=1.0`, `h=0.5`, `t=0.2`, `gamma in {0.1, 1.0}`, `N=4..10`, and the initial states `|1^N><1^N|`, `|0^N><0^N|`, `|+^N><+^N|`, and `I^{\otimes N}/2^N`.

### Heisenberg Robustness Benchmarks

The Heisenberg Hamiltonian is

```text
H = J sum_j (X_j X_{j+1} + Y_j Y_{j+1} + Z_j Z_{j+1}) + h sum_j X_j.
```

Two local noise mechanisms are currently generated:

- `amplitude_damping`: `L_nu = sqrt(gamma) |0><1|_nu`;
- `dephasing`: `L_nu = sqrt(gamma) Z_nu`.

The channel-error benchmark parameters are `J=1.0`, `h=0.5`, `t=0.15`, `gamma in {0.1, 1.0}`, `r in {4, 6, 8, 10}`, `N=4..10`, and initial states `0`, `1`, `+`, `I`. The corresponding data directories are:

- `data/noisy_heisenberg_amplitude_damping/`
- `data/noisy_heisenberg_dephasing/`

An exploratory Heisenberg amplitude-damping Richardson-extrapolation dataset is also available under `data/noisy_heisenberg_amplitude_damping_extrapolated/`. It uses `N=5`, `J=1.0`, `h=0.5`, `t=0.15`, `gamma in {0.1, 1.0}`, and initial states `|1^N><1^N|` and `|0^N><0^N|`. These figures are useful for diagnosis, but are not currently included in the manuscript.

## Main Files

- `Dissipative_TFIM_rho.py`: generate TFIM density-matrix channel-error data.
- `Dissipative_TFIM_diamond.py`: compare TFIM channels using the diamond norm.
- `extrapolation.py`: generate TFIM Richardson-extrapolation density matrices.
- `Noisy_Heisenberg_rho.py`: generate Heisenberg channel-error data for amplitude damping or dephasing.
- `Noisy_Heisenberg_extrapolation.py`: generate Heisenberg amplitude-damping extrapolation density matrices.
- `exp(tL)_Euler.py`: Euler-based utilities for `exp(tL)` experiments.
- `plot/plot_tfim_rho_with_grid.py`: regenerate TFIM channel-error PDFs with grid styling.
- `plot/plot_noisy_heisenberg_rho.ipynb`: plot the Heisenberg channel-error benchmark.
- `plot/plot_extrapolation_observable_scaling.py`: plot TFIM Richardson extrapolation for magnetization and nearest-neighbor `ZZ`.
- `plot/plot_heisenberg_extrapolation_observable_scaling.py`: plot exploratory Heisenberg extrapolation figures.
- `tests/test_simulation.py`: physics- and numerics-oriented regression tests.
- `codex_memory.md`: project-specific implementation notes and pitfalls.

## Repository Layout

- `data/`: generated density matrices and metadata.
- `jobs/`: job scripts and runtime logs.
- `plot/`: plotting notebooks, plotting scripts, PDFs, and summary JSON files.
- `overleaf/`: local Overleaf mirror for the PRX manuscript; ignored by this git repository.
- `archive/`: older experiments and reference implementations kept for comparison.
- `data_sigmam/`: archived datasets for an older operator convention.

Directory names encode initial states, for example `|000>initial`, `|111>initial`, `|+++>initial`, `|mixed>initial`, and their extrapolated analogues. Coupling strengths live in subdirectories such as `gamma_0.1` and `gamma_1.0`.

## Environment

Use the local conda environment named `Qubit` unless there is a specific reason not to:

```bash
conda run -n Qubit python --version
```

Typical dependencies are:

- `numpy`
- `scipy`
- `qutip`
- `cvxpy`
- `matplotlib`
- `pytest`
- `tqdm`

If you compute diamond norms with a specific backend, make sure the corresponding solver is installed and configured for `cvxpy` / QuTiP.

## Common Commands

Generate TFIM density-matrix data for one initial state:

```bash
conda run -n Qubit python Dissipative_TFIM_rho.py --N 4 --r 4 --t 0.2 --initial 0
```

Generate TFIM extrapolated data:

```bash
conda run -n Qubit python extrapolation.py --N 5 --t 0.2 --initial 1
```

Generate Heisenberg channel-error data:

```bash
conda run -n Qubit python Noisy_Heisenberg_rho.py --noise amplitude_damping
conda run -n Qubit python Noisy_Heisenberg_rho.py --noise dephasing
```

Generate exploratory Heisenberg extrapolation data:

```bash
conda run -n Qubit python Noisy_Heisenberg_extrapolation.py --noise amplitude_damping --initials 1 0
```

Regenerate script-based extrapolation figures:

```bash
conda run -n Qubit python plot/plot_extrapolation_observable_scaling.py
conda run -n Qubit python plot/plot_heisenberg_extrapolation_observable_scaling.py
```

Compile the PRX manuscript from the Overleaf mirror:

```bash
cd overleaf
latexcompile prx-q.tex
```

Upload selected Overleaf files after inspection:

```bash
cd overleaf
olcli upload prx-q.tex
olcli upload figs/<figure-name>.pdf
```

## Data and Plot Provenance

Heisenberg data scripts write `metadata.json` files containing the parameter values, vectorization convention, and sanity warnings. Current metadata files are:

- `data/noisy_heisenberg_amplitude_damping/metadata.json`
- `data/noisy_heisenberg_dephasing/metadata.json`
- `data/noisy_heisenberg_amplitude_damping_extrapolated/metadata.json`

The Heisenberg implementation uses column-major vectorization, `order="F"`, and the convention

```text
vec(A rho B) = (B^T \otimes A) vec(rho).
```

The local dissipative channels in `Noisy_Heisenberg_rho.py` are applied exactly with closed-form maps rather than Euler substeps or explicit local-superoperator exponentials.

## Testing

Run the regression tests from the repository root:

```bash
conda run -n Qubit pytest tests/test_simulation.py -v
```

The test suite checks basic operator structure, trace preservation, positivity, physical limits, and Trotter error convergence.

## Notes and Pitfalls

- Be careful with `sigmam` / `sigmap` naming. In this project the intended amplitude-damping jump is explicitly `|0><1|`, represented in `Noisy_Heisenberg_rho.py` by `np.array([[0, 1], [0, 0]])`.
- The original TFIM script applies local dissipators with many Euler substeps, so it can be much slower than the newer Heisenberg code, where local noise channels use exact closed forms.
- `overleaf/` is ignored by the project git repository. Git-based revert commands will not restore `overleaf/prx-q.tex`; use manual patches, local backups, or Overleaf history.
- Generated `.npy` datasets can be large. Prefer reading metadata and plotting summaries before re-running long jobs.
