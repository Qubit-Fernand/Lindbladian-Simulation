# Lindbladian Simulation

Numerical experiments for dissipative transverse-field Ising model (TFIM) dynamics based on Lindblad evolution and Suzuki-Trotter product formulas.

## Overview

This repository studies open-system spin-chain dynamics with:

- coherent TFIM evolution split into `H_x` and `H_z`,
- local Lindblad jump operators on each site,
- Trotterized propagation of density matrices and superoperators,
- comparisons against exact Liouvillian evolution.

The code is built around QuTiP objects and stores generated density matrices as `.npy` files for later analysis in notebooks and plots.

## Main Files

- `Dissipative_TFIM_rho.py`: evolve a density matrix directly and save `rho_N_*_r_*.npy`.
- `Dissipative_TFIM_diamond.py`: build the corresponding superoperator and compare against the exact channel with diamond norm.
- `extrapolation.py`: generate extrapolated datasets, including exact superoperator evolution and higher-`r` runs without Euler approximation for the dissipative part.
- `exp(tL)_Euler.py`: Euler-based evolution utilities for `exp(tL)` experiments.
- `Dissipative_TFIM.ipynb`, `precision_compare.ipynb`: notebook-based exploration and analysis.
- `tests/test_simulation.py`: physics- and numerics-oriented regression tests.

## Repository Layout

- `data/`: generated datasets for different initial states and dissipation strengths.
- `jobs/`: job-oriented scripts and intermediate outputs.
- `plot/`: plotting helpers and generated figure-related assets.
- `archive/`: older experiments and reference implementations kept for comparison.
- `data_sigmam/`: archived datasets for the `sigma_-` convention.

The directory names encode initial states such as `|000>initial`, `|111>initial`, `|+++>initial`, `|010>initial`, and `|mixed>initial`, with subdirectories like `gamma_0.1` and `gamma_1.0`.

## Environment

Python 3.10+ is recommended.

Typical dependencies used in this repository:

- `numpy`
- `qutip`
- `cvxpy`
- `pytest`
- `tqdm`

Example installation:

```bash
python -m venv .venv
source .venv/bin/activate
pip install numpy qutip cvxpy pytest tqdm
```

If you plan to compute diamond norms with a specific backend, make sure the corresponding solver is installed and configured for `cvxpy` / QuTiP.

## Quick Start

Generate density-matrix data for a chosen initial state:

```bash
python Dissipative_TFIM_rho.py --N 4 --r 4 --t 0.2 --initial 0
```

Generate extrapolated data:

```bash
python extrapolation.py --N 5 --t 0.2 --initial 1
```

Compute a channel-distance benchmark:

```bash
python Dissipative_TFIM_diamond.py --N 4 --gamma 0.1 --r 2
```

## Testing

Run the regression tests from the repository root:

```bash
pytest tests/test_simulation.py -v
```

The test suite checks basic operator structure, trace preservation, positivity, physical limits, and Trotter error convergence.

## Notes

- Several scripts loop over `gamma in [0.1, 1.0]` and write outputs into state-specific folders.
- The repository contains large precomputed datasets; generated `.npy` outputs under `jobs/` are ignored by git.
- Some archived files document older conventions or known issues and are preserved for research traceability rather than polished package design.
