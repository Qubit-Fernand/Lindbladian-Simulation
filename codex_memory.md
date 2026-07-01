# Codex Memory: Lindbladian Simulation

## Project Scope

This folder is the local numerical-experiment project for Lindbladian simulation. It is separate from the Overleaf manuscript folder. Put model/data-generation scripts at the project root, generated density-matrix data under `data/`, and plotting notebooks/figures under `plot/`.

## Current Heisenberg Robustness Benchmark

- Main data-generation script: `Noisy_Heisenberg_rho.py`.
- The script is intentionally data-only: it should save `.npy` density matrices and `metadata.json`, not generate figures directly.
- Output directories:
  - `data/noisy_heisenberg_amplitude_damping/`
  - `data/noisy_heisenberg_dephasing/`
- Each noise directory follows the existing TFIM state/gamma convention:
  - `|000>initial/gamma_0.1/`
  - `|111>initial/gamma_0.1/`
  - `|+++>initial/gamma_0.1/`
  - `|mixed>initial/gamma_0.1/`
  - and similarly for `gamma_1.0`.
- Saved files follow the existing rho naming style:
  - `rho_superexact_N_<N>.npy`
  - `rho_N_<N>_r_<r>.npy`
- Default Heisenberg parameters currently used:
  - `J=1.0`
  - `h=0.5`
  - `t=0.15`
  - `gamma in {0.1, 1.0}`
  - `r in {4, 6, 8, 10}`
  - initial states `0`, `1`, `+`, `I`

## Plotting

- Plotting is handled by `plot/plot_noisy_heisenberg_rho.ipynb`.
- Generated figures go under `plot/noisy_heisenberg/`.
- Follow the existing TFIM plotting notebook style:
  - Times New Roman
  - `plt.rcParams["text.usetex"] = True`
  - figure size around `(12.5, 6)` for a single initial-state plot
  - plot `log10(N)` on the x-axis but manually set tick labels to `4, 5, 6, ...`, not scientific notation
  - put the legend on the right with `bbox_to_anchor=(1.02, 0.5)`
- Do not make a combined four-panel figure unless explicitly requested; the current preference is one PDF per initial state.

## Operator Convention Pitfall

Be very careful with `sigmam` / `sigmap` naming. The confusing point in this project is that the symbol name can be opposite to the intended physical direction depending on basis convention and library notation.

For amplitude damping in the current Heisenberg code, do not rely on QuTiP naming. Use the explicit matrix

```python
AMP_DAMP_JUMP = np.array([[0, 1], [0, 0]], dtype=complex)
```

This is `|0><1|` in the column-vector computational basis:

- `AMP_DAMP_JUMP @ |0> = 0`
- `AMP_DAMP_JUMP @ |1> = |0>`

So the current `Noisy_Heisenberg_rho.py` amplitude-damping channel maps population from `|1>` to `|0>`, matching the intended `L = sqrt(gamma) |0><1|` convention. If changing the noise model or porting back to QuTiP, re-check the actual matrix action instead of trusting `sigmam` / `sigmap` names.

## Remote Mini Notes

The Mac mini has this folder synced through iCloud. As of the current run, the new Heisenberg data and plots had synced successfully. The mini did not have a `Qubit` conda environment; for `Noisy_Heisenberg_rho.py`, `/opt/anaconda3/bin/python3` was sufficient because the script only needs NumPy/SciPy/tqdm.

Long N=9 and N=10 Heisenberg data generation was started on the mini in tmux session `noisy_heis_N9_N10`, with log:

```text
jobs/noisy_heisenberg_N9_N10.log
```

Check with:

```bash
ssh mini 'tmux attach -t noisy_heis_N9_N10'
ssh mini 'tail -f "/Users/AntiEntropy/Documents/Research/Lindbladian Simulation/jobs/noisy_heisenberg_N9_N10.log"'
```
