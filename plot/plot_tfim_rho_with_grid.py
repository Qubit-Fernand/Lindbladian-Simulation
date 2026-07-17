from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({"font.size": 18})
plt.rcParams["text.usetex"] = True
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = PROJECT_ROOT / "data"
OUT_DIR = PROJECT_ROOT / "plot"

N_LIST = np.arange(4, 11)
R_LIST = np.array([1, 2, 3, 7])
GAMMAS = [0.1, 1.0]
STEP = 100000

INITIALS = ["111", "000", "+++", "mixed"]
LEGEND_TITLES = {
    "000": r"$\rho_0=|0^N\rangle\langle 0^N|$",
    "111": r"$\rho_0=|1^N\rangle\langle 1^N|$",
    "+++": r"$\rho_0=|+^N\rangle\langle +^N|$",
    "mixed": r"$\rho_0=I^{\otimes N}/2^N$",
}


def trace_norm_hermitian(mat):
    mat = 0.5 * (mat + mat.conj().T)
    return float(np.sum(np.abs(np.linalg.eigvalsh(mat))))


def load_rho(initial, gamma, filename):
    path = DATA_ROOT / f"|{initial}>initial" / f"gamma_{gamma}" / filename
    if not path.exists():
        raise FileNotFoundError(path)
    return np.load(path)


def make_plot(initial):
    fig, ax = plt.subplots(figsize=(12.5, 6))
    handles = []
    log_n = np.log10(N_LIST)

    for gamma in GAMMAS:
        marker = "o" if gamma == 0.1 else "d"
        linestyle = "-" if gamma == 0.1 else "--"
        for color, r in enumerate(R_LIST):
            errors = []
            for N in N_LIST:
                rho_trotter = load_rho(initial, gamma, f"rho_N_{N}_r_{r}.npy")
                rho_exact = load_rho(initial, gamma, f"rho_exact_N_{N}_step_{STEP}.npy")
                errors.append(trace_norm_hermitian(rho_trotter - rho_exact))

            log_e = np.log10(errors)
            slope, intercept = np.polyfit(log_n, log_e, 1)
            x_fit = np.linspace(log_n.min(), log_n.max(), 100)

            ax.plot(log_n, log_e, marker, color=f"C{color}", markersize=4)
            ax.plot(x_fit, slope * x_fit + intercept, linestyle=linestyle, color=f"C{color}", linewidth=1.3)
            handles.append(
                Line2D(
                    [0],
                    [0],
                    color=f"C{color}",
                    linestyle=linestyle,
                    marker=marker,
                    markersize=5,
                    label=rf"$\gamma={gamma:.1f}$, $r={r}$, slope={slope:.3f}",
                )
            )

    ax.legend(
        title=LEGEND_TITLES[initial],
        title_fontsize=20,
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=20,
    )
    ax.set_xticks(log_n)
    ax.set_xticklabels([str(n) for n in N_LIST])
    if initial != "mixed":
        ax.set_yticks([-3.5, -3, -2.5, -2, -1.5])
    ax.set_xlabel("qubit number $N$", fontsize=24)
    ax.set_ylabel(r"$\lg \|(e^{t\mathcal{L}}-\mathcal{S}(t/r)^r)\rho_0\|_1$", fontsize=21)
    ax.set_title(r"Lindbladian simulation error vs. system size $N$", fontsize=24)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()

    out_path = OUT_DIR / f"error-N-{initial}-sigmap.pdf"
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Wrote {out_path.relative_to(PROJECT_ROOT)}")


def main():
    for initial in INITIALS:
        make_plot(initial)


if __name__ == "__main__":
    main()
