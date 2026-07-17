import json
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
DATA_ROOT = PROJECT_ROOT / "data" / "noisy_heisenberg_amplitude_damping_extrapolated"
OUT_DIR = PROJECT_ROOT / "plot"

N = 5
GAMMAS = [0.1, 1.0]
SCALE_LIST = list(range(1, 10))
R_LIST = np.array([1, 2, 4])
RICHARDSON_COEFFS = np.array([1 / 45, -4 / 9, 64 / 45])
ERROR_FLOOR = 1e-18

INITIAL_DIRS = {
    "111": "|111>extrapolated",
    "000": "|000>extrapolated",
}

INITIAL_TITLES = {
    "111": r"\rho_0=|1^N\rangle\langle 1^N|",
    "000": r"\rho_0=|0^N\rangle\langle 0^N|",
}

GAMMA_LINESTYLES = {0.1: "-", 1.0: "--"}
GAMMA_MARKERS = {0.1: "o", 1.0: "d"}


def kron_all(ops):
    out = np.array([[1]], dtype=complex)
    for op in ops:
        out = np.kron(out, op)
    return out


def local_operator(site_ops):
    identity = np.eye(2, dtype=complex)
    return kron_all([site_ops.get(site, identity) for site in range(N)])


def observable_matrix(observable):
    z = np.array([[1, 0], [0, -1]], dtype=complex)
    if observable == "magnetization_z":
        return sum(local_operator({site: z}) for site in range(N))
    if observable == "nearest_neighbor_zz":
        return sum(local_operator({site: z, site + 1: z}) for site in range(N - 1))
    raise ValueError(f"Unknown observable: {observable}")


def observable_label(observable):
    if observable == "magnetization_z":
        return r"O=\sum\nolimits_{j=1}^{N} Z_j"
    if observable == "nearest_neighbor_zz":
        return r"O=\sum\nolimits_{j=1}^{N-1} Z_jZ_{j+1}"
    raise ValueError(f"Unknown observable: {observable}")


def legend_title(initial, observable):
    return rf"$N={N},\ {INITIAL_TITLES[initial]}$" + "\n" + rf"${observable_label(observable)}$"


def load_density(initial, gamma, filename):
    path = DATA_ROOT / INITIAL_DIRS[initial] / f"gamma_{gamma}" / filename
    if not path.exists():
        raise FileNotFoundError(path)
    return np.load(path)


def collect_plot_data(initial, observable):
    O = observable_matrix(observable)
    data = {}
    for gamma in GAMMAS:
        rho_exact = load_density(initial, gamma, f"rho_superexact_N_{N}_no_Euler.npy")
        exact_value = float(np.real(np.trace(O @ rho_exact)))

        scales = []
        extrapolated = []
        max_r = []
        for scale in SCALE_LIST:
            values = []
            missing = False
            for r_base in R_LIST:
                filename = f"rho_N_{N}_r_{scale * int(r_base)}_no_Euler.npy"
                path = DATA_ROOT / INITIAL_DIRS[initial] / f"gamma_{gamma}" / filename
                if not path.exists():
                    missing = True
                    break
                rho = np.load(path)
                values.append(float(np.real(np.trace(O @ rho))))
            if missing:
                continue

            scales.append(scale)
            max_r.append(values[-1])
            extrapolated.append(float(np.dot(RICHARDSON_COEFFS, values)))

        data[gamma] = {
            "scales": np.asarray(scales, dtype=float),
            "extrapolated": np.asarray(extrapolated, dtype=float),
            "max_r": np.asarray(max_r, dtype=float),
            "exact": exact_value,
        }
    return data


def safe_log10_error(values, exact):
    return np.log10(np.maximum(np.abs(values - exact), ERROR_FLOOR))


def make_plot(initial, observable, out_name):
    plot_data = collect_plot_data(initial, observable)
    fig, ax = plt.subplots(figsize=(12.5, 6))
    handles = []
    summary = {}

    for gamma in GAMMAS:
        scales = plot_data[gamma]["scales"]
        exact = plot_data[gamma]["exact"]
        extrapolated = plot_data[gamma]["extrapolated"]
        max_r = plot_data[gamma]["max_r"]
        marker = GAMMA_MARKERS[gamma]
        linestyle = GAMMA_LINESTYLES[gamma]

        x = np.log10(scales)
        x_fit = np.linspace(x.min(), x.max(), 200)

        y_extrapolated = safe_log10_error(extrapolated, exact)
        y_max_r = safe_log10_error(max_r, exact)

        k_extrapolated, b_extrapolated = np.polyfit(x, y_extrapolated, 1)
        k_max_r, b_max_r = np.polyfit(x, y_max_r, 1)

        ax.plot(x, y_extrapolated, color="C0", marker=marker, linestyle="None")
        ax.plot(x_fit, k_extrapolated * x_fit + b_extrapolated, color="C0", linestyle=linestyle)
        handles.append(
            Line2D(
                [0],
                [0],
                color="C0",
                linestyle=linestyle,
                marker=marker,
                markersize=6,
                label=rf"$\gamma={gamma:.1f}$, extrap., slope=${k_extrapolated:.3f}$",
            )
        )

        ax.plot(x, y_max_r, color="C1", marker=marker, linestyle="None")
        ax.plot(x_fit, k_max_r * x_fit + b_max_r, color="C1", linestyle=linestyle)
        handles.append(
            Line2D(
                [0],
                [0],
                color="C1",
                linestyle=linestyle,
                marker=marker,
                markersize=6,
                label=rf"$\gamma={gamma:.1f}$, original, slope=${k_max_r:.3f}$",
            )
        )

        summary[f"gamma={gamma:.1f}"] = {
            "available_scales": scales.astype(int).tolist(),
            "extrapolated_slope": float(k_extrapolated),
            "original_slope": float(k_max_r),
            "exact_observable": exact,
        }

    ax.set_xlabel(r"$r_{\mathrm{scale}}$", fontsize=22)
    ax.set_ylabel(
        r"$\lg|\mathrm{Tr}[O\mathcal{S}(t/r)^r\rho_0]-\mathrm{Tr}[Oe^{t\mathcal{L}}\rho_0]|$",
        fontsize=21,
    )
    ax.set_xticks(np.log10(SCALE_LIST))
    ax.set_xticklabels(SCALE_LIST)
    legend_order = [1, 3, 0, 2]
    legend = ax.legend(
        title=legend_title(initial, observable),
        title_fontsize=18,
        handles=[handles[i] for i in legend_order],
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=18,
    )
    legend.get_title().set_multialignment("center")
    legend.get_title().set_ha("center")
    ax.set_title("Trotter error scaling with Richardson extrapolation", fontsize=22)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()

    out_path = OUT_DIR / out_name
    fig.savefig(out_path)
    plt.close(fig)
    return out_path, summary


def main():
    jobs = [
        ("111", "magnetization_z", "heisenberg-extrapolation-5-111-amplitude_damping_scaling.pdf"),
        ("000", "magnetization_z", "heisenberg-extrapolation-5-000-amplitude_damping_scaling.pdf"),
        ("111", "nearest_neighbor_zz", "heisenberg-extrapolation-5-111-amplitude_damping_zz_scaling.pdf"),
        ("000", "nearest_neighbor_zz", "heisenberg-extrapolation-5-000-amplitude_damping_zz_scaling.pdf"),
    ]
    all_summary = {}
    for initial, observable, out_name in jobs:
        out_path, summary = make_plot(initial, observable, out_name)
        all_summary[out_name] = summary
        print(f"Wrote {out_path.relative_to(PROJECT_ROOT)}")

    summary_path = OUT_DIR / "heisenberg_extrapolation_observable_scaling_summary.json"
    summary_path.write_text(json.dumps(all_summary, indent=2))
    print(f"Wrote {summary_path.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
