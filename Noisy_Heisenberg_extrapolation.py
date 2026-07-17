import argparse
import json
from pathlib import Path

import numpy as np

from Noisy_Heisenberg_rho import (
    GAMMAS,
    INITIAL_DIRS,
    NOISE_LABELS,
    build_noisy_heisenberg,
    check_density,
    exact_rho,
    initial_rho,
    propagator_on_rho_0,
    trace_norm_hermitian,
)


R_BASES = (1, 2, 3, 4, 5, 6, 7)
SCALE_LIST = tuple(range(1, 11))
DEFAULT_INITIALS = ("1", "0")


def extrapolated_dir(data_root, noise, initial, gamma):
    return data_root / f"noisy_heisenberg_{noise}_extrapolated" / INITIAL_DIRS[initial].replace(
        "initial", "extrapolated"
    ) / f"gamma_{gamma}"


def save_density_outputs(out_dir, N, r, rho_exact, rho_t):
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / f"rho_superexact_N_{N}_no_Euler.npy", rho_exact)
    np.save(out_dir / f"rho_N_{N}_r_{r}_no_Euler.npy", rho_t)


def run_experiment(args):
    r_values = sorted({scale * base_r for scale in args.scales for base_r in args.r_bases})
    rows = []

    for gamma in args.gammas:
        H_terms, L_ops = build_noisy_heisenberg(args.N, args.J, args.h, gamma, args.noise)
        for initial in args.initials:
            rho_0 = initial_rho(args.N, initial)
            rho_exact = exact_rho(rho_0, args.N, H_terms, L_ops, args.t)
            check_density(f"exact N={args.N} initial={initial} gamma={gamma:g}", rho_exact)

            out_dir = extrapolated_dir(args.data_root, args.noise, initial, gamma)
            out_dir.mkdir(parents=True, exist_ok=True)
            np.save(out_dir / f"rho_superexact_N_{args.N}_no_Euler.npy", rho_exact)

            previous_error = None
            for r in r_values:
                rho_t = propagator_on_rho_0(
                    rho_0,
                    N=args.N,
                    H_terms=H_terms,
                    t=args.t,
                    r=r,
                    gamma=gamma,
                    noise=args.noise,
                )
                check_density(f"trotter N={args.N} initial={initial} gamma={gamma:g} r={r}", rho_t)
                np.save(out_dir / f"rho_N_{args.N}_r_{r}_no_Euler.npy", rho_t)
                error = max(trace_norm_hermitian(rho_exact - rho_t), 0.0)
                rows.append(
                    {
                        "model": "1D Heisenberg chain",
                        "noise": args.noise,
                        "gamma": gamma,
                        "N": args.N,
                        "r": r,
                        "initial": initial,
                        "trace_error": error,
                        "monotonicity_warning": bool(
                            previous_error is not None and error > previous_error * 1.02 + 1e-12
                        ),
                    }
                )
                previous_error = error

    return {
        "model": f"1D Heisenberg chain with {NOISE_LABELS[args.noise]}",
        "parameters": {
            "J": args.J,
            "h": args.h,
            "t": args.t,
            "noise": args.noise,
            "gammas": list(args.gammas),
            "N": args.N,
            "initials": list(args.initials),
            "r_values": r_values,
            "r_bases": list(args.r_bases),
            "scales": list(args.scales),
            "vectorization_order": "F",
            "liouvillian_convention": "vec(A rho B) = (B^T \\otimes A) vec(rho)",
            "output_note": "Density matrices are saved so observable choices can be made at plot time.",
        },
        "rows": rows,
        "sanity_warnings": [
            row
            for row in rows
            if row["monotonicity_warning"]
        ],
    }


def parse_args():
    parser = argparse.ArgumentParser(description="Heisenberg extrapolation density-matrix data")
    parser.add_argument("--noise", choices=("amplitude_damping", "dephasing"), default="amplitude_damping")
    parser.add_argument("--N", type=int, default=5)
    parser.add_argument("--J", type=float, default=1.0)
    parser.add_argument("--h", type=float, default=0.5)
    parser.add_argument("--t", type=float, default=0.15)
    parser.add_argument("--gammas", nargs="+", type=float, default=list(GAMMAS))
    parser.add_argument("--initials", nargs="+", choices=("0", "1", "+", "I"), default=list(DEFAULT_INITIALS))
    parser.add_argument("--r-bases", nargs="+", type=int, default=list(R_BASES))
    parser.add_argument("--scales", nargs="+", type=int, default=list(SCALE_LIST))
    parser.add_argument("--data-root", type=Path, default=Path("data"))
    return parser.parse_args()


def main():
    args = parse_args()
    args.gammas = tuple(args.gammas)
    args.initials = tuple(args.initials)
    args.r_bases = tuple(args.r_bases)
    args.scales = tuple(args.scales)

    results = run_experiment(args)
    metadata_dir = args.data_root / f"noisy_heisenberg_{args.noise}_extrapolated"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    summary_path = metadata_dir / "metadata.json"
    summary_path.write_text(json.dumps(results, indent=2), encoding="utf-8")

    if results["sanity_warnings"]:
        print("Sanity warnings:")
        for warning in results["sanity_warnings"][:20]:
            print(
                "  - "
                f"initial={warning['initial']}, gamma={warning['gamma']:g}, "
                f"r={warning['r']}, trace_error={warning['trace_error']:.3e}"
            )
        if len(results["sanity_warnings"]) > 20:
            print(f"  - ... {len(results['sanity_warnings']) - 20} more")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
