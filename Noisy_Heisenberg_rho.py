import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse import csc_matrix, eye, kron
from scipy.sparse.linalg import expm_multiply

try:
    from tqdm.auto import tqdm
except ImportError:  # pragma: no cover

    def tqdm(iterable, *args, **kwargs):
        return iterable


GAMMAS = (0.1, 1.0)
R_LIST = (4, 6, 8, 10)
N_LIST = tuple(range(4, 9))
INITIALS = ("0", "1", "+", "I")

INITIAL_DIRS = {
    "0": "|000>initial",
    "1": "|111>initial",
    "+": "|+++>initial",
    "I": "|mixed>initial",
}
NOISE_LABELS = {
    "amplitude_damping": "local amplitude damping",
    "dephasing": "local dephasing",
}


SX = csc_matrix(np.array([[0, 1], [1, 0]], dtype=complex))
SY = csc_matrix(np.array([[0, -1j], [1j, 0]], dtype=complex))
SZ = csc_matrix(np.array([[1, 0], [0, -1]], dtype=complex))
# Explicit |0><1| jump for amplitude damping. Avoid sigmam/sigmap-style names
# here because QuTiP's naming is easy to misread in this basis convention.
AMP_DAMP_JUMP = csc_matrix(np.array([[0, 1], [0, 0]], dtype=complex))
P0 = csc_matrix(np.array([[1, 0], [0, 0]], dtype=complex))
P1 = csc_matrix(np.array([[0, 0], [0, 1]], dtype=complex))
ID2 = eye(2, format="csc", dtype=complex)


def local_operator(N, site_ops):
    op = csc_matrix([[1]], dtype=complex)
    for site in range(N):
        op = kron(op, site_ops.get(site, ID2), format="csc")
    return op


def build_noisy_heisenberg(N, J, h, gamma, noise):
    dim = 2**N
    H_xx = csc_matrix((dim, dim), dtype=complex)
    H_yy = csc_matrix((dim, dim), dtype=complex)
    H_zz = csc_matrix((dim, dim), dtype=complex)
    H_x = csc_matrix((dim, dim), dtype=complex)

    for i in range(N - 1):
        H_xx += J * local_operator(N, {i: SX, i + 1: SX})
        H_yy += J * local_operator(N, {i: SY, i + 1: SY})
        H_zz += J * local_operator(N, {i: SZ, i + 1: SZ})

    for i in range(N):
        H_x += h * local_operator(N, {i: SX})

    L_ops = []
    for i in range(N):
        if noise == "amplitude_damping":
            jump = local_operator(N, {i: AMP_DAMP_JUMP})
        elif noise == "dephasing":
            jump = local_operator(N, {i: SZ})
        else:
            raise ValueError(f"Invalid noise model: {noise}")
        L_ops.append(np.sqrt(gamma) * jump)

    return (H_xx, H_yy, H_zz, H_x), L_ops


def build_liouvillian(N, H_terms, L_ops):
    dim = 2**N
    ident = eye(dim, format="csc", dtype=complex)
    H = sum(H_terms, csc_matrix((dim, dim), dtype=complex))
    L = -1j * (kron(ident, H, format="csc") - kron(H.T, ident, format="csc"))

    for jump in L_ops:
        number = jump.conj().T @ jump
        L += (
            kron(jump.conj(), jump, format="csc")
            - 0.5 * kron(ident, number, format="csc")
            - 0.5 * kron(number.T, ident, format="csc")
        )
    return L


def initial_rho(N, initial):
    dim = 2**N
    if initial == "0":
        ket = np.zeros(dim, dtype=complex)
        ket[0] = 1.0
        return np.outer(ket, ket.conj())
    if initial == "1":
        ket = np.zeros(dim, dtype=complex)
        ket[-1] = 1.0
        return np.outer(ket, ket.conj())
    if initial == "+":
        ket = np.ones(dim, dtype=complex) / np.sqrt(dim)
        return np.outer(ket, ket.conj())
    if initial == "I":
        return np.eye(dim, dtype=complex) / dim
    raise ValueError("Invalid initial state option.")


def dephasing_channel(rho, z_op, gamma, dt):
    decay = np.exp(-2.0 * gamma * dt)
    return 0.5 * ((1.0 + decay) * rho + (1.0 - decay) * (z_op @ rho @ z_op))


def amplitude_damping_kraus(N, gamma, dt):
    eta = np.exp(-gamma * dt)
    e0_local = P0 + np.sqrt(eta) * P1
    e1_local = np.sqrt(1.0 - eta) * AMP_DAMP_JUMP
    kraus = []
    for site in range(N):
        e0 = local_operator(N, {site: e0_local}).toarray()
        e1 = local_operator(N, {site: e1_local}).toarray()
        kraus.append((e0, e1))
    return kraus


def apply_noise(rho, noise, op, gamma, dt):
    if noise == "dephasing":
        return dephasing_channel(rho, op, gamma, dt)
    if noise == "amplitude_damping":
        e0, e1 = op
        return e0 @ rho @ e0.conj().T + e1 @ rho @ e1.conj().T
    raise ValueError(f"Invalid noise model: {noise}")


def propagator_on_rho_0(rho_0, N, H_terms, t, r=1, gamma=0.1, noise="amplitude_damping"):
    rho = rho_0.copy()
    s = t / r
    half_s = s / 2.0
    H_unitaries = [expm(-1j * H_term.toarray() * half_s) for H_term in H_terms]

    if noise == "dephasing":
        noise_ops = [local_operator(N, {site: SZ}).toarray() for site in range(N)]
    elif noise == "amplitude_damping":
        noise_ops = amplitude_damping_kraus(N, gamma, half_s)
    else:
        raise ValueError(f"Invalid noise model: {noise}")

    channels = [("H", U) for U in H_unitaries] + [(noise, op) for op in noise_ops]
    for _ in range(r):
        for kind, op in channels:
            if kind == "H":
                rho = op @ rho @ op.conj().T
            else:
                rho = apply_noise(rho, noise, op, gamma, half_s)
        for kind, op in reversed(channels):
            if kind == "H":
                rho = op @ rho @ op.conj().T
            else:
                rho = apply_noise(rho, noise, op, gamma, half_s)
    return rho


def exact_rho(rho_0, N, H_terms, L_ops, t):
    L = build_liouvillian(N, H_terms, L_ops)
    rho_vec = expm_multiply(t * L, rho_0.reshape(-1, order="F"))
    return rho_vec.reshape(rho_0.shape, order="F")


def trace_norm_hermitian(mat):
    mat = 0.5 * (mat + mat.conj().T)
    return float(np.sum(np.abs(np.linalg.eigvalsh(mat))))


def check_density(label, rho, tol=5e-8):
    hermiticity = np.linalg.norm(rho - rho.conj().T, ord="fro")
    trace_error = abs(np.trace(rho) - 1.0)
    if hermiticity > tol:
        raise ValueError(f"{label}: Hermiticity check failed ({hermiticity:.3e})")
    if trace_error > tol:
        raise ValueError(f"{label}: trace preservation check failed ({trace_error:.3e})")


def data_dir(data_root, noise, initial, gamma):
    return data_root / f"noisy_heisenberg_{noise}" / INITIAL_DIRS[initial] / f"gamma_{gamma}"


def save_density_outputs(out_dir, N, r, rho_exact, rho_t):
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / f"rho_superexact_N_{N}.npy", rho_exact)
    np.save(out_dir / f"rho_N_{N}_r_{r}.npy", rho_t)


def monotonicity_warnings(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["noise"], row["gamma"], row["N"], row["initial"])].append(row)

    warnings = []
    for key, group in grouped.items():
        ordered = sorted(group, key=lambda row: row["r"])
        for prev, curr in zip(ordered, ordered[1:]):
            if curr["trace_error"] > prev["trace_error"] * 1.02 + 1e-12:
                noise, gamma, N, initial = key
                warnings.append(
                    f"error increased with r for noise={noise}, gamma={gamma:g}, "
                    f"N={N}, initial={initial}: r={prev['r']} -> r={curr['r']}"
                )
    return warnings


def run_experiment(args):
    rows = []
    initials = INITIALS if args.initial == "all" else (args.initial,)

    for gamma in args.gammas:
        for N in tqdm(args.N_list, desc=f"{args.noise}, gamma={gamma:g}"):
            H_terms, L_ops = build_noisy_heisenberg(N, args.J, args.h, gamma, args.noise)
            for initial in initials:
                rho_0 = initial_rho(N, initial)
                rho_exact = exact_rho(rho_0, N, H_terms, L_ops, args.t)
                check_density(f"exact N={N} initial={initial} gamma={gamma:g}", rho_exact)
                out_dir = data_dir(args.data_root, args.noise, initial, gamma)

                for r in args.r_list:
                    rho_t = propagator_on_rho_0(
                        rho_0,
                        N=N,
                        H_terms=H_terms,
                        t=args.t,
                        r=r,
                        gamma=gamma,
                        noise=args.noise,
                    )
                    check_density(
                        f"trotter N={N} initial={initial} gamma={gamma:g} r={r}",
                        rho_t,
                    )
                    save_density_outputs(out_dir, N, r, rho_exact, rho_t)
                    error = trace_norm_hermitian(rho_exact - rho_t)
                    rows.append(
                        {
                            "model": "1D Heisenberg chain",
                            "noise": args.noise,
                            "gamma": gamma,
                            "N": N,
                            "r": r,
                            "initial": initial,
                            "trace_error": max(error, 0.0),
                        }
                    )

    results = {
        "model": f"1D Heisenberg chain with {NOISE_LABELS[args.noise]}",
        "parameters": {
            "J": args.J,
            "h": args.h,
            "t": args.t,
            "noise": args.noise,
            "gammas": list(args.gammas),
            "r_list": list(args.r_list),
            "N_list": list(args.N_list),
            "initials": list(initials),
            "vectorization_order": "F",
            "liouvillian_convention": "vec(A rho B) = (B^T \\otimes A) vec(rho)",
        },
        "rows": rows,
        "sanity_warnings": monotonicity_warnings(rows),
    }
    return results


def parse_args():
    parser = argparse.ArgumentParser(description="Noisy Heisenberg rho simulation")
    parser.add_argument("--noise", choices=("amplitude_damping", "dephasing"), default="amplitude_damping")
    parser.add_argument("--N", type=int, help="Single system size. Overrides --N-min/--N-max.")
    parser.add_argument("--N-min", type=int, default=min(N_LIST))
    parser.add_argument("--N-max", type=int, default=max(N_LIST))
    parser.add_argument("--J", type=float, default=1.0, help="Heisenberg coupling strength")
    parser.add_argument("--h", type=float, default=0.5, help="X-field strength")
    parser.add_argument("--t", type=float, default=0.15, help="Simulation time")
    parser.add_argument("--r", type=int, help="Single Trotter step count. Overrides --r-list.")
    parser.add_argument("--r-list", nargs="+", type=int, default=list(R_LIST))
    parser.add_argument("--gammas", nargs="+", type=float, default=list(GAMMAS))
    parser.add_argument("--initial", choices=("all",) + INITIALS, default="all")
    parser.add_argument("--data-root", type=Path, default=Path("data"))
    return parser.parse_args()


def main():
    args = parse_args()
    args.N_list = (args.N,) if args.N is not None else tuple(range(args.N_min, args.N_max + 1))
    args.r_list = (args.r,) if args.r is not None else tuple(args.r_list)
    args.gammas = tuple(args.gammas)

    results = run_experiment(args)
    metadata_dir = args.data_root / f"noisy_heisenberg_{args.noise}"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    summary_path = metadata_dir / "metadata.json"
    summary_path.write_text(json.dumps(results, indent=2), encoding="utf-8")

    if results["sanity_warnings"]:
        print("Sanity warnings:")
        for warning in results["sanity_warnings"]:
            print(f"  - {warning}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
