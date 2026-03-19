import argparse
from pathlib import Path

import numpy as np
import qutip as qt

try:
    from tqdm.auto import tqdm
except ImportError:  # pragma: no cover

    def tqdm(iterable, *args, **kwargs):
        return iterable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Exact evolution with first-order update for dissipative TFIM.")
    parser.add_argument("--N", type=int, default=5, help="Number of sites.")
    parser.add_argument("--J", type=float, default=1.0, help="Coupling strength J.")
    parser.add_argument("--h", type=float, default=0.5, help="Transverse field h.")
    parser.add_argument("--t", type=float, default=0.2, help="Total evolution time.")
    # parser.add_argument("--gamma", type=float, default=1.0, help="Dissipation strength.")
    parser.add_argument("--initial", type=str, default="0", help="Initial state")
    parser.add_argument("--step", type=int, default=10000, help="Number of steps.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    N = args.N
    J = args.J
    h = args.h
    t = args.t
    # gamma = args.gamma
    step = args.step

    # 构建 Hamiltonian
    si_x = [qt.tensor([qt.qeye(2)] * i + [qt.sigmax()] + [qt.qeye(2)] * (N - i - 1)) for i in range(N)]
    si_z = [qt.tensor([qt.qeye(2)] * i + [qt.sigmaz()] + [qt.qeye(2)] * (N - i - 1)) for i in range(N)]

    H_x = 0
    H_z = 0
    # 相互作用项
    for i in range(N - 1):
        H_z += -J * si_z[i] * si_z[i + 1]
    # 横场项
    for i in range(N):
        H_x += -h * si_x[i]

    for gamma in [0.1, 1.0]:
        # 密度矩阵 rho = |psi><psi|
        if args.initial == "0":
            psi = qt.tensor([qt.basis(2, 0) for _ in range(args.N)])
            rho_0 = qt.ket2dm(psi)
            gamma_dir = Path(f"./|000>initial/gamma_{gamma}")
        elif args.initial == "1":
            psi = qt.tensor([qt.basis(2, 1) for _ in range(args.N)])
            rho_0 = qt.ket2dm(psi)
            gamma_dir = Path(f"./|111>initial/gamma_{gamma}")
        elif args.initial == "+":
            plus = qt.basis(2, 0) + qt.basis(2, 1)  # .unit()
            psi = qt.tensor([plus for _ in range(args.N)])
            rho_0 = qt.ket2dm(psi)
            gamma_dir = Path(f"./|+++>initial/gamma_{gamma}")
        elif args.initial == "01":
            psi = qt.tensor([qt.basis(2, i % 2) for i in range(args.N)])
            rho_0 = qt.ket2dm(psi)
            gamma_dir = Path(f"./|010>initial/gamma_{gamma}")
        elif args.initial == "I":
            rho_0 = qt.qeye([2] * args.N)  # / (2**args.N)
            gamma_dir = Path(f"./|mixed>initial/gamma_{gamma}")

        else:
            raise ValueError("Invalid initial state option.")

        rho_exact = rho_0.copy()
        dt = t / step
        H = H_x + H_z
        L_ops = [np.sqrt(gamma) * qt.tensor([qt.qeye(2)] * i + [qt.sigmam()] + [qt.qeye(2)] * (N - i - 1)) for i in range(N)]
        for _ in tqdm(
            range(step),
            desc=f"Euler N={N} gamma={gamma} init={args.initial}",
        ):
            rho_delta = 0
            # Unitary part: -i [H, rho]
            rho_delta += -1j * (H * rho_exact - rho_exact * H) * dt
            for L_i in L_ops:
                rho_delta += (L_i @ rho_exact @ L_i.dag() - 1 / 2 * (L_i.dag() @ L_i @ rho_exact + rho_exact @ L_i.dag() @ L_i)) * dt

            rho_exact += rho_delta

        if args.initial in ["+", "I"]:
            rho_exact = rho_exact / (2**args.N)

        gamma_dir.mkdir(parents=True, exist_ok=True)
        np.save(gamma_dir / f"rho_exact_N_{N}_step_{step}.npy", rho_exact.full())


if __name__ == "__main__":
    main()
