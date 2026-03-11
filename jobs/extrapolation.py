"""
For N<=6, no Euler approximation
"""

import argparse
import numpy as np
import qutip as qt
from pathlib import Path

try:
    from tqdm.auto import tqdm
except ImportError:  # pragma: no cover

    def tqdm(iterable, *args, **kwargs):
        return iterable


def build_dissipative_tfim(N, J, h, gamma):
    # 1. 构建 Hamiltonian
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

    H = H_x + H_z

    # 2. 构建耗散算符 (Collapse Operators)
    # 这里使用 sigma_minus (自旋弛豫)
    L_ops = []
    First_version = 1
    for i in range(N):
        # if First_version:  # sigmam() |1><0|
        #     L_i = qt.tensor([qt.qeye(2)] * i + [qt.sigmam()] + [qt.qeye(2)] * (N - i - 1))
        # else:
        L_i = qt.tensor([qt.qeye(2)] * i + [qt.sigmap()] + [qt.qeye(2)] * (N - i - 1))
        L_ops.append(np.sqrt(gamma) * L_i)  # why square root here

    return H_x, H_z, L_ops


def propagator_on_rho_0(rho_0, N, H_x, H_z, t, r=1, L_ops=None):
    rho = rho_0.copy()
    s = t / r
    """
     Guarantee that the both e^{H} and dissipative are exact,
    """
    for _ in range(r):
        rho = (-1j * H_x * s / 2).expm() @ rho @ (1j * H_x * s / 2).expm()
        rho = (-1j * H_z * s / 2).expm() @ rho @ (1j * H_z * s / 2).expm()

        for L_i in L_ops:
            D_i = qt.lindblad_dissipator(L_i)
            rho = qt.propagator(D_i, s / 2)(rho)

        for L_i in reversed(L_ops):
            D_i = qt.lindblad_dissipator(L_i)
            rho = qt.propagator(D_i, s / 2)(rho)

        rho = (-1j * H_z * s / 2).expm() @ rho @ (1j * H_z * s / 2).expm()
        rho = (-1j * H_x * s / 2).expm() @ rho @ (1j * H_x * s / 2).expm()

    return rho


def parse_args():
    parser = argparse.ArgumentParser(description="Dissipative TFIM simulation")
    parser.add_argument("--N", type=int, default=4, help="Number of spins")
    parser.add_argument("--J", type=float, default=1, help="Coupling strength")
    parser.add_argument("--h", type=float, default=0.5, help="Transverse field")
    # parser.add_argument("--gamma", type=float, default=0.1, help="Dissipation rate")
    parser.add_argument("--t", type=float, default=0.2, help="Time step")
    # parser.add_argument("--r", type=int, default=1, help="Number of Trotter steps")
    parser.add_argument("--initial", type=str, default="1", help="Initial state")
    parser.add_argument(
        "--solver",
        type=str,
        help="Solver name for diamond norm (e.g. SCS, CVXOPT)",
    )
    # parser.add_argument("--scale", type=int, default="1", help="scale factor for Trotter steps")
    return parser.parse_args()


def main():
    args = parse_args()

    for gamma in [0.1, 1.0]:

        if args.initial == "0":
            psi = qt.tensor([qt.basis(2, 0) for _ in range(args.N)])
            rho_0 = qt.ket2dm(psi)
            gamma_dir = Path(f"./|000>initial/gamma_{gamma}")
        elif args.initial == "1":
            psi = qt.tensor([qt.basis(2, 1) for _ in range(args.N)])
            rho_0 = qt.ket2dm(psi)
            gamma_dir = Path(f"./|111>initial/gamma_{gamma}")
        elif args.initial == "+":
            plus = (qt.basis(2, 0) + qt.basis(2, 1)).unit()
            psi = qt.tensor([plus for _ in range(args.N)])
            rho_0 = qt.ket2dm(psi)
            gamma_dir = Path(f"./|+++>initial/gamma_{gamma}")
        elif args.initial == "01":
            psi = qt.tensor([qt.basis(2, i % 2) for i in range(args.N)])
            rho_0 = qt.ket2dm(psi)
            gamma_dir = Path(f"./|010>initial/gamma_{gamma}")
        elif args.initial == "I":
            rho_0 = qt.qeye([2] * args.N) / (2**args.N)
            gamma_dir = Path(f"./|mixed>initial/gamma_{gamma}")

        else:
            raise ValueError("Invalid initial state option.")

        H_x, H_z, L_ops = build_dissipative_tfim(args.N, J=args.J, h=args.h, gamma=gamma)
        L = qt.liouvillian(H_x + H_z, L_ops)
        exact = (L * args.t).expm()

        rho_superexact = exact * qt.operator_to_vector(rho_0)
        rho_superexact = qt.vector_to_operator(rho_superexact)

        gamma_dir.mkdir(parents=True, exist_ok=True)
        np.save(gamma_dir / f"rho_superexact_N_{args.N}_no_Euler.npy", rho_superexact.full())

        r_values = sorted({scale * base_r for scale in range(1, 11) for base_r in [1, 2, 3, 4, 5, 6, 7]})
        for r in r_values:
            rho_t = propagator_on_rho_0(rho_0, N=args.N, H_x=H_x, H_z=H_z, t=args.t, r=r, L_ops=L_ops)
            np.save(gamma_dir / f"rho_N_{args.N}_r_{r}_no_Euler.npy", rho_t.full())


if __name__ == "__main__":
    main()
