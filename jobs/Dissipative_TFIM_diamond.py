import argparse
import numpy as np
import qutip as qt
import cvxpy as cp
from pathlib import Path

try:
    from tqdm.auto import tqdm
except ImportError:  # pragma: no cover

    def tqdm(iterable, *args, **kwargs):
        return iterable


def build_dissipative_tfim(N, J, h, gamma):
    # 1. 构建 Hamiltonian
    si_x = [
        qt.tensor([qt.qeye(2)] * i + [qt.sigmax()] + [qt.qeye(2)] * (N - i - 1))
        for i in range(N)
    ]
    si_z = [
        qt.tensor([qt.qeye(2)] * i + [qt.sigmaz()] + [qt.qeye(2)] * (N - i - 1))
        for i in range(N)
    ]

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
    for i in range(N):
        # sigmam() 是降算符
        L_i = qt.tensor([qt.qeye(2)] * i + [qt.sigmam()] + [qt.qeye(2)] * (N - i - 1))
        L_ops.append(np.sqrt(gamma) * L_i)  # why square root here

    return H_x, H_z, L_ops


def build_superoperator(N, H_x, H_z, t, L_ops, exponential=False, r=1):
    PF = qt.to_super(qt.qeye([2] * N))
    s = t / r
    for _ in range(r):
        if exponential:
            PF = (qt.liouvillian(H_x, []) * s / 2).expm() * PF
            PF = (qt.liouvillian(H_z, []) * s / 2).expm() * PF
        else:
            PF = qt.to_super((-1j * H_x * s / 2).expm()) * PF
            PF = qt.to_super((-1j * H_z * s / 2).expm()) * PF

        for L_i in L_ops:
            D_i = qt.lindblad_dissipator(L_i)
            if exponential:
                PF = (D_i * (s / 2)).expm() * PF
            else:
                PF = qt.propagator(D_i, s / 2) * PF

        for L_i in reversed(L_ops):
            D_i = qt.lindblad_dissipator(L_i)
            if exponential:
                PF = (D_i * (s / 2)).expm() * PF
            else:
                PF = qt.propagator(D_i, s / 2) * PF

        if exponential:
            PF = (qt.liouvillian(H_z, []) * s / 2).expm() * PF
            PF = (qt.liouvillian(H_x, []) * s / 2).expm() * PF

        else:
            PF = qt.to_super((-1j * H_z * s / 2).expm()) * PF
            PF = qt.to_super((-1j * H_x * s / 2).expm()) * PF

    return PF


def parse_args():
    parser = argparse.ArgumentParser(description="Dissipative TFIM simulation")
    parser.add_argument("--N", type=int, default=4, help="Number of spins")
    parser.add_argument("--J", type=float, default=1, help="Coupling strength")
    parser.add_argument("--h", type=float, default=0.5, help="Transverse field")
    parser.add_argument("--gamma", type=float, default=0.1, help="Dissipation rate")
    parser.add_argument("--t", type=float, default=0.2, help="Time step")
    parser.add_argument("--r", type=int, default=2, help="Number of Trotter steps")

    parser.add_argument(
        "--solver",
        type=str,
        default="SCS",
        help="Solver name for diamond norm (e.g. SCS, CVXOPT)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    H_x, H_z, L_ops = build_dissipative_tfim(
        args.N, J=args.J, h=args.h, gamma=args.gamma
    )

    # Build propagator using Trotterization
    PF = build_superoperator(
        N=args.N, H_x=H_x, H_z=H_z, t=args.t, L_ops=L_ops, r=args.r
    )

    # Exact propagator
    L = qt.liouvillian(H_x + H_z, L_ops)
    exact = (L * args.t).expm()

    # Create output directory
    dir = Path(f"../data/diamond/gamma_{args.gamma}/txt")
    dir.mkdir(parents=True, exist_ok=True)

    # diamond norm distance
    dist = qt.dnorm(PF - exact, solver="MOSEK")  # args.solver
    print(f"Diamond Distance: {dist}")

    with open(dir / f"distance_N_{args.N}_r_{args.r}.txt", "w", encoding="utf-8") as f:
        f.write(f"{dist}\n")

    if 0:  # Compute trace norm for a sample initial state
        psi = qt.tensor([qt.basis(2, 0) for _ in range(args.N)])
        rho_0 = qt.ket2dm(psi)
        _ = rho_0  # kept for parity with the notebook; change initial state if needed

        rho_t = (PF - exact)(rho_0)  # do not need *, qutip sets it up
        trace_norm = rho_t.norm("tr")
        print(f"Trace norm of rho_t: {trace_norm}")

        with open(f"error_N_{args.N}_r_{args.r}.txt", "w", encoding="utf-8") as f:
            f.write(f"{trace_norm}\n")


if __name__ == "__main__":
    main()
