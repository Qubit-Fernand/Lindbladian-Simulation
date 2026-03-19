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
        if First_version:
            # sigmam() |1><0|
            L_i = qt.tensor([qt.qeye(2)] * i + [qt.sigmam()] + [qt.qeye(2)] * (N - i - 1))
        else:
            L_i = qt.tensor([qt.qeye(2)] * i + [qt.sigmap()] + [qt.qeye(2)] * (N - i - 1))
        L_ops.append(np.sqrt(gamma) * L_i)  # why square root here

    return H_x, H_z, L_ops


def propagator_on_rho_0(rho_0, N, H_x, H_z, t, r=1, L_ops=None, step=1000):
    rho = rho_0.copy()
    s = t / r
    """
     Guarantee that the e^{H} is exact, only dissipative part is approximated
    """
    Euler_H = False
    Euler_D = True
    if N >= 9:
        print("Use first-order Euler approximation for H", Euler_H)
        print("Use first-order Euler approximation for D/L_j", Euler_D)

    for _ in range(r):
        if Euler_H:
            for _ in range(step):
                delta = -1j * (H_x @ rho - rho @ H_x) * (s / (2 * step))
                rho += delta
            for _ in range(step):
                delta = -1j * (H_z @ rho - rho @ H_z) * (s / (2 * step))
                rho += delta
        else:
            rho = (-1j * H_x * s / 2).expm() @ rho @ (1j * H_x * s / 2).expm()
            rho = (-1j * H_z * s / 2).expm() @ rho @ (1j * H_z * s / 2).expm()

        # dissipative exponential always by 1st order Euler
        if Euler_D:
            for L_i in L_ops:
                for _ in tqdm(range(step)):
                    rho += (L_i @ rho @ L_i.dag() - 1 / 2 * (L_i.dag() @ L_i @ rho + rho @ L_i.dag() @ L_i)) * (s / (2 * step))

            for L_i in reversed(L_ops):
                for _ in range(step):
                    rho += (L_i @ rho @ L_i.dag() - 1 / 2 * (L_i.dag() @ L_i @ rho + rho @ L_i.dag() @ L_i)) * (s / (2 * step))

        if Euler_H:
            for _ in range(step):
                delta = -1j * (H_z @ rho - rho @ H_z) * (s / (2 * step))
                rho += delta
            for _ in range(step):
                delta = -1j * (H_x @ rho - rho @ H_x) * (s / (2 * step))
                rho += delta
        else:
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
    parser.add_argument("--r", type=int, default=1, help="Number of Trotter steps")
    parser.add_argument("--initial", type=str, default="0", help="Initial state")
    parser.add_argument(
        "--solver",
        type=str,
        help="Solver name for diamond norm (e.g. SCS, CVXOPT)",
    )
    parser.add_argument("--expH", type=bool, default=True, help="Use exponential for Hamiltonian")
    parser.add_argument("--step", type=int, default=1000, help="Euler steps for dissipator")
    return parser.parse_args()


def main():
    args = parse_args()
    for gamma in [0.1, 1.0]:
        H_x, H_z, L_ops = build_dissipative_tfim(args.N, J=args.J, h=args.h, gamma=gamma)

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

        rho_t = propagator_on_rho_0(
            rho_0,
            N=args.N,
            H_x=H_x,
            H_z=H_z,
            t=args.t,
            r=args.r,
            L_ops=L_ops,
            step=args.step,
        )

        if args.initial in ["+", "I"]:
            rho_t = rho_t / (2**args.N)

        gamma_dir.mkdir(parents=True, exist_ok=True)
        np.save(gamma_dir / f"rho_N_{args.N}_r_{args.r}.npy", rho_t.full())


if __name__ == "__main__":
    main()
