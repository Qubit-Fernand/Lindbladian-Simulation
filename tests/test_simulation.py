"""
test_simulation.py - 验证 Lindbladian 模拟的正确性

运行方法:
    cd jobs
    python -m pytest test_simulation.py -v
"""

import sys
from pathlib import Path

import numpy as np
import pytest
import qutip as qt

sys.path.insert(0, str(Path(__file__).parent.parent / "jobs"))
from archive.extrapolation_0to1 import build_dissipative_tfim, propagator_on_rho_0

# ─── 公共参数（N=2，跑得快） ───
N = 2
J, h, gamma, t = 1.0, 0.5, 0.1, 0.2


# ════════════════════════════════════════════
#  第一组：算符的基本结构
# ════════════════════════════════════════════


def test_hamiltonian_is_hermitian():
    """H_x 和 H_z 应该都是厄米算符（H = H†）"""
    H_x, H_z, _ = build_dissipative_tfim(N, J, h, gamma)
    assert H_x.isherm, "H_x 不是厄米的"
    assert H_z.isherm, "H_z 不是厄米的"


def test_collapse_operators_count():
    """L_ops 的数量应等于格点数 N"""
    _, _, L_ops = build_dissipative_tfim(N, J, h, gamma)
    assert len(L_ops) == N


def test_collapse_operator_scaling():
    """
    L_i = sqrt(gamma) * sigmam，所以对单格点 N=1：
        tr(L† L) = gamma * tr(sigmam† sigmam) = gamma * 1 = gamma
    注意：N>1 时 tensor product 展开后 tr(L_i† L_i) = gamma * 2^(N-1)，
    因为其他位点的单位矩阵 tr(I_2) = 2 会累乘。
    """
    # 用 N=1 验证最干净，直接对单格点
    _, _, L_ops_1 = build_dissipative_tfim(1, J=1.0, h=0.0, gamma=gamma)
    norm_sq = (L_ops_1[0].dag() * L_ops_1[0]).tr().real
    assert abs(norm_sq - gamma) < 1e-10, f"tr(L†L) = {norm_sq:.6f}，期望 {gamma}"


@pytest.mark.xfail(reason="已知 bug：N=1 时 H_z 循环不执行（range(0)），H_z 保持为 Python int 0，没有 .expm() 方法")
def test_hz_is_qobj_for_n1():
    """
    bug 检查：N=1 时 H_z 应是零 Qobj，而非 Python int 0。
    修复方案：将初始化改为 H_z = 0 * qt.qeye([2]*N) 或 qt.Qobj(...)。
    """
    H_x, H_z, _ = build_dissipative_tfim(1, J=1.0, h=0.5, gamma=gamma)
    assert isinstance(H_z, qt.Qobj), f"H_z 类型是 {type(H_z)}，而非 Qobj（N=1 时 range(N-1)=range(0) 不执行循环）"


# ════════════════════════════════════════════
#  第二组：密度矩阵基本性质
# ════════════════════════════════════════════


@pytest.fixture
def evolved_state():
    """N=2 的基本设置：初态 |00>, 演化参数 t=0.2, r=4"""
    psi = qt.tensor([qt.basis(2, 0)] * N)
    rho_0 = qt.ket2dm(psi)
    H_x, H_z, L_ops = build_dissipative_tfim(N, J, h, gamma)
    rho_t = propagator_on_rho_0(rho_0, N, H_x, H_z, t=t, r=4, L_ops=L_ops)
    return rho_0, rho_t, H_x, H_z, L_ops


def test_trace_preserved(evolved_state):
    """Lindblad 演化保迹：tr(ρ_t) = 1"""
    _, rho_t, _, _, _ = evolved_state
    tr = rho_t.tr().real
    assert abs(tr - 1.0) < 1e-8, f"迹 = {tr:.8f}（期望 1）"


def test_density_matrix_hermitian(evolved_state):
    """演化后密度矩阵仍是厄米矩阵：ρ = ρ†"""
    _, rho_t, _, _, _ = evolved_state
    diff = (rho_t - rho_t.dag()).norm()
    assert diff < 1e-10, f"ρ_t 不厄米，||ρ-ρ†|| = {diff:.2e}"


def test_density_matrix_positive_semidefinite(evolved_state):
    """密度矩阵应半正定：所有本征值 ≥ 0"""
    _, rho_t, _, _, _ = evolved_state
    eigvals = rho_t.eigenenergies()
    assert np.all(eigvals >= -1e-8), f"存在负本征值：min = {eigvals.min():.2e}"


# ════════════════════════════════════════════
#  第三组：物理极限检验
# ════════════════════════════════════════════


def test_zero_time_is_identity():
    """t=0 时演化等于恒等：ρ_t(t=0) = ρ_0"""
    psi = qt.tensor([qt.basis(2, 0)] * N)
    rho_0 = qt.ket2dm(psi)
    H_x, H_z, L_ops = build_dissipative_tfim(N, J, h, gamma)
    rho_t = propagator_on_rho_0(rho_0, N, H_x, H_z, t=0.0, r=1, L_ops=L_ops)
    diff = (rho_t - rho_0).norm()
    assert diff < 1e-12, f"t=0 时结果偏离初态，||ρ_t - ρ_0|| = {diff:.2e}"


def test_no_dissipation_preserves_purity():
    """gamma=0 时只有幺正演化，纯态应保持纯度 tr(ρ²) = 1"""
    psi = qt.tensor([qt.basis(2, 0)] * N)
    rho_0 = qt.ket2dm(psi)
    H_x, H_z, L_ops_zero = build_dissipative_tfim(N, J, h, gamma=0.0)
    rho_t = propagator_on_rho_0(rho_0, N, H_x, H_z, t=t, r=4, L_ops=L_ops_zero)
    purity = (rho_t * rho_t).tr().real
    assert abs(purity - 1.0) < 1e-8, f"gamma=0 时纯度 = {purity:.8f}（期望 1）"


def test_single_qubit_decay_rate_exact():
    """
    N=1, J=0, h=0, 初态 |0><0|（激发态）：纯自发辐射。
    精确解（用 qt.liouvillian 直接算，不依赖我们的 propagator）：
        激发态布居数 ρ_00(t) = e^{-gamma * t}

    QuTiP 约定：basis(2, 0) = |↑> = 激发态；
                 sigmam() = |1><0|，把激发态 |0> 衰减到基态 |1>。
    这验证了模型构建（build_dissipative_tfim）是正确的。
    注意：这里不调用 propagator_on_rho_0（N=1 时 H_z=int(0) 有 bug，见 test_hz_is_qobj_for_n1）。
    """
    N1 = 1
    rho_0 = qt.ket2dm(qt.basis(2, 0))  # |0><0| = 激发态（QuTiP 约定）
    H_x, _, L_ops = build_dissipative_tfim(N1, J=0.0, h=0.0, gamma=gamma)
    L = qt.liouvillian(H_x, L_ops)  # H=0（h=0），只有耗散项
    rho_t = qt.vector_to_operator((L * t).expm() * qt.operator_to_vector(rho_0))

    excited_pop = qt.expect(qt.ket2dm(qt.basis(2, 0)), rho_t)  # 测激发态 |0> 布居数
    expected = np.exp(-gamma * t)
    assert abs(excited_pop - expected) < 1e-8, f"激发态布居数：实际 {excited_pop:.8f}，精确解 {expected:.8f}"


def test_jump_operator_causes_decay_not_excitation():
    """
    extrapolation.py 用的是 sigmam（降算符 = 衰减方向）。
    从全激发态 |00><00| 出发，|00> 布居数应随时间减少。
    如果这个测试失败，说明跳跃算符方向搞反了（sigmam 换成了 sigmap）。

    QuTiP 约定：basis(2, 0) = 激发态，所以全激发 = |00> = index 0。
    """
    # 初态：|00><00|（全激发，QuTiP 约定 basis(2,0) = 激发态）
    psi = qt.tensor([qt.basis(2, 0)] * N)
    rho_0 = qt.ket2dm(psi)
    H_x, H_z, L_ops = build_dissipative_tfim(N, J=0.0, h=0.0, gamma=gamma)
    rho_t = propagator_on_rho_0(rho_0, N, H_x, H_z, t=0.5, r=4, L_ops=L_ops)

    pop_before = rho_0.full()[0, 0].real  # |00> 布居数 = 1.0
    pop_after = rho_t.full()[0, 0].real  # 演化后 |00> 布居数，应减少
    assert pop_after < pop_before, f"|00> 布居数未减少（{pop_before:.4f} → {pop_after:.4f}），" "跳跃算符方向可能错误（应是 sigmam，而非 sigmap）"


# ════════════════════════════════════════════
#  第四组：Trotter 误差收敛性
# ════════════════════════════════════════════


def _get_exact_rho(rho_0, H_x, H_z, L_ops):
    """用 qt.liouvillian 精确演化，作为参考基准"""
    L = qt.liouvillian(H_x + H_z, L_ops)
    return qt.vector_to_operator((L * t).expm() * qt.operator_to_vector(rho_0))


def test_trotter_error_decreases_with_r():
    """增大 Trotter 步数 r，与精确解的误差应单调递减"""
    psi = qt.tensor([qt.basis(2, 0)] * N)
    rho_0 = qt.ket2dm(psi)
    H_x, H_z, L_ops = build_dissipative_tfim(N, J, h, gamma)
    rho_exact = _get_exact_rho(rho_0, H_x, H_z, L_ops)

    r_list = [1, 2, 4, 8]
    errors = [(propagator_on_rho_0(rho_0, N, H_x, H_z, t=t, r=r, L_ops=L_ops) - rho_exact).norm() for r in r_list]

    for i in range(len(r_list) - 1):
        assert errors[i] > errors[i + 1], f"r={r_list[i]}→{r_list[i+1]}：误差未减小 " f"({errors[i]:.2e} → {errors[i+1]:.2e})"


def test_trotter_second_order_convergence():
    """
    对称 Suzuki-Trotter 分裂是 2 阶的：误差 ~ C * (t/r)²
    因此 error(r=2) / error(r=4) ≈ 4（允许范围 [2, 8]）
    """
    psi = qt.tensor([qt.basis(2, 0)] * N)
    rho_0 = qt.ket2dm(psi)
    H_x, H_z, L_ops = build_dissipative_tfim(N, J, h, gamma)
    rho_exact = _get_exact_rho(rho_0, H_x, H_z, L_ops)

    def err(r):
        return (propagator_on_rho_0(rho_0, N, H_x, H_z, t=t, r=r, L_ops=L_ops) - rho_exact).norm()

    ratio = err(2) / err(4)
    assert 2.0 < ratio < 8.0, f"误差比 = {ratio:.2f}，期望约 4（2 阶 Suzuki-Trotter 收敛）"


# ════════════════════════════════════════════
#  第五组：两种实现的一致性检验
# ════════════════════════════════════════════


def test_propagator_on_rho_vs_superoperator():
    """
    propagator_on_rho_0（直接作用于 ρ）
    与 build_superoperator（超算符矩阵）应给出相同结果。
    """
    from Dissipative_TFIM_diamond import (
        build_dissipative_tfim as build_diamond,
        build_superoperator,
    )

    psi = qt.tensor([qt.basis(2, 0)] * N)
    rho_0 = qt.ket2dm(psi)
    H_x, H_z, L_ops = build_diamond(N, J, h, gamma)  # diamond.py 也用 sigmam，一致

    # 方法一：直接在 ρ 上演化
    rho_direct = propagator_on_rho_0(rho_0, N, H_x, H_z, t=t, r=2, L_ops=L_ops)

    # 方法二：构建超算符再作用
    PF = build_superoperator(N, H_x, H_z, t=t, L_ops=L_ops, r=2, exponential=False)
    rho_super = qt.vector_to_operator(PF * qt.operator_to_vector(rho_0))

    diff = (rho_direct - rho_super).norm()
    assert diff < 1e-8, f"两种方法结果不一致，||差|| = {diff:.2e}"


@pytest.mark.xfail(reason="已知不一致：extrapolation_new.py 用 sigmap（激发），extrapolation.py 用 sigmam（衰减），物理方向相反")
def test_jump_operator_consistent_across_files():
    """
    extrapolation.py 和 extrapolation_new.py 的跳跃算符应相同。
    当前这个测试会失败——这说明两个文件的物理模型不一致，需要确认哪个是对的。
    """
    import archive.extrapolation_0to1 as old
    import extrapolation as new

    _, _, L_old = old.build_dissipative_tfim(N, J, h, gamma)
    _, _, L_new = new.build_dissipative_tfim(N, J, h, gamma)

    for i, (L_o, L_n) in enumerate(zip(L_old, L_new)):
        diff = (L_o - L_n).norm()
        assert diff < 1e-10, f"L_ops[{i}] 不一致：extrapolation.py 用 sigmam（衰减），" f"extrapolation_new.py 用 sigmap（激发），物理意义相反！"
