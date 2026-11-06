"""Verification of the scattering-matrix specification (smatrix_spec.py,
docs/PHASE2_METHOD_SELECTION.md section 2.3).

These tests use small TOY lattices (random matrices, a 1D chain). They
verify formulas and conventions, not physics. No physical solver result
is produced or implied.

Independent checks:
  * 1D chain: Fisher-Lee amplitudes (phases included) vs a transfer-matrix
    wave-function calculation that uses no Green functions;
  * continuum limit vs the analytic rectangular-barrier formula, which is
    itself checked against exact plane-wave matching.
Consistency checks on random multichannel toys: unitarity, S = S^T,
Caroli = Tr(t^dagger t), and the opposite-direction mode reciprocity.
"""
import numpy as np
import pytest

import smatrix_spec as sm
from analytic import rectangular_barrier_transmission


# ----------------------------------------------------------------------
# helpers: dense toy device (block tridiagonal), NOT the planned RGF solver
# ----------------------------------------------------------------------
def toy_green(E, Vs, E_perp, dx):
    """Retarded G of a toy device: slices j=1..N with real symmetric
    blocks Vs[j], hopping -t I, clean identical leads. Dense inversion."""
    t = sm.hopping(dx)
    N, ny = len(Vs), len(E_perp)
    modes = sm.lead_modes(E, E_perp, dx)
    H = np.zeros((N * ny, N * ny))
    for j, V in enumerate(Vs):
        s = slice(j * ny, (j + 1) * ny)
        H[s, s] = 2 * t * np.eye(ny) + np.diag(E_perp) + V
        if j + 1 < N:
            s2 = slice((j + 1) * ny, (j + 2) * ny)
            H[s, s2] = -t * np.eye(ny)
            H[s2, s] = -t * np.eye(ny)
    M = E * np.eye(N * ny, dtype=complex) - H
    M[:ny, :ny] -= np.diag(modes["sigma"])
    M[-ny:, -ny:] -= np.diag(modes["sigma"])
    G = np.linalg.inv(M)
    b = lambda i, j: G[i * ny:(i + 1) * ny, j * ny:(j + 1) * ny]
    return dict(G=G, G11=b(0, 0), GN1=b(N - 1, 0), G1N=b(0, N - 1), GNN=b(N - 1, N - 1), modes=modes)


def random_device(rng, N, ny, scale):
    Vs = []
    for _ in range(N):
        A = rng.normal(size=(ny, ny)) * scale
        Vs.append(0.5 * (A + A.T))                   # real symmetric (time-reversal invariant)
    return Vs


def toy_S(E, Vs, E_perp, dx):
    g = toy_green(E, Vs, E_perp, dx)
    m = g["modes"]
    amp = sm.fisher_lee(g["G11"], g["GN1"], g["G1N"], g["GNN"], m["gamma"], m["gamma"], m["open"], m["open"])
    return sm.assemble_S(amp["r_L"], amp["t_LR"], amp["t_RL"], amp["r_R"]), amp, g


E_PERP = 0.5 * (np.arange(1, 6) * np.pi) ** 2        # L_y = 1: 4.93, 19.7, 44.4, 79.0, 123.4
DX = 0.05


# ----------------------------------------------------------------------
# lead solutions
# ----------------------------------------------------------------------
def test_lead_dispersion_velocity_and_broadening():
    E = 60.0
    m = sm.lead_modes(E, E_PERP, DX)
    assert list(m["open"]) == [True, True, True, False, False]
    k = m["k"][m["open"]]
    assert np.allclose(E_PERP[m["open"]] + (1 - np.cos(k * DX)) / DX**2, E, atol=1e-10)   # discrete dispersion
    assert np.allclose(m["v"][m["open"]], np.sin(k * DX) / DX)
    assert np.allclose(m["gamma"], m["v"] / DX)                                         # Gamma = v / dx
    assert np.all(m["gamma"][~m["open"]] == 0.0)                                        # closed: no flux
    assert np.all(np.abs(m["lam"][~m["open"]]) < 1) and np.all(np.isreal(m["lam"][~m["open"]]))
    assert np.allclose(m["sigma"], -sm.hopping(DX) * m["lam"])


def test_retarded_root_selection():
    """lambda is the root of lambda + 1/lambda = 2c(E + i eta) with |lambda| < 1, eta -> 0+."""
    E, eta = 60.0, 1e-9
    m = sm.lead_modes(E, E_PERP, DX)
    for n, Ep in enumerate(E_PERP):
        c = 1 - ((E + 1j * eta) - Ep) * DX**2
        roots = np.roots([1.0, -2 * c, 1.0])
        lam = roots[np.argmin(np.abs(roots))]
        assert abs(lam - m["lam"][n]) < 1e-6
    assert np.all(np.imag(m["sigma"][m["open"]]) < 0)    # retarded: Im Sigma < 0


def test_continuum_limit_of_velocity():
    E = 60.0
    m = sm.lead_modes(E, E_PERP, 1e-4)
    vc = np.sqrt(2 * (E - E_PERP[m["open"]]))
    assert np.allclose(m["v"][m["open"]], vc, rtol=1e-6)


def test_threshold_energy_rejected():
    with pytest.raises(ValueError):
        sm.lead_modes(E_PERP[2], E_PERP, DX)


# ----------------------------------------------------------------------
# 1D chain: Fisher-Lee vs independent transfer-matrix wave function
# ----------------------------------------------------------------------
def transfer_matrix_1d(E, V, dx):
    """Scattering of a 1D chain (single channel, E_perp = 0) WITHOUT Green
    functions. Convention P: outgoing amplitude d = c_N = 1, c_{N+1} = lambda.
    Backward recursion through the device, then fit the left-lead solution
    c_j = a lambda^{j-1} + b lambda^{-(j-1)} to (c_0, c_1)."""
    t = sm.hopping(dx)
    lam = sm.lead_modes(E, [0.0], dx)["lam"][0]
    N = len(V)
    c = np.zeros(N + 2, dtype=complex)                # index j = 0..N+1
    c[N], c[N + 1] = 1.0, lam
    for j in range(N, 0, -1):
        eps = 2 * t + V[j - 1]
        c[j - 1] = (eps - E) / t * c[j] - c[j + 1]
    # c_1 = a + b ; c_0 = a / lam + b lam
    A = np.array([[1.0, 1.0], [1.0 / lam, lam]])
    a, b = np.linalg.solve(A, np.array([c[1], c[0]]))
    return 1.0 / a, b / a                             # t, r (same lead velocity: flux = amplitude ratio)


@pytest.mark.parametrize("E", [5.0, 60.0, 300.0])
def test_fisher_lee_1d_matches_transfer_matrix_including_phase(E):
    rng = np.random.default_rng(3)
    V = rng.normal(size=20) * 30.0
    Vs = [np.array([[v]]) for v in V]
    S, amp, _ = toy_S(E, Vs, np.array([0.0]), DX)
    t_tm, r_tm = transfer_matrix_1d(E, V, DX)
    assert abs(amp["t_RL"][0, 0] - t_tm) < 1e-10
    assert abs(amp["r_L"][0, 0] - r_tm) < 1e-10


def test_lattice_continuum_limit_rectangular_barrier():
    """1D lattice with a sharp barrier (sites inside [0, w)) converges to the
    analytic rectangular-barrier transmission as dx -> 0. This checks every
    dx, hopping and velocity factor of the Fisher-Lee normalization."""
    V0, w = 250.0, 0.1
    errs = []
    for eps in (150.0, 320.0):
        e_list = []
        for dx in (1 / 256, 1 / 512, 1 / 1024):
            nb = int(round(w / dx))
            V = np.concatenate([np.zeros(4), np.full(nb, V0), np.zeros(4)])
            S, amp, _ = toy_S(eps, [np.array([[v]]) for v in V], np.array([0.0]), dx)
            e_list.append(abs(abs(amp["t_RL"][0, 0]) ** 2 - rectangular_barrier_transmission(np.array([eps]), V0, w)[0]))
        assert e_list[2] < e_list[0]                  # converging
        assert e_list[2] < 1e-2
        errs.append(e_list)


def test_rectangular_barrier_formula_vs_matching():
    """Analytic formula vs exact matching of psi and psi' at x = 0 and x = w
    (continuum, piecewise constant potential)."""
    V0, w = 250.0, 0.1
    for eps in np.linspace(5.0, 600.0, 97):
        if abs(eps - V0) < 1e-6:
            continue
        k = np.sqrt(2 * eps + 0j)
        q = np.sqrt(2 * (eps - V0) + 0j)
        # unknowns: r, A, B, t ; psi = e^{ikx} + r e^{-ikx} (x<0), A e^{iqx} + B e^{-iqx} (0<x<w), t e^{ikx} (x>w)
        M = np.array([
            [-1, 1, 1, 0],
            [1j * k, 1j * q, -1j * q, 0],
            [0, np.exp(1j * q * w), np.exp(-1j * q * w), -np.exp(1j * k * w)],
            [0, 1j * q * np.exp(1j * q * w), -1j * q * np.exp(-1j * q * w), -1j * k * np.exp(1j * k * w)],
        ])
        rhs = np.array([1, 1j * k, 0, 0])
        r, A, B, t = np.linalg.solve(M, rhs)
        assert abs(abs(t) ** 2 + abs(r) ** 2 - 1) < 1e-10
        assert abs(abs(t) ** 2 - rectangular_barrier_transmission(np.array([eps]), V0, w)[0]) < 1e-10


# ----------------------------------------------------------------------
# multichannel: unitarity, reciprocity, Caroli
# ----------------------------------------------------------------------
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_multichannel_unitarity_reciprocity_caroli(seed):
    rng = np.random.default_rng(seed)
    Vs = random_device(rng, N=12, ny=5, scale=20.0)
    S, amp, g = toy_S(60.0, Vs, E_PERP, DX)
    assert S.shape == (6, 6)                          # 3 open channels per lead
    assert sm.unitarity_residual(S) < 1e-12
    assert sm.reciprocity_residual(S) < 1e-12       # S = S^T
    assert np.max(np.abs(amp["r_L"] - amp["r_L"].T)) < 1e-12
    assert np.max(np.abs(amp["r_R"] - amp["r_R"].T)) < 1e-12
    assert np.max(np.abs(amp["t_LR"] - amp["t_RL"].T)) < 1e-12
    assert sm.mode_reciprocity_residual(amp["t_RL"], amp["t_LR"]) < 1e-12
    T_caroli = sm.caroli_transmission(g["GN1"], g["modes"]["gamma"], g["modes"]["gamma"])
    assert abs(T_caroli - np.sum(np.abs(amp["t_RL"]) ** 2)) < 1e-12
    # total transmission is the same for both incidence directions (follows from unitarity)
    assert abs(np.sum(np.abs(amp["t_RL"]) ** 2) - np.sum(np.abs(amp["t_LR"]) ** 2)) < 1e-12
    # per incoming left channel: transmitted + reflected = 1
    assert np.allclose(np.sum(np.abs(amp["t_RL"]) ** 2, 0) + np.sum(np.abs(amp["r_L"]) ** 2, 0), 1, atol=1e-12)


def test_same_direction_symmetry_is_not_a_valid_criterion():
    """The old B3 criterion |T_{m<-n} - T_{n<-m}| (same direction) is NOT
    implied by time-reversal symmetry. For a generic (left-right
    asymmetric) device it is violated, while the correct opposite-direction
    relation holds to rounding."""
    rng = np.random.default_rng(7)
    Vs = random_device(rng, N=12, ny=5, scale=20.0)
    _, amp, _ = toy_S(60.0, Vs, E_PERP, DX)
    TLR_same = np.abs(amp["t_RL"]) ** 2
    assert np.max(np.abs(TLR_same - TLR_same.T)) > 1e-3
    assert sm.mode_reciprocity_residual(amp["t_RL"], amp["t_LR"]) < 1e-12


def test_y_uniform_potential_decouples_channels():
    """A potential proportional to the identity in mode space (uniform in y)
    gives a diagonal transmission matrix (no channel mixing)."""
    Vs = [np.eye(5) * v for v in np.linspace(0, 40, 10)]
    _, amp, _ = toy_S(60.0, Vs, E_PERP, DX)
    off = amp["t_RL"] - np.diag(np.diag(amp["t_RL"]))
    assert np.max(np.abs(off)) < 1e-12


def test_clean_guide_unit_transmission():
    N = 10
    Vs = [np.zeros((5, 5)) for _ in range(N)]
    S, amp, g = toy_S(60.0, Vs, E_PERP, DX)
    assert np.max(np.abs(np.abs(amp["t_RL"]) - np.eye(3))) < 1e-12
    assert np.max(np.abs(amp["r_L"])) < 1e-12 and np.max(np.abs(amp["r_R"])) < 1e-12
    assert sm.unitarity_residual(S) < 1e-12
    # convention P: free propagation from slice 1 to slice N gives t_nn = lambda_n^(N-1)
    lam = g["modes"]["lam"][g["modes"]["open"]]
    assert np.max(np.abs(np.diag(amp["t_RL"]) - lam ** (N - 1))) < 1e-12
    assert np.max(np.abs(np.diag(amp["t_LR"]) - lam ** (N - 1))) < 1e-12


def test_continuum_velocities_break_unitarity():
    """Using continuum velocities sqrt(2(E - E_perp)) instead of the discrete
    lattice velocities in the flux normalization violates unitarity at
    O((k dx)^2). The discrete Gamma = v_disc/dx is required."""
    rng = np.random.default_rng(0)
    Vs = random_device(rng, N=12, ny=5, scale=20.0)
    g = toy_green(60.0, Vs, E_PERP, DX)
    m = g["modes"]
    gam_cont = np.where(m["open"], np.sqrt(np.maximum(2 * (60.0 - E_PERP), 0)) / DX, 0.0)
    amp = sm.fisher_lee(g["G11"], g["GN1"], g["G1N"], g["GNN"], gam_cont, gam_cont, m["open"], m["open"])
    S = sm.assemble_S(amp["r_L"], amp["t_LR"], amp["t_RL"], amp["r_R"])
    assert sm.unitarity_residual(S) > 1e-3


def test_phase_convention_invariance():
    """Moving reference planes multiplies amplitudes by channel phases.
    Same phase for incoming and outgoing waves of a channel (S' = D S D):
    unitarity and S = S^T preserved. Different phases (S'' = D1 S D2):
    probabilities and unitarity preserved, but S'' != S''^T. So the
    reciprocity test requires a consistent channel phase convention."""
    rng = np.random.default_rng(1)
    S, _, _ = toy_S(60.0, random_device(rng, 12, 5, 20.0), E_PERP, DX)
    D = np.diag(np.exp(1j * rng.uniform(0, 2 * np.pi, 6)))
    D2 = np.diag(np.exp(1j * rng.uniform(0, 2 * np.pi, 6)))
    S1, S2 = D @ S @ D, D @ S @ D2
    for X in (S1, S2):
        assert sm.unitarity_residual(X) < 1e-12
        assert np.allclose(np.abs(X), np.abs(S))
    assert sm.reciprocity_residual(S1) < 1e-12
    assert sm.reciprocity_residual(S2) > 1e-3


# ----------------------------------------------------------------------
# block ordering and index mapping
# ----------------------------------------------------------------------
def test_block_ordering_and_index_mapping():
    nL, nR = 3, 2
    r_L = np.full((nL, nL), 1.0)
    t_LR = np.full((nL, nR), 2.0)
    t_RL = np.arange(nR * nL, dtype=float).reshape(nR, nL) + 10
    r_R = np.full((nR, nR), 4.0)
    S = sm.assemble_S(r_L, t_LR, t_RL, r_R)
    # incidence from left channel n -> outgoing right channel m sits at S[nL + m, n]
    for m in range(nR):
        for n in range(nL):
            assert S[nL + m, n] == t_RL[m, n]
    assert np.all(S[:nL, :nL] == 1) and np.all(S[:nL, nL:] == 2) and np.all(S[nL:, nL:] == 4)
    back = sm.split_S(S, nL)
    for k, v in dict(r_L=r_L, t_LR=t_LR, t_RL=t_RL, r_R=r_R).items():
        assert np.array_equal(back[k], v)
    with pytest.raises(ValueError):
        sm.assemble_S(r_L, t_RL, t_LR, r_R)          # swapped blocks have the wrong shapes


def test_mode_reciprocity_index_mapping():
    rng = np.random.default_rng(5)
    t_RL = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    assert sm.mode_reciprocity_residual(t_RL, t_RL.T) == 0.0    # correct pairing: T^{L->R}_{m<-n} = T^{R->L}_{n<-m}
    assert sm.mode_reciprocity_residual(t_RL, t_RL) > 1e-3      # wrong pairing (same direction) is detected
