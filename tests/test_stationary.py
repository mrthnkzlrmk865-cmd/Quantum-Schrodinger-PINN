"""Stationary RGF solver (stationary.py): construction, recursion, and the
B1 clean-waveguide benchmark of docs/PHASE2_METHOD_SELECTION.md §3.

B1 thresholds are those of the specification, unchanged:
  |T_nn - 1| <= 1e-12, sum_{m != n} T_mn <= 1e-12, R_n <= 1e-12,
  U = max|S^dagger S - I| (and |S S^dagger - I|) <= 1e-10,
  phase |t_nn - lambda_n^(N-1)| <= 1e-10 (provisional in the spec).
Reciprocity has no B1 threshold in the spec; the existing B2/B3 value
Q = max|S - S^T| <= 1e-8 is used (not a new or relaxed threshold).

The random-toy comparisons with a dense inverse are CODING checks of the
recursion, not the B3 benchmark (which needs a Potential disorder
realization and is not implemented here).
"""
import numpy as np
import pytest

import smatrix_spec as sm
from config import waveguide_config
from potential import Potential
from stationary import (Device, build_device, device_hamiltonian_blocks, potential_blocks, rgf_endpoint_blocks,
                        scattering_matrix, transverse_basis)

DX = 1 / 512
X0, NSTEP = 1.7, 717           # device x in [1.7, 1.7 + 717 dx]: slab 1.9-2.9 plus >= 0.2 buffers, N = 718


# ----------------------------------------------------------------------
# basis and potential blocks
# ----------------------------------------------------------------------
def test_transverse_energies_and_quadrature_orthonormality():
    for Ly, ny in [(1.0, 63), (1.3, 31)]:
        b = transverse_basis(Ly, ny)
        assert b.n_quad == 2 * ny + 1
        assert np.allclose(b.E_perp, 0.5 * (np.arange(1, ny + 1) * np.pi / Ly) ** 2, rtol=1e-15)
        assert np.max(np.abs(b.X.T @ b.X * b.weight - np.eye(ny))) < 1e-13


@pytest.mark.parametrize("kw", [dict(Ly=0.0, n_modes=5), dict(Ly=1.0, n_modes=0), dict(Ly=1.0, n_modes=5, n_quad=9),
                                dict(Ly=1.0, n_modes=2.5)])
def test_basis_input_validation(kw):
    with pytest.raises(ValueError):
        transverse_basis(**kw)


def test_potential_blocks_symmetric_and_y_uniform_diagonal():
    b = transverse_basis(1.0, 31)
    rng = np.random.default_rng(0)
    B = potential_blocks(rng.normal(size=(7, b.n_quad)) * 50, b)
    assert np.array_equal(B, np.transpose(B, (0, 2, 1)))               # exactly symmetric
    v = np.linspace(-3, 5, 7)
    U = potential_blocks(np.repeat(v[:, None], b.n_quad, axis=1), b)    # V independent of y
    for j in range(7):
        assert np.max(np.abs(U[j] - v[j] * np.eye(31))) <= 1e-13 * max(1.0, abs(v[j]))


def test_potential_matrix_element_analytic():
    """V(y) = cos(2 pi y / L_y): (V)_nm = 1/2 [|n - m| = 2] - 1/2 [n = m = 1]
    (from sin a sin b cos c identities); the DST-I quadrature is exact here."""
    Ly, ny = 1.0, 20
    b = transverse_basis(Ly, ny)
    B = potential_blocks(np.cos(2 * np.pi * b.y / Ly)[None, :], b)[0]
    n = np.arange(1, ny + 1)
    ref = 0.5 * (np.abs(n[:, None] - n[None, :]) == 2) - 0.5 * ((n[:, None] == 1) & (n[None, :] == 1))
    assert np.max(np.abs(B - ref)) < 1e-14


def test_build_device_samples_the_potential_class():
    cfg = waveguide_config(disorder_strength=100.0)
    pot = Potential(cfg)
    dev = build_device(pot.numpy, X0, X0 + NSTEP * DX, DX, 1.0, 31)
    assert dev.n_slices == NSTEP + 1 and np.isclose(dev.xs[-1] - dev.xs[0], NSTEP * DX)
    Xg, Yg = np.meshgrid(dev.xs, dev.basis.y, indexing="ij")
    assert np.array_equal(dev.V, potential_blocks(pot.numpy(Xg, Yg), dev.basis))   # same samples, same blocks
    assert dev.V_lead_max < 1e-6


def test_build_device_input_validation():
    pot = Potential(waveguide_config(disorder_strength=100.0))
    with pytest.raises(ValueError, match="integer"):
        build_device(pot.numpy, 1.7, 3.1, DX, 1.0, 15)                  # 1.4/dx = 716.8
    with pytest.raises(ValueError, match="clean-lead"):
        build_device(pot.numpy, 2.0, 2.5, DX, 1.0, 15)                  # inside the slab: leads not clean
    with pytest.raises(ValueError):
        build_device(pot.numpy, 1.7, 1.7 + 10 * DX, -DX, 1.0, 15)


def test_hamiltonian_blocks():
    b = transverse_basis(1.0, 5)
    V = np.zeros((3, 5, 5))
    V[1] = np.diag(np.arange(5.0))
    dev = Device(xs=np.arange(3) * 0.1, dx=0.1, basis=b, V=V, V_lead_max=0.0)
    H, t = device_hamiltonian_blocks(dev)
    assert t == sm.hopping(0.1) and np.isclose(t, 50.0, rtol=1e-14)    # 1/(2 dx^2); 0.1**2 is inexact
    assert np.allclose(H[0], 2 * t * np.eye(5) + np.diag(b.E_perp))
    assert np.allclose(H[1], 2 * t * np.eye(5) + np.diag(b.E_perp) + V[1])


# ----------------------------------------------------------------------
# recursion (coding checks on random toys; NOT the B3 benchmark)
# ----------------------------------------------------------------------
def _toy(rng, N=15, ny=6, dx=0.05, scale=20.0):
    b = transverse_basis(1.0, ny)
    A = rng.normal(size=(N, ny, ny)) * scale
    return Device(xs=np.arange(N) * dx, dx=dx, basis=b, V=0.5 * (A + np.transpose(A, (0, 2, 1))), V_lead_max=0.0)


def _dense_G(E, dev, sigma):
    H, t = device_hamiltonian_blocks(dev)
    N, ny = dev.n_slices, dev.basis.n_modes
    M = np.zeros((N * ny, N * ny), dtype=complex)
    for j in range(N):
        s = slice(j * ny, (j + 1) * ny)
        M[s, s] = E * np.eye(ny) - H[j]
        if j + 1 < N:
            s2 = slice((j + 1) * ny, (j + 2) * ny)
            M[s, s2] = M[s2, s] = t * np.eye(ny)
    M[:ny, :ny] -= np.diag(sigma)
    M[-ny:, -ny:] -= np.diag(sigma)
    G = np.linalg.inv(M)
    blk = lambda i, j: G[i * ny:(i + 1) * ny, j * ny:(j + 1) * ny]
    return dict(G11=blk(0, 0), GN1=blk(N - 1, 0), G1N=blk(0, N - 1), GNN=blk(N - 1, N - 1))


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_rgf_matches_dense_inverse(seed):
    dev = _toy(np.random.default_rng(seed))
    for E in (30.0, 60.0, 150.0):
        sig = sm.lead_modes(E, dev.basis.E_perp, dev.dx)["sigma"]
        G, D = rgf_endpoint_blocks(E, dev, sig), _dense_G(E, dev, sig)
        for k in ("G11", "GN1", "G1N", "GNN"):
            assert np.max(np.abs(G[k] - D[k])) <= 1e-12 * np.max(np.abs(D[k])), (k, E)
        # G_1N comes from its own recursion; it equals G_N1^T only through M = M^T
        assert np.max(np.abs(G["G1N"] - G["GN1"].T)) <= 1e-12 * np.max(np.abs(G["GN1"]))


def test_closed_channel_self_energy_is_included_and_matters():
    """All N_y channels enter Sigma; closed ones with a real, negative
    Sigma = -t lambda (0 < lambda < 1) and zero broadening. Dropping the
    closed-channel part changes G, so it cannot be omitted."""
    dev = _toy(np.random.default_rng(4))
    E = 60.0
    m = sm.lead_modes(E, dev.basis.E_perp, dev.dx)
    closed = ~m["open"]
    assert np.any(closed)
    assert np.all(np.isreal(m["sigma"][closed])) and np.all(m["sigma"][closed].real < 0)
    assert np.all(m["gamma"][closed] == 0)
    full = rgf_endpoint_blocks(E, dev, m["sigma"])
    assert np.max(np.abs(full["GN1"] - _dense_G(E, dev, m["sigma"])["GN1"])) < 1e-12 * np.max(np.abs(full["GN1"]))
    trunc = rgf_endpoint_blocks(E, dev, np.where(closed, 0.0, m["sigma"]))
    assert np.max(np.abs(full["GN1"] - trunc["GN1"])) > 1e-3 * np.max(np.abs(full["GN1"]))


def test_scattering_matrix_guards():
    dev = _toy(np.random.default_rng(0), ny=3)
    with pytest.raises(ValueError, match="closed channels"):
        scattering_matrix(200.0, dev)                     # all 3 basis modes open
    with pytest.raises(ValueError, match="no open channel"):
        scattering_matrix(1.0, dev)                       # below the first threshold
    with pytest.raises(ValueError, match="threshold"):
        scattering_matrix(float(dev.basis.E_perp[1]), dev)


def test_toy_scattering_matrix_unitary_and_reciprocal():
    dev = _toy(np.random.default_rng(1))
    r = scattering_matrix(60.0, dev)
    assert r["S"].shape == (2 * r["n_open"], 2 * r["n_open"])
    assert sm.unitarity_residual(r["S"]) < 1e-12 and sm.reciprocity_residual(r["S"]) < 1e-12
    T = np.abs(r["t_RL"]) ** 2
    assert np.max(np.abs(T - T.T)) > 1e-3        # generic device: same-direction symmetry does NOT hold
    assert sm.mode_reciprocity_residual(r["t_RL"], r["t_LR"]) < 1e-12
    assert abs(r["caroli"] - np.sum(T)) < 1e-12


# ----------------------------------------------------------------------
# B1: clean, matched waveguide (production geometry)
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def clean_device_63():
    pot = Potential(waveguide_config())            # no disorder, no barrier: V = 0
    return build_device(pot.numpy, X0, X0 + NSTEP * DX, DX, 1.0, 63)


def b1_errors(r, N):
    lam = r["modes"]["lam"][r["modes"]["open"]]
    S = r["S"]
    T, Tb = np.abs(r["t_RL"]) ** 2, np.abs(r["t_LR"]) ** 2
    I = np.eye(S.shape[0])
    return dict(
        T_diag=max(np.max(np.abs(np.diag(T) - 1)), np.max(np.abs(np.diag(Tb) - 1))),
        mixing=max(np.max(T.sum(0) - np.diag(T)), np.max(Tb.sum(0) - np.diag(Tb))),
        R=max(np.max((np.abs(r["r_L"]) ** 2).sum(0)), np.max((np.abs(r["r_R"]) ** 2).sum(0))),
        U=max(np.max(np.abs(S.conj().T @ S - I)), np.max(np.abs(S @ S.conj().T - I))),
        phase=max(np.max(np.abs(np.diag(r["t_RL"]) - lam ** (N - 1))),
                  np.max(np.abs(np.diag(r["t_LR"]) - lam ** (N - 1)))),
        Q=sm.reciprocity_residual(S))


B1_LIMITS = dict(T_diag=1e-12, mixing=1e-12, R=1e-12, U=1e-10, phase=1e-10, Q=1e-8)


def spec_b1_energies(E_perp):
    """50 energies spanning the 1st-6th channel windows, avoiding
    |E - E_perp_n| < 1e-3 E (spec B1)."""
    Es = np.linspace(E_perp[0], E_perp[6], 52)[1:-1]
    keep = [E for E in Es if np.all(np.abs(E - E_perp) >= 1e-3 * E)]
    return np.array(keep)


@pytest.mark.parametrize("E,n_open_expected", [(10.0, 1), (100.0, 4), (230.0, 6)])
def test_b1_selected_energies(clean_device_63, E, n_open_expected):
    """One open channel (E = 10), several (E = 100), and the six channels of
    the transport problem (E = 230); every case also has 57-62 closed
    channels in the basis."""
    r = scattering_matrix(E, clean_device_63)
    assert r["n_open"] == n_open_expected and clean_device_63.basis.n_modes - r["n_open"] >= 57
    err = b1_errors(r, clean_device_63.n_slices)
    for k, lim in B1_LIMITS.items():
        assert err[k] <= lim, (k, err[k])
    assert r["min_rcond"] > 1e-6


def test_b1_closed_channels_match_infinite_chain(clean_device_63):
    """Clean device + exact leads = infinite clean chain, diagonal in modes:
    G_11,nn = 1/(t (lambda - 1/lambda)),  G_N1,nn = lambda^(N-1) G_11,nn,
    for open AND closed channels (closed: real, decaying)."""
    dev, E = clean_device_63, 100.0
    r = scattering_matrix(E, dev)
    lam, t, N = r["modes"]["lam"], dev.t_x, dev.n_slices
    g0 = 1.0 / (t * (lam - 1.0 / lam))
    for key, ref in (("G11", g0), ("GN1", lam ** (N - 1) * g0)):
        G = r["G"][key]
        off = G - np.diag(np.diag(G))
        assert np.max(np.abs(off)) <= 1e-12 * np.max(np.abs(ref))
        closed = ~r["modes"]["open"]
        # closed channels with a non-negligible value: relative check
        sig = closed & (np.abs(ref) > 1e-200)
        assert np.max(np.abs(np.diag(G)[sig] / ref[sig] - 1)) <= 1e-9
        assert np.max(np.abs(np.diag(G)[~closed] - ref[~closed])) <= 1e-10 * np.max(np.abs(ref[~closed]))


def test_b1_spec_energy_grid(clean_device_63):
    Es = spec_b1_energies(clean_device_63.basis.E_perp)
    assert len(Es) == 50
    worst = {k: 0.0 for k in B1_LIMITS}
    for E in Es:
        err = b1_errors(scattering_matrix(E, clean_device_63), clean_device_63.n_slices)
        for k in worst:
            worst[k] = max(worst[k], err[k])
    for k, lim in B1_LIMITS.items():
        assert worst[k] <= lim, (k, worst[k])


def test_b1_independent_of_number_of_transverse_modes():
    """Clean guide: the open-channel S matrix must not depend on how many
    closed channels are kept (N_y = 15 ... 127)."""
    pot = Potential(waveguide_config())
    E = 100.0
    S = {}
    for ny in (15, 63, 127):
        dev = build_device(pot.numpy, X0, X0 + NSTEP * DX, DX, 1.0, ny)
        S[ny] = scattering_matrix(E, dev)["S"]
    assert np.max(np.abs(S[15] - S[63])) <= 1e-12 and np.max(np.abs(S[127] - S[63])) <= 1e-12


# ----------------------------------------------------------------------
# Longest device (dx = 1/1024, L_d = 3, N = 3483): measured behaviour
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def long_device():
    pot = Potential(waveguide_config(slab_length=3.0))
    dx = 1 / 1024
    n = int(np.ceil(3.4 / dx - 1e-9))
    return build_device(pot.numpy, 1.7, 1.7 + n * dx, dx, 1.0, 63)


def _long_E(dev):
    """Index 25 of the spec B1 grid (E = 125.692...): the worst energy for the
    phase on this device (channel 5 close to its threshold), taken from the
    grid itself rather than typed digits."""
    return spec_b1_energies(dev.basis.E_perp)[25]


def test_long_device_magnitudes_unitarity_reciprocity(long_device):
    r = scattering_matrix(_long_E(long_device), long_device)
    err = b1_errors(r, long_device.n_slices)
    for k in ("T_diag", "mixing", "R", "U", "Q"):
        assert err[k] <= B1_LIMITS[k], (k, err[k])


@pytest.mark.xfail(strict=True, reason="KNOWN, MEASURED: the provisional B1 phase criterion (1e-10) fails on the "
                   "3483-slice device (1.77e-10). Cause: FD conditioning near thresholds, error ~ N eps/sin(kdx); "
                   "see docs/RGF_B1_VALIDATION.md. Threshold intentionally NOT relaxed.")
def test_long_device_provisional_phase_criterion(long_device):
    r = scattering_matrix(_long_E(long_device), long_device)
    assert b1_errors(r, long_device.n_slices)["phase"] <= B1_LIMITS["phase"]
