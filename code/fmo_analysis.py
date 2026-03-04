#!/usr/bin/env python3
"""
FMO complex analysis for:
"Information Geometry of the Quantum-Classical Transition in Photosynthetic Exciton Transport"

Computes:
  - Lindblad dephasing + sink evolution of the FMO density matrix
  - Transport efficiency (ENAQT) via reaction-centre trapping
  - QFIM and effective dimensionality vs dephasing rate
  - Principal angles between diagonal/off-diagonal tangent subspaces

Usage:
  cd physics/70_decoherence_dimensional_collapse
  python3 code/fmo_analysis.py
"""

import numpy as np
from scipy.linalg import expm

# ── FMO Hamiltonian (Adolphs & Renger 2006, cm^-1 relative to 12210) ────────

H_FMO = np.array([
    [200,   -87.7,   5.5,  -5.9,   6.7, -13.7,  -9.9],
    [-87.7,  320,   30.8,   8.2,   0.7,  11.8,   4.3],
    [  5.5,  30.8,    0,  -53.5,  -2.2,  -9.6,   6.0],
    [ -5.9,   8.2, -53.5,  110,  -70.7, -17.0, -63.3],
    [  6.7,   0.7,  -2.2, -70.7,  270,   81.1,  -1.3],
    [-13.7,  11.8,  -9.6, -17.0,  81.1,  420,   39.7],
    [ -9.9,   4.3,   6.0, -63.3,  -1.3,  39.7,  230],
], dtype=np.float64)

D_FMO = 7  # Number of sites
COUPLING_THRESHOLD = 5.0  # cm^-1, below which coupling is "negligible"

# Convert Hamiltonian to angular frequency (rad/fs)
# 1 cm^-1 = 2*pi*c = 2*pi*3e10 cm/s = 1.884e11 rad/s = 1.884e-4 rad/fs
CM_TO_RADFS = 2.0 * np.pi * 2.998e-5  # cm^-1 to rad/fs


def fmo_coupling_graph(threshold=COUPLING_THRESHOLD):
    """Return adjacency matrix of significant couplings."""
    adj = np.zeros((D_FMO, D_FMO), dtype=bool)
    for i in range(D_FMO):
        for j in range(i + 1, D_FMO):
            if abs(H_FMO[i, j]) > threshold:
                adj[i, j] = True
                adj[j, i] = True
    return adj


def weak_coupling_pairs(threshold=COUPLING_THRESHOLD):
    """Find site pairs with coupling below threshold in the Hamiltonian."""
    pairs = []
    for i in range(D_FMO):
        for j in range(i + 1, D_FMO):
            if abs(H_FMO[i, j]) < threshold:
                pairs.append((i + 1, j + 1, H_FMO[i, j]))  # 1-indexed
    return pairs


def lindblad_dephasing_rhs(rho, H, gamma):
    """Right-hand side of Lindblad master equation with Haken-Strobl dephasing.

    drho/dt = -i[H, rho] + gamma * sum_k (|k><k| rho |k><k| - 0.5 {|k><k|, rho})
    """
    d = H.shape[0]
    # Unitary part
    commutator = -1j * (H @ rho - rho @ H)

    # Dephasing part: for site-basis dephasing,
    # sum_k |k><k| rho |k><k| = diag(rho)
    # sum_k {|k><k|, rho} = 2 * diag(diag(rho)) + off-diag terms...
    # Actually: |k><k| rho |k><k| = rho_kk |k><k|
    # {|k><k|, rho}_mn = delta_mk rho_kn + rho_mk delta_kn
    # So sum_k {|k><k|, rho}_mn = rho_mn (delta_mm + delta_nn) for m!=n -> rho_mn
    # and for m=n: sum_k {|k><k|, rho}_mm = 2*rho_mm
    # Net effect: off-diagonal elements decay at rate gamma, diagonal preserved.
    dephasing = np.zeros_like(rho)
    for k in range(d):
        proj = np.zeros((d, d), dtype=complex)
        proj[k, k] = 1.0
        dephasing += proj @ rho @ proj - 0.5 * (proj @ rho + rho @ proj)

    return commutator + gamma * dephasing


def evolve_lindblad(H, gamma, rho0, t_final, dt=0.1):
    """Evolve density matrix under Lindblad dephasing using RK4.

    Parameters:
        H: Hamiltonian in rad/fs
        gamma: dephasing rate in rad/fs
        rho0: initial density matrix
        t_final: evolution time in fs
        dt: time step in fs

    Returns:
        rho: final density matrix
    """
    rho = rho0.copy().astype(complex)
    t = 0.0
    n_steps = int(t_final / dt)

    for _ in range(n_steps):
        k1 = dt * lindblad_dephasing_rhs(rho, H, gamma)
        k2 = dt * lindblad_dephasing_rhs(rho + 0.5 * k1, H, gamma)
        k3 = dt * lindblad_dephasing_rhs(rho + 0.5 * k2, H, gamma)
        k4 = dt * lindblad_dephasing_rhs(rho + k3, H, gamma)
        rho = rho + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

        # Enforce Hermiticity and trace
        rho = 0.5 * (rho + rho.conj().T)
        rho /= np.trace(rho)

    return rho


def lindblad_dephasing_sink_rhs(rho, H, gamma, kappa, sink_site=2):
    """Lindblad RHS with Haken-Strobl dephasing + trapping sink.

    The sink at site `sink_site` (0-indexed, default 2 = site 3) is modeled as
    population loss: L_sink[ρ] = -κ/2 (P_s ρ + ρ P_s), where P_s = |s><s|.
    This makes the evolution non-trace-preserving; lost trace = trapped population.

    Parameters:
        rho: density matrix
        H: Hamiltonian in rad/fs
        gamma: dephasing rate in rad/fs
        kappa: trapping rate in rad/fs
        sink_site: 0-indexed site for the reaction-centre trap
    """
    # Start with dephasing dynamics
    drho = lindblad_dephasing_rhs(rho, H, gamma)

    # Add sink: -κ/2 (|s><s| ρ + ρ |s><s|)
    d = H.shape[0]
    proj_s = np.zeros((d, d), dtype=complex)
    proj_s[sink_site, sink_site] = 1.0
    drho -= 0.5 * kappa * (proj_s @ rho + rho @ proj_s)

    return drho


def compute_transport_efficiency(gamma_cm, kappa_trap_cm=5.3, t_final_fs=20000.0,
                                  dt_fs=1.0, initial_site=0, sink_site=2):
    """Compute FMO transport efficiency with Haken-Strobl dephasing + sink.

    Parameters:
        gamma_cm: dephasing rate in cm^-1
        kappa_trap_cm: trapping rate in cm^-1 (default 5.3 ≈ 1 ps^-1)
        t_final_fs: evolution time in fs
        dt_fs: time step in fs
        initial_site: 0-indexed initial excitation site
        sink_site: 0-indexed reaction-centre site

    Returns:
        efficiency: fraction of population trapped (0 to 1)
    """
    H = H_FMO * CM_TO_RADFS
    gamma = gamma_cm * CM_TO_RADFS
    kappa = kappa_trap_cm * CM_TO_RADFS
    d = D_FMO

    rho = np.zeros((d, d), dtype=complex)
    rho[initial_site, initial_site] = 1.0

    n_steps = int(t_final_fs / dt_fs)
    trapped = 0.0  # accumulated trapped population

    for _ in range(n_steps):
        # Accumulate trapping: η = κ ∫ p_sink(t) dt
        trapped += kappa * np.real(rho[sink_site, sink_site]) * dt_fs

        k1 = dt_fs * lindblad_dephasing_sink_rhs(rho, H, gamma, kappa, sink_site)
        k2 = dt_fs * lindblad_dephasing_sink_rhs(rho + 0.5 * k1, H, gamma, kappa, sink_site)
        k3 = dt_fs * lindblad_dephasing_sink_rhs(rho + 0.5 * k2, H, gamma, kappa, sink_site)
        k4 = dt_fs * lindblad_dephasing_sink_rhs(rho + k3, H, gamma, kappa, sink_site)
        rho = rho + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

        # Enforce Hermiticity (don't normalize trace — it should decrease)
        rho = 0.5 * (rho + rho.conj().T)

    return min(trapped, 1.0)


def scan_efficiency(gamma_values_cm, initial_site=0, **kwargs):
    """Scan dephasing rates and compute transport efficiency.

    Returns array of efficiencies, one per gamma value.
    """
    efficiencies = []
    for gamma_cm in gamma_values_cm:
        eta = compute_transport_efficiency(gamma_cm, initial_site=initial_site, **kwargs)
        efficiencies.append(eta)
        print(f"  gamma={gamma_cm:7.1f} cm^-1  eta={eta:.4f}  (site {initial_site+1})")
    return np.array(efficiencies)


def evolve_lindblad_sink(H, gamma, kappa, rho0, t_final, dt=1.0, sink_site=2):
    """Evolve density matrix under Lindblad dephasing + sink using RK4.

    Unlike evolve_lindblad, this does NOT renormalize trace during evolution,
    since the sink makes the dynamics non-trace-preserving.

    Parameters:
        H: Hamiltonian in rad/fs
        gamma: dephasing rate in rad/fs
        kappa: trapping rate in rad/fs
        rho0: initial density matrix
        t_final: evolution time in fs
        dt: time step in fs
        sink_site: 0-indexed reaction-centre site

    Returns:
        rho: final (sub-normalized) density matrix
    """
    rho = rho0.copy().astype(complex)
    n_steps = int(t_final / dt)

    for _ in range(n_steps):
        k1 = dt * lindblad_dephasing_sink_rhs(rho, H, gamma, kappa, sink_site)
        k2 = dt * lindblad_dephasing_sink_rhs(rho + 0.5 * k1, H, gamma, kappa, sink_site)
        k3 = dt * lindblad_dephasing_sink_rhs(rho + 0.5 * k2, H, gamma, kappa, sink_site)
        k4 = dt * lindblad_dephasing_sink_rhs(rho + k3, H, gamma, kappa, sink_site)
        rho = rho + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        rho = 0.5 * (rho + rho.conj().T)  # Enforce Hermiticity only

    return rho


def compute_qfim(rho):
    """Compute the QFIM for a full-rank density matrix in the pointer basis.

    Returns the (d^2-1) x (d^2-1) QFIM in the standard traceless Hermitian basis.
    """
    d = rho.shape[0]
    eigvals, U = np.linalg.eigh(rho)

    # Regularize small eigenvalues
    eigvals = np.maximum(eigvals, 1e-14)

    # Build traceless Hermitian basis
    # First d-1: diagonal (population) directions
    # Then d(d-1): off-diagonal (coherence) directions
    basis = []

    # Diagonal basis: e_a = diag with +1 at a, -1 at a+1
    for a in range(d - 1):
        D = np.zeros((d, d), dtype=complex)
        D[a, a] = 1.0
        D[a + 1, a + 1] = -1.0
        basis.append(D)

    # Off-diagonal basis: real and imaginary parts
    for k in range(d):
        for l in range(k + 1, d):
            R = np.zeros((d, d), dtype=complex)
            R[k, l] = 1.0
            R[l, k] = 1.0
            basis.append(R)

            I_mat = np.zeros((d, d), dtype=complex)
            I_mat[k, l] = 1j
            I_mat[l, k] = -1j
            basis.append(I_mat)

    n_params = len(basis)
    F_Q = np.zeros((n_params, n_params))

    # Compute SLDs in eigenbasis
    def sld(X):
        X_eig = U.conj().T @ X @ U
        L_eig = np.zeros((d, d), dtype=complex)
        for m in range(d):
            for n in range(d):
                denom = eigvals[m] + eigvals[n]
                if denom > 1e-14:
                    L_eig[m, n] = 2.0 * X_eig[m, n] / denom
        return U @ L_eig @ U.conj().T

    slds = [sld(b) for b in basis]

    for i in range(n_params):
        for j in range(i, n_params):
            val = 0.5 * np.real(np.trace(rho @ (slds[i] @ slds[j] + slds[j] @ slds[i])))
            F_Q[i, j] = val
            F_Q[j, i] = val

    return F_Q


def effective_dimensionality(F_Q):
    """Compute participation ratio from QFIM eigenvalues."""
    eigvals = np.linalg.eigvalsh(F_Q)
    eigvals = eigvals[eigvals > 1e-12]
    if len(eigvals) == 0:
        return 0.0
    return (eigvals.sum()) ** 2 / (eigvals ** 2).sum()


def principal_angles(rho):
    """Compute the d-1 principal angles between T^diag and T^off at state rho."""
    d = rho.shape[0]
    eigvals, U = np.linalg.eigh(rho)
    eigvals = np.maximum(eigvals, 1e-14)

    # Diagonal basis
    diag_basis = []
    for a in range(d - 1):
        D = np.zeros((d, d), dtype=complex)
        D[a, a] = 1.0
        D[a + 1, a + 1] = -1.0
        diag_basis.append(D)

    # Off-diagonal basis
    off_basis = []
    for k in range(d):
        for l in range(k + 1, d):
            R = np.zeros((d, d), dtype=complex)
            R[k, l] = 1.0
            R[l, k] = 1.0
            off_basis.append(R)
            I_mat = np.zeros((d, d), dtype=complex)
            I_mat[k, l] = 1j
            I_mat[l, k] = -1j
            off_basis.append(I_mat)

    def sld(X):
        X_eig = U.conj().T @ X @ U
        L_eig = np.zeros((d, d), dtype=complex)
        for m in range(d):
            for n in range(d):
                denom = eigvals[m] + eigvals[n]
                if denom > 1e-14:
                    L_eig[m, n] = 2.0 * X_eig[m, n] / denom
        return U @ L_eig @ U.conj().T

    def bures_inner(X, Y):
        LX = sld(X)
        LY = sld(Y)
        return 0.25 * np.real(np.trace(rho @ (LX @ LY + LY @ LX))) / 2.0

    n_diag = len(diag_basis)
    n_off = len(off_basis)

    G_DD = np.zeros((n_diag, n_diag))
    for i in range(n_diag):
        for j in range(i, n_diag):
            v = bures_inner(diag_basis[i], diag_basis[j])
            G_DD[i, j] = v
            G_DD[j, i] = v

    G_OO = np.zeros((n_off, n_off))
    for i in range(n_off):
        for j in range(i, n_off):
            v = bures_inner(off_basis[i], off_basis[j])
            G_OO[i, j] = v
            G_OO[j, i] = v

    G_DO = np.zeros((n_diag, n_off))
    for i in range(n_diag):
        for j in range(n_off):
            G_DO[i, j] = bures_inner(diag_basis[i], off_basis[j])

    def matrix_sqrt_inv(M):
        ev, V = np.linalg.eigh(M)
        ev = np.maximum(ev, 1e-12)
        return V @ np.diag(1.0 / np.sqrt(ev)) @ V.T

    M = matrix_sqrt_inv(G_DD) @ G_DO @ matrix_sqrt_inv(G_OO)
    svs = np.linalg.svd(M, compute_uv=False)
    cosines = np.clip(svs[:n_diag], 0.0, 1.0)
    return np.arccos(cosines)


def scan_dephasing_rates(gamma_values_cm, t_final_fs=10000.0, dt_fs=0.5):
    """Scan dephasing rates and compute D_eff and min principal angle.

    Parameters:
        gamma_values_cm: array of dephasing rates in cm^-1
        t_final_fs: evolution time in fs
        dt_fs: time step in fs

    Returns:
        dict with 'gamma_cm', 'D_eff', 'min_angle_deg'
    """
    H = H_FMO * CM_TO_RADFS
    rho0 = np.zeros((D_FMO, D_FMO), dtype=complex)
    rho0[0, 0] = 1.0  # Initial excitation at site 1

    D_effs = []
    min_angles = []

    for gamma_cm in gamma_values_cm:
        gamma = gamma_cm * CM_TO_RADFS
        rho = evolve_lindblad(H, gamma, rho0, t_final_fs, dt_fs)

        # Check that rho is sensible
        eigvals = np.linalg.eigvalsh(rho)
        if eigvals.min() < -1e-10:
            print(f"  WARNING: negative eigenvalue {eigvals.min():.2e} at gamma={gamma_cm}")

        # Make rho strictly full-rank for QFIM computation
        rho_reg = rho + 1e-10 * np.eye(D_FMO) / D_FMO
        rho_reg /= np.trace(rho_reg)

        F_Q = compute_qfim(rho_reg)
        D_eff = effective_dimensionality(F_Q)
        D_effs.append(D_eff)

        # Principal angles (expensive for d=7)
        angles = principal_angles(rho_reg)
        min_angles.append(np.degrees(np.min(angles)))

        print(f"  gamma={gamma_cm:8.1f} cm^-1  D_eff={D_eff:6.2f}  "
              f"min_angle={np.degrees(np.min(angles)):5.1f}°")

    return {
        'gamma_cm': np.array(gamma_values_cm),
        'D_eff': np.array(D_effs),
        'min_angle_deg': np.array(min_angles),
    }


if __name__ == "__main__":
    print("FMO Complex Analysis")
    print("=" * 60)

    # 1. Coupling graph
    adj = fmo_coupling_graph()
    print(f"\nWeak coupling pairs (|J| < {COUPLING_THRESHOLD} cm^-1):")
    for i, j, J in weak_coupling_pairs():
        print(f"  ({i},{j}): J = {J:.1f} cm^-1")

    # 2. Transport efficiency scan
    print("\n" + "=" * 60)
    print("Transport Efficiency (site 1 → site 3 sink)")
    print("=" * 60)

    gamma_values = np.concatenate([
        np.array([0.1, 0.5, 1.0, 2.0, 5.0]),
        np.arange(10, 110, 10),
        np.arange(150, 600, 50),
    ])
    gamma_values = np.sort(np.unique(gamma_values))

    etas = scan_efficiency(gamma_values, initial_site=0)

    idx_max = np.argmax(etas)
    print(f"\nPeak efficiency: eta={etas[idx_max]:.4f} at gamma={gamma_values[idx_max]:.1f} cm^-1")

    print("\nDone.")
