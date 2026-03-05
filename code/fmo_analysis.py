#!/usr/bin/env python3
"""
Photosynthetic complex analysis for:
"Information Geometry of the Quantum-Classical Transition: From Photosynthetic Excitons to Neural Oscillations"

Computes:
  - Lindblad dephasing + sink evolution of density matrices
  - Transport efficiency via target-manifold trapping
  - QFIM and effective dimensionality vs dephasing rate
  - Principal angles between diagonal/off-diagonal tangent subspaces

Supports multiple complexes via COMPLEXES registry (FMO, PE545, Chromatin, NeuralGamma, NeuralPD).

Usage:
  cd physics/70_decoherence_dimensional_collapse
  python3 code/fmo_analysis.py
"""

import numpy as np
from scipy.linalg import expm

# ── Unit conversion ──────────────────────────────────────────────────────────
# 1 cm^-1 = 2*pi*c = 2*pi*3e10 cm/s = 1.884e11 rad/s = 1.884e-4 rad/fs
CM_TO_RADFS = 2.0 * np.pi * 2.998e-5  # cm^-1 to rad/fs
PS_TO_FS = 1000.0  # 1 ps = 1000 fs

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

# ── PE545 Hamiltonian (Novoderezhkin et al. 2010, Model E, cm^-1) ────────────
# Site ordering: 0=PEB50/61C, 1=DBV19A, 2=DBV19B, 3=PEB82C,
#                4=PEB158C, 5=PEB50/61D, 6=PEB82D, 7=PEB158D

H_PE545 = np.array([
    [18532,    1,  -37,   37,   23,   92,  -16,   12],
    [    1, 18008,   4,  -11,   33,  -39,  -46,    3],
    [  -37,    4, 17973,  45,    3,    2,  -11,   34],
    [   37,  -11,   45, 18040,  -7,  -17,   -3,    6],
    [   23,   33,    3,   -7, 18711,  18,    7,    6],
    [   92,  -39,    2,  -17,   18, 19574,  40,   26],
    [  -16,  -46,  -11,   -3,    7,   40, 19050,   7],
    [   12,    3,   34,    6,    6,   26,    7, 18960],
], dtype=np.float64)

D_FMO = 7  # Number of sites
COUPLING_THRESHOLD = 5.0  # cm^-1, below which coupling is "negligible"

# ── Neural Hamiltonian: Gamma-band cortical microcircuit (d=10) ──────────────
# Synthetic model of a cortical column with 10 coupled oscillator populations
# spanning the gamma band (30–70 Hz). Coupling: ring + 2 long-range shortcuts
# (small-world topology). Site energies and couplings in Hz, then scaled to
# cm^-1 range for numerical compatibility (Lindblad scale invariance: principal
# angles depend on γ/J_max, not absolute energy).
#
# Scaling: alpha = 50 cm^-1 / Hz → entries O(10^3), matching FMO/PE545 range.
# Physical dephasing: γ_bio = kT/ℏ ≈ 4e13 Hz at 310K → γ_bio/J_max ≈ 5e12.

_NEURAL_GAMMA_FREQS_HZ = np.array([30, 33, 36, 40, 44, 48, 53, 58, 64, 70],
                                    dtype=np.float64)
_NEURAL_GAMMA_SCALE = 50.0  # cm^-1 per Hz

def _build_neural_gamma_hamiltonian():
    """Build 10×10 gamma-band cortical microcircuit Hamiltonian (cm^-1)."""
    d = 10
    H = np.diag(_NEURAL_GAMMA_FREQS_HZ * _NEURAL_GAMMA_SCALE)
    # Ring coupling: nearest-neighbor, J = 5 Hz → 250 cm^-1
    for i in range(d):
        j = (i + 1) % d
        H[i, j] = 5.0 * _NEURAL_GAMMA_SCALE
        H[j, i] = 5.0 * _NEURAL_GAMMA_SCALE
    # Long-range shortcuts (small-world): sites 1↔6, 3↔8, J = 2 Hz → 100 cm^-1
    for i, j in [(1, 6), (3, 8)]:
        H[i, j] = 2.0 * _NEURAL_GAMMA_SCALE
        H[j, i] = 2.0 * _NEURAL_GAMMA_SCALE
    return H

H_NEURAL_GAMMA = _build_neural_gamma_hamiltonian()

# Maximum coupling: J_max = 5 Hz (nearest-neighbor ring coupling).
NEURAL_GAMMA_JMAX_HZ = 5.0  # Hz (nearest-neighbor coupling)
NEURAL_GAMMA_JMAX_CM = NEURAL_GAMMA_JMAX_HZ * _NEURAL_GAMMA_SCALE  # cm^-1

# Biological dephasing at 310K: γ_bio = kT/ℏ ≈ 4.06e13 Hz
NEURAL_BIO_DEPHASING_HZ = 4.06e13  # Hz
NEURAL_BIO_GAMMA_OVER_JMAX = NEURAL_BIO_DEPHASING_HZ / NEURAL_GAMMA_JMAX_HZ

# ── Neural Hamiltonian: Potjans-Diesmann cortical column (d=8) ───────────────
# Empirical connectivity from Potjans & Diesmann 2014 (Cerebral Cortex).
# 4 layers × {E, I} = L2/3e, L2/3i, L4e, L4i, L5e, L5i, L6e, L6i.
# Connection probabilities from their Table 5, converted to effective coupling.
# Site energies: representative firing rates (Hz) scaled to cm^-1.

_PD_FREQS_HZ = np.array([3.0, 6.0, 5.0, 8.0, 7.0, 9.0, 1.0, 4.0],
                          dtype=np.float64)  # Approximate mean rates

def _build_neural_pd_hamiltonian():
    """Build 8×8 Potjans-Diesmann cortical column Hamiltonian (cm^-1).

    Connection probabilities from Table 5 of Potjans & Diesmann 2014.
    Converted to symmetric effective couplings: J_ij = sqrt(C_ij * C_ji) * w,
    where w is a scale factor. Excitatory-to-excitatory couplings are positive,
    involving inhibitory populations are negative (sign convention).
    """
    d = 8
    # Connection probabilities (from → to), Table 5 of Potjans & Diesmann 2014
    # Rows/cols: L2/3e, L2/3i, L4e, L4i, L5e, L5i, L6e, L6i
    C = np.array([
        [0.101, 0.169, 0.044, 0.082, 0.032, 0.0,   0.008, 0.0  ],
        [0.135, 0.137, 0.032, 0.052, 0.075, 0.0,   0.004, 0.0  ],
        [0.008, 0.006, 0.050, 0.135, 0.007, 0.0003,0.045, 0.0  ],
        [0.069, 0.003, 0.079, 0.160, 0.003, 0.0,   0.106, 0.0  ],
        [0.100, 0.062, 0.051, 0.006, 0.083, 0.373, 0.020, 0.0  ],
        [0.055, 0.027, 0.026, 0.002, 0.060, 0.316, 0.009, 0.0  ],
        [0.016, 0.007, 0.021, 0.017, 0.057, 0.020, 0.040, 0.225],
        [0.036, 0.001, 0.003, 0.001, 0.028, 0.008, 0.066, 0.144],
    ], dtype=np.float64)
    # Symmetric effective coupling
    J_raw = np.sqrt(C * C.T)
    # Sign: inhibitory populations (indices 1,3,5,7) get negative coupling
    signs = np.ones(d)
    signs[[1, 3, 5, 7]] = -1.0
    sign_matrix = np.outer(signs, signs)
    J_eff = J_raw * sign_matrix
    np.fill_diagonal(J_eff, 0.0)
    # Scale to cm^-1 range (coupling_scale makes max |J| ~ 200 cm^-1)
    coupling_scale = 200.0 / np.max(np.abs(J_eff))
    H = np.diag(_PD_FREQS_HZ * _NEURAL_GAMMA_SCALE) + J_eff * coupling_scale
    return H

H_NEURAL_PD = _build_neural_pd_hamiltonian()

NEURAL_PD_JMAX_CM = np.max(np.abs(H_NEURAL_PD - np.diag(np.diag(H_NEURAL_PD))))

# ── Chromatin Hamiltonian: nucleosome-remodeler complex (d=8) ─────────────
# Illustrative molecular-oscillator case study of a nucleosome-remodeler complex.
# Sites: 0-1 = DNA entry/exit loci, 2-3 = histone tail contacts,
#        4-5 = remodeler ATPase contacts, 6-7 = linker DNA modes.
# Site energies in cm^-1 (THz protein mode range, 20-65 cm^-1).
# Coupling: ring topology (J_nn = 25 cm^-1) + cross-link 0↔7 (J_lr = 10 cm^-1).
# Biological dephasing: γ_bio = 200 cm^-1 (mid-range of 2D-IR measured
# protein aqueous dephasing at 310K; Fleming/Hamm/Tokmakoff groups).

_CHROMATIN_SITE_ENERGIES_CM = np.array([20, 35, 50, 65, 40, 55, 30, 45],
                                        dtype=np.float64)
_CHROMATIN_J_NN = 25.0   # cm^-1, nearest-neighbor ring coupling
_CHROMATIN_J_LR = 10.0   # cm^-1, cross-link 0↔7

def _build_chromatin_hamiltonian():
    """Build 8×8 chromatin nucleosome-remodeler Hamiltonian (cm^-1).

    Open chain (0-1-...-7) with nearest-neighbor coupling J_nn = 25 cm^-1,
    plus a cross-link between sites 0 and 7 at J_lr = 10 cm^-1.
    """
    d = 8
    H = np.diag(_CHROMATIN_SITE_ENERGIES_CM)
    for i in range(d - 1):
        H[i, i + 1] = _CHROMATIN_J_NN
        H[i + 1, i] = _CHROMATIN_J_NN
    H[0, 7] = _CHROMATIN_J_LR
    H[7, 0] = _CHROMATIN_J_LR
    return H

H_CHROMATIN = _build_chromatin_hamiltonian()
CHROMATIN_JMAX_CM = 25.0   # J_max = max nearest-neighbor coupling
CHROMATIN_BIO_GAMMA_CM = 200.0  # cm^-1, protein aqueous dephasing at 310K
CHROMATIN_BIO_GAMMA_OVER_JMAX = CHROMATIN_BIO_GAMMA_CM / CHROMATIN_JMAX_CM  # = 8.0


# ── Ion channel selectivity filter (d=4) ────────────────────────────────────
# Illustrative KcsA-like potassium channel selectivity filter.
# Four binding sites (S1-S4) in a linear chain.
# Site energies: outer sites (S1,S4) slightly lower due to proximity to
# bulk water; inner sites (S2,S3) deeper wells from carbonyl coordination.
# Couplings: nearest-neighbor ion-ion Coulomb + backbone-mediated, 15-40 cm^-1.
# Dephasing: protein-water, 100-300 cm^-1 (2D-IR literature).
_ION_CHANNEL_SITE_ENERGIES_CM = np.array([30.0, 45.0, 45.0, 30.0])  # cm^-1
_ION_CHANNEL_J_NN = 30.0  # cm^-1, nearest-neighbor coupling

def _build_ion_channel_hamiltonian():
    """Build 4×4 ion channel selectivity filter Hamiltonian (cm^-1).

    Linear chain S1-S2-S3-S4 with nearest-neighbor coupling J_nn = 30 cm^-1.
    """
    d = 4
    H = np.diag(_ION_CHANNEL_SITE_ENERGIES_CM)
    for i in range(d - 1):
        H[i, i + 1] = _ION_CHANNEL_J_NN
        H[i + 1, i] = _ION_CHANNEL_J_NN
    return H

H_ION_CHANNEL = _build_ion_channel_hamiltonian()
ION_CHANNEL_JMAX_CM = 30.0   # J_max = nearest-neighbor coupling
ION_CHANNEL_BIO_GAMMA_CM = 200.0  # cm^-1, protein aqueous dephasing at 310K
ION_CHANNEL_BIO_GAMMA_OVER_JMAX = ION_CHANNEL_BIO_GAMMA_CM / ION_CHANNEL_JMAX_CM  # ≈ 6.7


# ── Protein mid-fold microdomain (d=6) ───────────────────────────────────────
# Illustrative coarse-grained folding intermediate network.
# Sites represent metastable conformational microstates along a local folding
# trajectory (pre-nucleus -> partial pack -> near-native pocket).
# Energies/couplings are in molecular vibrational range (cm^-1), not full-protein
# thermodynamic free-energy units.
_PROTEIN_MIDFOLD_H_CM = np.array([
    [20.0, 25.0, 12.0,  0.0,  0.0,  0.0],
    [25.0, 35.0, 18.0, 22.0,  0.0,  0.0],
    [12.0, 18.0, 42.0, 20.0, 10.0,  0.0],
    [ 0.0, 22.0, 20.0, 50.0, 24.0, 12.0],
    [ 0.0,  0.0, 10.0, 24.0, 38.0, 26.0],
    [ 0.0,  0.0,  0.0, 12.0, 26.0, 28.0],
], dtype=np.float64)

def _build_protein_midfold_hamiltonian():
    """Build 6×6 illustrative protein mid-fold Hamiltonian (cm^-1)."""
    return _PROTEIN_MIDFOLD_H_CM.copy()

H_PROTEIN_MIDFOLD = _build_protein_midfold_hamiltonian()
PROTEIN_MIDFOLD_JMAX_CM = np.max(
    np.abs(H_PROTEIN_MIDFOLD - np.diag(np.diag(H_PROTEIN_MIDFOLD))))
PROTEIN_MIDFOLD_BIO_GAMMA_CM = 200.0
PROTEIN_MIDFOLD_BIO_GAMMA_OVER_JMAX = (
    PROTEIN_MIDFOLD_BIO_GAMMA_CM / PROTEIN_MIDFOLD_JMAX_CM)


# ── Complex registry ─────────────────────────────────────────────────────────

COMPLEXES = {
    'FMO': {
        'name': 'FMO',
        'labels': ['BChl1', 'BChl2', 'BChl3', 'BChl4', 'BChl5', 'BChl6', 'BChl7'],
        'H_cm': H_FMO,
        'initial_sites': [0, 5],
        'target_sites': [2],
        'kappa_ps': [1.0],              # ps^-1 per target site
        'gamma_scan_cm': (1, 500),      # cm^-1
        'reference': 'Adolphs & Renger 2006',
    },
    'PE545': {
        'name': 'PE545',
        'labels': ['PEB50/61C', 'DBV19A', 'DBV19B', 'PEB82C',
                    'PEB158C', 'PEB50/61D', 'PEB82D', 'PEB158D'],
        'H_cm': H_PE545,
        'initial_sites': [0, 5],
        'target_sites': [1, 2],
        'kappa_ps': [0.5, 0.5],         # ps^-1 each → 1.0 total
        'gamma_scan_cm': (1, 1500),     # cm^-1 (peak ~725 cm^-1)
        'reference': 'Novoderezhkin et al. 2010',
    },
    'Chromatin': {
        'name': 'Chromatin',
        'labels': ['DNA-entry', 'DNA-exit', 'HistTail1', 'HistTail2',
                    'ATPase1', 'ATPase2', 'LinkerDNA1', 'LinkerDNA2'],
        'H_cm': H_CHROMATIN,
        'initial_sites': [0],
        'target_sites': [3, 4],
        'kappa_ps': [0.5, 0.5],
        'gamma_scan_cm': (1, 1000),
        'reference': 'Illustrative (nucleosome-remodeler)',
    },
    'IonChannel': {
        'name': 'IonChannel',
        'labels': ['S1', 'S2', 'S3', 'S4'],
        'H_cm': H_ION_CHANNEL,
        'initial_sites': [0],
        'target_sites': [3],
        'kappa_ps': [0.5],
        'gamma_scan_cm': (1, 1000),
        'reference': 'Illustrative (KcsA-like selectivity filter)',
    },
    'ProteinMidFold': {
        'name': 'ProteinMidFold',
        'labels': ['U0', 'I1', 'I2', 'I3', 'I4', 'N*'],
        'H_cm': H_PROTEIN_MIDFOLD,
        'initial_sites': [0],
        'target_sites': [5],
        'kappa_ps': [0.5],
        'gamma_scan_cm': (1, 1000),
        'reference': 'Illustrative (protein mid-fold microdomain)',
    },
    'NeuralGamma': {
        'name': 'NeuralGamma',
        'labels': ['L4in', 'L4a', 'L4b', 'L2/3a', 'L2/3b',
                    'L2/3c', 'L5a', 'L5b', 'L6a', 'L6b'],
        'H_cm': H_NEURAL_GAMMA,
        'initial_sites': [0],
        'target_sites': [8, 9],
        'kappa_ps': [0.5, 0.5],         # ps^-1 each → 1.0 total
        'gamma_scan_cm': (1, 5000),     # cm^-1
        'reference': 'Synthetic (gamma-band microcircuit)',
    },
    'NeuralPD': {
        'name': 'NeuralPD',
        'labels': ['L2/3e', 'L2/3i', 'L4e', 'L4i',
                    'L5e', 'L5i', 'L6e', 'L6i'],
        'H_cm': H_NEURAL_PD,
        'initial_sites': [2],           # L4e (thalamic input)
        'target_sites': [4, 6],         # L5e, L6e (output layers)
        'kappa_ps': [0.5, 0.5],
        'gamma_scan_cm': (1, 5000),
        'reference': 'Potjans & Diesmann 2014',
    },
}


def kappa_ps_to_cm(kappa_ps):
    """Convert trapping rate from ps^-1 to cm^-1.

    1 ps^-1 ≈ 5.3 cm^-1 (via 1/ps = 1e12/s, divided by 2*pi*c).
    """
    return kappa_ps * 5.3


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


def lindblad_dephasing_sink_rhs(rho, H, gamma, kappa, sink_site=2,
                                target_sites=None, kappa_list=None):
    """Lindblad RHS with Haken-Strobl dephasing + trapping sink(s).

    Supports both single-site (backward-compatible) and multisite target manifolds.
    L_sink[ρ] = -½ Σ_j κ_j {|j⟩⟨j|, ρ} for each target site j.

    Parameters:
        rho: density matrix
        H: Hamiltonian in rad/fs
        gamma: dephasing rate in rad/fs
        kappa: trapping rate in rad/fs (used if target_sites is None)
        sink_site: 0-indexed site for single-site trap (backward compat)
        target_sites: list of 0-indexed target sites (overrides sink_site)
        kappa_list: list of per-site trapping rates in rad/fs (same length as target_sites)
    """
    drho = lindblad_dephasing_rhs(rho, H, gamma)
    d = H.shape[0]

    if target_sites is not None:
        for j, kj in zip(target_sites, kappa_list):
            proj = np.zeros((d, d), dtype=complex)
            proj[j, j] = 1.0
            drho -= 0.5 * kj * (proj @ rho + rho @ proj)
    else:
        proj_s = np.zeros((d, d), dtype=complex)
        proj_s[sink_site, sink_site] = 1.0
        drho -= 0.5 * kappa * (proj_s @ rho + rho @ proj_s)

    return drho


def compute_transport_efficiency(gamma_cm, kappa_trap_cm=5.3, t_final_fs=15000.0,
                                  dt_fs=1.0, initial_site=0, sink_site=2):
    """Compute FMO transport efficiency with Haken-Strobl dephasing + sink.

    Legacy single-site wrapper. For multisite targets, use compute_transport_generic.
    """
    H = H_FMO * CM_TO_RADFS
    gamma = gamma_cm * CM_TO_RADFS
    kappa = kappa_trap_cm * CM_TO_RADFS
    d = D_FMO

    rho = np.zeros((d, d), dtype=complex)
    rho[initial_site, initial_site] = 1.0

    n_steps = int(t_final_fs / dt_fs)
    trapped = 0.0

    for _ in range(n_steps):
        trapped += kappa * np.real(rho[sink_site, sink_site]) * dt_fs

        k1 = dt_fs * lindblad_dephasing_sink_rhs(rho, H, gamma, kappa, sink_site)
        k2 = dt_fs * lindblad_dephasing_sink_rhs(rho + 0.5 * k1, H, gamma, kappa, sink_site)
        k3 = dt_fs * lindblad_dephasing_sink_rhs(rho + 0.5 * k2, H, gamma, kappa, sink_site)
        k4 = dt_fs * lindblad_dephasing_sink_rhs(rho + k3, H, gamma, kappa, sink_site)
        rho = rho + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

        rho = 0.5 * (rho + rho.conj().T)

    return min(trapped, 1.0)


def compute_transport_generic(H_cm, gamma_cm, target_sites, kappa_ps_list,
                               initial_site, t_final_ps=15.0, dt_fs=1.0):
    """Compute transport efficiency for any complex with multisite target manifold.

    Parameters:
        H_cm: Hamiltonian in cm^-1
        gamma_cm: dephasing rate in cm^-1
        target_sites: list of 0-indexed target sites
        kappa_ps_list: list of per-site trapping rates in ps^-1
        initial_site: 0-indexed initial excitation site
        t_final_ps: evolution time in ps (default 15)
        dt_fs: time step in fs (default 1.0)

    Returns:
        efficiency: fraction of population trapped (0 to 1)
    """
    H = H_cm * CM_TO_RADFS
    gamma = gamma_cm * CM_TO_RADFS
    kappa_radfs = [k * kappa_ps_to_cm(1.0) * CM_TO_RADFS for k in kappa_ps_list]
    d = H_cm.shape[0]
    t_final_fs = t_final_ps * PS_TO_FS

    rho = np.zeros((d, d), dtype=complex)
    rho[initial_site, initial_site] = 1.0

    n_steps = int(t_final_fs / dt_fs)
    trapped = 0.0

    for _ in range(n_steps):
        # Accumulate: η = Σ_j κ_j ∫ p_j(t) dt
        for j, kj in zip(target_sites, kappa_radfs):
            trapped += kj * np.real(rho[j, j]) * dt_fs

        k1 = dt_fs * lindblad_dephasing_sink_rhs(rho, H, gamma, None, None,
                                                   target_sites, kappa_radfs)
        k2 = dt_fs * lindblad_dephasing_sink_rhs(rho + 0.5 * k1, H, gamma, None, None,
                                                   target_sites, kappa_radfs)
        k3 = dt_fs * lindblad_dephasing_sink_rhs(rho + 0.5 * k2, H, gamma, None, None,
                                                   target_sites, kappa_radfs)
        k4 = dt_fs * lindblad_dephasing_sink_rhs(rho + k3, H, gamma, None, None,
                                                   target_sites, kappa_radfs)
        rho = rho + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        rho = 0.5 * (rho + rho.conj().T)

    return min(trapped, 1.0)


def evolve_lindblad_sink_generic(H_cm, gamma_cm, target_sites, kappa_ps_list,
                                  rho0, t_final_ps=5.0, dt_fs=1.0):
    """Evolve density matrix with multisite target trapping.

    Returns the final (sub-normalized) density matrix.
    """
    H = H_cm * CM_TO_RADFS
    gamma = gamma_cm * CM_TO_RADFS
    kappa_radfs = [k * kappa_ps_to_cm(1.0) * CM_TO_RADFS for k in kappa_ps_list]
    t_final_fs = t_final_ps * PS_TO_FS

    rho = rho0.copy().astype(complex)
    n_steps = int(t_final_fs / dt_fs)

    for _ in range(n_steps):
        k1 = dt_fs * lindblad_dephasing_sink_rhs(rho, H, gamma, None, None,
                                                   target_sites, kappa_radfs)
        k2 = dt_fs * lindblad_dephasing_sink_rhs(rho + 0.5 * k1, H, gamma, None, None,
                                                   target_sites, kappa_radfs)
        k3 = dt_fs * lindblad_dephasing_sink_rhs(rho + 0.5 * k2, H, gamma, None, None,
                                                   target_sites, kappa_radfs)
        k4 = dt_fs * lindblad_dephasing_sink_rhs(rho + k3, H, gamma, None, None,
                                                   target_sites, kappa_radfs)
        rho = rho + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        rho = 0.5 * (rho + rho.conj().T)

    return rho


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


def l1_coherence(rho):
    """L1-norm of coherence: sum of |rho_kl| for k != l."""
    d = rho.shape[0]
    return np.real(sum(abs(rho[k, l]) for k in range(d) for l in range(d) if k != l))


def purity(rho):
    """Purity Tr[rho^2]."""
    return np.real(np.trace(rho @ rho))


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


if __name__ == "__main__":
    print("Photosynthetic Complex Analysis")
    print("=" * 60)

    # Run both FMO and PE545 through the generic pipeline
    for cname in ['FMO', 'PE545']:
        cx = COMPLEXES[cname]
        H_cm = cx['H_cm']
        d = H_cm.shape[0]
        targets = cx['target_sites']
        kappas = cx['kappa_ps']

        print(f"\n{'=' * 60}")
        print(f"{cname} ({cx['reference']})")
        print(f"  d={d}, targets={targets}, kappa_ps={kappas}")
        print("=" * 60)

        gamma_values = np.concatenate([
            np.array([1.0, 5.0, 10.0, 50.0, 100.0]),
            np.arange(200, 1600, 200),
        ])

        for site_idx in cx['initial_sites']:
            label = cx['labels'][site_idx]
            print(f"\n  Initial: {label} (site {site_idx})")
            etas = [compute_transport_generic(H_cm, g, targets, kappas, site_idx,
                                              t_final_ps=30.0, dt_fs=2.0)
                    for g in gamma_values]
            etas = np.array(etas)
            idx_max = np.argmax(etas)
            print(f"  Peak: eta={etas[idx_max]:.4f} at gamma={gamma_values[idx_max]:.0f} cm^-1")

            # Principal angles at a few gamma values
            rho0 = np.zeros((d, d), dtype=complex)
            rho0[site_idx, site_idx] = 1.0
            for gamma_cm in [1.0, 100.0, gamma_values[idx_max]]:
                rho = evolve_lindblad_sink_generic(H_cm, gamma_cm, targets, kappas,
                                                    rho0, t_final_ps=5.0, dt_fs=1.0)
                tr = np.real(np.trace(rho))
                rho_n = rho / tr if tr > 1e-12 else np.eye(d, dtype=complex) / d
                rho_reg = rho_n + 1e-10 * np.eye(d) / d
                rho_reg /= np.trace(rho_reg)
                angles = np.sort(np.degrees(principal_angles(rho_reg)))
                print(f"    gamma={gamma_cm:7.0f} cm^-1  min_angle={angles[0]:5.1f}  "
                      f"max_angle={angles[-1]:5.1f}  tr={tr:.4f}")

    print("\nDone.")
