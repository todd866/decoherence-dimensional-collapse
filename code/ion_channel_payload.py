#!/usr/bin/env python3
"""Minimal ion-channel payload analysis for the rebuilt manuscript.

This script computes one concrete molecular payload anchor:

- 4-site KcsA-like selectivity filter Hamiltonian
- Lindblad dephasing + target trapping evolution
- Bures principal-angle spectrum at the biological operating point
- Non-classicality index chi
- Transport optimum over a dimensionless gamma/J scan

Outputs:
- results/ion_channel_summary.json
- results/ion_channel_scan.csv
- results/ion_channel_sensitivity.csv
- results/ion_channel_topology_sensitivity.csv
- results/ion_payload_benchmarks.csv
- results/threshold_scaling_scan.csv
- results/threshold_scaling_sensitivity.csv
- results/biological_anchor_points.csv
- figures/biological_anchor_map.pdf
- figures/biological_anchor_map.png
- figures/carrier_payload_schematic.pdf
- figures/carrier_payload_schematic.png
- figures/ion_channel_anchor.pdf
- figures/ion_channel_anchor.png
- figures/threshold_entrainment_scaling.pdf
- figures/threshold_entrainment_scaling.png
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# 1 cm^-1 = 2*pi*c = 2*pi*2.998e10 s^-1 = 2*pi*2.998e-5 fs^-1
CM_TO_RADFS = 2.0 * np.pi * 2.998e-5
PS_TO_FS = 1000.0
KAPPA_PS_TO_CM = 5.3

ION_CHANNEL_SITE_ENERGIES_CM = np.array([30.0, 45.0, 45.0, 30.0], dtype=np.float64)
ION_CHANNEL_J_NN_CM = 30.0
ION_CHANNEL_BIO_GAMMA_CM = 200.0
ION_CHANNEL_TARGET_SITES = [3]
ION_CHANNEL_INITIAL_SITE = 0
ION_CHANNEL_KAPPA_PS = [0.5]
GEOMETRY_READOUT_PS = 5.0
TRANSPORT_WINDOW_PS = 8.0
SCAN_POINTS = 25
DT_FS = 5.0
PHASE_WINDOW_S = 2.0
PHASE_WINDOW_N0 = 1.0e4
PHASE_WINDOW_REDUNDANCIES = [1.0, 10.0, 100.0]
PHASE_WINDOW_RAMP_OVER_R0 = 4.0
PHASE_WINDOW_SOFT_ETA = 2.0
PHASE_WINDOW_SENSITIVITY_S_VALUES = [1.0, 2.0, 3.0]
PHASE_WINDOW_SENSITIVITY_RAMP_VALUES = [2.0, 4.0, 8.0]
PHASE_WINDOW_SENSITIVITY_ETA_VALUES = [1.0, 2.0, 4.0]
PHASE_WINDOW_HETERO_SEED = 7
PHASE_WINDOW_HETERO_MODULES = 12000
PHASE_WINDOW_HETERO_DISC_RADIUS_OVER_R0 = 6.0
PHASE_WINDOW_HETERO_LOGSIGMA = 0.35
BENCHMARK_N_VALUES = [1.0e3, 1.0e4, 1.0e5]
ION_CHANNEL_SENSITIVITY_J_VALUES = [15.0, 20.0, 30.0, 40.0]
ION_CHANNEL_SENSITIVITY_GAMMA_VALUES = [100.0, 150.0, 200.0, 250.0, 300.0]
ION_CHANNEL_TOPOLOGY_DELTA_VALUES = [10.0, 15.0, 20.0]
ION_CHANNEL_TOPOLOGY_JNNN_VALUES = [0.0, 5.0, 10.0]

# Fallback legacy anchors, used only if the rebuilt photosynthetic summary
# has not yet been generated.
LEGACY_ANCHORS = [
    {
        "system": "FMO",
        "role": "proof_of_principle",
        "status": "legacy-validated",
        "dimension": 7,
        "collapse": "48->6",
        "gamma_bio_over_jmax": 100.0 / 87.7,
        "theta_min_bio_deg": 87.7,
        "theta_reference": "biological point",
        "provenance": "big-overhaul-checkpoint",
    },
    {
        "system": "PE545",
        "role": "proof_of_principle",
        "status": "legacy-validated",
        "dimension": 8,
        "collapse": "63->7",
        "gamma_bio_over_jmax": 725.0 / 92.0,
        "theta_min_bio_deg": 89.0,
        "theta_reference": "biological point",
        "provenance": "big-overhaul-checkpoint",
    },
]


def build_ion_channel_hamiltonian_with_params(
    j_nn_cm: float,
    inner_delta_cm: float = 15.0,
    j_nnn_cm: float = 0.0,
) -> np.ndarray:
    """4x4 selectivity-filter Hamiltonian with simple topology controls."""
    site_energies_cm = np.array([30.0, 30.0 + inner_delta_cm, 30.0 + inner_delta_cm, 30.0], dtype=np.float64)
    h_cm = np.diag(site_energies_cm)
    for idx in range(3):
        h_cm[idx, idx + 1] = j_nn_cm
        h_cm[idx + 1, idx] = j_nn_cm
    if j_nnn_cm > 0.0:
        h_cm[0, 2] = j_nnn_cm
        h_cm[2, 0] = j_nnn_cm
        h_cm[1, 3] = j_nnn_cm
        h_cm[3, 1] = j_nnn_cm
    return h_cm


def build_ion_channel_hamiltonian() -> np.ndarray:
    """Baseline 4x4 KcsA-like selectivity-filter Hamiltonian in cm^-1."""
    return build_ion_channel_hamiltonian_with_params(
        ION_CHANNEL_J_NN_CM,
        inner_delta_cm=15.0,
        j_nnn_cm=0.0,
    )


H_ION_CHANNEL = build_ion_channel_hamiltonian()
ION_CHANNEL_JMAX_CM = ION_CHANNEL_J_NN_CM


def build_ion_channel_hamiltonian_with_coupling(j_nn_cm: float) -> np.ndarray:
    """4x4 KcsA-like selectivity-filter Hamiltonian with variable coupling."""
    return build_ion_channel_hamiltonian_with_params(j_nn_cm, inner_delta_cm=15.0, j_nnn_cm=0.0)


def kappa_ps_to_radfs(kappa_ps: float) -> float:
    """Convert trapping rate from ps^-1 to rad/fs."""
    return kappa_ps * KAPPA_PS_TO_CM * CM_TO_RADFS


def lindblad_dephasing_rhs(rho: np.ndarray, h_radfs: np.ndarray, gamma_radfs: float) -> np.ndarray:
    """Haken-Strobl dephasing in the site basis."""
    commutator = -1j * (h_radfs @ rho - rho @ h_radfs)
    d = rho.shape[0]
    dephasing = np.zeros_like(rho)
    for k in range(d):
        proj = np.zeros((d, d), dtype=complex)
        proj[k, k] = 1.0
        dephasing += proj @ rho @ proj - 0.5 * (proj @ rho + rho @ proj)
    return commutator + gamma_radfs * dephasing


def lindblad_dephasing_sink_rhs(
    rho: np.ndarray,
    h_radfs: np.ndarray,
    gamma_radfs: float,
    target_sites: list[int],
    kappa_radfs: list[float],
) -> np.ndarray:
    """Dephasing plus irreversible trapping at target sites."""
    drho = lindblad_dephasing_rhs(rho, h_radfs, gamma_radfs)
    d = rho.shape[0]
    for site, rate in zip(target_sites, kappa_radfs):
        proj = np.zeros((d, d), dtype=complex)
        proj[site, site] = 1.0
        drho -= 0.5 * rate * (proj @ rho + rho @ proj)
    return drho


def evolve_lindblad_sink(
    h_cm: np.ndarray,
    gamma_cm: float,
    target_sites: list[int],
    kappa_ps: list[float],
    rho0: np.ndarray,
    t_final_ps: float,
    dt_fs: float = 1.0,
) -> np.ndarray:
    """Evolve the sub-normalized state under dephasing plus sink."""
    h_radfs = h_cm * CM_TO_RADFS
    gamma_radfs = gamma_cm * CM_TO_RADFS
    kappa_radfs = [kappa_ps_to_radfs(k) for k in kappa_ps]

    rho = rho0.astype(complex).copy()
    n_steps = int((t_final_ps * PS_TO_FS) / dt_fs)

    for _ in range(n_steps):
        k1 = dt_fs * lindblad_dephasing_sink_rhs(rho, h_radfs, gamma_radfs, target_sites, kappa_radfs)
        k2 = dt_fs * lindblad_dephasing_sink_rhs(rho + 0.5 * k1, h_radfs, gamma_radfs, target_sites, kappa_radfs)
        k3 = dt_fs * lindblad_dephasing_sink_rhs(rho + 0.5 * k2, h_radfs, gamma_radfs, target_sites, kappa_radfs)
        k4 = dt_fs * lindblad_dephasing_sink_rhs(rho + k3, h_radfs, gamma_radfs, target_sites, kappa_radfs)
        rho = rho + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        rho = 0.5 * (rho + rho.conj().T)

    return rho


def compute_transport_efficiency(
    h_cm: np.ndarray,
    gamma_cm: float,
    target_sites: list[int],
    kappa_ps: list[float],
    initial_site: int,
    t_final_ps: float,
    dt_fs: float = 1.0,
) -> float:
    """Integrated trapping yield."""
    h_radfs = h_cm * CM_TO_RADFS
    gamma_radfs = gamma_cm * CM_TO_RADFS
    kappa_radfs = [kappa_ps_to_radfs(k) for k in kappa_ps]

    d = h_cm.shape[0]
    rho = np.zeros((d, d), dtype=complex)
    rho[initial_site, initial_site] = 1.0
    trapped = 0.0
    n_steps = int((t_final_ps * PS_TO_FS) / dt_fs)

    for _ in range(n_steps):
        for site, rate in zip(target_sites, kappa_radfs):
            trapped += rate * float(np.real(rho[site, site])) * dt_fs

        k1 = dt_fs * lindblad_dephasing_sink_rhs(rho, h_radfs, gamma_radfs, target_sites, kappa_radfs)
        k2 = dt_fs * lindblad_dephasing_sink_rhs(rho + 0.5 * k1, h_radfs, gamma_radfs, target_sites, kappa_radfs)
        k3 = dt_fs * lindblad_dephasing_sink_rhs(rho + 0.5 * k2, h_radfs, gamma_radfs, target_sites, kappa_radfs)
        k4 = dt_fs * lindblad_dephasing_sink_rhs(rho + k3, h_radfs, gamma_radfs, target_sites, kappa_radfs)
        rho = rho + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        rho = 0.5 * (rho + rho.conj().T)

    return min(trapped, 1.0)


def regularize_density_matrix(rho: np.ndarray, eps: float = 1e-10) -> np.ndarray:
    """Normalize surviving state and keep it full rank for geometry."""
    tr = float(np.real(np.trace(rho)))
    d = rho.shape[0]
    if tr <= 1e-12 or not np.all(np.isfinite(rho)):
        return np.eye(d, dtype=complex) / d

    rho_n = rho / tr
    rho_reg = rho_n + eps * np.eye(d) / d
    rho_reg = 0.5 * (rho_reg + rho_reg.conj().T)
    rho_reg /= np.trace(rho_reg)
    return rho_reg


def principal_angles(rho: np.ndarray) -> np.ndarray:
    """Principal angles between diagonal and off-diagonal tangent spaces."""
    d = rho.shape[0]
    eigvals, eigvecs = np.linalg.eigh(rho)
    eigvals = np.maximum(np.real(eigvals), 1e-14)

    diag_basis = []
    for a in range(d - 1):
        basis_vec = np.zeros((d, d), dtype=complex)
        basis_vec[a, a] = 1.0
        basis_vec[a + 1, a + 1] = -1.0
        diag_basis.append(basis_vec)

    off_basis = []
    for k in range(d):
        for l in range(k + 1, d):
            real_part = np.zeros((d, d), dtype=complex)
            real_part[k, l] = 1.0
            real_part[l, k] = 1.0
            off_basis.append(real_part)

            imag_part = np.zeros((d, d), dtype=complex)
            imag_part[k, l] = 1j
            imag_part[l, k] = -1j
            off_basis.append(imag_part)

    def sld(direction: np.ndarray) -> np.ndarray:
        direction_eig = eigvecs.conj().T @ direction @ eigvecs
        sld_eig = np.zeros((d, d), dtype=complex)
        for m in range(d):
            for n in range(d):
                denom = eigvals[m] + eigvals[n]
                if denom > 1e-14:
                    sld_eig[m, n] = 2.0 * direction_eig[m, n] / denom
        return eigvecs @ sld_eig @ eigvecs.conj().T

    def bures_inner(x: np.ndarray, y: np.ndarray) -> float:
        lx = sld(x)
        ly = sld(y)
        return 0.25 * np.real(np.trace(rho @ (lx @ ly + ly @ lx))) / 2.0

    n_diag = len(diag_basis)
    n_off = len(off_basis)

    g_dd = np.zeros((n_diag, n_diag))
    g_oo = np.zeros((n_off, n_off))
    g_do = np.zeros((n_diag, n_off))

    for i in range(n_diag):
        for j in range(i, n_diag):
            val = bures_inner(diag_basis[i], diag_basis[j])
            g_dd[i, j] = val
            g_dd[j, i] = val

    for i in range(n_off):
        for j in range(i, n_off):
            val = bures_inner(off_basis[i], off_basis[j])
            g_oo[i, j] = val
            g_oo[j, i] = val

    for i in range(n_diag):
        for j in range(n_off):
            g_do[i, j] = bures_inner(diag_basis[i], off_basis[j])

    def matrix_sqrt_inv(mat: np.ndarray) -> np.ndarray:
        evals, evecs = np.linalg.eigh(mat)
        evals = np.maximum(evals, 1e-12)
        return evecs @ np.diag(1.0 / np.sqrt(evals)) @ evecs.T

    m = matrix_sqrt_inv(g_dd) @ g_do @ matrix_sqrt_inv(g_oo)
    singular_values = np.linalg.svd(m, compute_uv=False)
    cosines = np.clip(singular_values[:n_diag], 0.0, 1.0)
    return np.arccos(cosines)


def non_classicality_index(angles_rad: np.ndarray) -> float:
    """chi = mean(cos^2(theta_a))."""
    return float(np.mean(np.cos(angles_rad) ** 2))


def local_nonclassical_load(chi_value: float, dimension: int) -> float:
    """L = chi * (d^2 - d)."""
    return float(chi_value * (dimension * dimension - dimension))


def geometry_snapshot(gamma_cm: float) -> dict[str, float | list[float]]:
    """Compute the geometry of the surviving state at one gamma value."""
    rho0 = np.zeros((4, 4), dtype=complex)
    rho0[ION_CHANNEL_INITIAL_SITE, ION_CHANNEL_INITIAL_SITE] = 1.0
    rho = evolve_lindblad_sink(
        H_ION_CHANNEL,
        gamma_cm,
        ION_CHANNEL_TARGET_SITES,
        ION_CHANNEL_KAPPA_PS,
        rho0,
        t_final_ps=GEOMETRY_READOUT_PS,
        dt_fs=DT_FS,
    )
    rho_reg = regularize_density_matrix(rho)
    angles_rad = np.sort(principal_angles(rho_reg))
    angles_deg = [float(np.degrees(val)) for val in angles_rad]
    return {
        "theta_min_deg": angles_deg[0],
        "theta_spectrum_deg": angles_deg,
        "chi": non_classicality_index(angles_rad),
        "trace_after_readout": float(np.real(np.trace(rho))),
    }


def compute_ion_channel_sensitivity_rows() -> list[dict[str, float]]:
    """Sweep a small, physically reasonable ion-channel parameter band."""
    rows: list[dict[str, float]] = []
    for j_nn_cm in ION_CHANNEL_SENSITIVITY_J_VALUES:
        h_cm = build_ion_channel_hamiltonian_with_coupling(j_nn_cm)
        for gamma_cm in ION_CHANNEL_SENSITIVITY_GAMMA_VALUES:
            rho0 = np.zeros((4, 4), dtype=complex)
            rho0[ION_CHANNEL_INITIAL_SITE, ION_CHANNEL_INITIAL_SITE] = 1.0
            rho = evolve_lindblad_sink(
                h_cm,
                gamma_cm,
                ION_CHANNEL_TARGET_SITES,
                ION_CHANNEL_KAPPA_PS,
                rho0,
                t_final_ps=GEOMETRY_READOUT_PS,
                dt_fs=DT_FS,
            )
            rho_reg = regularize_density_matrix(rho)
            angles_rad = np.sort(principal_angles(rho_reg))
            theta_min_deg = float(np.degrees(angles_rad[0]))
            chi = float(non_classicality_index(angles_rad))
            rows.append(
                {
                    "j_nn_cm": float(j_nn_cm),
                    "gamma_cm": float(gamma_cm),
                    "gamma_over_j": float(gamma_cm / j_nn_cm),
                    "theta_min_deg": theta_min_deg,
                    "chi": chi,
                    "local_load": float(local_nonclassical_load(chi, 4)),
                }
            )
    return rows


def compute_ion_channel_topology_rows() -> list[dict[str, float]]:
    """Sweep a small topology band around the baseline ion-channel anchor."""
    rows: list[dict[str, float]] = []
    for inner_delta_cm in ION_CHANNEL_TOPOLOGY_DELTA_VALUES:
        for j_nnn_cm in ION_CHANNEL_TOPOLOGY_JNNN_VALUES:
            h_cm = build_ion_channel_hamiltonian_with_params(
                ION_CHANNEL_J_NN_CM,
                inner_delta_cm=inner_delta_cm,
                j_nnn_cm=j_nnn_cm,
            )
            rho0 = np.zeros((4, 4), dtype=complex)
            rho0[ION_CHANNEL_INITIAL_SITE, ION_CHANNEL_INITIAL_SITE] = 1.0
            rho = evolve_lindblad_sink(
                h_cm,
                ION_CHANNEL_BIO_GAMMA_CM,
                ION_CHANNEL_TARGET_SITES,
                ION_CHANNEL_KAPPA_PS,
                rho0,
                t_final_ps=GEOMETRY_READOUT_PS,
                dt_fs=DT_FS,
            )
            rho_reg = regularize_density_matrix(rho)
            angles_rad = np.sort(principal_angles(rho_reg))
            theta_min_deg = float(np.degrees(angles_rad[0]))
            chi = float(non_classicality_index(angles_rad))
            rows.append(
                {
                    "inner_delta_cm": float(inner_delta_cm),
                    "j_nnn_cm": float(j_nnn_cm),
                    "gamma_over_j": float(ION_CHANNEL_BIO_GAMMA_CM / ION_CHANNEL_JMAX_CM),
                    "theta_min_deg": theta_min_deg,
                    "chi": chi,
                    "local_load": float(local_nonclassical_load(chi, 4)),
                }
            )
    return rows


def summarize_ion_channel() -> tuple[dict, list[dict[str, float]]]:
    """Compute a reproducible ion-channel anchor summary."""
    gamma_over_j_values = np.logspace(-2, np.log10(30.0), SCAN_POINTS)
    scan_rows: list[dict[str, float]] = []
    efficiencies = []

    for gamma_over_j in gamma_over_j_values:
        gamma_cm = gamma_over_j * ION_CHANNEL_JMAX_CM
        eta = compute_transport_efficiency(
            H_ION_CHANNEL,
            gamma_cm,
            ION_CHANNEL_TARGET_SITES,
            ION_CHANNEL_KAPPA_PS,
            ION_CHANNEL_INITIAL_SITE,
            t_final_ps=TRANSPORT_WINDOW_PS,
            dt_fs=DT_FS,
        )
        eta_val = float(eta) if np.isfinite(eta) else float("nan")
        scan_rows.append(
            {
                "gamma_over_j": float(gamma_over_j),
                "gamma_cm": float(gamma_cm),
                "eta": eta_val,
            }
        )
        efficiencies.append(eta_val)

    idx_peak = int(np.nanargmax(efficiencies))
    peak_row = scan_rows[idx_peak]
    peak_is_interior = 0 < idx_peak < (len(scan_rows) - 1)
    bio_geometry = geometry_snapshot(ION_CHANNEL_BIO_GAMMA_CM)
    eta_bio = compute_transport_efficiency(
        H_ION_CHANNEL,
        ION_CHANNEL_BIO_GAMMA_CM,
        ION_CHANNEL_TARGET_SITES,
        ION_CHANNEL_KAPPA_PS,
        ION_CHANNEL_INITIAL_SITE,
        t_final_ps=TRANSPORT_WINDOW_PS,
        dt_fs=DT_FS,
    )
    sensitivity_rows = compute_ion_channel_sensitivity_rows()
    sensitivity_theta = [row["theta_min_deg"] for row in sensitivity_rows]
    sensitivity_chi = [row["chi"] for row in sensitivity_rows]
    topology_rows = compute_ion_channel_topology_rows()
    topology_theta = [row["theta_min_deg"] for row in topology_rows]
    topology_chi = [row["chi"] for row in topology_rows]

    summary = {
        "model": "KcsA-like selectivity filter",
        "dimension": 4,
        "initial_site": ION_CHANNEL_INITIAL_SITE,
        "target_sites": ION_CHANNEL_TARGET_SITES,
        "kappa_ps": ION_CHANNEL_KAPPA_PS,
        "jmax_cm": ION_CHANNEL_JMAX_CM,
        "gamma_bio_cm": ION_CHANNEL_BIO_GAMMA_CM,
        "gamma_bio_over_jmax": ION_CHANNEL_BIO_GAMMA_CM / ION_CHANNEL_JMAX_CM,
        "eta_bio": float(eta_bio),
        "theta_min_bio_deg": bio_geometry["theta_min_deg"],
        "theta_spectrum_bio_deg": bio_geometry["theta_spectrum_deg"],
        "chi_bio": bio_geometry["chi"],
        "trace_after_bio_readout": bio_geometry["trace_after_readout"],
        "transport_scan_peak_gamma_cm": peak_row["gamma_cm"],
        "transport_scan_peak_gamma_over_jmax": peak_row["gamma_over_j"],
        "transport_scan_peak_eta": peak_row["eta"],
        "transport_scan_peak_interior": peak_is_interior,
        "geometry_readout_ps": GEOMETRY_READOUT_PS,
        "transport_window_ps": TRANSPORT_WINDOW_PS,
        "scan_points": len(scan_rows),
        "sensitivity_j_values_cm": ION_CHANNEL_SENSITIVITY_J_VALUES,
        "sensitivity_gamma_values_cm": ION_CHANNEL_SENSITIVITY_GAMMA_VALUES,
        "sensitivity_theta_min_deg_min": float(min(sensitivity_theta)),
        "sensitivity_theta_min_deg_max": float(max(sensitivity_theta)),
        "sensitivity_chi_min": float(min(sensitivity_chi)),
        "sensitivity_chi_max": float(max(sensitivity_chi)),
        "topology_delta_values_cm": ION_CHANNEL_TOPOLOGY_DELTA_VALUES,
        "topology_j_nnn_values_cm": ION_CHANNEL_TOPOLOGY_JNNN_VALUES,
        "topology_theta_min_deg_min": float(min(topology_theta)),
        "topology_theta_min_deg_max": float(max(topology_theta)),
        "topology_chi_min": float(min(topology_chi)),
        "topology_chi_max": float(max(topology_chi)),
    }
    return summary, scan_rows


def write_outputs(summary: dict, scan_rows: list[dict[str, float]]) -> None:
    """Write JSON summary and CSV scan for the manuscript."""
    root = Path(__file__).resolve().parents[1]
    results_dir = root / "results"
    results_dir.mkdir(exist_ok=True)

    summary_path = results_dir / "ion_channel_summary.json"
    scan_path = results_dir / "ion_channel_scan.csv"

    with summary_path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
        fh.write("\n")

    with scan_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["gamma_over_j", "gamma_cm", "eta"],
        )
        writer.writeheader()
        writer.writerows(scan_rows)


def write_figure(summary: dict, scan_rows: list[dict[str, float]]) -> None:
    """Write a compact two-panel figure for the manuscript."""
    root = Path(__file__).resolve().parents[1]
    figures_dir = root / "figures"
    figures_dir.mkdir(exist_ok=True)

    gamma_over_j = np.array([row["gamma_over_j"] for row in scan_rows], dtype=float)
    eta = np.array([row["eta"] for row in scan_rows], dtype=float)
    theta_spectrum = np.array(summary["theta_spectrum_bio_deg"], dtype=float)

    fig, (ax_eta, ax_theta) = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)

    ax_eta.semilogx(gamma_over_j, eta, color="#1f4e79", lw=2.0)
    ax_eta.axvline(summary["gamma_bio_over_jmax"], color="#b03a2e", ls="--", lw=1.1)
    ax_eta.scatter(
        [summary["gamma_bio_over_jmax"]],
        [summary["eta_bio"]],
        color="#b03a2e",
        s=24,
        zorder=3,
    )
    ax_eta.set_xlabel(r"$\gamma/J_{\max}$")
    ax_eta.set_ylabel(r"Bounded yield $\eta$")
    ax_eta.set_title("Trap scan")
    ax_eta.text(
        0.05,
        0.10,
        "bio point",
        color="#b03a2e",
        transform=ax_eta.transAxes,
        fontsize=8,
    )

    x = np.arange(1, len(theta_spectrum) + 1)
    ax_theta.bar(x, theta_spectrum, color="#4c956c", width=0.65)
    ax_theta.axhline(90.0, color="black", ls=":", lw=1.0)
    ax_theta.set_xticks(x)
    ax_theta.set_xlabel(r"Angle index $a$")
    ax_theta.set_ylabel(r"$\theta_a$ (degrees)")
    ax_theta.set_ylim(86.5, 90.2)
    ax_theta.set_title("Biological-point geometry")
    ax_theta.text(
        0.05,
        0.10,
        rf"$\chi \approx {summary['chi_bio']:.2e}$",
        transform=ax_theta.transAxes,
        fontsize=8,
    )

    fig.savefig(figures_dir / "ion_channel_anchor.pdf", bbox_inches="tight")
    fig.savefig(figures_dir / "ion_channel_anchor.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def compute_soft_bridge_curves(
    *,
    s_dim: float,
    ramp_over_r0: float,
    eta_soft: float,
) -> dict[str, np.ndarray | float]:
    """Compute analytic and heterogeneous soft bridge curves for one parameter set."""
    freq_ratio = np.logspace(np.log10(0.25), np.log10(4.0), 300)

    r_phase_over_r0 = freq_ratio ** (-1.0)
    r_amp_over_r0 = np.full_like(freq_ratio, ramp_over_r0)
    inv_reff_eta = r_amp_over_r0 ** (-eta_soft) + r_phase_over_r0 ** (-eta_soft)
    r_eff_over_r0 = inv_reff_eta ** (-1.0 / eta_soft)
    r_eff_ref = (ramp_over_r0 ** (-eta_soft) + 1.0) ** (-1.0 / eta_soft)
    n_payload = PHASE_WINDOW_N0 * (r_eff_over_r0 / r_eff_ref) ** s_dim

    rng = np.random.default_rng(PHASE_WINDOW_HETERO_SEED)
    radii = PHASE_WINDOW_HETERO_DISC_RADIUS_OVER_R0 * np.sqrt(rng.random(PHASE_WINDOW_HETERO_MODULES))
    ramp_i = np.exp(
        rng.normal(
            loc=np.log(ramp_over_r0),
            scale=PHASE_WINDOW_HETERO_LOGSIGMA,
            size=PHASE_WINDOW_HETERO_MODULES,
        )
    )
    rphi0_i = np.exp(
        rng.normal(
            loc=0.0,
            scale=PHASE_WINDOW_HETERO_LOGSIGMA,
            size=PHASE_WINDOW_HETERO_MODULES,
        )
    )
    weights_f0 = np.exp(-((radii / ramp_i) ** eta_soft) - ((radii / rphi0_i) ** eta_soft))
    hetero_scale = PHASE_WINDOW_N0 / max(float(np.mean(weights_f0)), 1e-12)
    n_payload_hetero = np.zeros_like(freq_ratio)
    for idx, f_ratio in enumerate(freq_ratio):
        weights = np.exp(-((radii / ramp_i) ** eta_soft) - ((f_ratio * radii / rphi0_i) ** eta_soft))
        n_payload_hetero[idx] = hetero_scale * float(np.mean(weights))

    return {
        "freq_ratio": freq_ratio,
        "r_phase_over_r0": r_phase_over_r0,
        "r_eff_over_r0": r_eff_over_r0,
        "n_payload_analytic": n_payload,
        "n_payload_heterogeneous": n_payload_hetero,
        "analytic_gain_low_over_high": float(n_payload[0] / n_payload[-1]),
        "hetero_gain_low_over_high": float(n_payload_hetero[0] / n_payload_hetero[-1]),
        "monotone_analytic": float(np.all(np.diff(n_payload) <= 1e-9)),
        "monotone_heterogeneous": float(np.all(np.diff(n_payload_hetero) <= 1e-9)),
    }


def write_phase_window_sensitivity_outputs(summary: dict) -> None:
    """Write a compact parameter sweep for the soft carrier bridge."""
    root = Path(__file__).resolve().parents[1]
    results_dir = root / "results"
    results_dir.mkdir(exist_ok=True)

    local_load = local_nonclassical_load(summary["chi_bio"], summary["dimension"])
    csv_path = results_dir / "threshold_scaling_sensitivity.csv"

    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "s_dim",
                "eta_soft",
                "ramp_over_r0",
                "analytic_gain_low_over_high",
                "hetero_gain_low_over_high",
                "monotone_analytic",
                "monotone_heterogeneous",
                "leff_low_varrho_1",
                "leff_high_varrho_1",
            ]
        )
        for s_dim in PHASE_WINDOW_SENSITIVITY_S_VALUES:
            for eta_soft in PHASE_WINDOW_SENSITIVITY_ETA_VALUES:
                for ramp_over_r0 in PHASE_WINDOW_SENSITIVITY_RAMP_VALUES:
                    curves = compute_soft_bridge_curves(
                        s_dim=s_dim,
                        ramp_over_r0=ramp_over_r0,
                        eta_soft=eta_soft,
                    )
                    if abs(s_dim - PHASE_WINDOW_S) < 1e-12:
                        n_payload_hetero = curves["n_payload_heterogeneous"]
                        hetero_gain = curves["hetero_gain_low_over_high"]
                        monotone_heterogeneous = int(bool(curves["monotone_heterogeneous"]))
                        leff_low_varrho_1 = float(n_payload_hetero[0] * local_load / PHASE_WINDOW_REDUNDANCIES[0])
                        leff_high_varrho_1 = float(n_payload_hetero[-1] * local_load / PHASE_WINDOW_REDUNDANCIES[0])
                    else:
                        hetero_gain = ""
                        monotone_heterogeneous = ""
                        leff_low_varrho_1 = ""
                        leff_high_varrho_1 = ""
                    writer.writerow(
                        [
                            s_dim,
                            eta_soft,
                            ramp_over_r0,
                            curves["analytic_gain_low_over_high"],
                            hetero_gain,
                            int(bool(curves["monotone_analytic"])),
                            monotone_heterogeneous,
                            leff_low_varrho_1,
                            leff_high_varrho_1,
                        ]
                    )


def write_phase_window_scaling_outputs(summary: dict) -> None:
    """Write a soft phase-coupling bridge figure and CSV using the ion-channel anchor."""
    root = Path(__file__).resolve().parents[1]
    results_dir = root / "results"
    figures_dir = root / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)

    local_load = local_nonclassical_load(summary["chi_bio"], summary["dimension"])
    eta = PHASE_WINDOW_SOFT_ETA
    curves = compute_soft_bridge_curves(
        s_dim=PHASE_WINDOW_S,
        ramp_over_r0=PHASE_WINDOW_RAMP_OVER_R0,
        eta_soft=eta,
    )
    freq_ratio = curves["freq_ratio"]
    r_phase_over_r0 = curves["r_phase_over_r0"]
    r_amp_over_r0 = np.full_like(freq_ratio, PHASE_WINDOW_RAMP_OVER_R0)
    r_eff_over_r0 = curves["r_eff_over_r0"]
    n_payload = curves["n_payload_analytic"]
    n_payload_hetero = curves["n_payload_heterogeneous"]
    crossover_ratio = 1.0 / PHASE_WINDOW_RAMP_OVER_R0

    csv_path = results_dir / "threshold_scaling_scan.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "f_over_f0",
                "r_phase_over_r0",
                "r_amp_over_r0",
                "r_eff_over_r0",
                "n_payload_soft_analytic",
                "n_payload_soft_heterogeneous",
                "classical_baseline",
                "leff_varrho_1",
                "leff_varrho_10",
                "leff_varrho_100",
                "leff_hetero_varrho_1",
                "leff_hetero_varrho_10",
                "leff_hetero_varrho_100",
                "uplift_frac_varrho_1",
                "uplift_frac_varrho_10",
                "uplift_frac_varrho_100",
            ]
        )
        for idx, f_ratio in enumerate(freq_ratio):
            classical_baseline = float(n_payload[idx] * (summary["dimension"] - 1))
            writer.writerow(
                [
                    float(f_ratio),
                    float(r_phase_over_r0[idx]),
                    float(r_amp_over_r0[idx]),
                    float(r_eff_over_r0[idx]),
                    float(n_payload[idx]),
                    float(n_payload_hetero[idx]),
                    classical_baseline,
                    float(n_payload[idx] * local_load / PHASE_WINDOW_REDUNDANCIES[0]),
                    float(n_payload[idx] * local_load / PHASE_WINDOW_REDUNDANCIES[1]),
                    float(n_payload[idx] * local_load / PHASE_WINDOW_REDUNDANCIES[2]),
                    float(n_payload_hetero[idx] * local_load / PHASE_WINDOW_REDUNDANCIES[0]),
                    float(n_payload_hetero[idx] * local_load / PHASE_WINDOW_REDUNDANCIES[1]),
                    float(n_payload_hetero[idx] * local_load / PHASE_WINDOW_REDUNDANCIES[2]),
                    float(local_load / ((summary["dimension"] - 1) * PHASE_WINDOW_REDUNDANCIES[0])),
                    float(local_load / ((summary["dimension"] - 1) * PHASE_WINDOW_REDUNDANCIES[1])),
                    float(local_load / ((summary["dimension"] - 1) * PHASE_WINDOW_REDUNDANCIES[2])),
                ]
            )

    fig, (ax_field, ax_scale) = plt.subplots(1, 2, figsize=(7.2, 3.1), constrained_layout=True)

    ax_field.loglog(freq_ratio, r_phase_over_r0, color="#1f4e79", lw=2.0,
                    label=r"$R_{\phi}(f)/R_0 \propto f^{-1}$")
    ax_field.loglog(freq_ratio, r_amp_over_r0, color="#b03a2e", lw=2.0, ls="--",
                    label=r"$R_{\rm amp}/R_0$")
    ax_field.loglog(freq_ratio, r_eff_over_r0, color="#4c956c", lw=2.2,
                    label=rf"$R_{{\rm eff}}^{{(\eta)}}(f)/R_0,\ \eta={eta:.0f}$")
    ax_field.axvline(crossover_ratio, color="black", ls=":", lw=1.0)
    ax_field.set_xlabel(r"carrier frequency $f/f_0$")
    ax_field.set_ylabel(r"normalized coordination radius")
    ax_field.set_title("Soft phase-coupling radius")
    ax_field.invert_xaxis()
    ax_field.legend(frameon=False, fontsize=8, loc="lower left")
    ax_field.text(
        0.05,
        0.08,
        r"$R_{\rm eff}^{(\eta)}=(R_{\rm amp}^{-\eta}+R_{\phi}^{-\eta})^{-1/\eta}$" "\n"
        r"$R_{\phi}=v\Delta\phi_{\max}/(2\pi f)$",
        transform=ax_field.transAxes,
        fontsize=8,
    )

    colors = ["#4c956c", "#6c5ce7", "#b08968"]
    for color, redundancy in zip(colors, PHASE_WINDOW_REDUNDANCIES):
        leff = n_payload * local_load / redundancy
        leff_hetero = n_payload_hetero * local_load / redundancy
        ax_scale.loglog(
            freq_ratio,
            leff_hetero,
            color=color,
            lw=2.0,
            label=rf"$\varrho={int(redundancy)}$",
        )
        ax_scale.loglog(
            freq_ratio,
            leff,
            color=color,
            lw=1.4,
            ls="--",
            alpha=0.8,
        )
    ax_scale.set_xlabel(r"carrier frequency $f/f_0$")
    ax_scale.set_ylabel(r"effective payload load $\mathcal{L}_{\rm eff}$")
    ax_scale.set_title(r"Soft payload scaling: heterogeneous vs. analytic")
    ax_scale.invert_xaxis()
    ax_scale.text(
        0.05,
        0.08,
        rf"$L_{{\rm ion}}\approx {local_load:.2e}$" "\n"
        rf"$N_0={int(PHASE_WINDOW_N0):d},\ s={PHASE_WINDOW_S:.0f},\ R_{{\rm amp}}/R_0={PHASE_WINDOW_RAMP_OVER_R0:.0f},\ \eta={eta:.0f}$" "\n"
        r"solid: hetero soft field, dashed: soft analytic envelope",
        transform=ax_scale.transAxes,
        fontsize=8,
    )
    ax_scale.legend(frameon=False, fontsize=8, loc="upper right")

    fig.savefig(figures_dir / "threshold_entrainment_scaling.pdf", bbox_inches="tight")
    fig.savefig(figures_dir / "threshold_entrainment_scaling.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_ion_channel_sensitivity() -> None:
    """Write a compact robustness sweep for the ion-channel anchor."""
    root = Path(__file__).resolve().parents[1]
    results_dir = root / "results"
    results_dir.mkdir(exist_ok=True)

    rows = compute_ion_channel_sensitivity_rows()
    csv_path = results_dir / "ion_channel_sensitivity.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "j_nn_cm",
                "gamma_cm",
                "gamma_over_j",
                "theta_min_deg",
                "chi",
                "local_load",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def write_ion_channel_topology_sensitivity() -> None:
    """Write a small topology sweep around the baseline ion-channel anchor."""
    root = Path(__file__).resolve().parents[1]
    results_dir = root / "results"
    results_dir.mkdir(exist_ok=True)

    rows = compute_ion_channel_topology_rows()
    csv_path = results_dir / "ion_channel_topology_sensitivity.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "inner_delta_cm",
                "j_nnn_cm",
                "gamma_over_j",
                "theta_min_deg",
                "chi",
                "local_load",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def write_payload_benchmarks(summary: dict) -> None:
    """Write a compact benchmark table for concrete entrained site counts."""
    root = Path(__file__).resolve().parents[1]
    results_dir = root / "results"
    results_dir.mkdir(exist_ok=True)

    local_load = local_nonclassical_load(summary["chi_bio"], summary["dimension"])
    d_minus_1 = summary["dimension"] - 1
    csv_path = results_dir / "ion_payload_benchmarks.csv"

    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "n_payload",
                "varrho",
                "classical_baseline",
                "leff",
                "dcoord_total",
                "uplift_frac",
            ]
        )
        for n_payload in BENCHMARK_N_VALUES:
            classical_baseline = n_payload * d_minus_1
            for redundancy in PHASE_WINDOW_REDUNDANCIES:
                leff = n_payload * local_load / redundancy
                writer.writerow(
                    [
                        int(n_payload),
                        int(redundancy),
                        float(classical_baseline),
                        float(leff),
                        float(classical_baseline + leff),
                        float(local_load / (d_minus_1 * redundancy)),
                    ]
                    )


def load_photosynthetic_anchors(root: Path) -> list[dict[str, str | float | int]]:
    """Load rebuilt photosynthetic anchor points if available."""
    csv_path = root / "results" / "photosynthetic_anchor_points.csv"
    if not csv_path.exists():
        return [dict(row) for row in LEGACY_ANCHORS]

    rows: list[dict[str, str | float | int]] = []
    with csv_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(
                {
                    "system": row["system"],
                    "role": "proof_of_principle",
                    "status": row["source"],
                    "dimension": int(row["dimension"]),
                    "collapse": row["collapse"],
                    "gamma_bio_over_jmax": float(row["gamma_bio_over_jmax"]),
                    "theta_min_bio_deg": float(row["theta_min_bio_mean_deg"]),
                    "theta_reference": "biological point (site-averaged)",
                    "provenance": "rebuilt-repo",
                }
            )
    return rows


def load_protein_microdomain_anchor(root: Path) -> dict[str, str | float | int] | None:
    """Load the rebuilt protein microdomain anchor if present."""
    summary_path = root / "results" / "protein_microdomain_summary.json"
    if not summary_path.exists():
        return None

    with summary_path.open("r", encoding="utf-8") as fh:
        summary = json.load(fh)

    return {
        "system": "Protein mid-fold",
        "role": "secondary-payload",
        "status": "rebuilt-computation",
        "dimension": int(summary["dimension"]),
        "collapse": "35->5",
        "gamma_bio_over_jmax": float(summary["gamma_bio_over_jmax"]),
        "theta_min_bio_deg": float(summary["theta_min_bio_deg"]),
        "theta_reference": "biological point",
        "provenance": "rebuilt-repo",
    }


def load_neural_carrier_proxy_anchor(root: Path) -> dict[str, str | float | int] | None:
    """Load the rebuilt neural carrier proxy anchor if present."""
    csv_path = root / "results" / "neural_carrier_proxy.csv"
    if not csv_path.exists():
        return None

    with csv_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))

    if not rows:
        return None

    row = max(rows, key=lambda entry: float(entry["gamma_over_j"]))
    return {
        "system": "Carrier proxy",
        "role": "carrier-proxy",
        "status": "rebuilt-proxy",
        "dimension": 10,
        "collapse": "99->9",
        "gamma_bio_over_jmax": float(row["gamma_over_j"]),
        "theta_min_bio_deg": float(row["theta_min_deg"]),
        "theta_reference": "proxy point",
        "provenance": "rebuilt-repo",
    }


def anchor_rows(summary: dict) -> list[dict[str, str | float | int]]:
    """Build the cross-system anchor list used by the manuscript."""
    root = Path(__file__).resolve().parents[1]
    rows: list[dict[str, str | float | int]] = load_photosynthetic_anchors(root)
    rows.append(
        {
            "system": "Ion channel",
            "role": "primary-neural-payload",
            "status": "rebuilt-computation",
            "dimension": summary["dimension"],
            "collapse": "15->3",
            "gamma_bio_over_jmax": summary["gamma_bio_over_jmax"],
            "theta_min_bio_deg": summary["theta_min_bio_deg"],
            "theta_reference": "biological point",
            "provenance": "rebuilt-repo",
        }
    )
    protein_row = load_protein_microdomain_anchor(root)
    if protein_row is not None:
        rows.append(protein_row)
    carrier_row = load_neural_carrier_proxy_anchor(root)
    if carrier_row is not None:
        rows.append(carrier_row)
    return rows


def write_biological_anchor_outputs(summary: dict) -> None:
    """Write a compact anchor CSV and proof-of-principle figure."""
    root = Path(__file__).resolve().parents[1]
    results_dir = root / "results"
    figures_dir = root / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)

    rows = anchor_rows(summary)
    csv_path = results_dir / "biological_anchor_points.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "system",
                "role",
                "status",
                "dimension",
                "collapse",
                "gamma_bio_over_jmax",
                "theta_min_bio_deg",
                "theta_reference",
                "provenance",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(5.9, 3.6), constrained_layout=True)

    style = {
        "FMO": {"color": "#1f4e79", "marker": "o"},
        "PE545": {"color": "#d97904", "marker": "s"},
        "Ion channel": {"color": "#4c956c", "marker": "D"},
        "Protein mid-fold": {"color": "#6c757d", "marker": "^"},
        "Carrier proxy": {"color": "#7b2cbf", "marker": "P"},
    }

    ax.axhspan(87.0, 90.0, color="#e8f1f2", alpha=0.9, zorder=0)
    ax.axhline(90.0, color="black", ls=":", lw=1.0, zorder=1)

    for row in rows:
        system = str(row["system"])
        s = style[system]
        x_val = float(row["gamma_bio_over_jmax"])
        y_val = float(row["theta_min_bio_deg"])
        ax.scatter(
            [x_val],
            [y_val],
            s=58,
            marker=s["marker"],
            color=s["color"],
            edgecolor="white",
            linewidth=0.7,
            zorder=3,
        )
        ax.text(
            x_val * 1.06,
            y_val + 0.08,
            system,
            fontsize=8,
            color=s["color"],
        )

    ax.set_xscale("log")
    ax.set_xlim(0.8, 80.0)
    ax.set_ylim(86.9, 90.2)
    ax.set_xlabel(r"Biological ratio $\gamma_{\mathrm{bio}}/J_{\max}$")
    ax.set_ylabel(r"Minimum angle at bio point $\theta_{\min}^{\rm bio}$")
    ax.set_title("Payload anchors cluster left; carrier proxy sits deeper classical")
    ax.text(
        0.03,
        0.06,
        "Left cluster: rebuilt molecular payloads\nRight point: rebuilt neural carrier proxy",
        transform=ax.transAxes,
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#cccccc"},
    )

    fig.savefig(figures_dir / "biological_anchor_map.pdf", bbox_inches="tight")
    fig.savefig(figures_dir / "biological_anchor_map.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_carrier_payload_schematic() -> None:
    """Draw a conceptual carrier/payload schematic for the manuscript."""
    root = Path(__file__).resolve().parents[1]
    figures_dir = root / "figures"
    figures_dir.mkdir(exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.1), constrained_layout=True)

    x_positions = np.linspace(0.12, 0.88, 7)
    y_positions = np.linspace(0.16, 0.58, 4)
    payload_points = np.array([(x, y) for y in y_positions for x in x_positions], dtype=float)
    carrier_color = "#355070"
    payload_color = "#d17a22"
    muted = "#c8ccd0"

    configs = [
        {
            "ax": axes[0],
            "title": "High-$f$ local access",
            "subtitle": r"small $|E(f)|$, fast routing",
            "wave_cycles": 4.8,
            "radius": 0.18,
            "center": np.array([0.5, 0.38]),
        },
        {
            "ax": axes[1],
            "title": "Low-$f$ broad coordination",
            "subtitle": r"large $|E(f)|$, deep-classical carrier",
            "wave_cycles": 1.8,
            "radius": 0.34,
            "center": np.array([0.5, 0.38]),
        },
    ]

    for cfg in configs:
        ax = cfg["ax"]
        center = cfg["center"]
        radius = float(cfg["radius"])

        distances = np.linalg.norm(payload_points - center, axis=1)
        engaged = distances <= radius
        inactive = ~engaged

        ax.scatter(
            payload_points[inactive, 0],
            payload_points[inactive, 1],
            s=32,
            color=muted,
            edgecolor="white",
            linewidth=0.5,
            zorder=2,
        )
        ax.scatter(
            payload_points[engaged, 0],
            payload_points[engaged, 1],
            s=46,
            color=payload_color,
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )

        theta = np.linspace(0.0, 2.0 * np.pi, 300)
        ax.plot(
            center[0] + radius * np.cos(theta),
            center[1] + radius * np.sin(theta),
            color=carrier_color,
            ls="--",
            lw=1.2,
            alpha=0.7,
            zorder=1,
        )

        xs = np.linspace(0.02, 0.98, 500)
        ys = 0.84 + 0.055 * np.sin(2.0 * np.pi * cfg["wave_cycles"] * xs)
        ax.plot(xs, ys, color=carrier_color, lw=2.4, zorder=4)

        ax.annotate(
            "",
            xy=(center[0], center[1] + radius * 0.82),
            xytext=(center[0], 0.76),
            arrowprops={"arrowstyle": "->", "color": carrier_color, "lw": 1.5},
            zorder=4,
        )

        ax.text(0.04, 0.93, cfg["title"], transform=ax.transAxes, fontsize=11, weight="bold")
        ax.text(0.04, 0.86, cfg["subtitle"], transform=ax.transAxes, fontsize=9, color="#444444")
        ax.text(0.05, 0.77, r"carrier: $\chi \approx 0$", transform=ax.transAxes, fontsize=9, color=carrier_color)
        ax.text(0.05, 0.08, r"payload modules: $\chi > 0$", transform=ax.transAxes, fontsize=9, color=payload_color)
        ax.text(
            0.68,
            0.08,
            rf"$|E(f)|={int(engaged.sum())}$",
            transform=ax.transAxes,
            fontsize=9,
            color="#444444",
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "edgecolor": "#cccccc"},
        )

        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    fig.savefig(figures_dir / "carrier_payload_schematic.pdf", bbox_inches="tight")
    fig.savefig(figures_dir / "carrier_payload_schematic.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    summary, scan_rows = summarize_ion_channel()
    write_outputs(summary, scan_rows)
    write_figure(summary, scan_rows)
    write_ion_channel_sensitivity()
    write_ion_channel_topology_sensitivity()
    write_phase_window_scaling_outputs(summary)
    write_phase_window_sensitivity_outputs(summary)
    write_payload_benchmarks(summary)
    write_biological_anchor_outputs(summary)
    write_carrier_payload_schematic()

    print("Ion-channel payload summary")
    print("=" * 32)
    print(f"gamma_bio / Jmax : {summary['gamma_bio_over_jmax']:.3f}")
    print(f"eta(bio)         : {summary['eta_bio']:.6f}")
    print(f"theta_min(bio)   : {summary['theta_min_bio_deg']:.3f} deg")
    print(f"chi(bio)         : {summary['chi_bio']:.6e}")
    print(
        "sensitivity      : "
        f"theta in [{summary['sensitivity_theta_min_deg_min']:.3f}, {summary['sensitivity_theta_min_deg_max']:.3f}] deg, "
        f"chi in [{summary['sensitivity_chi_min']:.3e}, {summary['sensitivity_chi_max']:.3e}]"
    )
    print(
        "topology sweep   : "
        f"theta in [{summary['topology_theta_min_deg_min']:.3f}, {summary['topology_theta_min_deg_max']:.3f}] deg, "
        f"chi in [{summary['topology_chi_min']:.3e}, {summary['topology_chi_max']:.3e}]"
    )
    print(f"peak scan gamma/J: {summary['transport_scan_peak_gamma_over_jmax']:.3f}")
    print(f"peak scan eta    : {summary['transport_scan_peak_eta']:.6f}")
    print(f"peak interior?   : {summary['transport_scan_peak_interior']}")
    print("Wrote results/ion_channel_summary.json, results/ion_channel_scan.csv,")
    print("results/{ion_channel_sensitivity,ion_channel_topology_sensitivity,threshold_scaling_scan,threshold_scaling_sensitivity,ion_payload_benchmarks,biological_anchor_points}.csv,")
    print("and figures/{biological_anchor_map,carrier_payload_schematic,ion_channel_anchor,threshold_entrainment_scaling}.{pdf,png}")


if __name__ == "__main__":
    main()
