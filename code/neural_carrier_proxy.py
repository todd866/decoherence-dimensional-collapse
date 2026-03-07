#!/usr/bin/env python3
"""Reduced neural carrier proxies for the rebuilt manuscript.

This script computes only the carrier-side geometry the rebuilt paper needs:

- a synthetic gamma-band cortical microcircuit
- a reduced laminar E/I carrier proxy
- principal-angle geometry at a few proxy dephasing ratios
- explicit high-dephasing proxy points showing carrier classicality

Outputs:
- results/neural_carrier_proxy.csv
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from ion_channel_payload import (
    evolve_lindblad_sink,
    non_classicality_index,
    principal_angles,
    regularize_density_matrix,
)


NEURAL_GAMMA_FREQS_HZ = np.array([30, 33, 36, 40, 44, 48, 53, 58, 64, 70], dtype=np.float64)
NEURAL_GAMMA_SCALE = 50.0  # cm^-1 per Hz
NEURAL_J_NN_HZ = 5.0
NEURAL_J_SHORT_HZ = 2.0
NEURAL_TARGET_SITES = [8, 9]
NEURAL_KAPPA_PS = [0.5, 0.5]
NEURAL_INITIAL_SITE = 0
NEURAL_READOUT_PS = 5.0
NEURAL_DT_FS = 1.0
PROXY_GAMMA_OVER_J = [1.0, 10.0, 50.0]

LAMINAR_FREQS_HZ = np.array([48, 60, 44, 56, 34, 46, 28, 40], dtype=np.float64)
LAMINAR_SCALE = 50.0  # cm^-1 per Hz
LAMINAR_WITHIN_LAYER_HZ = [5.0, 5.0, 4.5, 4.0]
LAMINAR_FEEDFORWARD_HZ = {
    (2, 0): 3.5,  # L4E -> L23E
    (0, 4): 2.5,  # L23E -> L5E
    (4, 6): 3.0,  # L5E -> L6E
    (6, 2): 2.0,  # L6E -> L4E
}
LAMINAR_INHIBITORY_HZ = {
    (1, 3): 2.0,
    (3, 5): 2.0,
    (5, 7): 1.5,
}
LAMINAR_CROSS_EI_HZ = {
    (0, 3): 1.0,
    (2, 5): 1.0,
    (4, 7): 1.0,
    (6, 1): 0.8,
}
LAMINAR_TARGET_SITES = [6, 7]
LAMINAR_KAPPA_PS = [0.5, 0.5]
LAMINAR_INITIAL_SITE = 2


def build_neural_gamma_hamiltonian() -> np.ndarray:
    """10x10 synthetic gamma-band microcircuit Hamiltonian in cm^-1."""
    d = 10
    h_cm = np.diag(NEURAL_GAMMA_FREQS_HZ * NEURAL_GAMMA_SCALE)
    for idx in range(d):
        nxt = (idx + 1) % d
        h_cm[idx, nxt] = NEURAL_J_NN_HZ * NEURAL_GAMMA_SCALE
        h_cm[nxt, idx] = NEURAL_J_NN_HZ * NEURAL_GAMMA_SCALE
    for i, j in [(1, 6), (3, 8)]:
        h_cm[i, j] = NEURAL_J_SHORT_HZ * NEURAL_GAMMA_SCALE
        h_cm[j, i] = NEURAL_J_SHORT_HZ * NEURAL_GAMMA_SCALE
    return h_cm


def build_laminar_ei_hamiltonian() -> np.ndarray:
    """8x8 layered E/I carrier proxy in cm^-1."""
    h_cm = np.diag(LAMINAR_FREQS_HZ * LAMINAR_SCALE)

    for layer_idx, (exc_idx, inh_idx) in enumerate([(0, 1), (2, 3), (4, 5), (6, 7)]):
        coupling = LAMINAR_WITHIN_LAYER_HZ[layer_idx] * LAMINAR_SCALE
        h_cm[exc_idx, inh_idx] = coupling
        h_cm[inh_idx, exc_idx] = coupling

    for (src, dst), coupling_hz in LAMINAR_FEEDFORWARD_HZ.items():
        coupling = coupling_hz * LAMINAR_SCALE
        h_cm[src, dst] = coupling
        h_cm[dst, src] = coupling

    for (src, dst), coupling_hz in LAMINAR_INHIBITORY_HZ.items():
        coupling = coupling_hz * LAMINAR_SCALE
        h_cm[src, dst] = coupling
        h_cm[dst, src] = coupling

    for (src, dst), coupling_hz in LAMINAR_CROSS_EI_HZ.items():
        coupling = coupling_hz * LAMINAR_SCALE
        h_cm[src, dst] = coupling
        h_cm[dst, src] = coupling

    return h_cm


def compute_proxy_rows(
    *,
    model: str,
    h_cm: np.ndarray,
    target_sites: list[int],
    kappa_ps: list[float],
    initial_site: int,
) -> list[dict[str, float | str]]:
    """Compute carrier-proxy geometry rows for one model."""
    jmax_cm = float(np.max(np.abs(h_cm - np.diag(np.diag(h_cm)))))
    rho0 = np.zeros((h_cm.shape[0], h_cm.shape[0]), dtype=complex)
    rho0[initial_site, initial_site] = 1.0

    rows: list[dict[str, float | str]] = []
    for gamma_over_j in PROXY_GAMMA_OVER_J:
        gamma_cm = gamma_over_j * jmax_cm
        rho = evolve_lindblad_sink(
            h_cm,
            gamma_cm,
            target_sites,
            kappa_ps,
            rho0,
            t_final_ps=NEURAL_READOUT_PS,
            dt_fs=NEURAL_DT_FS,
        )
        rho_reg = regularize_density_matrix(rho)
        angles = np.sort(principal_angles(rho_reg))
        theta_min_deg = float(np.degrees(angles[0]))
        rows.append(
            {
                "model": model,
                "dimension": int(h_cm.shape[0]),
                "gamma_over_j": float(gamma_over_j),
                "gamma_cm": float(gamma_cm),
                "jmax_cm": jmax_cm,
                "theta_min_deg": theta_min_deg,
                "chi": float(non_classicality_index(angles)),
                "trace_after_readout": float(np.real(np.trace(rho))),
            }
        )
    return rows


def write_proxy_outputs() -> None:
    root = Path(__file__).resolve().parents[1]
    results_dir = root / "results"
    results_dir.mkdir(exist_ok=True)

    rows = []
    rows.extend(
        compute_proxy_rows(
            model="Synthetic gamma",
            h_cm=build_neural_gamma_hamiltonian(),
            target_sites=NEURAL_TARGET_SITES,
            kappa_ps=NEURAL_KAPPA_PS,
            initial_site=NEURAL_INITIAL_SITE,
        )
    )
    rows.extend(
        compute_proxy_rows(
            model="Laminar E/I",
            h_cm=build_laminar_ei_hamiltonian(),
            target_sites=LAMINAR_TARGET_SITES,
            kappa_ps=LAMINAR_KAPPA_PS,
            initial_site=LAMINAR_INITIAL_SITE,
        )
    )

    csv_path = results_dir / "neural_carrier_proxy.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "model",
                "dimension",
                "gamma_over_j",
                "gamma_cm",
                "jmax_cm",
                "theta_min_deg",
                "chi",
                "trace_after_readout",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print("Neural carrier proxy")
    print("=" * 20)
    for model in sorted({str(row["model"]) for row in rows}):
        print(model)
        for row in [r for r in rows if r["model"] == model]:
            print(
                f"  gamma/J={row['gamma_over_j']:.1f}: "
                f"theta_min={row['theta_min_deg']:.3f}°, chi={row['chi']:.3e}"
            )
    print("Wrote results/neural_carrier_proxy.csv")


if __name__ == "__main__":
    write_proxy_outputs()
