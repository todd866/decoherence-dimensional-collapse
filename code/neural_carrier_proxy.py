#!/usr/bin/env python3
"""Reduced neural carrier proxy for the rebuilt manuscript.

This script computes only the carrier-side geometry the rebuilt paper needs:

- a synthetic gamma-band cortical microcircuit
- principal-angle geometry at a few proxy dephasing ratios
- one explicit high-dephasing proxy point showing carrier classicality

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


def write_proxy_outputs() -> None:
    root = Path(__file__).resolve().parents[1]
    results_dir = root / "results"
    results_dir.mkdir(exist_ok=True)

    h_cm = build_neural_gamma_hamiltonian()
    jmax_cm = float(np.max(np.abs(h_cm - np.diag(np.diag(h_cm)))))
    rho0 = np.zeros((h_cm.shape[0], h_cm.shape[0]), dtype=complex)
    rho0[NEURAL_INITIAL_SITE, NEURAL_INITIAL_SITE] = 1.0

    rows = []
    for gamma_over_j in PROXY_GAMMA_OVER_J:
        gamma_cm = gamma_over_j * jmax_cm
        rho = evolve_lindblad_sink(
            h_cm,
            gamma_cm,
            NEURAL_TARGET_SITES,
            NEURAL_KAPPA_PS,
            rho0,
            t_final_ps=NEURAL_READOUT_PS,
            dt_fs=NEURAL_DT_FS,
        )
        rho_reg = regularize_density_matrix(rho)
        angles = np.sort(principal_angles(rho_reg))
        theta_min_deg = float(np.degrees(angles[0]))
        rows.append(
            {
                "gamma_over_j": float(gamma_over_j),
                "gamma_cm": float(gamma_cm),
                "jmax_cm": jmax_cm,
                "theta_min_deg": theta_min_deg,
                "chi": float(non_classicality_index(angles)),
                "trace_after_readout": float(np.real(np.trace(rho))),
            }
        )

    csv_path = results_dir / "neural_carrier_proxy.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
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
    for row in rows:
        print(
            f"gamma/J={row['gamma_over_j']:.1f}: "
            f"theta_min={row['theta_min_deg']:.3f}°, chi={row['chi']:.3e}"
        )
    print("Wrote results/neural_carrier_proxy.csv")


if __name__ == "__main__":
    write_proxy_outputs()
