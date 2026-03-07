#!/usr/bin/env python3
"""Reduced neural carrier proxies for the rebuilt manuscript.

This script computes only the carrier-side geometry the rebuilt paper needs:

- a synthetic gamma-band cortical microcircuit
- a reduced Potjans-Diesmann cortical-column carrier proxy
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

PD_FREQS_HZ = np.array([3.0, 6.0, 5.0, 8.0, 7.0, 9.0, 1.0, 4.0], dtype=np.float64)
PD_TARGET_SITES = [6, 7]
PD_KAPPA_PS = [0.5, 0.5]
PD_INITIAL_SITE = 2


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


def build_potjans_diesmann_hamiltonian() -> np.ndarray:
    """8x8 Potjans-Diesmann-inspired cortical column Hamiltonian in cm^-1.

    Rows/cols are L2/3e, L2/3i, L4e, L4i, L5e, L5i, L6e, L6i.
    Connection probabilities are taken from Table 5 of Potjans & Diesmann 2014
    and converted into a symmetric effective coupling matrix.
    """
    d = 8
    c = np.array([
        [0.101, 0.169, 0.044, 0.082, 0.032, 0.0,   0.008, 0.0  ],
        [0.135, 0.137, 0.032, 0.052, 0.075, 0.0,   0.004, 0.0  ],
        [0.008, 0.006, 0.050, 0.135, 0.007, 0.0003,0.045, 0.0  ],
        [0.069, 0.003, 0.079, 0.160, 0.003, 0.0,   0.106, 0.0  ],
        [0.100, 0.062, 0.051, 0.006, 0.083, 0.373, 0.020, 0.0  ],
        [0.055, 0.027, 0.026, 0.002, 0.060, 0.316, 0.009, 0.0  ],
        [0.016, 0.007, 0.021, 0.017, 0.057, 0.020, 0.040, 0.225],
        [0.036, 0.001, 0.003, 0.001, 0.028, 0.008, 0.066, 0.144],
    ], dtype=np.float64)
    j_raw = np.sqrt(c * c.T)
    signs = np.ones(d)
    signs[[1, 3, 5, 7]] = -1.0
    sign_matrix = np.outer(signs, signs)
    j_eff = j_raw * sign_matrix
    np.fill_diagonal(j_eff, 0.0)
    coupling_scale = 200.0 / np.max(np.abs(j_eff))
    return np.diag(PD_FREQS_HZ * NEURAL_GAMMA_SCALE) + j_eff * coupling_scale


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
            model="Potjans-Diesmann",
            h_cm=build_potjans_diesmann_hamiltonian(),
            target_sites=PD_TARGET_SITES,
            kappa_ps=PD_KAPPA_PS,
            initial_site=PD_INITIAL_SITE,
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
