#!/usr/bin/env python3
"""Illustrative protein microdomain payload anchor for the rebuilt manuscript.

Outputs:
- results/protein_microdomain_summary.json
- results/protein_microdomain_scan.csv
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

import ion_channel_payload as core


PROTEIN_MICRODOMAIN_H_CM = np.array(
    [
        [20.0, 25.0, 12.0, 0.0, 0.0, 0.0],
        [25.0, 35.0, 18.0, 22.0, 0.0, 0.0],
        [12.0, 18.0, 42.0, 20.0, 10.0, 0.0],
        [0.0, 22.0, 20.0, 50.0, 24.0, 12.0],
        [0.0, 0.0, 10.0, 24.0, 38.0, 26.0],
        [0.0, 0.0, 0.0, 12.0, 26.0, 28.0],
    ],
    dtype=np.float64,
)
PROTEIN_MICRODOMAIN_INITIAL_SITE = 0
PROTEIN_MICRODOMAIN_TARGET_SITES = [5]
PROTEIN_MICRODOMAIN_KAPPA_PS = [0.5]
PROTEIN_MICRODOMAIN_BIO_GAMMA_CM = 200.0
PROTEIN_MICRODOMAIN_JMAX_CM = float(
    np.max(np.abs(PROTEIN_MICRODOMAIN_H_CM - np.diag(np.diag(PROTEIN_MICRODOMAIN_H_CM))))
)


def geometry_snapshot(gamma_cm: float) -> dict[str, float | list[float]]:
    """Geometry of the protein microdomain at a given dephasing rate."""
    rho0 = np.zeros((6, 6), dtype=complex)
    rho0[PROTEIN_MICRODOMAIN_INITIAL_SITE, PROTEIN_MICRODOMAIN_INITIAL_SITE] = 1.0
    rho = core.evolve_lindblad_sink(
        PROTEIN_MICRODOMAIN_H_CM,
        gamma_cm,
        PROTEIN_MICRODOMAIN_TARGET_SITES,
        PROTEIN_MICRODOMAIN_KAPPA_PS,
        rho0,
        t_final_ps=core.GEOMETRY_READOUT_PS,
        dt_fs=core.DT_FS,
    )
    rho_reg = core.regularize_density_matrix(rho)
    angles_rad = np.sort(core.principal_angles(rho_reg))
    angles_deg = [float(np.degrees(val)) for val in angles_rad]
    return {
        "theta_min_deg": angles_deg[0],
        "theta_spectrum_deg": angles_deg,
        "chi": float(core.non_classicality_index(angles_rad)),
        "trace_after_readout": float(np.real(np.trace(rho))),
    }


def summarize_protein_microdomain() -> tuple[dict, list[dict[str, float]]]:
    """Compute a reproducible protein microdomain anchor summary."""
    gamma_over_j_values = np.logspace(-2, np.log10(30.0), core.SCAN_POINTS)
    scan_rows: list[dict[str, float]] = []
    efficiencies: list[float] = []

    for gamma_over_j in gamma_over_j_values:
        gamma_cm = gamma_over_j * PROTEIN_MICRODOMAIN_JMAX_CM
        eta = core.compute_transport_efficiency(
            PROTEIN_MICRODOMAIN_H_CM,
            gamma_cm,
            PROTEIN_MICRODOMAIN_TARGET_SITES,
            PROTEIN_MICRODOMAIN_KAPPA_PS,
            PROTEIN_MICRODOMAIN_INITIAL_SITE,
            t_final_ps=core.TRANSPORT_WINDOW_PS,
            dt_fs=core.DT_FS,
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
    bio_geometry = geometry_snapshot(PROTEIN_MICRODOMAIN_BIO_GAMMA_CM)
    eta_bio = core.compute_transport_efficiency(
        PROTEIN_MICRODOMAIN_H_CM,
        PROTEIN_MICRODOMAIN_BIO_GAMMA_CM,
        PROTEIN_MICRODOMAIN_TARGET_SITES,
        PROTEIN_MICRODOMAIN_KAPPA_PS,
        PROTEIN_MICRODOMAIN_INITIAL_SITE,
        t_final_ps=core.TRANSPORT_WINDOW_PS,
        dt_fs=core.DT_FS,
    )

    summary = {
        "model": "Illustrative protein mid-fold microdomain",
        "dimension": 6,
        "initial_site": PROTEIN_MICRODOMAIN_INITIAL_SITE,
        "target_sites": PROTEIN_MICRODOMAIN_TARGET_SITES,
        "kappa_ps": PROTEIN_MICRODOMAIN_KAPPA_PS,
        "jmax_cm": PROTEIN_MICRODOMAIN_JMAX_CM,
        "gamma_bio_cm": PROTEIN_MICRODOMAIN_BIO_GAMMA_CM,
        "gamma_bio_over_jmax": PROTEIN_MICRODOMAIN_BIO_GAMMA_CM / PROTEIN_MICRODOMAIN_JMAX_CM,
        "eta_bio": float(eta_bio),
        "theta_min_bio_deg": bio_geometry["theta_min_deg"],
        "theta_spectrum_bio_deg": bio_geometry["theta_spectrum_deg"],
        "chi_bio": bio_geometry["chi"],
        "local_load_bio": float(core.local_nonclassical_load(float(bio_geometry["chi"]), 6)),
        "trace_after_bio_readout": bio_geometry["trace_after_readout"],
        "transport_scan_peak_gamma_cm": peak_row["gamma_cm"],
        "transport_scan_peak_gamma_over_jmax": peak_row["gamma_over_j"],
        "transport_scan_peak_eta": peak_row["eta"],
        "transport_scan_peak_interior": peak_is_interior,
        "geometry_readout_ps": core.GEOMETRY_READOUT_PS,
        "transport_window_ps": core.TRANSPORT_WINDOW_PS,
        "scan_points": len(scan_rows),
    }
    return summary, scan_rows


def write_outputs(summary: dict, scan_rows: list[dict[str, float]]) -> None:
    """Write JSON summary and scan CSV."""
    root = Path(__file__).resolve().parents[1]
    results_dir = root / "results"
    results_dir.mkdir(exist_ok=True)

    summary_path = results_dir / "protein_microdomain_summary.json"
    scan_path = results_dir / "protein_microdomain_scan.csv"

    with summary_path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
        fh.write("\n")

    with scan_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["gamma_over_j", "gamma_cm", "eta"])
        writer.writeheader()
        writer.writerows(scan_rows)


def main() -> None:
    summary, scan_rows = summarize_protein_microdomain()
    write_outputs(summary, scan_rows)

    print("Protein microdomain summary")
    print("=" * 33)
    print(f"gamma_bio / Jmax : {summary['gamma_bio_over_jmax']:.3f}")
    print(f"eta(bio)         : {summary['eta_bio']:.6f}")
    print(f"theta_min(bio)   : {summary['theta_min_bio_deg']:.3f} deg")
    print(f"chi(bio)         : {summary['chi_bio']:.6e}")
    print(f"peak scan gamma/J: {summary['transport_scan_peak_gamma_over_jmax']:.3f}")
    print(f"peak scan eta    : {summary['transport_scan_peak_eta']:.6f}")
    print(f"peak interior?   : {summary['transport_scan_peak_interior']}")
    print("Wrote results/protein_microdomain_summary.json and results/protein_microdomain_scan.csv")


if __name__ == "__main__":
    main()
