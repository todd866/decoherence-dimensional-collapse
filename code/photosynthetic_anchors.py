#!/usr/bin/env python3
"""Reduced photosynthesis proof-of-principle computations for the rebuilt paper.

This script reintroduces only the anchor quantities the rebuilt draft needs:

- biological-point minimum principal angle for FMO and PE545
- transport-optimum dephasing ratio from the same sink-inclusive dynamics
- minimum principal angle at the transport optimum

Outputs:
- results/photosynthetic_anchor_points.csv
- results/photosynthetic_anchor_runs.csv
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from ion_channel_payload import (
    compute_transport_efficiency,
    evolve_lindblad_sink,
    local_nonclassical_load,
    non_classicality_index,
    principal_angles,
    regularize_density_matrix,
)

PHOTO_DT_FS = 5.0


H_FMO = np.array(
    [
        [200.0, -87.7, 5.5, -5.9, 6.7, -13.7, -9.9],
        [-87.7, 320.0, 30.8, 8.2, 0.7, 11.8, 4.3],
        [5.5, 30.8, 0.0, -53.5, -2.2, -9.6, 6.0],
        [-5.9, 8.2, -53.5, 110.0, -70.7, -17.0, -63.3],
        [6.7, 0.7, -2.2, -70.7, 270.0, 81.1, -1.3],
        [-13.7, 11.8, -9.6, -17.0, 81.1, 420.0, 39.7],
        [-9.9, 4.3, 6.0, -63.3, -1.3, 39.7, 230.0],
    ],
    dtype=np.float64,
)

H_PE545 = np.array(
    [
        [18532.0, 1.0, -37.0, 37.0, 23.0, 92.0, -16.0, 12.0],
        [1.0, 18008.0, 4.0, -11.0, 33.0, -39.0, -46.0, 3.0],
        [-37.0, 4.0, 17973.0, 45.0, 3.0, 2.0, -11.0, 34.0],
        [37.0, -11.0, 45.0, 18040.0, -7.0, -17.0, -3.0, 6.0],
        [23.0, 33.0, 3.0, -7.0, 18711.0, 18.0, 7.0, 6.0],
        [92.0, -39.0, 2.0, -17.0, 18.0, 19574.0, 40.0, 26.0],
        [-16.0, -46.0, -11.0, -3.0, 7.0, 40.0, 19050.0, 7.0],
        [12.0, 3.0, 34.0, 6.0, 6.0, 26.0, 7.0, 18960.0],
    ],
    dtype=np.float64,
)


SYSTEMS = [
    {
        "system": "FMO",
        "dimension": 7,
        "collapse": "48->6",
        "hamiltonian_cm": H_FMO,
        "initial_sites": [0, 5],
        "target_sites": [2],
        "kappa_ps": [1.0],
        "gamma_scan_cm": np.logspace(np.log10(5.0), np.log10(500.0), 15),
        "gamma_bio_cm": 100.0,
        "transport_window_ps": 15.0,
        "geometry_readout_ps": 5.0,
    },
    {
        "system": "PE545",
        "dimension": 8,
        "collapse": "63->7",
        "hamiltonian_cm": H_PE545,
        "initial_sites": [0, 5],
        "target_sites": [1, 2],
        "kappa_ps": [0.5, 0.5],
        "gamma_scan_cm": np.logspace(np.log10(50.0), np.log10(1500.0), 15),
        "gamma_bio_cm": 725.0,
        "transport_window_ps": 30.0,
        "geometry_readout_ps": 5.0,
    },
]


def jmax_cm(hamiltonian_cm: np.ndarray) -> float:
    off_diag = hamiltonian_cm - np.diag(np.diag(hamiltonian_cm))
    return float(np.max(np.abs(off_diag)))


def theta_snapshot(
    hamiltonian_cm: np.ndarray,
    gamma_cm: float,
    initial_site: int,
    target_sites: list[int],
    kappa_ps: list[float],
    geometry_readout_ps: float,
) -> tuple[float, float]:
    rho0 = np.zeros((hamiltonian_cm.shape[0], hamiltonian_cm.shape[0]), dtype=complex)
    rho0[initial_site, initial_site] = 1.0
    rho = evolve_lindblad_sink(
        hamiltonian_cm,
        gamma_cm,
        target_sites,
        kappa_ps,
        rho0,
        t_final_ps=geometry_readout_ps,
        dt_fs=PHOTO_DT_FS,
    )
    rho_reg = regularize_density_matrix(rho)
    angles = np.sort(principal_angles(rho_reg))
    return float(np.degrees(angles[0])), float(non_classicality_index(angles))


def evaluate_system(spec: dict) -> tuple[list[dict[str, float | int | str]], dict[str, float | str | int]]:
    hamiltonian_cm = spec["hamiltonian_cm"]
    gamma_scan_cm = spec["gamma_scan_cm"]
    jmax = jmax_cm(hamiltonian_cm)
    per_run_rows: list[dict[str, float | int | str]] = []

    bio_thetas = []
    bio_chis = []
    opt_thetas = []
    opt_chis = []
    opt_goj = []
    opt_eta = []

    for initial_site in spec["initial_sites"]:
        eta_scan = []
        for gamma_cm in gamma_scan_cm:
            eta_scan.append(
                compute_transport_efficiency(
                    hamiltonian_cm,
                    float(gamma_cm),
                    spec["target_sites"],
                    spec["kappa_ps"],
                    initial_site,
                    t_final_ps=spec["transport_window_ps"],
                    dt_fs=PHOTO_DT_FS,
                )
            )
        idx_peak = int(np.nanargmax(eta_scan))
        gamma_opt_cm = float(gamma_scan_cm[idx_peak])
        theta_bio_deg, chi_bio = theta_snapshot(
            hamiltonian_cm,
            spec["gamma_bio_cm"],
            initial_site,
            spec["target_sites"],
            spec["kappa_ps"],
            spec["geometry_readout_ps"],
        )
        theta_opt_deg, chi_opt = theta_snapshot(
            hamiltonian_cm,
            gamma_opt_cm,
            initial_site,
            spec["target_sites"],
            spec["kappa_ps"],
            spec["geometry_readout_ps"],
        )

        bio_thetas.append(theta_bio_deg)
        bio_chis.append(chi_bio)
        opt_thetas.append(theta_opt_deg)
        opt_chis.append(chi_opt)
        opt_goj.append(gamma_opt_cm / jmax)
        opt_eta.append(float(eta_scan[idx_peak]))

        per_run_rows.append(
            {
                "system": spec["system"],
                "initial_site": initial_site,
                "jmax_cm": jmax,
                "gamma_bio_cm": spec["gamma_bio_cm"],
                "gamma_bio_over_jmax": spec["gamma_bio_cm"] / jmax,
                "theta_min_bio_deg": theta_bio_deg,
                "chi_bio": chi_bio,
                "gamma_opt_cm": gamma_opt_cm,
                "gamma_opt_over_jmax": gamma_opt_cm / jmax,
                "theta_min_opt_deg": theta_opt_deg,
                "chi_opt": chi_opt,
                "eta_opt": float(eta_scan[idx_peak]),
            }
        )

    summary_row = {
        "system": spec["system"],
        "dimension": spec["dimension"],
        "collapse": spec["collapse"],
        "source": "rebuilt-validated",
        "n_initial_sites": len(spec["initial_sites"]),
        "jmax_cm": jmax,
        "gamma_bio_over_jmax": spec["gamma_bio_cm"] / jmax,
        "theta_min_bio_mean_deg": float(np.mean(bio_thetas)),
        "theta_min_bio_min_deg": float(np.min(bio_thetas)),
        "theta_min_bio_max_deg": float(np.max(bio_thetas)),
        "chi_bio_mean": float(np.mean(bio_chis)),
        "load_bio_mean": float(local_nonclassical_load(float(np.mean(bio_chis)), spec["dimension"])),
        "gamma_opt_over_jmax_mean": float(np.mean(opt_goj)),
        "gamma_opt_over_jmax_min": float(np.min(opt_goj)),
        "gamma_opt_over_jmax_max": float(np.max(opt_goj)),
        "theta_min_opt_mean_deg": float(np.mean(opt_thetas)),
        "theta_min_opt_min_deg": float(np.min(opt_thetas)),
        "theta_min_opt_max_deg": float(np.max(opt_thetas)),
        "chi_opt_mean": float(np.mean(opt_chis)),
        "load_opt_mean": float(local_nonclassical_load(float(np.mean(opt_chis)), spec["dimension"])),
        "eta_opt_mean": float(np.mean(opt_eta)),
        "eta_opt_min": float(np.min(opt_eta)),
        "eta_opt_max": float(np.max(opt_eta)),
    }
    return per_run_rows, summary_row


def write_outputs() -> None:
    root = Path(__file__).resolve().parents[1]
    results_dir = root / "results"
    results_dir.mkdir(exist_ok=True)

    all_runs: list[dict[str, float | int | str]] = []
    summaries: list[dict[str, float | int | str]] = []

    for spec in SYSTEMS:
        runs, summary = evaluate_system(spec)
        all_runs.extend(runs)
        summaries.append(summary)

    runs_path = results_dir / "photosynthetic_anchor_runs.csv"
    points_path = results_dir / "photosynthetic_anchor_points.csv"

    with runs_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "system",
                "initial_site",
                "jmax_cm",
                "gamma_bio_cm",
                "gamma_bio_over_jmax",
                "theta_min_bio_deg",
                "chi_bio",
                "gamma_opt_cm",
                "gamma_opt_over_jmax",
                "theta_min_opt_deg",
                "chi_opt",
                "eta_opt",
            ],
        )
        writer.writeheader()
        writer.writerows(all_runs)

    with points_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "system",
                "dimension",
                "collapse",
                "source",
                "n_initial_sites",
                "jmax_cm",
                "gamma_bio_over_jmax",
                "theta_min_bio_mean_deg",
                "theta_min_bio_min_deg",
                "theta_min_bio_max_deg",
                "chi_bio_mean",
                "load_bio_mean",
                "gamma_opt_over_jmax_mean",
                "gamma_opt_over_jmax_min",
                "gamma_opt_over_jmax_max",
                "theta_min_opt_mean_deg",
                "theta_min_opt_min_deg",
                "theta_min_opt_max_deg",
                "chi_opt_mean",
                "load_opt_mean",
                "eta_opt_mean",
                "eta_opt_min",
                "eta_opt_max",
            ],
        )
        writer.writeheader()
        writer.writerows(summaries)

    print("Photosynthetic anchors")
    print("=" * 24)
    for row in summaries:
        print(
            f"{row['system']}: gamma_bio/J={row['gamma_bio_over_jmax']:.3f}, "
            f"theta_bio={row['theta_min_bio_mean_deg']:.3f}°, "
            f"gamma_opt/J={row['gamma_opt_over_jmax_mean']:.3f}, "
            f"theta_opt={row['theta_min_opt_mean_deg']:.3f}°, "
            f"eta_opt={row['eta_opt_mean']:.3f}"
        )
    print("Wrote results/photosynthetic_anchor_points.csv and results/photosynthetic_anchor_runs.csv")


if __name__ == "__main__":
    write_outputs()
