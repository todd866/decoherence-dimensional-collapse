#!/usr/bin/env python3
"""
Generate all figures for:
"Information Geometry of the Quantum-Classical Transition: From Photosynthetic Excitons to Neural Oscillations"

Produces:
  figures/fig1_bloch_collapse.{pdf,png}     — Bloch sphere collapse
  figures/fig2_fmo_graph.{pdf,png}          — FMO coupling graph
  figures/fig3_bures_angle.{pdf,png}        — Bures angle across Bloch disk
  figures/fig4_fmo_collapse.{pdf,png}       — FMO transport + Bures angles (dual panel)
  figures/fig5_qudit_angles.{pdf,png}       — Qudit principal angles (d=3,5,7)
  figures/fig6_robustness.{pdf,png}         — Angle spectrum + sink-rate robustness
  figures/fig7_pe545.{pdf,png}              — PE545 transport + principal angles
  figures/fig8_neural_comparison.{pdf,png}  — Cross-scale quantum-classical boundary

Usage:
  cd physics/70_decoherence_dimensional_collapse
  python3 code/generate_figures.py
"""

import os
import sys
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import expm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import FancyArrowPatch, Patch
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import matplotlib.patches as mpatches

# Add code directory to path for FMO imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fmo_analysis import (
    H_FMO, D_FMO, CM_TO_RADFS, COUPLING_THRESHOLD,
    fmo_coupling_graph,
    evolve_lindblad, evolve_lindblad_sink, principal_angles, scan_efficiency,
    compute_transport_efficiency, l1_coherence, purity,
    H_PE545, compute_transport_generic, evolve_lindblad_sink_generic, COMPLEXES,
    H_NEURAL_GAMMA, H_NEURAL_PD, NEURAL_GAMMA_JMAX_CM, NEURAL_PD_JMAX_CM,
    NEURAL_BIO_GAMMA_OVER_JMAX, _NEURAL_GAMMA_SCALE,
    H_CHROMATIN, CHROMATIN_JMAX_CM, CHROMATIN_BIO_GAMMA_CM,
    CHROMATIN_BIO_GAMMA_OVER_JMAX,
    _CHROMATIN_SITE_ENERGIES_CM, _CHROMATIN_J_NN, _CHROMATIN_J_LR,
)

# ── Global style ─────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "text.usetex": False,
    "mathtext.fontset": "cm",
})

# Professional colour palette
PALETTE = {
    "blue":    "#2c5f8a",
    "orange":  "#d4762c",
    "green":   "#3a7d44",
    "red":     "#b83232",
    "purple":  "#7b4f9e",
    "grey":    "#6b6b6b",
}
PAL = list(PALETTE.values())

FIGDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")
os.makedirs(FIGDIR, exist_ok=True)


def save(fig, stem):
    """Save figure as PDF and PNG."""
    fig.savefig(os.path.join(FIGDIR, f"{stem}.pdf"), format="pdf")
    fig.savefig(os.path.join(FIGDIR, f"{stem}.png"), format="png", dpi=150)
    plt.close(fig)
    print(f"  saved {stem}.pdf / .png")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 1 — Bloch-sphere collapse under dephasing
# ══════════════════════════════════════════════════════════════════════════════

def fig1_bloch_collapse():
    print("Figure 1: Bloch-sphere collapse ...")

    gamma_vals = [0.0, 0.15, 0.35, 0.5]
    cmap = matplotlib.colormaps.get_cmap("coolwarm").resampled(len(gamma_vals))
    colours = [cmap(i) for i in range(len(gamma_vals))]

    fig = plt.figure(figsize=(7.0, 3.2))

    # Left panel: 3-D Bloch sphere
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(0, np.pi, 25)
    sx = np.outer(np.cos(u), np.sin(v))
    sy = np.outer(np.sin(u), np.sin(v))
    sz = np.outer(np.ones_like(u), np.cos(v))

    for idx, gamma in enumerate(gamma_vals):
        shrink = max(1.0 - 2.0 * gamma, 0.0)
        x = shrink * sx
        y = shrink * sy
        z = sz
        alpha = 0.18 if idx < len(gamma_vals) - 1 else 0.35
        ax3d.plot_surface(
            x, y, z, color=colours[idx], alpha=alpha,
            linewidth=0.3, edgecolor=(*colours[idx][:3], 0.25), shade=False,
        )

    ax3d.plot([0, 0], [0, 0], [-1.05, 1.05], color="k", lw=0.6, ls="--")
    ax3d.set_xlim(-1.1, 1.1); ax3d.set_ylim(-1.1, 1.1); ax3d.set_zlim(-1.1, 1.1)
    ax3d.set_xlabel(r"$r_x$", labelpad=1)
    ax3d.set_ylabel(r"$r_y$", labelpad=1)
    ax3d.set_zlabel(r"$r_z$", labelpad=1)
    ax3d.set_title("Accessible state space", fontsize=11, pad=2)
    ax3d.view_init(elev=20, azim=-55)
    ax3d.tick_params(axis="both", labelsize=7, pad=0)

    legend_patches = [
        Patch(facecolor=colours[i], alpha=0.5, label=rf"$\gamma = {gamma_vals[i]:.2f}$")
        for i in range(len(gamma_vals))
    ]
    ax3d.legend(handles=legend_patches, loc="upper left", fontsize=8,
                framealpha=0.85, handlelength=1.0, borderpad=0.3)

    # Right panel: D_eff trajectory
    ax = fig.add_subplot(1, 2, 2)
    gamma = np.linspace(0, 0.5, 500)
    shrink = (1.0 - 2.0 * gamma)
    f_phi = shrink ** 2
    Deff = (1.0 + f_phi) ** 2 / (1.0 + f_phi ** 2)
    ax.plot(gamma, Deff, color=PAL[0], lw=2.0)

    for idx, gv in enumerate(gamma_vals):
        s = max(1.0 - 2.0 * gv, 0.0)
        fp = s ** 2
        de = (1.0 + fp) ** 2 / (1.0 + fp ** 2)
        ax.plot(gv, de, "o", color=colours[idx], ms=7, zorder=5,
                markeredgecolor="white", markeredgewidth=0.8)

    ax.set_xlabel(r"Dephasing parameter $\gamma$")
    ax.set_ylabel(r"$D_{\mathrm{eff}}^Q$")
    ax.set_xlim(0, 0.5); ax.set_ylim(0.9, 2.05)
    ax.set_title(r"Effective dimension ($\theta = \pi/2$)", fontsize=11)
    ax.axhline(1.0, color=PAL[5], ls="--", lw=0.8, alpha=0.6)
    ax.text(0.42, 1.03, "Classical", fontsize=8, color=PAL[5])

    fig.tight_layout(w_pad=2.5)
    save(fig, "fig1_bloch_collapse")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 — FMO coupling graph
# ══════════════════════════════════════════════════════════════════════════════

def fig2_fmo_graph():
    print("Figure 2: FMO coupling graph ...")

    # Node positions (approximate spatial layout of FMO sites)
    positions = {
        1: (0.0,  1.5),
        2: (1.2,  1.0),
        3: (1.5, -0.3),
        4: (0.5, -1.2),
        5: (-0.8, -0.8),
        6: (-1.3,  0.3),
        7: (0.3, -0.2),
    }

    adj = fmo_coupling_graph()

    fig, ax = plt.subplots(figsize=(5.5, 5.0))

    # Draw significant coupling edges (|J| > 5 cm^-1)
    for i in range(D_FMO):
        for j in range(i + 1, D_FMO):
            x1, y1 = positions[i + 1]
            x2, y2 = positions[j + 1]
            coupling = abs(H_FMO[i, j])

            if adj[i, j]:
                lw = 0.5 + 2.5 * (coupling / 90.0)
                ax.plot([x1, x2], [y1, y2], color=PAL[0], lw=lw,
                        alpha=0.6, zorder=1)
            else:
                # Weak coupling edges (dotted)
                ax.plot([x1, x2], [y1, y2], color=PAL[5], lw=0.5,
                        ls=":", alpha=0.3, zorder=0)

    # Draw nodes
    for site, (x, y) in positions.items():
        circle = plt.Circle((x, y), 0.22, color=PAL[0], ec="white", lw=2,
                             zorder=3)
        ax.add_patch(circle)
        ax.text(x, y, str(site), ha="center", va="center", fontsize=11,
                fontweight="bold", color="white", zorder=4)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=PAL[0], lw=2,
               label=r"$|J_{kl}| > 5\,\mathrm{cm}^{-1}$"),
        Line2D([0], [0], color=PAL[5], lw=1, ls=":",
               label=r"$|J_{kl}| < 5\,\mathrm{cm}^{-1}$"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=9,
              framealpha=0.9)

    ax.set_xlim(-2.0, 2.3)
    ax.set_ylim(-1.8, 2.2)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.tight_layout()
    save(fig, "fig2_fmo_graph")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 3 — Bures separation angle across the Bloch disk
# ══════════════════════════════════════════════════════════════════════════════

def fig3_bures_angle():
    print("Figure 3: Bures separation angle ...")

    n = 500
    rperp = np.linspace(0, 1, n)
    rz = np.linspace(-1, 1, n)
    Rp, Rz = np.meshgrid(rperp, rz)
    R2 = Rp**2 + Rz**2

    mask = R2 >= 1.0
    denom = np.clip((1.0 - Rp**2) * (1.0 - Rz**2), 1e-15, None)
    sin2 = np.clip((1.0 - R2) / denom, 0.0, 1.0)
    theta = np.degrees(np.arcsin(np.sqrt(sin2)))
    theta[mask] = np.nan

    fig, ax = plt.subplots(figsize=(5.5, 4.5))

    im = ax.pcolormesh(
        rperp, rz, theta, cmap="RdYlBu", shading="auto",
        vmin=0, vmax=90, rasterized=True,
    )

    phi = np.linspace(-np.pi / 2, np.pi / 2, 300)
    ax.plot(np.cos(phi), np.sin(phi), "k-", lw=1.2, alpha=0.6)

    traj_rz = [0.3, 0.5, 0.7, 0.9]
    for rz_val in traj_rz:
        rp_max = np.sqrt(1.0 - rz_val**2) * 0.92
        ax.annotate("", xy=(0.02, rz_val), xytext=(rp_max, rz_val),
                     arrowprops=dict(arrowstyle="->", color="k", lw=1.0, linestyle="dashed"))
        ax.annotate("", xy=(0.02, -rz_val), xytext=(rp_max, -rz_val),
                     arrowprops=dict(arrowstyle="->", color="k", lw=1.0, linestyle="dashed"))

    cb = fig.colorbar(im, ax=ax, label=r"$\theta_B$ (degrees)", pad=0.02)
    cb.set_ticks([0, 15, 30, 45, 60, 75, 90])

    ax.set_xlabel(r"Coherence magnitude $r_\perp$")
    ax.set_ylabel(r"Population imbalance $r_z$")
    ax.set_xlim(0, 1.05); ax.set_ylim(-1.05, 1.05)
    ax.set_aspect("equal")
    ax.text(0.03, 0.88, r"$\theta_B = 90^\circ$", fontsize=8, color="k",
            fontweight="bold", alpha=0.8)
    ax.text(0.35, 0.06, r"$\theta_B = 90^\circ$", fontsize=8, color="k",
            fontweight="bold", alpha=0.8)

    fig.tight_layout()
    save(fig, "fig3_bures_angle")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 4 — FMO dimensional collapse trajectory (NEW)
# ══════════════════════════════════════════════════════════════════════════════

def fig4_fmo_collapse():
    """FMO transport efficiency and Bures principal angles vs dephasing.

    Two-panel figure:
      (a) Transport efficiency η(γ) with ENAQT peak (sink at site 3)
      (b) Minimum Bures principal angle θ_min(γ) showing geometry
    Both panels show site-1 and site-6 initial excitations.
    """
    print("Figure 4: FMO transport + geometry (dual panel) ...")

    H = H_FMO * CM_TO_RADFS
    d = D_FMO

    gamma_values_cm = np.concatenate([
        np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0]),
        np.arange(10, 110, 10),
        np.arange(150, 600, 50),
    ])
    gamma_values_cm = np.sort(np.unique(gamma_values_cm))

    # ── Panel (a): Transport efficiency ──
    print("  Computing transport efficiency (site 1)...")
    etas_s1 = scan_efficiency(gamma_values_cm, initial_site=0,
                               t_final_fs=15000.0, dt_fs=2.0)
    print("  Computing transport efficiency (site 6)...")
    etas_s6 = scan_efficiency(gamma_values_cm, initial_site=5,
                               t_final_fs=15000.0, dt_fs=2.0)

    # ── Panel (b): Principal angles (sink-inclusive, renormalized) ──
    kappa_trap_cm = 5.3  # 1 ps^-1
    kappa = kappa_trap_cm * CM_TO_RADFS

    def compute_angles_scan(initial_site):
        rho0 = np.zeros((d, d), dtype=complex)
        rho0[initial_site, initial_site] = 1.0
        min_angles = []
        for gamma_cm in gamma_values_cm:
            gamma = gamma_cm * CM_TO_RADFS
            # Evolve with same sink-inclusive dynamics as efficiency
            rho = evolve_lindblad_sink(H, gamma, kappa, rho0,
                                       t_final=5000.0, dt=1.0, sink_site=2)
            # Renormalize surviving state (conditional on not being trapped)
            tr = np.real(np.trace(rho))
            if tr > 1e-12:
                rho_renorm = rho / tr
            else:
                rho_renorm = np.eye(d, dtype=complex) / d
            # Regularize to full rank for QFIM
            rho_reg = rho_renorm + 1e-10 * np.eye(d) / d
            rho_reg /= np.trace(rho_reg)
            angles_deg = np.sort(np.degrees(principal_angles(rho_reg)))
            min_angles.append(angles_deg[0])
            print(f"    gamma={gamma_cm:7.1f}  min_angle={angles_deg[0]:5.1f}°"
                  f"  tr={tr:.4f}  (site {initial_site+1})")
        return np.array(min_angles)

    print("  Computing principal angles (site 1)...")
    angles_s1 = compute_angles_scan(0)
    print("  Computing principal angles (site 6)...")
    angles_s6 = compute_angles_scan(5)

    # ── Plot ──
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10.0, 4.0))

    # Panel (a): Transport efficiency
    ax_a.semilogx(gamma_values_cm, etas_s1, "o-", color=PAL[0], lw=1.5,
                   ms=3, markeredgecolor="white", markeredgewidth=0.4,
                   label="Site 1")
    ax_a.semilogx(gamma_values_cm, etas_s6, "s-", color=PAL[3], lw=1.5,
                   ms=3, markeredgecolor="white", markeredgewidth=0.4,
                   label="Site 6")

    # Mark peaks
    idx1 = np.argmax(etas_s1)
    idx6 = np.argmax(etas_s6)
    ax_a.axvline(gamma_values_cm[idx1], color=PAL[0], ls=":", lw=0.8, alpha=0.5)
    ax_a.axvline(gamma_values_cm[idx6], color=PAL[3], ls=":", lw=0.8, alpha=0.5)

    ax_a.set_xlabel(r"Dephasing rate $\gamma$ (cm$^{-1}$)")
    ax_a.set_ylabel(r"Transport efficiency $\eta$")
    ax_a.set_xlim(0.08, 600)
    ax_a.set_ylim(0, 1.0)
    ax_a.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax_a.set_title("(a)", loc="left", fontweight="bold", fontsize=11)

    # Panel (b): Min principal angle
    ax_b.semilogx(gamma_values_cm, angles_s1, "o-", color=PAL[0], lw=1.5,
                   ms=3, markeredgecolor="white", markeredgewidth=0.4,
                   label="Site 1")
    ax_b.semilogx(gamma_values_cm, angles_s6, "s-", color=PAL[3], lw=1.5,
                   ms=3, markeredgecolor="white", markeredgewidth=0.4,
                   label="Site 6")

    ax_b.axhline(90, color=PAL[5], ls="--", lw=0.8, alpha=0.5)
    ax_b.text(0.12, 91, r"$90°$", fontsize=8, color=PAL[5], alpha=0.7)

    # Mark the efficiency-peak γ values on the angle plot
    ax_b.axvline(gamma_values_cm[idx1], color=PAL[0], ls=":", lw=0.8, alpha=0.5)
    ax_b.axvline(gamma_values_cm[idx6], color=PAL[3], ls=":", lw=0.8, alpha=0.5)

    ax_b.set_xlabel(r"Dephasing rate $\gamma$ (cm$^{-1}$)")
    ax_b.set_ylabel(r"Min principal angle $\theta_{\min}$ (degrees)")
    ax_b.set_xlim(0.08, 600)
    ax_b.set_ylim(0, 95)
    ax_b.legend(loc="center right", fontsize=9, framealpha=0.9)
    ax_b.set_title("(b)", loc="left", fontweight="bold", fontsize=11)

    fig.tight_layout()
    save(fig, "fig4_fmo_collapse")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 5 — Qudit principal angles: numerical study (including d=7 for FMO)
# ══════════════════════════════════════════════════════════════════════════════

def _random_density_matrix(d, rng):
    A = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    rho = A @ A.conj().T
    rho /= np.trace(rho)
    return rho


def _sld_matrix(rho, X):
    eigvals, U = np.linalg.eigh(rho)
    d = len(eigvals)
    X_eig = U.conj().T @ X @ U
    L_eig = np.zeros((d, d), dtype=complex)
    for m in range(d):
        for n in range(d):
            denom = eigvals[m] + eigvals[n]
            if denom > 1e-14:
                L_eig[m, n] = 2.0 * X_eig[m, n] / denom
    return U @ L_eig @ U.conj().T


def _bures_inner(rho, X, Y):
    L_X = _sld_matrix(rho, X)
    L_Y = _sld_matrix(rho, Y)
    return 0.25 * np.real(np.trace(rho @ (L_X @ L_Y + L_Y @ L_X))) / 2.0


def _principal_angles(rho, d):
    diag_basis = []
    for a in range(d - 1):
        D = np.zeros((d, d), dtype=complex)
        D[a, a] = 1.0
        D[a + 1, a + 1] = -1.0
        diag_basis.append(D)

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

    n_diag = len(diag_basis)
    n_off = len(off_basis)

    G_DD = np.zeros((n_diag, n_diag))
    for i in range(n_diag):
        for j in range(i, n_diag):
            v = _bures_inner(rho, diag_basis[i], diag_basis[j])
            G_DD[i, j] = v
            G_DD[j, i] = v

    G_OO = np.zeros((n_off, n_off))
    for i in range(n_off):
        for j in range(i, n_off):
            v = _bures_inner(rho, off_basis[i], off_basis[j])
            G_OO[i, j] = v
            G_OO[j, i] = v

    G_DO = np.zeros((n_diag, n_off))
    for i in range(n_diag):
        for j in range(n_off):
            G_DO[i, j] = _bures_inner(rho, diag_basis[i], off_basis[j])

    def matrix_sqrt_inv(M):
        eigvals, V = np.linalg.eigh(M)
        eigvals = np.maximum(eigvals, 1e-12)
        return V @ np.diag(1.0 / np.sqrt(eigvals)) @ V.T

    M = matrix_sqrt_inv(G_DD) @ G_DO @ matrix_sqrt_inv(G_OO)
    svs = np.linalg.svd(M, compute_uv=False)
    cosines = np.clip(svs[:n_diag], 0.0, 1.0)
    return np.arccos(cosines)


def _params_to_rho(params, d):
    log_w = params[:d]
    diag = np.exp(log_w - np.max(log_w))
    diag /= diag.sum()

    n_off = d * (d - 1) // 2
    re_off = params[d : d + n_off]
    im_off = params[d + n_off : d + 2 * n_off]

    rho = np.diag(diag.astype(complex))
    idx = 0
    for k in range(d):
        for l in range(k + 1, d):
            max_mag = 0.95 * np.sqrt(diag[k] * diag[l])
            min_mag = 0.05 * np.sqrt(diag[k] * diag[l])
            raw = np.sqrt(re_off[idx]**2 + im_off[idx]**2)
            frac = min_mag + (max_mag - min_mag) * (np.tanh(raw) if raw > 0 else 0.0)
            if raw > 1e-15:
                phase = complex(re_off[idx], im_off[idx]) / raw
            else:
                phase = complex(1.0, 0.0)
            z = frac * phase
            rho[k, l] = z
            rho[l, k] = np.conj(z)
            idx += 1

    return rho


def _adversarial_search(d, n_starts, rng):
    n_off = d * (d - 1) // 2
    n_params = d + 2 * n_off

    best_angle = 0.0
    best_rho = None

    for _ in range(n_starts):
        x0 = np.zeros(n_params)
        x0[:d] = rng.standard_normal(d)
        x0[d:] = rng.standard_normal(2 * n_off) * 0.5

        def neg_min_angle(params):
            rho = _params_to_rho(params, d)
            eigvals = np.linalg.eigvalsh(rho)
            if eigvals.min() < 1e-10:
                return 0.0
            off_diag_norm = np.sqrt(sum(
                abs(rho[k, l])**2 for k in range(d) for l in range(k+1, d)
            ))
            if off_diag_norm < 1e-6:
                return 0.0
            angles = _principal_angles(rho, d)
            return -np.min(angles)

        result = minimize(neg_min_angle, x0, method="L-BFGS-B",
                          options={"maxiter": 200, "ftol": 1e-10})

        achieved = -result.fun
        if achieved > best_angle:
            best_angle = achieved
            best_rho = _params_to_rho(result.x, d)

    return np.degrees(best_angle), best_rho


def fig5_qudit_angles():
    print("Figure 5: Qudit principal angles (including d=7 for FMO) ...")

    rng = np.random.default_rng(42)
    dims = [3, 5, 7]  # Include d=7 for FMO
    n_samples_map = {3: 2000, 5: 1000, 7: 500}  # Fewer samples for d=7 (expensive)
    n_optim_map = {3: 100, 5: 50, 7: 20}

    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.8), sharey=True)

    for ax, d in zip(axes, dims):
        n_samples = n_samples_map[d]
        n_optim = n_optim_map[d]

        min_angles = []
        for i in range(n_samples):
            rho = _random_density_matrix(d, rng)
            angles = _principal_angles(rho, d)
            min_angles.append(np.degrees(np.min(angles)))
            if (i + 1) % 200 == 0:
                print(f"  d={d}: {i+1}/{n_samples} random states done")

        print(f"  d={d}: adversarial search ({n_optim} starts) ...", end="",
              flush=True)
        opt_angle, opt_rho = _adversarial_search(d, n_optim, rng)
        print(f" best = {opt_angle:.2f}")

        ax.hist(min_angles, bins=50, color=PAL[0], alpha=0.8,
                edgecolor="white", linewidth=0.3, density=True)
        ax.axvline(90, color=PAL[3], ls="--", lw=1.2, alpha=0.8)
        ax.axvline(opt_angle, color=PAL[1], ls="-", lw=1.5, alpha=0.9)
        ax.set_xlabel(r"$\min_a \theta_a$ (degrees)")
        label = rf"$d = {d}$"
        if d == 7:
            label += " (FMO)"
        ax.set_title(label, fontsize=11)
        ax.set_xlim(0, 95)

        median_val = np.median(min_angles)
        max_val = np.max(min_angles)
        ax.text(0.97, 0.95,
                f"median = {median_val:.1f}\u00b0\n"
                f"max = {max_val:.1f}\u00b0\n"
                f"optim = {opt_angle:.1f}\u00b0",
                transform=ax.transAxes, fontsize=7, ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor=PAL[5], alpha=0.8))

    axes[0].set_ylabel("Density")
    fig.tight_layout(w_pad=1.0)
    save(fig, "fig5_qudit_angles")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 6 — Angle spectrum + sink-rate robustness
# ══════════════════════════════════════════════════════════════════════════════

def fig6_robustness():
    """Full principal angle spectrum and robustness to trapping rate.

    Two-panel figure:
      (a) All 6 principal angles vs γ (site-1, Γ_trap = 1 ps⁻¹)
      (b) Min angle at ENAQT optimum vs Γ_trap (both sites)
    """
    print("Figure 6: Angle spectrum + sink-rate robustness ...")

    H = H_FMO * CM_TO_RADFS
    d = D_FMO

    # ── Panel (a): Full 6-angle spectrum for site-1 ──
    kappa_ref_cm = 5.3  # 1 ps⁻¹
    kappa_ref = kappa_ref_cm * CM_TO_RADFS

    gamma_values_cm = np.concatenate([
        np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0]),
        np.arange(10, 110, 10),
        np.arange(150, 600, 50),
    ])
    gamma_values_cm = np.sort(np.unique(gamma_values_cm))

    rho0_s1 = np.zeros((d, d), dtype=complex)
    rho0_s1[0, 0] = 1.0

    all_angles = []  # (n_gamma, 6)
    l1_vals = []
    for gamma_cm in gamma_values_cm:
        gamma = gamma_cm * CM_TO_RADFS
        rho = evolve_lindblad_sink(H, gamma, kappa_ref, rho0_s1,
                                    t_final=5000.0, dt=1.0, sink_site=2)
        tr = np.real(np.trace(rho))
        rho_n = rho / tr if tr > 1e-12 else np.eye(d, dtype=complex) / d
        rho_reg = rho_n + 1e-10 * np.eye(d) / d
        rho_reg /= np.trace(rho_reg)
        angles_deg = np.sort(np.degrees(principal_angles(rho_reg)))
        all_angles.append(angles_deg)
        l1_vals.append(l1_coherence(rho_reg))
        print(f"    gamma={gamma_cm:7.1f}  angles: {angles_deg[0]:.1f}..{angles_deg[-1]:.1f}")

    all_angles = np.array(all_angles)  # (n_gamma, 6)
    l1_vals = np.array(l1_vals)

    # Print scalar comparison values for tex remark
    for gc in [1.0, 10.0, 100.0, 150.0]:
        idx = np.argmin(np.abs(gamma_values_cm - gc))
        print(f"  ** gamma={gc}: l1={l1_vals[idx]:.4f}, "
              f"angles={np.round(all_angles[idx], 1)}")

    # ── Panel (b): Sink-rate robustness ──
    trap_rates_per_ps = np.array([0.2, 0.5, 1.0, 2.0, 3.0, 5.0])

    gamma_coarse = np.concatenate([
        np.array([0.1, 1.0, 5.0]),
        np.arange(10, 110, 10),
        np.arange(150, 600, 50),
    ])
    gamma_coarse = np.sort(np.unique(gamma_coarse))

    rho0_s6 = np.zeros((d, d), dtype=complex)
    rho0_s6[5, 5] = 1.0

    results_s1 = []
    results_s6 = []

    for Gamma_ps in trap_rates_per_ps:
        kappa_cm = Gamma_ps * 5.3
        kappa_val = kappa_cm * CM_TO_RADFS

        # Find ENAQT optimum for site-1
        etas_s1 = [compute_transport_efficiency(g, kappa_trap_cm=kappa_cm,
                    t_final_fs=15000.0, dt_fs=2.0, initial_site=0)
                   for g in gamma_coarse]
        idx_opt_s1 = np.argmax(etas_s1)
        g_opt_s1 = gamma_coarse[idx_opt_s1]

        # Compute angle at optimum for site-1
        rho = evolve_lindblad_sink(H, g_opt_s1 * CM_TO_RADFS, kappa_val, rho0_s1,
                                    t_final=5000.0, dt=1.0, sink_site=2)
        tr = np.real(np.trace(rho))
        rho_n = rho / tr if tr > 1e-12 else np.eye(d, dtype=complex) / d
        rho_reg = rho_n + 1e-10 * np.eye(d) / d
        rho_reg /= np.trace(rho_reg)
        min_ang_s1 = np.degrees(np.min(principal_angles(rho_reg)))
        results_s1.append((Gamma_ps, g_opt_s1, min_ang_s1))

        # Find ENAQT optimum for site-6
        etas_s6 = [compute_transport_efficiency(g, kappa_trap_cm=kappa_cm,
                    t_final_fs=15000.0, dt_fs=2.0, initial_site=5)
                   for g in gamma_coarse]
        idx_opt_s6 = np.argmax(etas_s6)
        g_opt_s6 = gamma_coarse[idx_opt_s6]

        # Compute angle at optimum for site-6
        rho = evolve_lindblad_sink(H, g_opt_s6 * CM_TO_RADFS, kappa_val, rho0_s6,
                                    t_final=5000.0, dt=1.0, sink_site=2)
        tr = np.real(np.trace(rho))
        rho_n = rho / tr if tr > 1e-12 else np.eye(d, dtype=complex) / d
        rho_reg = rho_n + 1e-10 * np.eye(d) / d
        rho_reg /= np.trace(rho_reg)
        min_ang_s6 = np.degrees(np.min(principal_angles(rho_reg)))
        results_s6.append((Gamma_ps, g_opt_s6, min_ang_s6))

        print(f"    Gamma={Gamma_ps:.1f} ps^-1: "
              f"s1 gamma_opt={g_opt_s1:.0f} angle={min_ang_s1:.1f}, "
              f"s6 gamma_opt={g_opt_s6:.0f} angle={min_ang_s6:.1f}")

    results_s1 = np.array(results_s1)
    results_s6 = np.array(results_s6)

    # ── Plot ──
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10.0, 4.0))

    # Panel (a): Full angle spectrum
    angle_colors = plt.cm.viridis(np.linspace(0.15, 0.85, 6))
    for i in range(6):
        label = rf"$\theta_{i+1}$" if i in [0, 5] else None
        ax_a.semilogx(gamma_values_cm, all_angles[:, i], "-",
                       color=angle_colors[i], lw=1.2, label=label)
    # Shaded band between min and max
    ax_a.fill_between(gamma_values_cm, all_angles[:, 0], all_angles[:, 5],
                       alpha=0.12, color=PAL[0])
    ax_a.axhline(90, color=PAL[5], ls="--", lw=0.8, alpha=0.5)
    ax_a.set_xlabel(r"Dephasing rate $\gamma$ (cm$^{-1}$)")
    ax_a.set_ylabel(r"Principal angles (degrees)")
    ax_a.set_xlim(0.08, 600)
    ax_a.set_ylim(0, 95)
    ax_a.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax_a.set_title("(a)", loc="left", fontweight="bold", fontsize=11)

    # Panel (b): Min angle vs trap rate
    ax_b.plot(results_s1[:, 0], results_s1[:, 2], "o-", color=PAL[0], lw=1.5,
              ms=6, markeredgecolor="white", markeredgewidth=0.8, label="Site 1")
    ax_b.plot(results_s6[:, 0], results_s6[:, 2], "s-", color=PAL[3], lw=1.5,
              ms=6, markeredgecolor="white", markeredgewidth=0.8, label="Site 6")
    ax_b.axvline(1.0, color=PAL[5], ls=":", lw=0.8, alpha=0.5)
    ax_b.text(1.05, 84.5, r"$\Gamma_{\rm trap} = 1\,{\rm ps}^{-1}$",
              fontsize=8, color=PAL[5])
    ax_b.axhline(85, color=PAL[1], ls="--", lw=0.6, alpha=0.4)
    ax_b.set_xlabel(r"Trapping rate $\Gamma_{\rm trap}$ (ps$^{-1}$)")
    ax_b.set_ylabel(r"Min angle at $\gamma_{\rm opt}$ (degrees)")
    ax_b.set_xlim(0, 5.5)
    ax_b.set_ylim(83, 91)
    ax_b.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax_b.set_title("(b)", loc="left", fontweight="bold", fontsize=11)

    fig.tight_layout()
    save(fig, "fig6_robustness")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 7 — PE545 transport + geometry (dual panel)
# ══════════════════════════════════════════════════════════════════════════════

def fig7_pe545():
    """PE545 transport efficiency and Bures principal angles vs dephasing.

    Two-panel figure matching fig4 style:
      (a) Transport efficiency η(γ) for both initial sites (PEB50/61C, PEB50/61D)
      (b) Seven Bures principal angles vs γ, with throughput-optimal dephasing marked
    """
    print("Figure 7: PE545 transport + geometry (dual panel) ...")

    d = 8
    target_sites = [1, 2]     # DBV19A, DBV19B (red pool)
    kappa_ps = [0.5, 0.5]     # 1.0 ps^-1 total

    gamma_values_cm = np.concatenate([
        np.array([1.0, 2.0, 5.0]),
        np.arange(10, 110, 10),
        np.arange(150, 1100, 50),
        np.arange(1200, 1600, 100),
    ])
    gamma_values_cm = np.sort(np.unique(gamma_values_cm))

    # ── Panel (a): Transport efficiency ──
    print("  Computing transport efficiency (PEB50/61C)...")
    etas_s0 = np.array([
        compute_transport_generic(H_PE545, g, target_sites, kappa_ps, 0, t_final_ps=30.0, dt_fs=2.0)
        for g in gamma_values_cm
    ])
    print("  Computing transport efficiency (PEB50/61D)...")
    etas_s5 = np.array([
        compute_transport_generic(H_PE545, g, target_sites, kappa_ps, 5, t_final_ps=30.0, dt_fs=2.0)
        for g in gamma_values_cm
    ])

    idx0 = np.argmax(etas_s0)
    idx5 = np.argmax(etas_s5)
    print(f"    PEB50/61C peak: eta={etas_s0[idx0]:.4f} at gamma={gamma_values_cm[idx0]:.0f}")
    print(f"    PEB50/61D peak: eta={etas_s5[idx5]:.4f} at gamma={gamma_values_cm[idx5]:.0f}")

    # ── Panel (b): All 7 principal angles for BOTH initial sites ──
    angle_gammas = np.concatenate([
        np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]),
        np.arange(200, 1100, 100),
        np.array([1500.0]),
    ])

    def _compute_pe545_angles(initial_site, targets, kappas, label):
        print(f"  Computing principal angles ({label})...")
        all_angles = []
        rho0 = np.zeros((d, d), dtype=complex)
        rho0[initial_site, initial_site] = 1.0
        for gamma_cm in angle_gammas:
            rho = evolve_lindblad_sink_generic(H_PE545, gamma_cm, targets, kappas,
                                                rho0, t_final_ps=5.0, dt_fs=1.0)
            tr = np.real(np.trace(rho))
            rho_n = rho / tr if tr > 1e-12 else np.eye(d, dtype=complex) / d
            rho_reg = rho_n + 1e-10 * np.eye(d) / d
            rho_reg /= np.trace(rho_reg)
            angles_deg = np.sort(np.degrees(principal_angles(rho_reg)))
            all_angles.append(angles_deg)
            print(f"    gamma={gamma_cm:7.0f}  angles: {angles_deg[0]:.1f}..{angles_deg[-1]:.1f}")
        return np.array(all_angles)

    all_angles_s0 = _compute_pe545_angles(0, target_sites, kappa_ps, "PEB50/61C")
    all_angles_s5 = _compute_pe545_angles(5, target_sites, kappa_ps, "PEB50/61D")

    # ── Plot ──
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10.0, 4.0))

    # Panel (a): Transport efficiency
    ax_a.semilogx(gamma_values_cm, etas_s0, "o-", color=PAL[0], lw=1.5,
                   ms=3, markeredgecolor="white", markeredgewidth=0.4,
                   label="PEB50/61C")
    ax_a.semilogx(gamma_values_cm, etas_s5, "s-", color=PAL[3], lw=1.5,
                   ms=3, markeredgecolor="white", markeredgewidth=0.4,
                   label="PEB50/61D")

    # Mark peaks
    ax_a.axvline(gamma_values_cm[idx0], color=PAL[0], ls=":", lw=0.8, alpha=0.5)
    ax_a.axvline(gamma_values_cm[idx5], color=PAL[3], ls=":", lw=0.8, alpha=0.5)

    ax_a.set_xlabel(r"Dephasing rate $\gamma$ (cm$^{-1}$)")
    ax_a.set_ylabel(r"Transport efficiency $\eta$")
    ax_a.set_xlim(0.8, 1600)
    ax_a.set_ylim(0, 1.0)
    ax_a.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax_a.set_title("(a)", loc="left", fontweight="bold", fontsize=11)

    # Panel (b): Full 7-angle spectrum (PEB50/61C) + min angle for PEB50/61D
    angle_colors = plt.cm.viridis(np.linspace(0.15, 0.85, 7))
    for i in range(7):
        label = rf"$\theta_{i+1}$" if i in [0, 6] else None
        ax_b.semilogx(angle_gammas, all_angles_s0[:, i], "-",
                       color=angle_colors[i], lw=1.2, label=label)

    ax_b.fill_between(angle_gammas, all_angles_s0[:, 0], all_angles_s0[:, 6],
                       alpha=0.12, color=PAL[0])

    # Overlay min angle for site 5 (dashed) to show both-site robustness
    ax_b.semilogx(angle_gammas, all_angles_s5[:, 0], "--",
                   color=PAL[3], lw=1.5, label=r"$\theta_{\min}$ (PEB50/61D)")

    ax_b.axhline(90, color=PAL[5], ls="--", lw=0.8, alpha=0.5)

    # Mark throughput-optimal γ
    ax_b.axvline(gamma_values_cm[idx0], color=PAL[0], ls=":", lw=0.8, alpha=0.5)

    ax_b.set_xlabel(r"Dephasing rate $\gamma$ (cm$^{-1}$)")
    ax_b.set_ylabel(r"Principal angles (degrees)")
    ax_b.set_xlim(0.8, 1600)
    ax_b.set_ylim(0, 95)
    ax_b.legend(loc="lower right", fontsize=8, framealpha=0.9)
    ax_b.set_title("(b)", loc="left", fontweight="bold", fontsize=11)

    fig.tight_layout()
    save(fig, "fig7_pe545")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 8 — Cross-scale quantum-classical boundary (neural oscillations)
# ══════════════════════════════════════════════════════════════════════════════

def fig8_neural_comparison():
    """Cross-scale quantum-classical boundary: photosynthesis vs neural oscillations.

    Two-panel figure:
      (a) θ_min vs dimensionless dephasing γ/J_max for FMO, PE545, neural (d=10)
          on log scale, with biological dephasing marked for each system
      (b) Detection threshold: θ_min vs hypothetical dephasing for neural system,
          with microenvironment shielding factor on secondary x-axis
    """
    print("Figure 8: Cross-scale quantum-classical boundary ...")

    # ── Compute θ_min vs γ/J_max for all three systems ──

    # FMO: J_max ≈ 87.7 cm^-1 (sites 1-2 coupling)
    fmo_jmax = np.max(np.abs(H_FMO - np.diag(np.diag(H_FMO))))
    # PE545: J_max ≈ 92 cm^-1
    pe545_jmax = np.max(np.abs(H_PE545 - np.diag(np.diag(H_PE545))))
    # Chromatin: J_max = 25 cm^-1
    chromatin_jmax = CHROMATIN_JMAX_CM
    # Neural: J_max = 250 cm^-1 (scaled)
    neural_jmax = NEURAL_GAMMA_JMAX_CM

    # Dimensionless γ/J_max scan: 0.01 to 10^6. Beyond γ/J_max ~ 100 the
    # system is fully decohered (θ_min = 90° to machine precision), so we
    # cap the Lindblad simulation range and set 90° analytically above that.
    gamma_over_j_values = np.logspace(-2, 6, 60)
    GAMMA_SIM_CAP = 500.0  # max γ/J_max for actual Lindblad simulation

    def compute_theta_min_scan(H_cm, initial_site, target_sites, kappa_ps,
                                jmax_cm, gamma_over_j_arr, label):
        """Compute minimum principal angle at each γ/J_max value.

        For γ/J_max > GAMMA_SIM_CAP, returns 90° as an asymptotic limit:
        at high dephasing the state approaches the classical simplex and
        all principal angles approach 90°.  This is a physical extrapolation,
        not a simulation result.  If simulation produces NaN/inf (numerical
        overflow at high γ), also returns 90° on the same asymptotic grounds.
        """
        d = H_cm.shape[0]
        min_angles = []
        rho0 = np.zeros((d, d), dtype=complex)
        rho0[initial_site, initial_site] = 1.0

        for goj in gamma_over_j_arr:
            if goj > GAMMA_SIM_CAP:
                min_angles.append(90.0)
                continue
            gamma_cm = goj * jmax_cm
            rho = evolve_lindblad_sink_generic(
                H_cm, gamma_cm, target_sites, kappa_ps,
                rho0, t_final_ps=5.0, dt_fs=1.0)
            # Check for numerical overflow
            if not np.all(np.isfinite(rho)):
                print(f"    WARNING: NaN/inf at γ/J_max={goj:.2f} for {label}"
                      " — returning 90° (fully decohered)")
                min_angles.append(90.0)
                continue
            tr = np.real(np.trace(rho))
            if tr < 1e-12:
                # Population fully trapped — state on simplex boundary
                min_angles.append(90.0)
                continue
            rho_n = rho / tr
            rho_reg = rho_n + 1e-10 * np.eye(d) / d
            rho_reg /= np.trace(rho_reg)
            angles_deg = np.sort(np.degrees(principal_angles(rho_reg)))
            min_angles.append(angles_deg[0])

        print(f"    {label}: {len(gamma_over_j_arr)} points computed")
        return np.array(min_angles)

    print("  Computing FMO angles...")
    fmo_angles = compute_theta_min_scan(
        H_FMO, 0, [2], [1.0], fmo_jmax, gamma_over_j_values, "FMO")

    print("  Computing PE545 angles...")
    pe545_angles = compute_theta_min_scan(
        H_PE545, 0, [1, 2], [0.5, 0.5], pe545_jmax, gamma_over_j_values, "PE545")

    print("  Computing Chromatin (d=8) angles...")
    chromatin_angles = compute_theta_min_scan(
        H_CHROMATIN, 0, [3, 4], [0.5, 0.5], chromatin_jmax,
        gamma_over_j_values, "Chromatin d=8")

    print("  Computing Neural (d=10) angles...")
    neural_angles = compute_theta_min_scan(
        H_NEURAL_GAMMA, 0, [8, 9], [0.5, 0.5], neural_jmax,
        gamma_over_j_values, "Neural d=10")

    # Biological dephasing ratios (dimensionless γ_bio/J_max)
    # FMO: γ_bio ≈ 100 cm^-1 at 300K → γ_bio/J_max ≈ 1.1
    fmo_bio_ratio = 100.0 / fmo_jmax
    # PE545: γ_bio ≈ 725 cm^-1 → γ_bio/J_max ≈ 7.9
    pe545_bio_ratio = 725.0 / pe545_jmax
    # Chromatin: γ_bio = 200 cm^-1 → γ_bio/J_max = 8.0
    chromatin_bio_ratio = CHROMATIN_BIO_GAMMA_OVER_JMAX
    # Neural: γ_bio/J_max ≈ 8.12e12
    neural_bio_ratio = NEURAL_BIO_GAMMA_OVER_JMAX

    # ── Panel (b): Detection threshold for neural system ──
    # Scan γ/J_max from 0.01 to 10^4 (simulation range), then extend
    # analytically to biological dephasing (90° plateau).
    detection_goj_sim = np.logspace(-2, np.log10(GAMMA_SIM_CAP), 50)
    # Add a few points in the plateau to show the biological dephasing mark
    detection_goj_plateau = np.logspace(np.log10(GAMMA_SIM_CAP) + 0.5,
                                         np.log10(neural_bio_ratio) + 0.5, 10)
    detection_goj = np.concatenate([detection_goj_sim, detection_goj_plateau])

    print("  Computing neural detection threshold...")
    detection_angles = compute_theta_min_scan(
        H_NEURAL_GAMMA, 0, [8, 9], [0.5, 0.5], neural_jmax,
        detection_goj, "Neural detection")

    # Find the ENAQT-equivalent regime (where θ_min first exceeds 85°)
    enaqt_idx = np.argmax(detection_angles > 85.0)
    if enaqt_idx > 0:
        enaqt_goj = detection_goj[enaqt_idx]
        shielding_needed = neural_bio_ratio / enaqt_goj
        print(f"    ENAQT-equivalent at γ/J_max ≈ {enaqt_goj:.2f}")
        print(f"    Shielding factor needed: {shielding_needed:.2e}")

    # ── Plot ──
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10.0, 4.2))

    # Panel (a): θ_min vs γ/J_max for all three systems
    ax_a.semilogx(gamma_over_j_values, fmo_angles, "o-", color=PAL[0], lw=1.5,
                   ms=3, markeredgecolor="white", markeredgewidth=0.4,
                   label=r"FMO ($d = 7$)")
    ax_a.semilogx(gamma_over_j_values, pe545_angles, "s-", color=PAL[1], lw=1.5,
                   ms=3, markeredgecolor="white", markeredgewidth=0.4,
                   label=r"PE545 ($d = 8$)")
    ax_a.semilogx(gamma_over_j_values, chromatin_angles, "D-", color=PAL[4], lw=1.5,
                   ms=3, markeredgecolor="white", markeredgewidth=0.4,
                   label=r"Chromatin ($d = 8$)")
    ax_a.semilogx(gamma_over_j_values, neural_angles, "^-", color=PAL[2], lw=1.5,
                   ms=3, markeredgecolor="white", markeredgewidth=0.4,
                   label=r"Neural ($d = 10$)")

    # Biological dephasing markers
    ax_a.axvline(fmo_bio_ratio, color=PAL[0], ls="--", lw=0.8, alpha=0.6)
    ax_a.axvline(pe545_bio_ratio, color=PAL[1], ls="--", lw=0.8, alpha=0.6)
    ax_a.axvline(chromatin_bio_ratio, color=PAL[4], ls="--", lw=0.8, alpha=0.6)

    # Neural bio dephasing is far off-scale; annotate with arrow
    ax_a.annotate(r"$\gamma_{\rm bio}/J_{\rm max} \sim 10^{13}$" "\n(neural)",
                  xy=(gamma_over_j_values[-1], 90), xytext=(1e4, 75),
                  fontsize=8, color=PAL[2], ha="center",
                  arrowprops=dict(arrowstyle="->", color=PAL[2], lw=1.0))

    ax_a.axhline(90, color=PAL[5], ls="--", lw=0.8, alpha=0.5)
    ax_a.text(0.015, 91.2, r"$90°$ (classical)", fontsize=8, color=PAL[5],
              alpha=0.7)

    ax_a.set_xlabel(r"Dimensionless dephasing $\gamma / J_{\max}$")
    ax_a.set_ylabel(r"Min principal angle $\theta_{\min}$ (degrees)")
    ax_a.set_xlim(0.008, 2e6)
    ax_a.set_ylim(0, 95)
    ax_a.legend(loc="center right", fontsize=9, framealpha=0.9)
    ax_a.set_title("(a)", loc="left", fontweight="bold", fontsize=11)

    # Panel (b): Detection threshold
    ax_b.semilogx(detection_goj, detection_angles, "-", color=PAL[2], lw=2.0)
    ax_b.axvline(neural_bio_ratio, color=PAL[3], ls="--", lw=1.2, alpha=0.8)
    ax_b.text(neural_bio_ratio * 1.5, 50,
              r"$\gamma_{\rm bio}$" "\n" r"($kT/\hbar$, 310 K)",
              fontsize=8, color=PAL[3], va="center")

    # Mark ENAQT-equivalent regime
    ax_b.axhline(85, color=PAL[5], ls=":", lw=0.8, alpha=0.5)
    ax_b.text(0.015, 86.5, r"$85°$", fontsize=8, color=PAL[5])

    ax_b.axhline(90, color=PAL[5], ls="--", lw=0.8, alpha=0.5)
    ax_b.text(0.015, 91.2, r"$90°$", fontsize=8, color=PAL[5])

    # Shielding factor axis (top)
    ax_b2 = ax_b.twiny()
    ax_b2.set_xscale("log")
    # shielding factor = γ_bio / γ_effective → as γ decreases, shielding increases
    # x-axis is γ/J_max; shielding = γ_bio_ratio / (γ/J_max)
    shielding_ticks = [1e0, 1e3, 1e6, 1e9, 1e12]
    shielding_positions = [neural_bio_ratio / s for s in shielding_ticks]
    ax_b2.set_xlim(ax_b.get_xlim())
    ax_b2.set_xticks(shielding_positions)
    ax_b2.set_xticklabels([r"$10^{0}$", r"$10^{3}$", r"$10^{6}$",
                            r"$10^{9}$", r"$10^{12}$"], fontsize=8)
    ax_b2.set_xlabel("Microenvironment shielding factor", fontsize=10)

    ax_b.set_xlabel(r"Dimensionless dephasing $\gamma / J_{\max}$")
    ax_b.set_ylabel(r"Min principal angle $\theta_{\min}$ (degrees)")
    ax_b.set_xlim(0.008, neural_bio_ratio * 5)
    ax_b.set_ylim(0, 95)
    ax_b.set_title("(b)", loc="left", fontweight="bold", fontsize=11)

    fig.tight_layout()
    save(fig, "fig8_neural_comparison")


def neural_robustness_check():
    """Run both neural models and print comparison table for Appendix C.

    Verifies:
      - Both give θ_min ≈ 90° at biological dephasing
      - ENAQT-equivalent regime at γ/J_max ~ O(1) for both
    Uses per-system parameters matching main figure computations.
    """
    print("Neural robustness check (Appendix C) ...")

    # Per-system parameters: (label, H, init_site, targets, kappas, jmax_cm,
    #                         t_final_ps for transport scan, bio_gamma_cm)
    all_models = [
        ("FMO (d=7)", H_FMO, 0, [2], [1.0],
         np.max(np.abs(H_FMO - np.diag(np.diag(H_FMO)))),
         15.0, 100.0),  # FMO: 15 ps, bio γ ≈ 100 cm⁻¹
        ("PE545 (d=8)", H_PE545, 0, [1, 2], [0.5, 0.5],
         np.max(np.abs(H_PE545 - np.diag(np.diag(H_PE545)))),
         30.0, 725.0),  # PE545: 30 ps, bio γ ≈ 725 cm⁻¹
        ("Chromatin (d=8)", H_CHROMATIN, 0, [3, 4], [0.5, 0.5],
         CHROMATIN_JMAX_CM,
         30.0, CHROMATIN_BIO_GAMMA_CM),  # Chromatin: 30 ps, bio γ = 200 cm⁻¹
        ("NeuralGamma (d=10)", H_NEURAL_GAMMA, 0, [8, 9], [0.5, 0.5],
         NEURAL_GAMMA_JMAX_CM,
         15.0, None),  # Neural: bio γ too large, use proxy
        ("NeuralPD (d=8)", H_NEURAL_PD, 2, [4, 6], [0.5, 0.5],
         NEURAL_PD_JMAX_CM,
         15.0, None),  # Neural: bio γ too large, use proxy
    ]

    print(f"\n  {'System':<22s}  {'d':>2s}  {'J_max':>8s}  "
          f"{'γ_opt/J':>8s}  {'θ_min@opt':>10s}  {'θ_min@bio':>10s}")
    print("  " + "-" * 72)

    for label, H_cm, init, targets, kappas, jmax, t_ps, bio_gamma in all_models:
        d = H_cm.shape[0]
        rho0 = np.zeros((d, d), dtype=complex)
        rho0[init, init] = 1.0

        # ENAQT scan: γ/J_max from 0.01 to 20 (well above known optima,
        # avoids edge-of-range artifacts)
        goj_scan = np.logspace(-2, np.log10(20), 60)
        etas = []
        for goj in goj_scan:
            gamma_cm = goj * jmax
            eta = compute_transport_generic(H_cm, gamma_cm, targets, kappas,
                                             init, t_final_ps=t_ps, dt_fs=2.0)
            etas.append(eta)
        etas = np.array(etas)
        idx_opt = np.argmax(etas)
        goj_opt = goj_scan[idx_opt]

        # Angle at optimum
        rho = evolve_lindblad_sink_generic(
            H_cm, goj_opt * jmax, targets, kappas,
            rho0, t_final_ps=5.0, dt_fs=1.0)
        tr = np.real(np.trace(rho))
        if not np.all(np.isfinite(rho)) or tr < 1e-12:
            angles_opt = np.array([90.0])
        else:
            rho_n = rho / tr
            rho_reg = rho_n + 1e-10 * np.eye(d) / d
            rho_reg /= np.trace(rho_reg)
            angles_opt = np.sort(np.degrees(principal_angles(rho_reg)))

        # Angle at biological dephasing
        if bio_gamma is None:
            # Neural: biological γ_bio/J_max ~ 10^12, far beyond sim range.
            # Use proxy γ/J_max = 50 to verify angles are near 90°.
            # Reported as proxy-sim value, NOT claimed as "90.00°".
            gamma_bio_cm = 50.0 * jmax
        else:
            gamma_bio_cm = bio_gamma

        rho0_bio = np.zeros((d, d), dtype=complex)
        rho0_bio[init, init] = 1.0
        rho_bio = evolve_lindblad_sink_generic(
            H_cm, gamma_bio_cm, targets, kappas,
            rho0_bio, t_final_ps=5.0, dt_fs=1.0)
        tr_bio = np.real(np.trace(rho_bio))
        if not np.all(np.isfinite(rho_bio)) or tr_bio < 1e-12:
            angles_bio = np.array([90.0])
        else:
            rho_n_bio = rho_bio / tr_bio
            rho_reg_bio = rho_n_bio + 1e-10 * np.eye(d) / d
            rho_reg_bio /= np.trace(rho_reg_bio)
            angles_bio = np.sort(np.degrees(principal_angles(rho_reg_bio)))

        print(f"  {label:<22s}  {d:2d}  {jmax:8.1f}  "
              f"{goj_opt:8.2f}  {angles_opt[0]:10.2f}°  {angles_bio[0]:10.2f}°")

    print()


# ══════════════════════════════════════════════════════════════════════════════
# PE545 robustness check — backs specific claims in manuscript Section 5.7
# ══════════════════════════════════════════════════════════════════════════════

def pe545_robustness_check():
    """Verify PE545 robustness claims: alt target manifold and halved trap rate.

    Manuscript claims (Section 5.7):
      - Including PEB82C as 3rd target site shifts γ_opt by <50 cm⁻¹
      - Including PEB82C changes minimum angles by <0.3°
      - Halving trapping rate preserves the qualitative picture
    """
    print("PE545 robustness check (backing manuscript claims) ...")

    d = 8
    gamma_fine = np.concatenate([
        np.arange(400, 1100, 25),
    ])

    configs = [
        ("default: tgt={1,2}, κ=1.0",   [1, 2],    [0.5, 0.5]),
        ("alt tgt: tgt={1,2,3}, κ=1.0",  [1, 2, 3], [0.33, 0.33, 0.34]),
        ("half κ: tgt={1,2}, κ=0.5",     [1, 2],    [0.25, 0.25]),
    ]

    results = {}
    for label, targets, kappas in configs:
        print(f"\n  {label}")
        for site_idx, site_name in [(0, "PEB50/61C"), (5, "PEB50/61D")]:
            etas = np.array([
                compute_transport_generic(H_PE545, g, targets, kappas,
                                           site_idx, t_final_ps=30.0, dt_fs=2.0)
                for g in gamma_fine
            ])
            idx = np.argmax(etas)
            g_opt = gamma_fine[idx]

            # Angles at optimum
            rho0 = np.zeros((d, d), dtype=complex)
            rho0[site_idx, site_idx] = 1.0
            rho = evolve_lindblad_sink_generic(H_PE545, g_opt, targets, kappas,
                                                rho0, t_final_ps=5.0, dt_fs=1.0)
            tr = np.real(np.trace(rho))
            rho_n = rho / tr if tr > 1e-12 else np.eye(d, dtype=complex) / d
            rho_reg = rho_n + 1e-10 * np.eye(d) / d
            rho_reg /= np.trace(rho_reg)
            angles = np.sort(np.degrees(principal_angles(rho_reg)))

            key = (label, site_name)
            results[key] = (g_opt, etas[idx], angles[0])
            print(f"    {site_name}: γ_opt={g_opt:.0f} cm⁻¹, η={etas[idx]:.4f}, "
                  f"θ_min={angles[0]:.2f}°")

    # Check claims
    print("\n  === Claim verification ===")
    for site_name in ["PEB50/61C", "PEB50/61D"]:
        g_def = results[("default: tgt={1,2}, κ=1.0", site_name)][0]
        g_alt = results[("alt tgt: tgt={1,2,3}, κ=1.0", site_name)][0]
        a_def = results[("default: tgt={1,2}, κ=1.0", site_name)][2]
        a_alt = results[("alt tgt: tgt={1,2,3}, κ=1.0", site_name)][2]
        a_half = results[("half κ: tgt={1,2}, κ=0.5", site_name)][2]

        dg = abs(g_def - g_alt)
        da = abs(a_def - a_alt)
        print(f"    {site_name}: Δγ_opt={dg:.0f} cm⁻¹ (<50? {'YES' if dg < 50 else 'NO'}), "
              f"Δθ_min={da:.2f}° (<0.3? {'YES' if da < 0.3 else 'NO'}), "
              f"half-κ θ_min={a_half:.1f}° (>89? {'YES' if a_half > 89 else 'NO'})")


# ══════════════════════════════════════════════════════════════════════════════
# Chromatin sensitivity sweep (mandatory, main text)
# ══════════════════════════════════════════════════════════════════════════════

def chromatin_sensitivity_sweep():
    """Sweep γ_bio × J_max for chromatin model; report θ_min grid.

    Kill criterion: if θ_min < 85° for ANY combination, chromatin goes to
    appendix only.

    Sweeps over physically reasonable parameter ranges:
      γ_bio ∈ [100, 150, 200, 250, 300] cm⁻¹
      J_max ∈ [10, 15, 20, 25, 30, 40, 50] cm⁻¹

    For each (γ_bio, J_max): rebuilds Hamiltonian with scaled couplings,
    computes θ_min.
    """
    print("Chromatin sensitivity sweep ...")

    gamma_bio_values = [100, 150, 200, 250, 300]
    jmax_values = [10, 15, 20, 25, 30, 40, 50]
    d = 8

    # Header
    header = f"  {'γ_bio':>6s}"
    for jmax in jmax_values:
        header += f"  J={jmax:>3.0f}"
    print(header)
    print("  " + "-" * (8 + 7 * len(jmax_values)))

    kill_triggered = False
    all_results = []

    for gamma_bio in gamma_bio_values:
        row = f"  {gamma_bio:6.0f}"
        for jmax in jmax_values:
            # Rebuild Hamiltonian with scaled couplings
            scale = jmax / 25.0  # scale relative to default J_nn=25
            H = np.diag(_CHROMATIN_SITE_ENERGIES_CM)
            for i in range(d - 1):
                H[i, i + 1] = _CHROMATIN_J_NN * scale
                H[i + 1, i] = _CHROMATIN_J_NN * scale
            H[0, 7] = _CHROMATIN_J_LR * scale
            H[7, 0] = _CHROMATIN_J_LR * scale

            # Compute θ_min at this (γ_bio, J_max)
            rho0 = np.zeros((d, d), dtype=complex)
            rho0[0, 0] = 1.0
            rho = evolve_lindblad_sink_generic(
                H, gamma_bio, [3, 4], [0.5, 0.5],
                rho0, t_final_ps=5.0, dt_fs=1.0)
            tr = np.real(np.trace(rho))
            if not np.all(np.isfinite(rho)) or tr < 1e-12:
                theta_min = 90.0
            else:
                rho_n = rho / tr
                rho_reg = rho_n + 1e-10 * np.eye(d) / d
                rho_reg /= np.trace(rho_reg)
                angles_deg = np.sort(np.degrees(principal_angles(rho_reg)))
                theta_min = angles_deg[0]

            if theta_min < 85.0:
                kill_triggered = True
            all_results.append((gamma_bio, jmax, theta_min))
            row += f"  {theta_min:5.1f}°"
        print(row)

    ratio_row = f"  {'γ/J':>6s}"
    for jmax in jmax_values:
        for gamma_bio in [200]:  # reference
            ratio_row += f"  {gamma_bio/jmax:5.1f} "
    print(f"\n  Reference ratios (γ_bio=200): ", end="")
    for jmax in jmax_values:
        print(f" {200/jmax:5.1f}", end="")
    print()

    if kill_triggered:
        print("\n  *** KILL CRITERION TRIGGERED: θ_min < 85° detected ***")
        print("  *** Chromatin should go to appendix only ***")
    else:
        print(f"\n  Kill criterion passed: all θ_min > 85°")
        min_overall = min(r[2] for r in all_results)
        print(f"  Minimum θ_min across all combinations: {min_overall:.1f}°")

    return not kill_triggered, all_results


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print(f"Output directory: {FIGDIR}\n")
    fig1_bloch_collapse()
    fig2_fmo_graph()
    fig3_bures_angle()
    fig4_fmo_collapse()
    fig5_qudit_angles()
    fig6_robustness()
    fig7_pe545()
    fig8_neural_comparison()
    pe545_robustness_check()
    neural_robustness_check()
    chromatin_sensitivity_sweep()
    print("\nAll figures generated.")
