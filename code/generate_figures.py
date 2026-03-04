#!/usr/bin/env python3
"""
Generate all figures for:
"Information Geometry of the Quantum-Classical Transition in Biological Systems"

Produces:
  figures/fig1_bloch_collapse.{pdf,png}     — Bloch sphere collapse (reused)
  figures/fig2_fmo_graph.{pdf,png}          — NEW: FMO coherence graph
  figures/fig3_bures_angle.{pdf,png}        — Bures angle across Bloch disk (reused)
  figures/fig4_fmo_collapse.{pdf,png}       — NEW: FMO dimensional collapse trajectory
  figures/fig5_two_stage.{pdf,png}          — NEW: Two-stage collapse schematic
  figures/fig6_qudit_angles.{pdf,png}       — Qudit principal angles (updated: d=7)

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
    evolve_lindblad, principal_angles,
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
# Figure 2 — FMO coherence graph with spectator structure (NEW)
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
    """FMO Bures principal angles as a function of dephasing rate.

    At each dephasing rate γ, evolve the FMO density matrix for 5 ps,
    then compute the d−1 = 6 principal angles between the diagonal
    (classical) and off-diagonal (quantum) tangent subspaces under the
    Bures metric.  At weak dephasing, angles are small (strong quantum-
    classical mixing); at strong dephasing, angles approach 90° (classical
    orthogonality, i.e. the quantum-classical split is clean).
    """
    print("Figure 4: FMO Bures principal angles vs dephasing ...")

    H = H_FMO * CM_TO_RADFS
    d = D_FMO
    rho0 = np.zeros((d, d), dtype=complex)
    rho0[0, 0] = 1.0

    gamma_values_cm = np.concatenate([
        np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0]),
        np.arange(10, 110, 10),
        np.arange(150, 600, 50),
    ])
    gamma_values_cm = np.sort(np.unique(gamma_values_cm))

    all_angles = []   # list of arrays, each shape (6,)

    for gamma_cm in gamma_values_cm:
        gamma = gamma_cm * CM_TO_RADFS
        rho = evolve_lindblad(H, gamma, rho0, t_final=5000.0, dt=1.0)
        rho_reg = rho + 1e-10 * np.eye(d) / d
        rho_reg /= np.trace(rho_reg)

        angles_rad = principal_angles(rho_reg)
        angles_deg = np.sort(np.degrees(angles_rad))
        all_angles.append(angles_deg)

        print(f"  gamma={gamma_cm:7.1f} cm^-1  "
              f"min={angles_deg[0]:5.1f}°  max={angles_deg[-1]:5.1f}°")

    all_angles = np.array(all_angles)   # shape (n_gamma, 6)
    min_angles = all_angles[:, 0]
    max_angles = all_angles[:, -1]

    fig, ax = plt.subplots(figsize=(6.0, 4.0))

    # Shaded band: range of all 6 principal angles
    ax.fill_between(gamma_values_cm, min_angles, max_angles,
                     alpha=0.2, color=PAL[0])

    # Min and max principal angle curves
    ax.semilogx(gamma_values_cm, min_angles, "o-", color=PAL[0], lw=1.5,
                ms=3, markeredgecolor="white", markeredgewidth=0.4,
                label=r"$\min_a\, \theta_a$")
    ax.semilogx(gamma_values_cm, max_angles, "s-", color=PAL[2], lw=1.0,
                ms=3, markeredgecolor="white", markeredgewidth=0.4,
                label=r"$\max_a\, \theta_a$", alpha=0.7)

    ax.set_xlabel(r"Dephasing rate $\gamma$ (cm$^{-1}$)")
    ax.set_ylabel(r"Principal angle $\theta$ (degrees)")
    ax.set_xlim(0.08, 600)
    ax.set_ylim(0, 95)

    # Classical orthogonality reference
    ax.axhline(90, color=PAL[5], ls="--", lw=0.8, alpha=0.5)
    ax.text(0.12, 91, r"$90°$ (classical orthogonality)", fontsize=8,
            color=PAL[5], alpha=0.7)

    # Physiological dephasing range
    ax.axvspan(50, 200, alpha=0.08, color=PAL[1])
    ax.text(75, 10, "Physiological\ndephasing", fontsize=8,
            color=PAL[1], style="italic")

    ax.legend(loc="center right", fontsize=9, framealpha=0.9)

    fig.tight_layout()
    save(fig, "fig4_fmo_collapse")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 5 — Two-stage dimensional collapse schematic (NEW)
# ══════════════════════════════════════════════════════════════════════════════

def fig5_two_stage():
    print("Figure 5: Two-stage collapse schematic ...")

    fig, ax = plt.subplots(figsize=(7.0, 3.0))

    # Stage 1: Quantum -> Classical simplex
    # Stage 2: Classical simplex -> Effective biological dimensions

    # Draw boxes
    boxes = [
        (0.5, 0.5, 1.8, 1.5, "Quantum\nstate space\n" + r"$d^2 - 1$ dim",
         PAL[4], 0.15),
        (3.5, 0.5, 1.8, 1.5, "Classical\nsimplex\n" + r"$d - 1$ dim",
         PAL[0], 0.15),
        (6.5, 0.5, 1.8, 1.5, "Effective\nbiological\n" + r"$D_{\mathrm{eff}}$ dim",
         PAL[2], 0.15),
    ]

    for x, y, w, h, label, color, alpha in boxes:
        rect = mpatches.FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.15",
            facecolor=color, alpha=alpha, edgecolor=color, lw=2
        )
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, label,
                ha="center", va="center", fontsize=10, color=color,
                fontweight="bold")

    # Arrows between stages
    arrow_props = dict(arrowstyle="-|>", color="k", lw=2, mutation_scale=20)
    ax.annotate("", xy=(3.4, 1.25), xytext=(2.4, 1.25),
                arrowprops=arrow_props)
    ax.annotate("", xy=(6.4, 1.25), xytext=(5.4, 1.25),
                arrowprops=arrow_props)

    # Stage labels
    ax.text(2.9, 2.2, "Stage 1", fontsize=10, ha="center",
            fontweight="bold", color=PAL[3])
    ax.text(2.9, 1.85, "Decoherence\n(fs)", fontsize=8, ha="center",
            color=PAL[3])

    ax.text(5.9, 2.2, "Stage 2", fontsize=10, ha="center",
            fontweight="bold", color=PAL[2])
    ax.text(5.9, 1.85, "Dissipation\n(ms\u2013s)", fontsize=8, ha="center",
            color=PAL[2])

    # Dimension annotations
    ax.text(1.4, 0.25, "e.g. FMO: 48 dim", fontsize=7, ha="center",
            color=PAL[5], style="italic")
    ax.text(4.4, 0.25, "e.g. FMO: 6 dim", fontsize=7, ha="center",
            color=PAL[5], style="italic")

    # Metric annotation
    ax.text(4.4, -0.1, "Bures " + r"$\to$" + " Fisher\u2013Rao",
            fontsize=8, ha="center", color=PAL[5])

    ax.set_xlim(0, 9)
    ax.set_ylim(-0.4, 2.7)
    ax.axis("off")

    fig.tight_layout()
    save(fig, "fig5_two_stage")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 6 — Qudit principal angles: numerical study (updated with d=7)
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


def fig6_qudit_angles():
    print("Figure 6: Qudit principal angles (including d=7 for FMO) ...")

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
    save(fig, "fig6_qudit_angles")


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print(f"Output directory: {FIGDIR}\n")
    fig1_bloch_collapse()
    fig2_fmo_graph()
    fig3_bures_angle()
    fig4_fmo_collapse()
    fig5_two_stage()
    fig6_qudit_angles()
    print("\nAll figures generated.")
