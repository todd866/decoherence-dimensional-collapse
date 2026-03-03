#!/usr/bin/env python3
"""
Generate all figures for:
"Decoherence as Dimensional Collapse in Quantum State Space"

Produces:
  figures/fig1_bloch_collapse.{pdf,png}
  figures/fig2_classical_fraction.{pdf,png}
  figures/fig3_deff_trajectories.{pdf,png}
  figures/fig4_quantum_dlb.{pdf,png}

Usage:
  cd physics/70_decoherence_dimensional_collapse
  python3 code/generate_figures.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

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
    # Sequential colour map: blue -> red
    cmap = matplotlib.colormaps.get_cmap("coolwarm").resampled(len(gamma_vals))
    colours = [cmap(i) for i in range(len(gamma_vals))]

    fig = plt.figure(figsize=(7.0, 3.2))

    # ── Left panel: 3-D Bloch sphere ──────────────────────────────────────
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")

    # Sphere mesh
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(0, np.pi, 25)
    sx = np.outer(np.cos(u), np.sin(v))
    sy = np.outer(np.sin(u), np.sin(v))
    sz = np.outer(np.ones_like(u), np.cos(v))

    for idx, gamma in enumerate(gamma_vals):
        shrink = max(1.0 - 2.0 * gamma, 0.0)
        x = shrink * sx
        y = shrink * sy
        z = sz  # z-component unaffected
        alpha = 0.18 if idx < len(gamma_vals) - 1 else 0.35
        ax3d.plot_surface(
            x, y, z,
            color=colours[idx],
            alpha=alpha,
            linewidth=0.3,
            edgecolor=(*colours[idx][:3], 0.25),
            shade=False,
        )

    # z-axis reference line
    ax3d.plot([0, 0], [0, 0], [-1.05, 1.05], color="k", lw=0.6, ls="--")

    ax3d.set_xlim(-1.1, 1.1)
    ax3d.set_ylim(-1.1, 1.1)
    ax3d.set_zlim(-1.1, 1.1)
    ax3d.set_xlabel(r"$r_x$", labelpad=1)
    ax3d.set_ylabel(r"$r_y$", labelpad=1)
    ax3d.set_zlabel(r"$r_z$", labelpad=1)
    ax3d.set_title("Accessible state space", fontsize=11, pad=2)
    ax3d.view_init(elev=20, azim=-55)
    ax3d.tick_params(axis="both", labelsize=7, pad=0)

    # Build legend manually
    from matplotlib.patches import Patch
    legend_patches = [
        Patch(facecolor=colours[i], alpha=0.5,
              label=rf"$\gamma = {gamma_vals[i]:.2f}$")
        for i in range(len(gamma_vals))
    ]
    ax3d.legend(handles=legend_patches, loc="upper left", fontsize=8,
                framealpha=0.85, handlelength=1.0, borderpad=0.3)

    # ── Right panel: D_eff(gamma) trajectory ──────────────────────────────
    ax = fig.add_subplot(1, 2, 2)

    gamma = np.linspace(0, 0.5, 500)
    shrink = (1.0 - 2.0 * gamma)
    f_phi = shrink ** 2  # at theta = pi/2
    # D_eff = (1 + f_phi)^2 / (1 + f_phi^2)
    Deff = (1.0 + f_phi) ** 2 / (1.0 + f_phi ** 2)

    ax.plot(gamma, Deff, color=PAL[0], lw=2.0)

    # Mark the 4 gamma values
    for idx, gv in enumerate(gamma_vals):
        s = max(1.0 - 2.0 * gv, 0.0)
        fp = s ** 2
        de = (1.0 + fp) ** 2 / (1.0 + fp ** 2)
        ax.plot(gv, de, "o", color=colours[idx], ms=7, zorder=5,
                markeredgecolor="white", markeredgewidth=0.8)

    ax.set_xlabel(r"Dephasing parameter $\gamma$")
    ax.set_ylabel(r"$D_{\mathrm{eff}}^Q$")
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0.9, 2.05)
    ax.set_title(r"Effective dimension ($\theta = \pi/2$)", fontsize=11)

    # Horizontal reference
    ax.axhline(1.0, color=PAL[5], ls="--", lw=0.8, alpha=0.6)
    ax.text(0.42, 1.03, "Classical", fontsize=8, color=PAL[5])

    fig.tight_layout(w_pad=2.5)
    save(fig, "fig1_bloch_collapse")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 — Classical dimensionality fraction
# ══════════════════════════════════════════════════════════════════════════════

def fig2_classical_fraction():
    print("Figure 2: Classical fraction ...")

    N = np.arange(1, 21)
    d = 2 ** N
    frac = (d - 1.0) / (d ** 2 - 1.0)  # = 1/(d+1)

    fig, ax = plt.subplots(figsize=(5.0, 3.5))

    ax.semilogy(N, frac, "o-", color=PAL[0], lw=1.5, ms=5,
                markeredgecolor="white", markeredgewidth=0.6)

    # Annotations for key values
    annot = {1: "33.3%", 2: "20%", 3: "11.1%", 10: "0.098%", 20: r"$\approx 10^{-6}$"}
    for n_val, label in annot.items():
        idx = n_val - 1
        yval = frac[idx]
        # Alternate offset direction for readability
        xyoff = (12, 8) if n_val <= 3 else (12, -10)
        ax.annotate(
            label, xy=(n_val, yval), xytext=xyoff,
            textcoords="offset points", fontsize=8,
            arrowprops=dict(arrowstyle="-", color=PAL[5], lw=0.6),
            color=PAL[3],
        )

    # 1% reference line
    ax.axhline(0.01, color=PAL[5], ls="--", lw=0.8, alpha=0.6)
    ax.text(17.5, 0.012, "1%", fontsize=8, color=PAL[5])

    ax.set_xlabel(r"Number of qubits $N$")
    ax.set_ylabel(r"Classical fraction $(d{-}1)/(d^2{-}1)$")
    ax.set_xlim(0.5, 20.5)
    ax.set_xticks([1, 5, 10, 15, 20])

    # Secondary x-axis: d = 2^N
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    sec_ticks = [1, 5, 10, 15, 20]
    ax2.set_xticks(sec_ticks)
    ax2.set_xticklabels([rf"$2^{{{n}}}$" for n in sec_ticks], fontsize=9)
    ax2.set_xlabel(r"Hilbert-space dimension $d = 2^N$", fontsize=11)

    fig.tight_layout()
    save(fig, "fig2_classical_fraction")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 3 — D_eff trajectories under spin-boson decoherence
# ══════════════════════════════════════════════════════════════════════════════

def fig3_deff_trajectories():
    print("Figure 3: Spin-boson D_eff trajectories ...")

    t = np.linspace(0, 5, 500)
    coupling_vals = [0.01, 0.1, 0.5, 2.0]
    colours = [PAL[0], PAL[2], PAL[1], PAL[3]]

    fig, ax = plt.subplots(figsize=(5.0, 3.5))

    for eta_kT, col in zip(coupling_vals, colours):
        Gamma_t = eta_kT * t
        f_phi = np.exp(-2.0 * Gamma_t)
        Deff = (1.0 + f_phi) ** 2 / (1.0 + f_phi ** 2)
        ax.plot(t, Deff, color=col, lw=1.8,
                label=rf"$\eta\, k_BT = {eta_kT}$")

    # Classical limit
    ax.axhline(1.0, color=PAL[5], ls="--", lw=0.8, alpha=0.6)
    ax.text(4.05, 1.03, "Classical limit", fontsize=8, color=PAL[5])

    ax.set_xlabel(r"Time $t$ (normalised)")
    ax.set_ylabel(r"$D_{\mathrm{eff}}^Q(t)$")
    ax.set_xlim(0, 5)
    ax.set_ylim(0.9, 2.05)
    ax.legend(loc="upper right", framealpha=0.9, borderpad=0.4)

    fig.tight_layout()
    save(fig, "fig3_deff_trajectories")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 4 — Thermodynamic cost of decoherence (quantum DLB)
# ══════════════════════════════════════════════════════════════════════════════

def fig4_quantum_dlb():
    print("Figure 4: Coherent free energy ...")

    gamma = np.linspace(0, 0.5, 500)

    def binary_entropy(p):
        """H_2(p) = -p ln p - (1-p) ln(1-p), safe at boundaries."""
        p = np.clip(p, 1e-15, 1.0 - 1e-15)
        return -p * np.log(p) - (1.0 - p) * np.log(1.0 - p)

    # ── Curve 1: |+> state (maximally coherent pure state) ────────────────
    # After dephasing by gamma:
    #   eigenvalues of rho(gamma) = (1 +/- (1-2gamma))/2
    #   S(rho) = H_2((1 + (1-2gamma))/2)
    #   S(Lambda_dec[rho]) = ln 2  (diagonal is always 1/2, 1/2 for |+>)
    #   C_r = ln 2 - S(rho)
    p_plus = (1.0 + (1.0 - 2.0 * gamma)) / 2.0  # = 1 - gamma
    S_rho_plus = binary_entropy(p_plus)
    Cr_plus = np.log(2) - S_rho_plus

    # ── Curve 2: partially coherent mixed state ───────────────────────────
    # rho_0 = 0.8 |+><+| + 0.2 * I/2
    # = 0.8 * (1/2)(I + sigma_x) + 0.2 * (1/2)I
    # = (1/2)I + 0.4 sigma_x
    # So Bloch vector = (0.4, 0, 0)
    # After dephasing: r_x -> (1-2gamma)*0.4
    # Eigenvalues: (1 +/- |(1-2gamma)*0.4|)/2
    r_mixed = 0.4
    coherence_mixed = (1.0 - 2.0 * gamma) * r_mixed
    p_mixed = (1.0 + np.abs(coherence_mixed)) / 2.0
    S_rho_mixed = binary_entropy(p_mixed)
    # Diagonal of rho_0 is (1/2, 1/2) regardless of off-diag, so S(Lambda_dec[rho]) = ln 2
    Cr_mixed = np.log(2) - S_rho_mixed

    # ── Curve 3: already-diagonal state I/2 ───────────────────────────────
    Cr_diag = np.zeros_like(gamma)

    fig, ax = plt.subplots(figsize=(5.0, 3.5))

    ax.plot(gamma, Cr_plus, color=PAL[0], lw=2.0,
            label=r"$|+\rangle$ (pure, max coherence)")
    ax.plot(gamma, Cr_mixed, color=PAL[1], lw=1.8, ls="--",
            label=r"$0.8\,|{+}\rangle\!\langle{+}| + 0.2\,I/2$")
    ax.plot(gamma, Cr_diag, color=PAL[5], lw=1.5, ls=":",
            label=r"$I/2$ (diagonal)")

    ax.set_xlabel(r"Dephasing parameter $\gamma$")
    ax.set_ylabel(r"$\Delta F\!/\,k_BT = C_r$")
    ax.set_xlim(0, 0.5)
    ax.set_ylim(-0.03, 0.75)
    ax.legend(loc="upper right", framealpha=0.9, borderpad=0.4)

    # Reference: ln 2
    ax.axhline(np.log(2), color=PAL[5], ls="--", lw=0.6, alpha=0.5)
    ax.text(0.01, np.log(2) + 0.02, r"$\ln 2$", fontsize=8, color=PAL[5])

    fig.tight_layout()
    save(fig, "fig4_quantum_dlb")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 5 — Bures separation angle across the Bloch disk
# ══════════════════════════════════════════════════════════════════════════════

def fig5_bures_angle():
    print("Figure 5: Bures separation angle ...")

    n = 500
    rperp = np.linspace(0, 1, n)
    rz = np.linspace(-1, 1, n)
    Rp, Rz = np.meshgrid(rperp, rz)
    R2 = Rp**2 + Rz**2

    # Mask out exterior of Bloch ball (|r| >= 1)
    mask = R2 >= 1.0

    # sin^2(theta_B) = (1 - |r|^2) / [(1 - r_perp^2)(1 - r_z^2)]
    # Avoid division by zero at r_perp=1 or r_z=±1
    denom = np.clip((1.0 - Rp**2) * (1.0 - Rz**2), 1e-15, None)
    sin2 = np.clip((1.0 - R2) / denom, 0.0, 1.0)
    theta = np.degrees(np.arcsin(np.sqrt(sin2)))
    theta[mask] = np.nan

    fig, ax = plt.subplots(figsize=(5.5, 4.5))

    # Heatmap
    im = ax.pcolormesh(
        rperp, rz, theta,
        cmap="RdYlBu",  # red (low angle) to blue (90°)
        shading="auto",
        vmin=0, vmax=90,
        rasterized=True,
    )

    # Bloch ball boundary arc
    phi = np.linspace(-np.pi / 2, np.pi / 2, 300)
    ax.plot(np.cos(phi), np.sin(phi), "k-", lw=1.2, alpha=0.6)

    # Dephasing trajectories: horizontal arrows at fixed r_z
    traj_rz = [0.3, 0.5, 0.7, 0.9]
    for rz_val in traj_rz:
        rp_max = np.sqrt(1.0 - rz_val**2) * 0.92  # start inside boundary
        ax.annotate(
            "", xy=(0.02, rz_val), xytext=(rp_max, rz_val),
            arrowprops=dict(
                arrowstyle="->", color="k", lw=1.0,
                linestyle="dashed",
            ),
        )
        # Mirror to negative r_z
        ax.annotate(
            "", xy=(0.02, -rz_val), xytext=(rp_max, -rz_val),
            arrowprops=dict(
                arrowstyle="->", color="k", lw=1.0,
                linestyle="dashed",
            ),
        )

    cb = fig.colorbar(im, ax=ax, label=r"$\theta_B$ (degrees)", pad=0.02)
    cb.set_ticks([0, 15, 30, 45, 60, 75, 90])

    ax.set_xlabel(r"Coherence magnitude $r_\perp$")
    ax.set_ylabel(r"Population imbalance $r_z$")
    ax.set_xlim(0, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_aspect("equal")

    # Mark the axes where theta = 90
    ax.text(0.03, 0.88, r"$\theta_B = 90^\circ$", fontsize=8, color="k",
            fontweight="bold", alpha=0.8)
    ax.text(0.35, 0.06, r"$\theta_B = 90^\circ$", fontsize=8, color="k",
            fontweight="bold", alpha=0.8)

    fig.tight_layout()
    save(fig, "fig5_bures_angle")


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print(f"Output directory: {FIGDIR}\n")
    fig1_bloch_collapse()
    fig2_classical_fraction()
    fig3_deff_trajectories()
    fig4_quantum_dlb()
    fig5_bures_angle()
    print("\nAll figures generated.")
