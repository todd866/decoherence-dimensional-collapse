# Information Geometry of Classical Carriers and Non-Classical Payloads

This repository is a fresh rebuild of the biological information-geometry project.

The new paper is organized around one claim:

$$
\text{dimensionality lives in the coordinated payload, not in the carrier frequency.}
$$

The carrier/payload split is:

$$
\chi_{\mathrm{carrier}} \approx 0,
\qquad
\chi_{\mathrm{payload}} > 0,
$$

where $\chi$ is the non-classicality index derived from Bures principal angles.

## Core idea

Decoherence maps a $d$-level quantum system from a $(d^2-1)$-dimensional Bures state manifold onto a $(d-1)$-dimensional Fisher-Rao simplex:

$$
d^2 - 1 \longrightarrow d - 1.
$$

This gives a continuous measure of non-classicality:

$$
\chi(\rho) = \frac{1}{d-1}\sum_{a=1}^{d-1}\cos^2\theta_a,
$$

with:

- $\chi = 0$ for geometrically classical states
- $\chi > 0$ when coherence directions remain geometrically accessible

The new manuscript uses that geometry to argue:

1. Molecular biological systems can sit near the quantum-classical boundary.
2. Macroscopic neural oscillations are deep classical at carrier scale.
3. Low-frequency waves can still support higher effective dimensionality by coordinating larger non-classical molecular payloads, with ion channels as the primary neural payload example.

## Ion-channel anchor

The rebuilt repo now includes one concrete neural payload computation in
[code/ion_channel_payload.py](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/code/ion_channel_payload.py).

For the minimal KcsA-like 4-site selectivity-filter model:

$$
\gamma_{\mathrm{bio}}/J_{\max} = 6.67,
\qquad
\theta_{\min} \approx 88.55^\circ,
\qquad
\chi \approx 2.76\times 10^{-4}.
$$

So the ion-channel payload is not strongly quantum and not fully classical. It is
boundary-near and slightly non-classical in exactly the sense the paper needs.

The generated summary files are:

- [results/ion_channel_summary.json](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/results/ion_channel_summary.json)
- [results/ion_channel_scan.csv](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/results/ion_channel_scan.csv)
- [results/ion_payload_benchmarks.csv](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/results/ion_payload_benchmarks.csv)
- [results/threshold_scaling_scan.csv](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/results/threshold_scaling_scan.csv)
- [results/biological_anchor_points.csv](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/results/biological_anchor_points.csv)
- [figures/biological_anchor_map.pdf](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/figures/biological_anchor_map.pdf)
- [figures/ion_channel_anchor.pdf](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/figures/ion_channel_anchor.pdf)
- [figures/threshold_entrainment_scaling.pdf](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/figures/threshold_entrainment_scaling.pdf)

![Biological anchor map](figures/biological_anchor_map.png)
![Ion-channel anchor](figures/ion_channel_anchor.png)

## Current anchor set

| System | Role | Source | `γ_bio/J_max` | `θ_min^bio` | Current use |
|---|---|---|---:|---:|---|
| FMO | Proof of principle | Legacy validated | 1.14 | 87.7° | Functional near-boundary photosynthesis anchor |
| PE545 | Proof of principle | Legacy validated | 7.88 | 89.0° | Independent photosynthetic confirmation of the same band |
| Ion channel | Primary neural payload | Rebuilt computation | 6.67 | 88.55° | Fresh rebuilt-repo molecular anchor |
| Protein microdomain | Secondary payload | Planned | — | — | Payload extension, not yet recomputed here |
| Neural oscillation | Carrier scale | Theory | `≫1` | `→90°` | Deep-classical coordination process |

The new anchor map is the compact visual summary: FMO, PE545, and the rebuilt ion-channel point all sit in the same near-boundary band. That is the concrete reason photosynthesis remains in the paper. It is the proof of principle that functional biology can inhabit this regime, and the ion-channel point then shows that a neural payload candidate can inhabit it too.

## Main scaling objects

For module $i$ with local dimension $d_i$ and non-classicality $\chi_i$, define the local non-classical load

$$
L_i = \chi_i(d_i^2-d_i).
$$

For an entrained set of modules $E$, define the redundancy-adjusted payload

$$
\mathcal{L}_{\mathrm{eff}}(E)
=
\frac{1}{\varrho_E}\sum_{i\in E} L_i,
\qquad
\varrho_E \ge 1,
$$

and the coordinated dimensionality

$$
D_{\mathrm{coord}}(E)
=
\sum_{i\in E}(d_i-1) + \mathcal{L}_{\mathrm{eff}}(E).
$$

This is the paper's working thesis:

$$
f \downarrow \;\Longrightarrow\; |E(f)| \uparrow
\;\Longrightarrow\;
D_{\mathrm{coord}}(E(f)) \uparrow,
$$

provided the payload statistics remain comparable and redundancy stays bounded.

The manuscript now also includes a minimal carrier-to-payload coupling bridge:

$$
A_f(r)=A_0 e^{-r/\ell(f)},
\qquad
\ell(f)=\ell_0\left(\frac{f_0}{f}\right)^\alpha,
$$

so that, in the threshold-entrainment toy model,

$$
D_{\mathrm{coord}}(E(f)) \propto \left(\frac{f_0}{f}\right)^{\alpha s}.
$$

This is not a fitted cortical law. It is the first explicit model in the rebuilt draft showing how lower-frequency carriers can yield higher coordinated dimensionality while remaining classical.

Using the actual ion-channel anchor, the local payload load is

$$
L_{\mathrm{ion}}=\chi(d^2-d)\approx 3.31\times 10^{-3}
$$

per channel at the biological point. The new threshold-scaling figure then plots

$$
\mathcal{L}_{\mathrm{eff}}(f)\approx \frac{N_0L_{\mathrm{ion}}}{\varrho}
\left(\frac{f_0}{f}\right)^\beta
$$

for illustrative redundancy assumptions, so the bridge figure now carries actual rebuilt-repo numbers rather than only a normalized schematic.

For homogeneous payload modules, the fractional uplift above the classical baseline is

$$
\Phi = \frac{\mathcal{L}_{\mathrm{eff}}}{N(d_m-1)}=\frac{\chi d_m}{\varrho}.
$$

For the ion-channel anchor \((d_m=4)\),

$$
\Phi_{\mathrm{ion}}\approx \frac{1.10\times 10^{-3}}{\varrho}.
$$

So the rebuilt draft is making a sober claim: the per-module uplift is small, but the absolute payload still scales linearly with entrained site count.

The manuscript now includes a compact benchmark table and the same numbers are written to:

- [results/ion_payload_benchmarks.csv](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/results/ion_payload_benchmarks.csv)

![Threshold-entrainment scaling](figures/threshold_entrainment_scaling.png)

## Files

- [paper.tex](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/paper.tex): main manuscript
- [references.bib](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/references.bib): bibliography
- [reader_guide.md](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/reader_guide.md): section-by-section reading guide
- [build.sh](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/build.sh): local LaTeX build script
- [code/ion_channel_payload.py](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/code/ion_channel_payload.py): minimal ion-channel payload analysis

## Build

```bash
./build.sh
```

`build.sh` regenerates the ion-channel anchor before compiling the manuscript.

This produces:

- `paper.pdf`
- `results/ion_channel_summary.json`
- `results/ion_channel_scan.csv`
- `results/ion_payload_benchmarks.csv`
- `results/threshold_scaling_scan.csv`
- `results/biological_anchor_points.csv`
- `figures/biological_anchor_map.pdf`
- `figures/ion_channel_anchor.pdf`
- `figures/threshold_entrainment_scaling.pdf`

## Current scope

This repo is intentionally narrow at the moment:

- manuscript-first rebuild
- theory and positioning first
- numerics/code to be reintroduced selectively from the tagged checkpoint, starting with ion-channel payload models

The old overhaul is preserved at git tag `big-overhaul-checkpoint`.
