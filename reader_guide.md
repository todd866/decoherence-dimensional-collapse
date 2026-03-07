# Reader Guide

## What this draft is trying to do

This is not a generic quantum-biology paper.

It is trying to make one argument precise:

> low-frequency waves can support high effective dimensionality because they are reliable classical carriers of large non-classical molecular payloads, with ion channels as the primary neural payload candidate.

## How to read it

### If you want the thesis fast

Read:

1. [paper.tex](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/paper.tex), Abstract
2. Section 1: Introduction
3. Section 3: Carrier-Payload Architecture
4. Section 5.3: minimal soft phase-coupling carrier model
5. Section 5.4: rival predictions table
6. Section 5: Low-Frequency Coordination and Intelligence
7. Section 6: Discussion

### If you want the geometry

Read:

1. Section 2.1: Dimensional collapse
2. Section 2.2: Non-classicality index
3. Section 3.1: Local payload load
4. Section 3.2: Extensive payload theorem
5. Section 3.3: Low-frequency coordination corollary

### If you want the biological anchors

Read:

1. Section 4.1: Photosynthesis
2. Section 4.2: Ion channels as the primary payload example
3. Section 4.3: Protein microdomains as secondary payloads
4. Table 1: anchor comparison
5. Figure 1: biological anchor map
6. Figure 2: ion-channel payload anchor
7. Figure 3: carrier/payload schematic
8. Figure 4: soft phase-coupling scaling schematic
9. Table 2: payload benchmark table
10. Section 5.1: Neural oscillations as carriers
11. Table 3: rival predictions
12. Table 4: operational measurement route

The ion-channel subsection now contains one actual rebuilt-repo anchor:

$$
\gamma_{\mathrm{bio}}/J_{\max}=6.67,\quad
\theta_{\min}\approx 88.55^\circ,\quad
\chi\approx 2.76\times 10^{-4}.
$$

It also includes a compact robustness sweep over
$$
J_{\max}\in\{15,20,30,40\}\,\mathrm{cm}^{-1},
\qquad
\gamma\in\{100,150,200,250,300\}\,\mathrm{cm}^{-1},
$$
which keeps the payload in the same boundary-near band:
$$
86.43^\circ \le \theta_{\min} \le 88.93^\circ,
\qquad
1.53\times10^{-4}\le \chi \le 1.55\times10^{-3}.
$$

There is now a second, tighter robustness check at fixed `J_max = 30 cm^-1` and `γ_bio = 200 cm^-1`: the repo varies the inner-well asymmetry `Δ ∈ {10,15,20} cm^-1` and a next-nearest coupling `J_2 ∈ {0,5,10} cm^-1` in
[results/ion_channel_topology_sensitivity.csv](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/results/ion_channel_topology_sensitivity.csv).
That keeps the payload point in a narrow band:
$$
88.55^\circ \le \theta_{\min} \le 88.77^\circ,
\qquad
2.20\times10^{-4}\le \chi \le 2.77\times10^{-4}.
$$
So the ion-channel anchor is not tied to one exact 4-site topology.

Figure 1 is the proof-of-principle bridge: it places the rebuilt molecular payload anchors on one biological operating map and shows the rebuilt neural carrier proxy farther right and closer to the classical limit. That makes the carrier/payload split visible at a glance instead of leaving it distributed across text and tables.

Section 4.3 is no longer only a placeholder. The rebuilt repo now includes a secondary protein payload anchor:
$$
\gamma_{\mathrm{bio}}/J_{\max}=7.69,\quad
\theta_{\min}\approx 88.21^\circ,\quad
\chi\approx 3.66\times 10^{-4}.
$$
That matters because it broadens the payload class beyond ion conduction without displacing ion channels from the main neural route.

Figure 2 is the cleanest fresh rebuilt-repo object: it shows the bounded trap scan only to locate the modeled operating point, and then shows the principal-angle spectrum that actually matters for the argument.

Figure 3 is the intuition object: it shows the disagreement with the frequency-as-dimensionality view in one glance. High frequency can still mark local access, while low frequency can coordinate a larger payload.

Figure 4 is the bridge object: it makes the low-frequency claim explicit in one soft coupling model instead of leaving it as pure verbal intuition.

It now does one extra job that matters: the solid curves are a heterogeneous soft-field simulation, while the dashed curves are the analytic soft envelope. The low-frequency advantage survives that heterogeneity, so the inversion is not being smuggled in by a single sharp threshold.

There is now also a parameter-band check behind that figure. The repo writes
[results/threshold_scaling_sensitivity.csv](/Users/iantodd/Projects/highdimensional/physics/70_decoherence_dimensional_collapse/results/threshold_scaling_sensitivity.csv), which shows that the analytic bridge stays monotone across `s ∈ {1,2,3}`, `η ∈ {1,2,4}`, and `R_amp/R0 ∈ {2,4,8}`, while the heterogeneous 2D field check stays monotone across the corresponding `η` and `R_amp/R0` ranges. That matters because the low-`f` inversion is now robust to a real parameter band, not just one default bridge setting.

Table 2 is the interpretation object: it shows that the ion-channel uplift fraction is small, but the absolute payload still grows linearly with entrained site count.

Table 3 is the falsifiability object: it states exactly how the present carrier-payload view differs from a frequency-as-dimensionality view. If dimensionality tracks carrier frequency after controlling for payload size and payload diversity, this paper loses.

Table 4 is the execution object: it translates the theory variables into things a neuroscience paper could actually measure. The current draft no longer hides behind direct molecular tomography as a prerequisite. The test is conditional and practical: carrier coherence, coordinated territory, payload diversity, and observed manifold dimensionality.

Right after that, the paper now gives the minimal regression form of the wager:
$$
D_{\mathrm{obs}}=\beta_0+\beta_1\log \widetilde N+\beta_2 H_{\mathrm{payload}}+\beta_3 C_{\mathrm{carrier}}+\beta_4\log f+\varepsilon.
$$
That matters because it reduces the disagreement to one coefficient test. If `β4` stays dominant after the payload and carrier controls, this paper loses.

Section 5.1 now also has one rebuilt carrier-side proxy: the synthetic gamma microcircuit gives `θ_min ≈ 89.78°` and `χ ≈ 3.72×10^-6` already at `γ/J_max = 50`, so the carrier/payload split is no longer just verbal.

Read Table 1 with that in mind: it now shows the `χ` gap directly. The payload anchors sit at `χ ~ 10^-4`, while the carrier proxy is already down at `χ ~ 10^-6`. That quantitative separation is the paper’s cleanest rebuilt-repo evidence for “classical carrier, non-classical payload.”

## Claim tiers

- `proved`: dimensional collapse, metric restriction, geometric definitions, extensive-payload result under stated assumptions
- `model-based`: rebuilt photosynthetic proof-of-principle anchors from validated Hamiltonians, one fresh ion-channel payload anchor from the rebuilt repo, and a soft phase-coupling bridge linking carrier scale to payload scale
- `hypothesis`: low-frequency/high-dimensional intelligence architecture

## What is deliberately not claimed

- not “high dimensionality automatically implies quantumness”
- not “neural oscillations are quantum”
- not “classical models are impossible”

The claim is narrower:

- large biological systems can carry an extensive non-classical payload if a nonzero fraction of their molecular modules sit near the boundary
- low-frequency carriers matter because they can coordinate more of that payload
