# Information Geometry of the Quantum-Classical Boundary in Biological Information Processing

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

This project introduces a geometric way to quantify **how non-classical** a system is, then applies it across biological scales.

## What this repo is doing

Standard coherence metrics usually compress dynamics to one scalar without resolving geometry.  
Here, the core observable is the set of Bures principal angles between:
- classical tangent directions (population/diagonal)
- coherence tangent directions (off-diagonal)

From those angles, we define a scalar non-classicality index:

$$
\chi(\rho) \;=\; \frac{1}{d-1}\sum_{a=1}^{d-1}\cos^2\theta_a,
$$

where:
- $\chi = 0$ means geometrically classical (orthogonal split)
- $\chi > 0$ means non-classical mixing is present
- larger $\chi$ means stronger classical-quantum geometric coupling

So the claim is not binary (“quantum vs classical”).  
It is geometric and continuous: **deep classical**, **near-classical but non-classical**, or **strongly mixed**.

## Core geometric statement

For dephasing $\Lambda_{\mathrm{dec}}$:

$$
\Lambda_{\mathrm{dec}}(\rho) = \sum_{k=1}^{d} |k\rangle\langle k|\,\rho\,|k\rangle\langle k|,
$$

the state-space collapse is:

$$
d^2 - 1 \;\xrightarrow{\;\Lambda_{\mathrm{dec}}\;}\; d - 1.
$$

The framework is scale-invariant: results are controlled by $\gamma/J_{\max}$.

## What is proved vs modeled vs hypothesized

### Proved (mathematics)
- Dimensional collapse in QFIM rank: $d^2-1 \to d-1$
- Bures $\to$ Fisher-Rao metric restriction on the classical simplex
- Exact qubit Bures angle formula
- Global orthogonality criterion for $d \ge 3$ (orthogonal iff diagonal)
- Non-classicality index definition from principal angles
- Free-energy identity under dephasing

### Empirical / model-based computations

| System | $d$ | Collapse | $\gamma_{\rm bio}/J_{\max}$ | $\theta_{\min}$ at bio | Interpretation |
|---|---:|---:|---:|---:|---|
| FMO (validated) | 7 | $48 \to 6$ | $\sim 1.1$ | $\sim 87.7^\circ$ | Near transition |
| PE545 (validated) | 8 | $63 \to 7$ | $\sim 7.9$ | $\sim 89.0^\circ$ | Near-classical transition zone |
| Ion channel (illustrative) | 4 | $15 \to 3$ | $\sim 6.7$ | $\sim 88.5^\circ$ | Non-classical ($\chi \approx 2.8\times10^{-4}$), near-classical |
| Neural gamma (model) | 10 | $99 \to 9$ | $\sim 8\times 10^{12}$ | $\to 90^\circ$ | Deep classical at oscillation scale |

Chromatin ($d=8$) is included as an additional illustrative molecular model (sensitivity-limited at weak coupling).

### Hypothesis layer (explicitly tagged in manuscript)

Low-frequency oscillations are deep classical carriers ($\chi_{\rm osc}\approx 0$), while coordinated molecular subsystems can remain non-classical ($\chi_m>0$).  
Coordinated dimensionality is defined as:

$$
D_{\mathrm{coord}} \;=\; N\left[(d_m-1)+\chi_m(d_m^2-d_m)\right].
$$

This formalizes “low frequency can coordinate high-dimensional non-classical content” through site count $N$, not oscillation frequency itself.

## Reader map

| Goal | Where to read |
|---|---|
| Core idea | §1 + Fig. 1 |
| Geometry framework | §2 |
| Main theorems | §3 (incl. $\chi$ definition) |
| Thermodynamics | §4 |
| Photosynthesis (FMO/PE545) | §5 + Figs. 4, 6, 7 |
| Neural + cross-scale | §6 + Fig. 8 |
| Low-F/high-D formalization | §6.4 (Def. 6.1, Prop. 6.2) |
| Ion channel model anchor | Appendix C.1 |
| Chromatin sensitivity | Appendix C |
| Numerics/methods | Appendix A |

## Build and reproduce

```bash
pip install numpy scipy matplotlib

cd code
python generate_figures.py

cd ..
pdflatex decoherence_biology.tex
bibtex decoherence_biology
pdflatex decoherence_biology.tex
pdflatex decoherence_biology.tex
```

## Repository layout

```text
decoherence_biology.tex          # main manuscript
reader_guide.tex                 # companion guide
cover_letter.tex                 # submission letter
references.bib                   # bibliography
code/fmo_analysis.py             # models + dynamics + geometry utilities
code/generate_figures.py         # figure generation + robustness/sweeps
figures/                         # generated figures
archive/FoP/                     # archived FoP version
archive/NJP/                     # archived NJP version
```

## Manuscript

- [Main paper](decoherence_biology.pdf)
- [Reader guide](reader_guide.pdf)

## Citation

```bibtex
@article{todd2026information,
  author  = {Todd, Ian},
  title   = {Information Geometry of the Quantum-Classical Boundary in Biological Information Processing},
  journal = {BioSystems},
  year    = {2026},
  note    = {Submitted}
}
```

## License

MIT
