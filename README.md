# Information Geometry of the Quantum-Classical Transition in Photosynthetic Exciton Transport

**Information-geometric characterization of decoherence as dimensional collapse, applied to the Fenna-Matthews-Olson photosynthetic complex.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

Decoherence maps a d-dimensional quantum system from a (d^2-1)-dimensional manifold with the Bures metric onto the (d-1)-dimensional classical simplex with the Fisher-Rao metric. We characterize this dimensional collapse exactly for qubits (closed-form Bures separation angle) and completely for d >= 3 (global orthogonality theorem: diagonal and off-diagonal tangent subspaces are Bures-orthogonal at every non-diagonal state). Applied to the 7-chromophore FMO complex, the collapse is 48 -> 6 dimensions. Using a Haken-Strobl model with reaction-centre trapping, we compute transport efficiency and the six Bures principal angles from the same sink-inclusive dynamics (geometry extracted from the renormalized conditional state). At the ENAQT optimum, all angles exceed 87 degrees: biologically optimal transport occurs in a regime that is geometrically close to classical.

## Key Results

- **Dimensional collapse theorem**: Dephasing reduces the QFIM rank from d^2-1 to at most d-1
- **Bures separation angle**: Exact closed-form for qubits; global orthogonality theorem and perturbative cross-Gram formula for d >= 3
- **FMO transport + geometry**: Same-model ENAQT efficiency peak and Bures principal angle analysis, both site-1 and site-6 initial conditions
- **Thermodynamic cost**: Coherences store free energy k_BT * C_r, destroyed irreversibly by decoherence

## Repository Structure

```
├── decoherence_biology.tex             # Main manuscript (LaTeX)
├── decoherence_biology.pdf             # Compiled paper
├── references.bib                      # Bibliography
├── cover_letter.tex                    # NJP cover letter
├── figures/                            # Generated publication figures
│   ├── fig1_bloch_collapse.pdf         # Bloch ball dimensional collapse
│   ├── fig2_fmo_graph.pdf              # FMO coupling graph
│   ├── fig3_bures_angle.pdf            # Bures angle across Bloch ball
│   ├── fig4_fmo_collapse.pdf           # Dual panel: transport efficiency + principal angles
│   ├── fig5_qudit_angles.pdf           # Qudit principal angles (d=3,5,7)
│   └── fig6_robustness.pdf            # Angle spectrum + sink-rate robustness
├── code/
│   ├── generate_figures.py             # Figure generation script
│   └── fmo_analysis.py                 # FMO Hamiltonian, Lindblad dynamics, transport efficiency
├── archive/fop_version/                # Original Foundations of Physics version
├── README.md
└── LICENSE
```

## Reproducing the Paper

### Generate figures

```bash
pip install numpy scipy matplotlib networkx
cd code
python generate_figures.py      # Outputs PDFs to ../figures/
```

### Compile manuscript

```bash
pdflatex decoherence_biology.tex
bibtex decoherence_biology
pdflatex decoherence_biology.tex
pdflatex decoherence_biology.tex
```

## Paper

**[Information Geometry of the Quantum-Classical Transition in Photosynthetic Exciton Transport](decoherence_biology.pdf)**

Todd, I. (2026). Target: *New Journal of Physics*.

## Citation

```bibtex
@article{todd2026information,
  author  = {Todd, Ian},
  title   = {Information Geometry of the Quantum-Classical Transition in Photosynthetic Exciton Transport},
  journal = {New Journal of Physics},
  year    = {2026},
  note    = {In preparation}
}
```

## License

MIT License
