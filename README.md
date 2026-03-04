# Information Geometry of the Quantum-Classical Transition in Biological Systems

**Information-geometric characterization of decoherence as dimensional collapse, applied to photosynthetic energy transfer in the FMO complex.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

Decoherence maps a d-dimensional quantum system from a (d²−1)-dimensional manifold with the Bures metric onto the (d−1)-dimensional classical simplex with the Fisher-Rao metric. We characterize this dimensional collapse exactly for qubits and partially for arbitrary dimension, then apply the framework to the 7-chromophore Fenna-Matthews-Olson (FMO) photosynthetic complex. At the ENAQT-optimal dephasing rate, the system occupies an intermediate position in the collapse from 48 quantum dimensions to 6 classical dimensions. The spectator theorem identifies specific chromophore pairs where the quantum-classical split is geometrically non-orthogonal during the transport window. After quantum decoherence completes (~fs), the surviving Fisher-Rao geometry governs classical biological coherence on ms-s timescales.

## Key Results

- **Dimensional collapse theorem**: Dephasing reduces the QFIM rank from d²−1 to at most d−1
- **Bures separation angle**: Exact closed-form for qubits; spectator criterion and perturbative formula for d ≥ 3
- **FMO application**: 26 spectator triples in the 7-site complex; coherence magnitude quantifies the quantum-classical transition under progressive dephasing
- **Two-stage collapse**: Quantum decoherence (fs) and classical dissipation (ms-s) connected through the same Fisher-Rao geometry
- **Thermodynamic cost**: Coherences store free energy k_BT · C_r, destroyed irreversibly by decoherence

## Repository Structure

```
├── decoherence_biology.tex             # Main manuscript (LaTeX)
├── decoherence_biology.pdf             # Compiled paper
├── references.bib                      # Bibliography
├── cover_letter.tex                    # NJP cover letter
├── figures/                            # Generated publication figures
│   ├── fig1_bloch_collapse.pdf         # Bloch ball dimensional collapse
│   ├── fig2_fmo_graph.pdf              # FMO coherence graph with spectators
│   ├── fig3_bures_angle.pdf            # Bures angle across Bloch ball
│   ├── fig4_fmo_collapse.pdf           # FMO coherence vs dephasing rate
│   ├── fig5_two_stage.pdf              # Two-stage collapse schematic
│   └── fig6_qudit_angles.pdf           # Qudit principal angles (d=3,5,7)
├── code/
│   ├── generate_figures.py             # Figure generation script
│   └── fmo_analysis.py                 # FMO-specific computations
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

**[Information Geometry of the Quantum-Classical Transition in Biological Systems](decoherence_biology.pdf)**

Todd, I. (2026). Target: *New Journal of Physics*.

## Citation

```bibtex
@article{todd2026information,
  author  = {Todd, Ian},
  title   = {Information Geometry of the Quantum-Classical Transition in Biological Systems},
  journal = {New Journal of Physics},
  year    = {2026},
  note    = {In preparation}
}
```

## License

MIT License
