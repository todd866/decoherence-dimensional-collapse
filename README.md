# Information Geometry of the Quantum-Classical Boundary in Biological Information Processing

**Decoherence is a dimensional collapse.** A $d$-level quantum system lives on a $(d^2{-}1)$-dimensional state manifold; dephasing projects it onto the $(d{-}1)$-dimensional classical simplex. This paper characterises the geometry of that projection exactly, then applies it across biological scales to show where different systems sit on the quantum-classical spectrum.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## The Idea

The quantum-classical transition is usually described by timescales (how fast does coherence decay?) or scalar measures (how much coherence remains?). Both compress a high-dimensional process into a single number.

We use **principal angles** between the diagonal (population) and off-diagonal (coherence) tangent subspaces under the Bures metric. This gives $d{-}1$ angles --- a full geometric characterisation of how "quantum" a state is. At $0°$, populations and coherences are maximally entangled in the information geometry; at $90°$, the system is indistinguishable from a classical probability distribution.

$$\Lambda_{\mathrm{dec}}(\rho) = \sum_{k=1}^{d} |k\rangle\langle k|\,\rho\,|k\rangle\langle k| \qquad \Longrightarrow \qquad d^2{-}1 \;\;\xrightarrow{\;\Lambda_{\mathrm{dec}}\;}\;\; d{-}1 \text{ dimensions}$$

The framework is scale-invariant: principal angles depend on the dimensionless ratio $\gamma/J_{\max}$ (dephasing to coupling), not absolute energy. Any $d$-level system --- photosynthetic antenna, molecular complex, neural oscillator assembly --- can be placed on the same axis.

## Results

### Theory (any $d$) --- **Proved**
- **Dimensional collapse**: QFIM rank drops from $d^2{-}1$ to $d{-}1$; collapse ratio $1/(d{+}1) \to 0$
- **Exact Bures angle** (qubits): closed-form expression for the separation angle
- **Global orthogonality theorem** ($d \geq 3$): diagonal and off-diagonal subspaces are Bures-orthogonal *if and only if* $\rho$ is diagonal --- any coherence breaks orthogonality
- **Non-classicality index**: $\chi(\rho) = \frac{1}{d-1}\sum_a \cos^2\theta_a$ --- continuous scalar measure of departure from classical simplex
- **Thermodynamic cost**: coherences store free energy $k_BT \cdot C_r(\rho)$; cost per collapsed dimension $\sim \ln d / d^2$

### Applications --- **Empirical**

| System | $d$ | Collapse | $\gamma_{\rm bio}/J_{\max}$ | $\theta_{\min}$ at bio | Regime |
|--------|-----|----------|-----------------------------|------------------------|--------|
| **FMO** (photosynthesis) | 7 | $48 \to 6$ | ~1.1 | $87°$ | Near transition --- quantum effects exploitable |
| **PE545** (photosynthesis) | 8 | $63 \to 7$ | ~7.9 | $89°$ | Just past optimum --- still in transition zone |
| **Cortical gamma** (neural) | 10 | $99 \to 9$ | ~$8 \times 10^{12}$ | $\to 90°$ | Deep classical --- 12 orders past the boundary |

An illustrative molecular-oscillator case study (chromatin, $d = 8$) is presented in Appendix C. At reference parameters ($\gamma_{\rm bio}/J_{\max} = 8.0$), $\theta_{\min} \approx 87°$, but the position is sensitive to coupling strength (see sensitivity sweep in Appendix C).

### Multi-scale architecture --- **Hypothesis**
Neural oscillations are deep classical ($\chi_{\rm osc} \approx 0$), but the molecular machinery they coordinate --- ion channels ($d \sim 4$), synaptic receptors --- operates at energy scales where $\chi_m > 0$. Coordinated dimensionality $D_{\rm coord} = N[(d_m{-}1) + \chi_m(d_m^2{-}d_m)]$ scales with entrained site count $N$, not oscillation frequency. Low frequency = high dimensionality.

### Cost-benefit interpretation --- **Hypothesis**
The ENAQT regime may represent a cost-benefit sweet spot where extra dimensional capacity from partial coherence comes at modest thermodynamic cost.

## Reader's Guide

| If you want to... | Read... |
|---|---|
| Understand the core idea | §1 Introduction + Figure 1 (Bloch sphere collapse) |
| See the mathematical framework | §2 (state manifold, QFIM, Bures metric, effective dimensionality) |
| Read the main theorems | §3: Thm 3.1 (collapse), Thm 3.4 (Bures angle), Thm 3.6 (global orthogonality), Def 3.11 ($\chi$) |
| See the thermodynamic cost | §4 (free energy of coherence, cost per collapsed dimension) |
| See the FMO application | §5.1--5.5 + Figures 4, 6 |
| See the PE545 second complex | §5.7 + Figure 7 |
| See the neural application | §6 + Figure 8 |
| See the multi-scale architecture | §6.4: Def 6.1 ($D_{\rm coord}$), Prop 6.2 (low-F = high-D) |
| See cross-scale comparison | §6.4 + Table in Appendix D |
| See the chromatin case study | Appendix C (illustrative model + sensitivity sweep) |
| Understand robustness | §5.6 + Figure 6b (trapping-rate robustness) |
| Check the numerics | Appendix A (RK4, SLD, principal angle algorithm) |

## Repository Structure

```
├── decoherence_biology.tex             # Main manuscript
├── decoherence_biology.pdf             # Compiled paper (8 figures)
├── references.bib                      # Bibliography
├── cover_letter.tex                    # BioSystems cover letter
├── reader_guide.tex                    # Companion reader's guide
├── figures/
│   ├── fig1_bloch_collapse.pdf         # Bloch ball dimensional collapse
│   ├── fig2_fmo_graph.pdf              # FMO coupling graph
│   ├── fig3_bures_angle.pdf            # Bures angle across Bloch disk
│   ├── fig4_fmo_collapse.pdf           # FMO: transport efficiency + principal angles
│   ├── fig5_qudit_angles.pdf           # Qudit principal angles (d=3,5,7)
│   ├── fig6_robustness.pdf             # Angle spectrum + sink-rate robustness
│   ├── fig7_pe545.pdf                  # PE545: transport efficiency + principal angles
│   └── fig8_neural_comparison.pdf      # Cross-scale quantum-classical boundary
├── code/
│   ├── fmo_analysis.py                 # Complex registry, Lindblad dynamics, transport, geometry
│   └── generate_figures.py             # All figure generation + sensitivity sweeps
├── archive/
│   ├── FoP/                            # Foundations of Physics submission version
│   └── NJP/                            # New Journal of Physics submission version
└── README.md
```

## Reproducing the Paper

```bash
# Generate all 8 figures + robustness checks + sensitivity sweep
pip install numpy scipy matplotlib
cd code
python generate_figures.py              # Outputs PDFs to ../figures/

# Compile manuscript
cd ..
pdflatex decoherence_biology.tex
bibtex decoherence_biology
pdflatex decoherence_biology.tex
pdflatex decoherence_biology.tex
```

## Paper

**[Information Geometry of the Quantum-Classical Boundary in Biological Information Processing](decoherence_biology.pdf)**

Todd, I. (2026). Target: *BioSystems*.

## Citation

```bibtex
@article{todd2026information,
  author  = {Todd, Ian},
  title   = {Information Geometry of the Quantum-Classical Boundary
             in Biological Information Processing},
  journal = {BioSystems},
  year    = {2026},
  note    = {Submitted}
}
```

## License

MIT License
