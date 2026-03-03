# Decoherence as Dimensional Collapse: The Geometry of the Quantum-Classical Transition

**Information-geometric characterization of decoherence as structured dimensional projection in quantum state space.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

Decoherence explains *that* the quantum-to-classical transition occurs, but its geometric structure has not been quantified. We show that dephasing is a dimensional collapse: a d-dimensional quantum system lives on a (d²−1)-dimensional manifold with the Bures metric, and dephasing projects it onto the (d−1)-dimensional classical simplex with the Fisher-Rao metric, annihilating the d²−d coherence directions. We compute the exact Bures principal angle between the classical and quantum tangent subspaces at arbitrary qubit states:

cos²θ_B = r_z² r_⊥² / [(1 − r_⊥²)(1 − r_z²)]

## Key Results

- **Dimensional collapse theorem**: Dephasing reduces the QFIM rank from d²−1 to at most d−1, annihilating specifically the coherence dimensions
- **Orthogonal decomposition**: At full-rank diagonal states, the differential of the dephasing map is the Bures-orthogonal projector onto the classical tangent subspace
- **Bures separation angle**: Closed-form formula for the angle between classical and quantum tangent subspaces at arbitrary qubit states; the split is oblique exactly when the state carries both coherence and population imbalance
- **Free energy decomposition**: Coherences store non-equilibrium free energy k_BT · C_r, which decoherence irreversibly destroys

## Running Simulations

```bash
cd code
python generate_figures.py      # Generate all publication figures
```

Requires: numpy, matplotlib

## Paper

**Decoherence as Dimensional Collapse: The Geometry of the Quantum-Classical Transition**

Todd, I. (2026). Target: *Foundations of Physics* (in preparation).

## Citation

```bibtex
@article{todd2026decoherence,
  author  = {Todd, Ian},
  title   = {Decoherence as Dimensional Collapse: The Geometry of the Quantum-Classical Transition},
  journal = {Foundations of Physics},
  year    = {2026},
  note    = {In preparation}
}
```

## License

MIT License
