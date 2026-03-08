# Information Geometry of Classical Carriers and Non-Classical Payloads in Biological Information Processing

**A geometric framework for the quantum-classical boundary in biology, with implications for hybrid simulation and computing architecture.**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

Decoherence maps a $d$-level quantum system from a $(d^2-1)$-dimensional Bures manifold onto a $(d-1)$-dimensional Fisher--Rao simplex. This dimensional collapse yields a continuous non-classicality index

$$
\chi(\rho) = \frac{1}{d-1}\sum_{a=1}^{d-1}\cos^2\theta_a \in [0,1],
$$

from Bures principal angles between classical and coherence tangent subspaces. Applied to biology, $\chi$ reveals a *carrier-payload architecture*: deep-classical macroscopic carriers ($\chi \approx 0$) coordinate an extensive near-boundary molecular substrate ($\chi > 0$).

## Key Results

- Photosynthetic complexes (FMO, PE545), ion-channel selectivity filters, and protein microdomains all inhabit a near-boundary band ($\chi \sim 10^{-4}$), while neural carrier proxies sit at $\chi \sim 10^{-5}$ to $10^{-6}$.
- The total non-classical payload grows linearly with entrained module count under nonzero-density and bounded-redundancy assumptions.
- A hybrid simulation protocol follows: classical carrier dynamics plus small Lindblad equations for near-boundary payload modules, scaling linearly when inter-module quantum correlations are negligible.
- The architecture suggests a speculative engineering direction in which classical oscillatory coordination replaces qubit isolation.

## Running Simulations

```bash
# Build PDF with cached anchor data
./build.sh

# Force full numerical recomputation
REBUILD_ANCHORS=1 ./build.sh
```

Individual analysis scripts:

```bash
python code/photosynthetic_anchors.py    # FMO/PE545 proof-of-principle
python code/ion_channel_payload.py       # Ion-channel selectivity-filter anchor
python code/protein_microdomain_payload.py  # Protein mid-fold payload anchor
python code/neural_carrier_proxy.py      # Reduced neural carrier proxies
```

## Paper

**Information Geometry of Classical Carriers and Non-Classical Payloads in Biological Information Processing**

Ian Todd, Sydney Medical School, University of Sydney

- **Target:** BioSystems
- **Status:** In preparation
- **Companion:** [github.com/todd866/decoherence-dimensional-collapse](https://github.com/todd866/decoherence-dimensional-collapse)

## Citation

```bibtex
@article{todd2026carrierPayload,
  author  = {Todd, Ian},
  title   = {Information Geometry of Classical Carriers and Non-Classical
             Payloads in Biological Information Processing},
  journal = {BioSystems},
  year    = {2026},
  note    = {In preparation}
}
```

## License

MIT
