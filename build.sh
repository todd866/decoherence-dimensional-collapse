#!/bin/zsh
set -euo pipefail

refresh_anchors=0

if [[ "${REBUILD_ANCHORS:-0}" == "1" || ! -f results/photosynthetic_anchor_points.csv || ! -f results/photosynthetic_anchor_runs.csv ]]; then
  python3 code/photosynthetic_anchors.py
  refresh_anchors=1
fi

if [[ "${REBUILD_ANCHORS:-0}" == "1" || "${refresh_anchors}" == "1" || ! -f results/ion_channel_summary.json || ! -f results/biological_anchor_points.csv ]]; then
  python3 code/ion_channel_payload.py
fi

pdflatex -interaction=nonstopmode paper.tex
bibtex paper
pdflatex -interaction=nonstopmode paper.tex
pdflatex -interaction=nonstopmode paper.tex
