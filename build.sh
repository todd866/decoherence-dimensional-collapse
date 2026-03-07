#!/bin/zsh
set -euo pipefail

refresh_anchors=0

needs_refresh() {
  local script_path="$1"
  shift
  local output_path
  for output_path in "$@"; do
    if [[ ! -f "${output_path}" || "${script_path}" -nt "${output_path}" ]]; then
      return 0
    fi
  done
  return 1
}

if [[ "${REBUILD_ANCHORS:-0}" == "1" ]] || needs_refresh code/photosynthetic_anchors.py results/photosynthetic_anchor_points.csv results/photosynthetic_anchor_runs.csv; then
  python3 code/photosynthetic_anchors.py
  refresh_anchors=1
fi

if [[ "${REBUILD_ANCHORS:-0}" == "1" ]] || needs_refresh code/protein_microdomain_payload.py results/protein_microdomain_summary.json results/protein_microdomain_scan.csv; then
  python3 code/protein_microdomain_payload.py
  refresh_anchors=1
fi

if [[ "${REBUILD_ANCHORS:-0}" == "1" ]] || needs_refresh code/neural_carrier_proxy.py results/neural_carrier_proxy.csv; then
  python3 code/neural_carrier_proxy.py
  refresh_anchors=1
fi

if [[ "${REBUILD_ANCHORS:-0}" == "1" || "${refresh_anchors}" == "1" ]] || needs_refresh code/ion_channel_payload.py results/ion_channel_summary.json results/biological_anchor_points.csv figures/ion_channel_anchor.pdf figures/threshold_entrainment_scaling.pdf figures/biological_anchor_map.pdf; then
  python3 code/ion_channel_payload.py
fi

pdflatex -interaction=nonstopmode paper.tex
bibtex paper
pdflatex -interaction=nonstopmode paper.tex
pdflatex -interaction=nonstopmode paper.tex
