#!/bin/zsh
set -euo pipefail

python3 code/ion_channel_payload.py
pdflatex -interaction=nonstopmode paper.tex
bibtex paper
pdflatex -interaction=nonstopmode paper.tex
pdflatex -interaction=nonstopmode paper.tex
