#!/bin/bash
# Build all PDFs and clean up build artifacts
set -e

for tex in decoherence_biology cover_letter reader_guide; do
    [ -f "${tex}.tex" ] || continue
    pdflatex -interaction=nonstopmode "${tex}.tex" > /dev/null 2>&1
    if [ "$tex" = "decoherence_biology" ]; then
        bibtex "$tex" > /dev/null 2>&1
        pdflatex -interaction=nonstopmode "${tex}.tex" > /dev/null 2>&1
    fi
    pdflatex -interaction=nonstopmode "${tex}.tex" > /dev/null 2>&1
done

# Clean up
rm -f *.aux *.log *.out *.bbl *.blg *.toc *.lot *.lof *.fls *.fdb_latexmk *.synctex.gz

echo "Built: $(ls -1 *.pdf 2>/dev/null | tr '\n' ' ')"
