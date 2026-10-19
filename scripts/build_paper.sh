#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
python3 analysis/make_paper_assets.py
if command -v tectonic >/dev/null 2>&1; then
  tectonic --keep-logs paper.tex
else
  pdflatex -interaction=nonstopmode -halt-on-error paper.tex
  bibtex paper
  pdflatex -interaction=nonstopmode -halt-on-error paper.tex
  pdflatex -interaction=nonstopmode -halt-on-error paper.tex
fi
