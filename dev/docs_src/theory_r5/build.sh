#!/bin/sh
# Build lsqvce_theory_r5.pdf from the LaTeX source (pdflatex, two passes for TOC/refs).
# Created 2026-09-01 (CEST, UTC+02:00) with assistance from Claude (Anthropic), model Claude Fable 5.1 -- REVIEW R5.
set -e
cd "$(dirname "$0")"
pdflatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null
mv main.pdf ../../../lsqvce_theory.pdf   # THE shipped theory PDF (R6 WP-101: one builder, one name)
rm -f main.aux main.log main.out main.toc
echo "built lsqvce_theory.pdf"
