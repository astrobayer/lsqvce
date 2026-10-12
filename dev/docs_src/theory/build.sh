#!/bin/sh
# Build lsqvce_theory.pdf -- the mathematical framework of lsqvce -- from the
# LaTeX source in this folder (pdflatex, three passes: the contents and the
# cross-references settle). Run from anywhere; writes the PDF into the toolbox root.
# Created 2026-10-12 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
# model Claude Opus 5.5 -- WP-124.
set -e
cd "$(dirname "$0")"
pdflatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null
mv main.pdf ../../../lsqvce_theory.pdf
rm -f main.aux main.log main.out main.toc
echo "built lsqvce_theory.pdf"
