# docs_src — documentation sources and generators

Everything here is **build-time only**; MATLAB never calls it. Run every command
from the toolbox root.

## The three user documents

| Document | Source | Builder | Output |
|---|---|---|---|
| **Manual** — the step-by-step guide: install, workflow, declaring a problem, the stochastic model, priors, reading the report, the examples, adapting the templates | `manual/*.md` | `python3 dev/docs_src/build_manual.py` | `html/lsqvce_manual.html` (`doc lsqvce`), `lsqvce_manual.pdf` |
| **API reference** — the `data` contract, `buildDay`, covariance specs, `Options`, `Results`, plus generated chapters (user-tier help blocks, the engine / kernel / util units, the identifier catalogue) | `api/*.md` | the same `build_manual.py` | `html/lsqvce_api.html`, `lsqvce_api.pdf` |
| **Theory** — the mathematical framework only | `theory/*.tex` | `sh dev/docs_src/theory/build.sh` (pdflatex) | `lsqvce_theory.pdf` |

`build_manual.py` needs pandoc (HTML, and PDF through xelatex); the theory needs
pdflatex. The rule for what goes where: the Manual says **how to use** the toolbox,
the API reference **what each interface accepts and returns**, the theory **what is
computed and why**. Version history lives in `RELEASE_NOTES.md`, development
decisions in `dev/HANDOVER.md`; none of the three documents carries either.

## Per-function pages

    python3 dev/docs_src/build_missing_pages.py all

Regenerates `html/lsqvce.<name>.html` from the help block of every public unit
(check 25 of the audit compares the pages with a regeneration).

## The audit

    python3 dev/docs_src/audit_tree.py       # exit 0 = clean, 1 = findings
    python3 dev/docs_src/audit_selftest.py   # fault injection: every check fires on a broken copy

Check 23 rebuilds the Manual and the API reference and compares both HTML files
byte for byte; check 29 holds both PDF titles to the version in `+lsqvce/Contents.m`.

## The archived R5 record

The LaTeX of the September 2026 review R5 (theory, the review of v4.0, the v5
redesign plan) is frozen under `dev/reviews/R5_2026-09-01/theory_r5/`. It is
history, not documentation; the current theory is `theory/`.

## The other scripts

`refactor_tiers.py`, `refactor_fluent.py`, `convert_banners.py`,
`flatten_tests.py` and `review_mfiles.py` are one-off or review helpers kept as
the record of a mechanical change; `error_ids.csv`, `tiers.csv` and the
`*_baseline.csv` / `long_*.csv` files are the audit's registries.
