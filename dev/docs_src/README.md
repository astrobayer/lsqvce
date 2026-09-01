# docs_src — documentation generators

Everything here is **build-time only**; nothing from this folder ships in a way
that executes, and MATLAB never calls it.

There are now **two independent documentation paths**:

| Path | Builds | Needs |
|---|---|---|
| `build_master_pdf.py` + `master_content_*.py` | `lsqvce_manual.pdf`, **the single master manual** | `reportlab`, Liberation fonts |
| `build_manuals.py` | the browsable HTML manuals | a LaTeX toolchain (equations become PNGs) |

**`build_manuals.py` is currently GUARDED and refuses to run** — it is behind its
own output (see `HANDOVER` §1.3b): the user-manual interface sections exist only
in the shipped HTML, and the source still holds pre-rename method names. Running
it would revert both. Back-port before clearing the guard.

No JavaScript, no Java, no Node anywhere in the doc set.

## Checking the tree

    python3 docs_src/audit_tree.py      # THE GATE: 14 checks, exit 1 on a finding
    python3 docs_src/review_mfiles.py   # a REVIEW of help-text completeness

`audit_tree.py` is a gate and must exit 0. `review_mfiles.py` is deliberately
NOT a gate: it reports what a help block is missing so a human can judge whether
it matters. Both are run from the toolbox root.

## HTML manuals — equations as static PNG images

    python3 docs_src/build_manuals.py

Writes `html/lsqvce_user_manual.html`, `html/lsqvce_reference_manual.html`,
`html/lsqvce_interpretation_guide.html`. Every `\(...\)` / `\[...\]` equation is
rendered from LaTeX to a tightly-cropped transparent PNG under `html/eq/` (cached
by content hash) and embedded as an `<img>`. The pages therefore contain **no
`<script>`, no SVG-as-math, and no external references** — pure static HTML that
renders offline in any browser or MATLAB's Help browser.

Requires `pdflatex`, `pdftoppm` (poppler), and Python `Pillow`.

## PDF manuals — typeset from LaTeX

### The master manual (the one you want)

    python3 docs_src/build_master_pdf.py

Writes `lsqvce_manual.pdf` into the toolbox root. **This is the single
document** and it needs only `reportlab` plus the Liberation TrueType fonts --
no LaTeX toolchain. Full instructions, including how to add an equation or a
diagram, are in `README_manual.md`.

### build_pdfs.py -- RETIRED

Built the three separate PDFs (user manual, reference manual, interpretation
guide) via HTML -> LaTeX -> pdflatex. Those PDFs were replaced by the master
manual on 2026-07-19 and deleted from the toolbox root; the script now refuses
to run. It is kept only as the record of how the old equations were marked up.

## History / rationale

The equations were originally images inside PDFs, then briefly a bundled MathJax
JavaScript runtime (removed as a supply-chain / execution risk in locked-down
environments), and are now: static PNG images in the HTML, and native LaTeX math
in the PDFs. No JavaScript is shipped or required at any point.

Created 2026-07-17 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 4.8.

---

## 3. Per-function pages for functions that have none

    python3 docs_src/build_missing_pages.py

Parses the standardised MATLAB help block of every `+lsqvce` function that has
no `html/lsqvce.<name>.html` page yet and writes one in the house style (same
CSS, nav, see-also and footer as the hand-written pages). It **only writes pages
that do not already exist**, so it can never clobber a curated page — safe to
re-run at any time.

Added 2026-07-18, when an audit found twelve functions (the Ellmer joint
machinery, `contribution`, `streamContribution`, `resolveNc`) shipping without a
page, so `doc lsqvce.<name>` failed for them. After running it, register the new
pages in `+lsqvce/Contents.m`, `helptoc.xml` and `html/index.html` — the script
writes the page but does not edit those three files.

## 4. Tree consistency audit

    python3 docs_src/audit_tree.py        # exit 0 = clean, 1 = findings

Read-only check that the tree is internally consistent: test-file registration,
removed names in executable code, entry-point signatures at every call site,
doc-page and Contents.m coverage, link/helptoc integrity, the JavaScript-free
requirement, and function/end balance. Run it after any rename, after adding a
function, and before shipping. See HANDOVER.md §4e-bis.

## theory_r5/ (REVIEW R5, 2026-09-01)

LaTeX source of `dev/lsqvce_theory_r5.pdf`: theory, practical considerations,
implementation, error handling (Parts I–III) and the R5 review + v5 redesign
plan (Parts IV–V). Build: `dev/docs_src/theory_r5/build.sh` (pdflatex, two
passes). Supersedes `theory_content.py` / `build_theory_pdf.py` for the theory
PDF once decision R5-D6 (HANDOVER §12.6) is taken; until then both exist.
