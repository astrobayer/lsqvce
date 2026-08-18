# Building the lsqvce master manual

`lsqvce_manual.pdf` is the single manual for the toolbox. It replaces the former
user manual, reference manual and interpretation guide.

## Build it

From the **toolbox root** (the folder holding `+lsqvce`):

```bash
python3 docs_src/build_master_pdf.py
```

Writes `lsqvce_manual.pdf` into the toolbox root. Build-time only: nothing here
ships as a runtime dependency of the toolbox, and MATLAB never calls it.

## Requirements

| Need | Why | If missing |
|---|---|---|
| `reportlab` >= 4.0 | the whole document | `pip install reportlab` |
| FreeSerif TrueType fonts | **equations, Greek and maths operators** | see below -- the build *refuses* rather than emitting black boxes |

The fonts are expected at `/usr/share/fonts/truetype/freefont/`
(`fonts-freefont-ttf` on Debian/Ubuntu). **Liberation Serif is NOT sufficient**:
it lacks U+2297, U+2308/9, U+2208 and U+221D, and a glyph the font lacks draws
*nothing* rather than a box, silently truncating the rest of the equation.
`_assert_coverage()` checks the font's cmap at import and refuses to build. To use a different Unicode serif face,
point `_TTF_DIR` in `pdf_graphics.py` at it; it must supply Regular, Italic and
Bold.

## The files

| File | Role |
|---|---|
| `build_master_pdf.py` | entry point: styles, page furniture, `P/B/N/code/callout/table/caption`, and `build()` |
| `pdf_graphics.py` | the drawing layer: `eq()`, `workflow()`, `decision_tree()`, `rule()` |
| `master_content_a.py` | front matter, Part 1 (orientation), Part 2 (theory) |
| `master_content_b.py` | Part 3 (entry points), Part 4 (template guide) |
| `master_content_c.py` | Part 5 (interpretation), Part 6 (reference), references |

Content is split across three modules only so no single file is unreviewable;
`build_master_pdf.py` concatenates them in order.

## Why ReportLab, and not LaTeX or HTML

Both were tried for this manual and rejected. LaTeX needed a `pdflatex`
toolchain at build time and still laid out the mixed prose/diagram material
poorly; HTML-to-PDF lost the equations to raster images. ReportLab has no
external toolchain and draws equations and diagrams as real vector content.

## Two things that will bite you

**1. Equations are DRAWN, and the graphics layer ignores the Symbol font.**

`String(..., fontName='Symbol')` inside a `Drawing` emits `.notdef` for every
glyph -- black boxes -- while the identical call on a *canvas* renders fine.
Equations must be Drawings (they have to be flowables so they page-break with
the text), so the Symbol font is unusable here. Every glyph is therefore written
as its **true Unicode codepoint** and drawn in a registered TrueType face. Both
`pdf_graphics.py` (equations) and `build_master_pdf.py` (body text) register the
Liberation faces; entities like `&#8968;` in a Paragraph fail the same way if
they do not.

`eq()` raises if the faces are not registered. That is deliberate: a silent
fallback to base-14 Times produces a manual whose mathematics is black
rectangles.

**2. Verify by extracting text, never by eyeballing a render.**

A `.notdef` box is a real character in the content stream, so it can be counted:

```bash
python3 - <<'EOF'
import subprocess
t = subprocess.run(['pdftotext', 'lsqvce_manual.pdf', '-'],
                   capture_output=True).stdout.decode('utf-8', 'replace')
print('boxes:', sum(t.count(c) for c in '\u25a0\u25a1\ufffd'))
EOF
```

**This must print 0.** A pixel heuristic on a rasterised page does *not* work:
it cannot distinguish a `.notdef` box from a table border, a rule or a filled
arrowhead, and it silently measures the wrong page as soon as the document
repaginates. That mistake was made, and it reported a clean result on a manual
full of boxes.

## Adding an equation

Equations are token lists, not LaTeX. This is deliberate -- there are about
twenty display equations, and a half-implemented parser mis-renders silently
whereas a token list is auditable.

```python
eq([('b', 'l'), ('sp', ''), ('r', '='), ('sp', ''),
    ('b', 'A'), ('b', 'x'), ('sp', ''), ('r', '+'), ('sp', ''), ('b', 'e')],
   width=BODY_W, number='2.1')
```

| Token | Renders as |
|---|---|
| `('r', txt)` | upright roman -- operators, numbers, words |
| `('v', txt)` | italic -- scalar variables |
| `('b', txt)` | bold -- matrices and vectors |
| `('g', name)` | a Greek letter or operator by NAME, e.g. `('g','sigma')`, `('g','minus')`, `('g','otimes')` -- see `GREEK` and `SYM_RAW` |
| `('sup', txt)` / `('sub', txt)` | raised / lowered, smaller |
| `('hat', '')` / `('tilde', '')` | an accent centred over the PRECEDING glyph, advancing zero width |
| `('sp', '')` | a thin space |

Use `('g', 'minus')` rather than a hyphen: U+2212 is the real minus sign, and a
hyphen is too short and sits at the wrong height in mathematics.

## Adding a diagram

Two routines, and picking the wrong one produces a *plausible but wrong*
picture:

* **`workflow(rows, ...)`** -- a PIPELINE. Rows are stacked and connected
  positionally, with an optional feedback loop up the right gutter
  (`loop=(from_row, to_row, 'label')`) for the iterative engines.

* **`decision_tree(spec, ...)`** -- a BRANCHING CHOICE, where **edges are
  declared explicitly**. Use this whenever a diamond has two outcomes.
  `workflow()` infers edges from position, so with two nodes above and two below
  it zips them in order and a "no" branch silently acquires an arrow to whatever
  happens to sit beneath it. Figure 1.1 was drawn wrongly this way.

```python
decision_tree({
    'nodes': {'q': {'t': 'is it known?', 'k': 'dec', 'row': 0,
                    'col': 2, 'span': 8}, ...},
    'edges': [('q', 'yes_target', 'yes'), ('q', 'no_target', 'no')],
}, width=BODY_W, title='Figure N -- ...')
```

Columns are a 12-unit grid: `col` is the left edge, `span` the width. Node kinds
are `io` (rounded, blue -- what the user supplies or receives), `proc` (plain --
a toolbox step), `dec` (diamond, orange -- a branch) and `note` (dashed grey --
an aside).

Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 4.8.
