#!/usr/bin/env python3
"""Reusable ReportLab drawing helpers for the lsqvce master manual.

Two families, both returning platypus-flowable Drawings so they page-break
correctly inside the document:

  eq(...)        a DISPLAY EQUATION drawn glyph-by-glyph. Equations are drawn
                 rather than typeset as text because ReportLab's text engine has
                 no math layout: superscripts, fraction bars and vector overbars
                 have to be POSITIONED, not marked up. Drawing them gives exact
                 control over baseline shifts and lets Symbol-font Greek sit on
                 the same baseline as Times italic variables.

  workflow(...)  a boxes-and-arrows diagram from a compact spec, so every entry
                 point (adjust / adjustBlock / adjustVCE / ...) is drawn by ONE
                 routine and they cannot drift apart stylistically.

WHY NOT LaTeX OR HTML: both were tried for this manual and rejected -- LaTeX
needed a pdflatex toolchain at build time and still produced poor layout for the
mixed prose/diagram material; HTML-to-PDF lost the equations to raster images.
This path has no external dependency beyond reportlab itself.

Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 4.8.
"""
from reportlab.graphics.shapes import (Drawing, String, Line, Rect, Polygon,
                                       Group, PolyLine)
from reportlab.lib.colors import HexColor, white
from reportlab.pdfbase import pdfmetrics

# ---- house palette, shared with the HTML doc set ---------------------------
ORANGE = HexColor('#d9541e')     # headings, accents
BLUE   = HexColor('#4a90d9')     # examples, secondary boxes
GREY   = HexColor('#666666')
DARK   = HexColor('#000000')
LIGHT  = HexColor('#f3f3f3')
BORDER = HexColor('#d3d3d3')
GREEN  = HexColor('#2e7d32')
RED    = HexColor('#c0392b')

SERIF   = 'Times-Roman'
ITALIC  = 'Times-Italic'
BOLD    = 'Times-Bold'
SYMBOL  = 'Symbol'
MONO    = 'Courier'
SANS    = 'Helvetica'
SANSB   = 'Helvetica-Bold'

# ---- Unicode fonts for the drawn equations ---------------------------------
# WHY NOT THE BUILT-IN Symbol FONT. ReportLab's GRAPHICS layer (shapes.String
# inside a Drawing, which is how every equation here is built) does NOT honour
# the Symbol font's custom encoding -- it emits .notdef for every glyph. Proved
# in isolation: String(..., fontName='Symbol') in a Drawing extracts as
# BLACK BOXES, while the identical call on a canvas renders correctly. The
# canvas path is not available to us because equations must be flowables so they
# page-break with the text.
#
# So the equations use a real UNICODE TrueType face instead, and every glyph is
# written as its true codepoint. Liberation Serif is metric-compatible with
# Times and carries Greek, the true minus U+2212, the radical, the Kronecker
# product and the ceiling brackets.
# FreeSerif, not Liberation. Liberation is metric-compatible with Times but its
# MATH coverage is incomplete: it has Greek and the true minus, but NOT the
# Kronecker product U+2297. A glyph the font lacks does not draw a box -- it
# draws NOTHING, and takes the rest of the string with it, so eq 2.12 silently
# lost both the operator and the I after it. That is worse than a box: a box is
# visible. `_assert_coverage` below now refuses to start if any glyph the
# equations use is missing.
_TTF_DIR = '/usr/share/fonts/truetype/freefont/'
_TTF = [('EqSerif',  'FreeSerif.ttf'),
        ('EqItalic', 'FreeSerifItalic.ttf'),
        ('EqBold',   'FreeSerifBold.ttf')]

EQ_SERIF, EQ_ITALIC, EQ_BOLD = SERIF, ITALIC, BOLD      # fallbacks
_TTF_OK = False


def _register_fonts():
    """Register the Unicode faces used by eq(); fall back to base-14 if absent.

    A missing font must not break the build, but it MUST be visible: the
    fallback cannot draw Greek, so eq() checks _TTF_OK and refuses rather than
    silently emitting the black boxes this whole mechanism exists to avoid.
    """
    global EQ_SERIF, EQ_ITALIC, EQ_BOLD, _TTF_OK
    import os
    from reportlab.pdfbase.ttfonts import TTFont
    try:
        for name, fn in _TTF:
            path = os.path.join(_TTF_DIR, fn)
            if not os.path.exists(path):
                return False
            pdfmetrics.registerFont(TTFont(name, path))
        EQ_SERIF, EQ_ITALIC, EQ_BOLD = 'EqSerif', 'EqItalic', 'EqBold'
        _TTF_OK = True
        return True
    except Exception:
        return False


_register_fonts()


def _w(text, font, size):
    """Width of TEXT in FONT at SIZE, in points."""
    return pdfmetrics.stringWidth(text, font, size)


# ===========================================================================
#  EQUATIONS
# ===========================================================================
# An equation is a list of TOKENS. Each token is (kind, text) where kind is:
#   'v'   variable      -> Times-Italic      (x, A, sigma via 'g')
#   'r'   roman/upright -> Times-Roman       (operators, numbers, 'diag')
#   'g'   greek/symbol  -> Symbol font       ('s' renders as sigma, see GREEK)
#   'b'   bold          -> Times-Bold        (matrices/vectors, when needed)
#   'sup' superscript   -> raised, smaller
#   'sub' subscript     -> lowered, smaller
#   'sp'  thin space
# This is deliberately NOT a LaTeX parser: the manual has ~20 display equations
# and a token list is auditable, whereas a half-implemented parser is a source of
# silent mis-rendering.

GREEK = {
    'alpha': '\u03b1', 'beta': '\u03b2', 'gamma': '\u03b3', 'delta': '\u03b4',
    'epsilon': '\u03b5', 'zeta': '\u03b6', 'eta': '\u03b7', 'theta': '\u03b8',
    'lambda': '\u03bb', 'mu': '\u03bc', 'nu': '\u03bd', 'xi': '\u03be',
    'pi': '\u03c0', 'rho': '\u03c1', 'sigma': '\u03c3', 'tau': '\u03c4',
    'phi': '\u03c6', 'chi': '\u03c7', 'psi': '\u03c8', 'omega': '\u03c9',
    'kappa': '\u03ba',
    'Delta': '\u0394', 'Sigma': '\u03a3', 'Omega': '\u03a9', 'Phi': '\u03a6',
    'Lambda': '\u039b', 'Gamma': '\u0393',
}
# Mathematical operators, also as true codepoints. Note 'minus' is U+2212, the
# real minus sign, NOT the ASCII hyphen: in running mathematics the hyphen is
# too short and sits at the wrong height.
SYM_RAW = {
    'otimes': '\u2297', 'approx': '\u2248', 'neq': '\u2260',
    'leq': '\u2264', 'geq': '\u2265', 'minus': '\u2212',
    'prop': '\u221d', 'infty': '\u221e', 'partial': '\u2202',
    'sum': '\u2211', 'sqrt': '\u221a', 'cdot': '\u00b7',
    'in': '\u2208', 'to': '\u2192', 'lceil': '\u2308', 'rceil': '\u2309',
    'ldots': '\u2026', 'times': '\u00d7', 'pm': '\u00b1',
}


def _assert_coverage():
    """Fail at import if the equation font lacks a glyph the manual uses.

    A MISSING glyph is more dangerous than a .notdef box: the box is visible and
    countable in the extracted text, whereas a glyph the font does not have
    simply draws NOTHING and silently truncates the rest of the string. That is
    how equation 2.12 lost its Kronecker product AND the identity matrix after
    it, through a build that reported zero boxes.
    """
    if not _TTF_OK:
        return
    try:
        from fontTools.ttLib import TTFont as _FT
        cmap = _FT(_TTF_DIR + _TTF[0][1], fontNumber=0).getBestCmap()
    except Exception:
        return                      # fontTools absent: skip, do not fail a build
    need = set(GREEK.values()) | set(SYM_RAW.values())
    missing = sorted(c for c in need if ord(c) not in cmap)
    if missing:
        raise RuntimeError(
            'equation font %s lacks %d glyph(s) the manual uses: %s. They would '
            'render as NOTHING, not as boxes. Pick a font with fuller maths '
            'coverage.' % (_TTF[0][1], len(missing),
                           ' '.join('U+%04X' % ord(c) for c in missing)))


_assert_coverage()


def _tokens_width(tokens, size):
    """Total advance width of a token list at base SIZE."""
    total = 0.0
    for kind, text in tokens:
        if kind == 'sp':
            total += size * 0.28
        elif kind in ('hat', 'tilde'):
            # An accent is drawn OVER the preceding glyph and advances nothing:
            # writing ('v','x'),('hat','') must occupy exactly the width of x.
            # Trailing the accent as a normal character instead renders "x ^",
            # which is what a naive token list produces and is wrong.
            total += 0.0
        elif kind in ('sup', 'sub'):
            total += _w(text, SERIF, size * 0.72)
        elif kind == 'g':
            total += _w(_greek(text), EQ_SERIF, size)
        elif kind == 'v':
            total += _w(text, EQ_ITALIC, size)
        elif kind == 'b':
            total += _w(text, EQ_BOLD, size)
        else:
            total += _w(text, EQ_SERIF, size)
    return total


def _greek(name):
    """Map a greek NAME (or raw symbol key) to its Symbol-font character."""
    if name in SYM_RAW:
        return SYM_RAW[name]
    return GREEK.get(name, name)


def eq(tokens, width=440, size=12, pad_top=6, pad_bottom=6, number=None):
    """Draw ONE display equation, centred, as a flowable Drawing.

    tokens  list of (kind, text) -- see the table above
    number  optional equation number, drawn flush right as (n)

    Returns a reportlab Drawing sized to the equation's height.
    """
    if not _TTF_OK:
        # Refuse rather than emit black boxes. The base-14 fallback cannot draw
        # Greek in a Drawing, and a manual whose equations are boxes is worse
        # than no manual -- this is the exact failure this module was rewritten
        # to remove, and it must be loud.
        raise RuntimeError(
            'Unicode equation fonts are not registered (looked in %s). '
            'Equations would render as black boxes. Install the Liberation '
            'fonts or point _TTF_DIR at a Unicode serif face.' % _TTF_DIR)
    h = size * 1.9 + pad_top + pad_bottom
    d = Drawing(width, h)
    total = _tokens_width(tokens, size)
    x = (width - total) / 2.0
    y = pad_bottom + size * 0.45
    prev_w = 0.0                 # width of the last glyph, for accent centring

    for kind, text in tokens:
        if kind == 'sp':
            x += size * 0.28
            continue
        if kind in ('hat', 'tilde'):
            # Centre the accent over the glyph just drawn, and do NOT advance.
            # prev_w is the width of that glyph, captured below.
            ch = '\u02c6' if kind == 'hat' else '\u02dc'
            ax = x - prev_w / 2.0
            d.add(String(ax, y + size * 0.52, ch, fontName=SERIF,
                         fontSize=size * 0.95, fillColor=DARK,
                         textAnchor='middle'))
            continue
        if kind == 'sup':
            s = size * 0.72
            d.add(String(x, y + size * 0.45, text, fontName=EQ_SERIF,
                         fontSize=s, fillColor=DARK))
            prev_w = _w(text, EQ_SERIF, s); x += prev_w
        elif kind == 'sub':
            s = size * 0.72
            d.add(String(x, y - size * 0.26, text, fontName=EQ_SERIF,
                         fontSize=s, fillColor=DARK))
            prev_w = _w(text, EQ_SERIF, s); x += prev_w
        elif kind == 'g':
            ch = _greek(text)
            d.add(String(x, y, ch, fontName=EQ_SERIF, fontSize=size,
                         fillColor=DARK))
            prev_w = _w(ch, EQ_SERIF, size); x += prev_w
        elif kind == 'v':
            d.add(String(x, y, text, fontName=EQ_ITALIC, fontSize=size,
                         fillColor=DARK))
            prev_w = _w(text, EQ_ITALIC, size); x += prev_w
        elif kind == 'b':
            d.add(String(x, y, text, fontName=EQ_BOLD, fontSize=size,
                         fillColor=DARK))
            prev_w = _w(text, EQ_BOLD, size); x += prev_w
        else:
            d.add(String(x, y, text, fontName=EQ_SERIF, fontSize=size,
                         fillColor=DARK))
            prev_w = _w(text, EQ_SERIF, size); x += prev_w

    if number is not None:
        lbl = '(%s)' % number
        d.add(String(width - _w(lbl, EQ_SERIF, size) - 2, y, lbl,
                     fontName=EQ_SERIF, fontSize=size, fillColor=GREY))
    return d


# ===========================================================================
#  WORKFLOW DIAGRAMS
# ===========================================================================
# A diagram is a list of ROWS; each row is a list of nodes. A node is a dict:
#   {'t': 'label text', 'k': kind}
# kind: 'io'    rounded  -- what the USER supplies or receives
#       'proc'  plain    -- a toolbox step
#       'dec'   diamond  -- a branch
#       'note'  dashed   -- an aside, not on the main path
# Rows are stacked top to bottom with a vertical arrow between them, so the
# spec stays declarative and every entry point renders identically.

NODE_H = 26
ROW_GAP = 20
FS = 8.0


def _fit(text, maxw, font=SANS, size=FS):
    """Shrink SIZE until TEXT fits MAXW; return the size actually used."""
    s = size
    while s > 5.0 and _w(text, font, s) > maxw - 8:
        s -= 0.25
    return s


def _node(g, x, y, w, label, kind):
    """Draw one node centred at (x + w/2, y + NODE_H/2)."""
    if kind == 'dec':
        cx, cy = x + w / 2.0, y + NODE_H / 2.0
        g.add(Polygon([cx, y + NODE_H, x + w, cy, cx, y, x, cy],
                      fillColor=white, strokeColor=ORANGE, strokeWidth=1))
        col = ORANGE
    elif kind == 'io':
        g.add(Rect(x, y, w, NODE_H, rx=6, ry=6, fillColor=LIGHT,
                   strokeColor=BLUE, strokeWidth=1))
        col = DARK
    elif kind == 'note':
        g.add(Rect(x, y, w, NODE_H, fillColor=white, strokeColor=GREY,
                   strokeWidth=0.6, strokeDashArray=[2, 2]))
        col = GREY
    else:
        g.add(Rect(x, y, w, NODE_H, fillColor=white, strokeColor=DARK,
                   strokeWidth=0.9))
        col = DARK

    # up to two lines of label
    lines = label.split('|')
    fs = min(_fit(l, w) for l in lines)
    if len(lines) == 1:
        g.add(String(x + w / 2.0, y + NODE_H / 2.0 - fs * 0.35, lines[0],
                     fontName=SANS, fontSize=fs, fillColor=col,
                     textAnchor='middle'))
    else:
        g.add(String(x + w / 2.0, y + NODE_H / 2.0 + 1.0, lines[0],
                     fontName=SANS, fontSize=fs, fillColor=col,
                     textAnchor='middle'))
        g.add(String(x + w / 2.0, y + NODE_H / 2.0 - fs - 0.5, lines[1],
                     fontName=SANS, fontSize=fs, fillColor=col,
                     textAnchor='middle'))


def _arrow(g, x0, y0, x1, y1, label=None, colour=DARK):
    """A line from (x0,y0) to (x1,y1) with a solid head at the far end."""
    g.add(Line(x0, y0, x1, y1, strokeColor=colour, strokeWidth=0.9))
    # head
    if abs(x1 - x0) < 0.01:                    # vertical
        dy = -4 if y1 < y0 else 4
        g.add(Polygon([x1, y1, x1 - 3, y1 - dy, x1 + 3, y1 - dy],
                      fillColor=colour, strokeColor=colour))
    else:                                      # horizontal
        dx = -4 if x1 < x0 else 4
        g.add(Polygon([x1, y1, x1 - dx, y1 - 3, x1 - dx, y1 + 3],
                      fillColor=colour, strokeColor=colour))
    if label:
        mx, my = (x0 + x1) / 2.0, (y0 + y1) / 2.0
        g.add(String(mx + 4, my - 2, label, fontName=SANS, fontSize=6.4,
                     fillColor=GREY))


def workflow(rows, width=460, title=None, loop=None):
    """Draw a stacked boxes-and-arrows workflow.

    rows   list of rows; each row a list of {'t': label, 'k': kind} dicts.
           A '|' in the label splits it over two lines.
    loop   optional (from_row, to_row, 'label') drawn as a feedback arrow up
           the right-hand side -- this is what makes the ITERATIVE engines
           (adjustVCE) legible as loops rather than straight pipelines.

    Returns a Drawing.
    """
    n = len(rows)
    title_h = 14 if title else 0
    h = n * NODE_H + (n - 1) * ROW_GAP + title_h + 8
    d = Drawing(width, h)
    g = Group()

    if title:
        g.add(String(0, h - 10, title, fontName=SANSB, fontSize=8.6,
                     fillColor=ORANGE))

    inner = width - (34 if loop else 0)        # leave a gutter for the loop
    centres = []
    for i, row in enumerate(rows):
        y = h - title_h - (i + 1) * NODE_H - i * ROW_GAP
        k = len(row)
        gap = 10
        wid = (inner - (k - 1) * gap) / float(k)
        row_c = []
        for j, node in enumerate(row):
            x = j * (wid + gap)
            _node(g, x, y, wid, node['t'], node.get('k', 'proc'))
            row_c.append((x + wid / 2.0, y, y + NODE_H))
        centres.append(row_c)

    # vertical arrows between consecutive rows
    for i in range(n - 1):
        top, bot = centres[i], centres[i + 1]
        if len(top) == len(bot):
            pairs = zip(top, bot)
        else:
            # fan in/out from the row's midpoint
            if len(top) == 1:
                pairs = [(top[0], b) for b in bot]
            elif len(bot) == 1:
                pairs = [(t, bot[0]) for t in top]
            else:
                pairs = zip(top, bot)
        for t, b in pairs:
            _arrow(g, t[0], t[1], b[0], b[2] + 1)

    # feedback loop up the right gutter
    if loop:
        r_from, r_to, lbl = loop
        yf = centres[r_from][0][1] + NODE_H / 2.0
        yt = centres[r_to][0][1] + NODE_H / 2.0
        xr = inner + 16
        g.add(PolyLine([inner + 1, yf, xr, yf, xr, yt],
                       strokeColor=ORANGE, strokeWidth=0.9))
        _arrow(g, xr, yt, inner + 1, yt, colour=ORANGE)
        g.add(String(xr + 3, (yf + yt) / 2.0, lbl, fontName=SANS,
                     fontSize=6.4, fillColor=ORANGE))

    d.add(g)
    return d


def decision_tree(spec, width=460, title=None):
    """Draw a decision tree where EDGES ARE DECLARED, not inferred.

    workflow() stacks rows and connects them positionally, which is right for a
    pipeline and WRONG for a branching choice: with two nodes above and two
    below it zips them in order, so a 'no' branch silently acquires an arrow to
    whatever happens to sit beneath it. A decision diagram has to name its edges.

    spec = {
      'nodes': {key: {'t': label, 'k': kind, 'row': r, 'col': c, 'span': s}},
      'edges': [(from_key, to_key, 'edge label'), ...],
    }
    Columns are a 12-unit grid; 'col' is the left edge and 'span' the width, so
    a node can be centred (col 3, span 6) or sit in a half (col 0, span 6).

    Returns a Drawing.
    """
    nodes = spec['nodes']
    edges = spec['edges']
    nrow = max(n['row'] for n in nodes.values()) + 1
    title_h = 14 if title else 0
    h = nrow * NODE_H + (nrow - 1) * ROW_GAP + title_h + 8
    d = Drawing(width, h)
    g = Group()

    if title:
        g.add(String(0, h - 10, title, fontName=SANSB, fontSize=8.6,
                     fillColor=ORANGE))

    # Reserve gutters so a multi-row edge routed around the outside stays
    # inside the drawing; without this it is clipped at the frame.
    GUT = 13.0
    unit = (width - 2 * GUT) / 12.0
    box = {}
    for key, n in nodes.items():
        x = GUT + n['col'] * unit
        w = n.get('span', 4) * unit
        y = h - title_h - (n['row'] + 1) * NODE_H - n['row'] * ROW_GAP
        _node(g, x, y, w - 6, n['t'], n.get('k', 'proc'))
        box[key] = (x, y, w - 6)

    # Rows a horizontal elbow must not cross. An edge spanning MORE than one
    # row cannot bend at the midpoint: that midpoint lands on an intervening
    # row, and the elbow is drawn straight through whatever box sits there.
    # (adjustVCE -> Results in figure 1.1 did exactly this, cutting through
    # adjustBlock.) Such edges are routed around the OUTSIDE instead.
    row_of = {k: n['row'] for k, n in nodes.items()}
    left_edge = GUT + min(n['col'] for n in nodes.values()) * unit
    right_edge = GUT + max((n['col'] + n.get('span', 4)) for n in nodes.values()) * unit

    for e in edges:
        a, b, lbl = (e + (None,))[:3] if len(e) == 2 else e
        ax, ay, aw = box[a]
        bx, by, bw = box[b]
        acx, bcx = ax + aw / 2.0, bx + bw / 2.0
        span = row_of[b] - row_of[a]

        if abs(acx - bcx) < 2.0:
            # straight down, nothing in the way
            _arrow(g, acx, ay, bcx, by + NODE_H + 1, lbl)
        elif span <= 1:
            # adjacent rows: a simple elbow at the midpoint is safe
            midy = (ay + by + NODE_H) / 2.0
            g.add(Line(acx, ay, acx, midy, strokeColor=DARK, strokeWidth=0.9))
            g.add(Line(acx, midy, bcx, midy, strokeColor=DARK, strokeWidth=0.9))
            _arrow(g, bcx, midy, bcx, by + NODE_H + 1)
            if lbl:
                lx = (acx + bcx) / 2.0
                g.add(String(lx, midy + 2.6, lbl, fontName=SANS, fontSize=6.4,
                             fillColor=GREY, textAnchor='middle'))
        else:
            # spans intervening rows: leave the source sideways, run down the
            # nearer outside gutter, then come back in above the target.
            go_right = acx > (left_edge + right_edge) / 2.0
            gut = (right_edge + 5) if go_right else (left_edge - 5)
            sx = (ax + aw) if go_right else ax          # leave from the side
            scy = ay + NODE_H / 2.0
            ty = by + NODE_H + 1
            g.add(PolyLine([sx, scy, gut, scy, gut, ty + 6, bcx, ty + 6],
                           strokeColor=DARK, strokeWidth=0.9))
            _arrow(g, bcx, ty + 6, bcx, ty)
            if lbl:
                g.add(String(gut + (3 if go_right else -3), (scy + ty) / 2.0,
                             lbl, fontName=SANS, fontSize=6.4, fillColor=GREY,
                             textAnchor='start' if go_right else 'end'))
    d.add(g)
    return d


def rule(width=460, colour=ORANGE, thickness=1.0, space=4):
    """A horizontal rule as a flowable."""
    d = Drawing(width, thickness + space * 2)
    d.add(Line(0, space, width, space, strokeColor=colour,
               strokeWidth=thickness))
    return d
