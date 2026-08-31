#!/usr/bin/env python3
"""Build lsqvce_theory.pdf -- theory, practice, implementation, error handling.

Reuses the master-manual pipeline (build_master_pdf: styles, drawn equations,
page decoration) so the two documents cannot drift apart stylistically; only
the content module and the output path differ. Run from the toolbox root:

    python3 docs_src/build_theory_pdf.py

Writes lsqvce_theory.pdf into the toolbox root. Build-time only.

VERIFY BY LOOKING (HANDOVER §6.2.2): rasterise a few pages after building --
"no exception" is not a verification.

Created 2026-09-01 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Fable 5 -- R4 Level L3, WP-59.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_master_pdf as bm                      # noqa: E402
from theory_content import theory_story           # noqa: E402


def _decorate(canv, doc):
    """Footer: 'lsqvce <version> theory' -- same layout as the manual's."""
    canv.saveState()
    n = canv.getPageNumber()
    if n > 1:
        canv.setFont('Helvetica', 7.4)
        canv.setFillColor(bm.GREY)
        canv.drawString(bm.MARGIN_L, bm.MARGIN_B - 9, 'lsqvce %s theory' % bm.VERSION)
        canv.drawRightString(bm.PAGE_W - bm.MARGIN_R, bm.MARGIN_B - 9, str(n))
        canv.setStrokeColor(bm.BORDER)
        canv.setLineWidth(0.4)
        canv.line(bm.MARGIN_L, bm.MARGIN_B - 4, bm.PAGE_W - bm.MARGIN_R, bm.MARGIN_B - 4)
    canv.restoreState()


def main():
    bm._decorate = _decorate
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    bm.OUT = os.path.join(root, 'lsqvce_theory.pdf')
    bm.build(theory_story())
    print('wrote', bm.OUT)


if __name__ == '__main__':
    main()
