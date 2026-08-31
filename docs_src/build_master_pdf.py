#!/usr/bin/env python3
"""Build lsqvce_manual.pdf -- the SINGLE master manual, with ReportLab.

Replaces the three separate manuals (user / reference / interpretation), which
were built via an HTML->LaTeX->pdflatex chain that needed an external toolchain
and still laid out the mixed prose/diagram material poorly. Here every equation
is DRAWN (docs_src/pdf_graphics.eq) and every workflow diagram comes from ONE
routine (pdf_graphics.workflow), so the entry points cannot drift apart
stylistically.

Run from the toolbox root:
    python3 docs_src/build_master_pdf.py
Writes lsqvce_manual.pdf into the toolbox root.

Build-time only: nothing here ships as a runtime dependency of the toolbox.

Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 4.8.
"""
import os
import sys
import re

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.lib import colors
from reportlab.platypus import (BaseDocTemplate, PageTemplate, Frame,
                                Paragraph, Spacer, Table, TableStyle,
                                PageBreak, KeepTogether, Preformatted)
from reportlab.pdfbase import pdfmetrics

from pdf_graphics import (eq, workflow, rule, ORANGE, BLUE, GREY, DARK,
                          LIGHT, BORDER, GREEN, RED)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'lsqvce_manual.pdf')

PAGE_W, PAGE_H = A4
MARGIN_L = 22 * mm
MARGIN_R = 18 * mm
MARGIN_T = 20 * mm
MARGIN_B = 18 * mm
BODY_W = PAGE_W - MARGIN_L - MARGIN_R

VERSION = '2.0'   # R4 L5 (WP-65); Contents.m is the ONE version source, this mirrors it
BUILT = '2026-07-19'

# ===========================================================================
#  STYLES
# ===========================================================================
# Body text uses the SAME Unicode faces as the equations. ReportLab's base-14
# Times has no Greek, no true minus, no ceiling brackets and no subscript
# digits, so an entity like &#8968; in a Paragraph renders as a black box -- the
# same failure the equations had, in different clothing. Registering the
# Liberation faces (metric-compatible with Times) fixes prose and equations
# together and keeps them visually identical.
from reportlab.pdfbase.ttfonts import TTFont as _TTFont
_BODY_SERIF, _BODY_ITALIC, _BODY_BOLD = 'Times-Roman', 'Times-Italic', 'Times-Bold'
try:
    _d = '/usr/share/fonts/truetype/liberation/'
    for _n, _f in [('BodySerif', 'LiberationSerif-Regular.ttf'),
                   ('BodyItalic', 'LiberationSerif-Italic.ttf'),
                   ('BodyBold', 'LiberationSerif-Bold.ttf'),
                   ('BodyBoldItalic', 'LiberationSerif-BoldItalic.ttf')]:
        pdfmetrics.registerFont(_TTFont(_n, _d + _f))
    from reportlab.pdfbase.pdfmetrics import registerFontFamily
    registerFontFamily('BodySerif', normal='BodySerif', bold='BodyBold',
                       italic='BodyItalic', boldItalic='BodyBoldItalic')
    _BODY_SERIF, _BODY_ITALIC, _BODY_BOLD = 'BodySerif', 'BodyItalic', 'BodyBold'
except Exception:
    pass

_ss = getSampleStyleSheet()

S = {}
S['body'] = ParagraphStyle(
    'body', parent=_ss['BodyText'], fontName=_BODY_SERIF, fontSize=9.6,
    leading=13.4, alignment=TA_JUSTIFY, spaceBefore=0, spaceAfter=5)
S['bodyc'] = ParagraphStyle('bodyc', parent=S['body'], spaceAfter=2)
S['h1'] = ParagraphStyle(
    'h1', parent=_ss['Heading1'], fontName='Helvetica-Bold', fontSize=17,
    leading=21, textColor=ORANGE, spaceBefore=0, spaceAfter=10)
S['h2'] = ParagraphStyle(
    'h2', parent=_ss['Heading2'], fontName='Helvetica-Bold', fontSize=12.5,
    leading=16, textColor=ORANGE, spaceBefore=14, spaceAfter=6)
S['h3'] = ParagraphStyle(
    'h3', parent=_ss['Heading3'], fontName='Helvetica-Bold', fontSize=10.2,
    leading=13.5, textColor=DARK, spaceBefore=10, spaceAfter=4)
S['h4'] = ParagraphStyle(
    'h4', parent=_ss['Heading4'], fontName='Helvetica-Oblique', fontSize=9.4,
    leading=12.5, textColor=GREY, spaceBefore=8, spaceAfter=3)
S['code'] = ParagraphStyle(
    'code', fontName='Courier', fontSize=7.9, leading=10.2,
    textColor=DARK, leftIndent=6, spaceBefore=3, spaceAfter=5)
S['bullet'] = ParagraphStyle(
    'bullet', parent=S['body'], leftIndent=13, bulletIndent=4, spaceAfter=2.5)
S['num'] = ParagraphStyle(
    'num', parent=S['body'], leftIndent=16, bulletIndent=4, spaceAfter=2.5)
S['note'] = ParagraphStyle(
    'note', parent=S['body'], leftIndent=8, rightIndent=6, spaceBefore=3,
    spaceAfter=3, fontSize=9.2, leading=12.6)
S['caption'] = ParagraphStyle(
    'caption', parent=S['body'], fontSize=8.2, leading=11, textColor=GREY,
    alignment=TA_LEFT, spaceBefore=1, spaceAfter=8)
S['toc1'] = ParagraphStyle(
    'toc1', parent=S['body'], fontName='Helvetica-Bold', fontSize=9.6,
    spaceBefore=5, spaceAfter=1, alignment=TA_LEFT)
S['toc2'] = ParagraphStyle(
    'toc2', parent=S['body'], fontSize=9.0, leftIndent=12, spaceAfter=0.5,
    alignment=TA_LEFT)
S['title'] = ParagraphStyle(
    'title', parent=S['body'], fontName='Helvetica-Bold', fontSize=25,
    leading=29, textColor=ORANGE, alignment=TA_LEFT, spaceAfter=4)
S['subtitle'] = ParagraphStyle(
    'subtitle', parent=S['body'], fontName='Helvetica', fontSize=12.5,
    leading=16, textColor=GREY, alignment=TA_LEFT, spaceAfter=16)
S['ref'] = ParagraphStyle(
    'ref', parent=S['body'], fontSize=8.6, leading=11.4, leftIndent=17,
    firstLineIndent=-17, spaceAfter=3.4, alignment=TA_LEFT)


def P(txt, style='body'):
    return Paragraph(txt, S[style])


def B(txt):
    """A bulleted body paragraph."""
    return Paragraph(txt, S['bullet'], bulletText='\u2022')


def N(txt, n):
    """A numbered body paragraph."""
    return Paragraph(txt, S['num'], bulletText='%d.' % n)


def code(txt):
    """A monospaced code block on a tinted background."""
    lines = [l for l in txt.strip('\n').split('\n')]
    t = Table([[Preformatted('\n'.join(lines), S['code'])]],
              colWidths=[BODY_W])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, -1), LIGHT),
        ('BOX', (0, 0), (-1, -1), 0.6, BORDER),
        ('LINEBEFORE', (0, 0), (0, -1), 2.4, BLUE),
        ('LEFTPADDING', (0, 0), (-1, -1), 7),
        ('RIGHTPADDING', (0, 0), (-1, -1), 5),
        ('TOPPADDING', (0, 0), (-1, -1), 5),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
    ]))
    return t


def callout(txt, kind='note'):
    """A tinted callout box: note (blue) / warn (amber) / danger (red)."""
    bg, bar = {
        'note':   (colors.HexColor('#eaf3fb'), BLUE),
        'warn':   (colors.HexColor('#fff4e5'), colors.HexColor('#e8912d')),
        'danger': (colors.HexColor('#fdecea'), RED),
        'good':   (colors.HexColor('#eaf7ea'), GREEN),
    }[kind]
    t = Table([[Paragraph(txt, S['note'])]], colWidths=[BODY_W])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, -1), bg),
        ('LINEBEFORE', (0, 0), (0, -1), 3.0, bar),
        ('LEFTPADDING', (0, 0), (-1, -1), 8),
        ('RIGHTPADDING', (0, 0), (-1, -1), 7),
        ('TOPPADDING', (0, 0), (-1, -1), 6),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
    ]))
    return t


def table(rows, widths=None, header=True, size=8.3, align=None):
    """A bordered table; ROWS[0] is the header when HEADER."""
    st = ParagraphStyle('tc', parent=S['body'], fontSize=size,
                        leading=size * 1.32, alignment=TA_LEFT,
                        spaceAfter=0, spaceBefore=0)
    sth = ParagraphStyle('th', parent=st, fontName='Helvetica-Bold')
    data = []
    for i, r in enumerate(rows):
        use = sth if (header and i == 0) else st
        data.append([Paragraph(str(c), use) for c in r])
    if widths is None:
        widths = [BODY_W / len(rows[0])] * len(rows[0])
    t = Table(data, colWidths=widths, repeatRows=1 if header else 0)
    style = [
        ('GRID', (0, 0), (-1, -1), 0.45, BORDER),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LEFTPADDING', (0, 0), (-1, -1), 5),
        ('RIGHTPADDING', (0, 0), (-1, -1), 5),
        ('TOPPADDING', (0, 0), (-1, -1), 3.6),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 3.6),
    ]
    if header:
        style += [('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#eeeeee'))]
    if align:
        style += align
    t.setStyle(TableStyle(style))
    return t


def caption(txt):
    return P(txt, 'caption')


# ===========================================================================
#  PAGE FURNITURE
# ===========================================================================
_SECTION = {'name': ''}


def _decorate(canv, doc):
    canv.saveState()
    n = canv.getPageNumber()
    if n > 1:
        canv.setFont('Helvetica', 7.4)
        canv.setFillColor(GREY)
        canv.drawString(MARGIN_L, MARGIN_B - 9,
                        'lsqvce %s manual' % VERSION)
        canv.drawRightString(PAGE_W - MARGIN_R, MARGIN_B - 9, str(n))
        canv.setStrokeColor(BORDER)
        canv.setLineWidth(0.4)
        canv.line(MARGIN_L, MARGIN_B - 4, PAGE_W - MARGIN_R, MARGIN_B - 4)
    canv.restoreState()


def build(story):
    doc = BaseDocTemplate(
        OUT, pagesize=A4,
        leftMargin=MARGIN_L, rightMargin=MARGIN_R,
        topMargin=MARGIN_T, bottomMargin=MARGIN_B,
        title='lsqvce %s -- Manual' % VERSION,
        author='Generated with assistance from Claude (Anthropic)',
        subject='Least-squares adjustment, stochastic error modelling, '
                'variance component estimation')
    frame = Frame(MARGIN_L, MARGIN_B, BODY_W,
                  PAGE_H - MARGIN_T - MARGIN_B, id='body',
                  leftPadding=0, rightPadding=0,
                  topPadding=0, bottomPadding=0)
    doc.addPageTemplates([PageTemplate(id='main', frames=[frame],
                                       onPage=_decorate)])
    doc.build(story)
    return OUT


def main():
    """Assemble every content module in order and write the PDF."""
    # Imported HERE, not at module top: the content modules import names from
    # this one, so a top-level import would be circular.
    from master_content_a import (front_matter, part1_orientation,
                                  part2_theory)
    from master_content_b import part3_entry_points, part4_template
    story = []
    story += front_matter()
    story += part1_orientation()
    story += part2_theory()
    story += part3_entry_points()
    story += part4_template()
    try:
        from master_content_c import (part5_interpretation, part6_reference,
                                      references)
        story += part5_interpretation()
        story += part6_reference()
        story += references()
    except ImportError:
        pass                       # parts 5-6 not written yet
    n = len(story)          # doc.build() CONSUMES the list, so count first
    path = build(story)
    print('wrote %s (%d flowables)' % (path, n))
    return path


if __name__ == '__main__':
    main()
