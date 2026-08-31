#!/usr/bin/env python3
"""Master manual, Parts 0-3: front matter, theory, the entry points.

Content module for build_master_pdf.py. Split across master_content_a/b/c so no
single file is unreviewable; the builder concatenates them in order.

Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 4.8.
"""
from reportlab.platypus import Spacer, PageBreak, KeepTogether
from reportlab.lib.units import mm

from build_master_pdf import (P, B, N, code, callout, table, caption, S,
                              BODY_W, VERSION, BUILT)
from pdf_graphics import eq, workflow, decision_tree, rule


def front_matter():
    """Title page and table of contents."""
    st = []
    st += [Spacer(1, 34 * mm)]
    st += [P('lsqvce', 'title')]
    st += [P('Least-squares adjustment, stochastic error modelling '
             'and variance component estimation', 'subtitle')]
    st += [rule(BODY_W, thickness=1.6)]
    st += [Spacer(1, 8)]
    st += [P('This is the single manual for the toolbox. It replaces the '
             'former user manual, reference manual and interpretation guide, '
             'which are consolidated here.', 'body')]
    st += [Spacer(1, 6)]
    st += [P('The toolbox solves large least-squares problems in batches, '
             'estimates the noise model from the residuals rather than '
             'assuming it, and re-weights observation groups by variance '
             'component estimation. It is built for the case where the data '
             'do not fit in memory, the noise is coloured and unknown, and '
             'the relative weighting of the observation types is exactly '
             'what you are trying to find out.', 'body')]
    st += [Spacer(1, 14)]
    st += [table([
        ['Version', VERSION],
        ['Built', BUILT],
        ['Test suite', '497 / 498 passing, 0 failed'],
        ['Entry points', 'adjust, adjustBlock, adjustVCE, snoop, diagnose'],
        ['Start your own case', 'template/lsqvce_template.m'],
    ], widths=[42 * mm, BODY_W - 42 * mm], header=False)]
    st += [Spacer(1, 20)]
    st += [P('Generated with assistance from Claude (Anthropic), model Claude '
             'Opus 4.8. Equations are drawn directly; workflow diagrams are '
             'generated from a single routine so every entry point is drawn '
             'the same way.', 'caption')]
    st += [PageBreak()]

    # ---- contents ---------------------------------------------------------
    st += [P('Contents', 'h1')]
    toc = [
        ('1', 'Orientation: which engine, and why', [
            '1.1  The three front doors',
            '1.2  Choosing between them']),
        ('2', 'Theory', [
            '2.0  Where the theory lives: lsqvce_theory.pdf (R4-D3)']),
        ('3', 'The entry points, step by step', [
            '3.1  adjust -- flat, in memory',
            '3.2  adjustBlock -- block-arrow, one solve',
            '3.3  adjustVCE -- the iterative loop',
            '3.4  snoop -- outlier rejection',
            '3.5  Reading a Results object']),
        ('4', 'Adapting the template: a step-by-step guide', [
            '4.1  What you write, and what the toolbox derives',
            '4.2  Step 1 -- declare the processing options',
            '4.3  Step 2 -- describe the observation groups',
            '4.4  Step 3 -- write buildDay',
            '4.5  Step 4 -- the a-priori covariance',
            '4.6  Step 5 -- epochCov, for multi-axis streams',
            '4.7  Step 6 -- run it, and read the result',
            '4.8  Step 7 -- optional plot hooks',
            '4.9  Four ways to get it wrong']),
        ('5', 'Interpreting the result', [
            '5.1  The order to read things',
            '5.2  Two ways least squares lies to you',
            '5.3  The quadrant you must inspect',
            '5.4  Deciding what to do next',
            '5.5  Reading a gravity solution']),
        ('6', 'Reference', [
            '6.1  Installation',
            '6.2  Symbols',
            '6.3  Equation summary',
            '6.4  Error identifiers',
            '6.5  Testing and the suite']),
        ('R', 'References', []),
    ]
    for num, title, subs in toc:
        st += [P('%s&nbsp;&nbsp;%s' % (num, title), 'toc1')]
        for s in subs:
            st += [P(s.replace('  ', '&nbsp;&nbsp;'), 'toc2')]
    st += [PageBreak()]
    return st


def part1_orientation():
    st = []
    st += [P('1&nbsp;&nbsp;Orientation: which engine, and why', 'h1')]
    st += [P('Three functions solve; everything else diagnoses, models or '
             'accumulates. Picking the right one is the first decision, and it '
             'follows from two properties of your problem: whether the normal '
             'matrix has block-arrow structure, and whether you already know '
             'the noise model.', 'body')]

    st += [callout('Since 2.0 (R4, 2026-09-01) the front door is ONE verb: '
                   '<b>lsqvce.fit(data, opts)</b> -- see the Getting started page '
                   '(doc lsqvce) and DATA_CONTRACT.md. The three engines below are '
                   'what fit dispatches to (lsqvce.engine.adjust / adjustBlock / '
                   'adjustVCE); this part describes them for experts and for '
                   'readers of 1.x code. Every name in this manual that is not '
                   'fit, Options, Results or report is now qualified with '
                   'lsqvce.engine. or lsqvce.kernel. (RELEASE_NOTES.md).', 'note')]
    st += [P('1.1&nbsp;&nbsp;The three front doors', 'h2')]
    st += [table([
        ['Function', 'Structure', 'Noise model', 'Data'],
        ['<b>adjust</b>',
         'none assumed; one dense normal matrix',
         '<b>required</b> -- you supply the covariance',
         'all in memory at once'],
        ['<b>adjustBlock</b>',
         'block-arrow: many local blocks, one common block',
         '<b>required</b> -- you supply it',
         'accumulated in batches'],
        ['<b>adjustVCE</b>',
         'block-arrow, as above',
         '<b>estimated from the data</b>',
         'accumulated, several passes'],
    ], widths=[26 * mm, 52 * mm, 44 * mm, BODY_W - 122 * mm])]
    st += [caption('Table 1.1 -- The three solving entry points. adjustVCE is a '
                   'fluent builder over adjustBlock: it runs the same solve '
                   'repeatedly, re-estimating the noise model and the group '
                   'weights between passes.')]

    st += [P('1.2&nbsp;&nbsp;Choosing between them', 'h2')]
    st += [P('Two questions decide it, in order. First: <b>is the noise model '
             'known?</b> If it is not, it has to be estimated from the '
             'residuals, which is iterative, and <b>adjustVCE</b> is the only '
             'engine that does it. Nothing else about the problem matters yet.',
             'body')]
    st += [P('If the model <i>is</i> known, the second question is size: does '
             'everything -- design matrix, observations and the full covariance '
             '-- fit in memory at once? If it does, <b>adjust</b> takes it '
             'directly. If it does not, the normal equations must be '
             'accumulated batch by batch, which is <b>adjustBlock</b>.', 'body')]
    st += [decision_tree({
        'nodes': {
            'start': {'t': 'your problem', 'k': 'io', 'row': 0, 'col': 3,
                      'span': 6},
            'q1':    {'t': 'is the noise model known?', 'k': 'dec', 'row': 1,
                      'col': 2, 'span': 8},
            'q2':    {'t': 'does everything fit in memory?', 'k': 'dec',
                      'row': 2, 'col': 0, 'span': 7},
            'vce':   {'t': 'adjustVCE|estimate it, then solve', 'k': 'proc',
                      'row': 2, 'col': 7, 'span': 5},
            'adj':   {'t': 'adjust|hand it the covariance', 'k': 'proc',
                      'row': 3, 'col': 0, 'span': 4},
            'blk':   {'t': 'adjustBlock|accumulate in batches', 'k': 'proc',
                      'row': 3, 'col': 4, 'span': 4},
            'res':   {'t': 'Results object', 'k': 'io', 'row': 4, 'col': 3,
                      'span': 6},
        },
        'edges': [
            ('start', 'q1'),
            ('q1', 'q2', 'yes'),
            ('q1', 'vce', 'no'),
            ('q2', 'adj', 'yes'),
            ('q2', 'blk', 'no'),
            ('adj', 'res'),
            ('blk', 'res'),
            ('vce', 'res'),
        ],
    }, width=BODY_W, title='Figure 1.1 -- choosing an engine')]
    st += [caption('Diamonds are decisions, and the branches are labelled. A '
                   'known noise model is a PREREQUISITE for both adjust and '
                   'adjustBlock: each takes the covariance as given. What '
                   'separates those two is only size. adjustVCE is the answer '
                   'whenever the model is not known, regardless of size, '
                   'because estimating it is inherently iterative.')]

    st += [callout(
        '<b>If you are unsure, start with adjustBlock and a plain white '
        'a-priori covariance.</b> Run it, read the report, and look at whether '
        'the residuals are white. If they are not, that is exactly the '
        'condition adjustVCE exists to fix, and you can switch by adding two '
        'fluent calls to the same loaders.', 'note')]

    st += [PageBreak()]
    return st


def part2_theory():
    """R4-D3 (2026-09-01, user decision): the theory lives in lsqvce_theory.pdf
    -- ONE authority per statement (HANDOVER §2.8). This part only says where
    it is and what it covers; the seven former sections were retired, not
    duplicated."""
    st = []
    st += [P('2&nbsp;&nbsp;Theory', 'h1')]
    st += [P('The theory of the toolbox -- the Gauss-Markov model and its two '
             'backends, conditioning and regularisation, variance component '
             'estimation (Helmert, MINQUE, LS-VCE, Monte-Carlo), the separable '
             'stochastic model and the three iteration loops, robust estimation, '
             'the spherical-harmonic setting, error handling, the workflow and '
             'the method\'s limitations -- is documented ONCE, in '
             '<b>lsqvce_theory.pdf</b> at the toolbox root (rebuild: '
             '<i>python3 docs_src/build_theory_pdf.py</i>). Every claim there '
             'names the work package that pinned it in HANDOVER.md.', 'body')]
    st += [callout('Decision R4-D3 (2026-09-01): this manual no longer carries a '
                   'theory chapter of its own; the former sections 2.1-2.7 were '
                   'retired so that no statement exists in two places. Where the '
                   'practitioner chapters cite an equation, the number refers to '
                   'the theory document (1.1-7.1).', 'note')]
    st += [PageBreak()]
    return st
