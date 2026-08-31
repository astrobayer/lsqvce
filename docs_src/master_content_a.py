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
        ['Test suite', 'green on every push (CI); counts in README, audit check 17'],
        ['User surface', 'lsqvce.fit, lsqvce.Options, lsqvce.Results, lsqvce.report'],
        ['Engines', 'lsqvce.engine.adjust / adjustBlock / adjustVCE / adjustVCEjoint'],
        ['Start your own case', 'examples/demo_fit.m and DATA_CONTRACT.md'],
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
        ('0', 'Start here: the 2.0 user surface', [
            '0.1  Three lines',
            '0.2  The three kinds of data',
            '0.3  buildDay in six rules',
            '0.4  Reading the answer']),
        ('1', 'Orientation: which engine, and why (ENGINE tier)', [
            '1.1  The three engines fit dispatches to',
            '1.2  Choosing between them']),
        ('2', 'Theory', [
            '2.0  Where the theory lives: lsqvce_theory.pdf (R4-D3)']),
        ('3', 'The entry points, step by step', [
            '3.1  adjust -- flat, in memory',
            '3.2  adjustBlock -- block-arrow, one solve',
            '3.3  adjustVCE -- the iterative loop',
            '3.4  snoop -- outlier rejection',
            '3.5  Reading a Results object']),
        ('4', 'Adapting the template: a step-by-step guide (ENGINE tier)', [
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


def part0_start_here():
    """Chapter 0 -- the USER tier (R4, 2.0). Derived from DATA_CONTRACT.md and
    html/lsqvce_getting_started.html; those two are the authorities, this
    chapter is the manual's copy of their content (L7, WP-68)."""
    st = []
    st += [P('0&nbsp;&nbsp;Start here: the 2.0 user surface', 'h1')]
    st += [P('Since 2.0 the toolbox has four user-facing names: '
             '<font face="Courier">lsqvce.fit</font> (the one verb), '
             '<font face="Courier">lsqvce.Options</font> (what to do), '
             '<font face="Courier">lsqvce.Results</font> (what came back) and '
             '<font face="Courier">lsqvce.report</font> (the verdict). You '
             'describe the problem in a plain struct, <b>data</b>. Everything '
             'that does the work lives one level down, in '
             '<font face="Courier">lsqvce.engine.*</font> and '
             '<font face="Courier">lsqvce.kernel.*</font>, and is what the rest '
             'of this manual describes for experts.', 'body')]
    st += [P('0.1&nbsp;&nbsp;Three lines', 'h2')]
    st += [code("""opts = lsqvce.Options(Model="vce");      % what to do: "gls" one solve, "vce" the loop
R    = lsqvce.fit(data, opts);           % the one verb: your data struct -> a Results
R.report                                 % the five-level verdict""")]
    st += [P('0.2&nbsp;&nbsp;The three kinds of data', 'h2')]
    st += [P('<font face="Courier">fit</font> tells them apart by the fields '
             'present; every rule is checked before any data is read, and an '
             'unknown field name is an error on purpose.', 'body')]
    st += [table([
        ['You have', 'The data struct', 'Engine fit runs'],
        ['A design matrix in memory, diagonal weights, no arc-wise unknowns',
         'data.A (m x u), data.l (m x 1, observed minus computed), optional data.p (weights)',
         'adjust; Robust="irls" -> robust'],
        ['The same, batched or nonlinear',
         'data.batchFcn, data.numUnknowns, data.numBatches; with data.x0 the callback returns the Jacobian at x',
         'adjust / robust / solveNonlinear'],
        ['Days, arcs, arc-wise unknowns, coloured noise (GRACE-FO)',
         'data.common (nc once), data.days, data.systems(k) with your buildDay, observationTypes, covariance, systemId, localsPerArc, arcs',
         'adjustBlock (gls), adjustVCE (vce, pooled), adjustVCEjoint (ellmer)'],
    ], widths=[110, 220, 120])]
    st += [code("""data.common.commonParams = commonTbl;      % nc stated ONCE (or .numCommon)
data.days = 1:30;
s = struct();
s.buildDay         = @buildDayC;           % the one function you write
s.observationTypes = {'kbr','pos'};
s.covariance.kbr   = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',1e-7);
s.covariance.pos   = struct('type','diag','sigma',0.02);
s.systemId         = 'C';
s.arcs             = @(day, obsType) arcSchedule('C', day);   % ragged: ids per day
s.localsPerArc     = 6;                    % fit numbers the local blocks for you
data.systems(1) = s;                       % add satellite D as data.systems(2)""")]
    st += [P('0.3&nbsp;&nbsp;buildDay in six rules', 'h2')]
    st += [P('<font face="Courier">buildDay(day, obsType, ...)</font> returns, for '
             'one day and one observation type, a struct with A (m x nc), B (m x '
             'nl, may be sparse), l (observed minus computed), t (seconds) and, '
             'when the day carries several arcs, arcIdx / arcId.', 'body')]
    st += [N('<b>Deterministic.</b> Same arguments, identical output -- the data is read twice per pass.', 1)]
    st += [N('<b>Raw, not whitened.</b> The toolbox whitens per arc from your covariance spec; the VCE loop reweights there.', 2)]
    st += [N('<b>l is O minus C</b>, and t is in seconds.', 3)]
    st += [N('<b>Arcs by id when ragged.</b> A resolver in arcs plus arcId in the day; a missing id is a data gap.', 4)]
    st += [N('<b>Local blocks from arcInfo.</b> Fill B(rows, off + (1:k)); never compute an offset yourself.', 5)]
    st += [N('<b>A stream called pos is 3-D by default.</b> A scalar stream with that name needs ComponentsPerEpoch = struct(\'pos\', 1).', 6)]
    st += [P('0.4&nbsp;&nbsp;Reading the answer', 'h2')]
    st += [code("""R.report                 % level-0..4 verdict with actionable findings
R.coefficients()         % estimates, formal errors, 95 % bounds as a table
R.x, R.sigma0            % the estimate and the a-posteriori sigma0
R.history, R.noise       % the VCE trajectory and the converged noise model
R.meta                   % provenance: version, release, UTC time, seeds, the options used""")]
    st += [callout('A sigma0 of NaN is a result, not a bug: the quadratic form '
                   'cancelled catastrophically and the toolbox refuses to invent '
                   'a plausible number (lsqvce_theory.pdf, section 7).', 'warn')]
    st += [P('Where to go next: DATA_CONTRACT.md (the binding specification with a '
             'two-satellite worked example), examples/demo_fit.m (the three lines on '
             'every data kind, copy it), doc lsqvce.Options (every knob), '
             'lsqvce_theory.pdf (the method).', 'body')]
    st += [PageBreak()]
    return st


def part1_orientation():
    st = []
    st += [P('1&nbsp;&nbsp;Orientation: which engine, and why (ENGINE tier)', 'h1')]
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
    st += [P('1.1&nbsp;&nbsp;The three engines fit dispatches to', 'h2')]
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
