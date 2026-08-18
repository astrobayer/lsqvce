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
            '2.1  The Gauss-Markov model and BLUE',
            '2.2  QR versus the normal equations',
            '2.3  Residuals, redundancy, reliability',
            '2.4  Variance component estimation',
            '2.5  Block-arrow systems and pre-elimination',
            '2.6  Coloured noise and whitening',
            '2.7  The separable model']),
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
    st = []
    st += [P('2&nbsp;&nbsp;Theory', 'h1')]
    st += [P('This part fixes notation and states the results the toolbox '
             'implements. It is deliberately compact: each equation is here so '
             'that later parts can refer to it by number, not to derive the '
             'subject from scratch. Standard treatments are Koch [3], '
             'Teunissen [4,5] and Grafarend &amp; Schaffrin [8].', 'body')]

    # --- 2.1 ---------------------------------------------------------------
    st += [P('2.1&nbsp;&nbsp;The Gauss-Markov model and BLUE', 'h2')]
    st += [P('The linear model relates <i>n</i> observations <b>l</b> to '
             '<i>u</i> unknowns <b>x</b> through a design matrix <b>A</b>:',
             'body')]
    st += [eq([('b', 'l'), ('sp', ''), ('r', '='), ('sp', ''),
               ('b', 'A'), ('b', 'x'), ('sp', ''), ('r', '+'), ('sp', ''),
               ('b', 'e'), ('r', ','), ('sp', ''), ('sp', ''),
               ('r', 'E{'), ('b', 'e'), ('r', '} = 0,'), ('sp', ''),
               ('r', 'D{'), ('b', 'e'), ('r', '} = '), ('g', 'sigma'),
               ('sub', '0'), ('sup', '2'), ('b', 'Q')],
              width=BODY_W, number='2.1')]
    st += [P('with <b>Q</b> the cofactor matrix and <b>P</b> = <b>Q</b>'
             '<sup>&minus;1</sup> the weight matrix. The best linear unbiased '
             'estimator is the weighted least-squares solution:', 'body')]
    st += [eq([('v', 'x'), ('hat', ''), ('sp', ''), ('r', '='), ('sp', ''),
               ('r', '('), ('b', 'A'), ('sup', 'T'), ('b', 'P'), ('b', 'A'),
               ('r', ')'), ('sup', '\u22121'), ('sp', ''),
               ('b', 'A'), ('sup', 'T'), ('b', 'P'), ('b', 'l')],
              width=BODY_W, number='2.2')]
    st += [P('and its dispersion, with the estimated variance factor, is',
             'body')]
    st += [eq([('r', 'D{'), ('v', 'x'), ('hat', ''), ('r', '} = '),
               ('g', 'sigma'), ('hat', ''), ('sub', '0'), ('sup', '2'),
               ('sp', ''), ('r', '('), ('b', 'A'), ('sup', 'T'), ('b', 'P'),
               ('b', 'A'), ('r', ')'), ('sup', '\u22121'), ('r', ','),
               ('sp', ''), ('sp', ''),
               ('g', 'sigma'), ('hat', ''), ('sub', '0'), ('sup', '2'),
               ('r', ' = '), ('b', 'v'), ('sup', 'T'), ('b', 'P'), ('b', 'v'),
               ('r', ' / '), ('v', 'r')],
              width=BODY_W, number='2.3')]
    st += [P('where <b>v</b> = <b>A</b><b>x</b>&#770; &minus; <b>l</b> are the '
             'residuals and <i>r</i> = <i>n</i> &minus; <i>u</i> is the total '
             'redundancy.', 'body')]

    st += [callout(
        '<b>Weights, not variances.</b> The toolbox takes <b>P</b>, a weight. '
        'Passing a variance instead silently inverts the weighting: your best '
        'observations are down-weighted and every downstream statistic is '
        'wrong while still looking plausible. It is dimensionally valid, so '
        'nothing can detect it for you.', 'danger')]

    # --- 2.2 ---------------------------------------------------------------
    st += [P('2.2&nbsp;&nbsp;QR versus the normal equations', 'h2')]
    st += [P('Forming <b>N</b> = <b>A</b><sup>T</sup><b>P</b><b>A</b> squares '
             'the condition number: a design matrix with condition '
             '&#954;(<b>A</b>) yields &#954;(<b>N</b>) = &#954;(<b>A</b>)'
             '<sup>2</sup>. Where that matters, factor the whitened design '
             'matrix directly [11]:', 'body')]
    st += [eq([('b', 'P'), ('sup', '1/2'), ('b', 'A'), ('sp', ''),
               ('r', '='), ('sp', ''), ('b', 'Q'), ('b', 'R'), ('r', ','),
               ('sp', ''), ('sp', ''),
               ('b', 'R'), ('v', 'x'), ('hat', ''), ('sp', ''),
               ('r', '='), ('sp', ''), ('b', 'Q'), ('sup', 'T'), ('b', 'P'),
               ('sup', '1/2'), ('b', 'l')],
              width=BODY_W, number='2.4')]
    st += [P('The toolbox accumulates in either form. Streaming makes this '
             'practical: <b>R</b> is updated batch by batch and the full '
             'design matrix never exists at once.', 'body')]

    # --- 2.3 ---------------------------------------------------------------
    st += [P('2.3&nbsp;&nbsp;Residuals, redundancy, reliability', 'h2')]
    st += [P('The residuals are a projection of the errors, not the errors '
             'themselves:', 'body')]
    st += [eq([('b', 'v'), ('sp', ''), ('r', '='), ('sp', ''),
               ('g', 'minus'), ('b', 'R'), ('b', 'e'), ('r', ','), ('sp', ''),
               ('sp', ''),
               ('b', 'R'), ('sp', ''), ('r', '='), ('sp', ''),
               ('b', 'I'), ('sp', ''), ('g', 'minus'), ('sp', ''),
               ('b', 'A'), ('r', '('), ('b', 'A'), ('sup', 'T'), ('b', 'P'),
               ('b', 'A'), ('r', ')'), ('sup', '\u22121'), ('b', 'A'),
               ('sup', 'T'), ('b', 'P')],
              width=BODY_W, number='2.5')]
    st += [P('The diagonal of <b>R</b><b>Q</b> gives the redundancy numbers '
             '<i>r<sub>i</sub></i> &#8712; [0,1], which sum to <i>n</i> &minus; '
             '<i>u</i>. They say what fraction of an error in observation '
             '<i>i</i> becomes visible in its own residual:', 'body')]
    st += [eq([('v', 'r'), ('sub', 'i'), ('sp', ''), ('r', '='), ('sp', ''),
               ('r', '('), ('b', 'R'), ('b', 'Q'), ('r', ')'), ('sub', 'ii'),
               ('r', ','), ('sp', ''), ('sp', ''),
               ('g', 'sigma'), ('sub', 'v'), ('sub', 'i'), ('sp', ''),
               ('r', '='), ('sp', ''),
               ('g', 'sigma'), ('sub', '0'), ('g', 'sqrt'),
               ('v', 'r'), ('sub', 'i'), ('g', 'sigma'), ('sub', 'i')],
              width=BODY_W, number='2.6')]
    st += [callout(
        '<b>A residual reveals only the fraction r<sub>i</sub> of the error '
        'present in observation i.</b> Where the redundancy number is near '
        'zero the observation is barely checked by the others, and a gross '
        'error there produces almost no residual. <b>A small residual is not '
        'evidence of a good observation.</b> This is the single most common '
        'misreading of a least-squares result; section 5.3 is devoted to '
        'reading it correctly.', 'danger')]

    st += [P('Baarda\'s data-snooping statistic [1] normalises a residual by '
             'its own standard deviation, and the minimal detectable bias '
             '(MDB) inverts the test to ask how large an error would have to '
             'be before it is found with probability &#946;<sub>0</sub>:',
             'body')]
    st += [eq([('v', 'w'), ('sub', 'i'), ('sp', ''), ('r', '='), ('sp', ''),
               ('v', 'v'), ('sub', 'i'), ('r', ' / '), ('g', 'sigma'),
               ('sub', 'v'), ('sub', 'i'), ('r', ','), ('sp', ''), ('sp', ''),
               ('r', 'MDB'), ('sub', 'i'), ('sp', ''), ('r', '='), ('sp', ''),
               ('g', 'delta'), ('sub', '0'), ('sp', ''),
               ('g', 'sigma'), ('sub', 'i'), ('r', ' / '), ('g', 'sqrt'),
               ('v', 'r'), ('sub', 'i')],
              width=BODY_W, number='2.7')]
    st += [callout(
        'As <i>r<sub>i</sub></i> &#8594; 0 the MDB grows without bound: an '
        'unchecked observation can carry an arbitrarily large error and still '
        'produce <i>w<sub>i</sub></i> &#8776; 0. Pope\'s tau test [2] replaces '
        'the normal quantile when &#963;<sub>0</sub> is estimated from the same '
        'data, but it goes blind at roughly 10&#37; contamination: the blunders '
        'inflate &#963;<sub>0</sub>, tau divides by &#963;<sub>0</sub>, and the '
        'two hide one another.', 'warn')]

    # --- 2.4 ---------------------------------------------------------------
    st += [P('2.4&nbsp;&nbsp;Variance component estimation', 'h2')]
    st += [P('Partition the observations into <i>g</i> groups whose relative '
             'weighting is unknown. Each group carries its own variance '
             'component &#963;<sub><i>g</i></sub><sup>2</sup>, and the combined '
             'cofactor matrix is', 'body')]
    st += [eq([('b', 'Q'), ('sp', ''), ('r', '='), ('sp', ''),
               ('g', 'sum'), ('sub', 'g'), ('sp', ''),
               ('g', 'sigma'), ('sub', 'g'), ('sup', '2'), ('sp', ''),
               ('b', 'Q'), ('sub', 'g')],
              width=BODY_W, number='2.8')]
    st += [P('Helmert\'s estimator [7] and its modern forms [3,6] solve for '
             'the components from the group residual quadratic forms and the '
             'group redundancies. The toolbox uses the fixed-point iteration',
             'body')]
    st += [eq([('g', 'sigma'), ('hat', ''), ('sub', 'g'), ('sup', '2'),
               ('sp', ''), ('r', '='), ('sp', ''),
               ('b', 'v'), ('sub', 'g'), ('sup', 'T'), ('b', 'P'), ('sub', 'g'),
               ('b', 'v'), ('sub', 'g'), ('r', ' / '), ('v', 'r'), ('sub', 'g'),
               ('r', ','), ('sp', ''), ('sp', ''),
               ('v', 'r'), ('sub', 'g'), ('sp', ''), ('r', '='), ('sp', ''),
               ('v', 'n'), ('sub', 'g'), ('sp', ''), ('g', 'minus'), ('sp', ''),
               ('r', 'tr('), ('b', 'N'), ('sup', '\u22121'), ('b', 'N'),
               ('sub', 'g'), ('r', ')')],
              width=BODY_W, number='2.9')]
    st += [P('which is repeated until the components stop moving. The trace '
             'term is the expensive part; for large systems it is estimated by '
             'Hutchinson probes rather than formed exactly, following Koch '
             '&amp; Kusche [9] and Kusche [10]. Only ratios are identifiable, '
             'so one group is held fixed or the set is normalised.', 'body')]

    st += [callout(
        '<b>Groups need redundancy.</b> A group with <i>r<sub>g</sub></i> '
        'near zero has no information about its own variance and its component '
        'will oscillate or diverge. If the variance factors will not settle, '
        'the usual cause is too fine a grouping -- coarsen it.', 'warn')]

    # --- 2.5 ---------------------------------------------------------------
    st += [P('2.5&nbsp;&nbsp;Block-arrow systems and pre-elimination', 'h2')]
    st += [P('Gravity recovery and similar problems have a characteristic '
             'sparsity: many <i>local</i> parameter blocks (one per arc, each '
             'seen only by its own observations) and one <i>common</i> block '
             '(the gravity field, seen by all). The normal matrix is a block '
             'arrow:', 'body')]
    st += [workflow([
        [{'t': 'day of observations', 'k': 'io'}],
        [{'t': 'arc 1|local block', 'k': 'proc'},
         {'t': 'arc 2|local block', 'k': 'proc'},
         {'t': 'arc k|local block', 'k': 'proc'}],
        [{'t': 'Schur complement: eliminate the locals', 'k': 'proc'}],
        [{'t': 'reduced common system (the field)', 'k': 'io'}],
    ], width=BODY_W, title='Figure 2.1 -- block-arrow pre-elimination')]
    st += [caption('Each arc contributes to the common block only through its '
                   'Schur complement. The local parameters never enter the '
                   'reduced system, so the memory cost is set by the common '
                   'block alone.')]
    st += [P('Writing <b>N</b><sub>cc</sub>, <b>N</b><sub>ll</sub> and '
             '<b>N</b><sub>cl</sub> for the common, local and coupling blocks, '
             'the reduced common system is', 'body')]
    st += [eq([('b', 'N'), ('tilde', ''), ('sub', 'cc'), ('sp', ''),
               ('r', '='), ('sp', ''),
               ('b', 'N'), ('sub', 'cc'), ('sp', ''), ('g', 'minus'),
               ('sp', ''),
               ('b', 'N'), ('sub', 'cl'), ('b', 'N'), ('sub', 'll'),
               ('sup', '\u22121'), ('b', 'N'), ('sub', 'lc')],
              width=BODY_W, number='2.10')]
    st += [callout(
        '<b>A near-singular local block poisons the field.</b> If an arc has '
        'too few observations to determine its own state, '
        '<b>N</b><sub>ll</sub><sup>&minus;1</sup> is unreliable and its '
        'contribution to equation 2.10 is noise. The accumulator flags such an '
        'arc rather than absorbing it silently.', 'danger')]

    # --- 2.6 ---------------------------------------------------------------
    st += [P('2.6&nbsp;&nbsp;Coloured noise and whitening', 'h2')]
    st += [P('Observations along a track are correlated in time. Model the '
             'noise as an autoregressive process of order <i>p</i>, fit its '
             'coefficients from the residual autocovariance by the '
             'Levinson-Durbin recursion [15,16], and apply the resulting '
             'filter to both sides of the model:', 'body')]
    st += [eq([('v', 'e'), ('sub', 't'), ('sp', ''), ('r', '+'), ('sp', ''),
               ('g', 'sum'), ('sub', 'k=1'), ('sup', 'p'), ('sp', ''),
               ('g', 'alpha'), ('sub', 'k'), ('sp', ''),
               ('v', 'e'), ('sub', 't\u2212k'), ('sp', ''), ('r', '='), ('sp', ''),
               ('g', 'eta'), ('sub', 't'), ('r', ','), ('sp', ''), ('sp', ''),
               ('g', 'eta'), ('sub', 't'), ('sp', ''), ('r', '~'), ('sp', ''),
               ('r', 'N(0, '), ('g', 'sigma'), ('sup', '2'), ('r', ')')],
              width=BODY_W, number='2.11')]
    st += [P('Whitening is then a filtering of <b>A</b> and <b>l</b>, after '
             'which ordinary least squares applies. Whether the filter '
             'succeeded is testable: the Ljung-Box statistic [17] on the '
             'filtered residuals should not reject whiteness.', 'body')]

    # --- 2.7 ---------------------------------------------------------------
    st += [P('2.7&nbsp;&nbsp;The separable model', 'h2')]
    st += [P('A multi-axis stream -- three position components per epoch, say '
             '-- has both a spatial covariance between axes and a temporal '
             'correlation along the track. When the spatial part is known '
             'per epoch (from a PPP solution, for instance) the covariance '
             'separates:', 'body')]
    st += [eq([('b', 'C'), ('sp', ''), ('r', '='), ('sp', ''),
               ('g', 'sigma'), ('sup', '2'), ('sp', ''),
               ('b', 'D'), ('sp', ''), ('r', '('), ('b', 'R'), ('sub', 't'),
               ('sp', ''), ('g', 'otimes'), ('sp', ''), ('b', 'I'),
               ('sub', 'd'), ('r', ')'), ('sp', ''), ('b', 'D'),
               ('sup', 'T')],
              width=BODY_W, number='2.12')]
    st += [P('with <b>D</b> the known per-epoch factor, <b>R</b><sub>t</sub> '
             'the temporal correlation the toolbox estimates, and <i>d</i> the '
             'number of axes. Supplying <b>D</b> lets the temporal law be '
             'recovered by exact de-scaling rather than inferred from whitened '
             'residuals. This is the model that gives the toolbox its name and '
             'the reason multi-axis streams need the epochCov hook of section '
             '4.6.', 'body')]

    st += [P('How <b>D</b> is actually stored', 'h3')]
    st += [P('<b>D</b> is never formed as the large sparse matrix that '
             'equation 2.12 suggests. It is block-diagonal by epoch, so '
             'only the blocks are stored, one per epoch, and the toolbox '
             'assembles the effect on the fly. Two layouts are accepted:',
             'body')]
    st += [table([
        ['Layout', 'Shape', 'Meaning', 'Use when'],
        ['Diagonal',
         '<font face="Courier">nEpoch &#215; d</font>',
         'per-axis standard deviations for each epoch; the axes are '
         'uncorrelated within an epoch',
         'the axes are independent, or you only have per-axis sigmas'],
        ['Full',
         '<font face="Courier">nEpoch &#215; d &#215; d</font>',
         'the complete <i>d</i>&#215;<i>d</i> covariance of each epoch, '
         'e.g. straight from a PPP solution',
         'the axes are correlated within an epoch (the realistic 3-D case)'],
    ], widths=[19 * mm, 27 * mm, 62 * mm, BODY_W - 108 * mm])]
    st += [caption('Table 2.1 -- the two accepted layouts of Cep. For 3-D '
                   'positions d = 3, so these are nEpoch&#215;3 and '
                   'nEpoch&#215;3&#215;3. Anything else raises '
                   'lsqvce:covModel:badCep.')]

    st += [callout(
        '<b>Not <font face="Courier">n&#215;9</font>, and not a sparse '
        '<font face="Courier">3n&#215;3n</font>.</b> A flattened '
        '<font face="Courier">n&#215;9</font> array is rejected: the third '
        'dimension carries the axis pairing, and flattening it loses which '
        'entry is which. Building the full sparse <b>D</b> yourself is possible '
        'but wasteful -- it is block-diagonal, and the whitening only ever '
        'needs one epoch\'s block at a time.', 'warn')]

    st += [P('<b>The indexing rule that ties this to row order.</b> Because '
             'rows are epoch-major, row <i>k</i> belongs to epoch '
             '&#8968;<i>k</i>/<i>d</i>&#8969;. The toolbox slices the epochs '
             'of a batch with exactly that expression, so the '
             '<i>e</i>-th slice of <font face="Courier">Cep</font> must '
             'describe the same epoch as rows '
             '<i>d</i>(<i>e</i>&minus;1)+1 &hellip; <i>de</i> of your '
             '<b>A</b>, <b>l</b> and <b>t</b>:', 'body')]
    st += [eq([('r', 'epoch('), ('v', 'k'), ('r', ') = '),
               ('g', 'lceil'), ('v', 'k'), ('r', ' / '), ('v', 'd'),
               ('g', 'rceil'), ('r', ','), ('sp', ''), ('sp', ''),
               ('r', 'rows('), ('v', 'e'), ('r', ') = '),
               ('v', 'd'), ('r', '('), ('v', 'e'), ('g', 'minus'),
               ('r', '1) + 1 '), ('g', 'ldots'), ('r', ' '), ('v', 'd'), ('v', 'e')],
              width=BODY_W, number='2.13')]
    st += [P('An <font face="Courier">epochCov</font> hook therefore returns '
             'both the per-epoch array and, for each arc, which day-epochs and '
             'which rows it covers -- so the two indexings can be checked '
             'against each other rather than assumed to agree. Section 4.6 '
             'gives the concrete return shape.', 'body')]

    st += [callout(
        '<b>Row order is load-bearing here.</b> The Kronecker structure in '
        'equation 2.12 assumes <b>epoch-major</b> rows, '
        '[x&#8321; y&#8321; z&#8321; | x&#8322; y&#8322; z&#8322; | &hellip;]. '
        'Stored axis-major, equation 2.13 pairs every row with the wrong epoch, '
        'the temporal law is stretched by the axis count, and nothing errors. '
        'See section 4.9.', 'danger')]
    st += [PageBreak()]
    return st
