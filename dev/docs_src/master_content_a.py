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
        ['Engines', 'lsqvce.engine.adjustBlock / adjustRobust / adjustNonlinear / adjustVCEjoint (v5: one accumulator)'],
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
        ('0', 'Start here: the 4.0 user surface', [
            '0.1  Three lines',
            '0.2  The three kinds of data',
            '0.3  buildDay in six rules',
            '0.4  Reading the answer']),
        ('D', 'The data contract (v3) -- the binding specification', [
            'D.1  The principle: state nothing the data already states',
            'D.2  The data struct, field by field',
            'D.3  The buildDay contract',
            'D.4  Arc schedules',
            'D.5  Covariance specs',
            'D.6  Local (arc-wise) unknowns',
            'D.7  Multi-axis streams and the epochCovariance hook',
            'D.8  Worked example: two satellites, ragged arcs',
            'D.9  The equivalent engine calls',
            'D.10 Checklist, and the two in-memory kinds']),
        ('1', 'Orientation: which engine, and why (ENGINE tier)', [
            '1.1  The three engines fit dispatches to',
            '1.2  Choosing between them']),
        ('2', 'Theory', [
            '2.0  Where the theory lives: lsqvce_theory.pdf (R4-D3)']),
        ('3', 'The entry points, step by step', [
            '3.1  adjust -- flat, in memory',
            '3.2  adjustBlock -- block-arrow, one solve',
            '3.3  The one loop -- Estimate="noise" (adjustVCEjoint)',
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
    st += [P('0&nbsp;&nbsp;Start here: the 4.0 user surface', 'h1')]
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
    st += [code("""opts = lsqvce.Options(Estimate="noise");  % what to do: "none" one solve, "noise" the loop
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
         'the block adapter (batchSystem); Robust="huber" -> adjustRobust'],
        ['The same, batched or nonlinear',
         'data.batchFcn, data.numBatches, data.numUnknowns (linear form); with data.x0 the callback returns the Jacobian at x and u = numel(x0)',
         'adjust / robust / solveNonlinear'],
        ['Days, arcs, arc-wise unknowns, coloured noise (GRACE-FO)',
         'data.commonParams (scoped template rows), data.days, data.systems(k) with your buildDay and covariance (its fieldnames are the observation types); arcs / scoped localParams as needed',
         'adjustBlock (Estimate none), adjustVCEjoint (Estimate noise: Pool / ArcFactors say which policy)'],
    ], widths=[110, 220, 120])]
    st += [code("""data.commonParams = [repmat({'common'}, nc, 1), commonTbl];  % scoped rows, nc derived
data.days = 1:30;
s = struct();
s.buildDay         = @buildDayC;           % the one function you write
s.covariance.kbr   = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',1e-7);
s.covariance.pos   = struct('type','diag','sigma',0.02);   % fieldnames = the types
s.id               = 'C';
s.arcs             = @(day, obsType) arcSchedule('C', day);   % ragged: ids per day
s.localParams      = [repmat({'arc'}, 6, 1), icLabels];    % 6 per (day, arc); 'day'/
data.systems(1) = s;                       % 'system' rows scale the same way""")]
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
    st += [P('Where to go next: chapter D (the data contract in full; '
             'DATA_CONTRACT.md in the repository is its source), '
             'examples/demo_fit.m (the three lines on '
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
                   'what fit dispatches to (lsqvce.engine.adjustBlock / adjustRobust / '
                   'adjustNonlinear / adjustVCEjoint); this part describes them for experts and for '
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
        ['<b>adjustVCEjoint</b>',
         'block-arrow, as above',
         '<b>estimated from the data</b>',
         'accumulated, several passes'],
    ], widths=[26 * mm, 52 * mm, 44 * mm, BODY_W - 122 * mm])]
    st += [caption('Table 1.1 -- The solving entry points (v5: every kind on the one accumulator). adjustVCEjoint is a '
                   'fluent builder over adjustBlock: it runs the same solve '
                   'repeatedly, re-estimating the noise model and the group '
                   'weights between passes.')]

    st += [P('1.2&nbsp;&nbsp;Choosing between them', 'h2')]
    st += [P('Two questions decide it, in order. First: <b>is the noise model '
             'known?</b> If it is not, it has to be estimated from the '
             'residuals, which is iterative, and <b>adjustVCEjoint</b> is the only '
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
            'vce':   {'t': 'Estimate="noise"|estimate it, then solve', 'k': 'proc',
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
                   'separates those two is only size. Estimate="noise" is the answer '
                   'whenever the model is not known, regardless of size, '
                   'because estimating it is inherently iterative.')]

    st += [callout(
        '<b>If you are unsure, start with adjustBlock and a plain white '
        'a-priori covariance.</b> Run it, read the report, and look at whether '
        'the residuals are white. If they are not, that is exactly the '
        'condition the variance-component loop exists to fix, and you can switch by adding two '
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
             '<i>python3 dev/docs_src/build_theory_pdf.py</i>). Every claim there '
             'names the work package that pinned it in HANDOVER.md.', 'body')]
    st += [callout('Decision R4-D3 (2026-09-01): this manual no longer carries a '
                   'theory chapter of its own; the former sections 2.1-2.7 were '
                   'retired so that no statement exists in two places. Where the '
                   'practitioner chapters cite an equation, the number refers to '
                   'the theory document (1.1-7.1).', 'note')]
    st += [PageBreak()]
    return st


def partD_data_contract():
    """Chapter D -- the data contract, v3 (L10 WP-73). The manual's copy of
    DATA_CONTRACT.md; that file is the binding source, this chapter mirrors it
    so the printed manual is complete without the repository."""
    st = []
    st += [P('D&nbsp;&nbsp;The data contract (v3) -- the binding specification', 'h1')]
    st += [P('This chapter is the manual\'s copy of DATA_CONTRACT.md (version 4, '
             'adopted 2026-08-31, Level L10): what you give '
             '<font face="Courier">lsqvce.fit</font>, and how to write '
             '<font face="Courier">buildDay</font>. The repository file is the '
             'binding source; the two are kept in step by the release process.', 'body')]

    st += [P('D.1&nbsp;&nbsp;The principle: state nothing the data already states', 'h2')]
    st += [P('The contract declares <b>structure</b>, never <b>counts</b>. Every '
             'count is derived from the one place it already lives:', 'body')]
    st += [table([
        ['Quantity', 'Derived from'],
        ['nc (common unknowns)', "the commonParams expansion: 'common' rows once + 'day' rows per calendar day with arcs"],
        ['nl (local unknowns)', "the localParams expansion: 'arc' rows per (day, arc) + 'day' rows per day with arcs + 'system' rows once; not declared -> nl = 0"],
        ['observation types', 'fieldnames(covariance), in declaration order'],
        ['u (batched nonlinear)', 'numel(x0)'],
    ], widths=[130, 320])]
    st += [P('The v2 count fields (numCommon, numLocal, numUnknowns-with-x0) and '
             'the duplicate type list (observationTypes) are gone. One source '
             'each -- nothing left to disagree. A v2 declaration fails loudly: '
             'the retired names are unknown fields now, and '
             '<font face="Courier">lsqvce:fit:unknownField</font> names the field.', 'body')]

    st += [P('D.2&nbsp;&nbsp;The data struct, field by field', 'h2')]
    st += [P('<font face="Courier">data</font> is a scalar struct; field names are '
             'case-sensitive lowerCamelCase, unknown names are errors on purpose.', 'body')]
    st += [table([
        ['Field', 'Type / size', 'Meaning'],
        ['commonParams', 'cell, R x C', 'REQUIRED. One row per COMMON unknown: nc = R. Column count and meaning are your convention; buildDay receives the table as its 3rd argument. Minimal: num2cell((1:nc).\')'],
        ['days', 'double, 1 x D', 'default day list (numeric integers; the cache is keyed on them); required unless every system carries its own'],
        ['systems', 'struct, 1 x M', 'one entry per observing system (below)'],
    ], widths=[80, 80, 290])]
    st += [P('Each <font face="Courier">data.systems(k)</font> -- two required '
             'fields, the rest optional (fill unused optionals with []):', 'body')]
    st += [table([
        ['Field', 'Default', 'Meaning'],
        ['buildDay', 'REQUIRED', 'your day builder (D.3)'],
        ['covariance', 'REQUIRED', 'struct, one field per observation type holding a covariance spec (D.5). The fieldnames ARE the observation types; their order is the processing order'],
        ['arcs', '1', 'which arcs each day carries (D.4); \'fromDay\' is ENGINE-tier only'],
        ['localParams', '{}', "R x C cell of LOCAL templates: column 1 = 'arc' | 'day' | 'system' (D.6), nl derived; buildDay's 4th argument. Not declared -> nl = 0"],
        ['days', 'data.days', 'per-system schedule override'],
        ['id', 'the index k', 'numeric scalar or single char, unique across systems (explicit ids and index defaults share one namespace)'],
        ['name', "''", 'human-readable label for reports'],
        ['processingOptions', 'struct()', 'anything buildDay needs; its 5th argument'],
        ['epochCovariance', '[]', '@(day,obsType) per-epoch D hook for multi-axis streams (D.7)'],
    ], widths=[75, 55, 320])]
    st += [callout('struct() with cell values expands into a struct ARRAY: '
                   's = struct(\'localParams\', {\'bias\',\'drift\'}) is 1 x 2. '
                   'Use double braces, or assign field by field.', 'warn')]

    st += [P('D.3&nbsp;&nbsp;The buildDay contract', 'h2')]
    st += [P('fit calls buildDay once per (day, observationType) and cuts the day '
             'into arcs itself. Every pass reads each day twice (normals, then '
             'residuals); the VCE loop multiplies that by the passes. Hence: '
             '<b>deterministic</b> (same arguments, byte-identical output -- no '
             'randomness, clocks, mutable globals), and cache expensive raw reads '
             'yourself (the toolbox caches one built day).', 'body')]
    st += [P('fit passes as many of these six arguments as the handle declares '
             '(nargin); a varargin builder receives all six:', 'body')]
    st += [code("""D = buildDay(day, obsType)                                  % classic
D = buildDay(day, obsType, common, local, procOpt)          % parameter-driven
D = buildDay(day, obsType, common, local, procOpt, arcInfo) % with local-block plan""")]
    st += [P('Return raw, unwhitened matrices -- whitening happens per arc inside '
             'the toolbox, and the VCE loop reweights there. Required fields: '
             'A (m x nc, common design), B (m x nl, local design, may be sparse), '
             'l (m x 1, OBSERVED MINUS COMPUTED), t (m x 1, seconds). Optional: '
             'arcIdx (cell of row-index vectors; absent = one arc), arcId (the '
             'identity of each arc for ragged schedules), nAxis and Cep for '
             'multi-axis streams (D.7).', 'body')]
    st += [P('Declaration errors are raised by fit BEFORE any data is read '
             '(lsqvce:fit:unknownField / noSystems / noCommon / missingField / '
             'badScope / noDays / badId / duplicateSystem / builderKind / '
             'badArcs, plus every makeLoader construction id relabelled under '
             'lsqvce:fit with the original as the cause). Reading errors keep '
             'their engine ids: lsqvce:makeLoader:buildDay / badDay / dayShape / '
             'badArc / badArcId, lsqvce:accum:arcBadNc / arcBadNl -- the derived '
             'nc and nl are enforced against every day.', 'body')]

    st += [P('D.4&nbsp;&nbsp;Arc schedules', 'h2')]
    st += [P('<b>Static</b> -- every day looks the same: arcs = 4 (ids 1..4), '
             'arcs = [1 2 4], or per type arcs.kbr = 2. Arcs are addressed by '
             'position in arcIdx; an id beyond numel(arcIdx) is an error, which '
             'catches a mis-declared schedule.', 'body')]
    st += [P('<b>Ragged</b> -- days differ (gaps, manoeuvres, short months): '
             'declare a resolver that returns the arc ids of a day without '
             'reading data, and let buildDay return arcId naming the arcs it '
             'actually carries. Arcs are then addressed by id; an id a day does '
             'not carry is a data gap and is skipped -- no padding, no empty arcs.', 'body')]
    st += [code("""s.arcs = @(day, obsType) myArcTable{day};   % 1 x A vector of ids for that day""")]

    st += [P('D.5&nbsp;&nbsp;Covariance specs', 'h2')]
    st += [P('One struct per observation type; type selects the model. sigma is '
             'always a STANDARD DEVIATION, never a variance; the VCE loop scales '
             'the specs at run time, you give the starting point. Full reference: '
             'doc lsqvce.engine.covModel.', 'body')]
    st += [table([
        ['type', 'Fields', 'Use for'],
        ["'diag'", 'sigma (scalar or m x 1 std)', 'uncorrelated noise'],
        ["'arma'", 'ar [1 a1 ... ap], ma (default 1), sigma', 'coloured noise with a known filter'],
        ["'acf'", 'acf @(dt,opt), sigma', 'correlation as a function of the lag'],
        ["'psd'", 'f (nf x 1, Hz), S (nf x 1)', 'noise given as a one-sided PSD'],
        ["'var'", 'ar (d x d x p), sigma (d x d)', 'vector-AR for multi-axis streams'],
    ], widths=[45, 175, 230])]

    st += [P('D.6&nbsp;&nbsp;Local unknowns: the scoped localParams table', 'h2')]
    st += [P('One mechanism: the rows are templates, column 1 says how often each '
             'instantiates. <b>\'arc\'</b>: once per (day, arc) of this system '
             '(ICs, per-arc empiricals). <b>\'day\'</b>: once per day whose '
             'schedule carries an arc (ACC scale/bias/drift). <b>\'system\'</b>: '
             'once over the whole adjustment (instrument constants; also the v3 '
             'free table under migration). nl is derived. Binding layout of B: '
             'per day [arc blocks in sorted-id order | the \'day\' block], then '
             'one trailing \'system\' block; a day without arcs gets no blocks '
             'at all -- offsets stay gap-safe on manoeuvre days. Fill B only '
             'from buildDay\'s 6th argument: info.offsetOf(arcId), '
             'info.dayOffset/.dayCount, info.sysOffset/.sysCount -- and A\'s day '
             'columns from info.commonDayOffset/.commonDayCount; never compute '
             'an offset yourself. <b>Not declared</b>: nl = 0 -- return '
             'B = zeros(m, 0); the engine skips the pre-elimination and the local '
             'diagnostics come back as consistent empty shapes. Heterogeneous '
             'per-arc counts (@(day,arcId)) stay ENGINE-tier on purpose: '
             'lsqvce.engine.makeLoader(..., \'LocalsPerArc\', fh).', 'body')]

    st += [P('D.7&nbsp;&nbsp;Multi-axis streams and the epochCovariance hook', 'h2')]
    st += [P('A kinematic position has d = 3 components per epoch sharing one '
             'temporal law. Return the rows EPOCH-MAJOR (x1 y1 z1 x2 y2 z2 ...), '
             'set nAxis = 3, and provide the per-epoch scaling D. The STANDARD '
             'source is Cep on the returned day (nEpoch x d, or nEpoch x d x d '
             'for full per-epoch blocks); a static scaling can ride on the '
             'covariance spec instead. A multi-axis stream with no D is an error.', 'body')]
    st += [P('The <b>epochCovariance hook</b> is optional and matters only under '
             'Estimate="noise": whitening never needs it, but the per-axis estimation of the '
             'VCE iteration (re-estimating the temporal AR law from residuals) '
             'does not read the day\'s Cep. With the hook it takes the exact '
             'separable de-scale path (lsqvce.kernel.residualAcv with the full '
             'per-epoch D, fetched per day WITHOUT rebuilding the day); without '
             'it, it falls back to descaling by the stream\'s current scalar '
             'sigma -- an approximation that biases the estimated AR law when D '
             'varies strongly across epochs (changing geometry / PDOP). Declare '
             'the hook when you run VCE on multi-axis streams and have D '
             'available independently of the day build.', 'body')]

    st += [P('D.8&nbsp;&nbsp;Worked example: the GRACE parameterization, three scopes, one table', 'h2')]
    st += [P('Arc-wise initial conditions, day-wise accelerometer scale/bias/'
             'drift, per-arc KBR empiricals -- the case v3 could only '
             'hand-number, and the driver of v4:', 'body')]
    st += [code("""function D = buildDayGrace(day, obsType, common, ~, procOpt, info)
    raw = readGraceDay(procOpt.dataDir, day, obsType);   % your reader, deterministic
    D.A = shPartials(raw, common);                       % m x info.nc, 'common' block
    D.l = raw.observed - raw.computed;                   % O - C
    D.t = raw.tSeconds;
    D.arcIdx = raw.arcRows;   D.arcId = raw.arcIds;
    D.B = sparse(numel(D.l), info.nl);
    for a = 1:numel(D.arcId)
        off = info.offsetOf(D.arcId(a));
        D.B(D.arcIdx{a}, off + (1:info.count(a))) = icAndKbrPartials(raw, a);
    end
    if info.dayCount > 0                                 % this day's ACC block
        D.B(:, info.dayOffset + (1:info.dayCount)) = accPartials(raw);
    end
end""")]
    st += [code("""data.commonParams = [repmat({'common'}, 3717, 1), shParameterTable(60)];
data.days = 1:31;
cov.kbr = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',1e-7);
cov.pos = struct('type','diag','sigma',0.02);  % fieldnames = the two types
s = struct();
s.buildDay = @buildDayGrace;   s.covariance = cov;
s.id = 'C';   s.name = 'GRACE-FO C';
s.arcs = @(day, obsType) arcSchedule('C', day);
s.localParams = [ ...
    [repmat({'arc'}, 6, 1), {'x0';'y0';'z0';'vx0';'vy0';'vz0'}];  % ICs per arc
    {'arc', 'kbrBias'; 'arc', 'kbrCos1cpr'};                      % KBR per arc
    [repmat({'day'}, 3, 1), {'accScale';'accBias';'accDrift'}]];  % ACC per day
s.processingOptions = struct('dataDir', '/data/gracefo');
data.systems(1) = s;
s.id = 'D';  s.name = 'GRACE-FO D';  s.arcs = @(day, obsType) arcSchedule('D', day);
data.systems(2) = s;
% nl per satellite = SUM over days-with-arcs of (arcs*8 + 3) -- derived
R = lsqvce.fit(data, Estimate="noise", Pool="obsType", Parallel=struct('Mode', "parfor"));
R.report""")]

    st += [P('D.9&nbsp;&nbsp;The equivalent engine calls', 'h2')]
    st += [P('fit builds one Loader per system through lsqvce.engine.makeLoader '
             '(field names map 1:1 to its Name-Value names, lowerCamelCase to '
             'UpperCamelCase, plus the D.1 derivations) and dispatches on '
             'opts.Estimate / opts.Robust / data.x0 to adjustBlock / adjustRobust / adjustNonlinear / adjustVCEjoint. These are '
             'ENGINE-tier calls; a user needs none of them. The ENGINE tier '
             'keeps the v2/v3 conveniences (numCommon/numLocal counts, free '
             'un-scoped tables via Scoped=false, LocalsPerArc incl. @(day,arcId), '
             '\'fromDay\', cache, plot hooks, raw loader closures) -- they left '
             'the USER surface, not the toolbox.', 'body')]

    st += [P('D.10&nbsp;&nbsp;Checklist, and the two in-memory kinds', 'h2')]
    st += [N('buildDay returns A, B, l, t with one row count; l is O - C; nothing is whitened.', 1)]
    st += [N('Calling buildDay twice with the same arguments gives identical output; t is in seconds.', 2)]
    st += [N("commonParams rows carry the scope token in column 1; size(A, 2) must equal info.nc, blocks placed by the info offsets.", 3)]
    st += [N('Every observation type is a field of covariance; the fieldnames are the type list.', 4)]
    st += [N('Ragged schedule -> resolver + arcId; static -> count/vector, no arcId.', 5)]
    st += [N("Locals: scoped localParams templates, B filled from the info offsets only; not declared -> B = zeros(m, 0).", 6)]
    st += [N('Multi-axis: epoch-major rows, nAxis, Cep on the day; the epochCovariance hook for the exact VCE Shape path.', 7)]
    st += [N('No unknown field names -- a typo is an error, and so are the v2 count fields and, since v4, localsPerArc.', 8)]
    st += [P('Not every problem streams. When A fits in memory and there are no '
             'arc-wise unknowns, data has no systems field at all; fit tells the '
             'kinds apart by the fields present. <b>Kind B</b>: A (m x u), l '
             '(m x 1, O - C), optional p (weights). <b>Kind C</b>: batchFcn with '
             '[A, l, p] = batchFcn(k), numBatches, and numUnknowns on the LINEAR '
             'form; with x0 the problem is nonlinear ([A, w, p] = batchFcn(k, x), '
             'Gauss-Newton) and u = numel(x0) -- numUnknowns must then be absent '
             '(lsqvce:fit:uCount). Model must be "gls" on both kinds; '
             'Robust="huber" (or igg3, bisquare) selects M-estimation on every kind. The same Results comes back; '
             'residuals() and contribution() refuse without an accumulator and '
             'say so.', 'body')]
    st += [callout('A stream whose observation type is called pos is treated as 3 '
                   'components per epoch by default, which selects the separable '
                   'model and needs a per-epoch scaling. A scalar stream that '
                   'happens to be called pos must declare '
                   'ComponentsPerEpoch = struct(\'pos\', 1).', 'note')]
    st += [PageBreak()]
    return st
