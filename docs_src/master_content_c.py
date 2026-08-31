#!/usr/bin/env python3
"""Master manual, Parts 5-6 and the references.

Part 5 (interpretation) and the symbols table in Part 6 are consolidated from
the former interpretation guide and reference manual rather than rewritten --
that material was already correct, and re-authoring it would only introduce
drift.

Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 4.8.
"""
from reportlab.platypus import Spacer, PageBreak
from reportlab.lib.units import mm

from build_master_pdf import (P, B, N, code, callout, table, caption, S,
                              BODY_W)
from pdf_graphics import eq, workflow, decision_tree, rule


def part5_interpretation():
    st = []
    st += [P('5&nbsp;&nbsp;Interpreting the result', 'h1')]
    st += [P('This part is about judgement rather than syntax: which number to '
             'look at first, what a red flag means, and the two ways least '
             'squares will quietly mislead you. The equations behind it are in '
             'Part 2.', 'body')]

    # ---- 5.1 --------------------------------------------------------------
    st += [P('5.1&nbsp;&nbsp;The order to read things', 'h2')]
    st += [P('<font face="Courier">R.report()</font> is organised in five '
             'levels. <b>Read it top-down and stop at the first failure</b> -- '
             'a break at one level makes every level below it meaningless.',
             'body')]
    st += [table([
        ['Level', 'Reads', 'If red', 'Why it comes first'],
        ['<b>0 Integrity</b>', 'the trace identity',
         'a data or weighting bug',
         'if &#8721;<i>r<sub>i</sub></i> &#8800; <i>n</i>&minus;<i>u</i> the '
         'geometry is wrong and nothing else is trustworthy'],
        ['<b>1 Global</b>', 'model consistency', 'inconsistent',
         '<i>s</i><sub>0</sub><sup>2</sup>/&#963;<sub>0</sub><sup>2</sup> too '
         'large means blunders or optimistic weights'],
        ['<b>2 Localise</b>', 'find the culprit', 'snoop',
         'but check the contamination gauge first'],
        ['<b>3 Reliability</b>', 'what you cannot check', 'you are blind here',
         'low <i>r<sub>i</sub></i> means unverifiable, whatever the residual'],
        ['<b>4 Parameters</b>', 'result fragility', 'results fragile',
         'high correlation or conditioning means weakly separated unknowns'],
    ], widths=[24 * mm, 26 * mm, 30 * mm, BODY_W - 80 * mm])]
    st += [caption('Table 5.1 -- the five levels of the verdict. A Level 0 red '
                   'is not a blunder; it is a bug in what you supplied.')]

    st += [P('Why Level 0 is an identity, not a test', 'h3')]
    st += [P('The redundancy numbers are the diagonal of a projector '
             '(equation 2.5), and the trace of a projector equals its rank. So '
             'their sum is not an estimate that might be close -- it is fixed '
             'by the geometry:', 'body')]
    st += [eq([('g', 'sum'), ('sub', 'i'), ('sp', ''),
               ('v', 'r'), ('sub', 'i'), ('sp', ''), ('r', '= tr('),
               ('b', 'R'), ('b', 'Q'), ('r', ') = '), ('v', 'n'), ('sp', ''),
               ('g', 'minus'), ('sp', ''), ('v', 'u')],
              width=BODY_W, number='5.1')]
    st += [P('Any departure beyond rounding means the weights, the design '
             'matrix or the rank are not what the toolbox was told. That is why '
             'it is read first: no later number can be salvaged from a system '
             'whose geometry is misreported.', 'body')]

    st += [P('What Level 1 actually tests', 'h3')]
    st += [P('The global test compares the a-posteriori variance factor with '
             'the a-priori one. Under the model, the weighted sum of squares is '
             'chi-square distributed with <i>r</i> degrees of freedom:', 'body')]
    st += [eq([('b', 'v'), ('sup', 'T'), ('b', 'P'), ('b', 'v'), ('sp', ''),
               ('r', '/ '), ('g', 'sigma'), ('sub', '0'), ('sup', '2'),
               ('sp', ''), ('r', '~'), ('sp', ''),
               ('g', 'chi'), ('sup', '2'), ('r', '('), ('v', 'r'), ('r', ')'),
               ('r', ','), ('sp', ''), ('sp', ''),
               ('v', 's'), ('sub', '0'), ('sup', '2'), ('r', ' / '),
               ('g', 'sigma'), ('sub', '0'), ('sup', '2'), ('sp', ''),
               ('r', '= '), ('b', 'v'), ('sup', 'T'), ('b', 'P'), ('b', 'v'),
               ('r', ' / ('), ('v', 'r'), ('sp', ''), ('g', 'sigma'),
               ('sub', '0'), ('sup', '2'), ('r', ')')],
              width=BODY_W, number='5.2')]
    st += [P('A ratio far above one has exactly three causes, and they are not '
             'distinguishable from this number alone: gross errors in the '
             'observations, weights that were too optimistic, or a model that '
             'is missing a term. That ambiguity is why Level 1 tells you to '
             'look further rather than what to do -- section 5.4 separates the '
             'three.', 'body')]
    st += [callout(
        '<b>A ratio far BELOW one is also a failure</b>, and it is the one '
        'people ignore. It means the weights were pessimistic: the solution '
        'fits better than the stated precision allows. The estimate is not '
        'wrong, but every covariance derived from it is too large, so formal '
        'errors and MDBs are overstated and the network looks weaker than it '
        'is.', 'warn')]

    # ---- 5.2 --------------------------------------------------------------
    st += [P('5.2&nbsp;&nbsp;Two ways least squares lies to you', 'h2')]
    st += [P('Lie 1 -- "the residual is small, so the observation is good"',
             'h3')]
    st += [P('False. A residual shows only the fraction <i>r<sub>i</sub></i> of '
             'the error present:', 'body')]
    st += [eq([('r', 'E{'), ('g', 'Delta'), ('v', 'v'), ('sub', 'i'),
               ('r', '} = '), ('g', 'minus'), ('v', 'r'), ('sub', 'i'),
               ('sp', ''), ('v', 'b'), ('sub', 'i')],
              width=BODY_W, number='5.3')]
    st += [P('At <i>r<sub>i</sub></i> &#8776; 0 even a huge blunder leaves a '
             'near-zero residual. <b>Always pair a residual with its redundancy '
             'number.</b> A small residual at low <i>r<sub>i</sub></i> is the '
             'most dangerous reading in the whole report, not the most '
             'reassuring.', 'body')]
    st += [P('The mechanism is leverage. In whitened form the redundancy number '
             'is the complement of the leverage <i>h<sub>i</sub></i>, the '
             'observation\'s own influence on its own fitted value:', 'body')]
    st += [eq([('v', 'r'), ('sub', 'i'), ('sp', ''), ('r', '= 1 '),
               ('g', 'minus'), ('r', ' '), ('v', 'h'), ('sub', 'i'),
               ('r', ','), ('sp', ''), ('sp', ''),
               ('v', 'h'), ('sub', 'i'), ('sp', ''), ('r', '= '),
               ('b', 'a'), ('sub', 'i'), ('sup', 'T'), ('b', 'N'),
               ('sup', '\u22121'), ('b', 'a'), ('sub', 'i')],
              width=BODY_W, number='5.4')]
    st += [P('An observation with <i>h<sub>i</sub></i> &#8594; 1 determines its '
             'own fitted value: the adjustment simply follows it, so it cannot '
             'disagree with itself and its residual goes to zero whatever its '
             'error. This is geometric, not statistical -- no amount of data '
             'quality control detects it, and no test can recover information '
             'the network never collected. It is why reliability (Level 3) is '
             'reported separately from residuals (Level 2).', 'body')]

    st += [P('Lie 2 -- "the tau test found no outliers, so the data are clean"',
             'h3')]
    st += [P('Baarda\'s <i>w</i> divides by a KNOWN &#963;<sub>0</sub> and is '
             'normally distributed. Pope\'s &#964; divides by the ESTIMATED '
             '<i>s</i><sub>0</sub> -- which is the honest thing to do when '
             '&#963;<sub>0</sub> is unknown, and is also the trap:', 'body')]
    st += [eq([('g', 'tau'), ('sub', 'i'), ('sp', ''), ('r', '= '),
               ('v', 'v'), ('sub', 'i'), ('r', ' / ('), ('v', 's'),
               ('sub', '0'), ('g', 'sqrt'), ('v', 'r'), ('sub', 'i'),
               ('g', 'sigma'), ('sub', 'i'), ('r', ')'), ('r', ','),
               ('sp', ''), ('sp', ''),
               ('v', 's'), ('sub', '0'), ('sup', '2'), ('sp', ''),
               ('r', '= '), ('b', 'v'), ('sup', 'T'), ('b', 'P'), ('b', 'v'),
               ('r', ' / '), ('v', 'r')],
              width=BODY_W, number='5.5')]
    st += [P('The blunders are in <b>v</b>, so they are in '
             '<i>s</i><sub>0</sub>, so they are in the denominator of their own '
             'test statistic. One large error inflates <i>s</i><sub>0</sub> '
             'enough to shrink its own &#964;; several inflate it enough to '
             'hide each other. Past a few per cent contamination the test goes '
             'blind, and a 100&#963; error can pass.', 'body')]
    st += [callout(
        '<b>This is also why rejection is ONE AT A TIME.</b> The usual '
        'justification -- that a blunder smears its residual onto neighbours '
        'through the correlation in <b>Q</b><sub>vv</sub> -- was tested '
        'directly and is not what dominates: with a Bonferroni threshold, '
        'simultaneous rejection produced no false positives. The real reason is '
        'MASKING. Remove one blunder, <i>s</i><sub>0</sub> drops, and the next '
        'becomes visible. Reject several at once and you are still testing '
        'against the inflated scale.', 'note')]
    st += [P('The defence is to estimate the scale robustly, so the blunders '
             'cannot set the threshold that is supposed to catch them. Read the '
             'contamination gauge first, and switch to '
             '<font face="Courier">\'test\',\'mad\'</font> when it is high.',
             'body')]

    # ---- 5.3 --------------------------------------------------------------
    st += [P('5.3&nbsp;&nbsp;The quadrant you must inspect', 'h2')]
    st += [table([
        ['', 'low |<i>w</i>| (quiet)', 'high |<i>w</i>| (loud)'],
        ['<b>high <i>r</i></b> (checked)', 'healthy',
         'plain outlier -- snoop it'],
        ['<b>low <i>r</i></b> (unchecked)',
         '<b>blind spot</b> -- do not relax',
         'leverage point -- a small <i>w</i> still hides a big error'],
    ], widths=[42 * mm, 48 * mm, BODY_W - 90 * mm])]
    st += [caption('Table 5.2 -- residual against redundancy. The two are only '
                   'interpretable together.')]
    st += [callout(
        '<b>The bottom row is where reputations are lost.</b> If a parameter '
        'you care about depends on bottom-row observations, no amount of clean '
        'residuals should reassure you. Strengthen the geometry or add '
        'independent data -- there is no statistical treatment that recovers '
        'information the network never collected.', 'danger')]

    st += [P('Quantifying the blind spot: the MDB', 'h3')]
    st += [P('The quadrant is qualitative; the minimal detectable bias makes it '
             'a number. Invert the test of equation 2.7: how large must an '
             'error be before it is found with probability '
             '&#946;<sub>0</sub>?', 'body')]
    st += [eq([('r', 'MDB'), ('sub', 'i'), ('sp', ''), ('r', '= '),
               ('g', 'delta'), ('sub', '0'), ('r', '('), ('g', 'alpha'),
               ('r', ', '), ('g', 'beta'), ('sub', '0'), ('r', ') '),
               ('g', 'sigma'), ('sub', 'i'), ('r', ' / '), ('g', 'sqrt'),
               ('v', 'r'), ('sub', 'i')],
              width=BODY_W, number='5.6')]
    st += [P('&#948;<sub>0</sub> is the non-centrality the chosen significance '
             'and power demand -- about 4.13 for &#945; = 0.001 and '
             '&#946;<sub>0</sub> = 0.80. Everything problem-specific sits in '
             '&#963;<i><sub>i</sub></i>/&#8730;<i>r<sub>i</sub></i>, and the '
             'square root is the whole story: at <i>r<sub>i</sub></i> = 0.5 the '
             'MDB is 1.4 &#963;<i><sub>i</sub></i>; at 0.05 it is 4.5; at 0.01 '
             'it is 10. <b>Halving the redundancy number does not halve what '
             'you can detect -- it multiplies it by &#8730;2.</b>', 'body')]
    st += [callout(
        '<b>Report MDBs, not just residuals.</b> An MDB says what the network '
        'CANNOT rule out, which is the question a reader of your results '
        'actually has. A table of small residuals answers a different and '
        'easier question. And compute them at the CONVERGED weights: an MDB '
        'from a mis-weighted solve is wrong by exactly the factor the variance '
        'components were about to correct.', 'warn')]

    # ---- 5.4 --------------------------------------------------------------
    st += [P('5.4&nbsp;&nbsp;Deciding what to do next', 'h2')]
    st += [P('Level 1 failed -- blunders, or the model?', 'h3')]
    st += [P('Run <font face="Courier">lsqvce.engine.systematic(D)</font> <i>before</i> '
             'snooping. If the residuals show runs, serial correlation or '
             'per-group offsets, the problem is the model -- an unmodelled '
             'signal, a wrong weight ratio, coloured noise treated as white. '
             'Snooping that away merely deletes good data from around a defect '
             'that remains. Fix the model, then re-solve.', 'body')]

    st += [P('Diffuse contamination -- reject or reweight?', 'h3')]
    st += [P('Isolated blunders call for iterative '
             '<font face="Courier">lsqvce.engine.snoop</font>. Many small '
             'contaminations call for robust M-estimation '
             '(<font face="Courier">lsqvce.engine.robust</font>), which degrades '
             'gracefully where snooping would chase its own tail. The '
             'contamination gauge tells you which regime you are in.', 'body')]

    st += [P('Group weights look wrong -- variance components', 'h3')]
    st += [P('If one data type or satellite is systematically over- or '
             'under-fit, read the variance components. A factor far from 1 '
             'means that group was mis-weighted going in.', 'body')]
    st += [callout(
        '<b>Re-solve at the converged weights before reading reliability.</b> '
        'MDBs computed at the wrong weights are wrong. Confirm the loop '
        'converged and that the minimum redundancy was met for every group '
        '-- <font face="Courier">R.summary()</font> names the streams, '
        '<font face="Courier">R.convergence()</font> shows whether they '
        'settled.', 'warn')]

    # ---- 5.5 --------------------------------------------------------------
    st += [P('5.5&nbsp;&nbsp;Reading a gravity solution', 'h2')]
    st += [P('For the block-arrow path, three habits keep you honest.', 'body')]
    st += [N('<b>Diagnose through the accumulator bridge</b> '
             '(<font face="Courier">lsqvce.engine.diagnoseAccum</font>), never the '
             'reduced common system alone. The common-parameter redundancies '
             'are only correct once the arc-local parameters have been counted; '
             'skip the bridge and every reliability number is inflated.', 1)]
    st += [N('<b>Map the redundancy of the field</b>, not just its residuals. '
             'A region carried by a single arc is a blind spot: its formal '
             'error can look fine while being uncheckable.', 2)]
    st += [N('<b>Watch for a near-singular arc.</b> An arc with too few '
             'observations to fix its own state poisons the Schur complement '
             '(equation 2.10). The accumulator flags it -- do not force it in.',
             3)]
    st += [callout(
        '<b>A clean bill of health means all of:</b> Level 0 exact, Level 1 '
        'consistent, no high-<i>w</i> observations left, contamination low, and '
        'the parameters you care about backed by high-redundancy observations '
        'with acceptable MDB. Any one of those missing turns "passed" into '
        '"passed, but&hellip;".', 'good')]
    st += [PageBreak()]
    return st


def part6_reference():
    st = []
    st += [P('6&nbsp;&nbsp;Reference', 'h1')]

    # ---- 6.0 installation -------------------------------------------------
    st += [P('6.1&nbsp;&nbsp;Installation', 'h2')]
    st += [P('The toolbox is plain MATLAB. There is nothing to compile, no MEX '
             'file, and no runtime dependency outside MATLAB itself -- the '
             'Python in <font face="Courier">docs_src/</font> builds this '
             'manual and is never called from MATLAB.', 'body')]
    st += [P('Steps', 'h3')]
    st += [N('<b>Unpack</b> the toolbox anywhere you have write access. The '
             'folder holding <font face="Courier">+lsqvce</font> is the '
             '<i>toolbox root</i>; every instruction below assumes MATLAB\'s '
             'working directory is that root.', 1)]
    st += [N('<b>Run the setup</b> once:', 2)]
    st += [code("""cd /path/to/lsqvce
lsqvce_setup""")]
    st += [P('This adds the root (making the '
             '<font face="Courier">+lsqvce</font> package visible), adds '
             '<font face="Courier">tests/</font>, '
             '<font face="Courier">examples/</font> and '
             '<font face="Courier">template/</font>, and rebuilds the Help '
             'browser index so <font face="Courier">doc lsqvce</font> works. '
             'It should run <b>without warnings</b>; an '
             '<font face="Courier">addpath</font> warning means a folder is '
             'missing and the unpack was incomplete.', 'body')]
    st += [N('<b>Verify</b>, before trusting anything:', 3)]
    st += [code("""run_all_lsqvce         % expect 497/498 passing, 0 failed, ~130 s""")]
    st += [P('One test is reported <i>incomplete</i>, not failed: it is filtered '
             'by assumption because it needs a parallel pool. <b>0 failed</b> '
             'is the number that matters.', 'body')]
    st += [N('<b>Save the path</b> if you want it to persist across MATLAB '
             'sessions -- <font face="Courier">lsqvce_setup</font> modifies the '
             'path for the current session only:', 4)]
    st += [code("""savepath               % optional, and it is a global change""")]

    st += [P('Requirements', 'h3')]
    st += [table([
        ['Component', 'Needed for', 'Required?'],
        ['MATLAB R2019b or newer',
         'the toolbox itself (string arrays, <font face="Courier">arguments</font> '
         'blocks are not used, so older releases may work but are untested)',
         'yes'],
        ['Parallel Computing Toolbox',
         'the parallel accumulation paths only',
         'no -- everything falls back to serial'],
        ['<font face="Courier">reportlab</font> + FreeSerif fonts',
         'rebuilding <i>this manual</i>',
         'no -- build-time only'],
    ], widths=[44 * mm, BODY_W - 74 * mm, 30 * mm])]

    st += [callout(
        '<b>Nothing in the toolbox writes outside its own folder</b> except '
        '<font face="Courier">savepath</font>, which you invoke deliberately, '
        'and the test log <font face="Courier">tests/last_run_all.txt</font>. '
        'There is no installer, no registry entry and no network access at any '
        'point.', 'note')]

    st += [P('First steps after installing', 'h3')]
    st += [table([
        ['Command', 'What it does'],
        ['<font face="Courier">doc lsqvce</font>',
         'the landing page in the Help browser'],
        ['<font face="Courier">demo_workflow</font>',
         'a guided walk-through, section by section'],
        ['<font face="Courier">demo_grace</font>',
         'the worked GRACE-FO case, end to end'],
        ['<font face="Courier">edit lsqvce_template</font>',
         'start your own case (Part 4)'],
        ['<font face="Courier">run_lsqvce(\'adjustvce\')</font>',
         'one test suite, while iterating'],
    ], widths=[58 * mm, BODY_W - 58 * mm])]


    # ---- 6.1 symbols ------------------------------------------------------
    st += [P('6.2&nbsp;&nbsp;Symbols', 'h2')]
    st += [table([
        ['Symbol', 'Meaning'],
        ['<i>n</i>, <i>u</i>, <i>r</i>',
         'observations, parameters, redundancy <i>r</i> = <i>n</i>&minus;<i>u</i>'],
        ['<b>A</b>, <b>l</b>, <b>e</b>',
         'design (Jacobian), observations, errors'],
        ['<b>C</b>, <b>P</b>',
         'observation cofactor and weight, <b>P</b> = <b>C</b><sup>&minus;1</sup>'],
        ['&#963;<sub>0</sub><sup>2</sup>, <i>s</i><sub>0</sub><sup>2</sup>',
         'a-priori and a-posteriori variance of unit weight'],
        ['<b>N</b>, <b>b</b>',
         'normal matrix <b>A</b><sup>T</sup><b>PA</b> and right side '
         '<b>A</b><sup>T</sup><b>Pl</b>'],
        ['<b>Q</b><sub>xx</sub> = <b>N</b><sup>&minus;1</sup>',
         'parameter cofactor'],
        ['<b>v</b>, <b>Q</b><sub>vv</sub>', 'residuals, residual cofactor'],
        ['<i>r<sub>i</sub></i>, <i>h<sub>i</sub></i>',
         'redundancy number and leverage, <i>r<sub>i</sub></i> = '
         '1 &minus; <i>h<sub>i</sub></i> (whitened)'],
        ['<i>w<sub>i</sub></i>, &#964;<sub><i>i</i></sub>',
         'Baarda standardised residual, Pope tau'],
        ['<i>n<sub>c</sub></i>, <i>n<sub>l</sub></i>, <i>K</i>',
         'common parameters, local parameters per arc, number of arcs'],
        ['<b>D</b>, <b>R</b><sub>t</sub>, <i>d</i>',
         'per-epoch spatial factor, temporal correlation, number of axes'],
    ], widths=[42 * mm, BODY_W - 42 * mm])]

    # ---- 6.2 equations ----------------------------------------------------
    st += [P('6.3&nbsp;&nbsp;Equation summary', 'h2')]
    st += [table([
        ['Eq', 'Subject', 'Section'],
        ['2.1&ndash;2.3', 'Gauss-Markov model, BLUE, dispersion', '2.1'],
        ['2.4', 'QR of the whitened design matrix', '2.2'],
        ['2.5&ndash;2.7', 'residuals, redundancy numbers, Baarda / MDB', '2.3'],
        ['2.8&ndash;2.9', 'variance components, Helmert iteration', '2.4'],
        ['2.10', 'Schur complement, reduced common system', '2.5'],
        ['2.11', 'autoregressive noise model', '2.6'],
        ['2.12&ndash;2.13', 'separable covariance, epoch indexing', '2.7'],
        ['5.1', 'expected residual under a bias', '5.2'],
    ], widths=[24 * mm, BODY_W - 54 * mm, 30 * mm])]

    # ---- 6.3 error ids ----------------------------------------------------
    st += [P('6.4&nbsp;&nbsp;Error identifiers you are likely to meet', 'h2')]
    st += [table([
        ['Identifier', 'Means'],
        ['<font face="Courier">lsqvce:makeLoader:badDay</font>',
         'your <font face="Courier">buildDay</font> omitted a required field'],
        ['<font face="Courier">lsqvce:makeLoader:dayShape</font>',
         '<b>A</b>, <b>B</b>, <b>l</b>, <b>t</b> do not share a row count'],
        ['<font face="Courier">lsqvce:covModel:badCep</font>',
         '<font face="Courier">Cep</font> is neither <i>nEpoch</i>&#215;<i>d</i> '
         'nor <i>nEpoch</i>&#215;<i>d</i>&#215;<i>d</i>'],
        ['<font face="Courier">lsqvce:covModel:badSigma</font>',
         'a separable model needs a positive scalar innovation std'],
        ['<font face="Courier">lsqvce:Results:noHistory</font>',
         'you called an iteration reader on a single-solve result'],
        ['<font face="Courier">lsqvce:verbosity:badLevel</font>',
         'verbosity outside 0&ndash;2 (level 3 was merged into 2)'],
        ['<font face="Courier">lsqvce_template:notImplemented</font>',
         'a template stub was called before you wrote it'],
    ], widths=[70 * mm, BODY_W - 70 * mm])]

    # ---- 6.4 testing ------------------------------------------------------
    st += [P('6.5&nbsp;&nbsp;Testing and the suite', 'h2')]
    st += [code("""run_all_lsqvce              % the whole suite, with a full log
run_lsqvce('adjustvce')     % one suite, while iterating
run_lsqvce('*degenerate*')  % one test, wherever it lives""")]
    st += [callout(
        '<b>Do not use <font face="Courier">runtests</font> on a folder.</b> It '
        'SKIPS any file that fails to load and still reports the rest as green, '
        'so the suite can look clean while an entire file never ran. '
        '<font face="Courier">run_all_lsqvce</font> iterates the registered '
        'list, counts a load failure as a failure, keeps the warning inventory '
        '(warnings are a signal, not noise) and records per-suite seconds.',
        'warn')]
    st += [P('The suite deliberately provokes warnings in its failure-mode '
             'tests, so a warning-free run would itself be suspicious. Compare '
             'the warning inventory and the per-suite timings between runs, not '
             'just the pass count.', 'body')]
    st += [PageBreak()]
    return st


def references():
    st = []
    st += [P('References', 'h1')]
    refs = [
        'Baarda, W. (1968). <i>A Testing Procedure for Use in Geodetic '
        'Networks.</i> Netherlands Geodetic Commission, Publications on '
        'Geodesy, New Series, 2(5). Delft.',
        'Pope, A. J. (1976). <i>The Statistics of Residuals and the Detection '
        'of Outliers.</i> NOAA Technical Report NOS 65 NGS 1. Rockville, MD.',
        'Koch, K.-R. (1999). <i>Parameter Estimation and Hypothesis Testing in '
        'Linear Models</i>, 2nd ed. Springer, Berlin.',
        'Teunissen, P. J. G. (2006). <i>Testing Theory: an Introduction</i>, '
        '2nd ed. VSSD, Delft.',
        'Teunissen, P. J. G. (2003). <i>Adjustment Theory: an Introduction.</i> '
        'VSSD, Delft.',
        'F&ouml;rstner, W. (1979). Ein Verfahren zur Sch&auml;tzung von '
        'Varianz- und Kovarianzkomponenten. <i>Allgemeine '
        'Vermessungs-Nachrichten</i>, 86, 446&ndash;453.',
        'Helmert, F. R. (1907). <i>Die Ausgleichungsrechnung nach der Methode '
        'der kleinsten Quadrate</i>, 2nd ed. Teubner, Leipzig.',
        'Grafarend, E. W. &amp; Schaffrin, B. (1993). <i>Ausgleichungsrechnung '
        'in linearen Modellen.</i> BI-Wissenschaftsverlag, Mannheim.',
        'Koch, K.-R. &amp; Kusche, J. (2002). Regularization of geopotential '
        'determination from satellite data by variance components. <i>Journal '
        'of Geodesy</i>, 76(5), 259&ndash;268.',
        'Kusche, J. (2003). A Monte-Carlo technique for weight estimation in '
        'satellite geodesy. <i>Journal of Geodesy</i>, 76(11&ndash;12), '
        '641&ndash;652.',
        'Golub, G. H. &amp; Van Loan, C. F. (2013). <i>Matrix Computations</i>, '
        '4th ed. Johns Hopkins University Press, Baltimore.',
        'Hutchinson, M. F. (1990). A stochastic estimator of the trace of the '
        'influence matrix for Laplacian smoothing splines. <i>Communications '
        'in Statistics -- Simulation and Computation</i>, 19(2), '
        '433&ndash;450.',
        'Golub, G. H., Heath, M. &amp; Wahba, G. (1979). Generalized '
        'cross-validation as a method for choosing a good ridge parameter. '
        '<i>Technometrics</i>, 21(2), 215&ndash;223.',
        'Huber, P. J. (1981). <i>Robust Statistics.</i> Wiley, New York.',
        'Levinson, N. (1947). The Wiener RMS error criterion in filter design '
        'and prediction. <i>Journal of Mathematics and Physics</i>, 25, '
        '261&ndash;278.',
        'Durbin, J. (1959). Efficient estimation of parameters in '
        'moving-average models. <i>Biometrika</i>, 46(3/4), 306&ndash;316.',
        'Ljung, G. M. &amp; Box, G. E. P. (1978). On a measure of lack of fit '
        'in time series models. <i>Biometrika</i>, 65(2), 297&ndash;303.',
        'Wen, H. Y., Kruizinga, G., Paik, M., et al. (2019). <i>GRACE-FO '
        'Level-1 Data Product User Handbook.</i> JPL D-56935, NASA Jet '
        'Propulsion Laboratory.',
        'Kornfeld, R. P., Arnold, B. W., Gross, M. A., et al. (2019). GRACE-FO: '
        'the Gravity Recovery and Climate Experiment Follow-On mission. '
        '<i>Journal of Spacecraft and Rockets</i>, 56(3), 931&ndash;951.',
        'Schaffrin, B. (1997). Reliability measures for correlated '
        'observations. <i>Journal of Surveying Engineering</i>, 123(3), '
        '126&ndash;137.',
        'Ellmer, M. (2018). <i>Contributions to GRACE Gravity Field Recovery: '
        'Improvements in Dynamic Orbit Integration, Stochastic Modelling of '
        'Observations, and Variance Component Estimation.</i> PhD thesis, TU '
        'Graz.',
        'Kvas, A. (2019). <i>Stochastic Modelling in Satellite Gravimetry.</i> '
        'PhD thesis, TU Graz.',
    ]
    for i, r in enumerate(refs, 1):
        st += [P('[%d]&nbsp;&nbsp;%s' % (i, r), 'ref')]
    st += [Spacer(1, 14)]
    st += [rule(BODY_W)]
    st += [P('lsqvce 1.0 -- consolidated diagnostics and streaming variance '
             'component estimation. This manual was generated with assistance '
             'from Claude (Anthropic), model Claude Opus 4.8. Build it with '
             '<font face="Courier">python3 docs_src/build_master_pdf.py</font>; '
             'see <font face="Courier">docs_src/README_manual.md</font>.',
             'caption')]
    return st
