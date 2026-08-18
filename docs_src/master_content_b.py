#!/usr/bin/env python3
"""Master manual, Parts 3-4: the entry points, and the template guide.

Content module for build_master_pdf.py. Part 4 is the step-by-step guide for
adapting template/lsqvce_template.m to a real case; the user asked for special
attention there, so every step states what you write, what the toolbox derives,
and what goes wrong if you get it subtly wrong.

Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 4.8.
"""
from reportlab.platypus import Spacer, PageBreak
from reportlab.lib.units import mm

from build_master_pdf import (P, B, N, code, callout, table, caption, S,
                              BODY_W)
from pdf_graphics import eq, workflow, decision_tree, rule


def part3_entry_points():
    st = []
    st += [P('3&nbsp;&nbsp;The entry points, step by step', 'h1')]
    st += [P('Each engine is shown as a workflow: what you supply, what the '
             'toolbox does with it, and what comes back. The shaded rounded '
             'boxes are yours; the plain boxes are the toolbox.', 'body')]

    # ---- 3.1 adjust -------------------------------------------------------
    st += [P('3.1&nbsp;&nbsp;adjust -- flat, in memory', 'h2')]
    st += [P('The direct path. You supply a batch function that hands over '
             'design rows, observations and weights; <b>adjust</b> accumulates '
             'the normal equations, solves, and returns the solution together '
             'with its diagnostics. One pass, no iteration, no assumed '
             'structure.', 'body')]
    st += [workflow([
        [{'t': 'batchFcn(k) -> A, l, p', 'k': 'io'}],
        [{'t': 'accumulate normal equations', 'k': 'proc'}],
        [{'t': 'solve', 'k': 'proc'}],
        [{'t': 'second pass: residuals + redundancy', 'k': 'proc'}],
        [{'t': 'sol, D, V', 'k': 'io'}],
    ], width=BODY_W, title='Figure 3.1 -- adjust')]
    st += [caption('batchFcn is called TWICE -- once to accumulate, once to '
                   'diagnose -- and must return identical data both times. A '
                   'mismatch is an error, not a silent misattribution.')]
    st += [code("""[sol, D, V] = lsqvce.adjust(u, batchFcn, nBatch);
lsqvce.report(V);                  % the level-0..4 verdict""")]
    st += [P('The returns are structs, not a Results object: <font '
             'face="Courier">sol</font> carries the estimate and its '
             'dispersion, <font face="Courier">D</font> the diagnostics, '
             '<font face="Courier">V</font> the verdict. Because there is no '
             'iteration there is no history, so the iteration readers of '
             'section 3.5 do not apply here.', 'body')]

    # ---- 3.2 adjustBlock --------------------------------------------------
    st += [P('3.2&nbsp;&nbsp;adjustBlock -- block-arrow, one solve', 'h2')]
    st += [P('The same single solve, but exploiting block-arrow structure and '
             'streaming the data. You supply <i>loaders</i> rather than a batch '
             'function; each loader knows how to build one day of one '
             'observation type. The accumulator slices arcs, applies the '
             'a-priori covariance, forms each arc\'s Schur complement and adds '
             'it to the reduced common system.', 'body')]
    st += [workflow([
        [{'t': 'loaders (makeLoader)', 'k': 'io'}],
        [{'t': 'for each day: buildDay -> A, B, l, t', 'k': 'proc'}],
        [{'t': 'slice arcs, whiten with the a-priori covariance', 'k': 'proc'}],
        [{'t': 'Schur complement per arc -> reduced common system', 'k': 'proc'}],
        [{'t': 'solve the common block, back-substitute the locals', 'k': 'proc'}],
        [{'t': 'Results', 'k': 'io'}],
    ], width=BODY_W, title='Figure 3.2 -- adjustBlock')]
    st += [code("""R = lsqvce.adjustBlock({Lc, Ld}, days);
R.report();""")]
    st += [callout(
        '<b>The noise model is taken as given.</b> Whatever covariance you '
        'declared on the loaders is used as-is. If the residuals come out '
        'coloured, adjustBlock will not tell you to fix it -- that is what '
        'adjustVCE is for.', 'note')]

    # ---- 3.3 adjustVCE ----------------------------------------------------
    st += [P('3.3&nbsp;&nbsp;adjustVCE -- the iterative loop', 'h2')]
    st += [P('The same solve as adjustBlock, run repeatedly, with two '
             'estimation steps between passes. After each solve the residuals '
             'are examined: a temporal noise model is fitted per stream (the '
             '<i>shape</i> of the noise), and the variance components are '
             're-estimated per group (its <i>scale</i>). Both feed back into '
             'the next solve. The loop stops when the variance factors are '
             'near one and the residuals test white.', 'body')]
    st += [workflow([
        [{'t': 'loaders + a-priori covariance', 'k': 'io'}],
        [{'t': 'solve (as adjustBlock)', 'k': 'proc'}],
        [{'t': 'fit the temporal model per stream (AR / ARMA)', 'k': 'proc'}],
        [{'t': 'estimate variance components per group', 'k': 'proc'}],
        [{'t': 'converged? sigma ~ 1 and residuals white', 'k': 'dec'}],
        [{'t': 'Results + history', 'k': 'io'}],
    ], width=BODY_W, title='Figure 3.3 -- adjustVCE',
        loop=(4, 1, 'not yet'))]
    st += [caption('The orange arrow is the outer loop. Each pass re-solves '
                   'with the noise model and weights discovered by the '
                   'previous one; R.history records every pass, which is what '
                   'the two readers in section 3.5 report on.')]
    st += [code("""R = lsqvce.adjustVCE(L, days) ...
      .groupComponentsBy('obsType') ...
      .run();
R.report();
R.summary();          % per stream: what it ended up as
R.convergence();      % per iteration: when it settled""")]
    st += [P('The fluent calls configure the loop before it runs. '
             '<font face="Courier">groupComponentsBy</font> sets the '
             'granularity of the variance components; '
             '<font face="Courier">stochasticModel</font> and '
             '<font face="Courier">varianceComponents</font> tune the two '
             'estimation steps. All have documented defaults -- start without '
             'them and change one only when a diagnostic tells you to.', 'body')]

    # ---- 3.4 snoop --------------------------------------------------------
    st += [P('3.4&nbsp;&nbsp;snoop -- outlier rejection', 'h2')]
    st += [callout(
        '<b>snoop is a SEPARATE step, and it works on the flat path only.</b> '
        'It is not part of any chain in 3.1-3.3: no engine calls it, and '
        '<font face="Courier">adjustBlock</font> and '
        '<font face="Courier">adjustVCE</font> have no equivalent. Its '
        'signature says so -- it takes <font face="Courier">sol</font> and '
        '<font face="Courier">D</font>, the structs '
        '<font face="Courier">adjust</font> returns, plus the same '
        '<font face="Courier">batchFcn</font>, because it must re-solve after '
        'every rejection.', 'warn')]
    st += [P('Data snooping tests each residual against its own standard '
             'deviation (equation 2.7), rejects the single worst, re-solves, '
             'and repeats until nothing trips the critical value.', 'body')]
    st += [workflow([
        [{'t': 'sol, D from adjust', 'k': 'io'}],
        [{'t': 'compute w_i for every surviving observation', 'k': 'proc'}],
        [{'t': 'any |w| above the critical value?', 'k': 'dec'}],
        [{'t': 'reject the worst ONE, re-solve', 'k': 'proc'}],
        [{'t': 'sol, D, REJ', 'k': 'io'}],
    ], width=BODY_W, title='Figure 3.4 -- snoop, an iteration around adjust',
        loop=(3, 1, 'repeat'))]
    st += [caption('One at a time, and the reason is MASKING: several blunders '
                   'inflate s0 together, so each hides the others -- every |w| '
                   'is divided by a variance factor the blunders themselves '
                   'raised. Remove one, s0 drops, and the next becomes '
                   'visible.')]
    st += [code("""[sol, D]      = lsqvce.adjust(u, batchFcn, nBatch);     % first solve
[sol, D, REJ] = lsqvce.snoop(sol, D, batchFcn, 'test', 'mad');""")]
    st += [callout(
        '<b>On the block-arrow path there is no snoop.</b> Gross-error handling '
        'there is a different problem: the observations are streamed and never '
        'all resident, and an arc whose local block is near-singular is a '
        'model defect rather than an outlier. Read the verdict from '
        '<font face="Courier">R.report()</font> and the per-stream state from '
        '<font face="Courier">R.summary()</font>, and treat a stream that will '
        'not whiten as a modelling question -- section 4.9.', 'note')]
    st += [callout(
        '<b>Localise before rejecting.</b> A structured error is a model '
        'defect, and snooping it away mutilates good data. Run '
        '<font face="Courier">lsqvce.systematic(D)</font> first: if the '
        'residuals show a pattern, fix the model, not the observations. And '
        'read section 5.3 before trusting any rejection at low redundancy.',
        'warn')]

    # ---- 3.5 Results ------------------------------------------------------
    st += [P('3.5&nbsp;&nbsp;Reading a Results object', 'h2')]
    st += [P('<b>adjustBlock</b> and <b>adjustVCE</b> both return a '
             '<font face="Courier">lsqvce.Results</font>. It carries the '
             'solution, the accumulator, and -- for the iterative engine -- the '
             'full history.', 'body')]
    st += [table([
        ['Call', 'Returns'],
        ['<font face="Courier">R.report()</font>',
         'the level-0..4 verdict; read it top-down and stop at the first red'],
        ['<font face="Courier">R.redundancy()</font>',
         'total redundancy, <i>n</i> &minus; <i>u</i>'],
        ['<font face="Courier">R.residuals(...)</font>',
         'whitened or raw residuals'],
        ['<font face="Courier">R.contribution(...)</font>',
         'per-group information share of the solution'],
        ['<font face="Courier">R.summary()</font>',
         'one row per <b>stream</b>: final sigma, AR order, whiteness'],
        ['<font face="Courier">R.convergence()</font>',
         'one row per <b>iteration</b>: when the loop settled'],
    ], widths=[46 * mm, BODY_W - 46 * mm])]
    st += [callout(
        '<b>The two iteration readers are complementary, not duplicates.</b> '
        '<font face="Courier">convergence()</font> reduces across streams (one '
        'row per pass); <font face="Courier">summary()</font> reduces across '
        'passes (one row per stream). A loop that converges can still hide one '
        'stream that never whitened -- the per-pass view reports the worst '
        'p-value without naming its owner. Both need a history, so both raise '
        '<font face="Courier">lsqvce:Results:noHistory</font> on a single '
        'solve.', 'note')]
    st += [PageBreak()]
    return st


def part4_template():
    st = []
    st += [P('4&nbsp;&nbsp;Adapting the template: a step-by-step guide', 'h1')]
    st += [P('Copy <font face="Courier">template/lsqvce_template.m</font>, '
             'rename it, and work through the seven steps below. The template '
             'is written so that the parts you must replace fail loudly until '
             'you do: the two data stubs raise '
             '<font face="Courier">lsqvce_template:notImplemented</font> rather '
             'than returning plausible fake data.', 'body')]
    st += [callout(
        '<b>The template has two registers.</b> The loader loop is real, '
        'runnable code -- it builds one loader per observation group and you '
        'edit it only if your topology differs. <font face="Courier">myBuildDay'
        '</font> and <font face="Courier">myEpochCov</font> are stubs that '
        'error until written. Running the file as shipped therefore gets you '
        'loaders that fail on first use; that is the expected state, not a '
        'defect.', 'note')]

    # ---- 4.1 --------------------------------------------------------------
    st += [P('4.1&nbsp;&nbsp;What you write, and what the toolbox derives',
             'h2')]
    st += [workflow([
        [{'t': 'procOpt|your processing configuration', 'k': 'io'}],
        [{'t': 'observationGroups -> one loader per system', 'k': 'proc'}],
        [{'t': 'buildDay  AND  epochCov|the day, and its per-epoch D',
          'k': 'io'}],
        [{'t': 'makeLoader: validate, cache, wire the hooks', 'k': 'proc'}],
        [{'t': 'adjustVCE / adjustBlock', 'k': 'proc'}],
        [{'t': 'Results', 'k': 'io'}],
    ], width=BODY_W, title='Figure 4.1 -- what you supply')]
    st += [caption('Rounded boxes are yours; everything below makeLoader is the '
                   'toolbox\'s job. buildDay and epochCov are supplied '
                   'TOGETHER, not as alternatives: buildDay gives the design '
                   'and the observations, epochCov gives the per-epoch spatial '
                   'covariance D for the same rows. A scalar stream is the one '
                   'case where epochCov is omitted -- there is no D to supply '
                   'when d = 1.')]
    st += [table([
        ['You state', 'Where', 'The toolbox derives'],
        ['the common parametrisation',
         '<font face="Courier">commonParams</font> (one row per unknown)',
         '<font face="Courier">nc</font>, from the row count'],
        ['the local parametrisation',
         '<font face="Courier">localParams</font>, per system',
         '<font face="Courier">nl</font>, from the row count'],
        ['which streams exist',
         '<font face="Courier">observationTypes</font>',
         'the variance-component groups'],
        ['the day list',
         '<font face="Courier">procOpt.days</font>',
         'the batches to stream'],
    ], widths=[38 * mm, 58 * mm, BODY_W - 96 * mm])]
    st += [callout(
        '<b>State each value once.</b> There is no '
        '<font face="Courier">nc</font> argument on any entry point: the '
        'common-parameter count is stated where the parametrisation is defined '
        'and derived everywhere else, so the two cannot drift apart. The same '
        'discipline applies to your own configuration -- feed the toolbox '
        '<i>from</i> it rather than restating values.', 'note')]

    # ---- 4.2 --------------------------------------------------------------
    st += [P('4.2&nbsp;&nbsp;Step 1 -- declare the processing options', 'h2')]
    st += [P('<font face="Courier">procOpt</font> is your configuration. The '
             'toolbox never reads a field of it: it is captured when the loader '
             'is built, piped unchanged to your <font face="Courier">buildDay'
             '</font> as its fifth argument, and surfaced to plot hooks as '
             '<font face="Courier">ctx.procOpt</font>. Group it into nested '
             'structs; only the fields your builder needs have to exist.',
             'body')]
    st += [code("""procOpt.days          = 1:31;
procOpt.loggingLevel  = 1;              % 0 silent | 1 progress | 2 detail
procOpt.data          = struct('pathToLoad', ..., 'maximumDegree', 60);
procOpt.epochs        = struct('mEpoch', 300, 'dt', 5, 'arcsPerDay', 2);
procOpt.forces        = ...;            % your own sub-structs, opaque to lsqvce
procOpt.observationGroups = { groupC, groupD, groupE };""")]
    st += [callout(
        '<b>procOpt is immutable input, not shared state.</b> It is frozen when '
        '<font face="Courier">makeLoader</font> is called. A live cursor such '
        'as an iteration index would stay at its initial value while the loop '
        'iterates around it, and your builder would silently read the wrong '
        'pass forever. The engine tracks these itself: use '
        '<font face="Courier">ctx.pass</font> and '
        '<font face="Courier">ctx.sysID</font>.', 'danger')]
    st += [P('Values your pipeline needs before any toolbox object exists stay '
             'in <font face="Courier">procOpt</font>, and the toolbox is fed '
             '<i>from</i> them:', 'body')]
    st += [code("""R = lsqvce.adjustVCE(L, procOpt.days, ...
        'verbose', procOpt.loggingLevel) ...
      .groupComponentsBy('obsType') ...
      .run();""")]
    st += [P('Note that <font face="Courier">loggingLevel</font> maps directly '
             'onto <font face="Courier">verbose</font>, which accepts 0, 1 or '
             '2. <b>Level 3 is rejected</b>, so a pipeline that can emit 3 must '
             'clamp before passing it.', 'body')]

    # ---- 4.3 --------------------------------------------------------------
    st += [P('4.3&nbsp;&nbsp;Step 2 -- describe the observation groups', 'h2')]
    st += [P('Each entry of <font face="Courier">observationGroups</font> '
             'describes one <i>system</i>: which streams it carries and how it '
             'is identified. The template loops over them and calls '
             '<font face="Courier">makeLoader</font> once per group. This is '
             'where the topology is fixed.', 'body')]
    st += [callout(
        '<b>When two satellites must share one system.</b> Merge them whenever '
        'an observation\'s design couples both satellites\' <i>local</i> '
        'parameters -- an inter-satellite range or range-rate does exactly '
        'this. The block-arrow model has one local block per system '
        '(equation 2.10), so a coupling observation cannot span two of them: '
        'give both satellites\' streams the same group and stack their locals. '
        'Streams that are merely <i>correlated</i> do not need merging; that is '
        'what the variance components are for.', 'warn')]
    st += [callout(
        '<b><font face="Courier">systemId</font> must be numeric or a single '
        'character.</b> The accumulator indexes systems by '
        '<font face="Courier">[store.sysID]</font>, which flattens a multi-'
        'character string into separate elements. Carry a readable label in '
        '<font face="Courier">\'name\'</font> instead.', 'danger')]

    # ---- 4.4 --------------------------------------------------------------
    st += [P('4.4&nbsp;&nbsp;Step 3 -- write buildDay', 'h2')]
    st += [P('This is the one function nobody can write for you. It reads one '
             'day of one observation type and returns it in the shape the '
             'accumulator expects.', 'body')]
    st += [code("""function D = myBuildDay(day, obsType, common, local, procOpt)
    nc = size(common, 1);          % derive, never hardcode
    nl = size(local, 1);
    ...
    D.A      = ...;   % m x nc   common design matrix, WHOLE day
    D.B      = ...;   % m x nl   local design matrix (may be sparse)
    D.l      = ...;   % m x 1    observations
    D.t      = ...;   % m x 1    epoch timestamps [s]
    D.arcIdx = {1:n1, n1+1:n2};    % row indices, one cell per arc
    D.nAxis  = 3;                  % axes per epoch; omit for scalar streams
end""")]
    st += [table([
        ['Field', 'Size', 'Required?'],
        ['<font face="Courier">.A</font>',
         '<i>m</i> &#215; <i>nc</i>', 'yes'],
        ['<font face="Courier">.B</font>',
         '<i>m</i> &#215; <i>nl</i>', 'yes'],
        ['<font face="Courier">.l</font>',
         '<i>m</i> &#215; 1', 'yes'],
        ['<font face="Courier">.t</font>',
         '<i>m</i> &#215; 1', 'yes'],
        ['<font face="Courier">.arcIdx</font>',
         '1 &#215; nArc cell',
         'optional -- omitted means the whole day is ONE arc'],
        ['<font face="Courier">.nAxis</font>',
         'scalar',
         'multi-axis streams only'],
    ], widths=[26 * mm, 30 * mm, BODY_W - 56 * mm])]
    st += [callout(
        '<b>It must be deterministic.</b> The accumulator reads every day a '
        '<i>second</i> time for the residual pass and must get identical data '
        'back. Seed any random number generator from '
        '(<font face="Courier">day</font>, '
        '<font face="Courier">obsType</font>), never from the clock -- and '
        'beware a file being rewritten while the run is in progress.',
        'danger')]
    st += [P('<b>You return a whole day, not an arc.</b> The accumulator slices '
             'arcs out of what you return, using '
             '<font face="Courier">.arcIdx</font>. Omitting it makes the whole '
             'day one arc, which is right for a single continuous pass and '
             'wrong for a day containing several -- one local block would then '
             'be asked to describe several disjoint states.', 'body')]

    # ---- 4.5 --------------------------------------------------------------
    st += [P('4.5&nbsp;&nbsp;Step 4 -- the a-priori covariance', 'h2')]
    st += [P('One covariance spec per observation type. If you are going to run '
             '<b>adjustVCE</b>, a deliberately plain start is honest: the loop '
             'exists to discover the truth, and beginning with unit-variance '
             'white noise lets you watch it work.', 'body')]
    st += [code("""cov.pos = struct('type','separable', 'ar',1, 'sigma',1, 'd',3);
cov.rr  = struct('type','diag',      'sigma',1);""")]
    st += [P('Scalar streams use <font face="Courier">\'diag\'</font>; '
             'multi-axis streams use <font face="Courier">\'separable\'</font>, '
             'the model of equation 2.12, and must also supply the per-epoch '
             '<b>D</b> described next.', 'body')]
    st += [callout(
        '<b><font face="Courier">\'separable\'</font> here is a covariance '
        'MODEL, not an engine name.</b> It is one of eight values '
        '<font face="Courier">type</font> can take -- alongside '
        '<font face="Courier">diag</font>, <font face="Courier">acf</font>, '
        '<font face="Courier">psd</font>, <font face="Courier">arma</font>, '
        '<font face="Courier">var</font>, <font face="Courier">toeplitz</font> '
        'and <font face="Courier">custom</font> -- and it names the structure '
        'of equation 2.12, C = &#963;<sup>2</sup>D(R<sub>t</sub>&#8855;I<sub>d'
        '</sub>)D. The <i>class</i> once called '
        '<font face="Courier">SeparableProblem</font> was renamed '
        '<font face="Courier">adjustVCE</font> in 2026; the covariance type was '
        'not touched, because it describes mathematics rather than an entry '
        'point.', 'note')]

    # ---- 4.6 --------------------------------------------------------------
    st += [P('4.6&nbsp;&nbsp;Step 5 -- epochCov, for multi-axis streams', 'h2')]
    st += [P('Only multi-axis streams need this. It supplies the known spatial '
             'half <b>D</b> of the separable model, per epoch, so the temporal '
             'law can be recovered by exact de-scaling rather than inferred '
             'from whitened residuals.', 'body')]
    st += [code("""function dd = myEpochCov(day, obsType, procOpt)
    dd.Cep = ...;                  % nEpochDay x d   OR  nEpochDay x d x d
    dd.arc(1).epochs = 1:300;      % which day-epochs this arc covers
    dd.arc(1).rows   = 1:900;      % and which ROWS  (R = E*d)
end""")]
    st += [P('The two accepted layouts of <font face="Courier">Cep</font> are '
             'given in table 2.1: <font face="Courier">nEpoch &#215; d</font> '
             'for per-axis standard deviations, or '
             '<font face="Courier">nEpoch &#215; d &#215; d</font> for the full '
             'per-epoch covariance. Anything else raises '
             '<font face="Courier">lsqvce:covModel:badCep</font>.', 'body')]
    st += [callout(
        'The <font face="Courier">.arc(i).epochs</font> and '
        '<font face="Courier">.arc(i).rows</font> fields exist so the two '
        'indexings can be <i>checked against each other</i> rather than assumed '
        'to agree. Row <i>k</i> belongs to epoch '
        '&#8968;<i>k</i>/<i>d</i>&#8969; (equation 2.13); if your rows and your '
        'epochs disagree, this is where it surfaces.', 'note')]

    # ---- 4.7 --------------------------------------------------------------
    st += [P('4.7&nbsp;&nbsp;Step 6 -- run it, and read the result', 'h2')]
    st += [code("""L = myCase(procOpt);                    % the template returns loaders
R = lsqvce.adjustVCE(L, procOpt.days, 'verbose', procOpt.loggingLevel) ...
      .groupComponentsBy('obsType') ...
      .run();

R.report();          % the verdict -- read top-down, stop at the first red
R.summary();         % per stream: sigma, AR order, whiteness
R.convergence();     % per iteration: did it settle, and when""")]
    st += [P('Read <font face="Courier">R.report()</font> first and stop at the '
             'first red. A Level&nbsp;0 red is a data or weighting bug, not a '
             'blunder -- fix it before reading anything else. Only once the '
             'verdict is clean are the variance factors worth interpreting.',
             'body')]

    # ---- 4.8 --------------------------------------------------------------
    st += [P('4.8&nbsp;&nbsp;Step 7 -- optional plot hooks', 'h2')]
    st += [P('A hook is a handle <i>you</i> write, '
             '<font face="Courier">@(ctx) -&gt; figureHandle</font>. It only '
             '<i>builds</i> a figure; the engine owns showing, exporting and '
             'closing it. A hook that throws is logged and skipped, so a broken '
             'plot never breaks a solve.', 'body')]
    st += [table([
        ['Stage', 'Fires', 'ctx additionally carries'],
        ['<font face="Courier">postSolve</font>',
         'after EACH solve, before the residual pass',
         'the current estimate'],
        ['<font face="Courier">postModel</font>',
         'ONCE, after the stochastic model is fitted',
         '<font face="Courier">.psd .psdFreq .acv .arOrder</font>'],
    ], widths=[26 * mm, 62 * mm, BODY_W - 88 * mm])]
    st += [callout(
        '<b>Hooks attach at construction.</b> Passing a bare handle to '
        '<font face="Courier">makeLoader</font> broadcasts it to every '
        'observation type. Assigning '
        '<font face="Courier">L.plotPostSolve</font> <i>after</i> the loader '
        'exists requires a <font face="Courier">containers.Map</font> keyed by '
        'obsType -- the engine skips anything else <b>silently</b>.', 'warn')]

    # ---- 4.9 --------------------------------------------------------------
    st += [P('4.9&nbsp;&nbsp;Four ways to get it wrong', 'h2')]
    st += [P('Two of these the toolbox catches for you. Two it cannot, and '
             'those are the ones worth memorising.', 'body')]

    st += [P('I -- a builder that is not deterministic', 'h3')]
    st += [code("""% WRONG -- fresh noise on every call
D.l = readObs(day) + 0.01*randn(m,1);

% RIGHT -- seeded from (day, obsType), so both passes agree
rng(day*1000 + double(obsType(1)));
D.l = readObs(day) + 0.01*randn(m,1);""")]
    st += [P('<b>How you notice:</b> you do not, from the numbers alone. '
             'Residuals are attributed to the wrong observations while every '
             'statistic still looks plausible.', 'body')]

    st += [P('II -- axis-major row order', 'h3')]
    st += [code("""% WRONG -- all x, then all y, then all z
D.A = [Ax; Ay; Az];

% RIGHT -- epoch-major: [x1 y1 z1 | x2 y2 z2 | ...]
D.A = reshape([Ax'; Ay'; Az'], nc, [])';""")]
    st += [callout(
        '<b>This is the silent one.</b> The temporal model strides by '
        '<font face="Courier">nAxis</font> and equation 2.13 pairs each row '
        'with an epoch. Stored axis-major, every row is paired with the wrong '
        'epoch: a correlation length of 60 epochs is read as 180. The fit '
        'succeeds, the AR order looks sensible, the variance factors converge '
        '-- and the stochastic model is wrong. <b>Check row order first '
        'whenever a multi-axis stream will not whiten.</b>', 'danger')]

    st += [P('III -- variances where weights are expected', 'h3')]
    st += [P('The weight <b>P</b> is the <i>inverse</i> of the covariance. '
             'Passing &#963;<sup>2</sup> instead of '
             '1/&#963;<sup>2</sup> silently inverts the weighting, so your best '
             'observations are down-weighted. It is dimensionally valid, so '
             'nothing can detect it for you.', 'body')]

    st += [P('IV -- mismatched row counts', 'h3')]
    st += [code("""D.A = zeros(m, nc);   D.l = zeros(m+1, 1);   % A and l disagree""")]
    st += [P('<b>How you notice:</b> immediately, and by name. '
             '<font face="Courier">lsqvce:makeLoader:dayShape</font> reports '
             'the day and obsType; its sibling '
             '<font face="Courier">lsqvce:makeLoader:badDay</font> fires when a '
             'required field is missing. These two are why a shape error never '
             'reaches the solver.', 'body')]

    st += [P('A checklist for a stream that will not whiten', 'h3')]
    st += [N('<b>Row order</b> -- epoch-major? (Failure II: silent, and the '
             'most common cause.)', 1)]
    st += [N('<b>Arcs</b> -- is <font face="Courier">.arcIdx</font> really '
             'splitting the day, or did you inherit the whole-day default?', 2)]
    st += [N('<b>AR order</b> -- is the fit pinned at '
             '<font face="Courier">MaxOrder</font>? Read '
             '<font face="Courier">R.summary()</font>; a stream at the ceiling '
             'needs a higher one, or ARMA.', 3)]
    st += [N('<b>Grouping</b> -- does each group have enough redundancy? '
             'Oscillating variance factors in '
             '<font face="Courier">R.convergence()</font> mean it does not; '
             'coarsen <font face="Courier">groupComponentsBy</font>.', 4)]
    st += [N('<b>Determinism</b> -- re-run and compare. Different numbers from '
             'the same input means failure I.', 5)]
    st += [PageBreak()]
    return st
