function [R, s] = demo_priors(varargin)
%DEMO_PRIORS  Regularisation written in the columns LSQVCE.LAYOUT announces.
%
%   [R, s] = DEMO_PRIORS()
%   [R, s] = DEMO_PRIORS('verbose', true)
%
%   THE PROBLEM THIS EXAMPLE EXISTS FOR.  In a lsqvce.TaskDefinition you state
%   TEMPLATES -- task.addCommon and sys.addLocal, each with a Scope -- and the
%   toolbox expands them over the calendar, so you do not know -- and must
%   never guess -- which column of A or of B a given unknown ended up in. But
%   regularisation is written in exactly those columns: a prior group is an
%   operator L on the parameter vector. task.layout() (LSQVCE.LAYOUT) is the
%   interface between the two. It announces the layout BEFORE fit runs, so you
%   can write L in it:
%
%       lay = task.layout();
%       c   = lay.colCommon(row, day);           % a column of A
%       c   = lay.systems(k).col(row, day, arc); % a column of that system's B
%       lay.common(j).name                       % the name behind column j of A
%
%   ROW is the TEMPLATE row: the position of the name in addCommon (for A) or
%   in the system's addLocal calls (for B), in the order they were made.
%
%   THE CASE (a calibration problem that is genuinely rank-deficient).  Four
%   harmonic signal coefficients are observed by one instrument over three
%   days, two arcs per day. The instrument carries a per-DAY offset, and each
%   ARC carries its own bias and drift. Offset and arc biases are separable
%   only up to a constant: add c to a day's offset, subtract c from both of
%   that day's arc biases, and every observation is unchanged. The normal
%   matrix therefore has ONE null direction PER DAY -- it is exactly singular,
%   and no amount of data will fix it. Something outside the data must choose
%   the datum. That something is a prior group, task.addPrior.
%
%   WHY THE LAYOUT IS UNAVOIDABLE HERE.  The datum this case needs is "the
%   offsets average to zero and vary smoothly": one row of ones over the
%   day-offset columns, and one differencing row per consecutive pair. Neither
%   is a template row, so neither can be named by the pattern form (On=, a
%   regular expression on the names) -- they exist only as columns, and only
%   the layout knows which. They are given as L=.
%
%   A TASK ONLY GROWS.  A prior group is added, never removed. Every variant
%   compared below (no priors, the datum alone, the pattern form, another
%   datum, the noise run) is therefore a task of its own: calibTask declares
%   the instrument, the add-calls of each step give it its groups.
%
%   WHAT IT SHOWS, IN ORDER
%     1. the layout census: which template row lives in which column, per day;
%     2. what fit does with the singular problem when no prior is declared
%        (measured, and quoted by error identifier);
%     3. the datum: two common groups (task.addPrior(..., L=...)) written
%        through lay.colCommon -- one tight, one loose;
%     4. a local group (sys.addPrior(..., L=...)) written through
%        lay.systems(k).col, and the SAME group written as a pattern on the
%        names (On="^drift$"): both forms give the same solution;
%     5. datum invariance, measured: the observable part of the solution (the
%        four signal coefficients) does not move when the datum does;
%     6. the prior group is DATA: under Estimate="noise" it rides as its own
%        stream prior_<name> and its sigma is estimated like any other level.
%
%   INPUT  (Name-Value, optional)
%     verbose   1x1 logical   print the engine log too (default false)
%   OUTPUT
%     R         1x1 lsqvce.Results   the regularised solve of steps 3 and 4
%     s         1x1 struct           the measured summary:
%       .nc .nl             the DERIVED counts
%       .nullDim            null directions of the unregularised design
%       .unregularised      what fit did without priors: an error identifier
%       .sigma0             sigma0 of the regularised solve
%       .zSignal            1x4  (estimate - truth) / sigma of the signal
%       .nPseudo            pseudo-observations the groups added
%       .formGap            max|x(L form) - x(On form)|
%       .datumGap           max|x_signal(datum A) - x_signal(datum B)|
%       .priorLevel         the estimated level of the local group's stream
%
%   ACCURACY  (MEASURED 2026-09-07 on R2026a. The fixture noise is drawn from
%   a stream seeded on the day, so every number below repeats -- except the
%   Estimate="noise" level, whose redundancies come from unseeded Hutchinson
%   probes and which read 0.747, 0.753, 0.756 and 0.757 over four runs.)
%     nc = 7 and nl = 13, both DERIVED; the design has 3 null directions, one
%       per day, and rank() confirms it on the densely assembled [A B].
%     Without priors fit REFUSES: lsqvce:qr:singularR (a pivot below
%       QRSystem.RankTol; until 6.8 the normal-form fold dropped the null
%       directions and the refusal read lsqvce:qr:rankDeficient).
%     Regularised: sigma0 = 1.017; the four signal coefficients sit at
%       z = +0.16, +0.00, -0.21, +0.04 of their formal errors; the three
%       groups add 9 pseudo-observations (1 datum + 2 smoothness + 6 drift).
%     The L form and the On form agree BITWISE (max|dx| = 0).
%     Moving the datum (the offsets average to 0.5, not 0) moves the day
%       offsets by 0.289 = 0.5/sqrt(3), exactly what the datum row prescribes,
%       and the signal by 1.6e-11 -- the observable part is untouched.
%     Under Estimate="noise" the streams are x1|prior_drift and x1|rng; the
%       drift group's level 0.75 x its declared sigma 0.05 = 0.037 against the
%       planted drift rms 0.031. With r = 6 rows the relative standard error of
%       a variance component is 0.5*sqrt(2/r) = 0.29, so that gap is noise.
%
%   ERRORS  none of its own. A prior group is checked where it is declared:
%     addPrior refuses a group without Sigma
%     (lsqvce:TaskDefinition:missingArgument), with both On and L or neither
%     (lsqvce:TaskDefinition:priorTarget), and an Estimate other than "none" or
%     "scale" (MATLAB:validators:mustBeMember). What only the columns can tell
%     is refused when fit runs: an L of the wrong width (lsqvce:fit:badPrior
%     for a common group, lsqvce:makeLoader:badPrior for a local one) and a
%     LOCAL row that couples two days (lsqvce:makeLoader:badPrior). A COMMON
%     group may couple days: the smoothness group here does, and must, since
%     that is what a smoothness statement is.
%
%   See also LSQVCE.TASKDEFINITION, LSQVCE.SYSTEMDEFINITION, LSQVCE.LAYOUT,
%   LSQVCE.FIT, LSQVCE.OPTIONS, DEMO_NOISEPOLICY, DEMO_WORKFLOW
%
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: a lsqvce.TaskDefinition; each prior variant is its own task. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-17 (CEST, UTC+02:00), WP-150: the drift group is declared
%   estimate "none" for the single solves and switched to "scale" for the one
%   run that estimates it. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-112 (R7-16): the refusal
%   without priors is lsqvce:qr:singularR, the rank decision on the pivots.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Created 2026-09-07 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5 -- R6 WP-96: lsqvce.layout had no worked
%   example (HANDOVER 12.19, carried; 12.20).

    p = inputParser;
    p.addParameter('verbose', false);
    p.parse(varargin{:});
    verbose = p.Results.verbose;

    cfg   = priorCfg();
    task0 = calibTask(cfg);                         % the instrument, no prior yet

    fprintf(['\n=== Regularisation in the announced layout ===\n' ...
             '    %d harmonic signal coefficients + one offset per day;\n' ...
             '    %d days x %d arcs, each arc its own bias and drift.\n' ...
             '    Offset and arc biases are separable only up to a constant.\n\n'], ...
            cfg.nSig, numel(cfg.days), cfg.arcs);

    % ---- 1. THE LAYOUT ---------------------------------------------------
    % Ask BEFORE the run. Nothing here is your bookkeeping: nc and nl are
    % DERIVED from the expansion of the templates, and every column carries
    % the name it was declared with.
    lay = task0.layout();
    s.nc = lay.nc;
    s.nl = lay.systems(1).nl;
    fprintf('lsqvce.layout says: nc = %d, nl = %d\n', s.nc, s.nl);
    fprintf('  A: ');
    for j = 1:lay.nc
        ec = lay.common(j);
        if isempty(ec.day)
            fprintf('%d=%s ', j, ec.name);
        else
            fprintf('%d=%s(d%d) ', j, ec.name, ec.day);
        end
    end
    fprintf('\n  B: ');
    loc = lay.systems(1).local;
    for j = 1:numel(loc)
        el = loc(j);
        if isempty(el.day)
            fprintf('%d=%s ', j, el.name);
        else
            fprintf('%d=%s(d%d,a%d) ', j, el.name, el.day, el.arc);
        end
    end
    fprintf('\n');

    % The columns the datum needs. They cannot be named by a pattern: they are
    % ONE template row instantiated once per day, and the operator couples them.
    dayCols = arrayfun(@(d) lay.colCommon(cfg.rowOffset, d), cfg.days);
    fprintf('the day-offset columns of A are %s\n', mat2str(dayCols));

    % ---- 2. WITHOUT A PRIOR ----------------------------------------------
    % MEASURED, not assumed: what a singular normal matrix actually does.
    s.nullDim = nullDimension(cfg, lay);
    fprintf('\nnull directions of the unregularised design: %d (one per day)\n', s.nullDim);
    try
        lsqvce.fit(task0, lsqvce.Options(Verbose=false));
        s.unregularised = 'accepted (no refusal)';
    catch ME
        s.unregularised = ME.identifier;
    end
    fprintf('fit without priors: %s\n', s.unregularised);

    % ---- 3. THE DATUM, WRITTEN IN THE LAYOUT -----------------------------
    % Two groups, because they are two statements with two strengths:
    %   "datum"  -- the offsets average to zero.        TIGHT: it is a datum.
    %   "smooth" -- consecutive offsets differ little.  LOOSE: it is a belief.
    % addDatum gives them to a task: task.addPrior(name, L=..., Value=..., Sigma=...).
    nD = numel(cfg.days);
    Lmean = zeros(1, lay.nc);
    Lmean(dayCols) = 1/sqrt(nD);
    Ldiff = zeros(nD-1, lay.nc);
    for k = 1:nD-1
        Ldiff(k, dayCols(k))   = -1;
        Ldiff(k, dayCols(k+1)) =  1;
    end

    % ---- 4. A LOCAL GROUP, AND THE SAME GROUP BY PATTERN -----------------
    % Damp the arc drifts. Written twice: through the layout (works for ANY
    % operator) and as a pattern on the NAMES (works when the group IS a whole
    % template row, as here). Both must give one solution.
    unit = @(c) double((1:s.nl) == c);
    Ld = zeros(numel(cfg.days) * cfg.arcs, s.nl);
    r  = 0;
    for d = cfg.days
        for a = 1:cfg.arcs
            r = r + 1;
            Ld(r, :) = unit(lay.systems(1).col(cfg.rowDrift, d, a));
        end
    end
    % Estimate "none": these single solves use the declared sigma, and say so
    % by the declaration (section 6 declares the group "scale" for the noise
    % estimator, the one run that estimates it -- WP-150: a "scale" group in a
    % single solve warns that it is held)
    [task, sys] = calibTask(cfg);
    addDatum(task, cfg, Lmean, Ldiff, 0);
    sys.addPrior("drift", L=Ld, Value=zeros(size(Ld, 1), 1), Sigma=cfg.driftSigma, Estimate="none");
    R  = lsqvce.fit(task, lsqvce.Options(Verbose=verbose));
    R0 = fitWithoutGroups(cfg, Lmean, Ldiff);

    [taskPat, sysPat] = calibTask(cfg);
    addDatum(taskPat, cfg, Lmean, Ldiff, 0);
    sysPat.addPrior("drift", On="^drift$", Value=0, Sigma=cfg.driftSigma, Estimate="none");
    Rp = lsqvce.fit(taskPat, lsqvce.Options(Verbose=false));

    s.formGap = max(abs(R.x - Rp.x));
    s.nPseudo = size(Lmean, 1) + size(Ldiff, 1) + size(Ld, 1);
    assert(R.diag.n - R0.diag.n == size(Ld, 1), ...
           'the drift group adds one pseudo-observation per row of L');
    s.sigma0  = R.sigma0;
    se = R.coefficients().se;                       % se already carries sigma0
    s.zSignal = ((R.x(1:cfg.nSig) - cfg.xSig(:)) ./ se(1:cfg.nSig)).';

    fprintf('\nregularised solve: sigma0 = %.3f, %d pseudo-observations from the groups\n', ...
            s.sigma0, s.nPseudo);
    fprintf('signal coefficients, (estimate - truth)/sigma: %s\n', ...
            num2str(s.zSignal, '%+.2f  '));
    fprintf('L form vs `On` form: max|dx| = %.1e (round-off)\n', s.formGap);

    % ---- 5. DATUM INVARIANCE, MEASURED -----------------------------------
    % The offsets are NOT observable, so the datum DECIDES them; the signal IS
    % observable, so the datum must NOT touch it. That is a claim about the
    % algebra, and a claim gets measured: move the datum itself -- state that
    % the offsets average to 0.5 instead of 0 -- and look at both halves.
    [taskB, sysB] = calibTask(cfg);
    addDatum(taskB, cfg, Lmean, Ldiff, 0.5);
    sysB.addPrior("drift", L=Ld, Value=zeros(size(Ld, 1), 1), Sigma=cfg.driftSigma, Estimate="none");
    RB = lsqvce.fit(taskB, lsqvce.Options(Verbose=false));
    s.datumGap = max(abs(R.x(1:cfg.nSig) - RB.x(1:cfg.nSig)));
    offsetMove = max(abs(R.x(dayCols) - RB.x(dayCols)));
    fprintf('\na different datum: the signal moves %.1e, the day offsets move %.3f\n', ...
            s.datumGap, offsetMove);

    % ---- 6. A PRIOR GROUP IS DATA ----------------------------------------
    % With Estimate="scale" the group is a stream of its own and its sigma is
    % the level the joint estimator computes for every stream (Koch-Kusche).
    % Levels only here (the instrument stream asks for "scale" too): six arcs
    % of 150 epochs are too thin to identify a spectral SHAPE, and a demo must
    % not ask for what its data cannot give -- see DEMO_NOISEPOLICY for the
    % per-stream policy itself.
    [taskN, sysN] = calibTask(cfg, "scale");
    addDatum(taskN, cfg, Lmean, Ldiff, 0);
    sysN.addPrior("drift", L=Ld, Value=zeros(size(Ld, 1), 1), Sigma=cfg.driftSigma, Estimate="scale");
    Rn = lsqvce.fit(taskN, lsqvce.Options(Estimate="noise", Pool="stream", ...
                                          ArcFactors=false, MaxPasses=4, ...
                                          Verbose=verbose));
    k = find(strcmp(Rn.noise.stream, 'x1|prior_drift'), 1);
    s.priorLevel = Rn.noise.model{k}.sigma;
    fprintf('\nEstimate="noise": the streams are %s\n', strjoin(Rn.history.streams, ', '));
    fprintf('the drift group''s estimated level: %.2f (declared 1)\n', s.priorLevel);
    % A level is only as good as the redundancy behind it: say what it implies,
    % what the truth was, and how well 6 rows can possibly pin a variance.
    fprintf(['  -> implied drift sigma %.4f against the planted rms %.4f;\n' ...
             '     r = %d rows, so the relative standard error is 0.5*sqrt(2/r) = %.2f\n'], ...
            cfg.driftSigma * s.priorLevel, sqrt(mean(cfg.drift(:).^2)), ...
            size(Ld, 1), 0.5*sqrt(2/size(Ld, 1)));
    fprintf(['\nTHE RULE (the API reference, chapter Priors): the parameters and their\n' ...
             'columns are yours; the regularisation MECHANISM is the toolbox''s; the\n' ...
             'interface between them is an operator in your columns, announced by\n' ...
             'lsqvce.layout and checked by width.\n']);
end

%==========================================================================
function cfg = priorCfg()
%PRIORCFG  The case: sizes, the planted truth, and the template rows the demo
%   addresses by index (nothing here is read by the toolbox). A template row
%   is the position of a name in the order calibTask declares it.
    cfg.nSig        = 4;                % harmonic signal coefficients (Scope "common")
    cfg.sigNames    = compose("sig%d", 1:cfg.nSig);
    cfg.days        = 1:3;
    cfg.arcs        = 2;
    cfg.mArc        = 150;              % epochs per arc
    cfg.dt          = 30;               % s
    cfg.sigma       = 0.02;             % the instrument's true noise (known here)
    cfg.rowOffset   = 5;                % common template row of "offset" (after sig1..sig4)
    cfg.rowBias     = 1;                % local template row of "bias"
    cfg.rowDrift    = 2;                % local template row of "drift"
    cfg.rowTemp     = 3;                % local template row of "tempScale" (Scope "system")
    cfg.datumSigma  = 1e-3;             % the datum row: tight, but see below
    cfg.smoothSigma = 0.20;
    cfg.driftSigma  = 0.05;
    cfg.xSig        = [0.80 -0.45 0.30 0.15];               % the truth (here only)
    cfg.offset      = [0.30 -0.20 0.50];                    % per day; NOT observable
    cfg.bias        = [0.10 -0.05; 0.02 0.08; -0.06 0.03];  % per (day, arc)
    cfg.drift       = [0.04 -0.03; 0.01 0.05; -0.02 0.02];  % per (day, arc), per hour
    cfg.tempScale   = 0.25;             % the "system" row: one for the whole run
end

function [task, sys] = calibTask(cfg, rngEstimate)
%CALIBTASK  The declaration of the instrument, without any prior group.
%   The order of the add-calls IS the template-row order priorCfg addresses:
%   commons sig1..sig4 (Scope "common"), then offset (Scope "day"); locals
%   bias, drift (Scope "arc"), then tempScale (Scope "system"). RNGESTIMATE
%   (optional, "none" | "scale" | "law") is the stream's noise policy; without
%   it the stream has none of its own (the default).
    task = lsqvce.TaskDefinition(Days=cfg.days);
    task.addCommon(cfg.sigNames);
    task.addCommon("offset", Scope="day");
    sys = task.addSystem("calib", Arcs=cfg.arcs, Id=1);
    sys.addLocal(["bias" "drift"], Scope="arc");
    sys.addLocal("tempScale", Scope="system");
    build = @(day, ctx) buildCalib(cfg, day, ctx);
    if nargin < 2
        sys.addStream("rng", Build=build, Noise=lsqvce.Noise.white(cfg.sigma));
    else
        sys.addStream("rng", Build=build, Noise=lsqvce.Noise.white(cfg.sigma), Estimate=rngEstimate);
    end
end

function addDatum(task, cfg, Lmean, Ldiff, datum)
%ADDDATUM  The two common groups of step 3, written in the layout's columns:
%   "datum" (the day offsets average to DATUM, tight) and "smooth"
%   (consecutive offsets differ little, loose). TASK is a handle: it grows.
    task.addPrior("datum", L=Lmean, Value=datum, Sigma=cfg.datumSigma, Estimate="none");
    task.addPrior("smooth", L=Ldiff, Value=zeros(size(Ldiff, 1), 1), Sigma=cfg.smoothSigma, ...
                  Estimate="none");
end

function out = buildCalib(cfg, day, ctx)
%BUILDCALIB  The stream's Build: one day of the instrument, two arcs. Every
%   column of A and B is asked of ctx.cols by name; priorArcs places them.
    col.nc     = ctx.NumCommon;
    col.nl     = ctx.NumLocal;
    col.sig    = ctx.cols(cfg.sigNames);
    col.offset = ctx.cols("offset");                % this day's offset
    col.temp   = ctx.cols("tempScale");             % the "system" row
    col.bias   = zeros(1, cfg.arcs);
    col.drift  = zeros(1, cfg.arcs);
    for a = 1:cfg.arcs
        col.bias(a)  = ctx.cols("bias", Arc=ctx.ArcIds(a));
        col.drift(a) = ctx.cols("drift", Arc=ctx.ArcIds(a));
    end
    D = priorArcs(cfg, day, col);
    out = ctx.newDay();
    for a = 1:cfg.arcs
        out.addArc(D.A{a}, D.l{a}, D.t{a}, B=D.B{a});   % the arc's B as one matrix
    end
end

function D = priorArcs(cfg, day, col)
%PRIORARCS  One day of the instrument, placed in the columns COL names (from
%   ctx.cols in Build, from the layout in nullDimension). DETERMINISTIC: fit
%   reads every day twice per pass, so the noise comes from a stream seeded on
%   the day. D.A, D.B, D.l, D.t are 1 x cfg.arcs cells, one per arc.
    m = cfg.arcs * cfg.mArc;
    t = (day - 1)*86400 + (0:m-1).' * cfg.dt;
    w = 2*pi / (0.5 * 86400);
    Ssig = [sin(w*t), cos(w*t), sin(2*w*t), cos(2*w*t)];
    temp = cos(2*pi * t / (1.7 * 86400));           % a known auxiliary series

    A = zeros(m, col.nc);
    A(:, col.sig) = Ssig;
    A(:, col.offset) = 1;                           % this day's offset
    B = zeros(m, col.nl);
    B(:, col.temp) = temp;                          % the "system" row
    l = Ssig * cfg.xSig(:) + cfg.offset(day) + cfg.tempScale * temp;
    rs  = RandStream('mt19937ar', 'Seed', 4200 + day);
    l = l + cfg.sigma * randn(rs, m, 1);
    D = struct('A', {cell(1, cfg.arcs)}, 'B', {cell(1, cfg.arcs)}, ...
               'l', {cell(1, cfg.arcs)}, 't', {cell(1, cfg.arcs)});
    for a = 1:cfg.arcs
        rows = (a-1)*cfg.mArc + (1:cfg.mArc);
        tau  = (t(rows) - t(rows(1))) / 3600;       % hours since the arc start
        Ba = B(rows, :);
        Ba(:, col.bias(a))  = 1;
        Ba(:, col.drift(a)) = tau;
        D.A{a} = A(rows, :);
        D.B{a} = Ba;
        D.l{a} = l(rows) + cfg.bias(day, a) + cfg.drift(day, a) * tau;
        D.t{a} = t(rows);
    end
end

function n = nullDimension(cfg, lay)
%NULLDIMENSION  The null space of the FULL (common + local) design, assembled
%   densely HERE only to show what the priors repair. The toolbox never forms
%   it: it eliminates the locals arc by arc.
    M = zeros(0, lay.nc + lay.systems(1).nl);
    for d = cfg.days
        D = priorArcs(cfg, d, layoutColumns(cfg, lay, d));
        M = [M; vertcat(D.A{:}), vertcat(D.B{:})]; %#ok<AGROW>   the cells are the arcs
    end
    n = size(M, 2) - rank(M);
end

function col = layoutColumns(cfg, lay, day)
%LAYOUTCOLUMNS  The columns of DAY as the layout announces them -- the ones
%   ctx.cols hands Build -- so the dense assembly above sees the columns fit
%   will. In your own Build you NEVER do this: ask ctx.cols. Here it keeps the
%   rank check honest.
    col.nc     = lay.nc;
    col.nl     = lay.systems(1).nl;
    col.sig    = arrayfun(@(row) lay.colCommon(row, []), 1:cfg.nSig);
    col.offset = lay.colCommon(cfg.rowOffset, day);
    col.temp   = lay.systems(1).col(cfg.rowTemp, [], []);
    col.bias   = arrayfun(@(a) lay.systems(1).col(cfg.rowBias, day, a), 1:cfg.arcs);
    col.drift  = arrayfun(@(a) lay.systems(1).col(cfg.rowDrift, day, a), 1:cfg.arcs);
end

function R = fitWithoutGroups(cfg, Lmean, Ldiff)
%FITWITHOUTGROUPS  The same solve with the DATUM groups only, so the pseudo-
%   observations of the drift group can be reported as a delta.
    task = calibTask(cfg);
    addDatum(task, cfg, Lmean, Ldiff, 0);
    R = lsqvce.fit(task, lsqvce.Options(Verbose=false));
end
