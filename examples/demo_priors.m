function [R, s] = demo_priors(varargin)
%DEMO_PRIORS  Regularisation written in the columns LSQVCE.LAYOUT announces.
%
%   [R, s] = DEMO_PRIORS()
%   [R, s] = DEMO_PRIORS('verbose', true)
%
%   THE PROBLEM THIS EXAMPLE EXISTS FOR.  In a scoped declaration you state
%   TEMPLATES and the toolbox expands them over the calendar, so you do not
%   know -- and must never guess -- which column of A or of B a given unknown
%   ended up in. But regularisation is written in exactly those columns: a
%   prior group is an operator L on the parameter vector. LSQVCE.LAYOUT is the
%   interface between the two. It announces the layout BEFORE fit runs, so you
%   can write L in it:
%
%       lay = lsqvce.layout(data);
%       c   = lay.colCommon(row, day);           % a column of A
%       c   = lay.systems(k).col(row, day, arc); % a column of that system's B
%
%   THE CASE (a calibration problem that is genuinely rank-deficient).  Four
%   harmonic signal coefficients are observed by one instrument over three
%   days, two arcs per day. The instrument carries a per-DAY offset, and each
%   ARC carries its own bias and drift. Offset and arc biases are separable
%   only up to a constant: add c to a day's offset, subtract c from both of
%   that day's arc biases, and every observation is unchanged. The normal
%   matrix therefore has ONE null direction PER DAY -- it is exactly singular,
%   and no amount of data will fix it. Something outside the data must choose
%   the datum. That something is data.priors.
%
%   WHY THE LAYOUT IS UNAVOIDABLE HERE.  The datum this case needs is "the
%   offsets average to zero and vary smoothly": one row of ones over the
%   day-offset columns, and one differencing row per consecutive pair. Neither
%   is a template row, so neither can be named by the `group` pattern form --
%   they exist only as columns, and only LSQVCE.LAYOUT knows which.
%
%   WHAT IT SHOWS, IN ORDER
%     1. the layout census: which template row lives in which column, per day;
%     2. what fit does with the singular problem when no prior is declared
%        (measured, and quoted by error identifier);
%     3. the datum: two common groups (data.priors) written through
%        lay.colCommon -- one tight, one loose;
%     4. a local group (data.systems(k).priors) written through
%        lay.systems(k).col, and the SAME group written as a `group` regular
%        expression: both forms give the same solution;
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
%       .formGap            max|x(L form) - x(group form)|
%       .datumGap           max|x_signal(datum A) - x_signal(datum B)|
%       .priorLevel         the estimated level of the local group's stream
%
%   ACCURACY  (MEASURED 2026-09-07 on R2026a. The fixture noise is drawn from
%   a stream seeded on the day, so every number below repeats -- except the
%   Estimate="noise" level, whose redundancies come from unseeded Hutchinson
%   probes and which read 0.747, 0.753, 0.756 and 0.757 over four runs.)
%     nc = 7 and nl = 13, both DERIVED; the design has 3 null directions, one
%       per day, and rank() confirms it on the densely assembled [A B].
%     Without priors fit REFUSES: lsqvce:qr:rankDeficient.
%     Regularised: sigma0 = 1.017; the four signal coefficients sit at
%       z = +0.16, +0.00, -0.21, +0.04 of their formal errors; the three
%       groups add 9 pseudo-observations (1 datum + 2 smoothness + 6 drift).
%     The L form and the `group` form agree BITWISE (max|dx| = 0).
%     Moving the datum (the offsets average to 0.5, not 0) moves the day
%       offsets by 0.289 = 0.5/sqrt(3), exactly what the datum row prescribes,
%       and the signal by 1.6e-11 -- the observable part is untouched.
%     Under Estimate="noise" the streams are x1|prior_drift and x1|rng; the
%       drift group's level 0.75 x its declared sigma 0.05 = 0.037 against the
%       planted drift rms 0.031. With r = 6 rows the relative standard error of
%       a variance component is 0.5*sqrt(2/r) = 0.29, so that gap is noise.
%
%   ERRORS  none of its own; declaration and prior errors keep their
%     lsqvce:fit: identifiers -- badPrior for an L of the wrong width, for
%     estimate = "law" on a LOCAL group, and for a local row that couples two
%     days. A COMMON group may couple days: the smoothness group here does,
%     and must, since that is what a smoothness statement is.
%
%   See also LSQVCE.LAYOUT, LSQVCE.FIT, LSQVCE.OPTIONS, DEMO_NOISEPOLICY,
%   DEMO_WORKFLOW
%
%   Created 2026-09-07 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5 -- R6 WP-96: lsqvce.layout had no worked
%   example (HANDOVER 12.19, carried; 12.20).

    p = inputParser;
    p.addParameter('verbose', false);
    p.parse(varargin{:});
    verbose = p.Results.verbose;

    cfg  = priorCfg();
    data = priorData(cfg);

    fprintf(['\n=== Regularisation in the announced layout ===\n' ...
             '    %d harmonic signal coefficients + one offset per day;\n' ...
             '    %d days x %d arcs, each arc its own bias and drift.\n' ...
             '    Offset and arc biases are separable only up to a constant.\n\n'], ...
            cfg.nSig, numel(cfg.days), cfg.arcs);

    % ---- 1. THE LAYOUT ---------------------------------------------------
    % Ask BEFORE the run. Nothing here is your bookkeeping: nc and nl are
    % DERIVED from the expansion of the two scoped tables.
    lay = lsqvce.layout(data);
    s.nc = lay.nc;
    s.nl = lay.systems(1).nl;
    fprintf('lsqvce.layout says: nc = %d, nl = %d\n', s.nc, s.nl);
    fprintf('  A: ');
    for j = 1:lay.nc
        ec = lay.common(j);
        if isempty(ec.day)
            fprintf('%d=%s ', j, data.commonParams{ec.row, 2});
        else
            fprintf('%d=%s(d%d) ', j, data.commonParams{ec.row, 2}, ec.day);
        end
    end
    fprintf('\n  B: ');
    loc = lay.systems(1).local;
    for j = 1:numel(loc)
        el = loc(j);
        nm = data.systems(1).localParams{el.row, 2};
        if isempty(el.day)
            fprintf('%d=%s ', j, nm);
        else
            fprintf('%d=%s(d%d,a%d) ', j, nm, el.day, el.arc);
        end
    end
    fprintf('\n');

    % The columns the datum needs. They cannot be named by a pattern: they are
    % ONE template row instantiated once per day, and the operator couples them.
    dayCols = arrayfun(@(d) lay.colCommon(cfg.rowOffset, d), cfg.days);
    fprintf('the day-offset columns of A are %s\n', mat2str(dayCols));

    % ---- 2. WITHOUT A PRIOR ----------------------------------------------
    % MEASURED, not assumed: what a singular normal matrix actually does.
    s.nullDim = nullDimension(cfg, lay, dayCols);
    fprintf('\nnull directions of the unregularised design: %d (one per day)\n', s.nullDim);
    try
        lsqvce.fit(data, lsqvce.Options(Verbose=false));
        s.unregularised = 'accepted (no refusal)';
    catch ME
        s.unregularised = ME.identifier;
    end
    fprintf('fit without priors: %s\n', s.unregularised);

    % ---- 3. THE DATUM, WRITTEN IN THE LAYOUT -----------------------------
    % Two groups, because they are two statements with two strengths:
    %   'datum'  -- the offsets average to zero.        TIGHT: it is a datum.
    %   'smooth' -- consecutive offsets differ little.  LOOSE: it is a belief.
    nD = numel(cfg.days);
    Lmean = zeros(1, lay.nc);
    Lmean(dayCols) = 1/sqrt(nD);
    Ldiff = zeros(nD-1, lay.nc);
    for k = 1:nD-1
        Ldiff(k, dayCols(k))   = -1;
        Ldiff(k, dayCols(k+1)) =  1;
    end
    data.priors = [struct('L', Lmean, 'x0', 0, 'sigma', cfg.datumSigma, ...
                          'name', 'datum', 'estimate', "none"), ...
                   struct('L', Ldiff, 'x0', zeros(nD-1, 1), 'sigma', cfg.smoothSigma, ...
                          'name', 'smooth', 'estimate', "none")];

    % ---- 4. A LOCAL GROUP, AND THE SAME GROUP BY PATTERN -----------------
    % Damp the arc drifts. Written twice: through the layout (works for ANY
    % operator) and as a `group` pattern on the NAME column (works when the
    % group IS a whole template row, as here). Both must give one solution.
    unit = @(c) double((1:s.nl) == c);
    Ld = zeros(numel(cfg.days) * cfg.arcs, s.nl);
    r  = 0;
    for d = cfg.days
        for a = 1:cfg.arcs
            r = r + 1;
            Ld(r, :) = unit(lay.systems(1).col(cfg.rowDrift, d, a));
        end
    end
    data.systems(1).priors = struct('L', Ld, 'x0', zeros(size(Ld, 1), 1), ...
                                    'sigma', cfg.driftSigma, 'name', 'drift', ...
                                    'estimate', "scale");
    R  = lsqvce.fit(data, lsqvce.Options(Verbose=verbose));
    R0 = fitWithoutGroups(data);

    dataPat = data;
    dataPat.systems(1).priors = struct('group', '^drift$', 'x0', 0, ...
                                       'sigma', cfg.driftSigma, 'name', 'drift', ...
                                       'estimate', "scale");
    Rp = lsqvce.fit(dataPat, lsqvce.Options(Verbose=false));

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
    fprintf('L form vs `group` form: max|dx| = %.1e (round-off)\n', s.formGap);

    % ---- 5. DATUM INVARIANCE, MEASURED -----------------------------------
    % The offsets are NOT observable, so the datum DECIDES them; the signal IS
    % observable, so the datum must NOT touch it. That is a claim about the
    % algebra, and a claim gets measured: move the datum itself -- state that
    % the offsets average to 0.5 instead of 0 -- and look at both halves.
    dataB = data;
    dataB.priors(1).x0 = 0.5;
    RB = lsqvce.fit(dataB, lsqvce.Options(Verbose=false));
    s.datumGap = max(abs(R.x(1:cfg.nSig) - RB.x(1:cfg.nSig)));
    offsetMove = max(abs(R.x(dayCols) - RB.x(dayCols)));
    fprintf('\na different datum: the signal moves %.1e, the day offsets move %.3f\n', ...
            s.datumGap, offsetMove);

    % ---- 6. A PRIOR GROUP IS DATA ----------------------------------------
    % With estimate = "scale" the group is a stream of its own and its sigma is
    % the level the joint estimator computes for every stream (Koch-Kusche).
    % Levels only here (the instrument stream asks for "scale" too): six arcs
    % of 150 epochs are too thin to identify a spectral SHAPE, and a demo must
    % not ask for what its data cannot give -- see DEMO_NOISEPOLICY for the
    % policy field itself.
    dataN = data;
    dataN.systems(1).covariance.rng.estimate = "scale";
    Rn = lsqvce.fit(dataN, lsqvce.Options(Estimate="noise", Pool="stream", ...
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
    fprintf(['\nTHE RULE (the Manual, chapter 3 section 1b): the parameters and their\n' ...
             'columns are yours; the regularisation MECHANISM is the toolbox''s; the\n' ...
             'interface between them is an operator in your columns, announced by\n' ...
             'lsqvce.layout and checked by width.\n']);
end

%==========================================================================
function cfg = priorCfg()
%PRIORCFG  The case: sizes, the planted truth, and the two template rows the
%   demo addresses by index (nothing here is read by the toolbox).
    cfg.nSig        = 4;                % harmonic signal coefficients ('common')
    cfg.days        = 1:3;
    cfg.arcs        = 2;
    cfg.mArc        = 150;              % epochs per arc
    cfg.dt          = 30;               % s
    cfg.sigma       = 0.02;             % the instrument's true noise (known here)
    cfg.rowOffset   = 5;                % commonParams row holding the day offset
    cfg.rowBias     = 1;                % localParams row holding the arc bias
    cfg.rowDrift    = 2;                % localParams row holding the arc drift
    cfg.rowTemp     = 3;                % localParams row holding the 'system' scale
    cfg.datumSigma  = 1e-3;             % the datum row: tight, but see below
    cfg.smoothSigma = 0.20;
    cfg.driftSigma  = 0.05;
    cfg.xSig        = [0.80 -0.45 0.30 0.15];               % the truth (here only)
    cfg.offset      = [0.30 -0.20 0.50];                    % per day; NOT observable
    cfg.bias        = [0.10 -0.05; 0.02 0.08; -0.06 0.03];  % per (day, arc)
    cfg.drift       = [0.04 -0.03; 0.01 0.05; -0.02 0.02];  % per (day, arc), per hour
    cfg.tempScale   = 0.25;             % the 'system' row: one for the whole run
end

function data = priorData(cfg)
%PRIORDATA  The declaration. Both tables are SCOPED: column 1 is the scope
%   token the toolbox reads, column 2 the name the `group` form matches.
    data = struct();
    data.commonParams = [repmat({'common'}, cfg.nSig, 1), ...
                         arrayfun(@(j) sprintf('sig%d', j), (1:cfg.nSig).', ...
                                  'UniformOutput', false)];
    data.commonParams(cfg.rowOffset, :) = {'day', 'offset'};
    data.days = cfg.days;

    s = struct();
    s.buildDay    = @(day, ot, common, local, po, info) priorDay(cfg, day, info);
    s.covariance  = struct('rng', struct('type', 'diag', 'sigma', cfg.sigma));
    s.arcs        = cfg.arcs;
    s.id          = 1;
    s.name        = 'calib';
    s.localParams = {'arc', 'bias'; 'arc', 'drift'; 'system', 'tempScale'};
    data.systems  = s;
end

function D = priorDay(cfg, day, info)
%PRIORDAY  One day of the instrument. DETERMINISTIC: fit reads every day twice
%   per pass, so the noise comes from a stream seeded on the day.
    m = cfg.arcs * cfg.mArc;
    t = (day - 1)*86400 + (0:m-1).' * cfg.dt;
    w = 2*pi / (0.5 * 86400);
    Ssig = [sin(w*t), cos(w*t), sin(2*w*t), cos(2*w*t)];
    temp = cos(2*pi * t / (1.7 * 86400));           % a known auxiliary series

    D.A = zeros(m, info.nc);
    D.A(:, 1:cfg.nSig) = Ssig;
    D.A(:, info.commonDayOffset + 1) = 1;           % this day's offset
    D.B = zeros(m, info.nl);
    D.B(:, info.sysOffset + 1) = temp;              % the 'system' row
    D.l = Ssig * cfg.xSig(:) + cfg.offset(day) + cfg.tempScale * temp;
    D.arcIdx = cell(1, cfg.arcs);
    for a = 1:cfg.arcs
        rows = (a-1)*cfg.mArc + (1:cfg.mArc);
        tau  = (t(rows) - t(rows(1))) / 3600;       % hours since the arc start
        col  = info.offsetOf(a);
        D.B(rows, col + cfg.rowBias)  = 1;
        D.B(rows, col + cfg.rowDrift) = tau;
        D.l(rows) = D.l(rows) + cfg.bias(day, a) + cfg.drift(day, a) * tau;
        D.arcIdx{a} = rows;
    end
    rs  = RandStream('mt19937ar', 'Seed', 4200 + day);
    D.l = D.l + cfg.sigma * randn(rs, m, 1);
    D.t = t;
end

function n = nullDimension(cfg, lay, dayCols)
%NULLDIMENSION  The null space of the FULL (common + local) design, assembled
%   densely HERE only to show what the priors repair. The toolbox never forms
%   it: it eliminates the locals arc by arc.
    M = zeros(0, lay.nc + lay.systems(1).nl);
    for d = cfg.days
        D = priorDay(cfg, d, layoutInfo(cfg, lay, d, dayCols));
        M = [M; D.A, D.B]; %#ok<AGROW>
    end
    n = size(M, 2) - rank(M);
end

function info = layoutInfo(cfg, lay, day, dayCols)
%LAYOUTINFO  buildDay's 6th argument, rebuilt from the layout so the dense
%   assembly above sees the columns fit will. In your own code you NEVER do
%   this -- fit hands you info. Here it keeps the rank check honest.
    info.nc = lay.nc;
    info.nl = lay.systems(1).nl;
    info.commonDayOffset = dayCols(cfg.days == day) - 1;
    info.sysOffset = lay.systems(1).col(cfg.rowTemp, [], []) - 1;
    info.offsetOf  = @(a) lay.systems(1).col(cfg.rowBias, day, a) - 1;
end

function R = fitWithoutGroups(data)
%FITWITHOUTGROUPS  The same solve with the DATUM groups only, so the pseudo-
%   observations of the drift group can be reported as a delta.
    d = data;
    d.systems(1).priors = [];
    R = lsqvce.fit(d, lsqvce.Options(Verbose=false));
end
