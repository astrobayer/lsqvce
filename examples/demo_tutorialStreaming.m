function [R, s] = demo_tutorialStreaming(varargin)
%DEMO_TUTORIALSTREAMING  Tutorial 3: a streaming adjustment, from declaration to verdict.
%
%   [R, s] = DEMO_TUTORIALSTREAMING()
%   [R, s] = DEMO_TUTORIALSTREAMING('verbose', true, 'plot', true)
%
%   THE EXAMPLE OF THE MANUAL'S TUTORIAL 3 (chapter "Tutorials", part 3).
%   Every line of the chapter is a line of this file, so the chapter runs
%   whenever the test suite runs this file. Read the two side by side.
%
%   THE CASE.  A platform watches a periodic signal of one hour period,
%       y(t) = a cos(2 pi t / 3600) + b sin(2 pi t / 3600),
%   with TWO instruments on one clock, a 'gauge' and a 'radar'. It records two
%   arcs a day of two hours each (240 epochs at 30 s), on five days. Every arc
%   starts with an unknown offset of its own, which both instruments see.
%     common parameters   a, b                 (2, what the run is for)
%     local parameters    one bias per arc     (10, nuisance, eliminated)
%     streams             'gauge'  white noise, sigma 0.05, declared right
%                         'radar'  AR(1) noise (phi 0.5), sigma 0.02, declared
%                                  WHITE -- the estimator has to find the colour
%   The epochs and the signal basis of a day are what both instruments share.
%   They are computed ONCE per day in the system's Prepare and reach both
%   streams' Build as ctx.Prepared.
%
%   WHAT IT DOES
%     1. declares the task: lsqvce.TaskDefinition, its parameters, one system
%        with a Prepare and two streams with a Build each;
%     2. runs ONE solve under the declared noise (Estimate = "none");
%     3. estimates the noise of both streams (Estimate = "noise") and solves again;
%     4. prints the report and the estimates against the planted truth.
%   A simulator stands in for your file reader: the noise of each (day,
%   stream) is drawn from a stream seeded on the pair, so every call returns
%   the same rows (Build must be deterministic) and the demo repeats.
%
%   INPUT  (Name-Value, all optional)
%     verbose  1x1 logical  print the engine's progress lines (default false)
%     plot     1x1 logical  draw R.plot("psd") and R.plot("residuals") (default false)
%   OUTPUT
%     R        1x1 lsqvce.Results   the noise-estimation run
%     s        1x1 struct           the measured summary:
%       .truth      2x1 double   the planted [a; b]
%       .first      1x1 struct   the declared-noise run: .x, .se, .sigma0
%       .final      1x1 struct   the noise run: .x, .se, .sigma0, .stopReason
%       .factor     1x2 double   the final sigma of [gauge radar] RELATIVE to
%                                the declared one (R.summary().sigma)
%       .verdict    1x1 string   the report's verdict ("PASS" here)
%
%   ACCURACY  (MEASURED 2026-10-29 on R2026a, v6.35.1, MATLAB MCP bridge.)
%     declared noise:  a = 0.8011 +- 0.0005, b = -0.3006 +- 0.0005, sigma0 0.993
%     estimated noise: a = 0.8012 +- 0.0008, b = -0.3004 +- 0.0008, sigma0 1.000,
%                      "converged" in 4 passes, factors 0.985 (gauge) and
%                      0.999 (radar), report verdict PASS.
%     The lesson: the declared run's sigma0 is fine, because the radar's
%     LEVEL was declared right -- but its formal errors are 1.6x too small,
%     because its COLOUR was not. Only the noise run's errors cover the truth
%     (0.8, -0.3) honestly.
%
%   ERRORS  none of its own; a declaration error is raised by the add-call
%           that makes it (lsqvce:TaskDefinition:*).
%
%   See also LSQVCE.TASKDEFINITION, LSQVCE.FIT, LSQVCE.OPTIONS, LSQVCE.RESULTS,
%   DEMO_TUTORIALMATRIX, DEMO_TUTORIALBATCH, DEMO_WORKFLOW, DEMO_FIT
%
%   Updated 2026-11-03 (CET, UTC+01:00), WP-163: renamed from demo_tutorial, the
%   third of three tutorials (matrix, batched, streaming). With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: declared as a
%   lsqvce.TaskDefinition -- Prepare on the system, one Build per stream, the
%   arc's bias by name. With assistance from Claude (Anthropic), model Claude
%   Opus 5.5.
%   Created 2026-10-29 (CET, UTC+01:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-159: the Manual had no runnable
%   first streaming case, and prepareDay no worked example outside
%   template/kbr_variational.

    p = inputParser;
    p.addParameter('verbose', false);
    p.addParameter('plot', false);
    p.parse(varargin{:});
    o = p.Results;

    % ---- Step 1: the task and its parameters -------------------------------
    task = lsqvce.TaskDefinition(Days=1:5);
    names = ["cosAmp" "sinAmp"];
    task.addCommon(names);                             % Scope="common" (default)

    % ---- Step 2: the system -- its Prepare, locals and streams --------------
    truth = [0.8; -0.3];                               % the simulator's truth
    sys = task.addSystem("platform", Prepare=@prepareTutorialDay, Arcs=2, ...   % arcs 1 and 2
                         Settings=struct('epochs', 240, 'step', 30, 'truth', truth));
    sys.addLocal("bias", Scope="arc");                 % one offset per arc
    sys.addStream("gauge", Build=@buildGauge, Noise=lsqvce.Noise.white(0.05));
    sys.addStream("radar", Build=@buildRadar, Noise=lsqvce.Noise.white(0.02));   % colour unknown

    % ---- Step 3: one solve under the declared noise ----------------------
    R1 = lsqvce.fit(task, Verbose=o.verbose);
    T1 = R1.coefficients();

    % ---- Step 4: estimate the noise of both streams and solve again ------
    opts = lsqvce.Options(Estimate="noise", Seed=1, Verbose=o.verbose);
    R = lsqvce.fit(task, opts);
    verdict = R.report();
    T = R.coefficients();
    S = R.summary();

    % ---- Step 5: compare with the truth -----------------------------------
    fprintf('\n%-8s %8s %20s %20s\n', 'param', 'truth', 'declared noise', 'estimated noise');
    for k = 1:2
        fprintf('%-8s %8.4f %10.4f +- %6.4f %10.4f +- %6.4f\n', names(k), ...
                truth(k), T1.estimate(k), T1.se(k), T.estimate(k), T.se(k));
    end
    fprintf('sigma0   %8s %20.3f %20.3f   (%s)\n', '', R1.sigma0, R.sigma0, R.stopReason);

    if o.plot
        R.plot("psd");
        R.plot("residuals");
    end

    s = struct('truth', truth, ...
               'first', struct('x', T1.estimate, 'se', T1.se, 'sigma0', R1.sigma0), ...
               'final', struct('x', T.estimate, 'se', T.se, 'sigma0', R.sigma0, ...
                               'stopReason', R.stopReason), ...
               'factor', S.sigma(:).', ...
               'verdict', string(verdict.verdict));
end

function P = prepareTutorialDay(day, ctx)
%PREPARETUTORIALDAY  What both streams of DAY share, computed once: the
%   epochs of every arc, the signal basis there (the design rows of the two
%   common parameters) and the true signal plus the arc's offset. In a real
%   case this is the expensive part -- an orbit integration, a file every
%   stream reads.
    cfg  = ctx.Settings;
    nArc = numel(ctx.ArcIds);
    P = struct('t', {cell(1, nArc)}, 'G', {cell(1, nArc)}, 'y', {cell(1, nArc)});
    for a = 1:nArc
        start  = (day - 1) * 86400 + (ctx.ArcIds(a) - 1) * 43200;  % arcs at 0 h and 12 h
        t      = start + (0:cfg.epochs - 1).' * cfg.step;
        w      = 2 * pi / 3600;
        G      = [cos(w * t), sin(w * t)];
        rs     = RandStream('mt19937ar', 'Seed', 100 * day + ctx.ArcIds(a));
        bias   = 0.5 * randn(rs);                      % this arc's true offset
        P.t{a} = t;
        P.G{a} = G;
        P.y{a} = G * cfg.truth + bias;
    end
end

function out = buildGauge(day, ctx)
%BUILDGAUGE  The gauge's rows on DAY: the shared part plus white noise (0.05).
    out = instrumentDay(day, ctx, @(rs, m) 0.05 * randn(rs, m, 1));
end

function out = buildRadar(day, ctx)
%BUILDRADAR  The radar's rows on DAY: the shared part plus AR(1) noise (phi 0.5,
%   sigma 0.02) -- declared white, so the estimator has to find the colour.
    phi = 0.5;
    out = instrumentDay(day, ctx, @(rs, m) arNoise(rs, m, phi, 0.02));
end

function out = instrumentDay(day, ctx, noise)
%INSTRUMENTDAY  One stream's arcs on DAY: the design from ctx.Prepared, the
%   instrument's noise on top, and the arc's bias by name.
    P   = ctx.Prepared;                                % prepareTutorialDay's result
    out = ctx.newDay();
    c   = char(ctx.Stream);
    for a = 1:numel(ctx.ArcIds)
        m  = numel(P.t{a});
        rs = RandStream('mt19937ar', 'Seed', 1000 * day + 10 * ctx.ArcIds(a) + double(c(1)));
        out.addArc(P.G{a}, P.y{a} + noise(rs, m), P.t{a}, ...    % A, observed - computed, t [s]
                   Local={"bias", ones(m, 1)});                  % this arc's offset
    end
end

function e = arNoise(rs, m, phi, sigma)
%ARNOISE  AR(1) noise of marginal standard deviation SIGMA, start-up dropped.
    e = filter(1, [1 -phi], sigma * sqrt(1 - phi^2) * randn(rs, m + 200, 1));
    e = e(201:end);
end
