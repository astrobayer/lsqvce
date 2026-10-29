function [R, s] = demo_tutorial(varargin)
%DEMO_TUTORIAL  Your first streaming adjustment, from declaration to verdict.
%
%   [R, s] = DEMO_TUTORIAL()
%   [R, s] = DEMO_TUTORIAL('verbose', true, 'plot', true)
%
%   THE EXAMPLE OF THE MANUAL'S TUTORIAL (chapter "Your first adjustment").
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
%   They are computed ONCE per day in the system's prepareDay and reach both
%   streams' buildDay as procOpt.prepared.
%
%   WHAT IT DOES
%     1. declares the problem as a data struct (commonParams, days, systems);
%     2. runs ONE solve under the declared noise (Estimate = "none");
%     3. estimates the noise of both streams (Estimate = "noise") and solves again;
%     4. prints the report and the estimates against the planted truth.
%   A simulator stands in for your file reader: the noise of each (day,
%   stream) is drawn from a stream seeded on the pair, so every call returns
%   the same rows (buildDay must be deterministic) and the demo repeats.
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
%   ERRORS  none of its own; every declaration error is fit's (lsqvce:fit:*).
%
%   See also LSQVCE.FIT, LSQVCE.OPTIONS, LSQVCE.RESULTS, DEMO_WORKFLOW,
%   DEMO_FIT
%
%   Created 2026-10-29 (CET, UTC+01:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-159: the Manual had no runnable
%   first streaming case, and prepareDay no worked example outside
%   template/kbr_variational.

    p = inputParser;
    p.addParameter('verbose', false);
    p.addParameter('plot', false);
    p.parse(varargin{:});
    o = p.Results;

    % ---- Step 1: the parameters ----------------------------------------
    data = struct();
    data.commonParams = {'common', 'cosAmp'; 'common', 'sinAmp'};
    data.days = 1:5;

    % ---- Step 2: the system: its locals, arcs, readers and noise ---------
    sys = struct();
    sys.name        = 'platform';
    sys.localParams = {'arc', 'bias'};                 % one offset per arc
    sys.arcs        = 2;                               % arcs 1 and 2 on every day
    sys.processingOptions = struct('epochs', 240, 'step', 30, ...
                                   'truth', [0.8; -0.3]);   % the simulator's truth
    sys.prepareDay  = @prepareTutorialDay;             % once per day
    sys.buildDay    = @buildTutorialDay;               % once per day and stream
    sys.covariance.gauge = struct('type', 'diag', 'sigma', 0.05);
    sys.covariance.radar = struct('type', 'diag', 'sigma', 0.02);   % colour unknown
    data.systems = sys;

    % ---- Step 3: one solve under the declared noise ----------------------
    R1 = lsqvce.fit(data, Verbose=o.verbose);
    T1 = R1.coefficients();

    % ---- Step 4: estimate the noise of both streams and solve again ------
    opts = lsqvce.Options(Estimate="noise", Seed=1, Verbose=o.verbose);
    R = lsqvce.fit(data, opts);
    verdict = R.report();
    T = R.coefficients();
    S = R.summary();

    % ---- Step 5: compare with the truth -----------------------------------
    truth = sys.processingOptions.truth;
    fprintf('\n%-8s %8s %20s %20s\n', 'param', 'truth', 'declared noise', 'estimated noise');
    for k = 1:2
        fprintf('%-8s %8.4f %10.4f +- %6.4f %10.4f +- %6.4f\n', data.commonParams{k, 2}, ...
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

function P = prepareTutorialDay(day, ~, ~, procOpt, info)
%PREPARETUTORIALDAY  What both streams of DAY share, computed once: the
%   epochs of every arc, the signal basis there (the design rows of the two
%   common parameters) and the true signal plus the arc's offset. In a real
%   case this is the expensive part -- an orbit integration, a file every
%   stream reads.
    nArc = numel(info.arcId);
    P = struct('t', {cell(1, nArc)}, 'G', {cell(1, nArc)}, 'y', {cell(1, nArc)});
    for a = 1:nArc
        start  = (day - 1) * 86400 + (info.arcId(a) - 1) * 43200;   % arcs at 0 h and 12 h
        t      = start + (0:procOpt.epochs - 1).' * procOpt.step;
        w      = 2 * pi / 3600;
        G      = [cos(w * t), sin(w * t)];
        rs     = RandStream('mt19937ar', 'Seed', 100 * day + info.arcId(a));
        bias   = 0.5 * randn(rs);                      % this arc's true offset
        P.t{a} = t;
        P.G{a} = G;
        P.y{a} = G * procOpt.truth + bias;
    end
end

function D = buildTutorialDay(day, obsType, ~, ~, procOpt, info)
%BUILDTUTORIALDAY  The rows of one stream on DAY, one cell per arc: the
%   shared part from procOpt.prepared, this instrument's noise on top.
    P    = procOpt.prepared;                           % prepareTutorialDay's result
    nArc = numel(info.arcId);
    D = struct('A', {cell(1, nArc)}, 'B', {cell(1, nArc)}, ...
               'l', {cell(1, nArc)}, 't', {cell(1, nArc)});
    for a = 1:nArc
        m  = numel(P.t{a});
        rs = RandStream('mt19937ar', 'Seed', 1000 * day + 10 * info.arcId(a) + double(obsType(1)));
        switch obsType
            case 'gauge'
                e = 0.05 * randn(rs, m, 1);                            % white
            case 'radar'
                phi = 0.5;                                             % AR(1), sigma 0.02
                e = filter(1, [1 -phi], 0.02 * sqrt(1 - phi^2) * randn(rs, m + 200, 1));
                e = e(201:end);                                        % drop the start-up
        end
        B = zeros(m, info.nl);
        B(:, info.offsetOf(info.arcId(a)) + (1:info.count(a))) = 1;   % this arc's bias
        D.A{a} = P.G{a};                                               % m x 2
        D.B{a} = B;                                                    % m x 10
        D.l{a} = P.y{a} + e;                                           % observed - computed (computed = 0)
        D.t{a} = P.t{a};                                               % seconds
    end
end
