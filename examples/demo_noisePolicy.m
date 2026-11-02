function [R, s] = demo_noisePolicy(varargin)
%DEMO_NOISEPOLICY  Per stream: estimate a shape, only a level, or nothing.
%
%   [R, s] = DEMO_NOISEPOLICY()
%   [R, s] = DEMO_NOISEPOLICY('maxPasses', 20, 'verbose', true)
%
%   THE PROBLEM THIS EXAMPLE EXISTS FOR.  Estimate="noise" does not have to
%   treat every stream alike, and it should not: what you know about each
%   instrument differs. A stream therefore carries a POLICY beside its noise
%   law, both given where the stream is declared
%   (LSQVCE.SYSTEMDEFINITION/ADDSTREAM; the API reference, chapter Covariance
%   specifications):
%
%       sys.addStream(name, Build=@build, Noise=lsqvce.Noise.<law>(...), ...
%                     Estimate=policy, ArcFactors=tf)
%       policy   "none" | "scale" | "law"    (default "law")
%       tf       true | false                (default: Options ArcFactors)
%
%     "none"   the declared covariance is used as it stands. For an instrument
%              whose noise you have calibrated and are not willing to have
%              re-negotiated by an adjustment that also has model errors in it.
%     "scale"  the SHAPE is yours, the LEVEL is estimated. The commonest case:
%              a filter response known from the hardware, an amplitude not.
%     "law"    both are estimated. For a stream whose colour is genuinely
%              unknown -- and only where there are enough arcs to identify one.
%
%   ARC FACTORS are the other half: with ArcFactors=true a stream gets one
%   variance factor PER ARC on top of its level, which is how one bad pass, a
%   manoeuvre or an eclipse season is caught instead of being smeared over the
%   whole stream. It costs one component per arc, so it is declared per
%   stream, not globally.
%
%   THE CASE.  One system, three days, three arcs a day, four streams, each in
%   a different state of knowledge and each with a PLANTED truth to recover:
%     gps   white, sigma calibrated and CORRECT       -> "none"
%     acc   shape known (a first-order filter, given  -> "scale"
%           to the toolbox as lsqvce.Noise.psd), level 2x wrong
%     rr    colour unknown (AR(1)), declared white    -> "law"
%     mag   white, but ARC 2 of every day is 3x worse -> "scale" + arcFactors
%
%   ALL FOUR IN ONE RUN, which is worth a sentence of history. Until WP-98
%   this demo had to run TWICE -- the three policies in one run, the arc
%   factors in another -- because a "law" stream and an arcFactors stream in
%   the same run drove the factors to 0.00 / 41 / 0.00 against that planted 9
%   (finding F-01). With eq. 6.5.28 multiplying the factor the weighting
%   actually used, the same four streams settle together: 0.97 / 8.90 / 0.92.
%   One run is now the honest demonstration; two would be a memorial.
%
%   INPUT  (Name-Value, all optional)
%     maxPasses 1x1 double   pass BUDGET of the loop; it settles in 3
%                            (default 10)
%     verbose   1x1 logical  print the engine log (default false)
%   OUTPUT
%     R         1x1 lsqvce.Results   the run, with R.history / R.noise attached
%     s         1x1 struct           the measured summary:
%       .stream       1x4 string   the streams, in R.noise order
%       .declared     1x4 double   the sigma each stream declared
%       .level        1x4 double   the level each stream ended at
%       .truth        1x4 double   the planted level, in the declared units
%       .noneFixed    1x1 logical  the "none" stream's level == its declared
%       .rho1         1x2 double   [estimated planted] lag-1 correlation of rr
%       .arcFactor    1x3 double   mean factor of arcs 1..3 of mag
%       .arcRatio     1x2 double   [estimated planted] bad-arc variance ratio
%       .passes       1x1 double   passes the loop took
%       .badEstimate  1x1 string   the identifier a misspelt policy raises
%
%   ACCURACY  (MEASURED 2026-10-16 on R2026a, v6.25.2, all four streams in ONE
%   run. The fixture noise is drawn from a stream seeded on (day, stream) and
%   two consecutive runs reproduced every figure below, so this demo repeats.)
%     Levels: gps 0.2000 against its declared 0.2000 EXACTLY -- "none" does not
%       move; acc 0.4913 against the planted 0.5000 ("scale" on a given psd
%       shape); rr 1.195 against the planted 1.1547 = 1/sqrt(1-phi^2), 3.5 %
%       high, which is what one AR(1) shape fitted from nine arcs of 256 epochs
%       is worth; mag 0.5728 against 0.5745 = 0.30*sqrt(mean([1 9 1])).
%     The shape rr found: lag-1 correlation 0.550 against the planted 0.500
%       (R.noise.Cgroups, 60 lags).
%     Arc factors of mag: 0.978, 9.03, 0.932 -- the bad arc is 9.45x the good
%       ones against a planted 9 (a VARIANCE ratio: arc 2 is 3x in sigma). The
%       level absorbs their mean, so only that RATIO is identifiable.
%     sigma0 = 0.9956; the loop CONVERGES in 3 passes (R.stopReason
%       "converged") with Noise.Smooth = 1. WP-145, measured: the plain 60-bin
%       law needed 14 passes, more than the 12 this demo spent, and settled on
%       the same figures to three digits (rr 1.199) -- the smoothing changes
%       when the loop may stop, not what it finds.
%
%   LIMITATIONS  (v5.3, MEASURED -- HANDOVER 12.20, 12.21, 12.22)
%     * Of a stream's per-arc factors only the RATIO is identifiable: the
%       level absorbs their mean, which is why the demo reports both and
%       compares only the ratio against its planted 9. (Finding F-01, fixed in
%       WP-98: eq. 6.5.28 multiplied by a factor the weighting never used, so
%       the ratio grew by a constant 8.4x every pass -- 9.4, 79, 671, 5626 --
%       instead of settling. It now settles: 9.4, 8.6, 8.6, 8.6.)
%     * "law" needs the arcs to identify a shape. On this fixture it recovers
%       AR(1) colour up to phi = 0.6; at phi = 0.7 the level runs away to
%       2.8e4. phi = 0.5 is used here for that reason.
%
%   ERRORS  none of its own. A policy that is not one of the three names is
%     refused where it is declared, by addStream, with
%     MATLAB:validators:mustBeMember (shown by the demo, s.badEstimate);
%     asking for "law" on a stream whose sigma is a per-row vector
%     (lsqvce.Noise.white of a column: a FIXED structure, which can carry a
%     level but not a shape) is refused by fit with
%     lsqvce:fit:lawOnFixedStructure.
%
%   See also LSQVCE.TASKDEFINITION, LSQVCE.SYSTEMDEFINITION, LSQVCE.NOISE,
%   LSQVCE.FIT, LSQVCE.OPTIONS, LSQVCE.RESULTS, DEMO_PRIORS, DEMO_GRACE
%
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: declared as a lsqvce.TaskDefinition, policies on addStream. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Created 2026-09-07 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5 -- R6 WP-96: the per-stream policy fields
%   and the 'psd' law had no worked example (HANDOVER 12.19, carried; 12.20).
%
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-124: references to the data
%   contract point at the API reference. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-16 (CEST, UTC+02:00), WP-145: the run converges
%   (Noise.Smooth 1, budget 10, settles in 3), ACCURACY re-measured. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.

    p = inputParser;
    p.addParameter('maxPasses', 10);
    p.addParameter('verbose', false);
    p.parse(varargin{:});
    o = p.Results;

    cfg = policyCfg();
    fprintf(['\n=== The per-stream noise policy ===\n' ...
             '    %d days x %d arcs; four streams, four states of knowledge.\n\n'], ...
            numel(cfg.days), cfg.arcs);
    fprintf('  %-5s %-8s %-11s %s\n', 'strm', 'estimate', 'arcFactors', 'what is known');
    for k = 1:numel(cfg.streams)
        fprintf('  %-5s %-8s %-11s %s\n', cfg.streams{k}, cfg.estimate{k}, ...
                string(cfg.arcFactors(k)), cfg.knowledge{k});
    end

    % ---- THE RUN: all four policies at once ------------------------------
    % ArcFactors=false in Options is the POLICY DEFAULT for the whole run; the
    % ONE stream that needs per-arc factors overrides it on its own addStream.
    % Noise.Smooth 1 (WP-145, measured): the plain 60-bin law of rr stopped
    % on the old budget of 12 passes with bins still creeping; it needed 14.
    % Averaging each bin with its neighbours leaves every figure below the
    % same to three digits and settles in 3.
    sel = 1:numel(cfg.streams);
    R = lsqvce.fit(policyTask(cfg, sel), ...
                   lsqvce.Options(Estimate="noise", Pool="stream", ...
                                  ArcFactors=false, MaxPasses=o.maxPasses, ...
                                  Noise=struct('Smooth', 1), Verbose=o.verbose));
    [s.stream, s.declared, s.level, s.truth] = levels(cfg, R);
    s.passes = R.history.nIter;
    % Look every stream up BY NAME. R.noise.stream follows the order of the
    % addStream calls, but a demo that indexed by position would report the wrong
    % row the day that order changed, and say nothing about it.
    at = @(nm) find(extractAfter(s.stream, '|') == string(nm), 1);
    fprintf('\n--- level per stream (the estimand of "scale") ---\n');
    fprintf('  %-10s %10s %10s %10s   %s\n', 'stream', 'declared', 'level', ...
            'planted', 'policy');
    for k = 1:numel(s.stream)
        j = find(strcmp(cfg.streams, extractAfter(s.stream(k), '|')), 1);
        fprintf('  %-10s %10.3f %10.3f %10.3f   %s\n', s.stream(k), s.declared(k), ...
                s.level(k), s.truth(k), cfg.estimate{j});
    end
    kNone = at(cfg.streams{1});
    s.noneFixed = s.level(kNone) == s.declared(kNone);
    fprintf('  the "none" stream kept its declared sigma EXACTLY: %d\n', s.noneFixed);

    % The shape the "law" stream found. One law per POOL GROUP, and with
    % Pool="stream" the group IS the stream: R.noise.Cgroups is keyed by it and
    % holds the fitted lag covariance, so C(2)/C(1) is its lag-1 correlation.
    kLaw = at(cfg.streams{3});
    Cg = R.noise.Cgroups(char(s.stream(kLaw)));
    s.rho1 = [Cg(2)/Cg(1), cfg.phiRr];
    fprintf('\n--- the shape "law" found (%d lags) ---\n', numel(Cg));
    fprintf('  %s: lag-1 correlation %.3f, planted %.3f\n', s.stream(kLaw), ...
            s.rho1(1), s.rho1(2));
    fprintf('  gps and acc kept the shape they were given -- a psd law for acc.\n');

    % ---- THE ARC FACTORS -------------------------------------------------
    kArc = at(cfg.streams{4});
    f2 = R.noise.model{kArc}.sigma2;                 % one factor per arc
    s.arcFactor = mean(reshape(f2, cfg.arcs, []), 2).';
    good = setdiff(1:cfg.arcs, cfg.badArc);
    s.arcRatio = [s.arcFactor(cfg.badArc) / mean(s.arcFactor(good)), cfg.badScale^2];
    fprintf('\n--- arc factors (only mag declared arcFactors = true) ---\n');
    fprintf('  %s, mean VARIANCE factor per arc: %s\n', s.stream(kArc), ...
            num2str(s.arcFactor, '%.2f  '));
    fprintf('  arc %d against the good arcs: %.1f, planted %d (%dx in sigma)\n', ...
            cfg.badArc, s.arcRatio(1), s.arcRatio(2), cfg.badScale);
    fprintf('  its level %.3f against %.3f: the level absorbs the MEAN factor,\n', ...
            s.level(kArc), s.truth(kArc));
    fprintf('  so only the RATIO between arcs is identifiable.\n');
    fprintf('  the other three streams got ONE level each, no per-arc split.\n');

    % ---- THE POLICY IS CHECKED BY NAME -----------------------------------
    % ... where it is declared: addStream refuses the misspelt name, before any
    % data is read (so the fit below is never reached).
    bad = cfg;
    bad.estimate{1} = "shape";                                       % not a policy
    try
        lsqvce.fit(policyTask(bad, sel), lsqvce.Options(Estimate="noise", Verbose=false));
        s.badEstimate = "accepted (no refusal)";
    catch ME
        s.badEstimate = string(ME.identifier);
    end
    fprintf('\na misspelt policy is refused by name: %s\n', s.badEstimate);

    fprintf('\n--- R.summary(): one row per stream ---\n');
    disp(R.summary());
    fprintf(['HOW TO CHOOSE: "none" when you trust a calibration more than this\n' ...
             'adjustment; "scale" when the shape comes from the hardware and the\n' ...
             'amplitude from the data; "law" only where the arcs can identify one\n' ...
             '(see LIMITATIONS in the help). arcFactors when a single pass can be\n' ...
             'bad on its own -- one component per arc is the price.\n']);
end

%==========================================================================
function [stream, declared, level, truth] = levels(cfg, R)
%LEVELS  What each stream of a run declared, ended at, and should have been.
    stream = string(R.noise.stream);
    n = numel(stream);
    declared = zeros(1, n);  level = zeros(1, n);  truth = zeros(1, n);
    for k = 1:n
        j = find(strcmp(cfg.streams, extractAfter(stream(k), '|')), 1);
        declared(k) = cfg.declared(j);
        truth(k)    = cfg.truth(j);
        level(k)    = R.noise.model{k}.sigma;
    end
end

function cfg = policyCfg()
%POLICYCFG  Sizes, the four streams with their policies, and the planted
%   truth. Nothing here is read by the toolbox; it is this demo's oracle.
    cfg.nc   = 4;                    % common harmonic coefficients
    cfg.days = 1:3;
    cfg.arcs = 3;
    cfg.mArc = 256;                  % epochs per arc
    cfg.dt   = 5;                    % s
    cfg.streams    = {'gps', 'acc', 'rr', 'mag'};
    cfg.estimate   = {"none", "scale", "law", "scale"};
    cfg.arcFactors = [false, false, false, true];
    cfg.knowledge  = {'calibrated sigma, trusted', ...
                      'shape known (psd), level 2x wrong', ...
                      'colour unknown, declared white', ...
                      'white, but one arc a day is 3x worse'};
    cfg.phiAcc   = 0.85;             % the KNOWN first-order filter of acc
    cfg.phiRr    = 0.50;             % the UNKNOWN colour of rr (see LIMITATIONS)
    cfg.badArc   = 2;                % the arc of mag that is degraded
    cfg.badScale = 3;                % ... by this factor in SIGMA
    cfg.xTrue    = [0.80 -0.45 0.30 0.15];
    cfg.bias     = 0.05;             % the per-arc local, same for every arc
    % What each stream DECLARES, what GENERATES the planted error (an innovation
    % std for the coloured streams), and the level that truth amounts to IN THE
    % DECLARED LAW's units -- which is what the estimator can possibly return:
    %   acc  the same shape, so the level is the innovation ratio;
    %   rr   declared white, so the level is the process std, 1/sqrt(1-phi^2);
    %   mag  the level absorbs the MEAN arc factor, hence sqrt(mean([1 9 1])).
    cfg.declared = [0.20, 1.00, 1.00, 0.30];
    cfg.gen      = [0.20, 0.50, 1.00, 0.30];
    arcF = ones(1, cfg.arcs);
    arcF(cfg.badArc) = cfg.badScale^2;
    cfg.truth = [cfg.gen(1), cfg.gen(2), ...
                 cfg.gen(3) / sqrt(1 - cfg.phiRr^2), ...
                 cfg.gen(4) * sqrt(mean(arcF))];
end

function task = policyTask(cfg, sel)
%POLICYTASK  The declaration: one system, the selected streams, one policy
%   each -- Estimate and ArcFactors are given to addStream beside the Noise
%   law. Every stream has a Build of its own (buildGps, buildAcc, buildRr,
%   buildMag); the streams are named by cfg.streams, in the order added.
    task = lsqvce.TaskDefinition(Days=cfg.days);
    task.addCommon(compose("sig%d", 1:cfg.nc));
    sys = task.addSystem("policy", Arcs=cfg.arcs, Id=1);
    sys.addLocal("bias", Scope="arc");
    builds = struct('gps', @buildGps, 'acc', @buildAcc, 'rr', @buildRr, 'mag', @buildMag);
    for k = sel
        if strcmp(cfg.streams{k}, 'acc')
            % THE 'psd' LAW: the shape comes from the hardware -- here the
            % response of the very filter the noise really has -- as a
            % one-sided power spectral density on a frequency grid. The
            % spectrum IS the declared covariance: its level is 1, which is
            % this stream's declared sigma (cfg.declared), and that level is
            % what "scale" estimates.
            [f, S] = firstOrderPsd(cfg.phiAcc, cfg.dt, 2048);
            noise = lsqvce.Noise.psd(f, S);
        else
            noise = lsqvce.Noise.white(cfg.declared(k));
        end
        build = builds.(cfg.streams{k});
        sys.addStream(cfg.streams{k}, Build=@(day, ctx) build(cfg, day, ctx), Noise=noise, ...
                      Estimate=cfg.estimate{k}, ArcFactors=cfg.arcFactors(k));
    end
end

function [f, S] = firstOrderPsd(phi, dt, nf)
%FIRSTORDERPSD  One-sided PSD of the AR(1) filter 1/(1 - phi z^-1) driven by
%   UNIT-variance white noise: S = 2 dt |H(f)|^2, f in (0, Nyquist] Hz.
    f = linspace(1e-9, 0.5/dt, nf).';
    H = 1 ./ (1 - phi * exp(-1i * 2*pi * f * dt));
    S = 2 * dt * abs(H).^2;
end

function out = buildGps(cfg, day, ctx)
%BUILDGPS  The Build of gps: white noise at its calibrated sigma.
    [rs, k] = noiseSource(cfg, day, ctx);
    e = cfg.gen(k) * randn(rs, cfg.arcs * cfg.mArc, 1);
    out = policyDay(cfg, day, ctx, e);
end

function out = buildAcc(cfg, day, ctx)
%BUILDACC  The Build of acc: the KNOWN first-order filter (phi cfg.phiAcc).
    [rs, k] = noiseSource(cfg, day, ctx);
    out = policyDay(cfg, day, ctx, arNoise(cfg, rs, cfg.phiAcc, cfg.gen(k)));
end

function out = buildRr(cfg, day, ctx)
%BUILDRR  The Build of rr: the UNKNOWN AR(1) colour (phi cfg.phiRr).
    [rs, k] = noiseSource(cfg, day, ctx);
    out = policyDay(cfg, day, ctx, arNoise(cfg, rs, cfg.phiRr, cfg.gen(k)));
end

function out = buildMag(cfg, day, ctx)
%BUILDMAG  The Build of mag: white noise, one degraded arc a day.
    [rs, k] = noiseSource(cfg, day, ctx);
    e = cfg.gen(k) * randn(rs, cfg.arcs * cfg.mArc, 1);
    rows = (cfg.badArc-1)*cfg.mArc + (1:cfg.mArc);
    e(rows) = e(rows) * cfg.badScale;
    out = policyDay(cfg, day, ctx, e);
end

function out = policyDay(cfg, day, ctx, e)
%POLICYDAY  One (day, stream) with the planted error E (m x 1, the demo's
%   oracle, drawn by the stream's Build): the signal, the arc bias, the arcs.
%   The local partial is placed by name, ctx.cols("bias", Arc=id), in a dense
%   B of ctx.NumLocal columns.
    m = cfg.arcs * cfg.mArc;
    t = (day - 1)*86400 + (0:m-1).' * cfg.dt;
    w = 2*pi / (0.5 * 86400);
    A = [sin(w*t), cos(w*t), sin(2*w*t), cos(2*w*t)];   % sig1..sig4, the addCommon order
    l = A * cfg.xTrue(:) + e;
    out = ctx.newDay();
    for a = 1:cfg.arcs
        rows = (a-1)*cfg.mArc + (1:cfg.mArc);
        B = zeros(cfg.mArc, ctx.NumLocal);
        B(:, ctx.cols("bias", Arc=ctx.ArcIds(a))) = 1;
        out.addArc(A(rows, :), l(rows) + cfg.bias, t(rows), B=B);
    end
end

function [rs, k] = noiseSource(cfg, day, ctx)
%NOISESOURCE  The random stream of one (day, stream) -- the demo's oracle.
%   DETERMINISTIC: fit reads every day twice per pass, so the noise comes from
%   a stream seeded on the day and on the stream's position k in cfg.streams.
    k  = find(strcmp(cfg.streams, char(ctx.Stream)), 1);
    rs = RandStream('mt19937ar', 'Seed', 5100 + 10*day + k);
end

function e = arNoise(cfg, rs, phi, sigma)
%ARNOISE  AR(1) noise of innovation standard deviation SIGMA, one day of all
%   arcs (cfg.arcs * cfg.mArc rows), its start-up (2000 samples) dropped.
    m  = cfg.arcs * cfg.mArc;
    nB = 2000;                                   % burn-in of the coloured cases
    e = filter(1, [1 -phi], sigma * randn(rs, m + nB, 1));
    e = e(nB+1:end);
end
