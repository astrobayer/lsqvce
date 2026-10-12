function [R, s] = demo_noisePolicy(varargin)
%DEMO_NOISEPOLICY  Per stream: estimate a shape, only a level, or nothing.
%
%   [R, s] = DEMO_NOISEPOLICY()
%   [R, s] = DEMO_NOISEPOLICY('maxPasses', 20, 'verbose', true)
%
%   THE PROBLEM THIS EXAMPLE EXISTS FOR.  Estimate="noise" does not have to
%   treat every stream alike, and it should not: what you know about each
%   instrument differs. The covariance spec of a stream therefore carries a
%   POLICY beside its law (the API reference, chapter Covariance specifications):
%
%       covariance.<type>.estimate   = "none" | "scale" | "law"  (default "law")
%       covariance.<type>.arcFactors = true | false          (default: Options)
%
%     "none"   the declared covariance is used as it stands. For an instrument
%              whose noise you have calibrated and are not willing to have
%              re-negotiated by an adjustment that also has model errors in it.
%     "scale"  the SHAPE is yours, the LEVEL is estimated. The commonest case:
%              a filter response known from the hardware, an amplitude not.
%     "law"    both are estimated. For a stream whose colour is genuinely
%              unknown -- and only where there are enough arcs to identify one.
%
%   ARC FACTORS are the other half: with arcFactors = true a stream gets one
%   variance factor PER ARC on top of its level, which is how one bad pass, a
%   manoeuvre or an eclipse season is caught instead of being smeared over the
%   whole stream. It costs one component per arc, so it is declared per
%   stream, not globally.
%
%   THE CASE.  One system, three days, three arcs a day, four streams, each in
%   a different state of knowledge and each with a PLANTED truth to recover:
%     gps   white, sigma calibrated and CORRECT       -> "none"
%     acc   shape known (a first-order filter, given  -> "scale"
%           to the toolbox as a 'psd' law), level 2x wrong
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
%     maxPasses 1x1 double   passes of the loop (default 12)
%     verbose   1x1 logical  print the engine log (default false)
%   OUTPUT
%     R         1x1 lsqvce.Results   the run, with R.history / R.noise attached
%     s         1x1 struct           the measured summary:
%       .stream       1x4 string   the streams, in R.noise order
%       .declared     1x4 double   the sigma each spec declared
%       .level        1x4 double   the level each stream ended at
%       .truth        1x4 double   the planted level, in the declared units
%       .noneFixed    1x1 logical  the "none" stream's level == its declared
%       .rho1         1x2 double   [estimated planted] lag-1 correlation of rr
%       .arcFactor    1x3 double   mean factor of arcs 1..3 of mag
%       .arcRatio     1x2 double   [estimated planted] bad-arc variance ratio
%       .passes       1x1 double   passes the loop took
%       .badEstimate  1x1 string   the identifier a misspelt policy raises
%
%   ACCURACY  (MEASURED 2026-09-08 on R2026a, all four streams in ONE run. The
%   fixture noise is drawn from a stream seeded on (day, obsType) and two
%   consecutive runs reproduced every figure below, so this demo repeats.)
%     Levels: gps 0.2000 against its declared 0.2000 EXACTLY -- "none" does not
%       move; acc 0.4906 against the planted 0.5000 ("scale" on a given psd
%       shape); rr 1.214 against the planted 1.1547 = 1/sqrt(1-phi^2), 5 % high,
%       which is what one AR(1) shape fitted from nine arcs of 256 epochs is
%       worth; mag 0.5691 against 0.5745 = 0.30*sqrt(mean([1 9 1])).
%     The shape rr found: lag-1 correlation 0.559 against the planted 0.500
%       (R.noise.Cgroups, 60 lags).
%     Arc factors of mag: 0.97, 8.90, 0.92 -- the bad arc is 9.4x the good ones
%       against a planted 9 (a VARIANCE ratio: arc 2 is 3x in sigma). The level
%       absorbs their mean, so only that RATIO is identifiable.
%     sigma0 = 0.997; the loop stops on maxPasses (12), the shape still moving
%       in its last digits.
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
%     refused by fit with lsqvce:fit:badEstimate; asking for "law" on a spec
%     whose sigma is a per-row vector (a FIXED structure, which can carry a
%     level but not a shape) is lsqvce:fit:lawOnFixedStructure.
%
%   See also LSQVCE.FIT, LSQVCE.OPTIONS, LSQVCE.RESULTS, DEMO_PRIORS,
%   DEMO_GRACE
%
%   Created 2026-09-07 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5 -- R6 WP-96: the per-stream policy fields
%   and the 'psd' law had no worked example (HANDOVER 12.19, carried; 12.20).
%
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-124: references to the data
%   contract point at the API reference. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.

    p = inputParser;
    p.addParameter('maxPasses', 12);
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
    % ONE stream that needs per-arc factors overrides it on its own spec.
    sel = 1:numel(cfg.streams);
    R = lsqvce.fit(policyData(cfg, sel), ...
                   lsqvce.Options(Estimate="noise", Pool="stream", ...
                                  ArcFactors=false, MaxPasses=o.maxPasses, ...
                                  Verbose=o.verbose));
    [s.stream, s.declared, s.level, s.truth] = levels(cfg, R);
    s.passes = R.history.nIter;
    % Look every stream up BY NAME. R.noise.stream follows the covariance
    % fieldnames, but a demo that indexed by position would report the wrong
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
    bad = policyData(cfg, sel);
    bad.systems(1).covariance.(cfg.streams{1}).estimate = "shape";   % not a policy
    try
        lsqvce.fit(bad, lsqvce.Options(Estimate="noise", Verbose=false));
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
    % What each spec DECLARES, what GENERATES the planted error (an innovation
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

function data = policyData(cfg, sel)
%POLICYDATA  The declaration: one system, the selected streams, one policy
%   each. The fieldnames of covariance ARE the streams (the API reference).
    data = struct();
    data.commonParams = [repmat({'common'}, cfg.nc, 1), ...
                         arrayfun(@(j) sprintf('sig%d', j), (1:cfg.nc).', ...
                                  'UniformOutput', false)];
    data.days = cfg.days;

    cov = struct();
    for k = sel
        if strcmp(cfg.streams{k}, 'acc')
            % THE 'psd' LAW: the shape comes from the hardware -- here the
            % response of the very filter the noise really has -- as a
            % one-sided power spectral density on a frequency grid. The LEVEL
            % rides on cm.sigma, and that is what "scale" estimates.
            [f, S] = firstOrderPsd(cfg.phiAcc, cfg.dt, 2048);
            cm = struct('type', 'psd', 'f', f, 'S', S, 'sigma', cfg.declared(k));
        else
            cm = struct('type', 'diag', 'sigma', cfg.declared(k));
        end
        cm.estimate   = cfg.estimate{k};
        cm.arcFactors = cfg.arcFactors(k);
        cov.(cfg.streams{k}) = cm;
    end

    s = struct();
    s.buildDay    = @(day, ot, common, local, po, info) policyDay(cfg, day, ot, info);
    s.covariance  = cov;
    s.arcs        = cfg.arcs;
    s.id          = 1;
    s.name        = 'policy';
    s.localParams = {'arc', 'bias'};
    data.systems  = s;
end

function [f, S] = firstOrderPsd(phi, dt, nf)
%FIRSTORDERPSD  One-sided PSD of the AR(1) filter 1/(1 - phi z^-1) driven by
%   UNIT-variance white noise: S = 2 dt |H(f)|^2, f in (0, Nyquist] Hz.
    f = linspace(1e-9, 0.5/dt, nf).';
    H = 1 ./ (1 - phi * exp(-1i * 2*pi * f * dt));
    S = 2 * dt * abs(H).^2;
end

function D = policyDay(cfg, day, obsType, info)
%POLICYDAY  One (day, stream). DETERMINISTIC: fit reads every day twice per
%   pass, so the noise comes from a stream seeded on (day, obsType).
    m = cfg.arcs * cfg.mArc;
    t = (day - 1)*86400 + (0:m-1).' * cfg.dt;
    w = 2*pi / (0.5 * 86400);
    A = [sin(w*t), cos(w*t), sin(2*w*t), cos(2*w*t)];
    l = A * cfg.xTrue(:) + policyNoise(cfg, day, obsType);
    % v6: THE CELLS ARE THE ARCS -- one cell per arc for A, B, l, t
    D = struct('A', {cell(1, cfg.arcs)}, 'B', {cell(1, cfg.arcs)}, ...
               'l', {cell(1, cfg.arcs)}, 't', {cell(1, cfg.arcs)});
    for a = 1:cfg.arcs
        rows = (a-1)*cfg.mArc + (1:cfg.mArc);
        B = zeros(cfg.mArc, info.nl);
        B(:, info.offsetOf(a) + 1) = 1;
        D.A{a} = A(rows, :);
        D.B{a} = B;
        D.l{a} = l(rows) + cfg.bias;
        D.t{a} = t(rows);
    end
end

function e = policyNoise(cfg, day, obsType)
%POLICYNOISE  The planted error of one (day, stream) -- the demo's oracle.
    k  = find(strcmp(cfg.streams, obsType), 1);
    rs = RandStream('mt19937ar', 'Seed', 5100 + 10*day + k);
    m  = cfg.arcs * cfg.mArc;
    nB = 2000;                                   % burn-in of the coloured cases
    switch obsType
        case 'acc'
            e = filter(1, [1 -cfg.phiAcc], cfg.gen(k) * randn(rs, m + nB, 1));
            e = e(nB+1:end);
        case 'rr'
            e = filter(1, [1 -cfg.phiRr], cfg.gen(k) * randn(rs, m + nB, 1));
            e = e(nB+1:end);
        otherwise
            e = cfg.gen(k) * randn(rs, m, 1);
    end
    if strcmp(obsType, 'mag')                    % one degraded arc a day
        rows = (cfg.badArc-1)*cfg.mArc + (1:cfg.mArc);
        e(rows) = e(rows) * cfg.badScale;
    end
end
