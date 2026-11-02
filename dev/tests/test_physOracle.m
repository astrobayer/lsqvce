function tests = test_physOracle
%TEST_PHYSORACLE  WP-48 tier (a): the physically structured closed-loop oracle.
%
%   The Gaussian fixtures prove the ALGEBRA; this suite proves the loop on
%   the DESIGN STRUCTURE production will meet (HANDOVER 10.9/10.10): the
%   design is built from real spherical-harmonic potential partials along a
%   deterministic near-polar circular orbit with Earth rotation
%   (examples/demo_gravity_design) -- order-block correlation, polar-gap
%   conditioning, (1/r)^n attenuation -- NO random design entries anywhere.
%   Truth is DETERMINISTIC: a closed-formula coefficient set with one
%   planted anomaly; the noise LAW is deterministic per stream (AR colour +
%   sigma), only its realisation is seeded-random. Two satellites (two
%   systems, phase-offset orbits) observe the same field with different
%   planted noise levels; the closed loop must recover coefficients within
%   the formal errors, the error degree variances must match the formal
%   ones, and VCE must recover the planted weight ratio under REALISTIC
%   column correlation.
%
%   Run:
%     results = runtests('tests/test_physOracle');  disp(table(results))
%
%   See also DEMO_GRAVITY_DESIGN (examples/), LSQVCE.ENGINE.ADJUSTVCE, HANDOVER 10.10
%
%   Created 2026-08-20 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.
%   Updated 2026-10-17 (CEST, UTC+02:00), WP-148: the template's flagged
%   epochs stay on the grid with weight 0 (no gap, the dof unchanged, their
%   values irrelevant). With assistance from Claude (Anthropic), model Claude
%   Opus 5.5.
%   Updated 2026-10-17 (CEST, UTC+02:00), WP-144: template/gravfire_data.m
%   runs on a generated bundle (A7); a local group's members must agree (A4).
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-16 (CEST, UTC+02:00), WP-145: demo_gravity runs as
%   shipped and must converge; its statistics are banded. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
    tests = functiontests(localfunctions);
end

%% ---------------------------------------------------------------- fixtures
function setupOnce(tc)
    pathBefore = path;   % R7-45: a standalone run (run_lsqvce / runtests) leaves the path as found
    tc.addTeardown(@() path(pathBefore));
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % dev/tests -> root (L8)
    if exist('lsqvce.engine.adjustVCEjoint', 'file') ~= 2 && ...
       exist('+lsqvce', 'dir') ~= 7
        addpath(root);                                       %#ok<MCAP>
    end
    % The SH design generator lives in examples/ so the demo and this suite
    % share ONE implementation (WP-51) -- same route TEST_ADJUSTVCE uses to
    % reach DEMO_GRACE_DATALOADER.
    % WP-122: template/kbr_variational (the variational-equation builder) is
    % tested here, against its own physical oracles
    addpath(fullfile(root, 'template'));                     %#ok<MCAP>
    exDir = fullfile(root, 'examples');
    if exist(fullfile(exDir, 'demo_gravity_design.m'), 'file') == 2
        addpath(exDir);                                      %#ok<MCAP>
    end
    cfg = oracleCfg();
    tc.TestData.cfg = cfg;                       % the geometry oracle (design, Legendre)
    % R5 M4 4c (decision R5-D4): the RECOVERY oracle runs on the USER surface
    % with the one estimator -- the gravity demo's planted field, planted
    % stream sigmas and ragged schedule, through lsqvce.fit(Estimate="noise", Pool="stream", ArcFactors=false).
    % ONE closed-loop solve shared by the recovery tests. (Before 4c:
    % lsqvce.engine.adjustVCE on the oracle loaders -- the separate mode,
    % retired with this level.)
    [cfgR, dataR] = demo_gravity_dataLoader();
    w = warning('off', 'lsqvce:adjustVCEjoint:notConverged'); cleaner = onCleanup(@() warning(w)); %#ok<NASGU>
    tc.TestData.cfgR = cfgR;
    tc.TestData.R = lsqvce.fit(dataR, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=6, Seed=1, Verbose=false));
end

%% ------------------------------------------------- the analytic anchors
function testLegendreAdditionTheorem(tc)
% categories: reference
% The addition-theorem identity Pbar_n0^2 + sum_m Pbar_nm^2 = 2n+1 holds for
% every degree at every latitude -- the ONE mechanical check that certifies
% a fully-normalised Legendre implementation as a whole (HANDOVER 6.1.4,
% authority "Analytic SH algebra"). Probed at nmax = 12, beyond the suite's
% working degree, including near-pole latitudes.
    t = [-0.999; -0.7; -0.3; 0; 0.25; 0.6; 0.95; 0.9999];
    orbit = struct('u', asin(t), 'lam0', zeros(size(t)), ...
                   'incl', pi/2, 'r', 1.07);
    [A, meta] = demo_gravity_design(12, orbit);                 %#ok<ASGLU>
    r = orbit.r;
    for n = 2:12
        cols = meta.n == n;
        % undo the (1/r)^(n+1) attenuation, re-square-sum: at fixed lam = 0
        % the C/S split re-assembles cos^2 + sin^2 exactly.
        S = sum((A(:, cols) * r^(n + 1)).^2, 2);
        % tol: eps-level growth of an O(n) recursion; measured <= 5e-13 at
        % n = 12 (ACCURACY block of demo_gravity_design)
        tc.verifyEqual(S, (2*n + 1) * ones(size(S)), 'AbsTol', 5e-13, ...
            sprintf('addition theorem fails at degree %d', n));
    end
end

function testLegendreClosedForms(tc)
% categories: reference
% Closed forms of the low-degree fully-normalised functions, the point
% anchors complementing the addition theorem's sum anchor.
    t = 0.37;
    orbit = struct('u', asin(t), 'lam0', 0, 'incl', pi/2, 'r', 2);
    [A, meta] = demo_gravity_design(3, orbit);
    r = orbit.r;
    P20 = A(1, meta.n == 2 & meta.m == 0) * r^3;
    P22c = A(1, meta.n == 2 & meta.m == 2 & meta.cs == 'C') * r^3;
    % tol: exact closed forms, 1e-14 ~ 100*eps
    tc.verifyEqual(P20, sqrt(5) * (3*t^2 - 1) / 2, 'AbsTol', 1e-14);
    tc.verifyEqual(P22c, sqrt(15)/2 * (1 - t^2) * cos(0), ...
                   'AbsTol', 1e-14);   % tol: same source
end

function testGravityDesignContract(tc)
% categories: contract
% Shapes, column bookkeeping, DETERMINISM (the design carries no random
% numbers: two calls are bitwise equal), and the documented rejections.
    orbit = struct('u', (0:49).' * 0.1, 'lam0', -(0:49).' * 0.004, ...
                   'incl', deg2rad(89), 'r', 1.07);
    [A1, meta] = demo_gravity_design(6, orbit);
    tc.verifySize(A1, [50, (6+1)^2 - 4]);
    tc.verifyEqual(numel(meta.n), size(A1, 2));
    tc.verifyEqual(sum(meta.cs == 'C') - sum(meta.cs == 'S'), 6 - 1);
    [A2, ~] = demo_gravity_design(6, orbit);
    tc.verifyEqual(A2, A1);                 % bitwise: deterministic design

    caught = MException.empty;
    try
        demo_gravity_design(1, orbit);
    catch e
        caught = e;
    end
    tc.verifyNotEmpty(caught);
    tc.verifyEqual(caught.identifier, 'lsqvce:gravityDesign:badNmax');
    caught = MException.empty;
    try
        demo_gravity_design(6, struct('u', 1, 'lam0', [1 2], 'incl', 1, 'r', 1.07));
    catch e
        caught = e;
    end
    tc.verifyEqual(caught.identifier, 'lsqvce:gravityDesign:badOrbit');
end

%% ------------------------------------ the fixture is NOT another Gaussian
function testDesignIsNotBenign(tc)
% categories: robustness
% The oracle's reason to exist: a randn design has near-orthogonal columns
% (max |corr| ~ 1/sqrt(M)); THIS design must carry the real structure. The
% floors below are deliberately far above anything a Gaussian design can
% reach at this size, so if a refactor ever degenerates the fixture into
% another benign one, this test -- not a silent loss of power -- goes red.
    cfg = tc.TestData.cfg;
    A = designOfDay(cfg, 'A', 1);
    N = A.' * A;
    d = sqrt(diag(N));
    C = abs(N ./ (d * d.'));  C(1:size(C,1)+1:end) = 0;
    % tol floors: MEASURED on this fixture 2026-08-20 (pf48_measure, user
    % machine): max|corr| = 0.956 and cond(N) = 4.05e18 for ONE day -- a
    % single satellite-day is near-singular, exactly like a real single-arc
    % gravity problem; the 4-day multi-orbit stack resolves it. A randn
    % design of the same size measures max|corr| = 0.246, cond = 7.99.
    % Floors at 0.45 / 1e3 sit between the two regimes with wide margin.
    tc.verifyGreaterThan(max(C(:)), 0.45, ...
        'the design has lost its physical column correlation');
    tc.verifyGreaterThan(cond(N), 1e3, ...
        'the design has lost its physical conditioning');
end

%% ------------------------------------------------ closed-loop recoveries
function testCoefficientRecoveryWithinFormalErrors(tc)
% categories: precision
% The closed loop: x was PLANTED (deterministic formula + one anomaly);
% the estimate must sit within its own formal errors. Two bright lines:
% no single normalised error beyond 5 sigma, and the mean squared
% normalised error inside the chi^2 band.
    R  = tc.TestData.R;
    cfg = tc.TestData.cfgR;                 % the demo's planted truth (4c)
    sig = sqrt(diag(R.sol.Cxx)) * R.sigma0;
    z = (R.x - cfg.xTrue) ./ sig;
    tc.verifyLessThan(max(abs(z)), 5, ...
        'a coefficient left its formal error band');   % tol: 5-sigma bright line, P(false alarm) < 77*6e-7
    % tol: mean z^2 ~ chi2_K/K, sd = sqrt(2/K) with K = 77 -> 0.16; band 4 sd
    tc.verifyEqual(mean(z.^2), 1, 'AbsTol', 0.65, ...
        'the formal errors do not describe the actual errors');
end

function testErrorDegreeVariancesMatchFormal(tc)
% categories: precision | tags: diagnostics
% Per degree: the realised error power against its formal expectation --
% the curve a gravity-field analyst actually reads. Ratio bands per degree
% are wide (few coefficients per degree), the SUM is tight.
    R = tc.TestData.R;
    cfg = tc.TestData.cfgR;                 % the demo's planted truth (4c)
    sig = sqrt(diag(R.sol.Cxx)) * R.sigma0;
    err = R.x - cfg.xTrue;
    ratios = zeros(1, cfg.nmax - 1);
    for n = 2:cfg.nmax
        k = cfg.meta.n == n;
        ratios(n - 1) = sum(err(k).^2) / sum(sig(k).^2);
    end
    % tol: each ratio ~ chi2_(2n+1)/(2n+1); the loosest band (n=2, 5 dof)
    % spans [0.08, 3.5] at 3 sigma -- asserted per degree
    tc.verifyGreaterThan(min(ratios), 0.05);
    tc.verifyLessThan(max(ratios), 4.0);
    % tol: pooled over K=77: chi2_K/K, 4*sqrt(2/K) band
    tc.verifyEqual(sum(ratios .* (2*(2:cfg.nmax) + 1)) / 77, 1, ...
                   'AbsTol', 0.65);
end

function testVceRecoversPlantedRatioOnStructuredDesign(tc)
% categories: precision, reference
% THE WP-48 question: does VCE still identify the planted weight ratio when
% the columns are correlated the way a real gravity design correlates them?
% (On the Gaussian fixtures this is proven at 0.20; structure inflates the
% estimator scatter, so the band here is measured on THIS fixture.)
    R = tc.TestData.R;
    cfg = tc.TestData.cfgR;                 % the demo's planted truth (4c)
    got = zeros(1, numel(R.noise.stream));
    exp = zeros(1, numel(R.noise.stream));
    for k = 1:numel(R.noise.stream)
        ot = extractAfter(R.noise.stream{k}, '|');
        got(k) = R.noise.model{k}.sigma;
        exp(k) = cfg.trueSigma.(ot);
    end
    gotR = got / got(1);   expR = exp / exp(1);
    % tol: planted-ratio recovery under structured correlation; 0.25 =
    % Gaussian fixture's 0.20 + headroom for correlation-inflated scatter,
    % verified green on the pinned seed (postflight 2026-08-20, 8/8)
    tc.verifyEqual(gotR, expR, 'RelTol', 0.25, ...
        'VCE lost the planted weighting under physical column correlation');
end

function testWhitenedResidualsPassOnStructuredDesign(tc)
% categories: robustness | tags: diagnostics
% The SHAPE half under structure: after the loop, the per-stream whiteness
% verdicts recorded in the history must accept the whitened residuals --
% colour that survives whitening on a structured design would bias exactly
% the low degrees an analyst trusts most.
    % R5 M4 4c: the joint estimator's convergence GATE is the per-frequency
    % PSD criterion (max|alpha2 - 1| < JointTol), which a per-stream law on
    % few arcs does not meet in 6 passes; the physically meaningful claim is
    % that the stream LEVELS settled -- the last two passes agree to 2 %
    R = tc.TestData.R;
    h = R.history;
    tc.verifyTrue(~isempty(h) && size(h.sigma, 1) >= 2, 'history must be recorded');
    % tol: 2e-2 relative -- level factors of O(1) streams after 6 passes
    tc.verifyEqual(h.sigma(end, :), h.sigma(end-1, :), 'RelTol', 2e-2, 'the stream levels settled');
    tc.verifyTrue(all(h.whitePval(end, :) >= 0 & h.whitePval(end, :) <= 1));
end

%% ---------------------------------------------------------------- helpers
function cfg = oracleCfg()
%ORACLECFG  Deterministic truth + orbit + noise laws for the oracle.
    cfg.nmax = 8;
    cfg.days = 1:4;
    cfg.arcsPerDay = 2;
    cfg.mEpoch = 150;
    cfg.dt = 5;
    cfg.r  = 1.07;                                   % ~450 km altitude
    cfg.incl = deg2rad(89);                          % near-polar: real gap
    % deterministic truth: closed formula, O(1) unitless scale, plus one
    % planted anomaly at (5,3,C)
    probe = struct('u', 0, 'lam0', 0, 'incl', cfg.incl, 'r', cfg.r);
    [~, meta] = demo_gravity_design(cfg.nmax, probe);
    cfg.meta = meta;
    xT = (cos(3*meta.n + meta.m) ./ meta.n.^2).';
    isS = meta.cs == 'S';
    xT(isS) = (sin(2*meta.n(isS) - meta.m(isS)) ./ meta.n(isS).^2).';
    anom = meta.n == 5 & meta.m == 3 & meta.cs == 'C';
    xT(anom) = xT(anom) + 0.15;                      % the planted anomaly
    cfg.xTrue = xT;
    cfg.nc = numel(xT);
    % per-arc locals: bias + drift, block-diagonal per arc
    cfg.nlPerArc = 2;
    cfg.nl = cfg.nlPerArc * cfg.arcsPerDay * numel(cfg.days);
    cfg.yTrue.A = 0.02 * sin(1:cfg.nl).';
    cfg.yTrue.B = 0.02 * cos(1:cfg.nl).';
    % noise laws (deterministic parameters; realisations seeded)
    cfg.trueSigma = struct('potA', 1.0, 'potB', 2.5);
    cfg.trueAR    = struct('potA', [1 -0.6 0.2], 'potB', [1 -0.3]);
    cfg.cov.potA  = struct('type', 'diag', 'sigma', 1);
    cfg.cov.potB  = struct('type', 'diag', 'sigma', 1);
    cfg.seed = struct('potA', 480001, 'potB', 480002);
end

function A = designOfDay(cfg, sat, day)
%DESIGNOFDAY  The deterministic day design for one satellite (no noise).
    A = demo_gravity_design(cfg.nmax, dayOrbit(cfg, sat, day));
end

function orbit = dayOrbit(cfg, sat, day)
%DAYORBIT  Deterministic orbit samples for one satellite-day: circular,
%   inclined, Earth-rotating node; phase offsets shift the ground track per
%   day and per satellite (satellite B trails A by half an arc).
    M = cfg.arcsPerDay * cfg.mEpoch;
    k = (0:M-1).';
    du = 2*pi * 1.15 / cfg.mEpoch;                   % 1.15 rev per arc
    phase = 0.7 * day + 1.9 * double(sat == 'B');
    orbit = struct('u',    phase + du * k, ...
                   'lam0', 0.4 * day - (du / 15.3) * k, ...
                   'incl', cfg.incl, 'r', cfg.r);
end

function L = oracleLoader(cfg, sat)
%ORACLELOADER  makeLoader wrapper for one satellite system.
    ot = ['pot' sat];
    cov = struct(ot, cfg.cov.(ot));
    L = lsqvce.engine.makeLoader(@(d, o) oracleBuildDay(cfg, sat, d, o), ...
            'observationTypes', {ot}, 'covariance', cov, ...
            'arcs', cfg.arcsPerDay, 'systemId', 1 + double(sat == 'B'), ...
            'numLocal', cfg.nl, 'numCommon', cfg.nc);
end

function D = oracleBuildDay(cfg, sat, day, obsType)
%ORACLEBUILDDAY  One satellite-day: deterministic physical design +
%   deterministic truth + seeded coloured-noise realisation.
    M  = cfg.arcsPerDay * cfg.mEpoch;
    A  = designOfDay(cfg, sat, day);
    B  = zeros(M, cfg.nl);
    tD = (0:M-1).' * cfg.dt;
    for a = 1:cfg.arcsPerDay
        rows = (a-1)*cfg.mEpoch + (1:cfg.mEpoch);
        arcGlobal = (day - 1) * cfg.arcsPerDay + a;
        c0 = (arcGlobal - 1) * cfg.nlPerArc;
        B(rows, c0 + 1) = 1;                                   % bias
        B(rows, c0 + 2) = (0:cfg.mEpoch-1).' / cfg.mEpoch;     % drift
    end
    rs = RandStream('mt19937ar', 'Seed', cfg.seed.(obsType) + 97 * day);
    e  = zeros(M, 1);
    for a = 1:cfg.arcsPerDay
        rows = (a-1)*cfg.mEpoch + (1:cfg.mEpoch);
        w = cfg.trueSigma.(obsType) * randn(rs, cfg.mEpoch, 1);
        e(rows) = filter(1, cfg.trueAR.(obsType), w);
    end
    D.A = A;  D.B = B;
    D.l = A * cfg.xTrue + B * cfg.yTrue.(sat) + e;
    D.t = tD;
    D.arcIdx = arrayfun(@(a) (a-1)*cfg.mEpoch + (1:cfg.mEpoch), ...
                        1:cfg.arcsPerDay, 'UniformOutput', false);
    D = arcCells(D);                                 % v6: the cells are the arcs
end

%% ------------------------------------------- WP-51: the demo is executable
function testDemoGravityFixtureStructure(tc)
% categories: contract
% WP-51: the three things demo_gravity_dataLoader exists to exercise, asserted
% STRUCTURALLY (exact by construction, so they can be pinned without a
% measurement): a different local count per system, a ragged arc schedule
% with arcs absent BY ID, and a different length for every arc.
    [cfg, L] = fixture_gravity();
    tc.verifyNumElements(L, 3);

    % 1. different local counts per system
    tc.verifyEqual([cfg.sys.nlPerArc], [2 3 1]);
    tc.verifyEqual([cfg.sys.nl], [2 3 1] * cfg.nArcTot);
    tc.verifyEqual(numel(unique([cfg.sys.nl])), 3, ...
        'the three systems must carry DIFFERENT local counts');

    % 2. ragged schedule: arcs missing by ID, not merely by count
    counts = cellfun(@numel, cfg.arcsOfDay);
    tc.verifyEqual(counts, [3 4 2 4 3]);
    tc.verifyEqual(cfg.arcsOfDay{3}, [1 3], ...
        'day 3 must be missing arc 2 BY ID');
    tc.verifyEqual(cfg.arcsOfDay{5}, [2 3 4]);
    tc.verifyEqual(cfg.nArcTot, sum(counts));

    % 3. every arc its own length
    lens = cfg.mOfArc(cfg.mOfArc > 0).';
    tc.verifyEqual(numel(lens), cfg.nArcTot);
    tc.verifyGreaterThan(numel(unique(lens)), 5, ...
        'arc lengths must genuinely vary');

    % the loaders answer the PER-DAY meta probe with that day's arcs -- this
    % is what keeps the engines from ever asking for an absent arc
    for s = 1:3
        m3 = L{s}.closure('meta', 3, [], struct());
        tc.verifyEqual(m3.arcs{1}, [1 3]);
        m5 = L{s}.closure('meta', 5, [], struct());
        tc.verifyEqual(m5.arcs{1}, [2 3 4]);
    end
end

function testDemoGravityRunsEndToEnd(tc)
% categories: robustness, reference
% WP-51: the driver must actually RUN -- a demo nobody executes rots. WP-145:
% and it must run AS SHIPPED, to convergence. This test used MaxPasses 1 for
% cost, so it could never see that the demo's own run stopped on its budget
% (the plain 60-bin law needed 72 passes; HANDOVER 12.78). The default run
% settles in 4 passes (~1.6 s); its statistics are checked against bands set
% by their DISTRIBUTIONS, not by the measured values.
    out = evalc('[Rg, sc] = demo_gravity();');
    echoWarnings(out);   % WP-145: its warnings reach the gate
    tc.verifyTrue(Rg.converged, sprintf('the shipped run must converge (stopReason %s)', Rg.stopReason));
    tc.verifyEqual(string(Rg.stopReason), "converged");
    % tol: mean z^2 over 77 coefficients has sd sqrt(2/77) = 0.161 if the
    % formal errors are right; 3 sd either side (measured 0.850)
    tc.verifyLessThan(abs(sc.zMeanSq - 1), 3 * sqrt(2/77), sprintf('mean z^2 = %.3f', sc.zMeanSq));
    % tol: P(max |z| > 4) over 77 independent N(0,1) is 0.5 %, and the help
    % calls 5 alarming (measured 2.32)
    tc.verifyLessThan(sc.zMax, 4, sprintf('max |z| = %.2f', sc.zMax));
    % tol: the planted anomaly within 3 of its own formal sigma (measured 1.10)
    tc.verifyLessThan(abs(sc.anomaly.z), 3, sprintf('anomaly z = %.2f', sc.anomaly.z));
    % WP-150 (F7): sol.Cxx is sigma0^2-scaled already; the demo multiplied the
    % formal error by sigma0 a second time
    cfg = demo_gravity_dataLoader();
    % tol: the same arithmetic on the same Results -- eps-level
    tc.verifyEqual(sc.zMax, max(abs((Rg.x - cfg.xTrue) ./ sqrt(diag(Rg.sol.Cxx)))), 'RelTol', 1e-12);
    tc.verifyTrue(contains(out, 'The problem'));
    tc.verifyTrue(contains(out, 'The weights'));
    tc.verifyClass(Rg, 'lsqvce.Results');
    tc.verifyTrue(all(isfinite(Rg.x)), 'the solution must be finite');
    tc.verifyTrue(isfinite(Rg.sigma0) && Rg.sigma0 > 0);
    tc.verifyNumElements(Rg.x, 77);

    % WP-53: the printed self-assessment is also RETURNED, so it can be
    % checked instead of read. The structural contract; the statistics are
    % banded above, their measured values quoted in the demo's ACCURACY block.
    tc.verifyClass(sc, 'struct');
    documented = {'sigma0','stream','sigmaPlanted','sigmaRecovered', ...
                  'sigmaRatio','zMax','zMeanSq','anomaly','degree','degRatio'};
    tc.verifyEmpty(setdiff(documented, fieldnames(sc)), ...
        'score must carry at least the documented fields');
    tc.verifyEqual(sc.sigma0, Rg.sigma0, ...
        'the score must report the SAME sigma0 the result carries');
    tc.verifyEqual(numel(sc.sigmaRatio), numel(sc.stream));
    tc.verifyEqual(sc.sigmaRatio, sc.sigmaRecovered ./ sc.sigmaPlanted, ...
        'AbsTol', 1e-12);   % tol: exact arithmetic, eps-level
    tc.verifyEqual(sc.degree, 2:8);
    tc.verifyTrue(all(isfinite(sc.degRatio)) && all(sc.degRatio >= 0));
    % tol: the same arithmetic, recomputed -- eps-level agreement
    tc.verifyEqual(sc.anomaly.z, ...
        (sc.anomaly.est - sc.anomaly.true) / sc.anomaly.sigma, 'AbsTol', 1e-12);
end

function testDemoGravityUsesTheNativeRaggedPath(tc)
% categories: contract | tags: provenance
% WP-52: the demo must ride the ENGINE's ragged support, not its own wiring.
% Two mechanical checks that cannot pass by accident: nl is DERIVED by
% makeLoader (the loader carries the plan, and it agrees with the fixture's
% independent arc-ordinal table), and the source no longer wraps .closure.
    exDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'examples');   % L8
    tc.assumeTrue(exist(fullfile(exDir, 'demo_gravity_dataLoader.m'), 'file') == 2);
    p0 = path;   % R7-45: restore, not rmpath -- setupOnce had put exDir there already
    addpath(exDir);
    tc.addTeardown(@() path(p0));

    [cfg, L] = fixture_gravity();
    for s = 1:numel(L)
        tc.verifyNotEmpty(L{s}.localPlan, ...
            'the loader must carry makeLoader''s local plan');
        tc.verifyEqual(L{s}.localPlan.nArcTot, cfg.nArcTot, ...
            'the plan must number the same arcs the fixture declares');
        tc.verifyEqual(L{s}.nl, cfg.sys(s).nlPerArc * cfg.nArcTot, ...
            'nl must be derived from the schedule');
        % the per-day arc list comes from the engine's own probe
        m = L{s}.closure('meta', cfg.days(3), [], struct());
        tc.verifyEqual(m.arcs{1}, cfg.arcsOfDay{3});
    end
    src = fileread(fullfile(exDir, 'demo_gravity_dataLoader.m'));
    tc.verifyEmpty(regexp(src, '\.closure\s*=', 'once'), ...
        'the demo must not wrap the loader closure any more');
end

%% ------------------------------- WP-122: the variational-equation range-rate template
function testKbrHarmonicsAccelerationAndGradient(tc)
% categories: reference, precision
% template/kbr_variational: the normalized solid harmonics of the Cartesian
% recursion against MATLAB's legendre (Condon-Shortley phase removed) --
% near the pole too; the acceleration of a C00 + C20 field against the
% closed-form J2 acceleration; the analytic gradient against central
% differences of the acceleration.
    kit = kbrKit();
    GM = 3.986004415e14;  R = 6378136.3;
    pts = [6.8e6 1.2e6 -2.0e6; 10 -20 6.87e6; -3.1e6 4.4e6 4.9e6].';
    L = 20;  worst = 0;
    for p = pts
        [V, W] = kit.harmonics(p, L, R);
        r = norm(p);  lam = atan2(p(2), p(1));
        for n = 0:L
            P = legendre(n, p(3) / r);
            for m = 0:n
                Pb = (-1)^m * P(m + 1) * sqrt((2 - (m == 0)) * (2*n + 1) * factorial(n - m) / factorial(n + m));
                ref = max(abs(V(:)));
                worst = max([worst, abs(V(n+1, m+1) - (R/r)^(n+1) * Pb * cos(m*lam)) / ref, ...
                                    abs(W(n+1, m+1) - (R/r)^(n+1) * Pb * sin(m*lam)) / ref]);
            end
        end
    end
    % tol: 1e-9 -- the recursion vs legendre + factorials (mirror vs scipy: 7e-11)
    tc.verifyLessThan(worst, 1e-9, 'harmonics vs legendre');
    [nn, mm, cs] = kit.coeffIndex(4);
    c = zeros(25, 1);  c(nn == 0) = 1;  C20 = -4.841694573e-4;  c(nn == 2 & mm == 0 & cs == 0) = C20;
    F = kit.makeField(c, 4, find(nn == 2 & mm == 0 & cs == 0), GM, R);
    for p = pts
        [a, G] = kit.fieldEval(F, p);
        r = norm(p);  J2 = -sqrt(5) * C20;  f = 1.5 * J2 * GM * R^2 / r^5;  z2 = 5 * p(3)^2 / r^2;
        aJ2 = -GM * p / r^3 + f * [p(1) * (z2 - 1); p(2) * (z2 - 1); p(3) * (z2 - 3)];
        % tol: 1e-13 relative -- two closed forms of one field (mirror: 2e-16)
        tc.verifyEqual(a, aJ2, 'RelTol', 1e-13, 'acceleration vs closed-form J2');
        Gfd = zeros(3);
        for j = 1:3
            e = zeros(3, 1);  e(j) = 1;
            Gfd(:, j) = (kit.fieldEval(F, p + e) - kit.fieldEval(F, p - e)) / 2;
        end
        % tol: 1e-7 of |G| -- central differences at 1 m (mirror: 8e-10)
        tc.verifyLessThan(max(abs(G(:) - Gfd(:))) / max(abs(G(:))), 1e-7, 'gradient vs differences');
    end
end

function testKbrIntegratorAgainstKepler(tc)
% categories: reference, precision
% The ABM8 PECE integrator of the template on a point-mass field, six hours
% at 5 s, against Kepler's equation.
    kit = kbrKit();
    GM = 3.986004415e14;  R = 6378136.3;
    [nn, ~, ~] = kit.coeffIndex(2);
    c = double(nn == 0);
    F = kit.makeField(c, 2, zeros(0, 1), GM, R);
    y0 = kbrLeo(GM, R, 0.4);
    st = kit.propagate(F, @(t) eye(3), @(t) deal(eye(3), zeros(3, 1)), zeros(3, 1), ones(3, 1), 0, y0, 5, 4320);
    rk = keplerPosition(GM, y0(1:3), y0(4:6), 5 * 4320);
    % tol: 1e-4 m after 6 h -- ABM8 at 5 s (mirror: 6e-6 m after a full day)
    tc.verifyLessThan(norm(st(end, 1:3).' - rk), 1e-4, 'integrator vs Kepler');
end

function testKbrPartialsAgainstReintegration(tc)
% categories: precision, metamorphic
% The variational equations: Psi times a small change of the initial state,
% of one coefficient, of the ACC bias and scale equals the change of the
% re-integrated state (1.5 h, degree 6, attitude and ACC signal on).
    kit = kbrKit();
    GM = 3.986004415e14;  R = 6378136.3;  N = 6;
    [nn, mm, cs] = kit.coeffIndex(N);
    c = zeros(numel(nn), 1);  c(nn == 0) = 1;
    c(nn == 2 & mm == 0 & cs == 0) = -4.841694573e-4;
    rs = RandStream('mt19937ar', 'Seed', 3);
    k3 = nn >= 3;  c(k3) = 1e-6 * randn(rs, nnz(k3), 1) ./ nn(k3).^2;
    est = find(nn >= 2);  np = numel(est);
    y0 = kbrLeo(GM, R, 0.4);
    run = @(cc, yy, b, s) kit.propagate(kit.makeField(cc, N, est, GM, R), @(t) eye(3), @kbrAcc, b, s, 0, yy, 5, 1080);
    [Y, Psi] = run(c, y0, zeros(3, 1), ones(3, 1));
    j = find(est == find(nn == 5 & mm == 3 & cs == 1));
    cases = {'IC x',    1,         @() run(c, y0 + [1; 0; 0; 0; 0; 0], zeros(3, 1), ones(3, 1)), 1;
             'IC vy',   5,         @() run(c, y0 + [0; 0; 0; 0; 1e-3; 0], zeros(3, 1), ones(3, 1)), 1e-3;
             'S53',     6 + j,     @() run(c + 1e-7 * ((1:numel(c)).' == est(j)), y0, zeros(3, 1), ones(3, 1)), 1e-7;
             'bias y',  6 + np + 2, @() run(c, y0, [0; 1e-8; 0], ones(3, 1)), 1e-8;
             'scale x', 9 + np + 1, @() run(c, y0, zeros(3, 1), [1.1; 1; 1]), 0.1};
    for k = 1:size(cases, 1)
        Y2 = cases{k, 3}();
        d = (Y2(end, :) - Y(end, :)).';
        col = cases{k, 2};             % Psi columns: [r0 v0 | SH | bias | scale]
        % tol: 1e-4 relative -- second-order terms and rounding (mirror: 2e-7..3e-6).
        % The steps move the final state by >= 5 cm: the S53 step was 1e-9 in CI
        % round 1, its 4.7 mm met the ~5e-7 m rounding floor at 1.09e-4 (WP-122)
        tc.verifyLessThan(norm(Psi(:, col) * cases{k, 4} - d) / norm(d), 1e-4, cases{k, 1});
    end
end

function testKbrClosedLoopThroughFit(tc)
% categories: metamorphic, contract
% The whole template on one synthetic day of GRACE-FO L1B files (the real
% ASCII layouts, written by fixture_kbr_l1b): a planted correction of
% degrees 2..6 and a 1-3 % ACC scale error, two 12 h arcs, KBR + both
% positions through lsqvce.fit -- the field comes back. And the template's
% own contract: a missing required field and a missing file are refused by name.
    dataDir = tempname;  mkdir(dataDir);
    tc.addTeardown(@() rmdir(dataDir, 's'));
    [cfg, truth] = fixture_kbr_l1b(dataDir);
    [task, kit, data] = kbr_variational(cfg);      % WP-161: the task, and its compiled struct
    tc.verifyEqual(size(data.commonParams, 1), numel(truth.dc), 'degrees 2..6: 45 coefficients');
    R = lsqvce.fit(task, Verbose=false);
    dc = kit.unitCommon .* R.x(:);
    % tol: 1e-3 of |dc| -- one linearization at the a-priori ACC calibration
    % (mirror, one day, degree 8, the same planted errors: 3.5e-5)
    tc.verifyLessThan(norm(dc - truth.dc) / norm(truth.dc), 1e-3, 'the planted field comes back');
    tc.verifyError(@() kbr_variational(rmfield(cfg, 'trf2crf')), 'kbrvar:badConfig');
    tc.verifyError(@() kit.readL1B(fullfile(dataDir, 'none.txt'), struct('fmt', '%f', 't', 1)), 'kbrvar:noFile');
end

function kit = kbrKit()
%KBRKIT  The template's physics kit (a declaration that reads no data file).
    gf = struct('GM', 3.986004415e14, 'R', 6378136.3, 'C', blkdiag(1, zeros(2)), 'S', zeros(3));
    [~, kit] = kbr_variational(struct('days', 60000, 'dataDir', tempdir, 'gravityField', gf, ...
                                      'trf2crf', @(t) eye(3), 'degree', [2 2]));
end

function y0 = kbrLeo(GM, R, u)
%KBRLEO  A circular 490 km, 89 deg orbit at argument of latitude u.
    a = R + 490e3;  vc = sqrt(GM / a);  i = deg2rad(89);
    y0 = [a * [cos(u); sin(u) * cos(i); sin(u) * sin(i)]; vc * [-sin(u); cos(u) * cos(i); cos(u) * sin(i)]];
end

function [Q, acc] = kbrAcc(t)
%KBRACC  A turning attitude and a drag-like ACC signal (SRF, m/s^2).
    th = 2 * pi / 5650 * t + 0.2;
    Q = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1] * [1 0 0; 0 cos(0.3) -sin(0.3); 0 sin(0.3) cos(0.3)];
    acc = [-1e-7 * (1 + 0.3 * cos(th)); 2e-8 * sin(th); 5e-8 * cos(2 * th)];
end

function r = keplerPosition(GM, r0, v0, dt)
%KEPLERPOSITION  The oracle: two-body position after dt by Kepler's equation.
    rn = norm(r0);  a = 1 / (2 / rn - (v0.' * v0) / GM);  n = sqrt(GM / a^3);
    ecosE = 1 - rn / a;  esinE = (r0.' * v0) / sqrt(GM * a);
    E0 = atan2(esinE, ecosE);  e = hypot(ecosE, esinE);
    M = E0 - esinE + n * dt;  E = M;
    for k = 1:50, E = E - (E - e * sin(E) - M) / (1 - e * cos(E)); end
    dE = E - E0;
    r = (1 - a / rn * (1 - cos(dE))) * r0 + (dt - (dE - sin(dE)) / n) * v0;
end

%% ------------------------------- WP-144: the GRAVFIRE template on a generated bundle
function testGravfireTemplateRunsOnAGeneratedBundle(tc)
% categories: contract, robustness
% WP-144 (A7): template/gravfire_data.m ran in no test -- its only caller,
% demo_gravfire, is given a path that cannot exist so that every tier runs the
% same suite, and the user's bundle is measurement data that never enters the
% repository. FIXTURE_GRAVFIRE_BUNDLE generates one with the real structure;
% the whole template runs on it: the tables, column 11 -> priors, the
% epoch-major rows with Cep, the flagged epoch, the cut into pieces, the fit.
    [c, m, s] = fixture_gravfire_bundle();
    [task, data] = gravfire_data(c, m, s, 'DryRun', true);   % WP-161: the task, and its compiled struct
    % WP-161: the compiled tables hold scope and name; the 14 GRAVFIRE columns
    % travel as each template's Info
    tc.verifySize(data.commonParams, [6 2], 'the scope token and the unique name');
    tc.verifyEqual(data.days, [59396 59397]);
    tc.verifyEqual(fieldnames(data.systems.covariance).', {'posA', 'posB', 'kbr'});
    tc.verifySize(data.systems.localParams, [6 2], 'one day''s templates: 2 satellites x 3');
    tc.verifyTrue(all(strcmp(data.systems.localParams(:, 1), 'day')));
    % column 11: the labelled common group 'stokes' (rows 1-2, one sigma
    % estimated, x0 per member) and the unlabelled row 3 (fixed)
    tc.verifyEqual({data.priors.name}, {'stokes', 's_2_1'});
    tc.verifyEqual(data.priors(1).L, [eye(2), zeros(2, 4)]);
    tc.verifyEqual(data.priors(1).x0, [0.5; 0], 'each member keeps its a-priori value');
    tc.verifyEqual(string({data.priors.estimate}), ["scale", "none"]);
    tc.verifyEqual(data.systems.priors(1).group, '^(GraceA_accscale|GraceB_accscale)$');
    % the labelled groups are "scale" and this single solve holds them: the
    % warning is the truth here (WP-150), fired on purpose
    muteWarn(tc, 'lsqvce:fit:priorScaleHeld');
    R = lsqvce.fit(task, lsqvce.Options(Verbose=0));
    % rows: 2 days x (posA 59 x 3 + posB 59 x 3 + kbr 59) -- 2 pieces of 30
    % epochs per (day, group), one flagged epoch at weight 0 (WP-148); unknowns: 6 commons
    % + 6 locals per day; prior rows: 'stokes' 2 + row 3 1 on the commons,
    % 'accscale' 2 members x 2 days on the locals
    tc.verifyEqual(R.sol.dof, 2*(59*3 + 59*3 + 59) - (6 + 2*6) + (2 + 1) + 2*2, ...
        'every row of every piece reached the accumulator');
    tc.verifyNumElements(R.x, 6);
end

function testGravfireFlaggedEpochsStayOnTheGrid(tc)
% categories: contract, robustness
% WP-148 (HANDOVER 12.80, observation 2): the template REMOVED the flagged
% epochs, which left a gap inside the piece; on a noise run the 'psd' law
% took it for one step (lsqvce:covModel:unevenGrid). The epochs stay on the
% grid with weight 0 now (D.excluded): the grid is even, the dof is the one
% the removal gave, and the flagged epoch's values -- NaN in a real bundle --
% do not reach the solution.
    [c, m, s] = fixture_gravfire_bundle();
    [task, data] = gravfire_data(c, m, s, 'DryRun', true);   % WP-161: the task, and its compiled struct
    po = data.systems.processingOptions;  info = struct('nc', 6, 'nl', 12);
    D = data.systems.buildDay(59396, 'posA', [], [], po, info);
    tc.verifyEqual(cellfun(@numel, D.l), [90 90], 'every epoch of both pieces, 3 rows each');
    tc.verifyEqual(cellfun(@nnz, D.excluded), [0 3], 'the flagged epoch (piece 2) weighted out');
    tc.verifyEqual(unique(diff(D.t{2}(1:3:end))), 5, 'no gap inside the piece');
    % its values do not matter: NaN observations and covariance give the same fit
    s2 = s;
    for k = 1:numel(s2{1}.observations)
        o = s2{1}.observations{k};  f = o.nanIndexOfPiece{1};
        o.lobs{1}(f, :) = NaN;
        if isfield(o, 'iSpp'), o.iSpp{1}(f, :) = NaN; end
        s2{1}.observations{k} = o;
    end
    muteWarn(tc, 'lsqvce:fit:priorScaleHeld');      % the labelled groups, held by a single solve
    R1 = lsqvce.fit(task, lsqvce.Options(Verbose=0));
    R2 = lsqvce.fit(gravfire_data(c, m, s2, 'DryRun', true), lsqvce.Options(Verbose=0));
    % tol: 0 -- the flagged rows are zeroed before any arithmetic; the same doubles
    tc.verifyEqual(R2.x, R1.x, 'the flagged values never reach the solution');
    tc.verifyEqual(R1.sol.dof, 2*(59*3 + 59*3 + 59) - (6 + 2*6) + (2 + 1) + 2*2, ...
        'the dof the removal gave: a weight-0 epoch is no observation');
    % the noise run: the grid is even, the gap warning does not fire
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');      % 3 passes on 2 days stop early
    lsqvce.util.log.inventory('reset');
    lsqvce.fit(task, lsqvce.Options(Verbose=0, Estimate="noise", MaxPasses=3));
    tc.verifyFalse(any(strcmp(lsqvce.util.log.inventory(), 'lsqvce:covModel:unevenGrid')), ...
        'no gap inside a piece');
    lsqvce.util.log.inventory('reset');
end

function testGravfireLocalGroupRefusesDisagreeingMembers(tc)
% categories: contract
% WP-144 (A4): the members of a labelled LOCAL Tikhonov group were checked
% for sigma only; a later member's fromPass and a-priori value were dropped
% and the group kept the first member's (reproduced on the generated bundle).
% Both are refused now, as the common rows refuse a differing fromPass.
    [c, m, s] = fixture_gravfire_bundle();
    P = s{1}.parameter;
    k = find(strcmp(P(:, 1), 'GraceB_accscale')).';
    s2 = s;  Q = P;
    for r = k, Q{r, 11} = {'Tikhonov', 3, 1e-6, 'accscale'}; end
    s2{1}.parameter = Q;
    tc.verifyError(@() gravfire_data(c, m, s2, 'DryRun', true), 'gravfire:badConstraint');
    s3 = s;  Q = P;
    for r = k, Q{r, 10} = 0.25; end
    s3{1}.parameter = Q;
    tc.verifyError(@() gravfire_data(c, m, s3, 'DryRun', true), 'gravfire:badConstraint');
end
