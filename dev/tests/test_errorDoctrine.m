function tests = test_errorDoctrine
%TEST_ERRORDOCTRINE  R5 M1 (WP-80): the error doctrine's mechanisms.
%   tests = test_errorDoctrine
%
%   Pins what the doctrine (dev/lsqvce_theory_r5.pdf ch. 7) turned into code:
%     - lsqvce.util.callUser: outputs counted ONCE, user errors wrapped with
%       the original as cause (finding H-08);
%     - lsqvce.util.log inventory / warnOnce: warnings are a RESULT
%       (R.meta.warnings), and an identifier fires once per run;
%     - Results.stopReason: convergence failure is a result, not an error;
%     - Options.SkipThreshold: a typed, range-checked run policy;
%     - BlockAccumulator.windowMax: the covariance grid is sized from a
%       MEASURED window (finding H-06), never probed;
%     - the joint (Ellmer) estimator REFUSES a multi-axis stream (H-15)
%       instead of computing the wrong estimator slowly.
%   The skip-rate refusal (lsqvce:adjustVCE:skipThreshold) has no provoking
%   fixture yet -- recorded as the M1 gap, closed with the M3 noise law.
%
%   See also LSQVCE.UTIL.CALLUSER, LSQVCE.UTIL.LOG, LSQVCE.RESULTS, TEST_CI_GATE
%
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M1, WP-80.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    addpath(root);
    addpath(fullfile(root, 'examples'));
    tc.TestData.root = root;
end

%% ============================================================ callUser
function testCallUserRejectsNonHandle(tc)
% categories: contract
    verifyError(tc, @() lsqvce.util.callUser(42, 1, 1, 'ctx', 1), 'lsqvce:callUser:badHandle');
end

function testCallUserTooFewOutputsDeclared(tc)
% categories: contract
    % a NAMED function with a known nargout: refused without any call
    verifyError(tc, @() lsqvce.util.callUser(@sin, 2, 3, 'ctx', 1), 'lsqvce:callUser:tooFew');
end

function testCallUserProbesOnceAndReturnsTheOutputs(tc)
% categories: contract
    % anonymous handles report nargout = -1, so the count is found by ONE
    % trial that must not be wasted: the outputs of that call come back
    f2 = @(k) deal(k, 2*k);
    f3 = @(k) deal(k, 2*k, 3*k);
    [n2, a, b, c] = lsqvce.util.callUser(f2, 2, 3, 'f2', 5);
    verifyEqual(tc, n2, 2);  verifyEqual(tc, [a b], [5 10]);  verifyEmpty(tc, c);
    [n3, a, b, c] = lsqvce.util.callUser(f3, 2, 3, 'f3', 5);
    verifyEqual(tc, n3, 3);  verifyEqual(tc, [a b c], [5 10 15]);
    % fixed count: min == max, no probing path
    [n, a, b] = lsqvce.util.callUser(f2, 2, 2, 'f2', 7);
    verifyEqual(tc, n, 2);  verifyEqual(tc, [a b], [7 14]);
end

function testCallUserZeroOutputHook(tc)
% categories: contract
    % a side-effect-only hook (assigns nothing): asked for one output it
    % raises TooManyOutputs ONCE, is then called for its effect, figH = []
    [n, figH] = lsqvce.util.callUser(@(ctx) disp(''), 0, 1, 'hook', struct());
    verifyEqual(tc, n, 0);
    verifyEmpty(tc, figH);
end

function testCallUserWrapsUserErrorWithCause(tc)
% categories: robustness
    % a NAMED function with an output that throws inside (an anonymous
    % handle whose body is error(...) never runs: MATLAB raises TooManyOutputs
    % for the missing output first -- that case is the tooFew test below)
    try
        lsqvce.util.callUser(@badBatch, 1, 1, 'batchFcn(k)', 3);
        verifyFail(tc, 'must throw');
    catch ME
        verifyEqual(tc, ME.identifier, 'lsqvce:callUser:threw');
        verifyTrue(tc, contains(ME.message, 'batchFcn(k)'), 'the context is named');
        verifyEqual(tc, numel(ME.cause), 1);
        verifyEqual(tc, ME.cause{1}.identifier, 'mine:boom', 'the original is the cause, unmodified');
        verifyEqual(tc, ME.cause{1}.message, 'the user said 3');
    end
end

function A = badBatch(k)
    A = [];  %#ok<NASGU>
    error('mine:boom', 'the user said %d', k);
end

function testCallUserArityErrorFromInsideIsNotMisread(tc)
% categories: robustness
    % a user function that ITSELF raises TooManyOutputs (a bug in it) must
    % not be misread as "fewer outputs": below minOut the probe refuses
    liar = @(k) error('MATLAB:TooManyOutputs', 'inner');
    verifyError(tc, @() lsqvce.util.callUser(liar, 1, 1, 'ctx', 1), 'lsqvce:callUser:tooFew');
end

%% ============================================================ inventory
function testInventoryCollectsEachIdOnce(tc)
% categories: contract | tags: provenance
    lsqvce.util.log.inventory('reset');
    verifyEmpty(tc, lsqvce.util.log.inventory());
    st = warning('off', 'lsqvce:internal:testProbe');  c1 = onCleanup(@() warning(st)); %#ok<NASGU>
    st2 = warning('off', 'lsqvce:internal:testProbe2'); c2 = onCleanup(@() warning(st2)); %#ok<NASGU>
    lsqvce.util.log.warn('lsqvce:internal:testProbe', 'one');
    lsqvce.util.log.warn('lsqvce:internal:testProbe', 'again');
    lsqvce.util.log.warn('lsqvce:internal:testProbe2', 'two');
    inv = lsqvce.util.log.inventory();
    verifyEqual(tc, inv, {'lsqvce:internal:testProbe', 'lsqvce:internal:testProbe2'}, ...
        'first-fired order, each once');
    lsqvce.util.log.inventory('reset');
    verifyEmpty(tc, lsqvce.util.log.inventory());
end

function testWarnOnceFiresOncePerRun(tc)
% categories: contract
    lsqvce.util.log.inventory('reset');
    st = warning('off', 'lsqvce:internal:testProbe');  c1 = onCleanup(@() warning(st)); %#ok<NASGU>
    fired = false(1, 5);
    for k = 1:5
        fired(k) = lsqvce.util.log.warnOnce('lsqvce:internal:testProbe', 'k = %d', k);
    end
    verifyEqual(tc, fired, [true false false false false]);
    verifyEqual(tc, numel(lsqvce.util.log.inventory()), 1);
    lsqvce.util.log.inventory('reset');
end

function testFitStoresTheInventoryInMeta(tc)
% categories: contract | tags: provenance
    % a clean single solve: fit resets the inventory and stores it (empty)
    rs = RandStream('mt19937ar', 'Seed', 3);
    A  = randn(rs, 60, 4);  x = [1; -2; 0.5; 3];
    l  = A*x + 0.01*randn(rs, 60, 1);
    R  = lsqvce.fit(struct('A', A, 'l', l), lsqvce.Options(Verbose=0));
    verifyTrue(tc, isfield(R.meta, 'warnings'));
    verifyTrue(tc, iscell(R.meta.warnings));
    verifyEmpty(tc, R.meta.warnings, 'a clean solve fires no lsqvce warning');
    verifyEqual(tc, R.stopReason, "single");
end

%% ============================================================ stopReason
function testStopReasonIsDerived(tc)
% categories: contract
    rs = RandStream('mt19937ar', 'Seed', 4);
    A  = randn(rs, 50, 3);  l = A*[1; 2; 3] + 0.01*randn(rs, 50, 1);
    R  = lsqvce.fit(struct('A', A, 'l', l), lsqvce.Options(Verbose=0));
    verifyEqual(tc, R.stopReason, "single");
    R.setIterative(struct('pass', 1), struct(), true);
    verifyEqual(tc, R.stopReason, "converged");
    R.setIterative(struct('pass', 1), struct(), false);
    verifyEqual(tc, R.stopReason, "maxPasses");
    R.setIterative(struct([]), struct(), false);
    verifyEqual(tc, R.stopReason, "single");
end

%% ============================================================ SkipThreshold
function testSkipThresholdOption(tc)
% categories: contract
    % R5 M5: the skip threshold lives in the Noise group; the default is 0.2
    % and it is forwarded to the joint driver, which range-checks it
    o = lsqvce.Options();
    verifyEqual(tc, o.Noise.SkipThreshold, 0.2);
    o = lsqvce.Options(Noise=struct('SkipThreshold', 0.5));
    a = o.engineArgs("adjustVCEjoint");
    verifyEqual(tc, a.nv{find(strcmp(a.nv, 'SkipThreshold')) + 1}, 0.5);
end

%% ============================================================ measured window
function testWindowMaxIsMeasuredFromTheAccumulation(tc)
% categories: contract, reference
    % nothing accumulated -> 0; after a system: min(BatchSize + overlap, rows
    % of the largest (day, obsType)) -- with a diag law the overlap is 0
    cfg = makeCfg(6, 3, 40, 2);                       % nc, arcs/day, rows/arc, nl
    f   = @(d, ot, a, uo) lsqvce.engine.synthLoader(cfg, 1, d, ot, a, uo);
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose', false, 'BatchSize', 70);
    verifyEqual(tc, acc.windowMax, 0);
    acc.addSystem(1, cfg.nl, f, 1:2);
    verifyEqual(tc, acc.windowMax, 70, 'BatchSize caps the window (day rows = 120)');
    acc2 = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose', false);   % default 20000
    acc2.addSystem(1, cfg.nl, f, 1:2);
    verifyEqual(tc, acc2.windowMax, 120, 'the day''s rows cap the window');
end

function cfg = makeCfg(nc, arcsPerDay, mArc, nl)
    rs = RandStream('mt19937ar', 'Seed', 42);
    cfg.nc = nc;  cfg.nl = nl;  cfg.sats = 1;
    cfg.x_c = randn(rs, nc, 1);
    cfg.y = containers.Map('KeyType', 'double', 'ValueType', 'any');
    cfg.y(1) = randn(rs, nl, 1);
    cfg.obsTypes = {'kbr'};
    cfg.arcsPerDay = arcsPerDay;  cfg.mArc = mArc;  cfg.dt = 5;  cfg.seed = 42;
    cfg.collinear = false;  cfg.collinearEps = 1e-7;
    cfg.cov.kbr = struct('type', 'diag', 'sigma', 0.2);
end

%% ============================================================ H-15: multi-axis on the joint path
function testEllmerRunsOnMultiAxisStreams(tc)
% categories: contract, robustness
    % R5 M3 (H-15): the demo's 3-axis position streams are declared
    % 'separable' (AR temporal law + per-epoch Cep). The joint (Ellmer)
    % estimator now descales every epoch and takes the per-lag statistics
    % PER AXIS against the temporal law over EPOCHS -- so it RUNS (M1 refused
    % it) and returns a pass record; the numbers are a smoke check here, the
    % parity against the single-axis reference is test_ellmer's job.
    muteWarn(tc, 'lsqvce:decorrFactor:notPD');
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');
    [~, data] = demo_grace_dataLoader();
    % R5 M3 3b (H-17): the joint law now leaves jointModelOpt in the PSD
    % domain and covModel resamples it to every arc's length, so the KBR arcs
    % longer than NumLags = 16 are no longer non-PD: NO arc is skipped at the
    % DEFAULT threshold (before 3b: 8 of 32, all 's1_rr_*'). The 24 position
    % arcs -- the H-15 case -- are processed per axis.
    opts = lsqvce.Options(Estimate="noise", Noise=struct('Grid', 16), MaxPasses=1, Verbose=0);
    R = lsqvce.fit(data, opts);
    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyEqual(R.history.kind, 'ellmer');
    tc.verifyTrue(all(isfinite(R.x)), 'a finite estimate');
    tc.verifyTrue(all(isfinite(R.noise.C)), 'a finite pooled covariance function');
    tc.verifyEmpty(R.history.skips, 'no arc is skipped: the PSD law is valid for every arc length');
    tc.verifyEqual(numel(R.noise.arcKeys), 32);
    % R5 M4 4a: the joint sweep now COLLECTS whiteness per stream per pass,
    % so the convergence table exists on an ellmer Results (it could not
    % even run before): one row per joint iteration, the streams named,
    % Ljung-Box p-values in [0, 1]
    T = R.convergence();
    tc.verifyEqual(height(T), R.history.nIter);
    tc.verifyEqual(numel(R.history.streams), 4, 'posC, posD, rr of system 1 and posE of system 2');
    tc.verifyTrue(all(contains(R.history.streams, '|')), 'streams are keyed satKey|obsType');
    p = R.history.whitePval;
    tc.verifyTrue(all(p(:) >= 0 & p(:) <= 1), 'p-values in [0, 1]');
    tc.verifyEqual(size(R.history.sigma), size(p));
end
%% ============================================================ R5 M2 2b: IRLS kernel
function testRobustWeightsContract(tc)
% categories: contract
    verifyError(tc, @() lsqvce.kernel.robustWeights([1 2 3], 'cauchy', []), 'lsqvce:robustWeights:badInfluence');
    verifyError(tc, @() lsqvce.kernel.robustWeights([1 2 3], 'igg3', 2.5), 'lsqvce:robustWeights:badTuning');
    verifyError(tc, @() lsqvce.kernel.robustWeights([1 2 3], 'huber', [1 2]), 'lsqvce:robustWeights:badTuning');
end

function testRobustWeightsAgainstTheFormulas(tc)
% categories: reference
    % residuals whose MAD scale is exactly 1: symmetric +-1 around 0 with a
    % large outlier; then each influence's closed form is checked at 3 sigma
    vw = [-1 1 -1 1 -1 1 -1 1 9].';                   % median 1; |vw - med| = 0,2,... -> MAD 2
    [facH, s, w] = lsqvce.kernel.robustWeights(vw, 'huber', []);
    % tol: 1e-14 absolute throughout -- closed forms on 9 numbers of order 1
    % re-evaluated here with the same double arithmetic (a few ulps)
    tc.verifyEqual(s, 1.4826 * median(abs(vw - median(vw))), 'AbsTol', 1e-14);
    tc.verifyEqual(w, vw / s, 'AbsTol', 1e-14);
    expH = ones(9, 1);  big = abs(w) > 1.345;  expH(big) = 1.345 ./ abs(w(big));
    % tol: 1e-14 absolute (see above)
    tc.verifyEqual(facH, expH, 'AbsTol', 1e-14);
    facB = lsqvce.kernel.robustWeights(vw, 'bisquare', []);
    expB = zeros(9, 1);  in = abs(w) < 4.685;  expB(in) = (1 - (w(in)/4.685).^2).^2;
    % tol: 1e-14 absolute (see above)
    tc.verifyEqual(facB, expB, 'AbsTol', 1e-14);
    facI = lsqvce.kernel.robustWeights(vw, 'igg3', [1.5 3]);
    tc.verifyEqual(facI(abs(w) <= 1.5), ones(nnz(abs(w) <= 1.5), 1));
    tc.verifyEqual(facI(abs(w) > 3), zeros(nnz(abs(w) > 3), 1));
end

function testRowWeightKeyFormat(tc)
% categories: contract
    tc.verifyEqual(lsqvce.util.rowWeightKey(1, 59396, 'kbr', 3), '1|59396|kbr|3');
    tc.verifyEqual(lsqvce.util.rowWeightKey('C', 2, "posA", 'a1'), 'C|2|posA|a1');
    tc.verifyNotEqual(lsqvce.util.rowWeightKey(1, 2, 'x', 3), lsqvce.util.rowWeightKey(1, 2, 'x', 4));
end

%% ============================================================ R5 M2 2c': the batch adapter
function testBatchSystemContract(tc)
% categories: contract
    verifyError(tc, @() lsqvce.engine.batchSystem(3, 42, 1), 'lsqvce:batchSystem:badSource');
    rs = RandStream('mt19937ar', 'Seed', 9);
    A = randn(rs, 12, 3);  l = randn(rs, 12, 1);  p = 0.5 + rand(rs, 12, 1);
    L = lsqvce.engine.batchSystem(3, @(k) deal(A, l, p), 1);
    tc.verifyEqual(L.nc, 3);  tc.verifyEqual(L.nl, 0);
    uo = struct('sysID', 1, 'nc', 3, 'nl', 0, 'BatchSize', 20000, 'obsType', 'obs', 'day', 1, 'arc', 1);
    arc = L.closure(1, 'obs', 1, uo);
    % rows are pre-whitened by sqrt(p): exact (diagonal-P GLS identity)
    tc.verifyEqual(arc.A, sqrt(p) .* A);  tc.verifyEqual(arc.l, sqrt(p) .* l);
    % a full weight matrix is a covariance law of the streaming kind: refused
    Lm = lsqvce.engine.batchSystem(3, @(k) deal(A, l, eye(12)), 1);
    verifyError(tc, @() Lm.closure(1, 'obs', 1, uo), 'lsqvce:batchSystem:weightMatrix');
    % a wrong A shape is named
    Lw = lsqvce.engine.batchSystem(4, @(k) deal(A, l), 1);
    verifyError(tc, @() Lw.closure(1, 'obs', 1, uo), 'lsqvce:batchSystem:badMatrix');
end
