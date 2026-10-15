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
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-137 (code review of WP-135):
%   Results.compare reads the solve's rank and redundancy. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (third review, lanes 1 and 4):
%   ftest takes the redundancy; the joint robust report's s0 is R.sigma0.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-134 (third review): robustWeights survives an extreme blunder; tuning validated.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-133: weights >= 0 and finite, a
%   zero weight / a day's .excluded marks its rows; the joint hard solution's
%   residual pass and verdict agree (review 3 findings 14, 17, 18). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M1, WP-80.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    pathBefore = path;   % R7-45: a standalone run (run_lsqvce / runtests) leaves the path as found
    tc.addTeardown(@() path(pathBefore));
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

function testKindBWarningReachesTheInventory(tc)
% categories: contract | tags: provenance
    % R7-23: kinds B / C reset the inventory AFTER runBatch's priorNotActive
    % warnOnce, so the id never reached R.meta.warnings; and once any run had
    % fired it, the next kind-B warnOnce was suppressed before that reset --
    % no warning at all. fit now resets once, before anything can warn.
    muteWarn(tc, 'lsqvce:fit:priorNotActive');   % provoked on purpose (gate rule 7)
    rs = RandStream('mt19937ar', 'Seed', 4);
    A  = randn(rs, 40, 3);  l = A*[1; 2; 3] + 0.01*randn(rs, 40, 1);
    data = struct('A', A, 'l', l, 'priors', struct('L', [0 0 1], 'sigma', 1, 'from', 3));
    R1 = lsqvce.fit(data, Verbose=false);
    verifyTrue(tc, any(strcmp(R1.meta.warnings, 'lsqvce:fit:priorNotActive')), ...
        'the kind-B warning is in its own run''s inventory');
    lastwarn('');
    R2 = lsqvce.fit(data, Verbose=false);         % the id fired in the previous run
    [~, wid] = lastwarn;
    verifyEqual(tc, wid, 'lsqvce:fit:priorNotActive', 'and the next run warns again, once');
    verifyTrue(tc, any(strcmp(R2.meta.warnings, 'lsqvce:fit:priorNotActive')));
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
    % with the per-epoch Cep on the day (R6 WP-90: a 'diag' temporal law + Cep). The joint (Ellmer)
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
function testJointRobustPassKeepsCleanMultiAxisEpochs(tc)
% categories: robustness, precision
    % R7-08 (user decision 2026-10-10: chi-calibrated). The joint loop's robust
    % pass reduced a 3-axis epoch to its RMS -- chi-shaped, median ~1 -- and
    % robustWeights standardises vw/s uncentred, so IGG3 zeroed ~24 % of CLEAN
    % 3-axis epochs (mirror: mirrors/wp115/chi_calibrated.py). The demo's
    % streams carry no blunders: with the calibrated statistic a clean epoch
    % is zeroed about as often as a clean scalar residual (~0.3 %).
    muteWarn(tc, 'lsqvce:decorrFactor:notPD');
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');
    [~, data] = demo_grace_dataLoader();
    R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Robust="igg3", ...
                   Noise=struct('Grid', 16), MaxPasses=2, Verbose=0));
    fac = R.noise.robust.factor;
    % tol: 5 % -- the calibrated statistic zeroes ~0.3 % of clean epochs
    % (scalar-equivalent); the RMS statistic zeroed ~24 % of the 3-axis ones,
    % which make up most of the demo's epochs
    tc.verifyLessThan(mean(fac == 0), 0.05, sprintf( ...
        'clean data: %.1f %% of the epochs zeroed (want < 5 %%)', 100 * mean(fac == 0)));
    % WP-133 (review 3 finding 18): the hard solution's sigma0 carries the
    % truncation correction and its rejected epochs leave the dof; the
    % residual pass applies the same (sol.truncation) and counts the same
    % rows -- it warned sigma0Mismatch on every igg3 / bisquare noise run
    res = R.residuals('Verbose', false);
    tc.verifyFalse(res.supersedesSol, 'the residual pass agrees with the hard solution');
    % tol: the residual pass's own mismatch line (1e-6), not met by rounding
    tc.verifyEqual(res.sigma0, R.sigma0, 'RelTol', 1e-6);
    % (14) the rejected epochs are "rejected" in the verdict, and its sum(r)
    % identity holds: the nugget whitener rotates a deleted epoch onto its own
    % rows (dev/validation/wp133_nugget_rotation_mirror.py)
    tc.verifyEqual(R.diag.quad.rejected > 0, any(fac == 0), 'a zeroed epoch is rejected in R.diag');
    % tol: the report's own integrity line (SUM_R > 1e-6 is INVALID)
    tc.verifyLessThan(R.diag.sumRRelErr, 1e-6);
    % WP-135 (review 3, lane 4): the REPORT's s0 is the hard solution's too --
    % D.s0 read vPv / r without the truncation (sqrt(tau) off R.sigma0, and the
    % chi^2 test and the scaled PAR_TABLE column with it)
    tc.verifyLessThan(R.diag.truncation, 1, 'an igg3 run truncates');
    % tol: the report's vPv integrity line (1e-6 relative; s0 takes its root)
    tc.verifyEqual(R.diag.s0, R.sigma0, 'RelTol', 1e-6);
end

%% ============================================================ R5 M2 2b: IRLS kernel
function testRobustWeightsContract(tc)
% categories: contract
    verifyError(tc, @() lsqvce.kernel.robustWeights([1 2 3], 'cauchy', []), 'lsqvce:robustWeights:badInfluence');
    verifyError(tc, @() lsqvce.kernel.robustWeights([1 2 3], 'igg3', 2.5), 'lsqvce:robustWeights:badTuning');
    verifyError(tc, @() lsqvce.kernel.robustWeights([1 2 3], 'huber', [1 2]), 'lsqvce:robustWeights:badTuning');
end

function testRobustWeightsScaleFallbackIsReported(tc)
% categories: robustness, contract
    % R7-21: the MAD -> std -> 1 fallback was silent, and a NEAR-tie MAD was
    % accepted. Half the residuals tied (zero-redundancy rows) -> std, said
    % once per run; a near-tie (1e-15 spread) counts as a tie; healthy data
    % stays on the MAD.
    muteWarn(tc, 'lsqvce:robustWeights:scaleFallback');   % provoked on purpose
    lsqvce.util.log.inventory('reset');
    vw = [zeros(51, 1); linspace(-1, 1, 44).'; 20; -20; 20; -20; 20];   % 51 % ties
    [~, s, ~, p] = lsqvce.kernel.robustWeights(vw, 'huber', []);
    tc.verifyEqual(p, 'std');
    % tol: 1e-14 relative -- the same std on the same 100 numbers
    tc.verifyEqual(s, std(vw), 'RelTol', 1e-14);
    tc.verifyTrue(any(strcmp(lsqvce.util.log.inventory(), 'lsqvce:robustWeights:scaleFallback')), ...
                  'the fallback is reported, not silent');
    near = vw;  near(1:51) = 1e-15 * (1:51).';                       % near-ties
    [facN, ~, ~, pN] = lsqvce.kernel.robustWeights(near, 'huber', []);
    tc.verifyEqual(pN, 'std', 'a 1e-15 MAD is a tie, not a scale');
    tc.verifyGreaterThan(median(facN(52:95)), 0.5, 'ordinary rows keep their weight');
    [~, ~, ~, pH] = lsqvce.kernel.robustWeights(linspace(-3, 3, 101).', 'huber', []);
    tc.verifyEqual(pH, 'mad');
    lsqvce.util.log.inventory('reset');
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
    % R7-06: %g kept six significant digits, so YYYYMMDD days and Julian Day
    % Numbers collided (30 September days -> ONE key) and robust factors of
    % one day were applied to another. Every day of a month must stay apart.
    yyyymmdd = 20260901:20260930;
    k = arrayfun(@(d) lsqvce.util.rowWeightKey(1, d, 'kbr', 1), yyyymmdd, 'UniformOutput', false);
    tc.verifyEqual(numel(unique(k)), 30, 'YYYYMMDD days must give distinct keys');
    tc.verifyEqual(lsqvce.util.rowWeightKey(1, 20260901, 'kbr', 1), '1|20260901|kbr|1');
    tc.verifyEqual(lsqvce.util.rowWeightKey(2460580, 2460581, 'p', 7), '2460580|2460581|p|7');
    tc.verifyEqual(lsqvce.util.rowWeightKey(1, 0.1, 'p', 1), '1|0.10000000000000001|p|1');   % non-integer: round-trip exact
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
    % a full SPD weight matrix (correlated rows) is whitened exactly by its
    % Cholesky factor (R5 M5 5b'); a non-SPD one is refused by name
    Lm = lsqvce.engine.batchSystem(3, @(k) deal(A, l, 2*eye(12)), 1);
    arcM = Lm.closure(1, 'obs', 1, uo);
    % tol: 1e-14 relative -- chol(2 I) = sqrt(2) I applied to O(1) numbers, a few ulps
    tc.verifyEqual(arcM.A, sqrt(2) * A, 'RelTol', 1e-14);
    Lb = lsqvce.engine.batchSystem(3, @(k) deal(A, l, -eye(12)), 1);
    verifyError(tc, @() Lb.closure(1, 'obs', 1, uo), 'lsqvce:batchSystem:weightMatrix');
    % a wrong A shape is named
    Lw = lsqvce.engine.batchSystem(4, @(k) deal(A, l), 1);
    verifyError(tc, @() Lw.closure(1, 'obs', 1, uo), 'lsqvce:batchSystem:badMatrix');
    % WP-133 (review 3 finding 17): m x 1 weights are >= 0 and finite (a
    % negative one made the rows complex: the QR and normal backends then
    % disagreed); a zero weight zeroes its row and marks it excluded -- kept
    % for its index, out of the dof
    pz = p;  pz(3) = 0;
    Lz = lsqvce.engine.batchSystem(3, @(k) deal(A, l, pz), 1);
    arcZ = Lz.closure(1, 'obs', 1, uo);
    ex = false(12, 1);  ex(3) = true;
    tc.verifyEqual(arcZ.excluded, ex, 'a zero weight marks its row excluded');
    tc.verifyEqual(arcZ.A(3, :), zeros(1, 3), 'and zeroes it');
    for bad = [-1, NaN, Inf]
        pb = p;  pb(5) = bad;
        Lbad = lsqvce.engine.batchSystem(3, @(k) deal(A, l, pb), 1);
        verifyError(tc, @() Lbad.closure(1, 'obs', 1, uo), 'lsqvce:batchSystem:weightMatrix');
    end
    % the day contract's .excluded marks rows that ARE zero; a live row is refused
    verifyError(tc, @() excludedDay(A, l, ex, uo), 'lsqvce:makeLoader:badDay');
    arcX = excludedDay(A .* ~ex, l .* ~ex, ex, uo);
    tc.verifyEqual(arcX.excluded, ex, 'a zeroed row marked excluded travels with the arc');
end

function arc = excludedDay(A, l, ex, uo)
%EXCLUDEDDAY  One arc of a day that declares rows excluded (WP-133).
    D = struct('A', {{A}}, 'B', {{zeros(size(A, 1), 0)}}, 'l', {{l}}, ...
               't', {{(1:numel(l)).'}}, 'excluded', {{ex}});
    cov.obs = struct('type', 'diag', 'sigma', 1);
    L = lsqvce.engine.makeLoader(@(day, ot) D, 'ObservationTypes', {'obs'}, 'Covariance', cov, ...
                                 'Arcs', 1, 'SystemId', 1, 'NumLocal', 0, 'NumCommon', size(A, 2), 'Days', 1);
    arc = L.closure(1, 'obs', 1, uo);
end

%% ============================================================ R7-27: contract points the audit had credited from comments
% Audit check 16 used to credit a unit from any `lsqvce.<name>` mention,
% comments included, and accepted 'contract' exemptions for units that
% validate their arguments. Counted from code, these eight had no contract
% point; each now rejects by its documented identifier.
function testAdjustNonlinearContract(tc)
% categories: contract
    bf = @(k, x) deal(ones(3, 2), zeros(3, 1));
    verifyError(tc, @() lsqvce.engine.adjustNonlinear(2, bf, 1, [1; 2; 3]), ...
                'lsqvce:adjustNonlinear:badX0');
    verifyError(tc, @() lsqvce.engine.adjustNonlinear(0, bf, 1, 1), ...
                'MATLAB:validators:mustBePositive');
    verifyError(tc, @() lsqvce.engine.adjustNonlinear(2, bf, 1, [0; 0], 'Backend', "lu"), ...
                'MATLAB:validators:mustBeMember');
end

function testDiagnoseContract(tc)
% categories: contract
    % no batch count on the solution and none given: refused by name
    verifyError(tc, @() lsqvce.engine.verdict.diagnose(struct(), @(k) deal(1, 1, 1)), ...
                'lsqvce:diagnose:nBatch');
    verifyError(tc, @() lsqvce.engine.verdict.diagnose(struct(), [], 'Test', "bogus"), ...
                'MATLAB:validators:mustBeMember');
    verifyError(tc, @() lsqvce.engine.verdict.diagnose(struct(), [], 'Alpha', 2), ...
                'MATLAB:validators:mustBeLessThan');
end

function testResolveNcContract(tc)
% categories: contract
    verifyError(tc, @() lsqvce.engine.resolveNc({}), 'lsqvce:resolveNc:empty');
    verifyError(tc, @() lsqvce.engine.resolveNc({struct('nl', 1)}), 'lsqvce:resolveNc:missing');
    verifyError(tc, @() lsqvce.engine.resolveNc({struct('isLoader', true)}), 'lsqvce:resolveNc:missing');
    verifyError(tc, @() lsqvce.engine.resolveNc({struct('nc', 2), struct('nc', 3)}), ...
                'lsqvce:resolveNc:conflict');
    [nc, info] = lsqvce.engine.resolveNc({struct('nc', 4), struct('nc', 4)});
    tc.verifyEqual(nc, 4);
    tc.verifyEqual(info.values, [4 4]);
end

function testRunMetaContract(tc)
% categories: contract
    verifyError(tc, @() lsqvce.engine.runMeta(1, struct()), 'lsqvce:runMeta:badArg');
    verifyError(tc, @() lsqvce.engine.runMeta(struct(), 'x'), 'lsqvce:runMeta:badArg');
    verifyError(tc, @() lsqvce.engine.runMeta(struct()), 'MATLAB:narginchk:notEnoughInputs');
    m = lsqvce.engine.runMeta(struct('probe', 7), struct('a', 1));
    tc.verifyEqual(m.seeds.probe, 7);
    tc.verifyEqual(m.toolboxVersion, lsqvce.util.toolboxVersion());
end

function testReducedContributionContract(tc)
% categories: contract
    verifyError(tc, @() lsqvce.engine.ReducedContribution(ones(2, 3), [1; 1], 0, 1), ...
                'lsqvce:reduced:notSquare');
    verifyError(tc, @() lsqvce.engine.ReducedContribution(eye(2), [1; 1; 1], 0, 1), ...
                'lsqvce:reduced:sizeMismatch');
    c = lsqvce.engine.ReducedContribution(eye(2), [1 2], 3, 4);
    tc.verifyEqual(c.n, [1; 2], 'n is stored as a column');
end

function testFtestContract(tc)
% categories: contract
    s = @(u, n, v) struct('rank', u, 'u', u, 'n', n, 'vPv', v);
    verifyError(tc, @() lsqvce.engine.verdict.ftest(s(3, 10, 2), s(2, 11, 3), 'Verbose', false), ...
                'lsqvce:ftest:differentData');
    verifyError(tc, @() lsqvce.engine.verdict.ftest(s(3, 10, 2), s(3, 10, 3), 'Verbose', false), ...
                'lsqvce:ftest:notNested');
    verifyError(tc, @() lsqvce.engine.verdict.ftest(s(10, 10, 2), s(2, 10, 3), 'Verbose', false), ...
                'lsqvce:ftest:noRedundancy');
    verifyError(tc, @() lsqvce.engine.verdict.ftest(s(3, 10, 0), s(2, 10, 0), 'Verbose', false), ...
                'lsqvce:ftest:zeroVariance');
    verifyError(tc, @() lsqvce.engine.verdict.ftest(s(3, 10, 2), s(2, 10, 3), 'Alpha', 2), ...
                'MATLAB:validators:mustBeLessThan');
    % WP-135 (review 3, lane 1): the full model's redundancy where it carries
    % one -- rows of weight 0 and a rank-deficient local block make it differ
    % from n - rank (compare() passes R.diag.r)
    sr = s(3, 10, 2);  sr.r = 6;
    F = lsqvce.engine.verdict.ftest(sr, s(2, 10, 3), 'Verbose', false);
    verifyEqual(tc, F.df2, 6);
    verifyEqual(tc, F.df1, 1);
end

function testAdjustRobustContract(tc)
% categories: contract
    % the argument block rejects before any system is read
    verifyError(tc, @() lsqvce.engine.adjustRobust({}, 1, 'Influence', "cauchy"), ...
                'MATLAB:validators:mustBeMember');
    verifyError(tc, @() lsqvce.engine.adjustRobust({}, 1, 'MaxIter', 0), ...
                'MATLAB:validators:mustBePositive');
    verifyError(tc, @() lsqvce.engine.adjustRobust({}, 1, 'Backend', "lu"), ...
                'MATLAB:validators:mustBeMember');
end

function testDeclarationContract(tc)
% categories: contract
    verifyError(tc, @() lsqvce.engine.declaration(struct('bogus', 1)), 'lsqvce:fit:unknownField');
    verifyError(tc, @() lsqvce.engine.declaration(struct('commonParams', {{'common', 1}})), ...
                'lsqvce:fit:noSystems');
    sys = struct('buildDay', @(d) 1, 'covariance', struct('obs', 1));
    verifyError(tc, @() lsqvce.engine.declaration(struct('systems', sys)), 'lsqvce:fit:noCommon');
    verifyError(tc, @() lsqvce.engine.declaration(struct('systems', sys, ...
                'commonParams', {{'arc', 1}})), 'lsqvce:fit:badScope');
end

function testResultsRefusalsByName(tc)
% categories: contract
    % R7-28: two USER-tier identifiers no test had fired. A Results built
    % without an accumulator (kinds B/C, or by hand) refuses the accumulator
    % methods by name, and compare() refuses a solution without a v'Pv.
    sol = struct('x', [1; 2], 'sigma0', 1);
    Rb = lsqvce.Results(sol, struct(), [], struct(), @() deal([], [], []), 0, []);
    verifyError(tc, @() Rb.residuals(), 'lsqvce:Results:noAccumulator');
    verifyError(tc, @() Rb.contribution(), 'lsqvce:Results:noAccumulator');
    verifyError(tc, @() Rb.compare(Rb), 'lsqvce:Results:noVPv');
end

function testCompareUsesTheSolvesRankAndRedundancy(tc)
% categories: contract
    % WP-137 (code review of WP-135, its lane-1 low "compare ignores D.rank"):
    % Results.compare handed the engine rank = u and r = n - u, so a fit whose
    % local block is not fully estimable (D.rank < u, the pinv fallback) got
    % df2 one short, s0^2 too large and F too small. testFtestContract pins
    % the engine; here compare passes D's numbers. Hand-built Results:
    % u = 6 with one local direction not estimable (rank 5, r = 45 of n = 50).
    DF = struct('u', 6, 'n', 50, 'rank', 5, 'r', 45);
    DR = struct('u', 5, 'n', 50, 'rank', 4, 'r', 46);
    RF = lsqvce.Results(struct('x', zeros(2, 1), 'vtv', 10), DF, [], struct(), @() deal([], [], []), 0, []);
    RR = lsqvce.Results(struct('x', zeros(1, 1), 'vtv', 14), DR, [], struct(), @() deal([], [], []), 0, []);
    F = RF.compare(RR);
    verifyEqual(tc, [F.df1, F.df2], [1, 45], 'df from the solve''s rank and redundancy');
    % tol: the closed form ((14 - 10) / 1) / (10 / 45) = 18, rounding only
    verifyEqual(tc, F.F, 18, 'RelTol', 1e-14);
end

function testRobustWeightsSurviveAnExtremeBlunder(tc)
% categories: robustness, precision
    % WP-134 (review 3 finding 9): the near-tie test compared the MAD with the
    % LARGEST deviation -- the worst blunder. One residual beyond 1/(1e3 eps)
    % ~ 4.5e12 MAD switched the scale to the std it inflates: s 1.6e11, the
    % Huber estimate moved by 2e8, 8-sigma rows kept factor 1
    % (dev/validation/wp134_kernel_mirror.py). The reference is the 75th
    % percentile of the deviations now, which no blunder in the upper quarter
    % moves; the R7-21 tie fixtures above still fall back.
    rs = RandStream('mt19937ar', 'Seed', 134);
    vw = [randn(rs, 999, 1); 1e13];
    lsqvce.util.log.inventory('reset');
    [fac, s, ~, p] = lsqvce.kernel.robustWeights(vw, 'huber', []);
    tc.verifyEqual(p, 'mad');
    % tol: 10 % -- the normalised MAD of 999 N(0,1) draws estimates 1 to ~3 %
    tc.verifyEqual(s, 1, 'RelTol', 0.1);
    % tol: Huber's c / |w| at 1e13 MAD is ~1.3e-13, well below 1e-11
    tc.verifyLessThan(fac(end), 1e-11);
    tc.verifyFalse(any(strcmp(lsqvce.util.log.inventory(), 'lsqvce:robustWeights:scaleFallback')));
    % tuning values are validated: a c <= 0 zeroed or negated every factor,
    % an igg3 with k0 >= k1 became a hard 0/1 rejection
    tc.verifyError(@() lsqvce.kernel.robustWeights(vw, 'huber', 0), 'lsqvce:robustWeights:badTuning');
    tc.verifyError(@() lsqvce.kernel.robustWeights(vw, 'bisquare', NaN), 'lsqvce:robustWeights:badTuning');
    tc.verifyError(@() lsqvce.kernel.robustWeights(vw, 'igg3', [3 1.5]), 'lsqvce:robustWeights:badTuning');
    lsqvce.util.log.inventory('reset');
end
