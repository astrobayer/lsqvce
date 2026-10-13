function tests = test_speedBudget
%TEST_SPEEDBUDGET  Measured speed budgets as one-line ratio assertions (R5 M7).
%   tests = test_speedBudget
%
%   Every budget here is a RATIO measured by dev/tests/private/timePair
%   (warm-up, alternation, per-repeat ratios, the median) against a bound
%   pinned from a measurement in the ledger -- never an absolute time, so it
%   holds on a laptop and on the hosted runner alike. The bounds are loose
%   (2x the measured value or worse) on purpose: a budget catches an
%   algorithmic regression (an O(n^2) creeping back), not a slow day.
%
%   Budgets and their measurements:
%     lag sums by FFT vs the per-lag loop     M3 3c: 2160 lags, 1 ms/lag/arc
%                                             -> flat in nLag; bound: the FFT
%                                             route is at most 0.5x the loop at
%                                             nLag = 512
%     Hutchinson trace vs dense inverse       measured exact/probe 3.02 at
%     at nc = 500                             nc = 500; bound: 1 <= ratio < 10
%     one accumulate + solve on the demo      an absolute sanity bound only for
%                                             the runner: < 20 s (it is ~1 s)
%
%   See also TIMEPAIR, TEST_INSTRUMENTS, dev/DEVELOPMENT.md (WP-79 table)
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-134 (third review): the armaACV point's cost model (the response follows its decay).
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M7, WP-89.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    pathBefore = path;   % R7-45: a standalone run (run_lsqvce / runtests) leaves the path as found
    tc.addTeardown(@() path(pathBefore));
    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    addpath(root);  addpath(fullfile(root, 'examples'));
    tc.TestData.root = root;
end

function testFftLagSumsBeatTheLoop(tc)
% categories: speed
    % lsqvce.kernel.redundancyLags' sub-diagonal sums are one FFT per column;
    % the loop form (kept here as the reference) is O(N nLag). Bound 0.5:
    % measured 0.05-0.15 at N = 2048, nLag = 512 (2026-09-03, PCWIN64 R2026a).
    N = 2048;  C = 40;  nL = 512;
    rs = RandStream('mt19937ar', 'Seed', 1);  X = randn(rs, N, C);
    fft = @() lagSumsFft(X, nL);
    loop = @() lagSumsLoop(X, nL);
    % tol: the two routes agree to rounding (a precision leg inside the budget)
    tc.verifyEqual(fft(), loop(), 'RelTol', 1e-10, 'AbsTol', 1e-9);
    % R7-42: Repeats back to timePair's default 15 (its help: 5 repeats ->
    % 31 % spread, 15 -> 10 %). NOTE: both sides are TEST-LOCAL copies, so this
    % budget pins the algorithm, not the kernel's columnLagSums -- moving it
    % onto lsqvce.kernel.redundancyLags needs a measured bound (tier A).
    rec = timePair(fft, loop, 'Inner', 1, 'Label', 'fft/loop lag sums');
    tc.verifyLessThan(rec.ratio, 0.5, sprintf('FFT lag sums must be < 0.5x the loop (measured %.3f)', rec.ratio));
end

function testProbeAndExactSweepsAreComparable(tc)
% categories: speed
    % the Hutchinson trace with 100 probes (lsqvce.engine.traceEstimate) vs the
    % dense inverse at nc = 500. MEASURED exact/probe = 3.02 (2026-09-03,
    % PCWIN64 R2026a: two triangular solves per probe beat inv(N)); the bound
    % pins the ORDER -- the probe route must not fall behind the dense inverse
    % (ratio >= 1) and must not be pathologically slow (ratio <= 10). (WP-79's
    % 1.04-1.19 was the FULL sweep incl. re-reads -- a different quantity; its
    % conclusion, the exact sweep as default, stands.)
    nc = 500;  rs = RandStream('mt19937ar', 'Seed', 2);
    M = randn(rs, nc, nc);  N = M.'*M + nc*eye(nc);  R = chol(N);
    exact = @() sum(diag(inv(N)));
    probe = @() lsqvce.engine.traceEstimate(@(Z) R \ (R.' \ Z), 'Size', nc, 'NumProbes', 100, 'Seed', 3);
    rec = timePair(exact, probe, 'Inner', 1, 'Label', 'exact/probe trace');   % R7-42: default 15 repeats
    tc.verifyGreaterThanOrEqual(rec.ratio, 1, sprintf('the probe route fell behind the dense inverse (ratio %.2f)', rec.ratio));
    tc.verifyLessThan(rec.ratio, 10, sprintf('the dense inverse cannot be 10x slower than 100 probes (ratio %.2f)', rec.ratio));
end

function testOneSolveOnTheDemoIsBounded(tc)
% categories: speed
    % an absolute sanity bound for the runner only: one accumulate + solve
    % (lsqvce.engine.adjustBlock through lsqvce.fit) on the GRACE demo
    % (3 days, 2 systems, 32 arcs) takes ~1 s on the user's machine; 20 s
    % would mean something quadratic came back
    [~, data] = demo_grace_dataLoader();
    t0 = tic;  R = lsqvce.fit(data, Verbose=false);  dt = toc(t0);
    tc.verifyTrue(all(isfinite(R.x)));
    tc.verifyLessThan(dt, 20, sprintf('one solve took %.1f s', dt));
end

%% ============================================================ R7-42: the speed gaps, gated
% Until R7-27 the matrix credited these units' speed from test_speed, which
% ci_gate never runs. Each relation below is the benchmark's own claim, now
% measured by timePair (both sides in one repeat, rotation, median of
% per-repeat ratios) and asserted in the gate. A bound copied from test_speed
% keeps that suite's provenance; a bound DERIVED from the algorithm says so.
% Repeats = 9 on the heavy pairs (timePair's help: 14 % spread) -- every bound
% below sits at least 2x away from its expected value.

function testFilterWhiteningBeatsDense(tc)
% categories: speed
    % lsqvce.kernel.armaWhiten is O(m p); the dense route (decorrFactor of the
    % Toeplitz covariance) is O(m^3). test_speed's claim: the filter wins
    % outright at realistic batch sizes. tol: ratio < 1 (test_speed bound; at
    % m = 1000 the cubic route is an order of magnitude behind)
    a = [1 -1.6 0.9 -0.25];  s = 0.7;  m = 1000;
    rs = RandStream('mt19937ar', 'Seed', 4);
    r = s^2 * lsqvce.kernel.armaACV(a, 1, m);
    X = randn(rs, m, 20);  z = zeros(m, 0);
    filt  = @() lsqvce.kernel.armaWhiten(a, 1, s, r, z, z, X, true);
    dense = @() lsqvce.kernel.decorrFactor(toeplitz(r)) * X;
    rec = timePair(filt, dense, 'Label', 'armaWhiten/dense whitening');
    tc.verifyLessThan(rec.ratio, 1, sprintf('filter/dense = %.3f [%.3f %.3f]', rec.ratio, rec.band));
end

function testArmaACVIsNearLinearInLags(tc)
% categories: speed
    % lsqvce.kernel.armaACV: an impulse response as long as its decay (WP-134:
    % it was max(4n, 512) taps, an 8n-point FFT -- 65 536 points at n = 8000,
    % out of cache on some machines: 12-17x there, HANDOVER 12.68) and one
    % FFT of nh + n points -- O(n log n). tol: 4x the lags may cost at most 8x
    % (DERIVED: n log n gives ~4.6x; a quadratic regression gives 16x)
    a = [1 -1.6 0.9 -0.25];
    rec = timePair(@() lsqvce.kernel.armaACV(a, 1, 8000), @() lsqvce.kernel.armaACV(a, 1, 2000), ...
                   'Inner', 5, 'Label', 'armaACV 4n/n');
    tc.verifyLessThan(rec.ratio, 8, sprintf('armaACV(4n)/armaACV(n) = %.2f [%.2f %.2f]', rec.ratio, rec.band));
end

function testLevinsonIsQuadraticInOrder(tc)
% categories: speed
    % lsqvce.kernel.levinsonDurbin is O(p^2). tol: 4x the order may cost at
    % most 4^2.6 = 36.8x (test_speed's exponent bound 2.6; O(p^3) gives 64x)
    r = lsqvce.kernel.armaACV([1 -1.6 0.9 -0.25], 1, 401);
    rec = timePair(@() lsqvce.kernel.levinsonDurbin(r, 400), @() lsqvce.kernel.levinsonDurbin(r, 100), ...
                   'Inner', 3, 'Label', 'levinsonDurbin 4p/p');
    tc.verifyLessThan(rec.ratio, 4^2.6, sprintf('levinson(4p)/levinson(p) = %.1f [%.1f %.1f]', rec.ratio, rec.band));
end

function testPsd2acvIsLinearInLags(tc)
% categories: speed
    % lsqvce.kernel.psd2acv: one quadrature per lag -- linear. tol: 4x the lags
    % may cost at most 8x (test_speed bound)
    f = linspace(0, 0.1, 4096).';  S = exp(-f/0.01);
    rec = timePair(@() lsqvce.kernel.psd2acv(f, S, 5, 400), @() lsqvce.kernel.psd2acv(f, S, 5, 100), ...
                   'Inner', 3, 'Label', 'psd2acv 4L/L');
    tc.verifyLessThan(rec.ratio, 8, sprintf('psd2acv(4L)/psd2acv(L) = %.2f [%.2f %.2f]', rec.ratio, rec.band));
end

function testQRSystemCostsAFewNormalEquations(tc)
% categories: speed
    % lsqvce.engine.QRSystem buys a squared condition number for ~2x the flops
    % of lsqvce.engine.NormalEqSystem; it must keep R truncated to dim+1 rows.
    % tol: ratio < 4 (test_speed's backend bound)
    dim = 12;  rs = RandStream('mt19937ar', 'Seed', 5);
    A = randn(rs, 500, dim);  l = randn(rs, 500, 1);
    rec = timePair(@() feedSystem(lsqvce.engine.QRSystem(dim), A, l, 40), ...
                   @() feedSystem(lsqvce.engine.NormalEqSystem(dim), A, l, 40), ...
                   'Repeats', 9, 'Label', 'QRSystem/NormalEqSystem');
    tc.verifyLessThan(rec.ratio, 4, sprintf('QR/normal = %.2f [%.2f %.2f]', rec.ratio, rec.band));
end

function testBlockAccumulatorIsLinearInDays(tc)
% categories: speed
    % lsqvce.engine.BlockAccumulator holds (nc+nl+1)^2 plus one arc; 4x the
    % days must cost ~4x. tol: < 4^1.4 = 7.0 (test_speed's exponent bound 1.4)
    cfg = speedCfg();
    rec = timePair(@() runAcc(cfg, 1:8), @() runAcc(cfg, 1:2), 'Repeats', 9, ...
                   'Label', 'BlockAccumulator 8 days / 2 days');
    tc.verifyLessThan(rec.ratio, 4^1.4, sprintf('8 days / 2 days = %.2f [%.2f %.2f]', rec.ratio, rec.band));
end

function testFitAddsLittleToTheBlockDriver(tc)
% categories: speed
    % lsqvce.fit = declaration (lsqvce.engine.makeLoader per system) + Options
    % mapping + lsqvce.engine.adjustBlock. Both sides build the loaders, so the
    % ratio is fit's own dispatch. tol: < 1.25 (timePair's example: expected
    % 1.00, measured 2026-08-31 L4 postflight, budget 1.25)
    [~, data] = demo_grace_dataLoader();
    opts = lsqvce.Options(Verbose=false);
    rec = timePair(@() lsqvce.fit(data, opts), @() directBlock(data), 'Repeats', 9, ...
                   'Label', 'fit/adjustBlock dispatch');
    tc.verifyLessThan(rec.ratio, 1.25, sprintf('fit/adjustBlock = %.3f [%.3f %.3f]', rec.ratio, rec.band));
end

function testKernelLagSumsAreFlatInLags(tc)
% categories: speed
    % The kernel budget the review asked for (R7-42): lsqvce.kernel.redundancyLags
    % computes every lag sum of its N x Z blocks by one FFT per column, so 32x
    % the lags adds only the O(N nLag) Toeplitz sub-diagonal loop. tol: < 3
    % (DERIVED, N = 1024, Z = 200: the fixed solve dominates, ~1.2x expected;
    % a per-lag loop regression in columnLagSums costs ~10x)
    N = 1024;  Z = 200;  P = 4;  rs = RandStream('mt19937ar', 'Seed', 6);
    W = eye(N);                              % Toeplitz: the structured path
    Abar = randn(rs, N, P);  UiZ = randn(rs, P, Z) / sqrt(Z);  e = randn(rs, N, 1);
    many = @() lsqvce.kernel.redundancyLags(Abar, W, UiZ, e, 512);
    few  = @() lsqvce.kernel.redundancyLags(Abar, W, UiZ, e, 16);
    rec = timePair(many, few, 'Repeats', 9, 'Label', 'redundancyLags 512/16 lags');
    tc.verifyLessThan(rec.ratio, 3, sprintf('512 lags / 16 lags = %.2f [%.2f %.2f]', rec.ratio, rec.band));
end

function s = feedSystem(s, A, l, k)
    for j = 1:k
        s.add(A, l);
    end
end

function R = directBlock(data)
    [spec, systems] = lsqvce.engine.declaration(data);
    R = lsqvce.engine.adjustBlock(systems, spec.days, 'Verbose', false);
end

function cfg = speedCfg()
% test_speed's synthLoader fixture (3 systems, 2 obs types, 2 arcs/day).
    rs = RandStream('mt19937ar', 'Seed', 42);
    cfg.nc = 6;  cfg.nl = 4;  cfg.sats = [1 2 3];
    cfg.x_c = randn(rs, cfg.nc, 1);
    cfg.y = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for s = cfg.sats, cfg.y(s) = randn(rs, cfg.nl, 1); end
    cfg.obsTypes = {'pos', 'kbr'};  cfg.arcsPerDay = 2;  cfg.mArc = 400;  cfg.dt = 5;
    cfg.seed = 42;  cfg.collinear = false;  cfg.collinearEps = 1e-7;
    cfg.cov.pos = struct('type', 'arma', 'ar', [1 -1.6 0.9 -0.25], 'ma', 1, 'sigma', 0.7);
    cfg.cov.kbr = struct('type', 'diag', 'sigma', 0.2);
end

function sol = runAcc(cfg, days)
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend', 'qr', 'Verbose', false);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, @(d, ot, a, uo) lsqvce.engine.synthLoader(cfg, s, d, ot, a, uo), days);
    end
    sol = acc.solve();
end

function d = lagSumsFft(X, nL)
    [N, ~] = size(X);  nfft = 2^nextpow2(2*N - 1);
    F = fft(X, nfft, 1);  r = ifft(sum(real(F .* conj(F)), 2), 'symmetric');
    d = r(1:nL);
end

function d = lagSumsLoop(X, nL)
    N = size(X, 1);  d = zeros(nL, 1);
    d(1) = sum(X(:).^2);
    for n = 1:nL-1
        idx = 1:(N-n);
        d(n+1) = sum(sum(X(idx+n, :) .* X(idx, :), 2));
    end
end
