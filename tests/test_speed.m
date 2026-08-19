function tests = test_speed
%TEST_SPEED  Performance, complexity scaling, and memory behaviour.
%
%   These tests assert SCALING and RELATIVE cost, not absolute wall-clock
%   times -- absolute timings are machine dependent and would make the suite
%   flaky. What must hold on any machine:
%
%     * the filter path beats the dense path, and by more as the batch grows
%     * the filter path scales O(m), the dense path O(b^3)
%     * QR costs no more than ~3x the normal equations
%     * memory stays bounded regardless of the number of days
%
%   Run:  results = runtests('tests/test_speed');  disp(table(results))
%
%   PROFILING TIP
%     If a real run is slow, the cause is almost always a dense covariance with
%     too large a BatchSize. Check acc-side with:
%         [~,~,info] = lsqvce.covModel(cfg.cov.pos, t, opt);
%         info.cost        % 'O(b*40)' good, 'O(b^3)' expensive
%
%   See also TEST_PRECISION, TEST_RELIABILITY

    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    tc.TestData.cfg = makeCfg();
    % Warm up the JIT so the first timed call is not penalised.
    runAcc(tc.TestData.cfg, 1, 'Backend','qr', 'Verbose',false);
end

%% ============================================ filter vs dense
function testFilterBeatsDenseWhitening(tc)
% categories: speed
% The central performance claim. The FIR filter is O(m*p); the dense Cholesky
% is O(b^3). At realistic batch sizes the filter must win outright.
    a = [1 -1.6 0.9 -0.25];   s = 0.7;   p = 3;
    nCol = 20;

    for n = [500 1000 2000]
        r = s^2 * lsqvce.armaACV(a, 1, n);
        X = randn(n, nCol);
        z = zeros(n, 0);

        tFilt = timeit(@() lsqvce.armaWhiten(a, 1, s, r, z, z, X, true));
        tDense = timeit(@() denseWhiten(toeplitz(r), X));

        tc.verifyLessThan(tFilt, tDense, ...
            sprintf('n=%d: the filter (%.2f ms) must beat dense (%.2f ms)', ...
                    n, tFilt*1e3, tDense*1e3));
    end
end

function testFilterAdvantageGrowsWithBatchSize(tc)
% categories: speed
% O(m*p) vs O(b^3): the speed ratio must GROW with the batch size. If it does
% not, the dense path is being used by accident somewhere.
    a = [1 -1.6 0.9 -0.25];   s = 0.7;
    ratio = zeros(1,2);
    nn    = [400 1600];

    for k = 1:2
        n = nn(k);
        r = s^2 * lsqvce.armaACV(a, 1, n);
        X = randn(n, 10);
        z = zeros(n, 0);

        tF = timeit(@() lsqvce.armaWhiten(a, 1, s, r, z, z, X, true));
        tD = timeit(@() denseWhiten(toeplitz(r), X));
        ratio(k) = tD / tF;
    end

    % 4x the batch size should widen the gap by well over 2x
    tc.verifyGreaterThan(ratio(2)/ratio(1), 2, ...
        sprintf('speed ratio went %.0fx -> %.0fx; expected superlinear growth', ...
                ratio(1), ratio(2)));
end

function testFilterScalesLinearly(tc)
% categories: speed
% Doubling the arc length must roughly double the filter cost, not quadruple
% it. Tolerance is generous because timing on a shared machine is noisy.
    a = [1 -1.6 0.9 -0.25];   s = 0.7;
    t = zeros(1,3);
    nn = [2000 4000 8000];

    for k = 1:3
        n = nn(k);
        r = s^2 * lsqvce.armaACV(a, 1, max(numel(a), 2));
        X = randn(n, 10);
        z = zeros(n, 0);
        t(k) = timeit(@() lsqvce.armaWhiten(a, 1, s, r, z, z, X, true));
    end

    % empirical exponent from the first and last point
    expo = log(t(3)/t(1)) / log(nn(3)/nn(1));
    tc.verifyLessThan(expo, 1.6, ...
        sprintf('filter scaling exponent %.2f -- should be ~1.0', expo));
end

%% ============================================ backend cost
function testQROverheadIsAcceptable(tc)
% categories: speed
% QR buys a squared condition number for roughly 2x the flops. If it costs
% much more than 3x, something is wrong (e.g. R is not being truncated).
    cfg = tc.TestData.cfg;

    tN = timeit(@() runAcc(cfg, 1:2, 'Backend','normal', 'Verbose',false));
    tQ = timeit(@() runAcc(cfg, 1:2, 'Backend','qr',     'Verbose',false));

    tc.verifyLessThan(tQ/tN, 4.0, ...
        sprintf('QR is %.1fx the cost of normal equations; expected < 4x', tQ/tN));
end

function testTriangularFactorStaysBounded(tc)
% categories: contract, speed
% The QR factor must be truncated to dim+1 rows after every batch. If it were
% not, memory would grow linearly with the number of observations.
    dim = 12;
    acc = lsqvce.QRSystem(dim);
    for k = 1:60
        acc.add(randn(500, dim), randn(500,1));
    end
    tc.verifyLessThanOrEqual(size(acc.R,1), dim+1);
    tc.verifyEqual(acc.m, 60*500);
end

%% ============================================ accumulator scaling
function testCostScalesLinearlyWithDays(tc)
% categories: speed
% Doubling the number of days must roughly double the runtime. Anything
% superlinear means something is accumulating that should not be.
    cfg = tc.TestData.cfg;

    t1 = timeit(@() runAcc(cfg, 1:2, 'Backend','qr', 'Verbose',false));
    t2 = timeit(@() runAcc(cfg, 1:4, 'Backend','qr', 'Verbose',false));
    t3 = timeit(@() runAcc(cfg, 1:8, 'Backend','qr', 'Verbose',false));

    expo = log(t3/t1) / log(4);
    tc.verifyLessThan(expo, 1.4, ...
        sprintf('runtime scaling exponent %.2f in days -- should be ~1.0', expo));
end

function testMemoryIsIndependentOfDays(tc)
% categories: contract, speed
% The whole point of the lazy loader: peak memory must NOT grow with the
% campaign length. The accumulator holds only (nc+nl+1)^2 plus one arc.
    cfg    = tc.TestData.cfg;
    nDays  = [1 8];
    nbytes = zeros(1, numel(nDays));

    for k = 1:numel(nDays)
        acc = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
        acc.addSystem(1, cfg.nl, binder(cfg,1), 1:nDays(k));

        w = whos('acc');
        nbytes(k) = w.bytes;
    end

    % The accumulator holds (nc+nl+1)^2 plus the Schur blocks -- neither grows
    % with the campaign length. Only dayList itself gets longer.
    tc.verifyLessThan(nbytes(2)/nbytes(1), 1.5, ...
        sprintf('accumulator grew %.2fx from 1 day to 8 days', nbytes(2)/nbytes(1)));
end

function testBatchSizeDoesNotChangeCostOrder(tc)
% categories: speed
% With a filter covariance, BatchSize should barely affect the runtime -- the
% cost is O(m*p) either way. A strong dependence means the dense path is being
% taken by mistake.
    cfg = tc.TestData.cfg;

    tSmall = timeit(@() runAcc(cfg, 1:2, 'Backend','qr', ...
                               'BatchSize',200, 'Verbose',false));
    tLarge = timeit(@() runAcc(cfg, 1:2, 'Backend','qr', ...
                               'BatchSize',1e9, 'Verbose',false));

    % small batches add loop overhead, but nothing like an order of magnitude
    tc.verifyLessThan(tSmall/tLarge, 6, ...
        sprintf('BatchSize 200 vs whole-arc: %.1fx -- filter path should be flat', ...
                tSmall/tLarge));
end

%% ============================================ primitive cost
function testLevinsonIsQuadratic(tc)
% categories: speed
% Levinson-Durbin is O(p^2). Doubling the order must roughly quadruple the
% cost -- and stay far below the O(p^3) of a direct Yule-Walker solve.
    r = lsqvce.armaACV([1 -1.6 0.9 -0.25], 1, 401);

    t1 = timeit(@() lsqvce.levinsonDurbin(r, 100));
    t2 = timeit(@() lsqvce.levinsonDurbin(r, 200));
    t3 = timeit(@() lsqvce.levinsonDurbin(r, 400));

    expo = log(t3/t1) / log(4);
    tc.verifyLessThan(expo, 2.6, ...
        sprintf('Levinson scaling exponent %.2f -- should be ~2.0', expo));
end

function testPSD2ACVIsLinearInLags(tc)
% categories: speed
    f = linspace(0, 0.1, 4096).';
    S = exp(-f/0.01);

    t1 = timeit(@() lsqvce.psd2acv(f, S, 5, 100));
    t2 = timeit(@() lsqvce.psd2acv(f, S, 5, 400));

    tc.verifyLessThan(t2/t1, 8, ...
        sprintf('psd2acv is %.1fx slower for 4x the lags; expected ~4x', t2/t1));
end

%% ============================================ realistic throughput
function testRealisticArcThroughput(tc)
% categories: speed
% A sanity floor, deliberately loose so it cannot go flaky: one day of 5 s
% data (17280 epochs) across 2 obs types must accumulate in reasonable time.
% The purpose is to catch an accidental O(m^2) regression, not to benchmark.
    cfg = tc.TestData.cfg;
    cfg.mArc       = 4320;      % 6 h arcs at 5 s
    cfg.arcsPerDay = 4;         % -> 17280 epochs/day/obsType

    tic;
    sol = runAcc(cfg, 1, 'Backend','qr', 'BatchSize',1e9, 'Verbose',false);
    el  = toc;

    nObs = cfg.mArc * cfg.arcsPerDay * numel(cfg.obsTypes) * numel(cfg.sats);

    tc.verifyTrue(all(isfinite(sol.x)));
    tc.verifyLessThan(el, 60, ...
        sprintf('%d observations took %.1f s -- suspiciously slow', nObs, el));

    fprintf('\n    throughput: %d obs in %.2f s (%.0f obs/s)\n', ...
            nObs, el, nObs/el);
end

%% ============================================ helpers
function W = denseWhiten(C, X)
%DENSEWHITEN  The O(b^3) reference path, for timing comparison.
    W = lsqvce.decorrFactor(C) * X;
end

function cfg = makeCfg()
    rs = RandStream('mt19937ar','Seed',42);
    cfg.nc   = 6;
    cfg.nl   = 4;
    cfg.sats = [1 2 3];
    cfg.x_c  = randn(rs, cfg.nc, 1);
    cfg.y    = containers.Map('KeyType','double','ValueType','any');
    for s = cfg.sats, cfg.y(s) = randn(rs, cfg.nl, 1); end
    cfg.obsTypes     = {'pos','kbr'};
    cfg.arcsPerDay   = 2;
    cfg.mArc         = 400;
    cfg.dt           = 5;
    cfg.seed         = 42;
    cfg.collinear    = false;
    cfg.collinearEps = 1e-7;
    cfg.cov.pos = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',0.7);
    cfg.cov.kbr = struct('type','diag','sigma',0.2);
end

function f = binder(cfg, s)
    f = @(d,ot,a,uo) lsqvce.synthLoader(cfg, s, d, ot, a, uo);
end

function sol = runAcc(cfg, days, varargin)
    acc = lsqvce.BlockAccumulator(cfg.nc, varargin{:});
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg,s), days);
    end
    sol = acc.solve();
end
