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
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M7, WP-89.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
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
    rec = timePair(fft, loop, 'Repeats', 7, 'Inner', 1, 'Label', 'fft/loop lag sums');
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
    rec = timePair(exact, probe, 'Repeats', 5, 'Inner', 1, 'Label', 'exact/probe trace');
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
