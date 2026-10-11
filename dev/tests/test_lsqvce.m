function tests = test_lsqvce
%TEST_LSQVCE  Unit tests for the lsqvce toolbox.
%
%   Run with:  results = runtests('tests/test_lsqvce.m')
%
%   Every numerical tolerance here was established by cross-validating an
%   independent implementation of the same algorithms against dense textbook
%   formulae and against Monte-Carlo sampling theory. They are not guesses.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

tests = functiontests(localfunctions);
end

%% ===================================================== fixtures
function setupOnce(tc)
    seedFor(tc, 42);   % R7-45: was a bare rng(42) that leaked to later suites
    tc.TestData.mk = @makeProblem;
end

function [A, l, p, xt] = makeProblem(n, u, seed, outliers)
    if nargin < 4, outliers = []; end
    rs = RandStream('mt19937ar', 'Seed', seed);
    A  = randn(rs, n, u);
    xt = randn(rs, u, 1);
    sg = 0.5 + rand(rs, n, 1);
    p  = 1 ./ sg.^2;
    l  = A * xt + sg .* randn(rs, n, 1);
    for k = 1:size(outliers, 1)
        i = outliers(k, 1);  m = outliers(k, 2);
        l(i) = l(i) + m * sg(i);
    end
end

function B = chop(A, l, p, nb)
    n = numel(l);
    e = round(linspace(0, n, nb + 1));
    B = cell(nb, 1);
    for k = 1:nb
        idx = (e(k)+1):e(k+1);
        B{k} = {A(idx,:), l(idx), p(idx)};
    end
end

function f = srcOf(B)
    f = @(k) deal(B{k}{1}, B{k}{2}, B{k}{3});
end

function [u, src, gsrc, nb] = vcSmall()
% A small, well-conditioned 2-group varComp problem shared by the TraceMethod
% parity tests: enough redundancy that the machinery check and the exact
% components are both meaningful, small enough to stay fast.
    rs = RandStream('mt19937ar', 'Seed', 23);
    u = 5;  ng = 400;  n = 2*ng;  nb = 8;
    A = randn(rs, n, u);  xt = randn(rs, u, 1);
    g = [ones(ng,1); 2*ones(ng,1)];
    p = ones(n, 1);
    l = A*xt + randn(rs, n, 1);
    B  = chop(A, l, p, nb);
    gB = chop(A, g, p, nb);
    src  = srcOf(B);
    gsrc = @(k) gB{k}{2};
end

%% ================================================ stat kernel
function test_normCrit_textbook(tc)
% categories: precision, reference
    % The classic Baarda threshold.
    % tol: closed-form literature value; kernel's own accuracy claim (help text) bounds the error, tol adds headroom
    verifyEqual(tc, lsqvce.kernel.stat.normCrit(0.001), 3.29053, 'AbsTol', 1e-4);
end

function test_delta0_textbook(tc)
% categories: precision, reference
    % The classic MDB non-centrality, alpha=0.001 beta=0.20.
    % tol: closed-form literature value; kernel's own accuracy claim (help text) bounds the error, tol adds headroom
    verifyEqual(tc, lsqvce.kernel.stat.delta0(0.001, 0.20), 4.1321, 'AbsTol', 1e-3);
end

function test_chi2Inv_roundtrip(tc)
% categories: precision
    % The only meaningful accuracy test of a quantile is that plugging it back
    % into the CDF returns the probability you asked for.
    for r = [1 2 5 10 100 1e4 1e6]
        for p = [1e-6 0.005 0.5 0.995 1-1e-6]
            x = lsqvce.kernel.stat.chi2Inv(p, r);
            back = lsqvce.kernel.stat.chi2Cdf(x, r);
            % tol: quantile/CDF round-trip through base-MATLAB special functions; their relative error is ~1e-14, bound leaves headroom
            verifyEqual(tc, back, p, 'RelTol', 1e-8, ...
                sprintf('chi2Inv round-trip failed at r=%g p=%g', r, p));
        end
    end
end

function test_stat_contract_rejections(tc)
% categories: contract
% WP-43 (audit check 16): the stat kernel documents three rejection ids
% (badProb, badDof, sizeMismatch) and none was exercised directly. Rejected
% calls run before the body -- microseconds each. The badProb message must
% name the calling quantile (the "offending field" promise of the contract).
    % try/catch capture (WP-39 pattern): verifyError's output-argument form
    % errored on the function-based TestCase in CI (round-1 red) -- assert
    % the id AND the caller-naming message through the caught exception.
    caught = MException.empty;
    try
        lsqvce.kernel.stat.chi2Inv(1.5, 5);
    catch e
        caught = e;
    end
    tc.verifyNotEmpty(caught, 'p outside (0,1) must be rejected');
    tc.verifyEqual(caught.identifier, 'lsqvce:stat:badProb');
    tc.verifySubstring(caught.message, 'chi2Inv');    % message names the caller
    tc.verifyError(@() lsqvce.kernel.stat.chi2Inv(NaN, 5),   'lsqvce:stat:badProb');
    tc.verifyError(@() lsqvce.kernel.stat.chi2Inv(0.5, 0),   'lsqvce:stat:badDof');
    tc.verifyError(@() lsqvce.kernel.stat.tInv(0.9, -1),     'lsqvce:stat:badDof');
    tc.verifyError(@() lsqvce.kernel.stat.chi2Inv([0.1 0.5], [1;2;3]), ...
                   'lsqvce:stat:sizeMismatch');
end

function test_chi2Inv_singular_density(tc)
% categories: precision, robustness
    % r = 1 makes the chi2 density diverge at the origin. A Newton iteration
    % in x stalls here; the log-space iteration must not.
    x = lsqvce.kernel.stat.chi2Inv(1e-6, 1);
    % tol: quantile/CDF round-trip through base-MATLAB special functions; their relative error is ~1e-14, bound leaves headroom
    verifyEqual(tc, lsqvce.kernel.stat.chi2Cdf(x, 1), 1e-6, 'RelTol', 1e-9);
end

function test_tau_is_bounded_by_sqrt_r(tc)
% categories: reference
    % Pope's tau cannot exceed sqrt(r). This is the formal statement that a
    % weak network cannot detect its own blunders, and the toolbox relies on
    % it to warn that snooping has become vacuous.
    for r = [3 5 10 30]
        tcv = lsqvce.kernel.stat.tauCrit(1e-8, r);
        verifyLessThan(tc, tcv, sqrt(r) + 1e-9);
    end
end

function test_tau_tends_to_normal(tc)
% categories: precision
    % tol: tau -> normal is an ASYMPTOTIC limit at r = 1e7, not an identity;
    % the O(1/r) remainder sets the 1e-3 bar
    verifyEqual(tc, lsqvce.kernel.stat.tauCrit(0.001, 1e7), ...
                    lsqvce.kernel.stat.normCrit(0.001), 'AbsTol', 1e-3);
end

%% ================================================ accumulator
function test_redundancy_sum_identity(tc)
% categories: precision
    % sum(r_i) == n - u, exactly. This identity is the toolbox's own integrity
    % check and must hold to machine precision.
    [A, l, p] = tc.TestData.mk(400, 5, 21);
    B = chop(A, l, p, 11);
    [~, D] = flatAdjust(5, srcOf(B), 11, 'verbose', false);
    % tol: exact-path identity (same numbers recombined); house round-trip bound 1e-12/1e-10 (HANDOVER 6.1.2)
    verifyEqual(tc, D.sumR, D.r, 'RelTol', 1e-10);
end

function test_redundancy_in_unit_interval(tc)
% categories: contract
    [A, l, p] = tc.TestData.mk(300, 6, 33);
    B = chop(A, l, p, 7);
    [~, D] = flatAdjust(6, srcOf(B), 7, 'verbose', false);
    verifyGreaterThanOrEqual(tc, min(D.rNum), 0);
    verifyLessThanOrEqual(tc, max(D.rNum), 1);
end

function test_outliers_are_found(tc)
% categories: contract
    [A, l, p] = tc.TestData.mk(800, 6, 77, [17 9; 400 -8; 700 7]);
    B = chop(A, l, p, 10);
    [~, D] = flatAdjust(6, srcOf(B), 10, 'verbose', false, ...
                            'bonferroni', false);
    verifyTrue(tc, all(ismember([17 400 700], D.flagIdx)), ...
        'planted outliers were not all flagged');
end

function test_uncontrolled_observation_no_nan(tc)
% categories: precision, robustness | tags: diagnostics
    % A parameter seen by exactly one observation makes that observation
    % completely uncontrolled: r_i = 0 and w_i = 0/0. The code must clamp
    % rather than emit NaN, and must REPORT the observation as untestable.
    rs = RandStream('mt19937ar','Seed',13);
    n = 60; u = 6;
    A = randn(rs,n,u);
    A(:,6) = 0;  A(7,6) = 1;            % param 6 <- obs 7 alone
    l = randn(rs,n,1);  p = ones(n,1);
    B = chop(A, l, p, 4);
    [~, D] = flatAdjust(u, srcOf(B), 4, 'verbose', false);

    verifyTrue(tc, all(isfinite(D.w)),    'w contains NaN/Inf');
    verifyTrue(tc, all(isfinite(D.rNum)), 'r contains NaN/Inf');
    verifyGreaterThanOrEqual(tc, min(D.rNum), 0);
    verifyLessThan(tc, D.rNum(7), 1e-8);
    verifyGreaterThanOrEqual(tc, D.nUncontrolled, 1);
    verifyTrue(tc, isinf(D.mdb(7)), 'MDB of an uncontrolled obs must be Inf');
    % and the identity still holds
    % tol: redundancy partition identity sum(r_i)=n-u; fp accumulation over n terms, house bound with headroom
    verifyEqual(tc, D.sumR, D.r, 'RelTol', 1e-10);
end

function test_external_reliability_is_exact(tc)
% categories: precision, reference
    % Inject exactly the MDB into the worst observation and check that the
    % parameter really moves by the predicted number of sigmas. This validates
    % the headline reliability claim against brute force.
    rs = RandStream('mt19937ar','Seed',7);
    n = 120; u = 5;
    A = randn(rs,n,u);  l = randn(rs,n,1);  p = ones(n,1);
    B = chop(A, l, p, 3);
    [~, D] = flatAdjust(u, srcOf(B), 3, 'verbose', false);

    i = D.extMaxIdx;  j = D.extMaxPar;
    x1 = (A'*A) \ (A'*l);
    l2 = l;  l2(i) = l2(i) + D.mdb(i);
    x2 = (A'*A) \ (A'*l2);
    realised = abs(x2(j) - x1(j)) / D.sigX(j);

    % tol: independent dense solve of the same system; difference ~ cond*eps, house bounds add headroom (HANDOVER 6.1.2)
    verifyEqual(tc, realised, D.extMax, 'RelTol', 1e-7, ...
        'predicted external reliability does not match brute-force injection');
end

function test_variance_factor_unbiased(tc)
% categories: precision
    % Under a correct model E[s0^2] = sigma0^2 and sd = sqrt(2/r).
    M = 120;  n = 200;  u = 6;
    rat = zeros(M,1);
    for m = 1:M
        [A, l, p] = tc.TestData.mk(n, u, 1000 + m);
        B = chop(A, l, p, 5);
        [~, D] = flatAdjust(u, srcOf(B), 5, 'verbose', false);
        rat(m) = D.varRatio;
    end
    % tol: statistical identity on a finite sample; scatter ~ 1/sqrt(dof) of the fixture, bar sized 4/sqrt(Z)-style (HANDOVER 6.1.2)
    verifyEqual(tc, mean(rat), 1, 'AbsTol', 0.03);
    verifyEqual(tc, std(rat), sqrt(2/(n-u)), 'RelTol', 0.25);
end

function test_bonferroni_suppresses_false_alarms(tc)
% categories: contract
    % With n = 2000 and alpha = 1e-3 uncorrected, ~2 good observations are
    % flagged per realisation by chance. Bonferroni must remove essentially
    % all of them.
    nRaw = 0;  nBonf = 0;  M = 40;
    for m = 1:M
        [A, l, p] = tc.TestData.mk(2000, 6, 5000 + m);
        B = chop(A, l, p, 8);
        [~, Dr] = flatAdjust(6, srcOf(B), 8, 'verbose', false, ...
                                 'bonferroni', false);
        [~, Db] = flatAdjust(6, srcOf(B), 8, 'verbose', false, ...
                                 'bonferroni', true);
        nRaw  = nRaw  + Dr.nFlag;
        nBonf = nBonf + Db.nFlag;
    end
    verifyGreaterThan(tc, nRaw, 0.3 * 0.001 * 2000 * M);
    verifyLessThan(tc, nBonf, 0.2 * nRaw);
end

function test_report_verdict_clean(tc)
% categories: contract | tags: diagnostics
    [A, l, p] = tc.TestData.mk(500, 5, 91);
    B = chop(A, l, p, 6);
    [~, ~, V] = flatAdjust(5, srcOf(B), 6, 'verbose', false);
    verifyTrue(tc, V.integrityOK);
    verifyTrue(tc, any(strcmp(V.verdict, ...
        {'PASS', 'PASS WITH WARNINGS'})), ...
        sprintf('clean data gave verdict "%s"', V.verdict));
end

function test_report_flags_uncontrolled(tc)
% categories: contract | tags: diagnostics
    rs = RandStream('mt19937ar','Seed',31);
    n = 80; u = 5;
    A = randn(rs,n,u);
    A(:,5) = 0;  A(3,5) = 1;
    l = randn(rs,n,1);  p = ones(n,1);
    B = chop(A, l, p, 4);
    [~, ~, V] = flatAdjust(u, srcOf(B), 4, 'verbose', false);
    codes = cellfun(@(f) f.code, V.findings, 'UniformOutput', false);
    verifyTrue(tc, ismember('REL_UNCONTROLLED', codes), ...
        'an uncontrolled observation was not reported');
end
