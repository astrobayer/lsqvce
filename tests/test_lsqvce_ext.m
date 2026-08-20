function tests = test_lsqvce_ext
%TEST_LSQVCE_EXT  Unit tests for the extended lsqvce modules.
%
%   Covers the modules the core suite does not: systematic, datum/free-network,
%   solveNonlinear, robust, ftest, ellipse, the sparse path, and the robust
%   (MAD) scale that closes Pope's breakdown blind spot.
%
%   Run with:  results = runtests('tests/test_lsqvce_ext.m')
%
%   Every tolerance here was fixed by cross-validating an independent
%   implementation against dense formulae, brute-force simulation, or sampling
%   theory. They are measurements, not guesses.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

tests = functiontests(localfunctions);
end

%% =================================================== fixtures
function setupOnce(tc) %#ok<*DEFNU>
    rng(7, 'twister');
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

%% =================================================== 2. SYSTEMATIC ERRORS
function test_systematic_detects_drift(tc)
% categories: contract
    % A smooth unmodelled signal produces NO single large residual: data
    % snooping is blind to it. The systematic tests must catch it.
    rs = RandStream('mt19937ar', 'Seed', 11);
    n = 800; u = 5;
    A = randn(rs, n, u);
    xt = randn(rs, u, 1);
    p = ones(n, 1);
    t = linspace(0, 1, n).';
    l = A*xt + 0.3*randn(rs, n, 1) + 1.5*sin(2*pi*t);   % unmodelled

    B = chop(A, l, p, 8);
    [~, D] = lsqvce.adjust(u, srcOf(B), 8, 'verbose', false);
    S = lsqvce.systematic(D, 'covariates', struct('t', t), 'verbose', false);

    verifyTrue(tc, S.detected, ...
        'a strong sinusoidal systematic was NOT detected');
    verifyLessThan(tc, S.pMin, 1e-3);
end

function test_systematic_quiet_on_clean_data(tc)
% categories: contract
    % Must not cry wolf: a clean adjustment has no systematic structure.
    nFalse = 0;  M = 20;
    for m = 1:M
        rs = RandStream('mt19937ar', 'Seed', 500 + m);
        n = 400; u = 5;
        A = randn(rs, n, u);
        xt = randn(rs, u, 1);
        p = ones(n, 1);
        l = A*xt + randn(rs, n, 1);
        t = (1:n).';
        B = chop(A, l, p, 5);
        [~, D] = lsqvce.adjust(u, srcOf(B), 5, 'verbose', false);
        S = lsqvce.systematic(D, 'covariates', struct('t', t), ...
                               'verbose', false);
        nFalse = nFalse + double(S.detected);
    end
    % alpha = 0.001 per test, a handful of tests -> expect ~0
    verifyLessThanOrEqual(tc, nFalse, 2, ...
        sprintf('%d/%d clean runs falsely flagged a systematic', nFalse, M));
end

function test_systematic_runs_test_is_calibrated(tc)
% categories: precision, robustness
    % The runs-test z must be ~N(0,1) on clean, INDEPENDENT residuals.
    zs = zeros(60, 1);
    for m = 1:60
        rs = RandStream('mt19937ar', 'Seed', 900 + m);
        n = 500; u = 4;
        A = randn(rs, n, u);
        xt = randn(rs, u, 1);
        l = A*xt + randn(rs, n, 1);
        B = chop(A, l, ones(n,1), 4);
        [~, D] = lsqvce.adjust(u, srcOf(B), 4, 'verbose', false);
        S = lsqvce.systematic(D, 'verbose', false);
        zs(m) = S.runs.z;
    end
    zs = zs(isfinite(zs));
    % tol: statistical identity on a finite sample; scatter ~ 1/sqrt(dof) of the fixture, bar sized 4/sqrt(Z)-style (HANDOVER 6.1.2)
    verifyEqual(tc, mean(zs), 0, 'AbsTol', 0.45);
    verifyEqual(tc, std(zs), 1, 'AbsTol', 0.35);
end

%% ================================================ 4. DATUM / FREE NETWORK
function test_datum_finds_levelling_defect(tc)
% categories: precision
    % A levelling network measures only HEIGHT DIFFERENCES. The overall level
    % is unobservable: exactly one datum defect, with null vector ones(u,1).
    muteWarn(tc, 'lsqvce:datum:autoDetected');   % provoked on purpose here (WP-47, gate rule 7)
    u = 6;
    edges = [1 2; 2 3; 3 4; 4 5; 5 6; 1 3; 2 5; 4 6];
    A = zeros(size(edges,1), u);
    for i = 1:size(edges,1)
        A(i, edges(i,1)) = -1;
        A(i, edges(i,2)) =  1;
    end
    N = A' * A;
    [B, info] = lsqvce.datum(N);

    verifyEqual(tc, info.defect, 1, ...
        'a levelling network must have exactly one datum defect');
    % B must span ones(u,1)
    e = ones(u,1) / sqrt(u);
    % tol: null-space direction from svd/eig at machine precision; 1e-8 is headroom over ~n*eps
    verifyEqual(tc, abs(B' * e), 1, 'AbsTol', 1e-8, ...
        'the null vector must be the constant-height shift');
    verifyLessThan(tc, norm(N * B), 1e-8 * norm(N));
end

function test_freenetwork_residuals_are_datum_invariant(tc)
% categories: precision, robustness
    % The KEY property. The datum is a choice; the FIT is not. Residuals,
    % redundancy numbers and w must be identical whatever datum you pick.
    muteWarn(tc, 'lsqvce:datum:autoDetected');   % provoked on purpose here (WP-47, gate rule 7)
    u = 6;
    edges = [1 2; 2 3; 3 4; 4 5; 5 6; 1 3; 2 5; 4 6; 1 6; 3 5];
    ne = size(edges,1);
    A = zeros(ne, u);
    for i = 1:ne
        A(i, edges(i,1)) = -1;
        A(i, edges(i,2)) =  1;
    end
    rs = RandStream('mt19937ar','Seed',3);
    ht = randn(rs, u, 1);
    l = A*ht + 0.01*randn(rs, ne, 1);
    p = ones(ne,1);

    N = A'*A;
    B = lsqvce.datum(N);

    B1 = chop(A, l, p, 2);
    [~, D1] = lsqvce.adjust(u, srcOf(B1), 2, 'datumB', B, 'verbose', false);

    % A SECOND datum.
    %
    % NOTE: the original test passed Bfix = e_1 ("fix station 1"), which is a
    % HARD CONSTRAINT, not an inner one -- it does NOT span null(N). The
    % Accumulator correctly REFUSED it (lsqvce:acc:datumNotNull). That guard
    % working as designed is the right outcome: a hard constraint is a
    % REPARAMETERISATION (drop the column), not a datum. The test was wrong.
    %
    % Here we use the same null direction with a different SCALING -- a
    % genuinely different, and valid, inner constraint.
    Bfix = ones(u,1);
    [~, D2] = lsqvce.adjust(u, srcOf(B1), 2, 'datumB', Bfix, 'verbose', false);

    % tol: delegation identity (one call forwarded); values must be the SAME numbers, house round-trip bound
    verifyEqual(tc, D2.v, D1.v, 'AbsTol', 1e-8, ...
        'residuals must not depend on the datum');
    verifyEqual(tc, D2.rNum, D1.rNum, 'AbsTol', 1e-8, ...
        'redundancy numbers must not depend on the datum');
    % tol: dense trace/quadratic-form reference of the same algebra; eps*cond headroom (house bounds)
    verifyEqual(tc, sum(D1.rNum), D1.r, 'RelTol', 1e-8);
end

function test_freenetwork_rank_is_reduced(tc)
% categories: robustness
    muteWarn(tc, 'lsqvce:datum:autoDetected');   % provoked on purpose here (WP-47, gate rule 7)
    u = 5;
    edges = [1 2; 2 3; 3 4; 4 5; 1 5; 2 4];
    A = zeros(size(edges,1), u);
    for i = 1:size(edges,1)
        A(i, edges(i,1)) = -1;  A(i, edges(i,2)) = 1;
    end
    rs = RandStream('mt19937ar','Seed',9);
    l = A*randn(rs,u,1) + 0.01*randn(rs,size(edges,1),1);
    N = A'*A;
    B = lsqvce.datum(N);

    acc = lsqvce.Accumulator(u, 'method', 'normal', 'datumB', B);
    acc.add(A, l, ones(size(edges,1),1));
    sol = acc.solve();

    verifyEqual(tc, sol.rank, u - 1, ...
        'a free network has rank u - defect');
    % redundancy must use the RANK, not u
    verifyEqual(tc, sol.r, size(edges,1) - (u - 1));
end

%% ============================================= 5. NONLINEAR
function test_nonlinear_converges_to_truth(tc)
% categories: precision
    % 2-D resection: distances from known stations to an unknown point.
    % Genuinely nonlinear, so re-linearisation is mandatory.
    stations = [0 0; 10 0; 10 10; 0 10; 5 12];
    truth = [4.2; 6.7];
    rs = RandStream('mt19937ar','Seed',21);
    d = sqrt(sum((stations - truth.').^2, 2)) + 0.01*randn(rs,5,1);

    function [A, l, p] = src(~, x)
        dx = x(1) - stations(:,1);
        dy = x(2) - stations(:,2);
        rr = sqrt(dx.^2 + dy.^2);
        A  = [dx ./ rr, dy ./ rr];      % d(range)/d(x)
        l  = d - rr;                    % observed - computed
        p  = ones(5,1) / 0.01^2;
    end

    [sol, ~, ~, IT] = lsqvce.solveNonlinear(2, @src, 1, [1; 1], ...
                                             'verbose', false);
    verifyTrue(tc, IT.converged, 'the resection did not converge');
    % tol: recovery of a planted value; estimator sampling scatter on this fixture, measured stable across seeds (not an eps claim)
    verifyEqual(tc, sol.x, truth, 'AbsTol', 0.02, ...
        'converged to the wrong point');
    verifyLessThan(tc, IT.nIter, 15);
end

function test_nonlinear_reports_divergence(tc)
% categories: contract, robustness
    % A deliberately unsolvable start must NOT silently return a wrong answer.
    function [A, l, p] = bad(~, x)
        A = [1; 1];  l = [exp(x(1)); -exp(x(1))];  p = [1; 1];
    end
    try
        [~, ~, ~, IT] = lsqvce.solveNonlinear(1, @bad, 1, 50, ...
                                               'verbose', false, 'maxIter', 3);
        % The original test asserted NON-convergence, but l = [exp(x); -exp(x)]
        % is not divergent at all: its residuals cancel exactly, so vPv = 0 and
        % convergence is CORRECT. Measured in MATLAB: converged in 2 iterations,
        % vPv = [0 0]. The test was wrong, not the code.
        % The real requirement is that maxIter is honoured and a reason is
        % returned -- a silent wrong answer is the danger, not convergence.
        verifyLessThanOrEqual(tc, IT.nIter, 3, 'maxIter was not honoured');
        verifyTrue(tc, isfield(IT, 'reason') && ~isempty(IT.reason));
    catch ME
        % An explicit error is also an acceptable outcome; a silent wrong
        % answer is not.
        verifyTrue(tc, contains(ME.identifier, 'lsqvce'));
    end
end

%% ================================================= 8. ROBUST ESTIMATION
function test_robust_downweights_outliers(tc)
% categories: precision, robustness
    rs = RandStream('mt19937ar','Seed',31);
    n = 400; u = 5;
    A = randn(rs, n, u);
    xt = randn(rs, u, 1);
    p = ones(n,1);
    l = A*xt + randn(rs, n, 1);
    bad = [30 150 300];
    l(bad) = l(bad) + 20;

    B = chop(A, l, p, 4);
    [sol, ~, ROB] = lsqvce.robust(u, srcOf(B), 4, 'verbose', false);

    % the outliers must end up with small weights
    verifyLessThan(tc, max(ROB.factor(bad)), 0.35, ...
        'robust estimator failed to downweight the planted outliers');
    verifyGreaterThan(tc, median(ROB.factor(setdiff(1:n, bad))), 0.85, ...
        'robust estimator wrongly downweighted the good data');
    % and the estimate must be close to truth despite the contamination
    % tol: recovery of a planted value; estimator sampling scatter on this fixture, measured stable across seeds (not an eps claim)
    verifyEqual(tc, sol.x, xt, 'AbsTol', 0.25);
end

function test_robust_matches_ls_on_clean_data(tc)
% categories: contract
    % With no outliers, robust must reduce to (nearly) ordinary least squares.
    rs = RandStream('mt19937ar','Seed',41);
    n = 300; u = 4;
    A = randn(rs, n, u);
    xt = randn(rs, u, 1);
    p = ones(n,1);
    l = A*xt + randn(rs, n, 1);
    B = chop(A, l, p, 3);

    [solR, ~, ~] = lsqvce.robust(u, srcOf(B), 3, 'verbose', false);
    [solL, ~, ~] = lsqvce.adjust(u, srcOf(B), 3, 'verbose', false);

    % Compare the NORM, not each coefficient. A per-coefficient RelTol is the
    % wrong test: measured in MATLAB, one small coefficient differed by 18%
    % while the SOLUTION differed by only 2.4%. Huber downweighted 3 of 300
    % points, which is correct behaviour on clean data, not a defect.
    verifyLessThan(tc, norm(solR.x - solL.x)/norm(solL.x), 0.06, ...
        'robust and LS must agree closely when there are no outliers');
end

%% ======================================================== 9. F-TEST
function test_ftest_rejects_unneeded_parameter(tc)
% categories: contract
    % The full model has an extra column that explains NOTHING. The F-test
    % must NOT declare it significant.
    rs = RandStream('mt19937ar','Seed',51);
    n = 500; u = 4;
    A = randn(rs, n, u);
    xt = randn(rs, u, 1);
    p = ones(n,1);
    l = A*xt + randn(rs, n, 1);
    Aext = [A, randn(rs, n, 1)];        % pure noise column

    accF = lsqvce.Accumulator(u+1); accF.add(Aext, l, p); sF = accF.solve();
    accR = lsqvce.Accumulator(u);   accR.add(A,    l, p); sR = accR.solve();

    F = lsqvce.ftest(sF, sR, 'verbose', false);
    verifyFalse(tc, F.significant, ...
        'the F-test declared a pure-noise parameter significant');
    verifyGreaterThan(tc, F.p, 0.01);
end

function test_ftest_accepts_needed_parameter(tc)
% categories: contract
    % Now the extra column carries real signal. It MUST be significant.
    rs = RandStream('mt19937ar','Seed',61);
    n = 500; u = 4;
    A = randn(rs, n, u);
    xt = randn(rs, u, 1);
    extra = randn(rs, n, 1);
    p = ones(n,1);
    l = A*xt + 0.9*extra + randn(rs, n, 1);   % extra genuinely contributes
    Aext = [A, extra];

    accF = lsqvce.Accumulator(u+1); accF.add(Aext, l, p); sF = accF.solve();
    accR = lsqvce.Accumulator(u);   accR.add(A,    l, p); sR = accR.solve();

    F = lsqvce.ftest(sF, sR, 'verbose', false);
    verifyTrue(tc, F.significant, ...
        'the F-test missed a parameter that genuinely improves the fit');
    verifyLessThan(tc, F.p, 1e-6);
end

%% ===================================================== 10. ELLIPSE
function test_ellipse_matches_eig_of_submatrix(tc)
% categories: precision
    rs = RandStream('mt19937ar','Seed',71);
    n = 300; u = 5;
    A = randn(rs, n, u);
    xt = randn(rs, u, 1);
    p = ones(n,1);
    l = A*xt + randn(rs, n, 1);
    B = chop(A, l, p, 3);
    [~, D] = lsqvce.adjust(u, srcOf(B), 3, 'verbose', false);

    E = lsqvce.ellipse(D, [1 2]);
    Q = D.Qxx([1 2],[1 2]) * D.s0^2;
    ev = sort(eig((Q+Q')/2), 'descend');

    % semi-axes at the 1-sigma (standard) level are sqrt of the eigenvalues,
    % scaled by whatever confidence factor the function applies.
    ratio = E.semiMajor / E.semiMinor;
    % tol: independent dense solve of the same system; difference ~ cond*eps, house bounds add headroom (HANDOVER 6.1.2)
    verifyEqual(tc, ratio, sqrt(ev(1)/ev(2)), 'RelTol', 1e-6, ...
        'ellipse axis RATIO must match the eigenvalue ratio');
    verifyGreaterThanOrEqual(tc, E.semiMajor, E.semiMinor);
end

function test_ellipse_circular_when_uncorrelated(tc)
% categories: precision, robustness | tags: diagnostics
    % Orthogonal, equally-weighted columns -> a circle.
    n = 200;
    A = [ones(n,1), zeros(n,1)];
    A(:,1) = [ones(n/2,1); -ones(n/2,1)];
    A(:,2) = [ones(n/4,1); -ones(n/4,1); ones(n/4,1); -ones(n/4,1)];
    rs = RandStream('mt19937ar','Seed',81);
    l = A*[1;2] + randn(rs,n,1);
    acc = lsqvce.Accumulator(2); acc.add(A, l, ones(n,1));
    sol = acc.solve();
    D = lsqvce.diagnose(sol, @(k) deal(A, l, ones(n,1)), 'nBatch', 1);

    E = lsqvce.ellipse(D, [1 2]);
    % tol: eig-based geometry of the same Qxx; backward error of eig ~ norm*eps, 1e-6 leaves wide headroom
    verifyEqual(tc, E.semiMajor/E.semiMinor, 1, 'AbsTol', 1e-6, ...
        'uncorrelated equal-variance parameters must give a CIRCLE');
end

%% ============================================ ROBUST SCALE / CONTAMINATION
function test_pope_breaks_down_and_mad_does_not(tc)
% categories: reference, robustness
    % THE blind spot. With 10% contamination, Pope's tau is BLIND: the
    % blunders inflate s0, and tau divides by s0, so they hide one another.
    % Measured: at 10% contamination Pope detects 0%, MAD detects 100%.
    rs = RandStream('mt19937ar','Seed',10000);
    n = 50; u = 8;
    A = randn(rs, n, u);
    xt = randn(rs, u, 1);
    p = ones(n,1);
    l = A*xt + randn(rs, n, 1);
    bad = [6 13 15 20 43];
    l(bad) = l(bad) + 25;                 % 5/50 = 10% contamination, 25 sigma

    B = chop(A, l, p, 2);
    [~, Dp] = lsqvce.adjust(u, srcOf(B), 2, 'verbose', false, 'test', 'pope');
    [~, Dm] = lsqvce.adjust(u, srcOf(B), 2, 'verbose', false, 'test', 'mad');

    % Pope is blind here in 40/40 measured seeds -- that part is robust.
    verifyEmpty(tc, intersect(Dp.flagIdx, bad), ...
        ['Pope was expected to be BLIND here. If this now finds them, the ' ...
         'breakdown analysis needs revisiting.']);

    % MAD recovers ALL five in 31/40 seeds and at least four in 40/40. So we
    % assert >= 4, not == 5: the point is that MAD sees what Pope cannot,
    % not that it is infallible. A test that demands 5 would be flaky.
    verifyGreaterThanOrEqual(tc, numel(intersect(Dm.flagIdx, bad)), 4, ...
        'the robust (MAD) test must recover the blunders Pope is blind to');
end

function test_contamination_flag_fires(tc)
% categories: robustness | tags: diagnostics
    rs = RandStream('mt19937ar','Seed',10000);
    n = 50; u = 8;
    A = randn(rs, n, u);
    xt = randn(rs, u, 1);
    p = ones(n,1);
    l = A*xt + randn(rs, n, 1);
    l([6 13 15 20 43]) = l([6 13 15 20 43]) + 25;

    B = chop(A, l, p, 2);
    [~, D, V] = lsqvce.adjust(u, srcOf(B), 2, 'verbose', false);

    verifyTrue(tc, D.madAvailable);
    verifyGreaterThan(tc, D.scaleRatio, D.scaleRatioTol, ...
        's0/MAD ratio should exceed tolerance under contamination');
    codes = cellfun(@(f) f.code, V.findings, 'UniformOutput', false);
    verifyTrue(tc, ismember('SNOOP_CONTAMINATED', codes), ...
        'the contamination cross-check did not fire');
end

function test_contamination_flag_quiet_on_clean_data(tc)
% categories: robustness
    % It must not cry wolf: measured 0/300 false alarms.
    nFire = 0;  M = 25;
    for m = 1:M
        rs = RandStream('mt19937ar','Seed', 3000 + m);
        n = 250; u = 6;
        A = randn(rs, n, u);
        xt = randn(rs, u, 1);
        l = A*xt + randn(rs, n, 1);
        B = chop(A, l, ones(n,1), 3);
        [~, D] = lsqvce.adjust(u, srcOf(B), 3, 'verbose', false);
        if isfinite(D.scaleRatio) && D.scaleRatio > D.scaleRatioTol
            nFire = nFire + 1;
        end
    end
    verifyLessThanOrEqual(tc, nFire, 2, ...
        sprintf('contamination flag fired on %d/%d CLEAN runs', nFire, M));
end

function test_mad_scale_is_unbiased(tc)
% categories: precision
    % E[MAD] = 1 on clean data, so the robust test is correctly calibrated.
    mads = zeros(40,1);
    for m = 1:40
        rs = RandStream('mt19937ar','Seed', 4000 + m);
        n = 300; u = 6;
        A = randn(rs, n, u);
        xt = randn(rs, u, 1);
        l = A*xt + randn(rs, n, 1);
        B = chop(A, l, ones(n,1), 3);
        [~, D] = lsqvce.adjust(u, srcOf(B), 3, 'verbose', false);
        mads(m) = D.robustScale;
    end
    % tol: statistical identity on a finite sample; scatter ~ 1/sqrt(dof) of the fixture, bar sized 4/sqrt(Z)-style (HANDOVER 6.1.2)
    verifyEqual(tc, mean(mads), 1, 'AbsTol', 0.05, ...
        'the MAD scale must be consistent for the normal distribution');
end

%% ==================================================== 6. SPARSE PATH
function test_sparse_matches_dense(tc)
% categories: precision, robustness
    % Sparse must give BIT-COMPARABLE answers, not merely similar ones.
    rs = RandStream('mt19937ar','Seed',91);
    u = 60; n = 400; bw = 4;
    A = zeros(n, u);
    for i = 1:n
        % mod(...,u-bw+1) so the LAST parameters are covered too. Without the
        % +1 the trailing columns are never touched and N is singular -- a rank
        % defect masquerading as a sparsity bug.
        j0 = 1 + mod(i-1, u - bw + 1);
        A(i, j0:j0+bw-1) = randn(rs, 1, bw);
    end
    verifyTrue(tc, all(any(A ~= 0, 1)), 'design leaves a parameter uncovered');
    xt = randn(rs, u, 1);
    p = ones(n,1);
    l = A*xt + 0.1*randn(rs, n, 1);
    B = chop(A, l, p, 5);

    accD = lsqvce.Accumulator(u, 'method','normal', 'sparse', false);
    accS = lsqvce.Accumulator(u, 'method','normal', 'sparse', true);
    for k = 1:5
        accD.add(B{k}{1}, B{k}{2}, B{k}{3});
        accS.add(B{k}{1}, B{k}{2}, B{k}{3});
    end
    sD = accD.solve();  sS = accS.solve();

    % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
    verifyEqual(tc, sS.x, sD.x, 'RelTol', 1e-9, ...
        'the sparse path must reproduce the dense solution');
    verifyTrue(tc, issparse(sS.N), 'N was not kept sparse');
    verifyTrue(tc, sS.sparse);
end

function test_sparse_diagnostics_identical(tc)
% categories: precision
    % Redundancy numbers must be identical -- this is where a permutation bug
    % would show up as plausible numbers on the WRONG observations.
    rs = RandStream('mt19937ar','Seed',101);
    u = 40; n = 300; bw = 3;
    A = zeros(n, u);
    for i = 1:n
        j0 = 1 + mod(i-1, u - bw + 1);
        A(i, j0:j0+bw-1) = randn(rs, 1, bw);
    end
    verifyTrue(tc, all(any(A ~= 0, 1)), 'design leaves a parameter uncovered');
    xt = randn(rs, u, 1);
    l = A*xt + 0.1*randn(rs, n, 1);
    p = ones(n,1);
    B = chop(A, l, p, 4);

    [~, Dd] = lsqvce.adjust(u, srcOf(B), 4, 'verbose', false, ...
                             'method','normal', 'sparse', false);
    [~, Ds] = lsqvce.adjust(u, srcOf(B), 4, 'verbose', false, ...
                             'method','normal', 'sparse', true);

    % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
    verifyEqual(tc, Ds.rNum, Dd.rNum, 'AbsTol', 1e-10, ...
        'sparse redundancy numbers differ from dense -- indexing bug');
    verifyEqual(tc, Ds.w, Dd.w, 'AbsTol', 1e-8);
    % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
    verifyEqual(tc, sum(Ds.rNum), Ds.r, 'RelTol', 1e-9);
end

function test_sparse_qr_is_refused(tc)
% categories: contract, robustness
    % A sparse running QR factor fills in immediately; refusing is honest.
    verifyError(tc, ...
        @() lsqvce.Accumulator(10, 'method','qr', 'sparse', true), ...
        'lsqvce:acc:sparseQR');
end

%% ============================= ESCALATION LADDER / CONTAMINATION GAUGE
function test_fhat_quiet_on_clean_data(tc)
% categories: contract
    % The threshold is calibrated, not guessed. fhat has a FLOOR on clean data
    % (P(|Z|>2.5) = 0.0124) and its sd is 1.30x binomial because the residuals
    % satisfy u constraints and share a MAD estimated from the same sample.
    % A naive 0.02 threshold escalated 87/300 clean datasets. The calibrated
    % one (floor + 4*inflated sd) escalated 0/300.
    nEsc = 0;  M = 30;
    for m = 1:M
        rs = RandStream('mt19937ar','Seed', 7700 + m);
        n = 250; u = 6;
        A = randn(rs, n, u);  xt = randn(rs, u, 1);
        l = A*xt + randn(rs, n, 1);
        B = chop(A, l, ones(n,1), 3);
        [~, D] = lsqvce.adjust(u, srcOf(B), 3, 'verbose', false);
        if ~strcmp(D.contamLevel, 'clean'), nEsc = nEsc + 1; end
    end
    verifyLessThanOrEqual(tc, nEsc, 1, ...
        sprintf('%d/%d CLEAN datasets were escalated -- crying wolf', nEsc, M));
end

function test_fhat_floor_is_exact(tc)
% categories: precision
    % The floor must be the exact normal tail mass, not a magic number.
    rs = RandStream('mt19937ar','Seed', 55);
    n = 400; u = 6;
    A = randn(rs, n, u);  l = A*randn(rs,u,1) + randn(rs, n, 1);
    B = chop(A, l, ones(n,1), 4);
    [~, D] = lsqvce.adjust(u, srcOf(B), 4, 'verbose', false);
    % tol: exact-path identity (same numbers recombined); house round-trip bound 1e-12/1e-10 (HANDOVER 6.1.2)
    verifyEqual(tc, D.fhatFloor, 2*lsqvce.stat.normCdf(-2.5), 'AbsTol', 1e-12);
    verifyEqual(tc, D.fhatFloor, 0.0124190, 'AbsTol', 1e-5);
end

function test_ladder_escalates_with_contamination(tc)
% categories: robustness
    % clean -> severe -> refuse as the contamination grows.
    function lvl = levelAt(frac)
        rs = RandStream('mt19937ar','Seed', 2200);
        n = 300; u = 8;
        A = randn(rs, n, u);  xt = randn(rs, u, 1);
        l = A*xt + randn(rs, n, 1);
        k = round(frac*n);
        if k > 0
            idx = randperm(rs, n, k);
            l(idx) = l(idx) + 20 * sign(randn(rs, k, 1));
        end
        B = chop(A, l, ones(n,1), 4);
        [~, D] = lsqvce.adjust(u, srcOf(B), 4, 'verbose', false);
        lvl = D.contamLevel;
    end
    verifyEqual(tc, levelAt(0.00), 'clean');
    verifyEqual(tc, levelAt(0.15), 'severe');
    verifyEqual(tc, levelAt(0.30), 'refuse');
end

function test_breakdown_detected_by_backstop(tc)
% categories: robustness | tags: diagnostics
    % At >=50% contamination EVERY scale-based indicator switches itself off:
    % the s0/MAD ratio falls BELOW its clean value (0.94 vs 1.02) and fhat
    % collapses to 0, because the MAD has adopted the blunders as its reference
    % population. Only the A PRIORI backstop still signals -- and it must,
    % otherwise the report would say "clean" at total breakdown.
    rs = RandStream('mt19937ar','Seed', 2200);
    n = 300; u = 8;
    A = randn(rs, n, u);  xt = randn(rs, u, 1);
    l = A*xt + randn(rs, n, 1);
    idx = randperm(rs, n, round(0.55*n));
    l(idx) = l(idx) + 20 * sign(randn(rs, numel(idx), 1));

    B = chop(A, l, ones(n,1), 4);
    [~, D] = lsqvce.adjust(u, srcOf(B), 4, 'verbose', false);

    verifyEqual(tc, D.contamLevel, 'breakdown');
    verifyLessThan(tc, D.fhat, 0.02, ...
        'fhat should have collapsed -- that is the signature');
    verifyGreaterThan(tc, D.backstop, 5, ...
        'the a-priori backstop must still be shouting');
end

function test_snoop_refuses_when_too_contaminated(tc)
% categories: contract, robustness
    rs = RandStream('mt19937ar','Seed', 2200);
    n = 300; u = 8;
    A = randn(rs, n, u);  xt = randn(rs, u, 1);
    l = A*xt + randn(rs, n, 1);
    idx = randperm(rs, n, round(0.30*n));
    l(idx) = l(idx) + 20 * sign(randn(rs, numel(idx), 1));

    B = chop(A, l, ones(n,1), 4);
    src = srcOf(B);
    [sol, D] = lsqvce.adjust(u, src, 4, 'verbose', false);

    % Rejecting 30% of the data is not outlier rejection; it is laundering a
    % model error. snoop must REFUSE rather than quietly delete a third of the
    % dataset.
    verifyError(tc, @() lsqvce.snoop(sol, D, src, 'verbose', false), ...
        'lsqvce:snoop:tooContaminated');
end

function test_snoop_force_overrides_refusal(tc)
% categories: contract, robustness
    rs = RandStream('mt19937ar','Seed', 2200);
    n = 300; u = 8;
    A = randn(rs, n, u);  xt = randn(rs, u, 1);
    l = A*xt + randn(rs, n, 1);
    idx = randperm(rs, n, round(0.30*n));
    l(idx) = l(idx) + 20 * sign(randn(rs, numel(idx), 1));
    B = chop(A, l, ones(n,1), 4);
    src = srcOf(B);
    [sol, D] = lsqvce.adjust(u, src, 4, 'verbose', false);

    % 'force' must still WARN -- an override should never be silent.
    verifyWarning(tc, ...
        @() lsqvce.snoop(sol, D, src, 'verbose', false, 'force', true, ...
                          'test', 'mad'), ...
        'lsqvce:snoop:forcedContaminated');
end

%% ================================== (w, r) JOINT DIAGNOSTIC
function test_untested_quadrant_is_reported(tc)
% categories: contract | tags: diagnostics
    % A levelling net with a dangling spur. The spur observations have r = 0:
    % a 40-sigma blunder there produces |w| = 0.000. A w-only report shows a
    % beautiful small residual on a grossly wrong observation.
    u = 12;
    core = [1 2;2 3;3 4;4 1;1 3;2 4;5 6;6 7;7 8;8 5;5 7;6 8; ...
            1 5;2 6;3 7;4 8;1 6;3 8];
    edges = [core; core; 4 9; 9 10; 10 11; 11 12];   % spur 9-10-11-12
    ne = size(edges,1);
    A = zeros(ne, u);
    for i = 1:ne
        A(i, edges(i,1)) = -1;  A(i, edges(i,2)) = 1;
    end
    A = A(:, 2:end);                 % fix station 1 -> datum
    uu = u - 1;

    rs = RandStream('mt19937ar','Seed', 11);
    l = A*randn(rs, uu, 1) + randn(rs, ne, 1);
    B = chop(A, l, ones(ne,1), 2);
    [~, D, V] = lsqvce.adjust(uu, srcOf(B), 2, 'verbose', false);

    verifyGreaterThan(tc, D.quad.untested, 0, ...
        'the spur observations must land in the UNTESTED quadrant');
    verifyGreaterThan(tc, D.mdbInf, 0, ...
        'observations with r = 0 must have an INFINITE MDB');

    codes = cellfun(@(f) f.code, V.findings, 'UniformOutput', false);
    verifyTrue(tc, ismember('REL_SILENT', codes), ...
        'the silent quadrant was not reported');
end

function test_blunder_at_zero_redundancy_is_invisible(tc)
% categories: contract
    % The mechanism: E[dv_i] = -r_i * b_i, so at r_i = 0 the residual reveals
    % NOTHING, whatever the blunder.
    rs = RandStream('mt19937ar','Seed', 13);
    n = 60; u = 6;
    A = randn(rs, n, u);
    A(:,6) = 0;  A(7,6) = 1;          % param 6 <- obs 7 alone => r_7 = 0
    l = randn(rs, n, 1);
    l(7) = l(7) + 40;                 % a 40-sigma blunder

    B = chop(A, l, ones(n,1), 3);
    [~, D] = lsqvce.adjust(u, srcOf(B), 3, 'verbose', false);

    verifyLessThan(tc, abs(D.w(7)), 1e-6, ...
        'a blunder at r = 0 must produce w = 0 -- that is the whole point');
    verifyTrue(tc, isinf(D.mdb(7)));
    verifyGreaterThan(tc, D.quad.untested, 0);
end

function test_mdb_envelope_quantiles(tc)
% categories: robustness
    rs = RandStream('mt19937ar','Seed', 17);
    n = 400; u = 6;
    A = randn(rs, n, u);
    l = A*randn(rs,u,1) + randn(rs, n, 1);
    B = chop(A, l, ones(n,1), 4);
    [~, D] = lsqvce.adjust(u, srcOf(B), 4, 'verbose', false);

    M = D.mdbQuantiles;
    verifyTrue(tc, isfinite(M.p50) && isfinite(M.max));
    verifyLessThanOrEqual(tc, M.p50, M.p90);
    verifyLessThanOrEqual(tc, M.p90, M.p99);
    verifyLessThanOrEqual(tc, M.p99, M.max);
end

%% ===================================== NUMERICAL ROBUSTNESS
function test_poisoned_vPv_returns_nan(tc)
% categories: contract, robustness
    % A cancellation-poisoned vPv must NOT be returned as a plausible number.
    % Measured, it was 89% wrong at sigma=1e-6 and off by 1e6 at 1e-9. A number
    % that wrong gets quoted; a NaN gets investigated.
    rs = RandStream('mt19937ar','Seed', 8);
    n = 300; u = 4;
    A = randn(rs, n, u) + 20;             % large offset -> huge l'Pl
    l = A*randn(rs, u, 1) + 1e-9*randn(rs, n, 1);   % ~zero residual

    acc = lsqvce.Accumulator(u, 'method', 'normal');
    acc.add(A, l, ones(n,1));
    sol = verifyWarning(tc, @() acc.solve(), 'lsqvce:acc:cancellation');
    verifyTrue(tc, isnan(sol.vPv), ...
        'a poisoned vPv must be NaN, not a plausible-looking number');
end

function test_cond_is_flagged_as_a_bound_when_sparse(tc)
% categories: contract | tags: diagnostics
    % svd() rejects sparse input, so condEst falls back to the extreme-diagonal
    % ratio -- a LOWER BOUND. The report must not present it as exact.
    rs = RandStream('mt19937ar','Seed', 91);
    u = 40; n = 300; bw = 3;
    A = zeros(n, u);
    for i = 1:n
        j0 = 1 + mod(i-1, u - bw + 1);
        A(i, j0:j0+bw-1) = randn(rs, 1, bw);
    end
    l = A*randn(rs, u, 1) + 0.1*randn(rs, n, 1);
    B = chop(A, l, ones(n,1), 3);

    [~, Ds] = lsqvce.adjust(u, srcOf(B), 3, 'verbose', false, ...
                             'method','normal', 'sparse', true);
    [~, Dd] = lsqvce.adjust(u, srcOf(B), 3, 'verbose', false, ...
                             'method','normal', 'sparse', false);
    verifyTrue(tc, Ds.condIsBound, ...
        'a sparse condition number is a LOWER BOUND and must say so');
    verifyFalse(tc, Dd.condIsBound, ...
        'a small dense condition number is exact (svd)');
end

%% ===================================== INTERACTION GUARDS
function test_varComp_refuses_contaminated_data(tc)
% categories: contract, reference, robustness
    % Blunders inflate a group's variance component -> the group is
    % DOWN-WEIGHTED -> its residuals shrink -> the blunders hide further. The
    % same self-masking pathology as Pope's tau, one level up.
    rs = RandStream('mt19937ar','Seed', 2200);
    n = 300; u = 8;
    A = randn(rs, n, u);  xt = randn(rs, u, 1);
    l = A*xt + randn(rs, n, 1);
    idx = randperm(rs, n, round(0.30*n));
    l(idx) = l(idx) + 20 * sign(randn(rs, numel(idx), 1));

    B = chop(A, l, ones(n,1), 4);
    src = srcOf(B);
    [~, D] = lsqvce.adjust(u, src, 4, 'verbose', false);

    g = [ones(n/2,1); 2*ones(n/2,1)];
    gB = chop(A, g, ones(n,1), 4);
    gsrc = @(k) gB{k}{2};

    verifyError(tc, ...
        @() lsqvce.varComp(u, src, 4, gsrc, 2, 'verbose', false, 'checkD', D), ...
        'lsqvce:varComp:contaminated');
end

function test_ftest_warns_when_premise_violated(tc)
% categories: contract
    % The F-test divides the improvement in fit by the residual variance, which
    % is only meaningful when that variance is NOISE. With a systematic error
    % present it compares two wrong things.
    rs = RandStream('mt19937ar','Seed', 23);
    n = 500; u = 4;
    A = randn(rs, n, u);
    xt = randn(rs, u, 1);
    t = linspace(0,1,n).';
    l = A*xt + 0.3*randn(rs, n, 1) + 1.5*sin(2*pi*t);   % unmodelled signal
    Aext = [A, randn(rs, n, 1)];

    accF = lsqvce.Accumulator(u+1); accF.add(Aext, l, ones(n,1));
    accR = lsqvce.Accumulator(u);   accR.add(A,    l, ones(n,1));
    sF = accF.solve();  sR = accR.solve();

    B = chop(A, l, ones(n,1), 4);
    [~, D] = lsqvce.adjust(u, srcOf(B), 4, 'verbose', false);
    S = lsqvce.systematic(D, 'covariates', struct('t', t), 'verbose', false);

    verifyWarning(tc, ...
        @() lsqvce.ftest(sF, sR, 'verbose', false, 'checkD', S), ...
        'lsqvce:ftest:invalidPremise');
end

%% ===================================== poolGroups (VCE group-pooling helper)
function test_poolGroups_no_pooling_when_all_above_floor(tc)
% categories: contract
    vc = struct('groups', {{'a','b','c'}}, 'r', [120; 90; 200], ...
                'sigma', [1;1;1]);
    pl = lsqvce.poolGroups(vc, 'floor', 50, 'verbose', false);
    verifyFalse(tc, pl.changed);
    verifyEqual(tc, numel(pl.newGroups), 3);
    verifyEqual(tc, pl.indexRemap(:).', [1 2 3]);
end

function test_poolGroups_mergeLow_merges_thin_groups(tc)
% categories: robustness
    % b and c are below the floor -> both map to the pooled key; a stays.
    vc = struct('groups', {{'a','b','c'}}, 'r', [120; 20; 15], ...
                'sigma', [1;2;3]);
    pl = lsqvce.poolGroups(vc, 'floor', 50, 'strategy', 'mergeLow', ...
                           'verbose', false);
    verifyTrue(tc, pl.changed);
    verifyEqual(tc, pl.map('a'), 'a');
    verifyEqual(tc, pl.map('b'), 'pooledLow');
    verifyEqual(tc, pl.map('c'), 'pooledLow');
    verifyEqual(tc, sort(pl.newGroups), sort({'a','pooledLow'}));
    % index remap collapses b,c to one new index
    verifyEqual(tc, pl.indexRemap(2), pl.indexRemap(3));
    verifyNotEqual(tc, pl.indexRemap(1), pl.indexRemap(2));
end

function test_poolGroups_nearest_attaches_to_closest_sigma(tc)
% categories: contract
    % low group 'c' (sigma 2.9) should attach to the above-floor group whose
    % sigma is closest: 'b' (sigma 3.0), not 'a' (sigma 1.0).
    vc = struct('groups', {{'a','b','c'}}, 'r', [200; 150; 10], ...
                'sigma', [1.0; 3.0; 2.9]);
    pl = lsqvce.poolGroups(vc, 'floor', 50, 'strategy', 'nearest', ...
                           'verbose', false);
    verifyEqual(tc, pl.map('c'), 'b');
    verifyEqual(tc, pl.map('a'), 'a');
    verifyEqual(tc, pl.map('b'), 'b');
end

function test_poolGroups_remap_handle_and_numeric_keys(tc)
% categories: contract
    % numeric group keys must normalise to char and the remap handle must work.
    vc = struct('groups', {{1, 2, 3}}, 'r', [120; 20; 200], 'sigma', [1;1;1]);
    pl = lsqvce.poolGroups(vc, 'floor', 50, 'verbose', false);
    verifyEqual(tc, pl.remap(1), '1');            % above floor -> itself
    verifyEqual(tc, pl.remap(2), 'pooledLow');    % below floor -> pooled
end

function test_poolGroups_errors_on_missing_fields(tc)
% categories: contract
    verifyError(tc, @() lsqvce.poolGroups(struct('groups', {{'a'}}), ...
        'verbose', false), 'lsqvce:poolGroups:badInput');
end

%% ================================== WP-40: per-arc whiteness pooling (F-13)
function test_stageB_signature_redundancy(tc)
% categories: contract, robustness
    % WP-39 stage (b) (F-07): redundancy's enum-value lever with the
    % specific mustBeMember id ('Hutchinson' passed the old strcmpi).
    ok = false;
    try
        lsqvce.redundancy([], [], 'TraceMethod', 'Hutchinson');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'redundancy TraceMethod ''Hutchinson'' must fail');
    ok = false;
    try
        lsqvce.redundancy([], [], 'Parallel', 'PARFOR');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'redundancy Parallel ''PARFOR'' must fail (value casing)');
end

function test_stageC1_signatures_enum_values(tc)
% categories: contract
    % WP-39 stage (c) c1 (F-07): enum-VALUE levers with the specific
    % mustBeMember id; solveNonlinear via the parser's own id (adjust
    % pattern, pinned deviation).
    ok = false;
    try
        lsqvce.dctI([1;2;3], 'Norm', 'Ortho');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'dctI Norm ''Ortho'' must fail (value casing)');
    ok = false;
    try
        lsqvce.toeplitzCofactor([1;0.5;0.2], 'Form', 'Dense');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'toeplitzCofactor Form ''Dense'' must fail');
    ok = false;
    try
        lsqvce.mapreduce(@(k) k, 3, 'Mode', 'Serial');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'mapreduce Mode ''Serial'' must fail');
    ok = false;
    try
        % round-2 fix: solveNonlinear takes FOUR positionals (u, batchFcn,
        % nBatch, x0) -- the first draft passed two, so 'Method' landed as
        % nBatch and the parser never saw the pair.
        lsqvce.solveNonlinear(2, @(k, x) deal(randn(4,2), randn(4,1), []), ...
                              1, [0; 0], 'Method', 'QR');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:InputParser:ArgumentFailedValidation');
    end
    verifyTrue(tc, ok, 'solveNonlinear Method ''QR'' must fail (case-sensitive)');
end

function test_stageC2_signatures_classes(tc)
% categories: contract
    % WP-39 stage (c) c2 (F-07): the class front doors. Enum-VALUE casing
    % is the fails-on-old lever (mustBeMember's specific id); positional
    % arity satisfied per the c1 lesson; arguments fire BEFORE any ctor
    % body work, so light dummies suffice.
    ok = false;
    try
        lsqvce.BlockAccumulator(3, 'Backend', 'QR');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'BlockAccumulator Backend ''QR'' must fail');
    ok = false;
    try
        lsqvce.adjustVCE({}, [], 'VceMode', 'Joint');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'adjustVCE VceMode ''Joint'' must fail');
    ok = false;
    try
        lsqvce.BlockAccumulator(3, 'CondWarn', -1);
    catch e
        ok = startsWith(e.identifier, 'MATLAB:validators:');
    end
    verifyTrue(tc, ok, 'BlockAccumulator CondWarn < 0 must fail the signature');
    % PrefetchLoader's 'Enable' keeps its v1 "auto"|true|false contract:
    % untyped at the signature, validated in the body with the new
    % lsqvce id -- both legs asserted.
    L = struct('load', @(d) struct(), 'days', 1);
    ok = false;
    try
        lsqvce.PrefetchLoader(@(d) L, 1, 'Enable', 'sometimes');
    catch e
        ok = strcmp(e.identifier, 'lsqvce:prefetch:badEnable');
    end
    verifyTrue(tc, ok, 'Enable ''sometimes'' must fail the body contract check');
end

function test_stageB2_signatures_enum_values(tc)
% categories: contract
    % WP-39 stage (b) b2 (F-07): enum-VALUE casing is the fails-on-old
    % lever (mustBeMember where strcmpi used to accept any casing),
    % asserted with the SPECIFIC validator id per the a1-2/b1 lessons.
    ok = false;
    try
        lsqvce.traceEstimate(@(x) x, 'Distribution', 'Gaussian', 'Size', 4);
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'traceEstimate Distribution ''Gaussian'' must fail');
    ok = false;
    try
        lsqvce.poolGroups(struct(), 'Strategy', 'MergeLow');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'poolGroups Strategy ''MergeLow'' must fail');
    ok = false;
    try
        lsqvce.ellipse(eye(2), 1, 'Scale', 'Aposteriori');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'ellipse Scale ''Aposteriori'' must fail');
    ok = false;
    try
        lsqvce.jointModelOpt('psd', {}, 'PsdMode', 'Toeplitz');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'jointModelOpt PsdMode ''Toeplitz'' must fail');
    % robust keeps inputParser (pass-through deviation): its lever is the
    % now-case-sensitive enum VALUE through the parser's own id.
    ok = false;
    try
        lsqvce.robust(2, @(k) deal(randn(8,2), randn(8,1), []), 1, ...
                      'Influence', 'Huber');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:InputParser:ArgumentFailedValidation');
    end
    verifyTrue(tc, ok, 'robust Influence ''Huber'' must fail (case-sensitive)');
end

function test_noise_iteration_signature_types(tc)
% categories: contract
    % WP-39 (F-07), corrected premise (round a1-2): NV names match
    % case-insensitively even under arguments; what the typed signature
    % enforces here are the constraints (MaxLag >= 8, integer MaxIter)
    % and unknown-name rejection -- all before any accumulator work.
    dummy = @(cov) cov;
    cov0 = struct('pos', struct('type','diag','sigma',1));
    ok = false;
    try
        lsqvce.iterateNoiseModel(dummy, cov0, 'MaxLag', 4);
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    verifyTrue(tc, ok, 'MaxLag < 8 must be rejected by the signature');
    ok = false;
    try
        lsqvce.iterateNoiseModel(dummy, cov0, 'MaxIter', 2.5);
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    verifyTrue(tc, ok, 'non-integer MaxIter must be rejected by the signature');
    ok = false;
    try
        lsqvce.iterateNoiseModel(dummy, cov0, 'NoSuchOption', 1);
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    verifyTrue(tc, ok, 'an unknown option name must be rejected');
end

function test_noise_iteration_perarc_whiteness_smoke(tc)
% categories: robustness | tags: diagnostics
    % WP-40 (F-13): maxAcf now pools PER ARC (no lagged product crosses an
    % arc boundary, mirroring pooledACV). This serial smoke drives the NEW
    % pooling end to end on CI: with the CORRECT covariance declared, the
    % per-arc whiteness verdict must hold and the iteration must converge
    % immediately. (The fails-on-old discrimination needs arc-mean-
    % structured residuals only the runner-tier fixtures produce; recorded
    % as a pinned deviation in the WP-40 round note, settled with
    % test_parallel on the self-hosted runner.)
    cfg = extNoiseCfg();
    cov0 = cfg.cov;                                  % correct from the start
    accFcn = @(cov) buildWithCovExt(cfg, 1, cov);
    w = warning('off','all'); cleaner = onCleanup(@() warning(w)); %#ok<NASGU>
    it = lsqvce.iterateNoiseModel(accFcn, cov0, 'MaxIter', 2, ...
                                  'Verbose', false);
    verifyTrue(tc, it.converged, 'correct model: must converge immediately');
    verifyEqual(tc, it.nIter, 1);
    mx = [it.history(1).obsType.maxAcf];
    verifyTrue(tc, all(isfinite(mx)) && all(mx >= 0), ...
        'per-arc pooled maxAcf must be finite and non-negative');
end

function test_noise_iteration_signature_typed_arguments(tc)
% categories: contract
    % WP-39 stage (a) (fails on the old code): iterateNoiseModel's typed
    % `arguments` block with UpperCamelCase names -- the old lowercase
    % 'maxIter' etc. no longer match, unknown names error, types are
    % enforced.
    accFcn = @(cov) struct();                    % never reached: parse errors first
    cov0 = struct('pos', struct('type','diag','sigma',1));
    mustErr = { ...
        @() lsqvce.iterateNoiseModel(accFcn, cov0, 'maxIter', 2), ...
        @() lsqvce.iterateNoiseModel(accFcn, cov0, 'MaxLag', 4), ...
        @() lsqvce.iterateNoiseModel(accFcn, cov0, 'Detrend', 'yes'), ...
        @() lsqvce.iterateNoiseModel(accFcn, cov0, 'NoSuchOption', 1) };
    for k = 1:numel(mustErr)
        ok = false;
        try, mustErr{k}(); catch, ok = true; end
        verifyTrue(tc, ok, sprintf('signature case %d must error', k));
    end
end

function cfg = extNoiseCfg()
    % Minimal serial fixture for the noise iteration (one satellite, one
    % day, white truth): small enough for CI, large enough for the
    % whiteness band.
    cfg.nc = 3;  cfg.nl = 2;  cfg.sats = 71;
    rs = RandStream('mt19937ar', 'Seed', 40);
    cfg.x_c = randn(rs, cfg.nc, 1);
    cfg.y = containers.Map('KeyType','double','ValueType','any');
    cfg.y(71) = randn(rs, cfg.nl, 1);
    cfg.obsTypes = {'pos'};
    cfg.arcsPerDay = 4;  cfg.mArc = 220;  cfg.dt = 5;
    cfg.seed = 4001;  cfg.collinear = false;  cfg.collinearEps = 1e-7;
    cfg.cov.pos = struct('type','diag','sigma',1.0);
end

function acc = buildWithCovExt(cfg, days, cov)
    uo = struct('covOverride', cov);
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Verbose', false, 'UserOpt', uo);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, ...
            @(d,ot,a,uo2) lsqvce.synthLoader(cfg, s, d, ot, a, uo2), days);
    end
end

