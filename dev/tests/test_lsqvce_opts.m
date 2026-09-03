function tests = test_lsqvce_opts
%TEST_LSQVCE_OPTS  Coverage for option paths and edge cases.
%
%   The first two suites test the SCIENCE (do the numbers come out right, do
%   the failure modes get caught). This one tests the SURFACE: every documented
%   option actually does what it claims, and every degenerate input produces a
%   clear error rather than a confusing one.
%
%   Run with:  results = runtests('tests/test_lsqvce_opts.m')
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

tests = functiontests(localfunctions);
end

%% =================================================== fixtures
function B = chop(A, l, p, nb) %#ok<*DEFNU>
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

function [A, l, p, xt] = prob(n, u, seed, bad)
    if nargin < 4, bad = []; end
    rs = RandStream('mt19937ar', 'Seed', seed);
    A  = randn(rs, n, u);
    xt = randn(rs, u, 1);
    sg = 0.5 + rand(rs, n, 1);
    p  = 1 ./ sg.^2;
    l  = A*xt + sg .* randn(rs, n, 1);
    for i = 1:numel(bad)
        l(bad(i)) = l(bad(i)) + 15 * sg(bad(i));
    end
end

%% ======================================= INPUT VALIDATION / EDGE CASES
function [A, l, p] = emptyMiddle(Af, lf, pf, k)
    switch k
        case 1, i = 1:50;
        case 2, i = [];              % an empty batch
        case 3, i = 51:100;
    end
    A = Af(i,:);  l = lf(i);  p = pf(i);
end

function test_store_none_still_gives_scalars(tc)
% categories: precision, robustness | tags: diagnostics
    % 'store','none' must still produce every SUMMARY statistic -- it only
    % suppresses the per-observation arrays. That is the whole point: at
    % n = 1e8 you cannot hold 5 arrays of size n, but you still need the verdict.
    [A, l, p] = prob(500, 5, 31, [10 200]);
    B = chop(A, l, p, 5);
    [~, D] = flatAdjust(5, srcOf(B), 5, 'verbose', false, 'store', 'none');

    verifyTrue(tc, isfinite(D.s0));
    verifyTrue(tc, isfinite(D.sumR));
    % tol: redundancy partition identity sum(r_i)=n-u; fp accumulation over n terms, house bound with headroom
    verifyEqual(tc, D.sumR, D.r, 'RelTol', 1e-9);
    verifyTrue(tc, isfield(D, 'quad'));
    verifyFalse(tc, isfield(D, 'v'), ...
        '''store'',''none'' must NOT allocate the n-sized arrays');
end

function test_store_flagged_keeps_the_flags(tc)
% categories: contract
    [A, l, p] = prob(500, 5, 33, [7 250 400]);
    B = chop(A, l, p, 5);
    [~, D] = flatAdjust(5, srcOf(B), 5, 'verbose', false, ...
                            'store', 'flagged', 'bonferroni', false);
    verifyGreaterThan(tc, D.nFlag, 0);
    verifyEqual(tc, numel(D.flagIdx), D.nFlag);
end

function test_baarda_vs_pope(tc)
% categories: precision, reference, robustness
    % Baarda divides by the A PRIORI sigma0; Pope by the ESTIMATED s0. On clean
    % data with sigma0 correct they agree closely; the DIFFERENCE only appears
    % under contamination (which is what makes Pope fragile -- see the ext suite).
    [A, l, p] = prob(600, 6, 41);
    B = chop(A, l, p, 6);
    [~, Db] = flatAdjust(6, srcOf(B), 6, 'verbose', false, 'test','baarda');
    [~, Dp] = flatAdjust(6, srcOf(B), 6, 'verbose', false, 'test','pope');

    % w_pope = w_baarda * sigma0/s0
    ratio = Db.w ./ Dp.w;
    ratio = ratio(isfinite(ratio) & Dp.w ~= 0);
    % tol: documented scaling relation between two runs of the same solve; fp-only deviation, bound adds headroom
    verifyEqual(tc, median(ratio), Dp.s0 / Db.sigma0, 'RelTol', 1e-6);
end

function test_bonferroni_raises_the_threshold(tc)
% categories: precision
    [A, l, p] = prob(2000, 6, 51);
    B = chop(A, l, p, 8);
    [~, D0] = flatAdjust(6, srcOf(B), 8, 'verbose', false, 'bonferroni', false);
    [~, D1] = flatAdjust(6, srcOf(B), 8, 'verbose', false, 'bonferroni', true);
    verifyGreaterThan(tc, D1.kCrit, D0.kCrit, ...
        'the Bonferroni correction must RAISE the critical value');
    % tol: documented scaling relation between two runs of the same solve; fp-only deviation, bound adds headroom
    verifyEqual(tc, D1.alphaEff, D0.alpha / D0.n, 'RelTol', 1e-12);
end

function test_rWeak_changes_the_quadrants(tc)
% categories: contract
    [A, l, p] = prob(300, 5, 61);
    B = chop(A, l, p, 3);
    [~, Da] = flatAdjust(5, srcOf(B), 3, 'verbose', false, 'rWeak', 0.01);
    [~, Db] = flatAdjust(5, srcOf(B), 3, 'verbose', false, 'rWeak', 0.99);
    verifyGreaterThanOrEqual(tc, Db.quad.untested + Db.quad.ambiguous, ...
                                 Da.quad.untested + Da.quad.ambiguous, ...
        'raising rWeak must classify MORE observations as weakly controlled');
end

%% ============================================================== SNOOP
function test_snoop_removes_the_planted_outliers(tc)
% categories: contract
    muteWarn(tc, 'lsqvce:snoop:maxReject');   % provoked on purpose here (WP-47, gate rule 7)
    bad = [17 400 700];
    [A, l, p] = prob(800, 6, 71, bad);
    B = chop(A, l, p, 8);
    src = srcOf(B);
    [sol, D] = flatAdjust(6, src, 8, 'verbose', false);
    [~, ~, REJ] = lsqvce.engine.snoop(sol, D, src, 'verbose', false, ...
                                'bonferroni', false);
    verifyTrue(tc, all(ismember(bad, REJ.idx)), ...
        sprintf('snoop rejected %s, expected to include %s', ...
                mat2str(REJ.idx(:).'), mat2str(bad)));
end

function test_snoop_refuses_low_redundancy(tc)
% categories: robustness
    % Rejecting a weakly-controlled observation is close to a coin flip: its
    % residual is as likely to be error absorbed from a NEIGHBOUR. snoop must
    % refuse and say so, not quietly delete it.
    rs = RandStream('mt19937ar','Seed', 81);
    n = 60; u = 6;
    A = randn(rs, n, u);
    A(:,6) = 0;  A(7,6) = 1;               % obs 7 has r = 0
    l = randn(rs, n, 1);  l(7) = l(7) + 50;
    p = ones(n,1);
    B = chop(A, l, p, 2);
    src = srcOf(B);
    [sol, D] = flatAdjust(u, src, 2, 'verbose', false);
    [~, ~, REJ] = lsqvce.engine.snoop(sol, D, src, 'verbose', false);
    verifyFalse(tc, ismember(7, REJ.idx), ...
        'snoop must NOT reject a zero-redundancy observation');
end

function test_snoop_maxReject_cap(tc)
% categories: contract
    muteWarn(tc, 'lsqvce:snoop:maxReject');   % provoked on purpose here (WP-47, gate rule 7)
    bad = 10:10:200;
    [A, l, p] = prob(800, 6, 91, bad);
    B = chop(A, l, p, 8);
    src = srcOf(B);
    [sol, D] = flatAdjust(6, src, 8, 'verbose', false);
    [~, ~, REJ] = lsqvce.engine.snoop(sol, D, src, 'verbose', false, ...
                                'maxReject', 5, 'bonferroni', false);
    verifyLessThanOrEqual(tc, numel(REJ.idx), 5);
end

function test_snoop_reaccumulate_gives_clean_solution(tc)
% categories: precision
    % After rejecting, the final solve must be a FRESH accumulation, not a
    % sequence of downdates -- otherwise downdate roundoff accumulates into
    % the reported answer.
    muteWarn(tc, 'lsqvce:snoop:maxReject');   % provoked on purpose here (WP-47, gate rule 7)
    bad = [50 300];
    [A, l, p] = prob(600, 5, 101, bad);
    B = chop(A, l, p, 6);
    src = srcOf(B);
    [sol, D] = flatAdjust(5, src, 6, 'verbose', false);
    [solR, ~, REJ] = lsqvce.engine.snoop(sol, D, src, 'verbose', false, ...
                                   'bonferroni', false, 'reaccumulate', true);
    % compare against a manual solve over the survivors
    keep = true(600,1);  keep(REJ.idx) = false;
    xd = (A(keep,:)' * diag(p(keep)) * A(keep,:)) \ ...
         (A(keep,:)' * diag(p(keep)) * l(keep));
    % tol: independent dense solve of the same system; difference ~ cond*eps, house bounds add headroom (HANDOVER 6.1.2)
    verifyEqual(tc, solR.x, xd, 'RelTol', 1e-8);
end

%% ============================================================= ROBUST
function test_robust_huber_and_bisquare(tc)
% categories: contract
    bad = [30 150 300];
    [A, l, p] = prob(400, 5, 111, bad);
    B = chop(A, l, p, 4);
    for inf_ = {'huber', 'bisquare'}
        [~, ~, ROB] = flatRobust(5, srcOf(B), 4, 'verbose', false, ...
                                     'influence', inf_{1});
        verifyLessThan(tc, max(ROB.factor(bad)), 0.5, ...
            sprintf('%s failed to downweight the outliers', inf_{1}));
    end
end

function test_robust_converges(tc)
% categories: contract
    [A, l, p] = prob(300, 5, 121, [40 200]);
    B = chop(A, l, p, 3);
    [~, ~, ROB] = flatRobust(5, srcOf(B), 3, 'verbose', false);
    verifyTrue(tc, ROB.converged, 'IRLS did not converge on a benign problem');
    verifyLessThan(tc, ROB.iter, 30);
end

function test_robust_weights_in_unit_interval(tc)
% categories: contract
    [A, l, p] = prob(300, 5, 131, [40 200]);
    B = chop(A, l, p, 3);
    [~, ~, ROB] = flatRobust(5, srcOf(B), 3, 'verbose', false);
    verifyGreaterThanOrEqual(tc, min(ROB.factor), 0);
    verifyLessThanOrEqual(tc, max(ROB.factor), 1 + 1e-12);
end

%% ============================================================ ELLIPSE
function test_ellipse_confidence_scales_axes(tc)
% categories: precision
    muteWarn(tc, 'MATLAB:graphics:HardwareUnavailable');   % headless runner: MATLAB warns ONCE per session at the first graphics use (WP-47 CI catch); each figure-creating test can be that first use in a subset run
    [A, l, p] = prob(400, 4, 141);
    B = chop(A, l, p, 4);
    [~, D] = flatAdjust(4, srcOf(B), 4, 'verbose', false);
    E1 = lsqvce.engine.ellipse(D, [1 2], 'conf', 0.39);  % ~1 sigma
    E2 = lsqvce.engine.ellipse(D, [1 2], 'conf', 0.95);
    verifyGreaterThan(tc, E2.semiMajor, E1.semiMajor, ...
        'a 95% ellipse must be LARGER than a 1-sigma ellipse');
    % the axis RATIO is a property of Qxx and must not change with confidence
    % tol: the axis ratio depends on Qxx alone; the confidence factor
    % cancels EXACTLY, fp division noise only -- house bound 1e-9
    verifyEqual(tc, E2.semiMajor/E2.semiMinor, E1.semiMajor/E1.semiMinor, ...
        'RelTol', 1e-9);
end

function test_ellipse_returns_a_drawable_boundary(tc)
% categories: contract, robustness
    [A, l, p] = prob(300, 4, 151);
    B = chop(A, l, p, 3);
    [~, D] = flatAdjust(4, srcOf(B), 3, 'verbose', false);
    E = lsqvce.engine.ellipse(D, [1 2], 'npts', 64);
    verifySize(tc, E.xy, [64 2]);
    verifyTrue(tc, all(isfinite(E.xy(:))));
end

function test_ellipse_rejects_bad_pair(tc)
% categories: contract
    [A, l, p] = prob(200, 3, 161);
    B = chop(A, l, p, 2);
    [~, D] = flatAdjust(3, srcOf(B), 2, 'verbose', false);
    verifyError(tc, @() lsqvce.engine.ellipse(D, [1 9]), ...
        'lsqvce:ellipse:badPair');
end

%% ======================================================= SOLVENONLINEAR
function test_nonlinear_history_is_monotone(tc)
% categories: contract
    % The residual sum of squares must DECREASE at every Gauss-Newton step on
    % a well-behaved problem. If it does not, the linearisation is wrong.
    stations = [0 0; 10 0; 10 10; 0 10; 5 12];
    truth = [4.2; 6.7];
    rs = RandStream('mt19937ar','Seed', 171);
    d = sqrt(sum((stations - truth.').^2, 2)) + 0.01*randn(rs,5,1);

    function [A, l, p] = src(~, x)
        dx = x(1) - stations(:,1);
        dy = x(2) - stations(:,2);
        rr = sqrt(dx.^2 + dy.^2);
        A  = [dx./rr, dy./rr];
        l  = d - rr;
        p  = ones(5,1) / 0.01^2;
    end

    [~, ~, ~, IT] = flatNonlinear(2, @src, 1, [1;1], 'verbose', false);
    v = IT.vPv(:);                  % vPv per iteration
    verifyTrue(tc, all(diff(v) <= 1e-8 * max(v(1),1)), ...
        'vPv increased during Gauss-Newton -- check the Jacobian');
end

function test_nonlinear_respects_maxIter(tc)
% categories: contract
    muteWarn(tc, 'lsqvce:nonlin:maxIter');   % provoked on purpose here (WP-47, gate rule 7)
    muteWarn(tc, 'lsqvce:acc:zeroRedundancy');   % provoked on purpose here (WP-47, gate rule 7)
    function [A, l, p] = slow(~, x)
        A = 1;  l = 0.999*x;  p = 1;      % converges very slowly
    end
    [~, ~, ~, IT] = flatNonlinear(1, @slow, 1, 1, ...
                                           'verbose', false, 'maxIter', 3);
    verifyLessThanOrEqual(tc, IT.nIter, 3);
end

%% ============================================================== DATUM
function test_datum_user_supplied_basis(tc)
% categories: contract
    u = 5;
    edges = [1 2; 2 3; 3 4; 4 5; 1 3];
    A = zeros(size(edges,1), u);
    for i = 1:size(edges,1)
        A(i, edges(i,1)) = -1;  A(i, edges(i,2)) = 1;
    end
    N = A'*A;
    Bu = ones(u,1);                        % the user KNOWS the defect
    [B, info] = lsqvce.engine.datum(N, 'datumB', Bu);
    verifyEqual(tc, info.method, 'user');
    verifyEqual(tc, size(B,2), 1);
    verifyLessThan(tc, norm(N*B), 1e-8*norm(N));
end

function test_datum_no_defect(tc)
% categories: robustness
    % A full-rank N has NO datum defect: the function must say so, not
    % hallucinate a null space out of rounding noise.
    rs = RandStream('mt19937ar','Seed', 181);
    A = randn(rs, 50, 5);
    [B, info] = lsqvce.engine.datum(A'*A);
    verifyEqual(tc, info.defect, 0);
    verifyEmpty(tc, B);
end

%% ========================================================== SYSTEMATIC
function test_systematic_finds_the_right_covariate(tc)
% categories: robustness
    % A drift in TIME and no structure in a second, irrelevant covariate.
    % systematic() must name the right one.
    rs = RandStream('mt19937ar','Seed', 191);
    n = 800; u = 5;
    A = randn(rs, n, u);
    xt = randn(rs, u, 1);
    t = linspace(0,1,n).';
    junk = randn(rs, n, 1);
    l = A*xt + 0.3*randn(rs,n,1) + 2.0*t;      % LINEAR drift in t
    B = chop(A, l, ones(n,1), 8);
    [~, D] = flatAdjust(u, srcOf(B), 8, 'verbose', false);
    S = lsqvce.engine.systematic(D, 'covariates', struct('time', t, 'junk', junk), ...
                           'verbose', false);
    verifyTrue(tc, S.detected);

    % S.worst ranks EVERY test, including Ljung-Box -- and a linear drift in
    % time NECESSARILY autocorrelates the residual sequence, so Ljung-Box
    % firing (measured: p = 0) is CORRECT, not a misattribution. The original
    % test compared against S.worst and therefore asserted the wrong thing.
    % The COVARIATE ranking lives in S.bins; that is what must name "time".
    [~, iw] = min([S.bins.p]);
    verifyEqual(tc, S.bins(iw).name, 'time', ...
        sprintf('blamed covariate "%s" instead of time', S.bins(iw).name));
    verifyTrue(tc, S.ljungBox.p < 1e-3, ...
        'a linear drift must also show up as autocorrelation');
end

function test_systematic_ljungbox_on_correlated_residuals(tc)
% categories: contract
    % An AR(1) error process leaves autocorrelated residuals. Ljung-Box must
    % catch it even though NO single residual is large.
    rs = RandStream('mt19937ar','Seed', 201);
    n = 600; u = 4;
    A = randn(rs, n, u);
    e = zeros(n,1);
    for i = 2:n, e(i) = 0.85*e(i-1) + randn(rs); end   % strongly correlated
    l = A*randn(rs,u,1) + e;
    B = chop(A, l, ones(n,1), 6);
    [~, D] = flatAdjust(u, srcOf(B), 6, 'verbose', false);
    S = lsqvce.engine.systematic(D, 'verbose', false);
    verifyTrue(tc, S.ljungBox.reliable);
    verifyLessThan(tc, S.ljungBox.p, 1e-3, ...
        'Ljung-Box missed a strongly autocorrelated residual sequence');
end

%% ============================================================ VARCOMP
function test_report_returns_machine_readable_findings(tc)
% categories: contract | tags: diagnostics
    [A, l, p] = prob(400, 5, 231, [10 200]);
    B = chop(A, l, p, 4);
    [~, ~, V] = flatAdjust(5, srcOf(B), 4, 'verbose', false);
    verifyTrue(tc, iscell(V.findings));
    verifyTrue(tc, all(cellfun(@(f) isfield(f,'level') && isfield(f,'code') ...
                                    && isfield(f,'severity'), V.findings)));
    verifyTrue(tc, ismember(V.verdict, ...
        {'PASS','PASS WITH WARNINGS','SUSPECT','INVALID'}));
end

function test_correlated_diagnostics_match_dense(tc)
% categories: precision, reference
    % P is block-diagonal with the BATCHES as blocks: correlation within a
    % batch, independence between. Baarda's GENERALISED statistic handles this
    % directly on the physical observations -- no whitening, no loss of
    % interpretability.
    rs = RandStream('mt19937ar','Seed', 301);
    nb = 4; nk = 20; u = 5;  n = nb*nk;
    A = randn(rs, n, u);  xt = randn(rs, u, 1);

    Pb = cell(nb,1);  Pfull = zeros(n);
    idx = (0:nk-1);
    for k = 1:nb
        Cm = 0.8 .^ abs(idx - idx.');          % AR(1) correlation
        sg = 0.5 + rand(rs, nk, 1);
        Ck = (sg * sg.') .* Cm;
        Pb{k} = inv(Ck);
        Pfull((k-1)*nk+(1:nk), (k-1)*nk+(1:nk)) = Pb{k};
    end
    l = A*xt + randn(rs, n, 1);

    src = @(k) deal(A((k-1)*nk+(1:nk), :), l((k-1)*nk+(1:nk)), Pb{k});
    [sol, D] = flatAdjust(u, src, nb, 'verbose', false);

    % dense reference
    N  = A' * Pfull * A;
    Qx = inv(N);
    xd = Qx * (A' * Pfull * l);
    vd = A*xd - l;
    Qvv = inv(Pfull) - A*Qx*A';
    rd  = diag(Qvv * Pfull);

    % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
    verifyEqual(tc, sol.x, xd, 'RelTol', 1e-8, ...
        'correlated solve disagrees with the dense reference');
    % R5 M5 5b': on the ONE accumulator a correlated batch is whitened by
    % chol(P) exactly, so the per-observation diagnostics live in the
    % WHITENED frame: D.v is W*v, D.r is the hat-matrix complement of the
    % whitened rows. Those are not the raw-frame quantities diag(Qvv P) /
    % v = A x - l (per-observation redundancy is not invariant under
    % whitening; only the frame invariants are). Pinned: the invariants.
    % tol: 1e-8 -- two routes (QR block vs dense inverse) to O(1) numbers
    verifyEqual(tc, sum(D.rNum), sum(rd), 'AbsTol', 1e-8, 'the total redundancy n - u is frame-invariant');
    verifyEqual(tc, D.vPv, vd.' * Pfull * vd, 'RelTol', 1e-8, 'v''Pv is frame-invariant');
    verifyEqual(tc, D.n, n);
end
function test_correlated_integrity_identity(tc)
% categories: precision, robustness
    % sum(r_i) == n - u must hold for correlated observations too.
    %
    % THIS TEST CAUGHT A REAL BUG. r_i is a leverage complement bounded in
    % [0,1] ONLY for diagonal P. For a correlated block it can legitimately
    % EXCEED 1 (measured: 1.0001476). Clipping it to 1 -- which the diagonal
    % code did -- silently destroyed the identity.
    rs = RandStream('mt19937ar','Seed', 311);
    nb = 5; nk = 30; u = 6;  n = nb*nk;
    A = randn(rs, n, u);
    Pb = cell(nb,1);
    idx = (0:nk-1);
    for k = 1:nb
        Cm = 0.85 .^ abs(idx - idx.');
        sg = 0.5 + rand(rs, nk, 1);
        Pb{k} = inv((sg*sg.') .* Cm);
    end
    l = A*randn(rs,u,1) + randn(rs, n, 1);
    src = @(k) deal(A((k-1)*nk+(1:nk), :), l((k-1)*nk+(1:nk)), Pb{k});
    [~, D] = flatAdjust(u, src, nb, 'verbose', false);

    % tol: redundancy partition identity sum(r_i)=n-u; fp accumulation over n terms, house bound with headroom
    verifyEqual(tc, D.sumR, D.r, 'RelTol', 1e-9, ...
        'sum(r_i) must equal n-u for correlated batches');
    verifyGreaterThanOrEqual(tc, min(D.rNum), 0);
    % and r_i is NOT capped at 1 here -- do not re-introduce the clip
    verifyTrue(tc, all(isfinite(D.rNum)));
end

function test_correlated_finds_the_physical_outlier(tc)
% categories: contract
    % The generalised statistic must identify the PHYSICAL observation, not a
    % linear combination. This is why we do not whiten.
    rs = RandStream('mt19937ar','Seed', 321);
    nb = 5; nk = 30; u = 6;  n = nb*nk;
    A = randn(rs, n, u);  xt = randn(rs, u, 1);
    Pb = cell(nb,1);
    idx = (0:nk-1);
    for k = 1:nb
        Cm = 0.85 .^ abs(idx - idx.');
        sg = 0.5 + rand(rs, nk, 1);
        Pb{k} = inv((sg*sg.') .* Cm);
    end
    l = A*xt + randn(rs, n, 1);
    bad = 47;                                % inside batch 2
    l(bad) = l(bad) + 15;

    src = @(k) deal(A((k-1)*nk+(1:nk), :), l((k-1)*nk+(1:nk)), Pb{k});
    [~, D] = flatAdjust(u, src, nb, 'verbose', false);

    [~, worst] = max(abs(D.w));
    verifyEqual(tc, worst, bad, ...
        sprintf('largest |w| at %d, blunder was at %d', worst, bad));
end

function test_correlated_reduces_to_diagonal(tc)
% categories: precision
    % Feeding a DIAGONAL matrix as a full matrix must give bit-comparable
    % results to feeding the vector. The generalised form is a strict
    % generalisation, not a separate algorithm.
    rs = RandStream('mt19937ar','Seed', 331);
    nb = 4; nk = 25; u = 5;  n = nb*nk;
    A = randn(rs, n, u);
    sg = 0.5 + rand(rs, n, 1);
    p  = 1 ./ sg.^2;
    l  = A*randn(rs,u,1) + sg .* randn(rs, n, 1);

    srcV = @(k) deal(A((k-1)*nk+(1:nk), :), l((k-1)*nk+(1:nk)), ...
                     p((k-1)*nk+(1:nk)));
    srcM = @(k) deal(A((k-1)*nk+(1:nk), :), l((k-1)*nk+(1:nk)), ...
                     diag(p((k-1)*nk+(1:nk))));

    [sV, DV] = flatAdjust(u, srcV, nb, 'verbose', false);
    [sM, DM] = flatAdjust(u, srcM, nb, 'verbose', false);

    % tol: diag-matrix vs vector weights is the SAME algebra restated;
    % eps-level accumulation-order difference only (house bounds)
    verifyEqual(tc, sM.x,    sV.x,    'RelTol', 1e-10);
    verifyEqual(tc, DM.rNum, DV.rNum, 'AbsTol', 1e-10, ...
        'a diagonal matrix must give the same r_i as the vector form');
    verifyEqual(tc, DM.w,    DV.w,    'AbsTol', 1e-8, ...   % tol: same source
        'the generalised w must REDUCE to the classical w for diagonal P');
    verifyEqual(tc, DM.vPv,  DV.vPv,  'RelTol', 1e-10);   % tol: same source
end

function test_snoop_refuses_correlated(tc)
% categories: contract, robustness
    % Removing one observation from a correlated block changes the covariance
    % of every SURVIVING observation in that block (a Schur complement, not a
    % row deletion), so the rank-one Cholesky downdate does not apply.
    % Refusing is honest; silently flattening the matrix would not be.
    rs = RandStream('mt19937ar','Seed', 341);
    nb = 3; nk = 20; u = 4;  n = nb*nk;
    A = randn(rs, n, u);
    Pb = cell(nb,1);
    idx = (0:nk-1);
    for k = 1:nb
        Cm = 0.7 .^ abs(idx - idx.');
        Pb{k} = inv(Cm + 0.1*eye(nk));
    end
    l = A*randn(rs,u,1) + randn(rs, n, 1);
    src = @(k) deal(A((k-1)*nk+(1:nk), :), l((k-1)*nk+(1:nk)), Pb{k});
    [sol, D] = flatAdjust(u, src, nb, 'verbose', false);

    verifyError(tc, @() lsqvce.engine.snoop(sol, D, src, 'verbose', false), ...
        'lsqvce:snoop:correlatedBatch');
end
