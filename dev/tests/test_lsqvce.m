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
    rng(42, 'twister');
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
    verifyEqual(tc, lsqvce.engine.stat.normCrit(0.001), 3.29053, 'AbsTol', 1e-4);
end

function test_delta0_textbook(tc)
% categories: precision, reference
    % The classic MDB non-centrality, alpha=0.001 beta=0.20.
    % tol: closed-form literature value; kernel's own accuracy claim (help text) bounds the error, tol adds headroom
    verifyEqual(tc, lsqvce.engine.stat.delta0(0.001, 0.20), 4.1321, 'AbsTol', 1e-3);
end

function test_chi2Inv_roundtrip(tc)
% categories: precision
    % The only meaningful accuracy test of a quantile is that plugging it back
    % into the CDF returns the probability you asked for.
    for r = [1 2 5 10 100 1e4 1e6]
        for p = [1e-6 0.005 0.5 0.995 1-1e-6]
            x = lsqvce.engine.stat.chi2Inv(p, r);
            back = lsqvce.engine.stat.chi2Cdf(x, r);
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
        lsqvce.engine.stat.chi2Inv(1.5, 5);
    catch e
        caught = e;
    end
    tc.verifyNotEmpty(caught, 'p outside (0,1) must be rejected');
    tc.verifyEqual(caught.identifier, 'lsqvce:stat:badProb');
    tc.verifySubstring(caught.message, 'chi2Inv');    % message names the caller
    tc.verifyError(@() lsqvce.engine.stat.chi2Inv(NaN, 5),   'lsqvce:stat:badProb');
    tc.verifyError(@() lsqvce.engine.stat.chi2Inv(0.5, 0),   'lsqvce:stat:badDof');
    tc.verifyError(@() lsqvce.engine.stat.tInv(0.9, -1),     'lsqvce:stat:badDof');
    tc.verifyError(@() lsqvce.engine.stat.chi2Inv([0.1 0.5], [1;2;3]), ...
                   'lsqvce:stat:sizeMismatch');
end

function test_chi2Inv_singular_density(tc)
% categories: precision, robustness
    % r = 1 makes the chi2 density diverge at the origin. A Newton iteration
    % in x stalls here; the log-space iteration must not.
    x = lsqvce.engine.stat.chi2Inv(1e-6, 1);
    % tol: quantile/CDF round-trip through base-MATLAB special functions; their relative error is ~1e-14, bound leaves headroom
    verifyEqual(tc, lsqvce.engine.stat.chi2Cdf(x, 1), 1e-6, 'RelTol', 1e-9);
end

function test_tau_is_bounded_by_sqrt_r(tc)
% categories: reference
    % Pope's tau cannot exceed sqrt(r). This is the formal statement that a
    % weak network cannot detect its own blunders, and the toolbox relies on
    % it to warn that snooping has become vacuous.
    for r = [3 5 10 30]
        tcv = lsqvce.engine.stat.tauCrit(1e-8, r);
        verifyLessThan(tc, tcv, sqrt(r) + 1e-9);
    end
end

function test_tau_tends_to_normal(tc)
% categories: precision
    % tol: tau -> normal is an ASYMPTOTIC limit at r = 1e7, not an identity;
    % the O(1/r) remainder sets the 1e-3 bar
    verifyEqual(tc, lsqvce.engine.stat.tauCrit(0.001, 1e7), ...
                    lsqvce.engine.stat.normCrit(0.001), 'AbsTol', 1e-3);
end

%% ================================================ accumulator
function test_matches_dense_solution(tc)
% categories: precision
    [A, l, p] = tc.TestData.mk(600, 7, 11);
    xd = (A' * diag(p) * A) \ (A' * diag(p) * l);

    for meth = {'qr', 'normal'}
        B = chop(A, l, p, 9);
        acc = lsqvce.engine.Accumulator(7, 'method', meth{1});
        for k = 1:9, acc.add(B{k}{1}, B{k}{2}, B{k}{3}); end
        sol = acc.solve();
        % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
        verifyEqual(tc, sol.x, xd, 'RelTol', 1e-9, ...
            sprintf('method %s disagrees with the dense solution', meth{1}));
    end
end

function test_invariant_to_batch_partition(tc)
% categories: vectorisation
    % The whole point of a streaming accumulator: the answer must not depend
    % on how the caller happened to chop up the data.
    [A, l, p] = tc.TestData.mk(500, 6, 3);
    xs = [];
    for nb = [1 2 7 53 500]
        B = chop(A, l, p, nb);
        acc = lsqvce.engine.Accumulator(6, 'method', 'qr');
        for k = 1:nb, acc.add(B{k}{1}, B{k}{2}, B{k}{3}); end
        xs = [xs, acc.solve().x]; %#ok<AGROW>
    end
    spread = max(max(abs(xs - xs(:,1))));
    verifyLessThan(tc, spread, 1e-10);
end

function test_vPv_cancellation_is_caught(tc)
% categories: contract
    % A tight fit makes l'Pl - b'x catastrophically cancel. The normal-equation
    % path must WARN, and the QR path must not be affected.
    rs = RandStream('mt19937ar', 'Seed', 8);
    n = 300; u = 4;
    A = randn(rs, n, u) + 20;          % large offset -> huge l'Pl
    xt = randn(rs, u, 1);
    l = A * xt + 1e-9 * randn(rs, n, 1);   % ~zero residual
    p = ones(n, 1);

    accN = lsqvce.engine.Accumulator(u, 'method', 'normal');
    accN.add(A, l, p);
    verifyWarning(tc, @() accN.solve(), 'lsqvce:acc:cancellation');

    % QR is immune and needs no warning.
    accQ = lsqvce.engine.Accumulator(u, 'method', 'qr');
    accQ.add(A, l, p);
    solQ = accQ.solve();
    verifyGreaterThanOrEqual(tc, solQ.vPv, 0);
end

function test_rank_deficiency_errors_clearly(tc)
% categories: contract, robustness | tags: diagnostics
    % Build an EXACTLY rank-deficient design (column 4 = col1 + col2). Use a
    % fixed stream so the test is deterministic; the solver must report the
    % rank defect via chol failure OR its rcond guard, never a hollow solve.
    rs = RandStream('mt19937ar', 'Seed', 4);
    A = randn(rs, 50, 4);
    A(:, 4) = A(:, 1) + A(:, 2);       % exact linear dependence
    l = randn(rs, 50, 1);
    acc = lsqvce.engine.Accumulator(4, 'method', 'normal');
    acc.add(A, l, ones(50,1));
    verifyError(tc, @() acc.solve(), 'lsqvce:acc:rankDeficient');
end

function test_negative_weight_rejected(tc)
% categories: contract
    acc = lsqvce.engine.Accumulator(3);
    p = ones(10,1);  p(4) = -1;
    verifyError(tc, @() acc.add(randn(10,3), randn(10,1), p), ...
                'lsqvce:acc:negativeWeight');
end

function test_nonfinite_rejected(tc)
% categories: contract, robustness
    acc = lsqvce.engine.Accumulator(3);
    A = randn(10,3);  A(2,2) = NaN;
    verifyError(tc, @() acc.add(A, randn(10,1), ones(10,1)), ...
                'lsqvce:acc:nonFiniteA');
end

function test_correlated_batch_supported(tc)
% categories: precision
    % A full SPD weight matrix for one batch must give the same answer as
    % feeding the whitened equivalent.
    rs = RandStream('mt19937ar','Seed',5);
    n = 40; u = 3;
    A = randn(rs,n,u); l = randn(rs,n,1);
    M = randn(rs,n,n);  P = M*M' + n*eye(n);      % SPD
    xd = (A'*P*A) \ (A'*P*l);

    acc = lsqvce.engine.Accumulator(u, 'method','qr');
    acc.add(A, l, P);
    % tol: independent dense solve of the same system; difference ~ cond*eps, house bounds add headroom (HANDOVER 6.1.2)
    verifyEqual(tc, acc.solve().x, xd, 'RelTol', 1e-8);
end

%% ================================================= diagnostics
function test_redundancy_sum_identity(tc)
% categories: precision
    % sum(r_i) == n - u, exactly. This identity is the toolbox's own integrity
    % check and must hold to machine precision.
    [A, l, p] = tc.TestData.mk(400, 5, 21);
    B = chop(A, l, p, 11);
    [~, D] = lsqvce.engine.adjust(5, srcOf(B), 11, 'verbose', false);
    % tol: exact-path identity (same numbers recombined); house round-trip bound 1e-12/1e-10 (HANDOVER 6.1.2)
    verifyEqual(tc, D.sumR, D.r, 'RelTol', 1e-10);
end

function test_redundancy_in_unit_interval(tc)
% categories: contract
    [A, l, p] = tc.TestData.mk(300, 6, 33);
    B = chop(A, l, p, 7);
    [~, D] = lsqvce.engine.adjust(6, srcOf(B), 7, 'verbose', false);
    verifyGreaterThanOrEqual(tc, min(D.rNum), 0);
    verifyLessThanOrEqual(tc, max(D.rNum), 1);
end

function test_outliers_are_found(tc)
% categories: contract
    [A, l, p] = tc.TestData.mk(800, 6, 77, [17 9; 400 -8; 700 7]);
    B = chop(A, l, p, 10);
    [~, D] = lsqvce.engine.adjust(6, srcOf(B), 10, 'verbose', false, ...
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
    [~, D] = lsqvce.engine.adjust(u, srcOf(B), 4, 'verbose', false);

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
    [~, D] = lsqvce.engine.adjust(u, srcOf(B), 3, 'verbose', false);

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
        [~, D] = lsqvce.engine.adjust(u, srcOf(B), 5, 'verbose', false);
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
        [~, Dr] = lsqvce.engine.adjust(6, srcOf(B), 8, 'verbose', false, ...
                                 'bonferroni', false);
        [~, Db] = lsqvce.engine.adjust(6, srcOf(B), 8, 'verbose', false, ...
                                 'bonferroni', true);
        nRaw  = nRaw  + Dr.nFlag;
        nBonf = nBonf + Db.nFlag;
    end
    verifyGreaterThan(tc, nRaw, 0.3 * 0.001 * 2000 * M);
    verifyLessThan(tc, nBonf, 0.2 * nRaw);
end

function test_second_pass_mismatch_detected(tc)
% categories: contract | tags: diagnostics
    % If the caller's second-pass batch function does not reproduce the first
    % pass, every diagnostic would be silently attributed to the wrong
    % observation. That must be an error, not a surprise.
    [A, l, p] = tc.TestData.mk(200, 4, 55);
    B = chop(A, l, p, 5);
    acc = lsqvce.engine.Accumulator(4);
    for k = 1:5, acc.add(B{k}{1}, B{k}{2}, B{k}{3}); end
    sol = acc.solve();

    bad = @(k) deal(B{k}{1}(1:end-1,:), B{k}{2}(1:end-1), B{k}{3}(1:end-1));
    verifyError(tc, @() lsqvce.engine.diagnose(sol, bad, 'nBatch', 5), ...
                'lsqvce:diagnose:batchMismatch');
end

function test_zero_redundancy_is_flagged(tc)
% categories: contract
    % n == u: the model fits exactly and nothing can be tested. The toolbox
    % must say so rather than reporting a meaningless perfect result.
    A = randn(5,5);  l = randn(5,1);
    acc = lsqvce.engine.Accumulator(5);
    acc.add(A, l, ones(5,1));
    verifyWarning(tc, @() acc.solve(), 'lsqvce:acc:zeroRedundancy');
end

%% =========================================== variance components
function test_helmert_recovers_known_weight_error(tc)
% categories: precision
    % Group 2's true sigma is 3x what we claim. Helmert must find that factor.
    rs = RandStream('mt19937ar','Seed',4);
    u = 6;  ng = 600;  nG = 2;  n = ng*nG;
    A = randn(rs,n,u);  xt = randn(rs,u,1);
    g = [ones(ng,1); 2*ones(ng,1)];
    sigTrue = [ones(ng,1); 3*ones(ng,1)];
    p = ones(n,1);                       % we CLAIM sigma = 1 for both
    l = A*xt + sigTrue .* randn(rs,n,1);

    B = chop(A, l, p, 8);
    gB = chop(A, g, p, 8);               % reuse the chopper to split g
    src = srcOf(B);
    gsrc = @(k) gB{k}{2};

    VC = lsqvce.engine.varComp(u, src, 8, gsrc, 2, 'verbose', false);
    ratio = VC.sigma(2) / VC.sigma(1);
    % tol: recovery of a planted value; estimator sampling scatter on this fixture, measured stable across seeds (not an eps claim)
    verifyEqual(tc, ratio, 3, 'RelTol', 0.12);
end

function test_adjust_typed_validators_and_enum_values(tc)
% categories: contract | tags: diagnostics
    % WP-39 round a2 (F-07): adjust keeps inputParser BY DESIGN (open
    % pass-through to diagnose) but its five own options got typed
    % validators and canonical names. Fails-on-old levers, checked with
    % the SPECIFIC parser id (a loose 'MATLAB:' prefix is satisfied by
    % body crashes -- the a1-2 lesson): the old code had NO validators,
    % so these calls sailed into the body there.
    src_ = @(k) deal(randn(8,2), randn(8,1), []);
    ok = false;
    try
        lsqvce.engine.adjust(2, src_, 1, 'Sigma0', 'abc');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:InputParser:ArgumentFailedValidation');
    end
    verifyTrue(tc, ok, 'non-numeric Sigma0 must fail the typed validator');
    ok = false;
    try
        lsqvce.engine.adjust(2, src_, 1, 'Method', 'QR');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:InputParser:ArgumentFailedValidation');
    end
    verifyTrue(tc, ok, 'Method ''QR'' must be rejected (case-sensitive enum value)');
    % and the pass-through contract survives: lowercase casing + a
    % forwarded diagnose option must still run end to end.
    w = warning('off','all'); cleaner = onCleanup(@() warning(w)); %#ok<NASGU>
    [sol, D, V] = lsqvce.engine.adjust(2, src_, 1, 'verbose', false, 'store', 'none'); %#ok<ASGLU>
    verifyTrue(tc, isstruct(sol) || isobject(sol), ...
        'lowercase casing and forwarded options must keep working');
end

function test_stageB_signature_varComp(tc)
% categories: contract
    % WP-39 stage (b) (F-07): the enum-value lever with the specific
    % mustBeMember id (the old strcmpi accepted 'Rigorous').
    b = @(k) deal(randn(8,2), randn(8,1), ones(8,1));
    g = @(k) ones(8,1);
    ok = false;
    try
        lsqvce.engine.varComp(2, b, 1, g, 1, 'Method', 'Rigorous');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'varComp Method ''Rigorous'' must fail mustBeMember');
    ok = false;
    try
        lsqvce.engine.varComp(2.5, b, 1, g, 1);
    catch e
        ok = startsWith(e.identifier, 'MATLAB:validators:');
    end
    verifyTrue(tc, ok, 'non-integer u must fail the typed positional');
end

function test_varComp_negative_component_freezes_not_nukes(tc)
% categories: contract
    % WP-40 (fixes F-11; fails on the old code, whose s(bad) = 1e-6 divided
    % the group's CUMULATIVE variance by 1e6): a two-group rig (group 2:
    % tiny true noise, 6 obs on 4 shared params) drives the rigorous
    % Helmert component of group 2 negative. The seed is searched
    % deterministically (twister; the first firing seed is stable across
    % runs), the pathology is detected via the lsqvce:varComp:negative
    % warning, and the treatment is asserted: the component is FROZEN --
    % every cumulative scale stays O(1), never scale*1e-6.
    w = warning('off','all'); cleaner = onCleanup(@() warning(w)); %#ok<NASGU>
    fired = false;
    % Round-1 lesson: with m2 = 6 the group sat below varComp's minObs
    % default (10) and was excluded from estimation entirely -- the
    % negative branch can only fire for ESTIMABLE groups. m2 = 12 with a
    % tiny true variance goes negative in ~1/3 of seeds (Python frequency
    % check in the round note), so the deterministic seed search is safe.
    for seed = 1:40
        rs = RandStream('mt19937ar', 'Seed', seed);
        u = 4;  m1 = 80;  m2 = 12;
        A1 = randn(rs, m1, u);  A2 = randn(rs, m2, u);
        xt = randn(rs, u, 1);
        l1 = A1*xt + 1.00*randn(rs, m1, 1);
        l2 = A2*xt + sqrt(0.03)*randn(rs, m2, 1);
        A  = [A1; A2];  l = [l1; l2];  p = ones(m1+m2, 1);
        g  = [ones(m1,1); 2*ones(m2,1)];
        batchFcn = @(k) deal(A, l, p);               % one batch, all obs
        groupFcn = @(k) g;
        lastwarn('', '');
        VC = lsqvce.engine.varComp(u, batchFcn, 1, groupFcn, 2, ...
                            'maxIter', 1, 'verbose', false);
        [~, wid] = lastwarn;
        if strcmp(wid, 'lsqvce:varComp:negative')
            fired = true;
            verifyTrue(tc, all(VC.scale > 1e-3), ...
                ['a negative component must FREEZE its group (scale stays ' ...
                 'O(1)) -- the old 1e-6 clamp nuked the cumulative scale']);
            verifyTrue(tc, any(VC.s == 1), ...
                'the frozen group must carry s = 1 for the pass');
            break
        end
    end
    verifyTrue(tc, fired, ...
        'fixture sanity: no seed in 1:40 produced a negative component');
end

function test_helmert_no_false_correction(tc)
% categories: precision
    % When the weights are already right, the estimator must not invent a
    % correction beyond its own sampling noise.
    rs = RandStream('mt19937ar','Seed',17);
    u = 5;  ng = 800;  nG = 2;  n = ng*nG;
    A = randn(rs,n,u);  xt = randn(rs,u,1);
    g = [ones(ng,1); 2*ones(ng,1)];
    p = ones(n,1);
    l = A*xt + randn(rs,n,1);

    B  = chop(A, l, p, 6);
    gB = chop(A, g, p, 6);
    VC = lsqvce.engine.varComp(u, srcOf(B), 6, @(k) gB{k}{2}, 2, 'verbose', false);

    ratio = VC.sigma(2) / VC.sigma(1);
    % Resolution limit is ~0.5*sqrt(2/r_g) per group; allow 3 sigma of it.
    tol = 3 * sqrt(sum(VC.resolution.^2));
    % tol: 3-sigma of the estimator's own reported resolution/SE -- self-scaling statistical bar
    verifyEqual(tc, ratio, 1, 'AbsTol', max(tol, 0.05));
end

function test_helmert_redundancy_partition(tc)
% categories: precision
    rs = RandStream('mt19937ar','Seed',9);
    u = 4;  ng = 300;  n = 2*ng;
    A = randn(rs,n,u);  l = randn(rs,n,1);  p = ones(n,1);
    g = [ones(ng,1); 2*ones(ng,1)];
    B  = chop(A, l, p, 4);
    gB = chop(A, g, p, 4);
    VC = lsqvce.engine.varComp(u, srcOf(B), 4, @(k) gB{k}{2}, 2, ...
                         'verbose', false, 'maxIter', 1);
    % tol: redundancy partition identity sum(r_i)=n-u; fp accumulation over n terms, house bound with headroom
    verifyEqual(tc, sum(VC.redund), n - u, 'RelTol', 1e-8);
end

%% ============================= varComp TraceMethod parity (Hutchinson)
function test_varComp_exact_has_no_traceEstimate(tc)
% categories: contract
    % Default 'exact' must NOT attach the estimator struct: it is opt-in.
    [u, src, gsrc, nb] = vcSmall();
    VC = lsqvce.engine.varComp(u, src, nb, gsrc, 2, 'verbose', false);
    verifyFalse(tc, isfield(VC, 'traceEstimate'), ...
        'exact TraceMethod must not attach VC.traceEstimate');
end

function test_varComp_hutchinson_check_returns_u(tc)
% categories: precision
    % Machinery check: tr(N^-1 N) = u exactly; the estimate must return ~u.
    [u, src, gsrc, nb] = vcSmall();
    VC = lsqvce.engine.varComp(u, src, nb, gsrc, 2, 'verbose', false, ...
             'TraceMethod', 'hutchinson', 'Probes', 200, 'TraceSeed', 7);
    verifyTrue(tc, isfield(VC, 'traceEstimate'));
    te = VC.traceEstimate;
    verifyEqual(tc, te.expected, u);
    % Unbiased estimator; with Z=200 probes, ~few*relSE band around u.
    % tol: Monte-Carlo probe noise, sigma ~ 1/sqrt(Probes); bar per HANDOVER 6.1.2 (4/sqrt(Z)-style, not containment)
    verifyEqual(tc, te.value, u, 'AbsTol', max(6*te.se, 0.5));
    verifyEqual(tc, te.probes, 200);
end

function test_varComp_hutchinson_does_not_perturb_components(tc)
% categories: reference
    % The exact per-group components are authoritative: turning the estimator
    % on must leave VC.s / VC.sigma / VC.redund bit-for-bit unchanged.
    [u, src, gsrc, nb] = vcSmall();
    VCe = lsqvce.engine.varComp(u, src, nb, gsrc, 2, 'verbose', false);
    VCh = lsqvce.engine.varComp(u, src, nb, gsrc, 2, 'verbose', false, ...
             'TraceMethod', 'hutchinson', 'TraceSeed', 3);
    verifyEqual(tc, VCh.s,      VCe.s);
    verifyEqual(tc, VCh.sigma,  VCe.sigma);
    verifyEqual(tc, VCh.redund, VCe.redund);
end

%% ===================================== vceHutchinson (strict-Ellmer VCE path)
function test_vceHutchinson_shape(tc)
% categories: reference
    % Returns the documented struct, marked as the estimated path, and NO
    % per-observation fields (strict Ellmer).
    [u, src, gsrc, nb] = vcSmall();
    VC = lsqvce.engine.vceHutchinson(u, src, nb, gsrc, 2, 'verbose', false, ...
             'probes', 400, 'seed', 5);
    verifyEqual(tc, VC.estimator, 'hutchinson');
    verifyEqual(tc, numel(VC.s), 2);
    verifyEqual(tc, VC.probes, 400);
    verifyTrue(tc, isfield(VC, 'redundSE'));       % estimated -> carries an SE
    verifyFalse(tc, isfield(VC, 'w'));             % no per-obs snooping outputs
    verifyFalse(tc, isfield(VC, 'r_i'));
end

function test_vceHutchinson_matches_varComp(tc)
% categories: precision
    % TIER T3 validation: on the same problem, the estimated components and
    % redundancies must agree with EXACT varComp within the sampling band. With
    % a healthy probe count the agreement is tight.
    [u, src, gsrc, nb] = vcSmall();
    VCx = lsqvce.engine.varComp(u, src, nb, gsrc, 2, 'verbose', false);
    VCh = lsqvce.engine.vceHutchinson(u, src, nb, gsrc, 2, 'verbose', false, ...
             'probes', 4000, 'seed', 11);
    % sigma scales agree to ~1% with Z=4000 on this well-conditioned problem
    % tol: Monte-Carlo probe scatter ~ 1/sqrt(Probes) on this fixture; 4/sqrt(Z)-style bar (HANDOVER 6.1.2), not containment
    verifyEqual(tc, VCh.sigma, VCx.sigma, 'RelTol', 0.05);
    % estimated redundancies agree with exact within a few * their own SE
    band = max(4 * VCh.redundSE, 0.05 * abs(VCx.redund));
    verifyTrue(tc, all(abs(VCh.redund - VCx.redund) <= band), ...
        'estimated r_g outside the sampling band around exact r_g');
end

function test_hutchinson_recovers_planted_ratio_like_demos(tc)
% categories: precision
    % Guards the CLAIM the demos make (demo_lsqvce case 5, demo_workflow §4):
    % a 2-group problem where group 2 is 3x noisier but we CLAIM equal weights;
    % both the exact and the 'backend','hutchinson' paths must recover the true
    % sigma ratio of ~3. Since L10 (WP-74) the demos run the recovery through
    % lsqvce.fit(Model="vce") and only NAME the Hutchinson backend; the
    % exact-vs-Hutchinson agreement they cite is guarded HERE, at engine tier.
    rs = RandStream('mt19937ar', 'Seed', 2026);
    u = 8;  ng = 500;  n = 2*ng;  nb = 10;
    A  = randn(rs, n, u);  xt = randn(rs, u, 1);
    g  = [ones(ng,1); 2*ones(ng,1)];
    sg = ones(n, 1);  sg(g==2) = 3;            % group 2 is 3x noisier (truth)
    p  = ones(n, 1);                            % but we CLAIM equal weights
    l  = A*xt + sg .* randn(rs, n, 1);
    e  = round(linspace(0, n, nb+1));
    src  = @(k) deal(A(e(k)+1:e(k+1), :), l(e(k)+1:e(k+1)), p(e(k)+1:e(k+1)));
    gsrc = @(k) g(e(k)+1:e(k+1));

    VCe = lsqvce.engine.varComp(u, src, nb, gsrc, 2, 'verbose', false);
    VCh = lsqvce.engine.varComp(u, src, nb, gsrc, 2, 'verbose', false, ...
             'backend', 'hutchinson', 'probes', 4000, 'seed', 2026);

    % tol: Monte-Carlo probe scatter ~ 1/sqrt(Probes) on this fixture; 4/sqrt(Z)-style bar (HANDOVER 6.1.2), not containment
    verifyEqual(tc, VCe.sigma(2)/VCe.sigma(1), 3.0, 'RelTol', 0.15, ...
        'exact backend must recover the planted 3x ratio');
    % tol: Monte-Carlo probe scatter ~ 1/sqrt(Z), Z = 4000 -> ~1.6e-2;
    % 0.05 = ~3 sigma (HANDOVER 6.1.2), not containment
    verifyEqual(tc, VCh.sigma(2)/VCh.sigma(1), VCe.sigma(2)/VCe.sigma(1), ...
        'RelTol', 0.05, ...
        'Hutchinson backend (Z=4000) must match exact ratio the demos print');
end

function test_vceHutchinson_simplified_matches(tc)
% categories: precision
    % The simplified estimator (s_g = q_g / r_g) must also track exact varComp's
    % simplified path within the band.
    [u, src, gsrc, nb] = vcSmall();
    VCx = lsqvce.engine.varComp(u, src, nb, gsrc, 2, 'verbose', false, 'method', 'simplified');
    VCh = lsqvce.engine.vceHutchinson(u, src, nb, gsrc, 2, 'verbose', false, ...
             'method', 'simplified', 'probes', 4000, 'seed', 2);
    % tol: Monte-Carlo probe scatter ~ 1/sqrt(Probes) on this fixture; 4/sqrt(Z)-style bar (HANDOVER 6.1.2), not containment
    verifyEqual(tc, VCh.sigma, VCx.sigma, 'RelTol', 0.05);
end

function test_vceHutchinson_deterministic_with_seed(tc)
% categories: contract
    % Same seed + same probe count -> identical result (reproducible).
    [u, src, gsrc, nb] = vcSmall();
    V1 = lsqvce.engine.vceHutchinson(u, src, nb, gsrc, 2, 'verbose', false, 'probes', 300, 'seed', 9);
    V2 = lsqvce.engine.vceHutchinson(u, src, nb, gsrc, 2, 'verbose', false, 'probes', 300, 'seed', 9);
    verifyEqual(tc, V1.s,      V2.s);
    verifyEqual(tc, V1.redund, V2.redund);
end

%% ============================ varComp backend dispatch (item 0 tier T2)
function test_varComp_backend_hutchinson_delegates(tc)
% categories: contract
    % 'backend','hutchinson' hands the whole job to vceHutchinson: the returned
    % struct is the estimated one (marked, with a redundSE), not the exact path.
    [u, src, gsrc, nb] = vcSmall();
    VC = lsqvce.engine.varComp(u, src, nb, gsrc, 2, 'verbose', false, ...
             'backend', 'hutchinson', 'probes', 500, 'seed', 4);
    verifyEqual(tc, VC.estimator, 'hutchinson');
    verifyTrue(tc, isfield(VC, 'redundSE'));
end

function test_varComp_backend_matches_exact(tc)
% categories: precision
    % The delegated result must agree with the exact backend within the band.
    [u, src, gsrc, nb] = vcSmall();
    VCx = lsqvce.engine.varComp(u, src, nb, gsrc, 2, 'verbose', false);                % exact
    VCh = lsqvce.engine.varComp(u, src, nb, gsrc, 2, 'verbose', false, ...
             'backend', 'hutchinson', 'probes', 4000, 'seed', 6);
    % tol: Monte-Carlo probe scatter ~ 1/sqrt(Probes) on this fixture; 4/sqrt(Z)-style bar (HANDOVER 6.1.2), not containment
    verifyEqual(tc, VCh.sigma, VCx.sigma, 'RelTol', 0.05);
end

function test_varComp_backend_hutchinson_checkD_fallsback(tc)
% categories: robustness | tags: diagnostics
    % Supplying 'checkD' must force a fallback to the EXACT backend (with a
    % warning), because the contamination guard needs the per-obs path. The
    % result is then the exact struct (no .estimator field).
    [u, src, gsrc, nb] = vcSmall();
    ws = warning('off', 'lsqvce:varComp:hutchNoCheckD');
    c  = onCleanup(@() warning(ws)); %#ok<NASGU>
    cd = struct('contamLevel', 'clean', 'fhat', 0);   % benign guard struct
    VC = lsqvce.engine.varComp(u, src, nb, gsrc, 2, 'verbose', false, ...
             'backend', 'hutchinson', 'checkD', cd);
    verifyFalse(tc, isfield(VC, 'estimator'));         % exact path, not delegated
end

%% ==================================================== report
function test_report_verdict_clean(tc)
% categories: contract | tags: diagnostics
    [A, l, p] = tc.TestData.mk(500, 5, 91);
    B = chop(A, l, p, 6);
    [~, ~, V] = lsqvce.engine.adjust(5, srcOf(B), 6, 'verbose', false);
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
    [~, ~, V] = lsqvce.engine.adjust(u, srcOf(B), 4, 'verbose', false);
    codes = cellfun(@(f) f.code, V.findings, 'UniformOutput', false);
    verifyTrue(tc, ismember('REL_UNCONTROLLED', codes), ...
        'an uncontrolled observation was not reported');
end
