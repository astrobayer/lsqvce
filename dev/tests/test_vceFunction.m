function tests = test_vceFunction
%TEST_VCEFUNCTION  The two functions that had no dedicated test: lsqvce.engine.vce (Helmert
%   variance components, called end to end through its accFcn contract) and
%   lsqvce.engine.verdict.source (the bridge to the diagnostics layer).
%
%   test_vceMath already checks the VCE *mathematics* by hand; this checks the VCE
%   *function* -- the accFcn convention, the convergence loop, the output
%   struct -- and the bridge's block-arrow reconstruction.
%
%   Proven first in validation/v7verify.py (6/6). Run:
%     results = runtests('tests/test_vceFunction');  disp(table(results))
%
%   See also TEST_V3, TEST_V4, LSQVCE.ENGINE.VCE, LSQVCE.ENGINE.VERDICT.SOURCE,
%   LSQVCE.ENGINE.BLOCKCOFACTOR
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-129: the downdate of a local
%   block kept on its QR factor is pinned against a fresh QR of the
%   survivors. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 3: systemOf, and
%   snoop's clean re-solve keeps the block structure. With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-86 part 2: the bridge's
%   cofactor comes as blocks; BlockCofactor is pinned against the dense
%   inverse, and snoop's block downdate against its clean re-solve. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    tc.TestData.cfg = makeCfg();
end

%% ===================================================== lsqvce.engine.vce (the function)
function testVceFunctionRecoversComponents(tc)
% categories: precision
% Call lsqvce.engine.vce end to end. Two obsTypes with DIFFERENT true sigma; start
% from a deliberately wrong (all-ones) weighting; the recovered ratio must
% match the truth. Only the RATIO is identifiable (global sigma0 is free).
    cfg = tc.TestData.cfg;
    days = 1:4;

    accFcn = @(sigma) buildVceAcc(cfg, days, sigma);
    vc = lsqvce.engine.vce(accFcn, 'GroupBy','obsType', ...
                    'Groups', {'pos','kbr'}, 'Sigma0', [1 1], ...
                    'MinRedundancy', 10, 'Verbose', false);

    tc.verifyTrue(vc.converged);
    ratioEst  = vc.sigma(strcmp(vc.groups,'pos')) / vc.sigma(strcmp(vc.groups,'kbr'));
    ratioTrue = cfg.trueSigma.pos / cfg.trueSigma.kbr;
    % tol: recovery of a planted value; estimator sampling scatter on this fixture, measured stable across seeds (not an eps claim)
    tc.verifyEqual(ratioEst/ratioTrue, 1, 'AbsTol', 0.15);
end

function testVceIdentityAtConvergence(tc)
% categories: precision
% At the fixed point the whitened per-group vtv/r must be ~1 (the update is a
% no-op), and the global redundancy identity sum(r)=m-u must hold.
    cfg = tc.TestData.cfg;
    days = 1:4;
    accFcn = @(sigma) buildVceAcc(cfg, days, sigma);
    vc = lsqvce.engine.vce(accFcn, 'GroupBy','obsType', 'Groups', {'pos','kbr'}, ...
                    'Sigma0', [2 0.5], 'MinRedundancy', 10, 'Verbose', false);

    % re-whiten with the converged sigma, check vtv_g/r_g ~ 1 per group
    acc = accFcn(vc.sigma);
    sol = acc.solve();
    red = lsqvce.engine.redundancy(acc, sol, 'GroupBy','obsType', 'Verbose', false);
    for k = 1:numel(red.group)
        % tol: statistical identity on a finite sample; scatter ~ 1/sqrt(dof) of the fixture, bar sized 4/sqrt(Z)-style (HANDOVER 6.1.2)
        tc.verifyEqual(red.group(k).s0sq, 1, 'AbsTol', 0.10);
    end
    tc.verifyLessThan(red.check, 1e-6 * max(red.rTotal,1));
end

function testVceConvergesFast(tc)
% categories: contract | tags: diagnostics
    cfg = tc.TestData.cfg;
    accFcn = @(sigma) buildVceAcc(cfg, 1:3, sigma);
    vc = lsqvce.engine.vce(accFcn, 'GroupBy','obsType', 'Groups', {'pos','kbr'}, ...
                    'Sigma0', [1 1], 'MinRedundancy', 10, 'Verbose', false);
    tc.verifyTrue(vc.converged);
    tc.verifyLessThanOrEqual(vc.nIter, 5);
end

%% ============================================= lsqvce.engine.verdict.source (the bridge)
function testFromAccumSolContract(tc)
% categories: contract | tags: diagnostics
% The bridge must return a sol carrying exactly the fields lsqvce.engine.verdict.diagnose
% needs, with the FULL parameter count and the redundancy identity. The
% cofactor comes as blocks (WP-86 part 2): no u_full x u_full factor.
    cfg = tc.TestData.cfg;
    days = 1:3;
    acc = buildVceAcc(cfg, days, [1 1]);
    sol = acc.solve();

    [dsol, batchFcn, nBatch] = lsqvce.engine.verdict.source(acc, sol);

    tc.verifyTrue(all(isfield(dsol, {'x','Q','n','u','r','rank','defect'})));
    tc.verifyFalse(isfield(dsol, 'R'), 'no dense full-width factor (WP-86 part 2)');
    uFull = cfg.nc + numel(cfg.sats)*cfg.nl;
    tc.verifyEqual(dsol.u, uFull);
    tc.verifyEqual(dsol.r, dsol.n - dsol.u);
    tc.verifyClass(dsol.Q, 'lsqvce.engine.BlockCofactor');
    tc.verifyEqual(dsol.Q.u, uFull);
    tc.verifyEqual(numel(dsol.x), uFull);
    tc.verifyGreaterThan(nBatch, 0);
    tc.verifyTrue(isa(batchFcn, 'function_handle'));
end

function testFromAccumCommonMatchesReduced(tc)
% categories: precision | tags: diagnostics
% The full block-arrow x the bridge builds must reproduce the Schur-reduced
% common solution to machine precision (v7verify.py test 4).
    cfg = tc.TestData.cfg;
    days = 1:3;
    acc = buildVceAcc(cfg, days, [1 1]);
    sol = acc.solve();
    dsol = lsqvce.engine.verdict.source(acc, sol);
    % tol: recomputation of the same quadratic form; house exact-identity bound (HANDOVER 6.1.2)
    tc.verifyEqual(dsol.x(1:cfg.nc), sol.x, 'RelTol', 1e-9);
end

function testFromAccumBatchesReproduceNormals(tc)
% categories: precision, reference | tags: diagnostics
% Streaming the bridge's batches through the normal equations gives the full
% normal matrix N; the block cofactor must BE its inverse (the dense
% reference), and sum(r_i) = n - u over the replayed rows.
    cfg = tc.TestData.cfg;
    days = 1:2;
    acc = buildVceAcc(cfg, days, [1 1]);
    sol = acc.solve();
    [dsol, batchFcn, nBatch] = lsqvce.engine.verdict.source(acc, sol);

    uFull = dsol.u;
    N = zeros(uFull);  nObs = 0;  sumh = 0;
    for k = 1:nBatch
        [A, ~, p] = batchFcn(k);
        N = N + (A.*p).' * A;          % unit weight p=1 here
        nObs = nObs + size(A,1);
        sumh = sumh + sum(sum(dsol.Q.halfRows(A).^2, 1));
    end
    Qd = inv(N);
    % tol: block inverse vs a dense inverse of the same normals; ~cond*eps,
    % 1e-8 relative house headroom (HANDOVER 6.1.2)
    tc.verifyEqual(dsol.Q.block(1:uFull), Qd, 'AbsTol', 1e-8 * max(abs(Qd(:))));
    % tol: same comparison, the diagonal route (one triangular inverse)
    tc.verifyEqual(dsol.Q.diagonal(), diag(Qd), 'AbsTol', 1e-8 * max(abs(Qd(:))));
    tc.verifyEqual(nObs, dsol.n);
    % tol: redundancy partition sum(h_i) = u over all batches; fp accumulation
    % over n terms, 1e-6 leaves wide headroom
    tc.verifyEqual(sumh, uFull, 'AbsTol', 1e-6);
end

%% ============================================= lsqvce.engine.BlockCofactor (WP-86 part 2)
function testBlockCofactorMatchesDenseInverse(tc)
% categories: reference, precision | tags: diagnostics
% Every answer the diagnostics take from the cofactor equals the dense
% inverse of the block-arrow normals, prior rows included (mirror:
% dev/reviews/R7_2026-10-09/mirrors/wp86/blockcofactor.py).
    [Rc, blocks, N, A] = blockArrow(86);
    Q = lsqvce.engine.BlockCofactor(Rc, blocks);
    Qd = inv(N);
    sc = max(abs(Qd(:)));
    W = Q.halfRows(A);
    % tol: exact identities of the block-arrow inverse at ~1e-16 in the mirror;
    % 1e-10 relative is house headroom for BLAS order (HANDOVER 6.1.2)
    tc.verifyEqual(W.' * W, A * Qd * A.', 'AbsTol', 1e-10 * max(1, sc));
    % tol: as above
    tc.verifyEqual(Q.qxxTimes(A), Qd * A.', 'AbsTol', 1e-10 * max(1, sc));
    % tol: as above
    tc.verifyEqual(Q.diagonal(), diag(Qd), 'AbsTol', 1e-10 * sc);
    idx = [1 6 9 Q.u];                              % commons and three systems
    % tol: as above
    tc.verifyEqual(Q.block(idx), Qd(idx, idx), 'AbsTol', 1e-10 * sc);

    sd = sqrt(diag(Qd));  C = Qd ./ (sd * sd.');  C(1:Q.u+1:end) = 0;
    [cmax, pair] = Q.maxCorrelation();
    % tol: as above
    tc.verifyEqual(cmax, max(abs(C(:))), 'AbsTol', 1e-10);
    % tol: as above (the pair attains the maximum)
    tc.verifyEqual(abs(C(pair(1), pair(2))), cmax, 'AbsTol', 1e-10);

    k = Q.condLowerBound();
    ev = eig((N + N.')/2);
    tc.verifyLessThanOrEqual(k, max(ev) / min(ev) * (1 + 1e-10), 'a LOWER bound');
    tc.verifyGreaterThan(k, 0.5 * max(ev) / min(ev), '30 power steps converge on 14 parameters');
end

function testBlockCofactorDowndateMatchesDense(tc)
% categories: precision | tags: diagnostics
% snoop's step: removing one row (one of a system with locals, one of a
% local-free system) leaves the cofactor of N - p a a', and the closed form
% x + p v Qxx~ a is the re-solve without the row.
    [Rc, blocks, N, A, sysOf, Af] = blockArrow(87);
    Q = lsqvce.engine.BlockCofactor(Rc, blocks);
    rs = RandStream('mt19937ar', 'Seed', 87);
    lf = Af * randn(rs, Q.u, 1) + randn(rs, size(Af, 1), 1);
    x = N \ (Af.' * lf);
    keep = true(size(Af, 1), 1);
    for i = [find(sysOf == 2, 1), find(sysOf == 3, 1)]
        a = A(i, :);  p = 1.7;
        v = a * x - lf(i);
        [Q, ok] = Q.downdate(sqrt(1/p) * a, p);    % weight p on the row sqrt(1/p)*a
        tc.verifyTrue(ok);
        N = N - a.' * a;  keep(i) = false;
        Qd = inv(N);
        % tol: rank-one downdate vs a dense re-inverse; mirror 1e-16, house headroom
        tc.verifyEqual(Q.qxxTimes(A), Qd * A.', 'AbsTol', 1e-9 * max(1, max(abs(Qd(:)))));
        x = x + p * (v / sqrt(p)) * Q.qxxTimes(sqrt(1/p) * a);
        xd = N \ (Af(keep, :).' * lf(keep));
        % tol: Sherman-Morrison update vs a fresh solve; ~cond*eps
        tc.verifyEqual(x, xd, 'AbsTol', 1e-9 * max(1, max(abs(xd))));
    end
end

function testBlockCofactorDowndateOfAnIllConditionedLocalBlockIsExact(tc)
% categories: precision, robustness | tags: diagnostics
% WP-129: a local block the QR solve kept on its factor (a quadratic in raw
% seconds over 1e4 s, cond(B) = 1.3e8) is downdated ON that factor. The old
% refactor of the squared Nll took the eigen branch and lost the weak
% direction: after one rejection the leverage sum over the survivors was 6,
% not u = 7, and Qll was 100 % off (dev/validation/wp129_local_downdate_mirror.py).
    rs = RandStream('mt19937ar', 'Seed', 129);
    m = 200;  nc = 4;  t = linspace(0, 1e4, m).';
    A = randn(rs, m, nc);  B = [ones(m, 1), t, t.^2];
    l = A*randn(rs, nc, 1) + B*[1; 1e-3; 1e-7] + 1e-3*randn(rs, m, 1);
    function out = poly(d, ot, a, uo)                          %#ok<INUSD>
        if ischar(d) || isstring(d), out = struct('obsTypes', {{'x'}}, 'arcs', {{1}}); return, end
        out = struct('A', A, 'B', B, 'l', l);
    end
    acc = lsqvce.engine.BlockAccumulator(nc, 'Backend', 'qr', 'KeepLocal', true, 'Verbose', false);
    acc.addSystem(1, 3, @poly, 1);
    acc.solve();
    Q = lsqvce.engine.BlockCofactor.fromAccumulator(acc);
    out = 57;  keep = (1:m).' ~= out;
    [Q, ok] = Q.downdate([A(out, :), B(out, :)], 1);
    tc.assertTrue(ok, 'the row is not the only information on a local');
    h = Q.leverage(1, A(keep, :), B(keep, :));
    % tol: 1e-8 -- the identity holds to rounding on the factor; the refactor lost a whole unit
    tc.verifyEqual(sum(h), nc + 3, 'AbsTol', 1e-8, 'sum h over the survivors = u');
    [~, Rf] = qr([A(keep, :), B(keep, :)], 0);  Ri = Rf \ eye(nc + 3);   % the reference, on a factor
    % tol: 1e-6 relative -- two factorisations of a design at cond 1.3e8:
    % eps * cond ~ 3e-8 on these entries, with headroom
    tc.verifyEqual(Q.diagonal(), sum(Ri.^2, 2), 'RelTol', 1e-6);
end

function testBlockCofactorScaleEquivariance(tc)
% categories: metamorphic
% Scaling every observation by c scales N by c^2 and the cofactor by 1/c^2:
% the blocks must carry the scale through every answer.
    [Rc, blocks, ~, A] = blockArrow(88);
    c = 3;
    Q1 = lsqvce.engine.BlockCofactor(Rc, blocks);
    for s = 1:numel(blocks)
        blocks(s).Nll = c^2 * blocks(s).Nll;  blocks(s).Nlc = c^2 * blocks(s).Nlc;
    end
    Q2 = lsqvce.engine.BlockCofactor(c * Rc, blocks);
    % tol: exact scaling up to rounding of the scaled factors
    tc.verifyEqual(Q2.diagonal(), Q1.diagonal() / c^2, 'RelTol', 1e-12);
    % tol: as above
    tc.verifyEqual(Q2.qxxTimes(A), Q1.qxxTimes(A) / c^2, 'AbsTol', 1e-12 * max(1, max(max(abs(Q1.qxxTimes(A))))));
    % tol: correlation is scale-free
    tc.verifyEqual(Q2.maxCorrelation(), Q1.maxCorrelation(), 'AbsTol', 1e-12);
end

function testBlockCofactorSingularLocalBlock(tc)
% categories: robustness
% A local block the accumulator could only pre-eliminate by its pinv
% fallback (two equal local columns): the cofactor uses the SAME
% pseudo-inverse -- finite, and the hat matrix is A*pinv(N)*A'.
    [~, blocks, ~, A, sysOf, Af] = blockArrow(89);
    c1 = blocks(1).cols;
    A(:, c1(2)) = A(:, c1(1));  Af(:, c1(2)) = Af(:, c1(1));
    nc = c1(1) - 1;
    Kr = Af(size(A, 1)+1:end, 1:nc);                 % the prior rows
    Nred = Kr.' * Kr;
    for t = 1:numel(blocks)
        r = sysOf == t;  ct = blocks(t).cols;
        Ac = A(r, 1:nc);  B = A(r, ct);
        blocks(t).Nll = B.' * B;  blocks(t).Nlc = B.' * Ac;
        Nred = Nred + Ac.' * Ac - blocks(t).Nlc.' * pinv(blocks(t).Nll) * blocks(t).Nlc;
    end
    Q = lsqvce.engine.BlockCofactor(chol((Nred + Nred.')/2), blocks);
    W = Q.halfRows(A);
    tc.verifyTrue(all(isfinite(W(:))), 'a singular local block gives finite answers');
    % tol: the pseudo-inverse is formed from an eigen-decomposition of a
    % rank-deficient block; 1e-8 relative headroom
    tc.verifyEqual(W.' * W, A * pinv(Af.' * Af) * A.', 'AbsTol', 1e-8 * max(1, max(abs(A(:)))^2));
end

function testBlockCofactorRowsAreIndependent(tc)
% categories: vectorisation
% A batch is answered column by column: the batched call equals the
% row-at-a-time calls (the diagnostics' batches are whole arcs).
    [Rc, blocks, ~, A] = blockArrow(90);
    Q = lsqvce.engine.BlockCofactor(Rc, blocks);
    W = Q.halfRows(A);  Y = Q.qxxTimes(A);
    for i = [1 7 23 size(A, 1)]
        % tol: the same arithmetic on one column; BLAS blocking may reorder sums
        tc.verifyEqual(W(:, i), Q.halfRows(A(i, :)), 'AbsTol', 1e-13 * max(1, max(abs(W(:)))));
        % tol: as above
        tc.verifyEqual(Y(:, i), Q.qxxTimes(A(i, :)), 'AbsTol', 1e-13 * max(1, max(abs(Y(:)))));
    end
end

function testBlockCofactorContract(tc)
% categories: contract
    [Rc, blocks, ~, A] = blockArrow(91);
    bad = blocks;  bad(1).Nlc = bad(1).Nlc(:, 1:end-1);
    tc.verifyError(@() lsqvce.engine.BlockCofactor(Rc, bad), 'lsqvce:BlockCofactor:badBlocks');
    Q = lsqvce.engine.BlockCofactor(Rc, blocks);
    two = A(1, :);  two(blocks(2).cols) = 1;        % a row on the locals of systems 1 and 2
    tc.verifyError(@() Q.downdate(two, 1), 'lsqvce:BlockCofactor:multiSystemRow');
    tc.verifyError(@() Q.systemOf(two), 'lsqvce:BlockCofactor:multiSystemRow');
    tc.verifyEqual(Q.systemOf(A(1, :)), 1, 'a row of system 1');
    tc.verifyEqual(Q.systemOf([A(1, 1:5), zeros(1, Q.u - 5)]), 0, 'a commons-only row');
    flat = lsqvce.engine.BlockCofactor(Rc, []);      % local-free: the commons alone
    tc.verifyEqual(flat.u, size(Rc, 1));
end

function testDiagnoseReadsTheDataOnce(tc)
% categories: metamorphic, contract | tags: diagnostics
% R7 WP-121: diagnose's robust check (the MAD scale, the contamination gauge)
% is formed in the main pass for every test but 'mad' -- each batch is read
% ONCE (before: twice, with the O(nc^2) halfRows per row computed twice). The
% scale is the same either way: 'mad' still sweeps first, from the same raw
% scores, so the two routes agree bitwise.
    cfg = tc.TestData.cfg;
    acc = buildVceAcc(cfg, 1:2, [1 2.5]);
    sol = acc.solve();
    [dsol, bf, nB] = lsqvce.engine.verdict.source(acc, sol);
    cnt = containers.Map({'n'}, {0});
    D1 = lsqvce.engine.verdict.diagnose(dsol, @(k) countedBatch(bf, k, cnt), 'nBatch', nB);
    tc.verifyEqual(cnt('n'), nB, 'one read per batch');
    tc.assertTrue(D1.madAvailable, 'the robust check ran');
    cnt('n') = 0;
    D2 = lsqvce.engine.verdict.diagnose(dsol, @(k) countedBatch(bf, k, cnt), 'nBatch', nB, 'Test', "mad");
    tc.verifyEqual(cnt('n'), 2 * nB, '''mad'' needs its scale first: two reads');
    % bitwise: the same raw scores, the same median
    tc.verifyEqual(D1.robustScale, D2.robustScale, 'the fused scale is the swept one');
    tc.verifyEqual(D1.fhat, D2.fhat, 'the fused gauge is the swept one');
end

function testSnoopBlockDowndateMatchesCleanResolve(tc)
% categories: precision | tags: diagnostics
% snoop on a fit WITH locals: its downdated estimate (block cofactor + closed
% form, Reaccumulate=false) is the clean re-solve over the survivors
% (Reaccumulate=true) -- the same rejections, the same x; and the clean
% re-solve is block-arrow again (part 3), not all-commons.
    muteWarn(tc, 'lsqvce:snoop:maxReject');   % provoked on purpose here (WP-47, gate rule 7)
    cfg = tc.TestData.cfg;
    acc = buildVceAcc(cfg, 1:2, [1 2.5]);
    sol = acc.solve();
    [dsol, bf, nB] = lsqvce.engine.verdict.source(acc, sol);
    D = lsqvce.engine.verdict.diagnose(dsol, bf, 'nBatch', nB);
    args = {'Verbose', false, 'Alpha', 0.05, 'Bonferroni', false, 'MaxReject', 3};
    [solD, ~, R1] = lsqvce.engine.verdict.snoop(dsol, D, bf, args{:}, 'Reaccumulate', false);
    [solR, ~, R2] = lsqvce.engine.verdict.snoop(dsol, D, bf, args{:}, 'Reaccumulate', true);
    tc.assertEqual(R1.idx, R2.idx, 'the same rejections either way');
    % WP-86 part 3: the clean re-solve rebuilds the block-arrow system (commons
    % + each system's locals, the original layout), not an all-commons one
    tc.verifyEqual(solR.Q.nc, cfg.nc, 'the clean re-solve keeps the commons / locals split');
    tc.verifyEqual(solR.Q.u, dsol.u);
    tc.assertNotEmpty(R1.idx, 'the fixture rejects something at alpha 0.05');
    % tol: three rank-one downdates + closed-form updates vs a fresh QR solve;
    % ~cond*eps per step, 1e-8 relative headroom
    tc.verifyEqual(solD.x, solR.x, 'AbsTol', 1e-8 * max(abs(solR.x)));
end

function testFromAccumRequiresKeepLocal(tc)
% categories: contract | tags: diagnostics
% Without KeepLocal the bridge cannot back-substitute the locals -> it must
% error, not silently return a wrong full x.
    cfg = tc.TestData.cfg;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose',false, 'KeepLocal',false);
    acc.addSystem(cfg.sats(1), cfg.nl, binder(cfg, cfg.sats(1)), 1:2);
    sol = acc.solve();
    tc.verifyError(@() lsqvce.engine.verdict.source(acc, sol), ?MException);
end

%% =================================================================== helpers
function [Rc, blocks, N, A, sysOf, Af] = blockArrow(seed)
%BLOCKARROW  A random block-arrow problem (5 commons; systems with 2, 3, 0 and
%   4 locals) with a 2-row common prior: the Schur blocks the accumulator
%   keeps, the commons factor, and the dense normals N of the rows Af = [A;
%   prior] as the reference.
    rs = RandStream('mt19937ar', 'Seed', seed);
    nc = 5;  nls = [2 3 0 4];  u = nc + sum(nls);
    offs = nc + [0 cumsum(nls(1:end-1))];
    A = zeros(0, u);  sysOf = zeros(0, 1);
    for s = 1:numel(nls)
        m = 12 + 3 * nls(s);
        As = zeros(m, u);
        As(:, 1:nc) = randn(rs, m, nc);
        As(:, offs(s) + (1:nls(s))) = randn(rs, m, nls(s));
        A = [A; As];  sysOf = [sysOf; repmat(s, m, 1)];   %#ok<AGROW> fixture
    end
    K = randn(rs, 2, nc) / 0.5;
    Af = [A; K, zeros(2, u - nc)];
    N = Af.' * Af;
    Nred = K.' * K;
    blocks = struct('cols', {}, 'Nll', {}, 'Nlc', {});
    for s = 1:numel(nls)
        r = sysOf == s;  cols = offs(s) + (1:nls(s));
        Ac = A(r, 1:nc);  B = A(r, cols);
        Nll = B.' * B;  Nlc = B.' * Ac;
        Nred = Nred + Ac.' * Ac;
        if nls(s) > 0, Nred = Nred - Nlc.' * (Nll \ Nlc); end
        blocks(end+1) = struct('cols', cols, 'Nll', Nll, 'Nlc', Nlc);   %#ok<AGROW> fixture
    end
    Rc = chol((Nred + Nred.') / 2);
end

function acc = buildVceAcc(cfg, days, sigma)
% accFcn for lsqvce.engine.vce: build+accumulate with the given per-group sigma, via
% the loader's covOverride. sigma is indexed [pos kbr] (the Groups order).
    cov.pos = struct('type','diag', 'sigma', sigma(1));
    cov.kbr = struct('type','diag', 'sigma', sigma(2));
    uo  = struct('covOverride', cov);
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false, ...
                                  'KeepLocal',true, 'UserOpt', uo);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg, s), days);
    end
end

function f = binder(cfg, s)
    f = @(d, ot, a, uo) lsqvce.engine.synthLoader(cfg, s, d, ot, a, uo);
end

function cfg = makeCfg()
    rs = RandStream('mt19937ar', 'Seed', 71);
    cfg.nc   = 8;
    cfg.nl   = 4;
    cfg.sats = [1 2];
    cfg.x_c  = randn(rs, cfg.nc, 1);
    cfg.y    = containers.Map('KeyType','double','ValueType','any');
    for s = cfg.sats, cfg.y(s) = randn(rs, cfg.nl, 1); end
    cfg.obsTypes     = {'pos','kbr'};
    cfg.arcsPerDay   = 3;
    cfg.mArc         = 300;
    cfg.dt           = 5;
    cfg.seed         = 71;
    cfg.collinear    = false;
    cfg.collinearEps = 1e-7;
    % TRUE per-group noise: white, different sigma. The synthLoader draws noise
    % from cfg.cov.(ot); make those the truth, and let covOverride re-weight.
    cfg.trueSigma.pos = 1.0;
    cfg.trueSigma.kbr = 2.5;
    cfg.cov.pos = struct('type','diag', 'sigma', cfg.trueSigma.pos);
    cfg.cov.kbr = struct('type','diag', 'sigma', cfg.trueSigma.kbr);
end

function [A, l, p] = countedBatch(bf, k, cnt)
%COUNTEDBATCH  bf(k), counting the reads in the handle cnt (WP-121).
    cnt('n') = cnt('n') + 1;
    [A, l, p] = bf(k);
end
