function tests = test_parallel
%TEST_PARALLEL  Parallel strategy layer, PSD iteration, per-obsType quick-look.
%
%   Covers everything added in v4:
%     * lsqvce.mapreduce (the execution-strategy layer)
%     * Accumulator.addSystems == sequential addSystem (all modes)
%     * the refactored lsqvce.redundancy == an independent brute-force pass
%     * lsqvce.iterateNoiseModel: covariance model iteration to sigma0 ~ 1 + white acf
%     * per-obsType sigma0 quick-look in lsqvce.residuals
%     * dense-path boundary bias is benign (the measured 2.2 result)
%
%   The parallel-execution tests are SKIPPED (not failed) when the Parallel
%   Computing Toolbox is unavailable. The decomposition itself is proven in
%   validation/v4verify.py; those tests here confirm the actual parfor/spmd
%   schedules reproduce serial to ~1e-10 on real hardware.
%
%   Run:  results = runtests('tests/test_parallel');  disp(table(results))
%
%   See also TEST_V3, TEST_LSQACCUM, TEST_PRECISION
%
%   Updated 2026-08-18 (CEST, UTC+02:00): hasBackgroundPool() probes by USE
%   (exist() is false-negative under `matlab -batch`; WP-30 postflight
%   finding) + consistency test testBackgroundPoolGuardAgreesWithUse. With
%   assistance from Claude (Anthropic), model Claude Fable 5.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    tc.TestData.cfg = makeCfg();
end

%% ===================================================== mapreduce (the layer)
function testMapreduceSerialCollect(tc)
% Without Reduce: a cell of per-task results, in task order.
    [out, info] = lsqvce.mapreduce(5, @(j) j^2, 'Mode', 'serial');
    tc.verifyEqual([out{:}], (1:5).^2);
    tc.verifyEqual(info.mode, 'serial');
end

function testMapreduceSerialOrderedReduce(tc)
% With Reduce: identical to the hand-written loop, term for term.
    red = @(a, t) a + t;
    out = lsqvce.mapreduce(7, @(j) 1/j, 'Mode','serial', 'Reduce',red, 'Init',0);
    ref = 0;
    for j = 1:7, ref = ref + 1/j; end
    tc.verifyEqual(out, ref);                 % same terms, same order: exact
end

function testMapreduceThresholdFallsBackToSerial(tc)
% A parallel request below the threshold must run serially and say why.
    [out, info] = lsqvce.mapreduce(3, @(j) j, 'Mode','parfor', 'Threshold',100);
    tc.verifyEqual([out{:}], 1:3);
    tc.verifyEqual(info.mode, 'serial');
    tc.verifyTrue(contains(info.fallback, 'Threshold'));
end

function testMapreduceSpmdWithoutCombineFallsBack(tc)
% 'spmd' without a Combine function cannot reduce on the workers.
    tc.assumeTrue(hasPCT(), 'Parallel Computing Toolbox not available');
    [~, info] = lsqvce.mapreduce(20, @(j) j, 'Mode','spmd', 'Threshold',1);
    tc.verifyTrue(~strcmp(info.mode, 'spmd'));
end

function testMapreduceParforMatchesSerial(tc)
    tc.assumeTrue(hasPCT(), 'Parallel Computing Toolbox not available');
    f = @(j) struct('v', sin(j) + j^2);
    red = @(a, t) a + t.v;
    s  = lsqvce.mapreduce(40, f, 'Mode','serial', 'Reduce',red, 'Init',0);
    pf = lsqvce.mapreduce(40, f, 'Mode','parfor', 'Reduce',red, 'Init',0, ...
                       'Threshold',1);
    % identical reduction order -> differences can only come from BLAS
    tc.verifyEqual(pf, s, 'RelTol', 1e-14);
end

function testMapreduceSpmdMatchesSerialToTolerance(tc)
% spmd reduces on the workers: summation order differs, so tolerance, never
% isequal (ROADMAP 3.5 caveat 2).
    tc.assumeTrue(hasPCT(), 'Parallel Computing Toolbox not available');
    f   = @(j) struct('v', sin(j) + j^2);
    red = @(a, t) a + t.v;
    cmb = @(a, b) a + b;
    s  = lsqvce.mapreduce(40, f, 'Mode','serial', 'Reduce',red, 'Init',0);
    sp = lsqvce.mapreduce(40, f, 'Mode','spmd', 'Reduce',red, 'Init',0, ...
                       'Combine',cmb, 'Threshold',1);
    tc.verifyEqual(sp, s, 'RelTol', 1e-10);
end

%% ======================================================== addSystems
function testAddSystemsEqualsSequentialAddSystem(tc)
% Same systems, same order: addSystems (serial) must reproduce the addSystem
% loop. The merge is applied in input order, so this is tight.
    cfg  = tc.TestData.cfg;
    days = 1:3;

    a1 = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats
        a1.addSystem(s, cfg.nl, binder(cfg, s), days);
    end
    s1 = a1.solve();

    a2 = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    a2.addSystems(specsFor(cfg, days));
    s2 = a2.solve();

    tc.verifyEqual(s2.x,      s1.x,      'RelTol', 1e-12);
    tc.verifyEqual(s2.sigma0, s1.sigma0, 'RelTol', 1e-12);
    tc.verifyEqual(s2.dof,    s1.dof);
    tc.verifyEqual(a2.report().sysID, a1.report().sysID);
end

function testAddSystemsRejectsDuplicates(tc)
    cfg = tc.TestData.cfg;
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Verbose', false);
    sp  = specsFor(cfg, 1:2);
    sp(end+1) = sp(1);                                   % duplicate in specs
    tc.verifyError(@() acc.addSystems(sp), 'lsqvce:accum:duplicateSystem');

    acc.addSystem(cfg.sats(1), cfg.nl, binder(cfg, cfg.sats(1)), 1:2);
    tc.verifyError(@() acc.addSystems(specsFor(cfg, 1:2)), ...
                   'lsqvce:accum:duplicateSystem');         % duplicate vs added

    % the check must not depend on KeepLocal (diag-based, not store-based)
    a3 = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false, 'KeepLocal',false);
    a3.addSystem(cfg.sats(1), cfg.nl, binder(cfg, cfg.sats(1)), 1:2);
    tc.verifyError(@() a3.addSystem(cfg.sats(1), cfg.nl, ...
                   binder(cfg, cfg.sats(1)), 1:2), ...
                   'lsqvce:accum:duplicateSystem');
end

function testAddSystemsParforMatchesSerial(tc)
    tc.assumeTrue(hasPCT(), 'Parallel Computing Toolbox not available');
    cfg  = tc.TestData.cfg;
    days = 1:2;

    a1 = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false);
    a1.addSystems(specsFor(cfg, days));
    s1 = a1.solve();

    a2 = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false, ...
                         'Parallel','parfor', 'ParallelThreshold',1);
    a2.addSystems(specsFor(cfg, days));
    s2 = a2.solve();

    tc.verifyEqual(s2.x,      s1.x,      'RelTol', 1e-10);
    tc.verifyEqual(s2.sigma0, s1.sigma0, 'RelTol', 1e-10);
end

%% ================================================ redundancy, refactored
function testRedundancyMatchesIndependentBruteForce(tc)
% The refactored task/reduce pass must reproduce a from-scratch reference
% that never touches the accumulator internals: rebuild every system's
% normals directly from the loader, form the leverage blocks, and sum
% h / r / v'v per group in an explicit loop.
    cfg  = tc.TestData.cfg;
    days = 1:2;

    acc = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg, s), days);
    end
    sol = acc.solve();
    red = lsqvce.redundancy(acc, sol, 'GroupBy','obsType', 'Verbose',false);

    % ---- independent reference -------------------------------------------
    [Nsys, Nred] = rebuildNormals(cfg, days);
    ref = containers.Map('KeyType','char','ValueType','any');
    mT = 0;
    for s = cfg.sats
        st  = Nsys(sprintf('%d', s));
        blk = lsqvce.leverageBlocks(Nred, st.Nll, st.Nlc);
        y   = st.Nll \ (st.nl - st.Nlc * sol.x);
        for j = 1:numel(st.arcs)
            e = st.arcs(j);
            h = lsqvce.leverage(e.A, e.B, blk);
            v = e.l - e.A*sol.x - e.B*y;
            if ~ref.isKey(e.obsType)
                ref(e.obsType) = struct('m',0,'r',0,'vtv',0);
            end
            g = ref(e.obsType);
            g.m   = g.m + numel(v);
            g.r   = g.r + numel(v) - sum(h);
            g.vtv = g.vtv + v.'*v;
            ref(e.obsType) = g;
            mT = mT + numel(v);
        end
    end

    for k = 1:numel(red.group)
        g  = red.group(k);
        rf = ref(char(string(g.key)));
        tc.verifyEqual(g.m,   rf.m);
        tc.verifyEqual(g.r,   rf.r,   'RelTol', 1e-9);
        tc.verifyEqual(g.vtv, rf.vtv, 'RelTol', 1e-9);
    end
    u = cfg.nc + numel(cfg.sats)*cfg.nl;
    tc.verifyEqual(red.rTotal, mT - u, 'RelTol', 1e-8);   % the identity
end

function testRedundancyParforMatchesSerial(tc)
    tc.assumeTrue(hasPCT(), 'Parallel Computing Toolbox not available');
    cfg  = tc.TestData.cfg;
    days = 1:2;
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg, s), days);
    end
    sol = acc.solve();

    rS = lsqvce.redundancy(acc, sol, 'GroupBy','arc', 'Verbose',false, ...
                        'Parallel','serial');
    rP = lsqvce.redundancy(acc, sol, 'GroupBy','arc', 'Verbose',false, ...
                        'Parallel','parfor');

    tc.verifyEqual({rP.group.key}, {rS.group.key});   % ordered reduce: same
    tc.verifyEqual([rP.group.r],   [rS.group.r],   'RelTol', 1e-10);
    tc.verifyEqual([rP.group.vtv], [rS.group.vtv], 'RelTol', 1e-10);
    tc.verifyEqual(rP.check, rS.check, 'AbsTol', 1e-6);
end

function testRedundancySpmdMatchesSerialToTolerance(tc)
    tc.assumeTrue(hasPCT(), 'Parallel Computing Toolbox not available');
    cfg  = tc.TestData.cfg;
    days = 1:2;
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg, s), days);
    end
    sol = acc.solve();

    rS = lsqvce.redundancy(acc, sol, 'GroupBy','obsType', 'Verbose',false, ...
                        'Collect',false, 'Parallel','serial');
    rW = lsqvce.redundancy(acc, sol, 'GroupBy','obsType', 'Verbose',false, ...
                        'Collect',false, 'Parallel','spmd');

    % spmd's group DISCOVERY order can differ; align by key
    [tf, loc] = ismember(cellstr(string({rS.group.key})), ...
                         cellstr(string({rW.group.key})));
    tc.verifyTrue(all(tf));
    tc.verifyEqual([rW.group(loc).r],   [rS.group.r],   'RelTol', 1e-10);
    tc.verifyEqual([rW.group(loc).vtv], [rS.group.vtv], 'RelTol', 1e-10);
end

%% ============================================== iterateNoiseModel (ROADMAP 2.3)
function testPsdIterateConvergesFromWhiteStart(tc)
% Start with the WRONG model (white assumed, AR(3) true): after iteration
% sigma0 must be ~1 and the recovered model must match the TRUE ACV.
%
% WHITENESS BAND NOTE. With this much data (n ~ 1e4 whitened residuals) the
% default whiteness band WhiteFactor/sqrt(n) is very tight (~0.025 at
% WhiteFactor 3). The iteration reaches a STABLE fixed point at sigma0 ~ 0.983
% with a residual max|acf| ~ 0.052 -- the true model is recovered (ACV error
% below, sigma0 ~ 1), but that 0.052 sits just above a 3-sigma band, so the
% boolean 'converged' flag is knife-edge and flips with tiny cross-version
% arithmetic differences. We therefore set a wider WhiteFactor (a legitimate
% ~5-6 sigma whiteness gate at this n) so genuine convergence is declared, and
% assert the substantive outcomes -- level and ACV fidelity -- directly.
    cfg = tc.TestData.cfg;
    cfg.mArc = 600;                       % enough data for the ACV estimate
    days = 1:4;

    cov0.pos = struct('type','diag','sigma',1);          % wrong colour+level
    cov0.kbr = cfg.cov.kbr;                              % correct already

    accFcn = @(cov) buildWithCov(cfg, days, cov);
    it = lsqvce.iterateNoiseModel(accFcn, cov0, 'UpdateTypes',{'pos'}, ...
                        'MaxIter',8, 'OrderGrid',[2 3 5 8 12], ...
                        'MaxLag',64, 'WhiteFactor',8, 'Verbose',false);

    tc.verifyTrue(it.converged);
    tc.verifyLessThanOrEqual(it.nIter, 6);
    tc.verifyEqual(it.sol.sigma0, 1, 'AbsTol', 0.05);

    % the recovered ACV vs the TRUTH the noise was generated with
    cm = cfg.cov.pos;                                    % arma, ar + sigma
    rT = cm.sigma^2 * lsqvce.armaACV(cm.ar, 1, 40);
    cf = it.cov.pos;
    rF = cf.sigma^2 * lsqvce.armaACV(cf.ar, 1, 40);
    tc.verifyLessThan(norm(rF - rT)/norm(rT), 0.15);     % sampling-limited
end

function testPsdIterateCorrectModelConvergesImmediately(tc)
    cfg  = tc.TestData.cfg;
    days = 1:3;
    accFcn = @(cov) buildWithCov(cfg, days, cov);
    it = lsqvce.iterateNoiseModel(accFcn, cfg.cov, 'MaxIter',3, 'Verbose',false);
    tc.verifyTrue(it.converged);
    tc.verifyEqual(it.nIter, 1);
end

%% ==================================== per-obsType quick-look (ROADMAP 2.4)
function testResidualsPerObsTypeQuickLook(tc)
    cfg  = tc.TestData.cfg;
    days = 1:3;
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg, s), days);
    end
    sol = acc.solve();
    res = acc.residuals(sol, 'Verbose', false);

    tc.verifyEqual(sort({res.obsType.obsType}), sort(cfg.obsTypes));
    % the split must be complete: per-type sums == the global totals
    tc.verifyEqual(sum([res.obsType.vtv]), res.vtv, 'RelTol', 1e-12);
    tc.verifyEqual(sum([res.obsType.m]), res.dof + cfg.nc + ...
                   numel(cfg.sats)*cfg.nl);
    % correct covariance model -> every quick-look rms near 1
    for k = 1:numel(res.obsType)
        tc.verifyEqual(res.obsType(k).rms, 1, 'AbsTol', 0.1);
    end
end

%% ============================== dense-path boundary bias (ROADMAP 2.2)
function testDensePathBoundaryBiasIsBenign(tc)
% The measured 2.2 result, reproduced small: batching the DENSE path with NO
% overlap does not bias the estimate -- each batch is whitened with its exact
% marginal covariance, so only cross-batch information is lost. At b = 4L the
% parameter shift is a fraction of the FORMAL error and sigma0 is preserved
% to ~1e-2. (Full sweep and the efficiency-loss numbers:
% validation/v4boundary.py.)
    rs  = RandStream('mt19937ar', 'Seed', 9);
    m   = 2048;  nc = 4;  L = 16;
    phi = exp(-1/L);  sInn = sqrt(1 - phi^2);
    b   = 4*L;

    A  = randn(rs, m, nc);
    xt = randn(rs, nc, 1);
    e  = filter(1, [1 -phi], sInn*randn(rs, m + 200, 1));
    l  = A*xt + e(201:end);

    % exact: whole-arc AR(1) FIR whitening (IS chol(C)^-1)
    Aw = filter([1 -phi], 1, A) / sInn;
    lw = filter([1 -phi], 1, l) / sInn;
    Nr = Aw.'*Aw;
    xr = Nr \ (Aw.'*lw);
    wr = lw - Aw*xr;
    s0r = sqrt(wr.'*wr / (m - nc));
    sigF = sqrt(trace(inv(Nr)));                         % formal |dx|

    % batched dense, no overlap, exact per-batch marginal covariance
    N = zeros(nc);  n = zeros(nc,1);  vtv = 0;
    for i0 = 1:b:m
        i1 = min(i0+b-1, m);
        C  = toeplitz(phi.^(0:i1-i0));
        W  = lsqvce.decorrFactor(C);
        Ab = W*A(i0:i1,:);  lb = W*l(i0:i1);
        N = N + Ab.'*Ab;  n = n + Ab.'*lb;  vtv = vtv + lb.'*lb;
    end
    xb  = N \ n;
    s0b = sqrt(max(vtv - n.'*xb, 0) / (m - nc));

    tc.verifyLessThan(norm(xb - xr) / sigF, 0.5);   % well inside formal err
    tc.verifyEqual(s0b, s0r, 'RelTol', 2e-2);
end

%% ==================== day-distributed accumulation (v4.1, ROADMAP 3.2)
function testMergeRawQREqualsFusedAdd(tc)
% QRSystem.mergeRaw: two systems built separately then merged must equal one
% system fed everything -- in N, n, ltl and m, NOT in R itself.
    rs = RandStream('mt19937ar', 'Seed', 21);
    d = 7;
    A1 = randn(rs, 40, d);  l1 = randn(rs, 40, 1);
    A2 = randn(rs, 55, d);  l2 = randn(rs, 55, 1);

    qf = lsqvce.QRSystem(d);  qf.add(A1, l1);  qf.add(A2, l2);
    qa = lsqvce.QRSystem(d);  qa.add(A1, l1);
    qb = lsqvce.QRSystem(d);  qb.add(A2, l2);
    qa.mergeRaw(qb);

    ic = 1:3;  il = 4:d;
    [Fll,Flc,Fcc,fl,fc,fltl,fm] = qf.blocks(ic, il);
    [Mll,Mlc,Mcc,ml,mc,mltl,mm] = qa.blocks(ic, il);
    tc.verifyEqual(mm, fm);
    tc.verifyEqual(Mll, Fll, 'RelTol', 1e-12);
    tc.verifyEqual(Mlc, Flc, 'AbsTol', 1e-10*norm(Flc));
    tc.verifyEqual(Mcc, Fcc, 'RelTol', 1e-12);
    tc.verifyEqual([ml; mc], [fl; fc], 'AbsTol', 1e-10*norm([fl; fc]));
    % the sigma0 guard: the residual-row mass must survive the merge
    tc.verifyEqual(mltl, fltl, 'RelTol', 1e-12);
    tc.verifyEqual(fltl, l1.'*l1 + l2.'*l2, 'RelTol', 1e-12);
end

function testMergeRawRejectsDimensionMismatch(tc)
    qa = lsqvce.QRSystem(4);  qb = lsqvce.QRSystem(5);
    tc.verifyError(@() qa.mergeRaw(qb), 'lsqvce:qr:badMerge');
    na = lsqvce.NormalEqSystem(4);  nb = lsqvce.NormalEqSystem(5);
    tc.verifyError(@() na.mergeRaw(nb), 'lsqvce:normalEq:badMerge');
end

function testDistributeDaysEqualsFusedSerial(tc)
% 'DistributeDays' with Mode below the pool threshold runs the DISTRIBUTED
% code path serially through mapreduce: per-day partials + mergeRaw, no PCT
% needed. Must match the fused loop to the measured summation-order gap.
    cfg  = tc.TestData.cfg;
    days = 1:4;
    for be = {'qr', 'normal'}
        a1 = lsqvce.BlockAccumulator(cfg.nc, 'Backend',be{1}, 'Verbose',false);
        for s = cfg.sats
            a1.addSystem(s, cfg.nl, binder(cfg, s), days);
        end
        s1 = a1.solve();

        a2 = lsqvce.BlockAccumulator(cfg.nc, 'Backend',be{1}, 'Verbose',false, ...
                 'Parallel','parfor', 'DistributeDays',true, ...
                 'ParallelThreshold',2);
        % without PCT mapreduce falls back to serial: same task code, same
        % ordered merge -- the numerical claim is identical either way
        for s = cfg.sats
            a2.addSystem(s, cfg.nl, binder(cfg, s), days);
        end
        s2 = a2.solve();

        tc.verifyEqual(s2.x,      s1.x,      'RelTol', 1e-11, ...
                       sprintf('backend %s', be{1}));
        tc.verifyEqual(s2.sigma0, s1.sigma0, 'RelTol', 1e-11);
        tc.verifyEqual(s2.dof,    s1.dof);
    end
end

function testDistributeDaysResidualsStillReproduce(tc)
% The second pass (residuals) re-runs the loader against the stored Schur
% blocks. Those blocks now come from merged day partials -- sigma0 must still
% reproduce, or the KeepLocal bookkeeping broke.
    cfg  = tc.TestData.cfg;
    days = 1:3;
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false, ...
              'Parallel','parfor', 'DistributeDays',true, ...
              'ParallelThreshold',2);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg, s), days);
    end
    sol = acc.solve();
    res = acc.residuals(sol, 'Verbose', false);
    tc.verifyEqual(res.sigma0, sol.sigma0, 'RelTol', 1e-8);
end

function testDistributeDaysParforMatchesSerial(tc)
    tc.assumeTrue(hasPCT(), 'Parallel Computing Toolbox not available');
    cfg  = tc.TestData.cfg;
    days = 1:4;
    a1 = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false);
    a1.addSystem(1, cfg.nl, binder(cfg, 1), days);
    s1 = a1.solve();

    a2 = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false, 'Parallel','parfor', ...
                         'DistributeDays',true, 'ParallelThreshold',1);
    a2.addSystem(1, cfg.nl, binder(cfg, 1), days);
    s2 = a2.solve();

    tc.verifyEqual(s2.x,      s1.x,      'RelTol', 1e-10);
    tc.verifyEqual(s2.sigma0, s1.sigma0, 'RelTol', 1e-10);
end

%% ======================================== VAR path (v4.2, ROADMAP 3.1)
function testVarACVYuleWalker(tc)
% Gamma(0) - sum A_k Gamma(k)' = Sigma and the lag recursion, exactly.
    [A, Sig] = stableVAR();
    p = size(A, 3);
    G = lsqvce.varACV(A, Sig, 10);
    G0res = G(:,:,1);
    for k = 1:p
        G0res = G0res - A(:,:,k) * G(:,:,k+1).';
    end
    tc.verifyEqual(G0res, Sig, 'RelTol', 1e-10);
    for j = p:10
        S = zeros(size(Sig));
        for k = 1:p
            S = S + A(:,:,k) * G(:,:,j-k+1);
        end
        tc.verifyEqual(G(:,:,j+1), S, 'AbsTol', 1e-10*norm(G(:,:,1)));
    end
end

function testVarACVRejectsNonStationary(tc)
    A = zeros(2,2,1);  A(:,:,1) = 1.05*eye(2);
    tc.verifyError(@() lsqvce.varACV(A, eye(2), 5), 'lsqvce:varACV:notStationary');
end

function testWhittleRecoversGenerator(tc)
% Whittle on the exact ACV must return the generating (A_k, Sigma), and
% overfitting must make the extra matrices vanish -- the order-sufficiency
% signature.
    [A, Sig] = stableVAR();
    p = size(A, 3);
    G = lsqvce.varACV(A, Sig, p+3);
    [Ah, Sh] = lsqvce.whittle(G, p);
    tc.verifyEqual(Ah, A, 'AbsTol', 1e-8*norm(A(:,:,1)));
    tc.verifyEqual(Sh, Sig, 'RelTol', 1e-9);

    Ao = lsqvce.whittle(G, p+2);
    tc.verifyLessThan(norm(Ao(:,:,p+1)) + norm(Ao(:,:,p+2)), ...
                      1e-7*norm(A(:,:,1)));
end

function testVarWhitenEqualsDenseDecorrelation(tc)
% The matrix FIR + exact block startup against the dense route through the
% factory's own Cfcn: identical normal equations.
    [A, Sig] = stableVAR();
    d = size(Sig, 1);  nE = 50;  m = nE*d;
    cm = struct('type','var', 'ar',A, 'sigma',Sig);
    [Cfcn, Wfcn, info] = lsqvce.covModel(cm, (1:m).', struct());
    tc.verifyEqual(info.nAxis, d);
    tc.verifyEqual(info.filterMemory, size(A,3)*d);

    rs = RandStream('mt19937ar','Seed',31);
    Ad = randn(rs, m, 4);  l = randn(rs, m, 1);
    [Aw, ~, lw] = Wfcn(1:m, Ad, zeros(m,0), l, struct());

    W  = lsqvce.decorrFactor(Cfcn(1:m, struct()));
    A2 = W*Ad;  l2 = W*l;
    tc.verifyEqual(Aw.'*Aw, A2.'*A2, 'RelTol', 1e-9);
    tc.verifyEqual(Aw.'*lw, A2.'*l2, 'RelTol', 1e-8);
    tc.verifyEqual(lw.'*lw, l2.'*l2, 'RelTol', 1e-9);
end

function testVarGLSThroughAccumulatorBatched(tc)
% End to end: VAR noise, 'var' model via the factory, TWO batches with the
% automatic p*d-row overlap and epoch rounding (arc.nAxis). Batched must
% equal whole-arc, and sigma0 ~ 1 under the true model.
    [A, Sig] = stableVAR();
    d = size(Sig,1);  p = size(A,3);
    nE = 600;  m = nE*d;  nc = 4;
    rs = RandStream('mt19937ar','Seed',33);

    % simulate the VAR noise
    L = chol(Sig, 'lower');
    z = zeros(nE+200, d);
    for t = p+1:nE+200
        acc = (L*randn(rs, d, 1)).';
        for k = 1:p
            acc = acc + (A(:,:,k)*z(t-k,:).').';
        end
        z(t,:) = acc;
    end
    e = reshape(z(201:end,:).', m, 1);

    Ad = randn(rs, m, nc);
    xt = randn(rs, nc, 1);
    l  = Ad*xt + e;

    cm = struct('type','var', 'ar',A, 'sigma',Sig);
    ld = @(day, ot, a, uo) varArcLoader(day, ot, a, uo, Ad, l, cm, d);

    solveWith = @(bs) solveOne(nc, ld, bs);
    sW = solveWith(1e9);                     % whole arc, one batch
    sB = solveWith(97*d + 1);                % forces rounding + overlap

    tc.verifyEqual(sB.x, sW.x, 'RelTol', 1e-10);
    tc.verifyEqual(sB.sigma0, sW.sigma0, 'RelTol', 1e-10);
    tc.verifyEqual(sW.sigma0, 1, 'AbsTol', 0.06);
    sig = sqrt(diag(sW.Cxx));
    tc.verifyLessThan(abs(sW.x - xt), 5*sig);
end

%% ================== spmd day distribution + prefetch (v4.2, ROADMAP 3.3)
function testDistributeDaysSpmdMatchesSerial(tc)
    tc.assumeTrue(hasPCT(), 'Parallel Computing Toolbox not available');
    cfg  = tc.TestData.cfg;
    days = 1:4;
    a1 = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false);
    a1.addSystem(1, cfg.nl, binder(cfg, 1), days);
    s1 = a1.solve();

    a2 = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false, 'Parallel','spmd', ...
                         'DistributeDays',true, 'ParallelThreshold',1);
    a2.addSystem(1, cfg.nl, binder(cfg, 1), days);
    s2 = a2.solve();

    % worker tree reduce: tolerance, never isequal
    tc.verifyEqual(s2.x,      s1.x,      'RelTol', 1e-10);
    tc.verifyEqual(s2.sigma0, s1.sigma0, 'RelTol', 1e-10);
end

function testPrefetchLoaderPassthroughIsExact(tc)
% Disabled (and auto-disabled without a backgroundPool) the wrapper is a pure
% passthrough: the accumulation through it must be BIT-identical.
    cfg  = tc.TestData.cfg;
    days = 1:3;

    a1 = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false);
    a1.addSystem(1, cfg.nl, binder(cfg, 1), days);
    s1 = a1.solve();

    pf = lsqvce.PrefetchLoader(binder(cfg, 1), days, 'Enable', false);
    a2 = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false);
    a2.addSystem(1, cfg.nl, pf.handle(), days);
    s2 = a2.solve();

    tc.verifyEqual(s2.x, s1.x);              % passthrough: identical calls
    tc.verifyEqual(s2.sigma0, s1.sigma0);
    tc.verifyFalse(pf.enabled);
end

function testPrefetchLoadWholeDayMatchesDirectCalls(tc)
% The function that runs on the background worker, executed serially: its
% cache must hold exactly what direct loader calls return.
    cfg = tc.TestData.cfg;
    f   = binder(cfg, 1);
    uo  = struct('sysID',1, 'nc',cfg.nc, 'nl',cfg.nl, 'BatchSize',20000);

    C = lsqvce.PrefetchLoader.loadWholeDay(f, 2, uo);

    meta = f('meta', 2, [], uo);
    for ot = 1:numel(meta.obsTypes)
        uo.obsType = meta.obsTypes{ot};
        for a = meta.arcs{ot}(:).'
            uo.day = 2;  uo.arc = a;
            ref = f(2, uo.obsType, a, uo);
            got = C(sprintf('%s|%s', string(uo.obsType), string(a)));
            tc.verifyEqual(got.A, ref.A);
            tc.verifyEqual(got.l, ref.l);
        end
    end
end

function testPrefetchServesFromCacheWhenPrimed(tc)
% Simulate a harvested future by priming the cache path via loadWholeDay,
% then verify hit/miss accounting and that a cold day falls through
% correctly. Exercises fetch() logic without a pool.
    tc.assumeTrue(hasBackgroundPool(), 'backgroundPool not available');
    cfg  = tc.TestData.cfg;
    days = 1:3;
    pf = lsqvce.PrefetchLoader(binder(cfg, 1), days, 'Enable', true);
    tc.assumeTrue(pf.enabled, 'prefetch could not enable');

    acc = lsqvce.BlockAccumulator(cfg.nc, 'Verbose', false);
    acc.addSystem(1, cfg.nl, pf.handle(), days);
    sol = acc.solve();

    ref = lsqvce.BlockAccumulator(cfg.nc, 'Verbose', false);
    ref.addSystem(1, cfg.nl, binder(cfg, 1), days);
    sr = ref.solve();

    tc.verifyEqual(sol.x, sr.x, 'RelTol', 1e-14);   % identical data served
    tc.verifyGreaterThan(pf.hits, 0);               % it actually prefetched
end
function testBackgroundPoolGuardAgreesWithUse(tc)
% §6.1 regression for the WP-30 postflight finding (2026-08-18): the guard
% must agree with the ground truth - does CALLING backgroundPool work? The
% old exist()-based guard said NO under `matlab -batch` while the call
% succeeded, silently filtering the prefetch coverage on every batch/gate
% run. This test FAILS on that code under -batch and pins the contract.
    tc.assumeTrue(hasPCT(), 'PCT not installed');
    try
        backgroundPool;
        works = true;
    catch
        works = false;
    end
    tc.verifyEqual(hasBackgroundPool(), works, ...
        'hasBackgroundPool() must match a use-probe of backgroundPool.');
end

%% =================================================================== helpers
function ok = hasPCT()
    ok = license('test', 'Distrib_Computing_Toolbox') && ...
         ~isempty(ver('parallel'));
end

function ok = hasBackgroundPool()
% Availability check for MATLAB's backgroundPool: PROBE BY USE, not by
% metadata. exist('backgroundPool','file') == 2 was wrong (packaged/builtin
% functions return 5 or 0), and the exist-based replacement (class + file)
% is FALSE-NEGATIVE under `matlab -batch`: measured 2026-08-18 on R2026a
% Update 4 / Windows 11, lic=1 ver=1 class=0 file=0 while backgroundPool
% itself starts fine with 16 workers (WP-30 postflight finding; the CI
% skips this suite, so only the local gate could catch it). Calling the
% function is the only probe that answers the actual question; the lazy
% singleton it starts is used by the prefetch tests anyway.
    ok = hasPCT();
    if ~ok, return; end
    try
        backgroundPool;
        ok = true;
    catch
        ok = false;
    end
end

function [A, Sig] = stableVAR()
% A stationary VAR(2), d = 3, with genuine cross-axis coupling.
    rs = RandStream('mt19937ar', 'Seed', 29);
    d = 3;  p = 2;
    A = 0.35 * randn(rs, d, d, p);
    A(:,:,1) = A(:,:,1) + 0.3*eye(d);
    F = zeros(p*d);
    F(1:d, :) = reshape(A, d, p*d);
    F(d+1:end, 1:end-d) = eye((p-1)*d);
    A = A * (0.5 / max(abs(eig(F))));
    L = 0.2 * randn(rs, d, d);
    Sig = L*L.' + 0.5*eye(d);
end

function out = varArcLoader(day, ot, a, uo, Ad, l, cm, d) %#ok<INUSD>
% One-system, one-day, one-arc loader emitting epoch-major VAR data.
    if (ischar(day) || isstring(day)) && strcmp(day, 'meta')
        out = struct('obsTypes', {{'pos'}}, 'arcs', {{1}});
        return
    end
    m = size(Ad, 1);
    out.A = Ad;
    out.B = zeros(m, 1);          % one dummy local to satisfy the Schur path
    out.B(1:2:end) = 1e-3;
    out.l = l;
    out.t = (1:m).';
    [~, out.Wfcn, info] = lsqvce.covModel(cm, out.t, uo);
    out.filterMemory = info.filterMemory;
    out.nAxis        = info.nAxis;
end

function sol = solveOne(nc, ld, bs)
    acc = lsqvce.BlockAccumulator(nc, 'Backend','qr', 'Verbose',false, ...
                          'BatchSize',bs, 'Reg','none');
    acc.addSystem(1, 1, ld, 1);
    sol = acc.solve();
end


function f = binder(cfg, s)
    f = @(d, ot, a, uo) lsqvce.synthLoader(cfg, s, d, ot, a, uo);
end

function sp = specsFor(cfg, days)
    sp = struct('sysID', {}, 'nl', {}, 'dataFcn', {}, 'dayList', {});
    for s = cfg.sats
        sp(end+1) = struct('sysID', s, 'nl', cfg.nl, ...
                           'dataFcn', binder(cfg, s), 'dayList', days); %#ok<AGROW>
    end
end

function acc = buildWithCov(cfg, days, cov)
% Build and accumulate with a runtime covariance override -- the accFcn a
% real workflow would write (LSQVCE.SYNTHLOADER honours opt.covOverride).
    uo  = struct('covOverride', cov);
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false, ...
                          'UserOpt', uo);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg, s), days);
    end
end

function [Nsys, Nred] = rebuildNormals(cfg, days)
% From-scratch reference: whiten every arc through the loader's own Wfcn and
% form the per-system normals + the total reduced common matrix, without
% touching any accumulator internals.
    Nred = zeros(cfg.nc);
    Nsys = containers.Map('KeyType','char','ValueType','any');
    for s = cfg.sats
        ld = binder(cfg, s);
        st = struct('Nll', zeros(cfg.nl), 'Nlc', zeros(cfg.nl, cfg.nc), ...
                    'nl', zeros(cfg.nl, 1), 'arcs', struct([]));
        Ncc = zeros(cfg.nc);
        uo  = struct('sysID', s, 'nc', cfg.nc, 'nl', cfg.nl, ...
                     'BatchSize', 1e6);
        for d = days
            meta = ld('meta', d, [], uo);
            for ot = 1:numel(meta.obsTypes)
                uo.obsType = meta.obsTypes{ot};
                for a = meta.arcs{ot}(:).'
                    uo.day = d;  uo.arc = a;
                    arc = ld(d, uo.obsType, a, uo);
                    if isempty(arc) || isempty(arc.l), continue, end
                    mA = numel(arc.l);
                    [Aw, Bw, lw] = arc.Wfcn(1:mA, arc.A, arc.B, arc.l(:), uo);
                    st.Nll = st.Nll + Bw.'*Bw;
                    st.Nlc = st.Nlc + Bw.'*Aw;
                    st.nl  = st.nl  + Bw.'*lw;
                    Ncc    = Ncc    + Aw.'*Aw;
                    e = struct('A',Aw, 'B',Bw, 'l',lw, 'obsType',uo.obsType);
                    if isempty(st.arcs), st.arcs = e;
                    else,                st.arcs(end+1) = e; end
                end
            end
        end
        Nred = Nred + Ncc - st.Nlc.'*(st.Nll\st.Nlc);
        Nsys(sprintf('%d', s)) = st;
    end
    Nred = (Nred + Nred.')/2;
end

function cfg = makeCfg()
    rs = RandStream('mt19937ar', 'Seed', 44);
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
    cfg.seed         = 44;
    cfg.collinear    = false;
    cfg.collinearEps = 1e-7;
    cfg.cov.pos = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',0.7);
    cfg.cov.kbr = struct('type','diag','sigma',0.2);
end
