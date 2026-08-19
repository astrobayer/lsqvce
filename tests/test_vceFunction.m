function tests = test_vceFunction
%TEST_VCEFUNCTION  The two functions that had no dedicated test: lsqvce.vce (Helmert
%   variance components, called end to end through its accFcn contract) and
%   lsqvce.diagnoseAccum (the bridge to the diagnostics layer).
%
%   test_vceMath already checks the VCE *mathematics* by hand; this checks the VCE
%   *function* -- the accFcn convention, the convergence loop, the output
%   struct -- and the bridge's block-arrow reconstruction.
%
%   Proven first in validation/v7verify.py (6/6). Run:
%     results = runtests('tests/test_vceFunction');  disp(table(results))
%
%   See also TEST_V3, TEST_V4, LSQVCE.VCE, LSQVCE.DIAGNOSEACCUM
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    tc.TestData.cfg = makeCfg();
end

%% ===================================================== lsqvce.vce (the function)
function testVceFunctionRecoversComponents(tc)
% categories: precision
% Call lsqvce.vce end to end. Two obsTypes with DIFFERENT true sigma; start
% from a deliberately wrong (all-ones) weighting; the recovered ratio must
% match the truth. Only the RATIO is identifiable (global sigma0 is free).
    cfg = tc.TestData.cfg;
    days = 1:4;

    accFcn = @(sigma) buildVceAcc(cfg, days, sigma);
    vc = lsqvce.vce(accFcn, 'GroupBy','obsType', ...
                    'Groups', {'pos','kbr'}, 'Sigma0', [1 1], ...
                    'MinRedundancy', 10, 'Verbose', false);

    tc.verifyTrue(vc.converged);
    ratioEst  = vc.sigma(strcmp(vc.groups,'pos')) / vc.sigma(strcmp(vc.groups,'kbr'));
    ratioTrue = cfg.trueSigma.pos / cfg.trueSigma.kbr;
    tc.verifyEqual(ratioEst/ratioTrue, 1, 'AbsTol', 0.15);
end

function testVceIdentityAtConvergence(tc)
% categories: precision
% At the fixed point the whitened per-group vtv/r must be ~1 (the update is a
% no-op), and the global redundancy identity sum(r)=m-u must hold.
    cfg = tc.TestData.cfg;
    days = 1:4;
    accFcn = @(sigma) buildVceAcc(cfg, days, sigma);
    vc = lsqvce.vce(accFcn, 'GroupBy','obsType', 'Groups', {'pos','kbr'}, ...
                    'Sigma0', [2 0.5], 'MinRedundancy', 10, 'Verbose', false);

    % re-whiten with the converged sigma, check vtv_g/r_g ~ 1 per group
    acc = accFcn(vc.sigma);
    sol = acc.solve();
    red = lsqvce.redundancy(acc, sol, 'GroupBy','obsType', 'Verbose', false);
    for k = 1:numel(red.group)
        tc.verifyEqual(red.group(k).s0sq, 1, 'AbsTol', 0.10);
    end
    tc.verifyLessThan(red.check, 1e-6 * max(red.rTotal,1));
end

function testVceConvergesFast(tc)
% categories: contract | tags: diagnostics
    cfg = tc.TestData.cfg;
    accFcn = @(sigma) buildVceAcc(cfg, 1:3, sigma);
    vc = lsqvce.vce(accFcn, 'GroupBy','obsType', 'Groups', {'pos','kbr'}, ...
                    'Sigma0', [1 1], 'MinRedundancy', 10, 'Verbose', false);
    tc.verifyTrue(vc.converged);
    tc.verifyLessThanOrEqual(vc.nIter, 5);
end

%% ============================================= lsqvce.diagnoseAccum (the bridge)
function testFromAccumSolContract(tc)
% categories: contract | tags: diagnostics
% The bridge must return a sol carrying exactly the fields lsqvce.diagnose
% needs, with the FULL parameter count and the redundancy identity.
    cfg = tc.TestData.cfg;
    days = 1:3;
    acc = buildVceAcc(cfg, days, [1 1]);
    sol = acc.solve();

    [dsol, batchFcn, nBatch] = lsqvce.diagnoseAccum(acc, sol);

    tc.verifyTrue(all(isfield(dsol, {'x','R','n','u','r'})));
    uFull = cfg.nc + numel(cfg.sats)*cfg.nl;
    tc.verifyEqual(dsol.u, uFull);
    tc.verifyEqual(dsol.r, dsol.n - dsol.u);
    tc.verifyEqual(size(dsol.R,2), uFull);
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
    dsol = lsqvce.diagnoseAccum(acc, sol);
    tc.verifyEqual(dsol.x(1:cfg.nc), sol.x, 'RelTol', 1e-9);
end

function testFromAccumBatchesReproduceNormals(tc)
% categories: precision | tags: diagnostics
% Streaming the bridge's batches through the normal equations must rebuild the
% same full-width R (up to sign): N = R'R and sum(r_i) = n - u.
    cfg = tc.TestData.cfg;
    days = 1:2;
    acc = buildVceAcc(cfg, days, [1 1]);
    sol = acc.solve();
    [dsol, batchFcn, nBatch] = lsqvce.diagnoseAccum(acc, sol);

    uFull = dsol.u;
    N = zeros(uFull);  nObs = 0;
    for k = 1:nBatch
        [A, ~, p] = batchFcn(k);
        N = N + (A.*p).' * A;          % unit weight p=1 here
        nObs = nObs + size(A,1);
    end
    tc.verifyEqual(N, dsol.R.'*dsol.R, 'RelTol', 1e-8, 'AbsTol', 1e-10);
    tc.verifyEqual(nObs, dsol.n);

    % redundancy from the streamed R must satisfy the identity
    sumr = 0;
    for k = 1:nBatch
        A = batchFcn(k);
        W = dsol.R.' \ A.';
        sumr = sumr + sum(1 - sum(W.^2, 1));
    end
    tc.verifyEqual(sumr, dsol.n - dsol.u, 'AbsTol', 1e-6);
end

function testFromAccumRequiresKeepLocal(tc)
% categories: contract | tags: diagnostics
% Without KeepLocal the bridge cannot back-substitute the locals -> it must
% error, not silently return a wrong full x.
    cfg = tc.TestData.cfg;
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false, 'KeepLocal',false);
    acc.addSystem(cfg.sats(1), cfg.nl, binder(cfg, cfg.sats(1)), 1:2);
    sol = acc.solve();
    tc.verifyError(@() lsqvce.diagnoseAccum(acc, sol), ?MException);
end

%% =================================================================== helpers
function acc = buildVceAcc(cfg, days, sigma)
% accFcn for lsqvce.vce: build+accumulate with the given per-group sigma, via
% the loader's covOverride. sigma is indexed [pos kbr] (the Groups order).
    cov.pos = struct('type','diag', 'sigma', sigma(1));
    cov.kbr = struct('type','diag', 'sigma', sigma(2));
    uo  = struct('covOverride', cov);
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false, ...
                                  'KeepLocal',true, 'UserOpt', uo);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg, s), days);
    end
end

function f = binder(cfg, s)
    f = @(d, ot, a, uo) lsqvce.synthLoader(cfg, s, d, ot, a, uo);
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
