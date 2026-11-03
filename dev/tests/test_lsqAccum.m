function tests = test_lsqAccum
%TEST_LSQACCUM  Unit tests for the lsq accumulation toolbox.
%   Run:  results = runtests('tests/test_lsqAccum');  disp(table(results))
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-133: tsvd / ridge dof, a coupled
%   LocalBlockSize, a prior without redundancy, a rank-0 local block, a
%   sparse B through snoop (review 3 findings 3, 6 and the lows). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-131: residuals 'Whitened',
%   false keeps the whitened scalars. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
    tests = functiontests(localfunctions);
end

%% ---------------------------------------------------------------- fixtures
function setupOnce(tc)
    tc.TestData.cfg = baseCfg();
end

function cfg = baseCfg()
    rs = RandStream('mt19937ar','Seed',42);
    cfg.nc  = 6;
    cfg.nl  = 4;
    cfg.sats = [1 2 3];
    cfg.x_c = randn(rs, cfg.nc, 1);
    cfg.y   = containers.Map('KeyType','double','ValueType','any');
    for s = cfg.sats, cfg.y(s) = randn(rs, cfg.nl, 1); end
    cfg.obsTypes   = {'pos','kbr'};
    cfg.arcsPerDay = 2;
    cfg.mArc       = 400;
    cfg.dt         = 5;
    cfg.seed       = 42;
    cfg.collinear  = false;
    cfg.collinearEps = 1e-7;
    cfg.cov.pos = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',0.7);
    cfg.cov.kbr = struct('type','diag','sigma',0.2);
end

function f = binder(cfg, s)
    f = @(d,ot,a,uo) lsqvce.engine.synthLoader(cfg, s, d, ot, a, uo);
end

function sol = runIt(cfg, days, varargin)
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, varargin{:});
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg,s), days);
    end
    sol = acc.solve();
end

%% ---------------------------------------------------- signal-processing core
function testPSD2ACVMatchesAnalyticAR(tc)
% categories: precision, reference
% PSD of a known AR(3) -> ACV must match the analytic ARMA ACV.
    a = [1 -1.6 0.9 -0.25];  s = 0.7;  dt = 5;
    r_ana = s^2 * lsqvce.kernel.armaACV(a, 1, 64);
    f = linspace(1e-9, 0.5/dt, 8192).';
    H = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * a(:));
    S = s^2 * dt * 2 * abs(H).^2;
    r_num = lsqvce.kernel.psd2acv(f, S, dt, 64);
    tc.verifyLessThan(norm(r_num-r_ana)/norm(r_ana), 2e-3);
end

function testLevinsonRecoversAR(tc)
% categories: precision
% Levinson-Durbin on the exact ACV must return the generating coefficients.
    a = [1 -1.6 0.9 -0.25];  s = 0.7;
    r = s^2 * lsqvce.kernel.armaACV(a, 1, 8);
    [ah, s2] = lsqvce.kernel.levinsonDurbin(r, 3);
    % tol: Levinson on the EXACT generator ACV; recursion error growth over p steps ~ 1e3*eps, bound adds headroom
    tc.verifyEqual(ah, a, 'AbsTol', 1e-8);
    tc.verifyEqual(sqrt(s2), s, 'RelTol', 1e-8);
end

function testFilterEqualsCholeskyWhitening(tc)
% categories: precision
    seedFor(tc, 73);   % R7-45: drew from the global stream unseeded (suite-order dependent)
% THE key invariant: the FIR whitening filter and the dense Cholesky
% decorrelation give the SAME Y'*Y (they differ by an orthogonal rotation).
    a = [1 -1.6 0.9 -0.25];  s = 0.7;  n = 300;
    r = s^2 * lsqvce.kernel.armaACV(a, 1, n);
    X = randn(n, 5);
    W  = lsqvce.kernel.decorrFactor(toeplitz(r));
    Y1 = W*X;
    z  = zeros(n,0);
    [~,~,Y2] = lsqvce.kernel.armaWhiten(a, 1, s, r, z, z, X, true);
    % tol: dense reconstruction of the same algebra; eps-level accumulation, house bounds (HANDOVER 6.1.2)
    tc.verifyEqual(Y2.'*Y2, Y1.'*Y1, 'RelTol', 1e-9);
end

function testARMAACVAgainstBruteForce(tc)
% categories: precision, reference
% ARMA ACV via the FFT/Wold expansion vs a direct long simulation moment.
    a = [1 -0.9];  b = [1 0.5];
    r = lsqvce.kernel.armaACV(a, b, 5);
    % analytic ARMA(1,1) with unit innovation variance:
    phi = 0.9; th = 0.5;
    r0 = (1 + 2*phi*th + th^2)/(1-phi^2);
    r1 = phi*r0 + th;
    % tol: closed-form literature value; kernel's own accuracy claim (help text) bounds the error, tol adds headroom
    tc.verifyEqual(r(1), r0, 'RelTol', 1e-6);
    tc.verifyEqual(r(2), r1, 'RelTol', 1e-6);
end

%% ---------------------------------------------------------- core estimation
function testRecoversTruthWhiteNoise(tc)
% categories: precision
    cfg = tc.TestData.cfg;
    cfg.cov.pos = struct('type','diag','sigma',0.05);
    cfg.cov.kbr = struct('type','diag','sigma',0.05);
    sol = runIt(cfg, 1:3, 'Backend','qr', 'Verbose',false);
    % tol: STATISTICAL, not numerical (R7-31: the comment here had been copied
    % from a backend-agreement bound) -- x_c recovered from sigma-0.05 noise; 5e-3
    % is an estimation-error band whose margin in formal sigma is not recorded
    tc.verifyEqual(sol.x, cfg.x_c, 'AbsTol', 5e-3);
end

function testSigma0IsUnbiased(tc)
% categories: precision
% With a correctly specified covariance, sigma0 must be ~1.
    cfg = tc.TestData.cfg;
    sol = runIt(cfg, 1:5, 'Backend','qr', 'Verbose',false);
    % tol: statistical identity on a finite sample; scatter ~ 1/sqrt(dof) of the fixture, bar sized 4/sqrt(Z)-style (HANDOVER 6.1.2)
    tc.verifyEqual(sol.sigma0, 1, 'AbsTol', 0.06);
end

function testBackendsAgree(tc)
% categories: precision
% Well-conditioned problem: normal and QR must give the same answer.
    cfg = tc.TestData.cfg;
    s1 = runIt(cfg, 1:3, 'Backend','normal', 'Verbose',false);
    s2 = runIt(cfg, 1:3, 'Backend','qr',     'Verbose',false);
    % tol: two backends on one well-conditioned problem; difference ~ cond*eps, bound with release-drift headroom
    tc.verifyEqual(s2.x, s1.x, 'RelTol', 1e-7);
    tc.verifyEqual(s2.sigma0, s1.sigma0, 'RelTol', 1e-6);
end

function testPreeliminationEqualsJointSolve(tc)
% categories: precision
% Pre-elimination must be algebraically exact: the reduced solution equals
% the monolithic joint solve for [x_c; y_1; y_2; y_3].
    cfg = tc.TestData.cfg;
    days = 1:2;
    sol  = runIt(cfg, days, 'Backend','qr', 'Verbose',false);

    % assemble the full whitened system explicitly
    [Af, lf] = assembleFull(cfg, days);
    xfull = Af \ lf;
    % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
    tc.verifyEqual(sol.x, xfull(1:cfg.nc), 'RelTol', 1e-8);
end

function testBackSubstitutionRecoversLocals(tc)
% categories: precision
    cfg = tc.TestData.cfg;
    days = 1:2;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), days); end
    sol = acc.solve();

    [Af, lf] = assembleFull(cfg, days);
    xfull = Af \ lf;
    for k = 1:numel(cfg.sats)
        y = acc.backSubstitute(cfg.sats(k), sol.x);
        yref = xfull(cfg.nc + (k-1)*cfg.nl + (1:cfg.nl));
        % tol: independent dense solve of the same system; difference ~ cond*eps, house bounds add headroom (HANDOVER 6.1.2)
        tc.verifyEqual(y, yref, 'RelTol', 1e-7);
    end
end

%% -------------------------------------------------------------- batching
function testBatchOverlapExactness(tc)
% categories: vectorisation
% Batching with overlap >= filter memory must reproduce the unbatched result.
% With zero overlap the filter transient leaks in and the solution degrades.
    cfg = tc.TestData.cfg;
    days = 1:2;
    ref = runIt(cfg, days, 'Backend','qr','BatchSize',1e9, 'Verbose',false);

    ok  = runIt(cfg, days, 'Backend','qr','BatchSize',100, ...
                'BatchOverlap',3, 'Verbose',false);
    bad = runIt(cfg, days, 'Backend','qr','BatchSize',100, ...
                'BatchOverlap',0, 'Verbose',false);

    eOK  = norm(ok.x  - ref.x)/norm(ref.x);
    eBAD = norm(bad.x - ref.x)/norm(ref.x);

    tc.verifyLessThan(eOK, 1e-6);
    tc.verifyGreaterThan(eBAD, 10*eOK);     % zero overlap is demonstrably worse
end

function testAutoOverlapFromFilterMemory(tc)
% categories: precision, vectorisation
% With BatchOverlap unset, the loader's filterMemory must be picked up.
    cfg = tc.TestData.cfg;
    ref = runIt(cfg, 1:2, 'Backend','qr','BatchSize',1e9, 'Verbose',false);
    aut = runIt(cfg, 1:2, 'Backend','qr','BatchSize',100, 'Verbose',false);
    % tol: two backends on one well-conditioned problem; difference ~ cond*eps, bound with release-drift headroom
    tc.verifyEqual(aut.x, ref.x, 'RelTol', 1e-6);
end

function testBatchSizeInvariance(tc)
% categories: precision, robustness, vectorisation
    cfg = tc.TestData.cfg;
    x = [];
    for bs = [50 137 400 1e9]
        s = runIt(cfg, 1:2, 'Backend','qr', 'BatchSize',bs, ...
                  'BatchOverlap',3, 'Verbose',false);
        % tol: two backends on one well-conditioned problem; difference ~ cond*eps, bound with release-drift headroom
        if isempty(x), x = s.x; else, tc.verifyEqual(s.x, x, 'RelTol',1e-6); end
    end
end

function testArmaWithMaIsBatchInvariant(tc)
% categories: precision, vectorisation
% WP-125: an MA part makes the whitener A(z)/B(z) an IIR filter. The overlap
% used to count the AR order alone, so the restart transient of 1/B(z) leaked
% into every batch after the first. With the overlap taken from the decay of
% 1/B(z) (covModel filterMemory) small batches reproduce the whole arc.
    cfg = tc.TestData.cfg;
    cfg.cov.pos = struct('type','arma','ar',[1 -0.5],'ma',[1 0.4],'sigma',0.7);
    ref = runIt(cfg, 1:2, 'Backend','qr','BatchSize',1e9, 'Verbose',false);
    bat = runIt(cfg, 1:2, 'Backend','qr','BatchSize',120, 'Verbose',false);
    % tol: the same exact whitening in batches vs whole; summation order only, 1e-9 relative
    tc.verifyEqual(bat.x, ref.x, 'RelTol', 1e-9);
end

function testDenseLawIsBatchInvariant(tc)
% categories: precision, vectorisation
% WP-125 (user decision 2026-10-12): a dense law cannot be joined across
% batches by an overlap; it is now factorised over the whole arc, so the
% BatchSize no longer changes the answer (before: one Cholesky per batch,
% the correlation between batches dropped).
    cfg = tc.TestData.cfg;
    a = [1 -1.6 0.9 -0.25];  s = 0.7;  dt = cfg.dt;
    f = linspace(1e-9, 0.5/dt, 4096).';
    H = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * a(:));
    cfg.cov.pos = struct('type','psd', 'f',f, 'S', s^2*dt*2*abs(H).^2);
    oTP = struct('psdMode','toeplitz');
    ref = runIt(cfg, 1:2, 'Backend','qr','BatchSize',1e9, 'UserOpt',oTP, 'Verbose',false);
    bat = runIt(cfg, 1:2, 'Backend','qr','BatchSize',50,  'UserOpt',oTP, 'Verbose',false);
    % tol: identical whole-arc factorisation; summation order only, 1e-9 relative
    tc.verifyEqual(bat.x, ref.x, 'RelTol', 1e-9);
end

%% ------------------------------------------------------ PSD -> AR pathway
function testPSDModeMatchesToeplitz(tc)
% categories: precision, vectorisation
% The AR reduction of a PSD must give the same normals as the dense Toeplitz.
    cfg = tc.TestData.cfg;
    a = [1 -1.6 0.9 -0.25];  s = 0.7;  dt = cfg.dt;
    f = linspace(1e-9, 0.5/dt, 4096).';
    H = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * a(:));
    cfg.cov.pos = struct('type','psd', 'f',f, 'S', s^2*dt*2*abs(H).^2);
    cfg.cov.kbr = struct('type','diag','sigma',0.2);

    oAR = struct('psdMode','ar',      'arOrder',20);
    oTP = struct('psdMode','toeplitz');

    sAR = runIt(cfg, 1:2, 'Backend','qr','BatchSize',400,'BatchOverlap',20, ...
                'UserOpt',oAR, 'Verbose',false);
    sTP = runIt(cfg, 1:2, 'Backend','qr','BatchSize',400,'BatchOverlap',20, ...
                'UserOpt',oTP, 'Verbose',false);

    % tol: MODEL APPROXIMATION, not cond*eps (R7-31) -- an AR(20) fit with a
    % 20-row batch overlap against the exact Toeplitz whitener of the same law;
    % the estimates differ by that approximation, bounded at 1e-3 relative
    tc.verifyEqual(sAR.x, sTP.x, 'RelTol', 1e-3);
end

function testCheckPSD(tc)
% categories: contract
    dt = 5;  a = [1 -1.6 0.9 -0.25];  s = 0.7;
    f = linspace(1e-9, 0.5/dt, 4096).';
    H = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * a(:));
    cm = struct('type','psd','f',f,'S', s^2*dt*2*abs(H).^2);
    chk = lsqvce.engine.checkPSD(cm, dt, struct('arOrder',20), false);
    tc.verifyLessThan(chk.acvRelErr, 1e-2);
    tc.verifyTrue(chk.stable);
end

%% ------------------------------------------- conditioning & regularisation
function testQRBeatsNormalWhenIllConditioned(tc)
% categories: robustness
% QR avoids squaring the condition number. Measured on the backend in ISOLATION:
% in the full pipeline the effect is masked once regularisation engages, because
% then the regularisation sets the accuracy for both backends alike.
% See test_precision for the detailed treatment.
    muteWarn(tc, 'lsqvce:accum:pinvFallback');   % provoked on purpose here (WP-47, gate rule 7)
    muteWarn(tc, 'lsqvce:normalEq:cancellation');   % tight-fit / ill-conditioned regime provoker (WP-49 unified id)
    rs = RandStream('mt19937ar','Seed',5);
    m = 800;  n = 8;

    [U,~] = qr(randn(rs, m, n), 0);
    [V,~] = qr(randn(rs, n, n), 0);
    A = U * diag(logspace(0, -6, n)) * V.';      % cond(A) = 1e6
    l = A*randn(rs, n, 1) + 1e-10*randn(rs, m, 1);
    xr = A \ l;

    accN = lsqvce.engine.NormalEqSystem(n);  accN.add(A, l);
    accQ = lsqvce.engine.QRSystem(n);        accQ.add(A, l);

    eN = norm(accN.solve() - xr)/norm(xr);
    eQ = norm(accQ.solve() - xr)/norm(xr);

    tc.verifyLessThan(eQ, eN);
    tc.verifyLessThan(eQ, 1e-10);
end

function testDiagnosticsFlagIllConditioning(tc)
% categories: robustness
    muteWarn(tc, 'lsqvce:accum:pinvFallback');   % provoked on purpose here (WP-47, gate rule 7)
    cfg = tc.TestData.cfg;
    cfg.collinear = true;  cfg.collinearEps = 1e-9;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr','Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:2); end
    d = acc.diag;
    tc.verifyGreaterThan(d(1).condNll, 1e8);
    tc.verifyNotEqual(d(1).status, 'ok');
end

function testTikhonovReducesLocalCondition(tc)
% categories: robustness
    cfg = tc.TestData.cfg;
    cfg.collinear = true;  cfg.collinearEps = 1e-10;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Reg','tikhonov', ...
                          'RegLambda',1e-3, 'Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:2); end
    sol = acc.solve();
    tc.verifyTrue(all(isfinite(sol.x)));
    tc.verifyEqual(acc.diag(1).regUsed, 'tikhonov');
end

function testAutoRegSelectsByStatus(tc)
% categories: robustness
    cfg = tc.TestData.cfg;
    % well conditioned -> no regularisation
    a1 = lsqvce.engine.BlockAccumulator(cfg.nc,'Reg','auto','Verbose',false);
    a1.addSystem(1, cfg.nl, binder(cfg,1), 1:2);
    tc.verifyEqual(a1.diag(1).regUsed, 'none');

    % ill conditioned -> tikhonov or tsvd
    cfg.collinear = true;  cfg.collinearEps = 1e-11;
    a2 = lsqvce.engine.BlockAccumulator(cfg.nc,'Reg','auto','Verbose',false);
    a2.addSystem(1, cfg.nl, binder(cfg,1), 1:2);
    tc.verifyNotEqual(a2.diag(1).regUsed, 'none');
end

function testTSVDHandlesExactRankDeficiency(tc)
% categories: robustness, contract
% Exactly duplicated local column -> singular Nll. TSVD must still solve.
% WP-133 (review 3 finding 6): "tsvd" TRUNCATES -- the help's "equivalent to
% pinv" -- so the dof counts the kept rank and the verdict's identity holds;
% it floored the null direction at 1e-10 max and counted nl (sum(r) off by
% one per system: SUM_R ERROR, INVALID). No pinvFallback warning: truncation
% is what was asked for.
    cfg = tc.TestData.cfg;
    cfg.collinear = true;  cfg.collinearEps = 0;      % exact duplication
    acc = lsqvce.engine.BlockAccumulator(cfg.nc,'Backend','qr','Reg','tsvd','Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:2); end
    sol = acc.solve();
    tc.verifyTrue(all(isfinite(sol.x)));
    tc.verifyGreaterThanOrEqual(acc.diag(1).rankDef, 1);
    D = lsqvce.engine.verdict.fromAccumulator(acc, sol);
    kept = numel(cfg.sats) * cfg.nl - sum([acc.diag.rankDef]);
    tc.verifyEqual(sol.dof, D.n - cfg.nc - kept, 'the dof counts the kept rank');
    tc.verifyEqual(D.r, sol.dof, 'the verdict counts the same');
    % tol: the report's own integrity line (report.m: SUM_R > 1e-6 is INVALID)
    tc.verifyLessThan(D.sumRRelErr, 1e-6, 'sum(r) = n - (nc + kept rank)');
end

function testRidgeCountsItsEffectiveRedundancy(tc)
% categories: precision, metamorphic
% WP-133 (review 3 finding 6): a ridge on a local block ("tikhonov", "lm")
% takes tr(Q Lambda) out of the hat matrix's trace, so the redundancy is
% m - u + sum tr(Q_ll Lambda), not m - u: the verdict's sum(r) identity read
% SUM_R ERROR / INVALID for every regularised solve. The solve, the residual
% pass and the verdict count the same effective dof, and s0^2 = vtv / dof
% (mirror: dev/validation/wp133_regularised_dof_mirror.py, 1e-12).
    cfg = tc.TestData.cfg;
    cfg.collinear = true;  cfg.collinearEps = 1e-10;
    for reg = {'tikhonov', 'lm'}
        acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend', 'qr', 'Reg', reg{1}, ...
                                             'RegLambda', 1e-3, 'Verbose', false);
        for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:2); end
        sol = acc.solve();
        tc.verifyGreaterThan(sol.regDefect, 0, [reg{1} ': the ridge has a defect']);
        D = lsqvce.engine.verdict.fromAccumulator(acc, sol);
        tc.verifyGreaterThan(sol.dof, D.n - D.u, [reg{1} ': more than m - u']);
        % tol: 1e-12 relative -- the same defect summed by the solve and the bridge
        tc.verifyEqual(D.r, sol.dof, 'RelTol', 1e-12);
        % tol: 1e-12 relative -- s0^2 is vtv / dof by construction, two products
        tc.verifyEqual(sol.s0sq * sol.dof, sol.vtv, 'RelTol', 1e-12);
        % tol: the report's own integrity line (SUM_R > 1e-6 is INVALID)
        tc.verifyLessThan(D.sumRRelErr, 1e-6, [reg{1} ': sum(r) = m - u + tr(Q Lambda)']);
        res = tc.verifyWarningFree(@() lsqvce.engine.residuals(acc, sol, 'Verbose', false));
        % tol: 1e-12 relative -- the residual pass counts the same dof
        tc.verifyEqual(res.dof, sol.dof, 'RelTol', 1e-12);
    end
end

%% ----------------------------------------------------------------- residuals
function testResidualsReproduceSigma0(tc)
% categories: precision
% The second pass must recompute exactly the same sigma0 as accumulation.
    cfg = tc.TestData.cfg;
    days = 1:3;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc,'Backend','qr','Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), days); end
    sol = acc.solve();
    res = acc.residuals(sol, 'Verbose', false);
    % tol: two backends on one well-conditioned problem; difference ~ cond*eps, bound with release-drift headroom
    tc.verifyEqual(sqrt(res.s0sq), sol.sigma0, 'RelTol', 1e-8);
    % WP-131: 'Whitened', false only leaves e.w out -- it made vtv / sigma0 the
    % RAW sums, warned sigma0Mismatch and set supersedesSol (coloured arcs here)
    res2 = tc.verifyWarningFree(@() acc.residuals(sol, 'Verbose', false, 'Whitened', false));
    % tol: 1e-12 relative -- the same whitened sum, formed twice
    tc.verifyEqual(res2.s0sq, res.s0sq, 'RelTol', 1e-12, 'the scalars are the whitened ones');
end

function testResidualsAreWhite(tc)
% categories: precision
% With a correctly specified covariance the whitened residuals must have
% unit variance and negligible lag-1 autocorrelation.
    cfg = tc.TestData.cfg;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc,'Backend','qr','Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:5); end
    sol = acc.solve();
    res = acc.residuals(sol, 'Verbose', false);

    w = vertcat(res.sys(1).obs.w);
    % tol: statistical identity on a finite sample; scatter ~ 1/sqrt(dof) of the fixture, bar sized 4/sqrt(Z)-style (HANDOVER 6.1.2)
    tc.verifyEqual(std(w), 1, 'AbsTol', 0.12);
    rho1 = (w(1:end-1).'*w(2:end)) / (w.'*w);
    tc.verifyLessThan(abs(rho1), 0.1);
end

function testLocalCovarianceIsPD(tc)
% categories: precision
    cfg = tc.TestData.cfg;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc,'Backend','qr','Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:2); end
    sol = acc.solve();
    Cyy = acc.localCovariance(cfg.sats(1), sol.Cxx, sol.s0sq);
    % tol: two backends on one well-conditioned problem; difference ~ cond*eps, bound with release-drift headroom
    tc.verifyEqual(Cyy, Cyy.', 'AbsTol', 1e-12);
    tc.verifyGreaterThan(min(eig(Cyy)), 0);
end

%% ---------------------------------------------------------------- plumbing
function testUserOptReachesLoader(tc)
% categories: contract
% opt must be forwarded to the loader on every call.
    cfg = tc.TestData.cfg;
    seen = false;
    function out = probe(d, ot, a, uo)
        if ~(isfield(uo,'myFlag') && uo.myFlag == 77)
            error('opt not forwarded');
        end
        if ~(isfield(uo,'sysID') && isfield(uo,'day') || ischar(d))
            error('context not injected');
        end
        seen = true;
        out = lsqvce.engine.synthLoader(cfg, 1, d, ot, a, uo);
    end
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose',false, ...
                          'UserOpt', struct('myFlag',77));
    acc.addSystem(1, cfg.nl, @probe, 1);
    tc.verifyTrue(seen);
end

function testEmptyArcIsSkipped(tc)
% categories: robustness, vectorisation
    cfg = tc.TestData.cfg;
    function out = gappy(d, ot, a, uo)
        out = lsqvce.engine.synthLoader(cfg, 1, d, ot, a, uo);
        if ~ischar(d) && a == 2
            out.A = zeros(0,cfg.nc); out.B = zeros(0,cfg.nl); out.l = [];
        end
    end
    acc = lsqvce.engine.BlockAccumulator(cfg.nc,'Verbose',false);
    acc.addSystem(1, cfg.nl, @gappy, 1:2);
    sol = acc.solve();
    tc.verifyTrue(all(isfinite(sol.x)));
end

%% ------------------------------------------------------------------ helper
function [Af, lf] = assembleFull(cfg, days)
%ASSEMBLEFULL  Explicit whitened design matrix for [x_c; y_1; ...; y_S].
%   Reference implementation used only by the tests.
    ns = numel(cfg.sats);
    Af = [];  lf = [];
    uo = struct('BatchSize',1e9);
    for k = 1:ns
        s = cfg.sats(k);
        f = binder(cfg, s);
        for d = days(:).'
            meta = f('meta', d, [], uo);
            for ot = 1:numel(meta.obsTypes)
                for a = meta.arcs{ot}
                    arc = f(d, meta.obsTypes{ot}, a, uo);
                    m = size(arc.A,1);
                    if isfield(arc,'Wfcn') && ~isempty(arc.Wfcn)
                        [Aw,Bw,lw] = arc.Wfcn(1:m, arc.A, arc.B, arc.l, uo);
                    else
                        W = lsqvce.kernel.decorrFactor(arc.Cfcn(1:m, uo));
                        Aw = W*arc.A;  Bw = W*arc.B;  lw = W*arc.l;
                    end
                    blk = zeros(m, cfg.nc + ns*cfg.nl);
                    blk(:,1:cfg.nc) = Aw;
                    blk(:, cfg.nc + (k-1)*cfg.nl + (1:cfg.nl)) = Bw;
                    Af = [Af; blk];   %#ok<AGROW>
                    lf = [lf; lw];    %#ok<AGROW>
                end
            end
        end
    end
end

%% ============================================ verbosity is inherited, not reset
function testResidualsInheritTheAccumulatorVerbosity(tc)
% categories: robustness
    % REGRESSION GUARD (2026-07-18). LSQVCE.ENGINE.RESIDUALS defaults its own 'verbose'
    % to TRUE -- right for a direct call, wrong for one made THROUGH an
    % accumulator that was already told to be quiet. Because the forwarder did
    % not pass the level on, an adjustVCE run with verbose OFF still emitted a
    % per-day residual line for every system and day: 544 lines in the
    % 2026-07-18 suite log, and the ONLY output that "silent" run produced.
    cfg  = baseCfg();
    days = 1:2;
    acc  = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg,s), days);
    end
    sol = acc.solve();

    f = [tempname, '.log'];  fid = fopen(f, 'w');
    c = onCleanup(@() cleanupResidLog(fid, f));
    lsqvce.util.log.set('fid', fid, 'screen', false);
    acc.residuals(sol);                 % NO 'Verbose' -> must inherit silence
    lsqvce.util.log.reset();

    tc.verifyEmpty(strtrim(fileread(f)), ...
        'a silent accumulator must produce a silent residual pass');
end

function testResidualsExplicitVerboseStillWins(tc)
% categories: robustness
    % Inheriting must not take the choice away: an explicit 'Verbose' passed by
    % the caller overrides the accumulator's level in both directions.
    cfg  = baseCfg();
    acc  = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg,s), 1:2);
    end
    sol = acc.solve();

    f = [tempname, '.log'];  fid = fopen(f, 'w');
    c = onCleanup(@() cleanupResidLog(fid, f));
    lsqvce.util.log.set('fid', fid, 'screen', false);
    acc.residuals(sol, 'Verbose', true);      % opt IN despite a quiet acc
    lsqvce.util.log.reset();

    tc.verifyNotEmpty(strtrim(fileread(f)), ...
        'an explicit Verbose must override the inherited level');
end

function cleanupResidLog(fid, f)
%CLEANUPRESIDLOG  Restore the global log state and remove the temp file.
    lsqvce.util.log.reset();
    if fid > 0, fclose(fid); end
    if exist(f, 'file'), delete(f); end
end

%% ---------------------------------- R6 E-09 (WP-91): whole epochs in EVERY batch loop
function testMultiAxisArcLongerThanBatchInEverySweep(tc)
% categories: contract, robustness, vectorisation
    % A 3-axis arc (60 epochs = 180 rows) with a per-epoch Cep, whitened in
    % batches of 100 ROWS. The accumulator rounded its batch down to whole
    % epochs; the residual sweep and the diagnostic re-read had copied the
    % loop WITHOUT the rounding (finding E-09) and refused
    % lsqvce:covModel:notWholeEpochs on the split epoch. One iterator now:
    % batched == single batch (the same rows through the same filter, the
    % overlap = the filter memory in rows), eps-level.
    d = 3;  nE = 60;  m = nE*d;  nc = 3;  nl = 2;
    rs = RandStream('mt19937ar', 'Seed', 9);
    Cep = 0.5 + rand(rs, nE, d);
    D = struct('A', randn(rs, m, nc), 'B', randn(rs, m, nl), 'l', randn(rs, m, 1), ...
               't', reshape(repmat((0:nE-1)*5, d, 1), [], 1), 'Cep', Cep);
    D = arcCells(D);                                 % v6: the cells are the arcs
    cov.pos = struct('type', 'arma', 'ar', [1 -0.6], 'sigma', 0.8);
    L = lsqvce.engine.makeLoader(@(day, ot) D, 'ObservationTypes', {'pos'}, 'Covariance', cov, ...
                                 'Arcs', 1, 'SystemId', 1, 'NumLocal', nl, 'NumCommon', nc, 'Days', 1);
    [xRef, wRef, Aref] = batchedSweeps(L, nc, nl, 1e9);
    [x100, w100, A100] = batchedSweeps(L, nc, nl, 100);
    % tol: the same rows through the same AR(1) filter with the overlap equal
    % to the filter memory -- identical arithmetic up to summation order (1e-12 house)
    tc.verifyEqual(x100, xRef, 'RelTol', 1e-12);
    tc.verifyEqual(w100, wRef, 'RelTol', 1e-12);   % tol: same source
    tc.verifyEqual(A100, Aref, 'RelTol', 1e-12);   % tol: same source
end

function [x, w, Aw] = batchedSweeps(L, nc, nl, bs)
%BATCHEDSWEEPS  Solve, residual sweep and diagnostic re-read at one BatchSize.
    acc = lsqvce.engine.BlockAccumulator(nc, 'Backend', 'qr', 'BatchSize', bs, ...
                                         'KeepLocal', true, 'Verbose', false);
    acc.addSystem(1, nl, L.closure, 1);
    sol = acc.solve();
    x = sol.x;
    res = lsqvce.engine.residuals(acc, sol, 'verbose', false);
    w = res.sys(1).obs(1).w;
    [~, batchFcn, nB] = lsqvce.engine.verdict.source(acc, sol);
    Aw = [];
    for k = 1:nB
        Ak = batchFcn(k);
        Aw = [Aw; Ak];                                            %#ok<AGROW>
    end
end

%% ---------------------------------- WP-133 (review 3): the solve core's edges
function L = rowsLoader(A, B, l)
%ROWSLOADER  One system, one day, one arc of unit-weight rows [A | B] -> l.
    m = numel(l);
    D = struct('A', {{A}}, 'B', {{B}}, 'l', {{l}}, 't', {{(1:m).'}});
    cov.obs = struct('type', 'diag', 'sigma', 1);
    L = lsqvce.engine.makeLoader(@(day, ot) D, 'ObservationTypes', {'obs'}, 'Covariance', cov, ...
                                 'Arcs', 1, 'SystemId', 1, 'NumLocal', size(B, 2), ...
                                 'NumCommon', size(A, 2), 'Days', 1);
end

function testLocalBlockSizeRefusesACoupledBlock(tc)
% categories: robustness, metamorphic
% WP-133 (review 3): LocalBlockSize = k solves the local block as nl/k
% independent blocks and IGNORES whatever couples them -- a local column
% spanning two blocks put x 27 formal sd off without a word. Declared
% structure is now verified: a coupled block is warned and solved whole
% (bitwise the undeclared solve); a truly block-diagonal one keeps the fast
% path and the dense answer.
    rs = RandStream('mt19937ar', 'Seed', 1332);
    nc = 3;  m = 40;
    A = {randn(rs, m, nc), randn(rs, m, nc)};  l = {randn(rs, m, 1), randn(rs, m, 1)};
    Bdiag = {[randn(rs, m, 2), zeros(m, 2)], [zeros(m, 2), randn(rs, m, 2)]};
    Bfull = {randn(rs, m, 4), randn(rs, m, 4)};
    for coupled = [false true]
        Bk = Bdiag;  if coupled, Bk = Bfull; end
        L = lsqvce.engine.batchSystem(nc, @(k) deal([A{k}, Bk{k}], l{k}), 2, 'LocalColumns', nc + (1:4));
        x = cell(1, 2);  lbs = [0 2];                 % 0: undeclared (the dense Schur)
        for blk = 1:2
            acc = lsqvce.engine.BlockAccumulator(nc, 'Backend', 'normal', 'Verbose', false, ...
                                                 'LocalBlockSize', lbs(blk));
            if blk == 2 && coupled
                tc.verifyWarning(@() acc.addSystem(L.sysID, 4, L.closure, 1:2), ...
                                 'lsqvce:accum:notBlockDiagonal');
            else
                acc.addSystem(L.sysID, 4, L.closure, 1:2);
            end
            sol = acc.solve();  x{blk} = sol.x;
        end
        if coupled
            tc.verifyEqual(x{2}, x{1}, 'a coupled block is solved whole, as if undeclared');
        else
            % tol: 1e-12 relative -- block vs dense inverse of the same
            % well-conditioned 4 x 4 (cond ~ 10), rounding only
            tc.verifyEqual(x{2}, x{1}, 'RelTol', 1e-12, 'a block-diagonal block keeps its answer');
        end
    end
end

function testRankDecisionIsUnitFree(tc)
% categories: precision, metamorphic
% WP-151 (review 12.84, S1). A local block in mixed units -- a bias and a
% curvature in seconds squared, B = [1, t^2] -- is well determined, but the
% raw rank tests read its scale as a defect: the normal backend's rcond(Nll)
% from t ~ 2e4 s, the QR backend's pivots against the largest from 2e7 s.
% Both cut the t^2 direction (pinvFallback): x 10.7 formal sd off the
% least-squares solution, dof m - 5 instead of m - 6, sigma0 1.97. The tests
% now read the Jacobi-equilibrated block. THE ORACLE: the column-equilibrated
% least-squares solution, exact in any units. Silent, both backends, every
% range; and the diagnostics (through fit) say 'ok' and close sum(r).
    m = 2000;
    rs = RandStream('mt19937ar', 'Seed', 11);
    A = randn(rs, m, 4) + 0.5;  e = randn(rs, m, 1);
    for tmax = [2e4 2e8]
        t = linspace(0, tmax, m).';
        B = [ones(m, 1), t.^2];
        l = A * [1; -2; 0.5; 3] + B * [3; 4 / tmax^2] + e;
        M = [A B];  c = sqrt(sum(M.^2));  ref = ((M ./ c) \ l) ./ c.';
        for be = {'normal', 'qr'}
            acc = lsqvce.engine.BlockAccumulator(4, 'Backend', be{1}, 'Verbose', false);
            L = rowsLoader(A, B, l);
            tc.verifyWarningFree(@() acc.addSystem(1, 2, L.closure, 1), ...
                sprintf('%s, t to %g s: no pinvFallback, no nearlySingular', be{1}, tmax));
            sol = tc.verifyWarningFree(@() acc.solve());
            z = (sol.x - ref(1:4)) ./ sqrt(diag(sol.Cxx));
            % tol: 1e-8 formal sd -- measured 1.6e-13 to 3.5e-13; was 10.7
            tc.verifyLessThan(max(abs(z)), 1e-8, sprintf('%s, t to %g s: x is the least-squares x', be{1}, tmax));
            tc.verifyEqual(sol.dof, m - 6, sprintf('%s, t to %g s: both locals estimated', be{1}, tmax));
        end
    end
    dg = lsqvce.engine.localDiagnostics(B.' * B, struct());
    tc.verifyEqual(dg.status, 'ok', 'the equilibrated block is well conditioned');
    % tol: the equilibrated condition of [1, t^2] on a uniform grid is ~ 7
    tc.verifyLessThan(dg.condNll, 100);
    data = struct('commonParams', {[repmat({'common'}, 4, 1), {'c1'; 'c2'; 'c3'; 'c4'}]}, 'days', 1);
    data.systems = struct('buildDay', @(day, ot) struct('A', {{A}}, 'B', {{B}}, 'l', {{l}}, 't', {{t}}), ...
        'covariance', struct('o', struct('type', 'diag', 'sigma', 1)), 'arcs', 1, 'id', 1, ...
        'name', 'one', 'localParams', {{'arc', 'b1'; 'arc', 'b2'}});
    R = tc.verifyWarningFree(@() lsqvce.engine.fitDeclaration(data, lsqvce.Options(Backend="normal", Verbose=0)));
    % tol: the verdict's own sum(r) integrity line is 1e-6; the cofactor
    % reads the same equilibrated block (BlockCofactor)
    tc.verifyLessThan(R.diag.sumRRelErr, 1e-9, 'the cofactor counts both locals');
end

function testPriorWithoutRedundancyIsRefused(tc)
% categories: contract, robustness
% WP-133 (review 3): with a prior the dof is m + nk - nc; at 0 the backends
% themselves refuse (noRedundancy), but the prior paths returned
% s0sq = 0 / max(dof, 1) = 0 and Cxx = 0 -- formal errors of exactly zero.
    A = [1 0 0.3; 0 1 -0.2];  l = [1; 2];               % 2 rows, 3 commons
    L = lsqvce.engine.batchSystem(3, @(k) deal(A, l), 1);
    for be = {'qr', 'normal'}
        acc = lsqvce.engine.BlockAccumulator(3, 'Backend', be{1}, 'Verbose', false);
        acc.addSystem(L.sysID, 0, L.closure, 1);
        acc.setPrior(lsqvce.engine.CommonPrior([0 0 1], 0.5, 1));   % 2 + 1 - 3 = 0
        tc.verifyError(@() acc.solve(), 'lsqvce:accum:noRedundancy', be{1});
    end
end

function testUnobservedLocalBlockIsRankZero(tc)
% categories: robustness, precision
% WP-133 (review 3): a local block the data never see reduces on the factor
% to rank 0 -- Rll 0 x 0 and Vll nl x 0, both EMPTY -- and every reader's
% isempty test sent it to the singular Nll_r: NaN locals, NaN residuals of
% the system. Rank 0 is a zero block: y = 0 (minimum norm), no dof spent,
% the commons those of the data alone.
    muteWarn(tc, 'lsqvce:accum:pinvFallback');   % the rank-0 block IS the pinv case (WP-47, gate rule 7)
    rs = RandStream('mt19937ar', 'Seed', 1334);
    A = randn(rs, 30, 3);  l = A * [1; 2; 3] + 0.1 * randn(rs, 30, 1);
    L = rowsLoader(A, zeros(30, 2), l);
    acc = lsqvce.engine.BlockAccumulator(3, 'Backend', 'qr', 'Verbose', false);
    acc.addSystem(1, 2, L.closure, 1);
    sol = acc.solve();
    % tol: 1e-12 relative -- one QR of the same 30 x 3 rows, cond ~ 2
    tc.verifyEqual(sol.x, A \ l, 'RelTol', 1e-12, 'the commons of the data alone');
    tc.verifyEqual(sol.dof, 30 - 3, 'an unseen block spends no dof');
    tc.verifyEqual(acc.backSubstitute(1, sol.x), zeros(2, 1), 'the minimum-norm locals');
    res = lsqvce.engine.residuals(acc, sol, 'Verbose', false);
    tc.verifyTrue(all(isfinite(res.sys(1).obs(1).v)), 'the residuals are finite');
    D = lsqvce.engine.verdict.fromAccumulator(acc, sol);
    tc.verifyEqual(D.r, 27);
    % tol: the leverages of 30 rows on a cond ~ 2 design, rounding only
    tc.verifyLessThan(D.sumRRelErr, 1e-12);
end

function testSparseLocalDesignSolvesAndSnoops(tc)
% categories: robustness, metamorphic
% WP-133 (review 3 finding 3): buildDay may return a SPARSE B (the API
% reference says so), but the QR factor then stays sparse, and the rank
% decision's svd and snoop's cholupdate take full matrices only -- both
% died. A sparse B gives the full B's answer: a full-rank block through
% snoop's downdate, a rank-deficient one through the factor's SVD.
    rs = RandStream('mt19937ar', 'Seed', 1333);
    m = 60;  A = randn(rs, m, 3);  B = randn(rs, m, 2);
    l = A * [1; -2; 0.5] + B * [0.3; -0.7] + randn(rs, m, 1);
    l(17) = l(17) + 20;                                 % one blunder for snoop to remove
    out = cell(1, 2);
    for k = 1:2
        Bk = B;  if k == 2, Bk = sparse(B); end
        L = rowsLoader(A, Bk, l);
        acc = lsqvce.engine.BlockAccumulator(3, 'Backend', 'qr', 'Verbose', false);
        acc.addSystem(1, 2, L.closure, 1);
        sol = acc.solve();
        [D, dsol, bf] = lsqvce.engine.verdict.fromAccumulator(acc, sol);
        [solS, ~, REJ] = lsqvce.engine.verdict.snoop(dsol, D, bf, 'Verbose', false);
        out{k} = struct('x', sol.x, 'xs', solS.x, 'idx', REJ.idx);
    end
    tc.verifyEqual(out{1}.idx(1), 17, 'the planted blunder goes first');
    tc.verifyEqual(out{2}.idx, out{1}.idx, 'sparse B: snoop rejects the same observations');
    % tol: 1e-10 relative -- the same rows, sparse vs full storage of B
    tc.verifyEqual(out{2}.x, out{1}.x, 'RelTol', 1e-10);
    tc.verifyEqual(out{2}.xs, out{1}.xs, 'RelTol', 1e-10);   % tol: as above
    % a duplicated local column: the rank is decided on the (sparse) factor
    muteWarn(tc, 'lsqvce:accum:pinvFallback');   % the rank-deficient block IS the pinv case
    Bd = [B(:, 1), B(:, 1)];
    xd = cell(1, 2);
    for k = 1:2
        Bk = Bd;  if k == 2, Bk = sparse(Bd); end
        L = rowsLoader(A, Bk, l);
        acc = lsqvce.engine.BlockAccumulator(3, 'Backend', 'qr', 'Verbose', false);
        acc.addSystem(1, 2, L.closure, 1);
        sol = acc.solve();  xd{k} = sol.x;
        tc.verifyEqual(sol.dof, m - 3 - 1, 'the duplicated column costs one dof, not two');
    end
    tc.verifyEqual(xd{2}, xd{1}, 'RelTol', 1e-10);   % tol: as above
end

%% ================================================== WP-158: day contributions (engine)
function testDayContributionsNeedTheStoredLoaders(tc)
% categories: contract
    % dayContributions re-streams the days through the stored loaders, so an
    % accumulator built without KeepLocal has nothing to stream from.
    cfg = tc.TestData.cfg;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend', 'qr', 'KeepLocal', false, 'Verbose', false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg, s), 1:2); end
    tc.verifyError(@() acc.dayContributions(), 'lsqvce:accum:noContributionStore');
end

function testDayContributionsRefuseARegularisedLocalBlock(tc)
% categories: contract, robustness
    % A local block the accumulation regularised no longer is the block the
    % days' contributions would reproduce: refused, never approximated.
    cfg = tc.TestData.cfg;
    cfg.collinear = true;  cfg.collinearEps = 1e-10;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend', 'qr', 'Reg', 'tikhonov', ...
                                         'RegLambda', 1e-3, 'Verbose', false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg, s), 1:2); end
    tc.assertEqual(acc.diag(1).regUsed, 'tikhonov', 'the fixture reaches the regularisation');
    tc.verifyError(@() acc.dayContributions(), 'lsqvce:accum:regularisedLocals');
end

function testReduceDayRefusesUndeterminedLocals(tc)
% categories: contract, precision
    % A day whose own locals are not determined by it (two identical local
    % columns) has no reduced contribution, on either backend; a column of R
    % that carries data cannot be dropped by reducedRowsOf.
    rs = RandStream('mt19937ar', 'Seed', 4);
    m = 40;  nc = 2;  A = randn(rs, m, nc);  b = randn(rs, m, 1);  l = randn(rs, m, 1);
    for be = {'qr', 'normal'}
        if strcmp(be{1}, 'qr'), loc = lsqvce.engine.QRSystem(nc + 3);
        else,                   loc = lsqvce.engine.NormalEqSystem(nc + 3); end
        loc.add([A, b, b, zeros(m, 1)], l);              % locals 1 and 2 identical, 3 unused
        tc.verifyError(@() lsqvce.engine.BlockAccumulator.reduceDay(loc, nc, 3, be{1}, 1, 7), ...
                       'lsqvce:accum:dayLocalsRankDeficient', be{1});
    end
    q = lsqvce.engine.QRSystem(nc + 2);
    q.add([A, b, zeros(m, 1)], l);
    tc.verifyError(@() q.reducedRowsOf([], 1:nc), 'lsqvce:qr:droppedNonzero');
    % the unused column drops exactly: reducedRowsOf over the used one equals
    % reducedRows of a system that never had the unused column
    q2 = lsqvce.engine.QRSystem(nc + 1);  q2.add([A, b], l);
    r1 = q.reducedRowsOf(nc + 1, 1:nc);  r2 = q2.reducedRows(1);
    % tol: 1e-12 relative -- one re-triangularisation each, compared as normals
    tc.verifyEqual(r1.' * r1, r2.' * r2, 'RelTol', 1e-12);
end
