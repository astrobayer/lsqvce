function tests = test_lsqAccum
%TEST_LSQACCUM  Unit tests for the lsq accumulation toolbox.
%   Run:  results = runtests('tests/test_lsqAccum');  disp(table(results))
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
    f = @(d,ot,a,uo) lsqvce.synthLoader(cfg, s, d, ot, a, uo);
end

function sol = runIt(cfg, days, varargin)
    acc = lsqvce.BlockAccumulator(cfg.nc, varargin{:});
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
    r_ana = s^2 * lsqvce.armaACV(a, 1, 64);
    f = linspace(1e-9, 0.5/dt, 8192).';
    H = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * a(:));
    S = s^2 * dt * 2 * abs(H).^2;
    r_num = lsqvce.psd2acv(f, S, dt, 64);
    tc.verifyLessThan(norm(r_num-r_ana)/norm(r_ana), 2e-3);
end

function testLevinsonRecoversAR(tc)
% categories: precision
% Levinson-Durbin on the exact ACV must return the generating coefficients.
    a = [1 -1.6 0.9 -0.25];  s = 0.7;
    r = s^2 * lsqvce.armaACV(a, 1, 8);
    [ah, s2] = lsqvce.levinsonDurbin(r, 3);
    % tol: Levinson on the EXACT generator ACV; recursion error growth over p steps ~ 1e3*eps, bound adds headroom
    tc.verifyEqual(ah, a, 'AbsTol', 1e-8);
    tc.verifyEqual(sqrt(s2), s, 'RelTol', 1e-8);
end

function testFilterEqualsCholeskyWhitening(tc)
% categories: precision
% THE key invariant: the FIR whitening filter and the dense Cholesky
% decorrelation give the SAME Y'*Y (they differ by an orthogonal rotation).
    a = [1 -1.6 0.9 -0.25];  s = 0.7;  n = 300;
    r = s^2 * lsqvce.armaACV(a, 1, n);
    X = randn(n, 5);
    W  = lsqvce.decorrFactor(toeplitz(r));
    Y1 = W*X;
    z  = zeros(n,0);
    [~,~,Y2] = lsqvce.armaWhiten(a, 1, s, r, z, z, X, true);
    % tol: dense reconstruction of the same algebra; eps-level accumulation, house bounds (HANDOVER 6.1.2)
    tc.verifyEqual(Y2.'*Y2, Y1.'*Y1, 'RelTol', 1e-9);
end

function testARMAACVAgainstBruteForce(tc)
% categories: precision, reference
% ARMA ACV via the FFT/Wold expansion vs a direct long simulation moment.
    a = [1 -0.9];  b = [1 0.5];
    r = lsqvce.armaACV(a, b, 5);
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
    % tol: two backends on one well-conditioned problem; difference ~ cond*eps, bound with release-drift headroom
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
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
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

    % tol: two backends on one well-conditioned problem; difference ~ cond*eps, bound with release-drift headroom
    tc.verifyEqual(sAR.x, sTP.x, 'RelTol', 1e-3);
end

function testCheckPSD(tc)
% categories: contract
    dt = 5;  a = [1 -1.6 0.9 -0.25];  s = 0.7;
    f = linspace(1e-9, 0.5/dt, 4096).';
    H = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * a(:));
    cm = struct('type','psd','f',f,'S', s^2*dt*2*abs(H).^2);
    chk = lsqvce.checkPSD(cm, dt, struct('arOrder',20), false);
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
    muteWarn(tc, 'lsqvce:normalEq:negativeVtv');   % provoked on purpose here (WP-47, gate rule 7)
    rs = RandStream('mt19937ar','Seed',5);
    m = 800;  n = 8;

    [U,~] = qr(randn(rs, m, n), 0);
    [V,~] = qr(randn(rs, n, n), 0);
    A = U * diag(logspace(0, -6, n)) * V.';      % cond(A) = 1e6
    l = A*randn(rs, n, 1) + 1e-10*randn(rs, m, 1);
    xr = A \ l;

    accN = lsqvce.NormalEqSystem(n);  accN.add(A, l);
    accQ = lsqvce.QRSystem(n);        accQ.add(A, l);

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
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr','Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:2); end
    d = acc.diag;
    tc.verifyGreaterThan(d(1).condNll, 1e8);
    tc.verifyNotEqual(d(1).status, 'ok');
end

function testTikhonovReducesLocalCondition(tc)
% categories: robustness
    cfg = tc.TestData.cfg;
    cfg.collinear = true;  cfg.collinearEps = 1e-10;
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr', 'Reg','tikhonov', ...
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
    a1 = lsqvce.BlockAccumulator(cfg.nc,'Reg','auto','Verbose',false);
    a1.addSystem(1, cfg.nl, binder(cfg,1), 1:2);
    tc.verifyEqual(a1.diag(1).regUsed, 'none');

    % ill conditioned -> tikhonov or tsvd
    cfg.collinear = true;  cfg.collinearEps = 1e-11;
    a2 = lsqvce.BlockAccumulator(cfg.nc,'Reg','auto','Verbose',false);
    a2.addSystem(1, cfg.nl, binder(cfg,1), 1:2);
    tc.verifyNotEqual(a2.diag(1).regUsed, 'none');
end

function testTSVDHandlesExactRankDeficiency(tc)
% categories: robustness
% Exactly duplicated local column -> singular Nll. TSVD must still solve.
    muteWarn(tc, 'lsqvce:accum:pinvFallback');   % provoked on purpose here (WP-47, gate rule 7)
    cfg = tc.TestData.cfg;
    cfg.collinear = true;  cfg.collinearEps = 0;      % exact duplication
    acc = lsqvce.BlockAccumulator(cfg.nc,'Backend','qr','Reg','tsvd','Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:2); end
    sol = acc.solve();
    tc.verifyTrue(all(isfinite(sol.x)));
    tc.verifyGreaterThanOrEqual(acc.diag(1).rankDef, 1);
end

%% ----------------------------------------------------------------- residuals
function testResidualsReproduceSigma0(tc)
% categories: precision
% The second pass must recompute exactly the same sigma0 as accumulation.
    cfg = tc.TestData.cfg;
    days = 1:3;
    acc = lsqvce.BlockAccumulator(cfg.nc,'Backend','qr','Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), days); end
    sol = acc.solve();
    res = acc.residuals(sol, 'Verbose', false);
    % tol: two backends on one well-conditioned problem; difference ~ cond*eps, bound with release-drift headroom
    tc.verifyEqual(sqrt(res.s0sq), sol.sigma0, 'RelTol', 1e-8);
end

function testResidualsAreWhite(tc)
% categories: precision
% With a correctly specified covariance the whitened residuals must have
% unit variance and negligible lag-1 autocorrelation.
    cfg = tc.TestData.cfg;
    acc = lsqvce.BlockAccumulator(cfg.nc,'Backend','qr','Verbose',false);
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
    acc = lsqvce.BlockAccumulator(cfg.nc,'Backend','qr','Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:2); end
    sol = acc.solve();
    Cyy = acc.localCovariance(cfg.sats(1), sol.Cxx);
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
        out = lsqvce.synthLoader(cfg, 1, d, ot, a, uo);
    end
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Verbose',false, ...
                          'UserOpt', struct('myFlag',77));
    acc.addSystem(1, cfg.nl, @probe, 1);
    tc.verifyTrue(seen);
end

function testEmptyArcIsSkipped(tc)
% categories: robustness, vectorisation
    cfg = tc.TestData.cfg;
    function out = gappy(d, ot, a, uo)
        out = lsqvce.synthLoader(cfg, 1, d, ot, a, uo);
        if ~ischar(d) && a == 2
            out.A = zeros(0,cfg.nc); out.B = zeros(0,cfg.nl); out.l = [];
        end
    end
    acc = lsqvce.BlockAccumulator(cfg.nc,'Verbose',false);
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
                        W = lsqvce.decorrFactor(arc.Cfcn(1:m, uo));
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
    % REGRESSION GUARD (2026-07-18). LSQVCE.RESIDUALS defaults its own 'verbose'
    % to TRUE -- right for a direct call, wrong for one made THROUGH an
    % accumulator that was already told to be quiet. Because the forwarder did
    % not pass the level on, an adjustVCE run with verbose OFF still emitted a
    % per-day residual line for every system and day: 544 lines in the
    % 2026-07-18 suite log, and the ONLY output that "silent" run produced.
    cfg  = baseCfg();
    days = 1:2;
    acc  = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
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
    acc  = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
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
