function tests = test_precision
%TEST_PRECISION  Numerical accuracy and exactness of the lsq toolbox.
%
%   Every test here asserts a NUMBER, not just "it ran". The tolerances are
%   the specification: if one of these loosens, something regressed.
%
%   Run:  results = runtests('tests/test_precision');  disp(table(results))
%
%   See also TEST_RELIABILITY, TEST_SPEED, TEST_LSQACCUM
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    tc.TestData.cfg = makeCfg();
end

%% ================================================== whitening exactness
function testWhitenerIsExactInverseSquareRoot(tc)
% THE central claim of the toolbox: the FIR whitening filter F satisfies
% F'*F = inv(C) EXACTLY, not approximately. Everything else rests on this.
    a = [1 -1.6 0.9 -0.25];  s = 0.7;
    for n = [20 60 300]
        r  = s^2 * lsqvce.armaACV(a, 1, n);
        C  = toeplitz(r);
        Ci = C \ eye(n);

        z = zeros(n, 0);
        [~, ~, F] = lsqvce.armaWhiten(a, 1, s, r, z, z, eye(n), true);

        relErr = norm(F.'*F - Ci) / norm(Ci);
        tc.verifyLessThan(relErr, 1e-12, ...
            sprintf('n=%d: F''F deviates from inv(C) by %.2e', n, relErr));
    end
end

function testWhitenerIsTheCholeskyFactorItself(tc)
% The FIR filter is not merely SOME square root of inv(C) -- it IS the Cholesky
% inverse factor. Both are lower triangular, and the triangular square root is
% unique up to sign, so they must coincide row for row.
%
% The steady-state rows are the AR coefficients reversed and divided by sigma;
% the leading p rows are the exact stationary startup block.
    a = [1 -1.6 0.9 -0.25];  s = 0.7;  n = 200;
    r = s^2 * lsqvce.armaACV(a, 1, n);

    W = lsqvce.decorrFactor(toeplitz(r));          % dense Cholesky inverse
    z = zeros(n, 0);
    [~, ~, F] = lsqvce.armaWhiten(a, 1, s, r, z, z, eye(n), true);

    tc.verifyEqual(F, W, 'AbsTol', 1e-9, ...
        'the FIR filter must reproduce chol(C)^-1 exactly');

    % and therefore the normal equations are identical too
    rs = RandStream('mt19937ar','Seed',7);
    X  = randn(rs, n, 6);
    l  = randn(rs, n, 1);
    [~, ~, Xf] = lsqvce.armaWhiten(a, 1, s, r, z, z, X, true);
    [~, ~, lf] = lsqvce.armaWhiten(a, 1, s, r, z, z, l, true);

    tc.verifyEqual(Xf.'*Xf, (W*X).'*(W*X), 'RelTol', 1e-11);
    tc.verifyEqual(Xf.'*lf, (W*X).'*(W*l), 'RelTol', 1e-11);
end

function testSteadyStateRowsAreTheARCoefficients(tc)
% Explicit structural check: away from the startup block, row k of the whitener
% is [a(p+1) ... a(1)]/sigma placed at columns k-p..k. If this ever breaks, the
% filter argument order has been swapped (numerator/denominator).
    a = [1 -1.6 0.9 -0.25];  s = 0.7;  n = 20;  p = 3;
    r = s^2 * lsqvce.armaACV(a, 1, n);

    z = zeros(n, 0);
    [~, ~, F] = lsqvce.armaWhiten(a, 1, s, r, z, z, eye(n), true);

    expected = fliplr(a) / s;                    % [a3 a2 a1 1]/s
    for k = (p+1):n
        tc.verifyEqual(F(k, k-p:k), expected, 'AbsTol', 1e-12, ...
            sprintf('row %d is not the AR whitening stencil', k));
    end
end

function testDecorrFactorIsExact(tc)
    r = lsqvce.armaACV([1 -0.9], 1, 120);
    C = toeplitz(r);
    W = lsqvce.decorrFactor(C);
    tc.verifyEqual(W.'*W, C \ eye(120), 'RelTol', 1e-10, 'AbsTol', 1e-9);
    tc.verifyEqual(W*C*W.', eye(120), 'AbsTol', 1e-9);
end

function testDecorrFactorRescuePathIsUsable(tc)
% An indefinite covariance must still yield a factor with W*C*W' ~ I. Flooring
% at eps (the old bug) gives cond ~ 1/eps and W*C*W' off by ~1e-1.
    rs = RandStream('mt19937ar','Seed',3);
    n  = 40;
    Q  = randn(rs, n, n);
    C  = Q*Q.';
    [V, d] = eig((C+C.')/2, 'vector');
    d(1:3) = -1e-9;                             % force indefiniteness
    C = V*diag(d)*V.';
    C = (C+C.')/2;

    W = tc.verifyWarning(@() lsqvce.decorrFactor(C), 'lsqvce:decorrFactor:notPD');

    % reconstruct the floored covariance and check W is its inverse factor
    [V2, d2] = eig((C+C.')/2, 'vector');
    d2 = max(real(d2), 1e-10*max(real(d2)));
    Cf = real(V2)*diag(d2)*real(V2).';
    tc.verifyEqual(W*Cf*W.', eye(n), 'AbsTol', 1e-6);
end

%% ================================================== signal processing
function testLevinsonExactOnKnownAR(tc)
    a0 = [1 -1.6 0.9 -0.25];  s0 = 0.7;
    r  = s0^2 * lsqvce.armaACV(a0, 1, 8);
    [a, s2] = lsqvce.levinsonDurbin(r, 3);
    tc.verifyEqual(a, a0, 'AbsTol', 1e-11);
    tc.verifyEqual(sqrt(s2), s0, 'RelTol', 1e-11);
end

function testLevinsonReflectionCoefficients(tc)
% |kappa| < 1 for every order iff the autocovariance is positive definite.
    r = lsqvce.armaACV([1 -1.6 0.9 -0.25], 1, 21);
    [~, ~, refl] = lsqvce.levinsonDurbin(r, 20);
    tc.verifyLessThan(max(abs(refl)), 1);
end

function testARMAACVAgainstClosedForm(tc)
    phi = 0.9;  th = 0.5;
    r  = lsqvce.armaACV([1 -phi], [1 th], 4);
    r0 = (1 + 2*phi*th + th^2)/(1 - phi^2);
    r1 = phi*r0 + th;
    tc.verifyEqual(r(1), r0, 'RelTol', 1e-9);
    tc.verifyEqual(r(2), r1, 'RelTol', 1e-9);
    tc.verifyEqual(r(3), phi*r1, 'RelTol', 1e-9);   % r_k = phi*r_{k-1}, k > q
end

function testPSD2ACVSatisfiesParsevalAtZeroLag(tc)
% r(0) = trapz(f, S) exactly -- an identity, not an approximation.
    dt = 5;
    f  = linspace(0, 0.5/dt, 8192).';
    H  = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * [1 -1.6 0.9 -0.25].');
    S  = 0.7^2 * dt * 2 * abs(H).^2;

    r = lsqvce.psd2acv(f, S, dt, 4);
    tc.verifyEqual(r(1), trapz(f, S), 'RelTol', 1e-12);
end

function testPSDToACVRoundTrip(tc)
% PSD of a known AR -> ACV must match the analytic ACV of that AR.
    dt = 5;  a = [1 -1.6 0.9 -0.25];  s = 0.7;
    f  = linspace(0, 0.5/dt, 8192).';
    H  = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * a(:));
    S  = s^2 * dt * 2 * abs(H).^2;

    r_num = lsqvce.psd2acv(f, S, dt, 64);
    r_ana = s^2 * lsqvce.armaACV(a, 1, 64);
    tc.verifyLessThan(norm(r_num - r_ana)/norm(r_ana), 2e-3);
end

function testCheckPSDVerdicts(tc)
    dt = 5;  a = [1 -1.6 0.9 -0.25];
    f  = linspace(0, 0.5/dt, 4096).';
    H  = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * a(:));
    cm = struct('type','psd', 'f',f, 'S', 0.7^2*dt*2*abs(H).^2);

    chk = lsqvce.checkPSD(cm, dt, struct('arOrder',20), false);
    tc.verifyEqual(chk.verdict, 'good');
    tc.verifyTrue(chk.stable);
    tc.verifyLessThan(chk.acvRelErr, 1e-3);

    % AR(1) is far too low an order for an AR(3) spectrum
    chkBad = tc.verifyWarning(@() lsqvce.checkPSD(cm, dt, struct('arOrder',1), false), ...
                              'lsqvce:checkPSD:badFit');
    tc.verifyEqual(chkBad.verdict, 'bad');
end

%% ================================================== estimation exactness
function testPreeliminationEqualsJointSolve(tc)
% Pre-elimination is algebraically exact: the reduced solution must equal the
% monolithic joint solve for [x_c; y_1; ...; y_S] to machine precision.
    cfg  = tc.TestData.cfg;
    days = 1:2;

    sol = runAcc(cfg, days, 'Backend','qr', 'Verbose',false);
    [Af, lf] = assembleFull(cfg, days);
    xf = Af \ lf;

    relErr = norm(sol.x - xf(1:cfg.nc)) / norm(xf(1:cfg.nc));
    tc.verifyLessThan(relErr, 1e-9, ...
        sprintf('Schur complement is not exact: %.2e', relErr));
end

function testBackSubstitutionExact(tc)
    cfg  = tc.TestData.cfg;
    days = 1:2;

    acc = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), days); end
    sol = acc.solve();

    [Af, lf] = assembleFull(cfg, days);
    xf = Af \ lf;

    for k = 1:numel(cfg.sats)
        y    = acc.backSubstitute(cfg.sats(k), sol.x);
        yref = xf(cfg.nc + (k-1)*cfg.nl + (1:cfg.nl));
        tc.verifyEqual(y, yref, 'RelTol', 1e-8);
    end
end

function testBatchOverlapGivesMachinePrecision(tc)
% With overlap = filter memory, batching is EXACT -- not merely accurate.
% With overlap = 0 the filter transient leaks in and costs ~12 orders.
    cfg  = tc.TestData.cfg;
    days = 1:2;

    ref = runAcc(cfg, days, 'Backend','qr','BatchSize',1e9,'Verbose',false);
    ok  = runAcc(cfg, days, 'Backend','qr','BatchSize',100, ...
                 'BatchOverlap',3, 'Verbose',false);
    bad = runAcc(cfg, days, 'Backend','qr','BatchSize',100, ...
                 'BatchOverlap',0, 'Verbose',false);

    eOK  = norm(ok.x  - ref.x)/norm(ref.x);
    eBAD = norm(bad.x - ref.x)/norm(ref.x);

    tc.verifyLessThan(eOK, 1e-12, ...
        sprintf('overlap=p should be exact, got %.2e', eOK));
    tc.verifyGreaterThan(eBAD/eOK, 1e6, ...
        'zero overlap should be catastrophically worse');
end

function testBatchSizeInvariance(tc)
% The answer must not depend on how the arc was chopped up.
    cfg = tc.TestData.cfg;
    xRef = [];
    for bs = [37 100 400 1e9]
        s = runAcc(cfg, 1:2, 'Backend','qr', 'BatchSize',bs, ...
                   'BatchOverlap',3, 'Verbose',false);
        if isempty(xRef)
            xRef = s.x;
        else
            tc.verifyEqual(s.x, xRef, 'RelTol', 1e-10, ...
                sprintf('BatchSize = %g changed the answer', bs));
        end
    end
end

function testPSDARMatchesDenseToeplitz(tc)
% The AR filter path and the dense Cholesky path must give the same normals.
    cfg = tc.TestData.cfg;
    dt  = cfg.dt;  a = [1 -1.6 0.9 -0.25];
    f   = linspace(0, 0.5/dt, 4096).';
    H   = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * a(:));
    cfg.cov.pos = struct('type','psd', 'f',f, 'S', 0.7^2*dt*2*abs(H).^2);

    sAR = runAcc(cfg, 1:2, 'Backend','qr','BatchSize',400,'BatchOverlap',20, ...
                 'UserOpt', struct('psdMode','ar','arOrder',20), 'Verbose',false);
    sTP = runAcc(cfg, 1:2, 'Backend','qr','BatchSize',400,'BatchOverlap',20, ...
                 'UserOpt', struct('psdMode','toeplitz'), 'Verbose',false);

    tc.verifyEqual(sAR.x, sTP.x, 'RelTol', 1e-6);
end

function testWhitenModeCholMatchesFilter(tc)
% Forcing the dense path must reproduce the filter path.
    cfg = tc.TestData.cfg;
    sF = runAcc(cfg, 1:2, 'Backend','qr','BatchSize',400,'BatchOverlap',10, ...
                'UserOpt', struct('whitenMode','filter'), 'Verbose',false);
    sC = runAcc(cfg, 1:2, 'Backend','qr','BatchSize',400,'BatchOverlap',10, ...
                'UserOpt', struct('whitenMode','chol'), 'Verbose',false);
    tc.verifyEqual(sF.x, sC.x, 'RelTol', 1e-6);
end

%% ================================================== conditioning
function testQRBeatsNormalOnIllConditioned(tc)
% The whole justification for the QR backend: cond(A) versus cond(A)^2.
%
% Measured on the BACKEND IN ISOLATION, because in the full pipeline the effect
% is masked -- once the local block is ill conditioned enough to matter, the
% regularisation kicks in and sets the accuracy for both backends alike (see
% testRegularisationMasksBackendChoice). The backend difference is real, but it
% only shows where the design matrix is ill conditioned AND unregularised.
    rs = RandStream('mt19937ar','Seed',5);
    m  = 800;   n = 8;

    for k = [4 6 8]                                  % cond(A) = 1e4, 1e6, 1e8
        [U,~] = qr(randn(rs, m, n), 0);
        [V,~] = qr(randn(rs, n, n), 0);
        A = U * diag(logspace(0, -k, n)) * V.';

        xt = randn(rs, n, 1);
        l  = A*xt + 1e-10*randn(rs, m, 1);
        xr = A \ l;                                   % LAPACK reference (QR based)

        accN = lsqvce.NormalEqSystem(n);   accN.add(A, l);
        accQ = lsqvce.QRSystem(n);         accQ.add(A, l);

        xN = accN.solve();
        xQ = accQ.solve();

        eN = norm(xN - xr)/norm(xr);
        eQ = norm(xQ - xr)/norm(xr);

        tc.verifyLessThan(eQ, eN, ...
            sprintf('cond = 1e%d: QR (%.2e) must beat normal (%.2e)', k, eQ, eN));
        tc.verifyLessThan(eQ, 1e-8, ...
            sprintf('cond = 1e%d: QR error %.2e is too large', k, eQ));
    end
end

function testNormalEquationsLoseTheAnswerAtHighCondition(tc)
% The failure mode the QR backend exists to prevent. At cond(A) = 1e8 the normal
% matrix has cond(A)^2 = 1e16 > 1/eps: the information is destroyed by forming
% A'*A and NO solver can recover it afterwards.
    rs = RandStream('mt19937ar','Seed',11);
    m  = 800;   n = 8;

    [U,~] = qr(randn(rs, m, n), 0);
    [V,~] = qr(randn(rs, n, n), 0);
    A = U * diag(logspace(0, -8, n)) * V.';
    l = A*randn(rs, n, 1) + 1e-10*randn(rs, m, 1);
    xr = A \ l;

    accN = lsqvce.NormalEqSystem(n);   accN.add(A, l);
    accQ = lsqvce.QRSystem(n);         accQ.add(A, l);

    eN = norm(accN.solve() - xr)/norm(xr);
    eQ = norm(accQ.solve() - xr)/norm(xr);

    tc.verifyGreaterThan(eN, 1e-3, 'the normal equations should fail badly here');
    tc.verifyLessThan(eQ, 1e-7,   'QR should still be accurate here');
    tc.verifyGreaterThan(eN/eQ, 1e4);
end

function testRegularisationMasksBackendChoice(tc)
% Worth knowing, and easy to misread: once the local block is singular enough
% that TSVD engages, BOTH backends give the same answer -- the regularisation,
% not the factorisation, sets the accuracy. Do not conclude from a full-pipeline
% test that the backend does not matter; test the backend in isolation.
    cfg = tc.TestData.cfg;
    cfg.collinear = true;   cfg.collinearEps = 1e-9;

    sN = runAcc(cfg, 1:2, 'Backend','normal', 'Reg','tsvd', 'Verbose',false);
    sQ = runAcc(cfg, 1:2, 'Backend','qr',     'Reg','tsvd', 'Verbose',false);

    tc.verifyEqual(sQ.x, sN.x, 'RelTol', 1e-6, ...
        'with TSVD engaged the two backends should agree');
end

function testBackendsAgreeWhenWellConditioned(tc)
    cfg = tc.TestData.cfg;
    sN = runAcc(cfg, 1:3, 'Backend','normal', 'Verbose',false);
    sQ = runAcc(cfg, 1:3, 'Backend','qr',     'Verbose',false);
    tc.verifyEqual(sQ.x,      sN.x,      'RelTol', 1e-8);
    tc.verifyEqual(sQ.sigma0, sN.sigma0, 'RelTol', 1e-7);
end

%% ================================================== statistical fidelity
function testSigma0IsUnbiased(tc)
% With a correctly specified covariance, sigma0 -> 1. This is the end-to-end
% check that the whitening, the batching and the ltl bookkeeping all agree.
    cfg = tc.TestData.cfg;
    sol = runAcc(cfg, 1:6, 'Backend','qr', 'Verbose',false);
    tc.verifyEqual(sol.sigma0, 1, 'AbsTol', 0.05);
end

function testResidualPassReproducesSigma0(tc)
% Two independent routes to sigma0 must agree to machine precision.
    cfg = tc.TestData.cfg;
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:3); end
    sol = acc.solve();
    res = acc.residuals(sol, 'Verbose', false);
    tc.verifyEqual(res.sigma0, sol.sigma0, 'RelTol', 1e-10);
end

function testWhitenedResidualsAreWhite(tc)
% The stochastic model is right iff the whitened residuals are white with unit
% variance. This tests the noise COLOUR, not just its level.
    cfg = tc.TestData.cfg;
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:6); end
    sol = acc.solve();
    res = acc.residuals(sol, 'Verbose', false);

    w = vertcat(res.sys(1).obs.w);
    n = numel(w);

    tc.verifyEqual(std(w), 1, 'AbsTol', 0.08);

    % lag-1..5 autocorrelation inside the 95% band of white noise
    band = 2.5/sqrt(n);
    for k = 1:5
        rho = (w(1:end-k).'*w(k+1:end)) / (w.'*w);
        tc.verifyLessThan(abs(rho), band, ...
            sprintf('lag-%d autocorrelation %.4f exceeds %.4f', k, rho, band));
    end
end

function testCovarianceIsPositiveDefinite(tc)
    cfg = tc.TestData.cfg;
    acc = lsqvce.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:3); end
    sol = acc.solve();

    tc.verifyEqual(sol.Cxx, sol.Cxx.', 'AbsTol', 1e-14);
    tc.verifyGreaterThan(min(eig(sol.Cxx)), 0);

    Cyy = acc.localCovariance(cfg.sats(1), sol.Cxx);
    tc.verifyEqual(Cyy, Cyy.', 'AbsTol', 1e-14);
    tc.verifyGreaterThan(min(eig(Cyy)), 0);
end

function testFormalErrorsMatchEmpiricalScatter(tc)
% The strongest statistical test: repeat the estimation over many noise
% realisations and check that the empirical scatter of x matches the formal
% sigma from Cxx. This validates the covariance propagation end to end.
    nRep = 60;                                  % 40 was too few for a max-of-6 test
    cfg  = tc.TestData.cfg;
    cfg.mArc = 200;   cfg.arcsPerDay = 2;

    X  = zeros(cfg.nc, nRep);
    sF = zeros(cfg.nc, nRep);
    for k = 1:nRep
        c = cfg;  c.seed = 1000 + k;             % new noise realisation
        s = runAcc(c, 1:2, 'Backend','qr', 'Verbose',false);
        X(:,k)  = s.x - cfg.x_c;                 % error w.r.t. truth
        sF(:,k) = sqrt(diag(s.Cxx));
    end

    empirical = std(X, 0, 2);
    formal    = mean(sF, 2);
    ratio     = empirical ./ formal;

    % The MEAN ratio is the stable, meaningful check -- it must be ~1, which is
    % the actual claim (formal errors are calibrated). The MAX over the 6
    % parameters has a fat tail with finite nRep (the worst of 6 sample-std
    % ratios), so it gets a looser bound. Sample scatter of each std is
    % 1/sqrt(2*(nRep-1)) ~ 9% at nRep = 60.
    tc.verifyEqual(mean(ratio), 1, 'AbsTol', 0.10, ...
        'formal errors are not calibrated on average');
    tc.verifyLessThan(max(abs(ratio - 1)), 0.35, ...
        sprintf('formal/empirical sigma ratio spans [%.2f, %.2f]', ...
                min(ratio), max(ratio)));
end

%% ================================================== helpers
function cfg = makeCfg()
    rs = RandStream('mt19937ar','Seed',42);
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
    cfg.seed         = 42;
    cfg.collinear    = false;
    cfg.collinearEps = 1e-7;
    cfg.cov.pos = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',0.7);
    cfg.cov.kbr = struct('type','diag','sigma',0.2);
end

function f = binder(cfg, s)
    f = @(d,ot,a,uo) lsqvce.synthLoader(cfg, s, d, ot, a, uo);
end

function sol = runAcc(cfg, days, varargin)
    acc = lsqvce.BlockAccumulator(cfg.nc, varargin{:});
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg,s), days);
    end
    sol = acc.solve();
end

function [Af, lf] = assembleFull(cfg, days)
%ASSEMBLEFULL  Explicit whitened design matrix for [x_c; y_1; ...; y_S].
%   Reference implementation, used only by the tests. This is what the
%   accumulator must reproduce without ever building it.
    ns = numel(cfg.sats);
    Af = [];  lf = [];
    uo = struct('BatchSize', 1e9);

    for k = 1:ns
        f = binder(cfg, cfg.sats(k));
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
                    blk(:, 1:cfg.nc) = Aw;
                    blk(:, cfg.nc + (k-1)*cfg.nl + (1:cfg.nl)) = Bw;
                    Af = [Af; blk];   %#ok<AGROW>
                    lf = [lf; lw];    %#ok<AGROW>
                end
            end
        end
    end
end
