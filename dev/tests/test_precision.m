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
% categories: contract
% THE central claim of the toolbox: the FIR whitening filter F satisfies
% F'*F = inv(C) EXACTLY, not approximately. Everything else rests on this.
    a = [1 -1.6 0.9 -0.25];  s = 0.7;
    for n = [20 60 300]
        r  = s^2 * lsqvce.kernel.armaACV(a, 1, n);
        C  = toeplitz(r);
        Ci = C \ eye(n);

        z = zeros(n, 0);
        [~, ~, F] = lsqvce.kernel.armaWhiten(a, 1, s, r, z, z, eye(n), true);

        relErr = norm(F.'*F - Ci) / norm(Ci);
        tc.verifyLessThan(relErr, 1e-12, ...
            sprintf('n=%d: F''F deviates from inv(C) by %.2e', n, relErr));
    end
end

function testWhitenerIsTheCholeskyFactorItself(tc)
% categories: precision
% The FIR filter is not merely SOME square root of inv(C) -- it IS the Cholesky
% inverse factor. Both are lower triangular, and the triangular square root is
% unique up to sign, so they must coincide row for row.
%
% The steady-state rows are the AR coefficients reversed and divided by sigma;
% the leading p rows are the exact stationary startup block.
    a = [1 -1.6 0.9 -0.25];  s = 0.7;  n = 200;
    r = s^2 * lsqvce.kernel.armaACV(a, 1, n);

    W = lsqvce.kernel.decorrFactor(toeplitz(r));          % dense Cholesky inverse
    z = zeros(n, 0);
    [~, ~, F] = lsqvce.kernel.armaWhiten(a, 1, s, r, z, z, eye(n), true);

    % tol: both sides are Cholesky factors of the same Toeplitz C (cond ~1e3);
    % eps*cond*|W| ~ 1e-12; 1e-9 adds 3 decades of headroom
    tc.verifyEqual(F, W, 'AbsTol', 1e-9, ...
        'the FIR filter must reproduce chol(C)^-1 exactly');

    % and therefore the normal equations are identical too
    rs = RandStream('mt19937ar','Seed',7);
    X  = randn(rs, n, 6);
    l  = randn(rs, n, 1);
    [~, ~, Xf] = lsqvce.kernel.armaWhiten(a, 1, s, r, z, z, X, true);
    [~, ~, lf] = lsqvce.kernel.armaWhiten(a, 1, s, r, z, z, l, true);

    % tol: quadratic forms of the same whitener, eps-level; 1e-11 = 1e4*eps
    tc.verifyEqual(Xf.'*Xf, (W*X).'*(W*X), 'RelTol', 1e-11);
    tc.verifyEqual(Xf.'*lf, (W*X).'*(W*l), 'RelTol', 1e-11);  % tol: same source
end

function testSteadyStateRowsAreTheARCoefficients(tc)
% categories: precision
% Explicit structural check: away from the startup block, row k of the whitener
% is [a(p+1) ... a(1)]/sigma placed at columns k-p..k. If this ever breaks, the
% filter argument order has been swapped (numerator/denominator).
    a = [1 -1.6 0.9 -0.25];  s = 0.7;  n = 20;  p = 3;
    r = s^2 * lsqvce.kernel.armaACV(a, 1, n);

    z = zeros(n, 0);
    [~, ~, F] = lsqvce.kernel.armaWhiten(a, 1, s, r, z, z, eye(n), true);

    expected = fliplr(a) / s;                    % [a3 a2 a1 1]/s
    for k = (p+1):n
        % tol: steady-state rows are the AR stencil EXACTLY; 1e-12 = house
        % factorisation round-trip bound (HANDOVER 6.1.2)
        tc.verifyEqual(F(k, k-p:k), expected, 'AbsTol', 1e-12, ...
            sprintf('row %d is not the AR whitening stencil', k));
    end
end

function testDecorrFactorIsExact(tc)
% categories: precision
    r = lsqvce.kernel.armaACV([1 -0.9], 1, 120);
    C = toeplitz(r);
    W = lsqvce.kernel.decorrFactor(C);
    % tol: exact identity through a solve with cond(C) ~ 20 for AR(1) phi=0.9;
    % house exact-identity bound 1e-10 rel (HANDOVER 6.1.2)
    tc.verifyEqual(W.'*W, C \ eye(120), 'RelTol', 1e-10, 'AbsTol', 1e-9);
    tc.verifyEqual(W*C*W.', eye(120), 'AbsTol', 1e-9);   % tol: same source
end

function testDecorrFactorRescuePathIsUsable(tc)
% categories: contract, precision
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

    W = tc.verifyWarning(@() lsqvce.kernel.decorrFactor(C), 'lsqvce:decorrFactor:notPD');

    % reconstruct the floored covariance and check W is its inverse factor
    [V2, d2] = eig((C+C.')/2, 'vector');
    d2 = max(real(d2), 1e-10*max(real(d2)));
    Cf = real(V2)*diag(d2)*real(V2).';
    % tol: the floored covariance is reconstructed through an eig at the
    % 1e-10 eigenvalue floor; 1e-6 = floor * cond headroom, measured stable
    tc.verifyEqual(W*Cf*W.', eye(n), 'AbsTol', 1e-6);
end

%% ================================================== signal processing
function testLevinsonExactOnKnownAR(tc)
% categories: precision
    a0 = [1 -1.6 0.9 -0.25];  s0 = 0.7;
    r  = s0^2 * lsqvce.kernel.armaACV(a0, 1, 8);
    [a, s2] = lsqvce.kernel.levinsonDurbin(r, 3);
    % tol: Levinson on the EXACT ACV of the generator is algebraically exact;
    % 1e-11 = 1e4*eps for the recursion's error growth over p=3 steps
    tc.verifyEqual(a, a0, 'AbsTol', 1e-11);
    tc.verifyEqual(sqrt(s2), s0, 'RelTol', 1e-11);   % tol: same source
end

function testLevinsonReflectionCoefficients(tc)
% categories: precision
% |kappa| < 1 for every order iff the autocovariance is positive definite.
    r = lsqvce.kernel.armaACV([1 -1.6 0.9 -0.25], 1, 21);
    [~, ~, refl] = lsqvce.kernel.levinsonDurbin(r, 20);
    tc.verifyLessThan(max(abs(refl)), 1);
end

function testARMAACVAgainstClosedForm(tc)
% categories: precision
    phi = 0.9;  th = 0.5;
    r  = lsqvce.kernel.armaACV([1 -phi], [1 th], 4);
    r0 = (1 + 2*phi*th + th^2)/(1 - phi^2);
    r1 = phi*r0 + th;
    % tol: closed-form ARMA(1,1) ACV (textbook); FFT truncation of the
    % impulse response is cut at eps relative -- 1e-9 is documented in the
    % armaACV example (help text claims < 1e-8)
    tc.verifyEqual(r(1), r0, 'RelTol', 1e-9);
    tc.verifyEqual(r(2), r1, 'RelTol', 1e-9);   % tol: same source
    tc.verifyEqual(r(3), phi*r1, 'RelTol', 1e-9);   % tol: same; r_k = phi*r_{k-1}, k > q
end

function testPSD2ACVSatisfiesParsevalAtZeroLag(tc)
% categories: precision
% r(0) = trapz(f, S) exactly -- an identity, not an approximation.
    dt = 5;
    f  = linspace(0, 0.5/dt, 8192).';
    H  = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * [1 -1.6 0.9 -0.25].');
    S  = 0.7^2 * dt * 2 * abs(H).^2;

    r = lsqvce.kernel.psd2acv(f, S, dt, 4);
    % tol: r(0) = trapz(f,S) is the quadrature IDENTITY (same sum, same order);
    % 1e-12 = house factorisation round-trip bound
    tc.verifyEqual(r(1), trapz(f, S), 'RelTol', 1e-12);
end

function testPSDToACVRoundTrip(tc)
% categories: precision, reference
% PSD of a known AR -> ACV must match the analytic ACV of that AR.
    dt = 5;  a = [1 -1.6 0.9 -0.25];  s = 0.7;
    f  = linspace(0, 0.5/dt, 8192).';
    H  = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * a(:));
    S  = s^2 * dt * 2 * abs(H).^2;

    r_num = lsqvce.kernel.psd2acv(f, S, dt, 64);
    r_ana = s^2 * lsqvce.kernel.armaACV(a, 1, 64);
    tc.verifyLessThan(norm(r_num - r_ana)/norm(r_ana), 2e-3);
end

function testPSD2ACVHasNoLagPeriod(tc)
% categories: precision, reference
% R7-01 (WP-111 part 2): a COARSE law (60 bins, the joint estimator's Grid)
% must give its own ACV at every lag of a long arc. Trapezoidal quadrature
% of S cos on nf points is periodic in the lag with period 2(nf-1): on the
% uniform grid covModel used to resample to, r(m-1) came back as r(1) (even
% m) or r(0) (odd m) -- a circulant covariance (Toeplitz cond 3.3e8, measured
% on the runner). The exact integral of the piecewise-linear S has no period.
    L = 60;  phi = 0.5;  dt = 5;
    fN = (0:L-1).' / (2 * (L - 1));
    S  = 2 * dt ./ abs(1 - phi * exp(-2i * pi * fN)).^2;      % one-sided, physical
    f  = fN / dt;
    for m = [255 256 2111]
        r = lsqvce.kernel.psd2acv(f, S, dt, m);
        % tol: the AR(1) correlation at lag m-1 >= 254 is 0.5^254 ~ 0; what the
        % linear interpolant adds is O(1e-8) (measured 1e-8 .. 2e-8); a
        % circulant gives 0.5 or 1
        tc.verifyLessThan(abs(r(m) / r(1)), 1e-6, sprintf('no wrap at m = %d', m));
        % tol: the interpolant of a 60-bin AR(1) law reproduces lag 1 to 1e-3
        tc.verifyEqual(r(2) / r(1), phi, 'AbsTol', 1e-3, sprintf('lag 1 at m = %d', m));
        % tol: cond of the exact law's Toeplitz is ~9 at every size (AR(1) 0.5)
        tc.verifyLessThan(cond(toeplitz(r)), 20, sprintf('well conditioned at m = %d', m));
    end
    % tol: the identity r(0) = trapz(f, S) still holds (the same linear interpolant)
    tc.verifyEqual(r(1), trapz(f, S), 'RelTol', 1e-12);
end

function testCovModelToeplitzLawIsNotCirculant(tc)
% categories: contract, precision
% R7-01 at the seam: covModel 'psd' 'toeplitz' no longer resamples the law
% onto floor(m/2)+1 points; its covariance is psd2acv's exact ACV, so an odd
% arc is not a circulant (the shipped 6.6.1 gave corr(1, m) = 1.000).
    L = 60;  phi = 0.5;  dt = 5;  m = 255;
    fN = (0:L-1).' / (2 * (L - 1));
    cm = struct('type', 'psd', 'fNorm', fN, 'S', 2 ./ abs(1 - phi * exp(-2i * pi * fN)).^2);
    Cfcn = lsqvce.engine.covModel(cm, (0:m-1).' * dt, struct('psdMode', 'toeplitz'));
    C = Cfcn(1:m, struct());
    % tol: as testPSD2ACVHasNoLagPeriod -- ~1e-8 interpolant, a circulant gives 1
    tc.verifyLessThan(abs(C(m, 1) / C(1, 1)), 1e-6, 'epoch 1 and epoch m are not neighbours');
    % tol: cond ~9 for the exact AR(1) 0.5 law, 3.3e8 for the circulant
    tc.verifyLessThan(cond(C), 20);
end

function testCheckPSDVerdicts(tc)
% categories: contract | tags: diagnostics
    dt = 5;  a = [1 -1.6 0.9 -0.25];
    f  = linspace(0, 0.5/dt, 4096).';
    H  = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * a(:));
    cm = struct('type','psd', 'f',f, 'S', 0.7^2*dt*2*abs(H).^2);

    chk = lsqvce.engine.checkPSD(cm, dt, struct('arOrder',20), false);
    tc.verifyEqual(chk.verdict, 'good');
    tc.verifyTrue(chk.stable);
    tc.verifyLessThan(chk.acvRelErr, 1e-3);

    % AR(1) is far too low an order for an AR(3) spectrum
    chkBad = tc.verifyWarning(@() lsqvce.engine.checkPSD(cm, dt, struct('arOrder',1), false), ...
                              'lsqvce:checkPSD:badFit');
    tc.verifyEqual(chkBad.verdict, 'bad');
end

%% ================================================== estimation exactness
function testPreeliminationEqualsJointSolve(tc)
% categories: contract
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
% categories: precision
    cfg  = tc.TestData.cfg;
    days = 1:2;

    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), days); end
    sol = acc.solve();

    [Af, lf] = assembleFull(cfg, days);
    xf = Af \ lf;

    for k = 1:numel(cfg.sats)
        y    = acc.backSubstitute(cfg.sats(k), sol.x);
        yref = xf(cfg.nc + (k-1)*cfg.nl + (1:cfg.nl));
        tc.verifyEqual(y, yref, 'RelTol', 1e-8);   % tol: dense-vs-batched same algebra; 1e-8 = whitened-path house bound
    end
end

function testLocalCovarianceEqualsTheJointSolve(tc)
% categories: reference, precision
% R7-40: the local block of the monolithic joint covariance s0^2 (A'A)^-1 is
% what localCovariance(sysID, Cxx, s0sq) must return -- with the SOLVE's s0sq.
% It used the store's s0sq_hint (l'l / dof with the commons unsolved), which
% put the verbose locals line ~300x too high (R7 mirror m6_schur_cyy.py).
    cfg  = tc.TestData.cfg;
    days = 1:2;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), days); end
    sol = acc.solve();
    [Af, lf] = assembleFull(cfg, days);
    xf  = Af \ lf;
    s0f = sum((lf - Af*xf).^2) / (size(Af,1) - size(Af,2));
    Cf  = s0f * inv(Af.' * Af);                          %#ok<MINV> the reference, dense on purpose
    for k = 1:numel(cfg.sats)
        iy  = cfg.nc + (k-1)*cfg.nl + (1:cfg.nl);
        Cyy = acc.localCovariance(cfg.sats(k), sol.Cxx, sol.s0sq);
        % tol: dense-vs-batched same algebra; 1e-8 = the whitened-path house bound
        tc.verifyEqual(Cyy, Cf(iy, iy), 'RelTol', 1e-8, ...
            'the local block of the joint covariance, with the solve''s s0^2');
    end
end

function testBatchOverlapGivesMachinePrecision(tc)
% categories: vectorisation
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
% categories: precision, robustness, vectorisation
% The answer must not depend on how the arc was chopped up.
    cfg = tc.TestData.cfg;
    xRef = [];
    for bs = [37 100 400 1e9]
        s = runAcc(cfg, 1:2, 'Backend','qr', 'BatchSize',bs, ...
                   'BatchOverlap',3, 'Verbose',false);
        if isempty(xRef)
            xRef = s.x;
        else
            % tol: exact-identity claim (overlap makes batching exact);
            % house exact-identity bound 1e-10 rel (HANDOVER 6.1.2)
            tc.verifyEqual(s.x, xRef, 'RelTol', 1e-10, ...
                sprintf('BatchSize = %g changed the answer', bs));
        end
    end
end

function testPSDARMatchesDenseToeplitz(tc)
% categories: precision, vectorisation
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

    % tol: AR(20) fit vs exact Toeplitz are DIFFERENT models of one PSD;
    % 1e-6 is the measured model-agreement level (v10review.py), not eps
    tc.verifyEqual(sAR.x, sTP.x, 'RelTol', 1e-6);
end

function testWhitenModeCholMatchesFilter(tc)
% categories: precision, vectorisation
% Forcing the dense path must reproduce the filter path.
    cfg = tc.TestData.cfg;
    sF = runAcc(cfg, 1:2, 'Backend','qr','BatchSize',400,'BatchOverlap',10, ...
                'UserOpt', struct('whitenMode','filter'), 'Verbose',false);
    sC = runAcc(cfg, 1:2, 'Backend','qr','BatchSize',400,'BatchOverlap',10, ...
                'UserOpt', struct('whitenMode','chol'), 'Verbose',false);
    % tol: filter vs chol whitening agree at the whitened-path level;
    % 1e-6 measured across seeds (v10review.py), boundary rows differ
    tc.verifyEqual(sF.x, sC.x, 'RelTol', 1e-6);
end

%% ================================================== conditioning
function testQRBeatsNormalOnIllConditioned(tc)
% categories: robustness
% The whole justification for the QR backend: cond(A) versus cond(A)^2.
%
% Measured on the BACKEND IN ISOLATION; the accumulator path is pinned by
% testAccumulatorQRStaysSquareRoot (R7 WP-112). The masking this test used to
% blame on regularisation was structural (review R7-03): until 6.7 the
% accumulator squared the factor before the Schur step.
    muteWarn(tc, 'lsqvce:normalEq:cancellation');   % tight-fit / ill-conditioned regime provoker (WP-49 unified id)
    rs = RandStream('mt19937ar','Seed',5);
    m  = 800;   n = 8;

    for k = [4 6 8]                                  % cond(A) = 1e4, 1e6, 1e8
        [U,~] = qr(randn(rs, m, n), 0);
        [V,~] = qr(randn(rs, n, n), 0);
        A = U * diag(logspace(0, -k, n)) * V.';

        xt = randn(rs, n, 1);
        l  = A*xt + 1e-10*randn(rs, m, 1);
        xr = A \ l;                                   % LAPACK reference (QR based)

        accN = lsqvce.engine.NormalEqSystem(n);   accN.add(A, l);
        accQ = lsqvce.engine.QRSystem(n);         accQ.add(A, l);

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

function testAccumulatorQRStaysSquareRoot(tc)
% categories: precision, robustness
% R7 WP-112 (review R7-03): the ACCUMULATOR path of the QR backend, not the
% backend in isolation. Commons of cond(A) = 1e7 .. 1e9, with and without
% local parameters, one system through BlockAccumulator: x_c must equal the
% dense [A B] \ l to a tiny fraction of its formal sd, and no full-rank
% problem may be refused. Until 6.7 the locals were Schur-eliminated in
% normal form and merge took chol of the reduced N: median 2-3 formal sd off
% at cond 1e7-1e9, and 4-11 of 20 full-rank problems refused at 1e8.5-1e9
% (mirror mirrors/wp112/square_root_accumulator.py; the square-root path
% measured <= 1.6e-4 sd, mostly 1e-8).
    m = 800;  nc = 8;
    for k = [7 8.5 9]
        for nl = [0 3]
            rs = RandStream('mt19937ar', 'Seed', round(10*k) + nl);
            for seed = 1:4
                [U, ~] = qr(randn(rs, m, nc), 0);  [V, ~] = qr(randn(rs, nc, nc), 0);
                A = U * diag(logspace(0, -k, nc)) * V.';
                B = randn(rs, m, nl);
                l = [A B] * randn(rs, nc + nl, 1) + 1e-4 * randn(rs, m, 1);
                [Q, Rd] = qr([A B], 0);
                xr = Rd \ (Q.' * l);
                v  = l - [A B] * xr;
                Ri = Rd \ eye(nc + nl);
                sd = sqrt(v.' * v / (m - nc - nl)) * sqrt(sum(Ri.^2, 2));
                sol = oneSystemAcc(A, B, l);
                err = max(abs(sol.x - xr(1:nc)) ./ sd(1:nc));
                tc.verifyLessThan(err, 1e-3, sprintf(['cond 1e%g, nl = %d, seed %d: x_c is %.2e ' ...
                    'formal sd off the dense solve'], k, nl, seed, err));
            end
        end
    end
end

function testFreeNetworkRefusalIsDeterministic(tc)
% categories: contract, robustness
% R7 WP-112 (review R7-16): the rank decision is QRSystem.RankTol on the
% pivots, not LAPACK's chol. Random connected levelling networks: WITHOUT a
% datum every one is refused (lsqvce:qr:singularR), WITH one fixed height
% every one is solved. The review measured the shipped path at 61 % refused /
% 39 % solved with an arbitrary datum (sd ~1e7) on 2000 such networks; the
% mirror's square-root path refuses 2000 / 2000 and solves 2000 / 2000.
    rs = RandStream('mt19937ar', 'Seed', 112);
    for trial = 1:40
        u = randi(rs, [4 11]);
        pairs = [(1:u).', [2:u, 1].'];
        for c = 1:randi(rs, [1 u-1])
            pairs(end+1, :) = randperm(rs, u, 2);   %#ok<AGROW>
        end
        A = zeros(size(pairs, 1), u);
        for r = 1:size(pairs, 1), A(r, pairs(r, 1)) = -1;  A(r, pairs(r, 2)) = 1; end
        A = A * (0.5 + 2.5 * rand(rs));
        l = A * (10 * randn(rs, u, 1)) + 0.01 * randn(rs, size(A, 1), 1);
        id = '';
        try
            lsqvce.fit(struct('A', A, 'l', l), Verbose=false);
        catch err
            id = err.identifier;
        end
        tc.verifyEqual(id, 'lsqvce:qr:singularR', sprintf('free network %d (u = %d) must be refused', trial, u));
        Ad = [A; [1, zeros(1, u-1)]];  ld = [l; 0];
        R = lsqvce.fit(struct('A', Ad, 'l', ld), Verbose=false);
        tc.verifyTrue(all(isfinite(R.x)), sprintf('network %d with a datum must be solved', trial));
    end
end

function testNormalEquationsLoseTheAnswerAtHighCondition(tc)
% categories: contract
% The failure mode the QR backend exists to prevent. At cond(A) = 1e8 the normal
% matrix has cond(A)^2 = 1e16 > 1/eps: the information is destroyed by forming
% A'*A and NO solver can recover it afterwards.
    muteWarn(tc, 'lsqvce:normalEq:cancellation');   % tight-fit / ill-conditioned regime provoker (WP-49 unified id)
    rs = RandStream('mt19937ar','Seed',11);
    m  = 800;   n = 8;

    [U,~] = qr(randn(rs, m, n), 0);
    [V,~] = qr(randn(rs, n, n), 0);
    A = U * diag(logspace(0, -8, n)) * V.';
    l = A*randn(rs, n, 1) + 1e-10*randn(rs, m, 1);
    xr = A \ l;

    accN = lsqvce.engine.NormalEqSystem(n);   accN.add(A, l);
    accQ = lsqvce.engine.QRSystem(n);         accQ.add(A, l);

    eN = norm(accN.solve() - xr)/norm(xr);
    eQ = norm(accQ.solve() - xr)/norm(xr);

    tc.verifyGreaterThan(eN, 1e-3, 'the normal equations should fail badly here');
    tc.verifyLessThan(eQ, 1e-7,   'QR should still be accurate here');
    tc.verifyGreaterThan(eN/eQ, 1e4);
end

function testRegularisationMasksBackendChoice(tc)
% categories: precision, robustness
% Worth knowing, and easy to misread: once the local block is singular enough
% that TSVD engages, BOTH backends give the same answer -- the regularisation,
% not the factorisation, sets the accuracy. Do not conclude from a full-pipeline
% test that the backend does not matter; test the backend in isolation.
    muteWarn(tc, 'lsqvce:normalEq:cancellation');   % tight-fit / ill-conditioned regime provoker (WP-49 unified id)
    cfg = tc.TestData.cfg;
    cfg.collinear = true;   cfg.collinearEps = 1e-9;

    sN = runAcc(cfg, 1:2, 'Backend','normal', 'Reg','tsvd', 'Verbose',false);
    sQ = runAcc(cfg, 1:2, 'Backend','qr',     'Reg','tsvd', 'Verbose',false);

    % tol: TSVD-regularised backends at collinearEps 1e-9; agreement is set
    % by the truncation, 1e-6 measured stable across seeds
    tc.verifyEqual(sQ.x, sN.x, 'RelTol', 1e-6, ...
        'with TSVD engaged the two backends should agree');
end

function testBackendsAgreeWhenWellConditioned(tc)
% categories: precision
    cfg = tc.TestData.cfg;
    sN = runAcc(cfg, 1:3, 'Backend','normal', 'Verbose',false);
    sQ = runAcc(cfg, 1:3, 'Backend','qr',     'Verbose',false);
    % tol: well-conditioned problem, QR vs normal lose ~cond(A)^2*eps ~ 1e-10;
    % 1e-8/1e-7 add headroom for the release-to-release LAPACK drift
    tc.verifyEqual(sQ.x,      sN.x,      'RelTol', 1e-8);
    tc.verifyEqual(sQ.sigma0, sN.sigma0, 'RelTol', 1e-7);   % tol: same source
end

%% ================================================== statistical fidelity
function testSigma0IsUnbiased(tc)
% categories: precision
% With a correctly specified covariance, sigma0 -> 1. This is the end-to-end
% check that the whitening, the batching and the ltl bookkeeping all agree.
    cfg = tc.TestData.cfg;
    sol = runAcc(cfg, 1:6, 'Backend','qr', 'Verbose',false);
    % tol: sigma0^2 ~ chi2_r/r with r ~ 2e4 => sd(sigma0) ~ 1/sqrt(2r) ~ 5e-3;
    % 0.05 = 10 sigma, a bright-line unbiasedness check, not containment
    tc.verifyEqual(sol.sigma0, 1, 'AbsTol', 0.05);
end

function testResidualPassReproducesSigma0(tc)
% categories: precision
% Two independent routes to sigma0 must agree to machine precision.
    cfg = tc.TestData.cfg;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:3); end
    sol = acc.solve();
    res = acc.residuals(sol, 'Verbose', false);
    % tol: the residual pass recomputes the SAME quadratic form; house
    % exact-identity bound 1e-10 rel (HANDOVER 6.1.2)
    tc.verifyEqual(res.sigma0, sol.sigma0, 'RelTol', 1e-10);
end

function testWhitenedResidualsAreWhite(tc)
% categories: precision
% The stochastic model is right iff the whitened residuals are white with unit
% variance. This tests the noise COLOUR, not just its level.
    cfg = tc.TestData.cfg;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:6); end
    sol = acc.solve();
    res = acc.residuals(sol, 'Verbose', false);

    w = vertcat(res.sys(1).obs.w);
    n = numel(w);

    % tol: sd(std) ~ 1/sqrt(2m) over the pooled residuals; 0.08 ~ 4/sqrt(Z)
    % Monte-Carlo bar (HANDOVER 6.1.2), a 2.93-sigma test, not containment
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
% categories: precision
    cfg = tc.TestData.cfg;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:3); end
    sol = acc.solve();

    % tol: symmetry is enforced by construction; 1e-14 ~ 100*eps guards the
    % accumulation order only
    tc.verifyEqual(sol.Cxx, sol.Cxx.', 'AbsTol', 1e-14);
    tc.verifyGreaterThan(min(eig(sol.Cxx)), 0);

    Cyy = acc.localCovariance(cfg.sats(1), sol.Cxx, sol.s0sq);
    tc.verifyEqual(Cyy, Cyy.', 'AbsTol', 1e-14);   % tol: same source
    tc.verifyGreaterThan(min(eig(Cyy)), 0);
end

function testFormalErrorsMatchEmpiricalScatter(tc)
% categories: precision, robustness, vectorisation
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
    % tol: formal-vs-empirical scatter ratio over Z repeats; 0.10 ~ 4/sqrt(Z)
    % Monte-Carlo bar (HANDOVER 6.1.2)
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
    f = @(d,ot,a,uo) lsqvce.engine.synthLoader(cfg, s, d, ot, a, uo);
end

function sol = oneSystemAcc(A, B, l)
%ONESYSTEMACC  One system, one day, one white arc through BlockAccumulator
%   (Backend "qr"): the accumulator path R7 WP-112 is about.
    nc = size(A, 2);  nl = size(B, 2);
    build = @(d, ot) struct('A', {{A}}, 'B', {{B}}, 'l', {{l}}, 't', {{(1:numel(l)).'}});
    L = lsqvce.engine.makeLoader(build, 'ObservationTypes', {'x'}, ...
            'Covariance', struct('x', struct('type', 'diag', 'sigma', 1)), ...
            'NumLocal', nl, 'NumCommon', nc);
    acc = lsqvce.engine.BlockAccumulator(nc, 'Backend', 'qr', 'Verbose', false);
    acc.addSystem(1, nl, L.closure, 1);
    sol = acc.solve();
end

function sol = runAcc(cfg, days, varargin)
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, varargin{:});
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
                        W = lsqvce.kernel.decorrFactor(arc.Cfcn(1:m, uo));
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

%% ===================================================== vectorisation
function testArmaWhitenBatchedEqualsColumnAtATime(tc)
% categories: vectorisation
% WP-43: the single cheapest reference test in a vectorised language, and the
% one most often assumed -- whitening an n x k block in ONE call must equal
% whitening the k columns one at a time. Both paths run filter() per column
% and the same triangular startup solve, so the expectation is BITWISE for
% the filter part; asserted at 10*eps relative only to stay robust to
% multi-RHS BLAS kernels in the startup solve, which may order the flops
% differently per column.                                    % tol: 10*eps
    a = [1 -1.6 0.9 -0.25];   s = 0.7;   n = 300;   k = 7;
    r = s^2 * lsqvce.kernel.armaACV(a, 1, n);
    z = zeros(n, 0);
    rs = RandStream('mt19937ar', 'Seed', 11);
    X  = randn(rs, n, k);

    [~, ~, Xw] = lsqvce.kernel.armaWhiten(a, 1, s, r, z, z, X, true);

    Xcol = zeros(n, k);
    for j = 1:k
        [~, ~, Xcol(:, j)] = lsqvce.kernel.armaWhiten(a, 1, s, r, z, z, X(:, j), true);
    end

    % tol: 10*eps -- see the header of this test
    tc.verifyEqual(Xw, Xcol, 'RelTol', 10*eps, ...
        'batched whitening must equal column-at-a-time whitening');
end

%% ================================================== WP-49: honest vPv
function testNormalBackendPoisonedVtvIsNaN(tc)
% categories: contract, robustness | tags: diagnostics
% WP-49: a cancellation-destroyed vtv must leave the normal backend as NaN,
% never as a plausible-looking number (the flat path's measured argument:
% 89% wrong at sigma = 1e-6, factor 1e6 at 1e-9). The tight-fit recipe
% drives l'Pl / vtv far past the 8-digit guard deterministically; x itself
% is NOT poisoned and must still match the LAPACK reference.
    rs = RandStream('mt19937ar', 'Seed', 4901);
    m = 600;  n = 6;
    A  = randn(rs, m, n) + 20;                 % large offset -> huge l'Pl
    xt = randn(rs, n, 1);
    l  = A*xt + 1e-9*randn(rs, m, 1);          % ~zero residual

    accN = lsqvce.engine.NormalEqSystem(n);  accN.add(A, l);
    tc.verifyWarning(@() accN.solve(), 'lsqvce:normalEq:cancellation');

    muteWarn(tc, 'lsqvce:normalEq:cancellation');  % second call below re-fires
    [x, Cxx, s0sq, dof, vtv] = accN.solve();
    tc.verifyTrue(isnan(vtv)  && isnan(s0sq), ...
        'a poisoned vtv/s0sq must be NaN, not a number');
    tc.verifyTrue(all(isnan(Cxx(:))), ...
        'a covariance scaled by a destroyed variance factor is unknown');
    tc.verifyGreaterThan(dof, 0);
    % tol: x is solved BEFORE the vtv subtraction and is unaffected; benign
    % conditioning here, so backslash-reference agreement at 1e-8
    tc.verifyEqual(x, A \ l, 'RelTol', 1e-8);
end

function testResidualPassSupersedesPoisonedSigma0(tc)
% categories: contract, robustness | tags: diagnostics, provenance
% WP-49, the promise end to end on the BLOCK path: normal backend +
% tight fit -> sol.sigma0 = NaN (refused, not wrong); the residual pass
% recomputes it honestly, says so with its own identifier, flags
% res.supersedesSol, and its number agrees with the QR backend's stable
% route. The QR run on the same data is the second, independent route.
    cfg = tc.TestData.cfg;
    % signal/noise chosen so the NORMAL route is destroyed (l'Pl/vtv ~ 1e11,
    % far past the 8-digit guard) while QR's dynamic range (~3e5) stays well
    % inside double precision -- QR must remain the valid second route. The
    % first probe at 3e4/1e-5 broke QR TOO (sigma0 = 0; observation recorded
    % in 10.11); this is the measured sweet spot.
    cfg.x_c = 30 * (1:cfg.nc).';
    cfg.cov.pos.sigma = 1e-4;
    cfg.cov.kbr.sigma = 1e-4;

    accN = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','normal', ...
                                   'KeepLocal',true, 'Verbose',false);
    for s = cfg.sats
        accN.addSystem(s, cfg.nl, binder(cfg,s), 1:2);
    end
    solN = verifyWarning(tc, @() accN.solve(), 'lsqvce:normalEq:cancellation');
    tc.verifyTrue(isnan(solN.sigma0), ...
        'the block sol must carry the refused NaN, not a wrong number');

    res = verifyWarning(tc, @() accN.residuals(solN, 'verbose', false), ...
                        'lsqvce:residuals:supersedesSol');
    tc.verifyTrue(res.supersedesSol);
    tc.verifyTrue(isfinite(res.sigma0) && res.sigma0 > 0);

    solQ = runAcc(cfg, 1:2, 'Backend','qr', 'Verbose',false);
    % tol: MEASURED 8.0e-5 here (3.2e-3 on a harder 1-sat probe), 2026-08-20:
    % the pass inherits the normal-path x-hat (2nd-order penalty); well
    % below sigma0's sampling scatter 1/sqrt(2*dof) ~ 1.4e-2
    tc.verifyEqual(res.sigma0, solQ.sigma0, 'RelTol', 1e-3, ...
        'the superseding number must match the stable backend');
end

%% ================================================== WP-50: honest vtv, QR
function testQrMergeFoldKeepsTheResidualMass(tc)
% categories: contract, robustness | tags: diagnostics
% WP-50 asked the block-QR path to REFUSE (NaN) when the normal-form fold of a
% pre-eliminated contribution destroyed the residual mass (rho2 = ltl_red -
% |rhs|^2; at the WP-49 deep-regime probe it had produced sigma0 = EXACTLY 0).
% R7 WP-112 (R7-03) removed the subtraction: the contribution is merged as
% square-root rows, the residual row rides along, so at the SAME regime the
% accumulated sigma0 is now honest -- no warning at the fold, finite, ~1, and
% equal to the residual pass's.
    cfg = tc.TestData.cfg;
    cfg.x_c = 3e4 * (1:cfg.nc).';              % the exact observed regime
    cfg.cov.pos.sigma = 1e-5;
    cfg.cov.kbr.sigma = 1e-5;

    accQ = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', ...
                                   'KeepLocal',true, 'Verbose',false);
    tc.verifyWarningFree(@() addAllSystems(accQ, cfg, 1:2));
    solQ = tc.verifyWarningFree(@() accQ.solve());
    % tol: the declared law is exact, so sigma0 ~ 1 with sampling scatter
    % 1/sqrt(2*dof) ~ 1.4e-2 -- band 5 sigma (the WP-49 residual pass measured 1.0072)
    tc.verifyEqual(solQ.sigma0, 1, 'AbsTol', 0.07, ...
        'the square-root fold keeps the residual mass: an honest sigma0');
    % the residual pass recomputes l - A x - B y at a 1e9 dynamic range and may
    % differ from the square root past its 1e-6 consistency bar
    muteWarn(tc, 'lsqvce:residuals:sigma0Mismatch');
    res = accQ.residuals(solQ, 'verbose', false);
    % tol: 1e-3 relative, below sigma0's sampling scatter 1/sqrt(2*dof) ~ 1.4e-2
    tc.verifyEqual(res.sigma0, solQ.sigma0, 'RelTol', 1e-3);
end
function testQrHealthyPathUnchangedByTheGuards(tc)
% categories: contract
% Regression lock for WP-50: on a healthy fixture the guards must be
% invisible -- finite vtv, sigma0 agreeing with the residual pass at the
% established 1e-10, and NO cancellation warning (asserted warning-free
% so a guard that misfires on healthy data cannot land silently).
    cfg = tc.TestData.cfg;
    accQ = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', ...
                                   'KeepLocal',true, 'Verbose',false);
    for s = cfg.sats
        accQ.addSystem(s, cfg.nl, binder(cfg,s), 1:2);
    end
    solQ = tc.verifyWarningFree(@() accQ.solve());
    tc.verifyTrue(isfinite(solQ.sigma0) && solQ.sigma0 > 0);
    res = accQ.residuals(solQ, 'verbose', false);
    tc.verifyFalse(res.supersedesSol);
    % tol: the established healthy-route agreement (regression-locked)
    tc.verifyEqual(res.sigma0, solQ.sigma0, 'RelTol', 1e-10);
end

function addAllSystems(acc, cfg, days)
%ADDALLSYSTEMS  Plain-loop addSystem over cfg.sats (WP-50): arrayfun over a
%   void method errors on output collection, which broke the r1 warning
%   assertion and leaked the id into the inventory.
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg, s), days);
    end
end
