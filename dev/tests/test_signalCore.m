function tests = test_signalCore
%TEST_SIGNALCORE  Direct tests for three functions that previously ran only INSIDE
%   other tests, with no assertion on their own behaviour (roadmap item 2):
%
%     lsqvce.engine.outliers   -- data snooping w_i = |v_i|/sqrt(r_i) flags blunders
%     lsqvce.kernel.varWhiten  -- vector-AR whitener == dense block-Cholesky GLS
%     lsqvce.kernel.ar2arma    -- Durbin reduction, ACV round-trip of an ARMA(p,q)
%
%   Reference properties and tolerances established first in
%   validation/v8verify.py (16/16). These tests assert the same properties
%   through the real MATLAB signatures.
%
%   Run:
%     results = runtests('tests/test_signalCore');  disp(table(results))
%
%   See also TEST_V7, LSQVCE.ENGINE.OUTLIERS, LSQVCE.KERNEL.VARWHITEN, LSQVCE.KERNEL.AR2ARMA
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-134 (third review): covModel takes a VAR(1), varACV a VAR(0); ar2arma normalises aHi(1).
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    tc.TestData.cfg = makeCfg();
end

%% ===================================================== lsqvce.engine.outliers
function testOutliersFlagsInjectedBlunder(tc)
% categories: contract
% Wrap the synthetic loader so ONE arc carries a large spike, then confirm
% lsqvce.engine.outliers flags an observation in that arc and reports the correct
% redundancy identity. The blunder is +15 sigma at a known epoch.
    cfg  = tc.TestData.cfg;
    days = 1:3;
    tgt  = struct('sat', cfg.sats(1), 'obsType', 'pos', 'day', 2, 'arc', 1, ...
                  'row', 50, 'size', 15);

    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false, ...
                                  'KeepLocal',true);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, spikedBinder(cfg, s, tgt), days);
    end
    sol = acc.solve();

    out = lsqvce.engine.outliers(acc, sol, 'Threshold', 4.0, 'Verbose', false);

    tc.verifyGreaterThanOrEqual(out.nFlagged, 1);

    % at least one flag must sit in the target arc
    hit = false;
    for j = 1:numel(out.flagged)
        f = out.flagged(j);
        if f.sysID == tgt.sat && strcmp(char(f.obsType), tgt.obsType) && ...
           f.day == tgt.day && f.arc == tgt.arc
            hit = true;  break
        end
    end
    tc.verifyTrue(hit, 'the injected blunder was not flagged in its own arc');
end

function testOutliersCleanDataNoFalseStorm(tc)
% categories: reference
% On clean data the flag count must be a small fraction of n (Baarda at
% threshold 4.0), and the redundancy identity must hold exactly.
    cfg  = tc.TestData.cfg;
    days = 1:3;
    acc  = buildCleanAcc(cfg, days);
    sol  = acc.solve();
    out  = lsqvce.engine.outliers(acc, sol, 'Threshold', 4.0, 'Verbose', false);

    red  = lsqvce.engine.redundancy(acc, sol, 'Verbose', false);
    tc.verifyLessThan(red.check, 1e-6 * max(red.rTotal, 1));   % sum(r)=m-u
    tc.verifyLessThan(out.nFlagged, 0.05 * red.rTotal);        % no false storm
end

%% ===================================================== lsqvce.kernel.varWhiten
function testVarWhitenMatchesDenseGLS(tc)
% categories: precision, robustness, reference
% The vector-AR FIR whitener with the exact stationary startup must reproduce
% the dense block-Cholesky GLS solve, row for row (docstring: rel diff ~1e-10;
% v8verify.py: 2e-16). Everything is derived from the toolbox's own varACV so
% the test is convention-independent.
    d  = 2;  nE = 40;  u = 3;
    A1 = cat(3, [0.4 0.1; -0.2 0.5]);          % VAR(1), stable
    Lsig  = [1.0 0.0; 0.3 0.8];
    Sigma = Lsig * Lsig.';

    G = lsqvce.kernel.varACV(A1, Sigma, nE);          % Gamma(0..nE), toolbox convention
    C = blockToeplitzCov(G, d, nE);            % dense m x m covariance
    m = d * nE;

    rs = RandStream('mt19937ar', 'Seed', 88);
    Adesign = randn(rs, m, u);
    xTrue   = randn(rs, u, 1);
    Lc      = chol(C, 'lower');
    l       = Adesign * xTrue + Lc * randn(rs, m, 1);

    % dense GLS reference
    Cinv    = C \ eye(m);
    xDense  = (Adesign.' * Cinv * Adesign) \ (Adesign.' * Cinv * l);
    quadRef = l.' * Cinv * l;

    % whiten via varWhiten (B empty), first batch -> exact startup
    Bempty  = zeros(m, 0);
    [Aw, ~, lw] = lsqvce.kernel.varWhiten(A1, Sigma, G, d, Adesign, Bempty, l, true);
    xWhite  = (Aw.' * Aw) \ (Aw.' * lw);

    % tol: v8verify.py measures 2e-16 rel; 1e-8 leaves 8 decades for BLAS drift
    tc.verifyEqual(xWhite, xDense, 'RelTol', 1e-8);
    tc.verifyEqual(lw.' * lw, quadRef, 'RelTol', 1e-8);   % tol: same source
end

%% ===================================================== lsqvce.kernel.ar2arma
function testAr2armaRoundTrip(tc)
% categories: precision
% A resonant ARMA(2,1) -- a sharp spectral line, the 1/rev regime ar2arma is
% for. Build its exact ACV, fit a long AR, reduce to ARMA(2,1), and require the
% TRUE polynomials back, with no fallback. WP-126: the Hankel fit ran the
% recursion backwards (reciprocal roots, 0.94 -> 1.064), so the stationarity
% guard fired here and the test passed on the fallback AR(60) -- the muted
% 'unstable' warning was the defect. Reference: Durbin's method on the exact
% impulse response; dev/validation/wp126_ar2arma_mirror.py: 8e-16 on a,
% 2e-14 on b, ACV relErr 2e-14.
    rad = 0.94;  ang = 2*pi*0.12;
    aTrue = real(poly([rad*exp(1i*ang), rad*exp(-1i*ang)]));   % [1 a1 a2]
    bTrue = [1.0 0.5];
    nACV  = 256;
    rTarget = lsqvce.kernel.armaACV(aTrue, bTrue, nACV);              % s2 = 1

    [aHi, s2Hi] = lsqvce.kernel.levinsonDurbin(rTarget, 60);
    opt = struct('armaP', 2, 'armaN', 4096);
    [a, b, s2] = verifyWarningFree(tc, @() lsqvce.kernel.ar2arma(aHi, s2Hi, 1, opt));

    % tol: 1e-9 -- the AR(60) truncates the ARMA's 1/B(z) tail at 0.5^60
    % (1e-18), so what remains is the Levinson / least-squares rounding the
    % mirror measures at 1e-14; four orders of headroom across platforms
    tc.verifyEqual(a, aTrue, 'AbsTol', 1e-9);
    tc.verifyEqual(b, bTrue, 'AbsTol', 1e-9);               % tol: same source
    rArma  = s2 * lsqvce.kernel.armaACV(a, b, nACV);
    relErr = norm(rArma - rTarget) / norm(rTarget);
    tc.verifyLessThan(relErr, 1e-9);                         % tol: same source
end

function testAr2armaNonInvertibleMaFallsBackToTheAr(tc)
% categories: robustness, contract
% WP-126: an ARMA(1,2) whose MA part is reduced to q = 1 is fitted with a
% NON-INVERTIBLE B(z) (mirror: b = [1 1.160], root 1.16; the AR part 0.940 is
% stationary, so only the new guard catches it). The whitener A(z)/B(z)
% would be unstable, and COVMODEL refuses such a law -- the reduction must
% fall back to the input AR with lsqvce:ar2arma:nonInvertibleMA, not fail.
    rTarget = lsqvce.kernel.armaACV([1 -0.9], [1 1.2 0.81], 256);
    [aHi, s2Hi] = lsqvce.kernel.levinsonDurbin(rTarget, 60);
    opt = struct('armaP', 1, 'armaN', 1024);
    tc.verifyWarning(@() lsqvce.kernel.ar2arma(aHi, s2Hi, 1, opt), 'lsqvce:ar2arma:nonInvertibleMA');
    muteWarn(tc, 'lsqvce:ar2arma:nonInvertibleMA');   % provoked on purpose here (gate rule 7)
    [a, b, s2] = lsqvce.kernel.ar2arma(aHi, s2Hi, 1, opt);
    tc.verifyEqual(a, aHi(:).', 'the fallback is the input AR');
    tc.verifyEqual(b, 1);
    tc.verifyEqual(s2, s2Hi);
end

function testAr2armaSmoothFloorStaysStable(tc)
% categories: precision, robustness
% On a smooth coloured floor (a plain AR(2), no lines) an ARMA(6,1) is
% over-parameterised: the Hankel system is rank-deficient and the fit may
% come back non-stationary or non-invertible. Whatever comes back must be
% usable as a whitener: normalised, a stationary AR part and an invertible MA
% part -- the fitted model, or the input AR after a warning.
    muteWarn(tc, 'lsqvce:ar2arma:unstable');          % provoked on purpose here (WP-47, gate rule 7)
    muteWarn(tc, 'lsqvce:ar2arma:nonInvertibleMA');   % provoked on purpose here (WP-126)
    muteWarn(tc, 'MATLAB:rankDeficientMatrix');       % provoked on purpose here (WP-47, gate rule 7)
    aAR = [1 -0.6 0.1];                        % smooth AR(2)
    rTarget = lsqvce.kernel.armaACV(aAR, 1, 128);
    [aHi, s2Hi] = lsqvce.kernel.levinsonDurbin(rTarget, 40);
    opt = struct('armaP', 6, 'armaN', 2048);
    [a, b, ~] = lsqvce.kernel.ar2arma(aHi, s2Hi, 1, opt);
    % tol: normalisation exact by construction; house round-trip bound 1e-12
    tc.verifyEqual(a(1), 1, 'AbsTol', 1e-12);
    tc.verifyEqual(b(1), 1, 'AbsTol', 1e-12);               % tol: same source
    tc.verifyLessThan(max(abs(roots(a))), 1);
    if numel(b) > 1, tc.verifyLessThan(max(abs(roots(b))), 1); end
end

%% =================================================================== contract
function testSignalKernelContractRejections(tc)
% categories: contract
% WP-43 (audit check 16): ar2arma, varWhiten and whittle validate their
% arguments with documented ids, yet none of those rejections was exercised
% by a CI-run suite -- varWhiten/whittle coverage lived ONLY in the CI-skipped
% test_parallel. Rejected calls run before the body, so each check here costs
% microseconds; contract may not be exempted for a validating function.
    d = 2;
    A1    = cat(3, [0.4 0.1; -0.2 0.5]);
    Sigma = eye(d);
    G     = lsqvce.kernel.varACV(A1, Sigma, 8);

    % ar2arma: armaN too short for the requested ARMA(p,q) (< 4*(p+q+1))
    tc.verifyError(@() lsqvce.kernel.ar2arma([1 -0.5], 1.0, 1, ...
                       struct('armaP', 6, 'armaN', 8)), ...
                   'lsqvce:ar2arma:tooShort');

    % varWhiten: wrong coefficient shape / wrong Sigma size / partial epoch
    m = d * 8;  X = zeros(m, 1);  Z = zeros(m, 0);
    tc.verifyError(@() lsqvce.kernel.varWhiten(ones(3,2,1), Sigma, G, d, X, Z, X, true), ...
                   'lsqvce:varWhiten:badA');
    tc.verifyError(@() lsqvce.kernel.varWhiten(A1, eye(3), G, d, X, Z, X, true), ...
                   'lsqvce:varWhiten:badSigma');       % wrong size (3x3, d=2)
    tc.verifyError(@() lsqvce.kernel.varWhiten(A1, [1 2; 2 1], G, d, X, Z, X, true), ...
                   'lsqvce:varWhiten:badSigma');       % right size, not PD
    Xodd = zeros(m + 1, 1);
    tc.verifyError(@() lsqvce.kernel.varWhiten(A1, Sigma, G, d, Xodd, zeros(m+1,0), ...
                                        Xodd, true), ...
                   'lsqvce:varWhiten:notWholeEpochs');

    % whittle: malformed cross-ACV argument
    tc.verifyError(@() lsqvce.kernel.whittle(ones(2,3,4), 1), 'lsqvce:whittle:badG');
end

%% =================================================================== helpers
function C = blockToeplitzCov(G, d, nE)
% Dense block-Toeplitz covariance from cross-ACV G(:,:,k+1)=Gamma(k),
% epoch-major rows. Gamma(-k) = Gamma(k)'.
    m = d * nE;  C = zeros(m);
    for i = 1:nE
        for j = 1:nE
            k = i - j;
            if k >= 0, blk = G(:,:,k+1);
            else,      blk = G(:,:,-k+1).'; end
            C((i-1)*d+(1:d), (j-1)*d+(1:d)) = blk;
        end
    end
    C = 0.5 * (C + C.');
end

function acc = buildCleanAcc(cfg, days)
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false, ...
                                  'KeepLocal',true);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, @(d,ot,a,uo) lsqvce.engine.synthLoader(cfg,s,d,ot,a,uo), days);
    end
end

function f = spikedBinder(cfg, s, tgt)
% Loader wrapper that adds a spike to arc.l for one (sat,obsType,day,arc).
    f = @(d, ot, a, uo) spikedLoader(cfg, s, d, ot, a, uo, tgt);
end

function out = spikedLoader(cfg, s, day, obsType, arcID, uo, tgt)
    out = lsqvce.engine.synthLoader(cfg, s, day, obsType, arcID, uo);
    if ischar(obsType) || isstring(obsType)                    % data call
        if s == tgt.sat && strcmp(char(obsType), tgt.obsType) && ...
           day == tgt.day && arcID == tgt.arc && ...
           isfield(out,'l') && ~isempty(out.l)
            row = min(tgt.row, numel(out.l));
            sig = std(out.l);  if sig == 0, sig = 1; end
            out.l(row) = out.l(row) + tgt.size * sig;
        end
    end
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
    cfg.trueSigma.pos = 1.0;
    cfg.trueSigma.kbr = 2.5;
    cfg.cov.pos = struct('type','diag', 'sigma', cfg.trueSigma.pos);
    cfg.cov.kbr = struct('type','diag', 'sigma', cfg.trueSigma.kbr);
end

%% ============================================ WP-134: a VAR(1), a VAR(0)
function testCovModelTakesAVar1(tc)
% categories: contract, precision
% WP-134 (review 3 finding 10): MATLAB drops a trailing singleton, so a
% VAR(1)'s d x d x 1 coefficient array has ndims 2 -- covModel tested
% ndims ~= 3 and refused every VAR(1), whittle(G, 1)'s own output included.
% Reference: the factory's whitener against the dense decorrelation of its own
% Cfcn (as test_parallel does for a VAR(2)). A VAR(0) is white noise in varACV
% (the companion form of order 0 raised an unnamed eig error).
    A1 = 0.5 * eye(3);  Sig = [1 0.2 0; 0.2 1 0.1; 0 0.1 1];
    d = 3;  nE = 30;  m = nE * d;
    [Cfcn, Wfcn, info] = lsqvce.engine.covModel(struct('type','var', 'ar',A1, 'sigma',Sig), ...
                                                (1:m).', struct());
    tc.verifyEqual(info.arOrder, 1);
    rs = RandStream('mt19937ar', 'Seed', 134);
    Ad = randn(rs, m, 3);
    [Aw, ~, ~] = Wfcn(1:m, Ad, zeros(m, 0), zeros(m, 1), struct());
    W = lsqvce.kernel.decorrFactor(Cfcn(1:m, struct()));
    % tol: as testVarWhitenEqualsDenseDecorrelation -- quadratic forms agree
    % at eps*cond; 1e-9 house headroom
    tc.verifyEqual(Aw.' * Aw, (W * Ad).' * (W * Ad), 'RelTol', 1e-9);
    [Ah, Sh] = lsqvce.kernel.whittle(lsqvce.kernel.varACV(A1, Sig, 3), 1);
    [~, W2] = lsqvce.engine.covModel(struct('type','var', 'ar',Ah, 'sigma',Sh), (1:m).', struct());
    tc.verifyNotEmpty(W2, 'whittle(G, 1)''s output is accepted');
    G0 = lsqvce.kernel.varACV(zeros(3, 3, 0), Sig, 4);
    tc.verifyEqual(G0(:, :, 1), Sig);
    tc.verifyEqual(G0(:, :, 2:end), zeros(3, 3, 4));
end

function testAr2armaNormalisesTheLeadingCoefficient(tc)
% categories: metamorphic, contract
% WP-134 (review 3 low): b and s2 carried 1/aHi(1) and its square when
% aHi(1) ~= 1. c aHi with innovation variance c^2 s2Hi is the same process,
% so the reduction must return the same (a, b, s2) -- bitwise for c = 2,
% whose normalisation is exact. aHi(1) = 0 is refused by name.
    rad = 0.94;  ang = 2*pi*0.12;
    aTrue = real(poly([rad*exp(1i*ang), rad*exp(-1i*ang)]));
    r = lsqvce.kernel.armaACV(aTrue, [1 0.5], 256);
    [aHi, s2Hi] = lsqvce.kernel.levinsonDurbin(r, 60);
    opt = struct('armaP', 2, 'armaN', 4096);
    [a1, b1, s1] = lsqvce.kernel.ar2arma(aHi, s2Hi, 1, opt);
    [a2, b2, s2] = lsqvce.kernel.ar2arma(2 * aHi, 4 * s2Hi, 1, opt);
    tc.verifyEqual([a2, b2, s2], [a1, b1, s1]);
    tc.verifyError(@() lsqvce.kernel.ar2arma([0, aHi(2:end)], s2Hi, 1, opt), 'lsqvce:ar2arma:degenerate');
end
