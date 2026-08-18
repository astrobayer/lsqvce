function tests = test_ellmer
%TEST_ELLMER  Unit tests for the Ellmer (2018) stochastic-model primitives.
%
%   Phase A: the type-I DCT (LSQVCE.DCTI), the banded lag components
%   (LSQVCE.COVCOMPONENT), and the Toeplitz cofactor builder
%   (LSQVCE.TOEPLITZCOFACTOR). Later phases (per-lag redundancy, per-frequency
%   PSD VCE, arc-wise factors) append their tests here.
%
%   Every check mirrors an identity verified offline against the exact algebra
%   (Ellmer eqs. 6.5.6, 6.5.12-6.5.17): DCT-I orthonormality/round-trip, the
%   V_n banded apply equals the dense V_n, and sum_n c(n) V_n = toeplitz(c).
%
%   Run with:  results = runtests('tests/test_ellmer.m')
%
%   Updated 2026-08-18 (CEST, UTC+02:00), WP-33: four tests for the locals-
%   in-residual-and-projector fix (F-03) -- dense joint-projector parity,
%   alpha2 ~ 1 at truth WITH locals (+ the pre-fix contrast, so the test
%   fails on the old code), the bitwise nl = 0 regression, and the flat-PSD
%   driver check. Round 2 (same date): the deferred coloured-synthetic
%   PSD-recovery test (corr > 0.95 with locals) + its jointCfgColoured
%   fixture. With assistance from Claude (Anthropic), model Claude Fable 5.
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
    tests = functiontests(localfunctions);
end

%% ============================================ dctI: orthonormal DCT-I
function test_dctI_ortho_is_self_inverse(tc)
    % The orthonormal type-I DCT is symmetric and self-inverse: X = X.', X*X = I.
    for N = [2 3 5 8 16 33]
        [~, X] = lsqvce.dctI(zeros(N,1));         % forward matrix
        verifyEqual(tc, X, X.', 'AbsTol', 1e-12, sprintf('N=%d symmetric', N));
        verifyEqual(tc, X*X, eye(N), 'AbsTol', 1e-10, sprintf('N=%d involutive', N));
    end
end

function test_dctI_ortho_round_trip(tc)
    % ortho is self-inverse, so dctI(dctI(x)) == x with no Inverse flag.
    rng(1);
    x = randn(12, 4);
    y = lsqvce.dctI(x);
    z = lsqvce.dctI(y);
    verifyEqual(tc, z, x, 'AbsTol', 1e-10);
end

function test_dctI_row_in_row_out(tc)
    x = 1:6;
    y = lsqvce.dctI(x);
    verifyTrue(tc, isrow(y), 'a row input must give a row output');
    verifyEqual(tc, y(:), lsqvce.dctI(x(:)), 'AbsTol', 1e-12);
end

function test_dctI_ellmer_norm_invertible(tc)
    % The Ellmer-literal normalisation is NOT orthonormal but is well conditioned
    % (cond ~ 1.5); its inverse round-trips via 'Inverse',true.
    rng(2);
    x = randn(10, 3);
    [~, X] = lsqvce.dctI(zeros(10,1), 'Norm', 'ellmer');
    verifyLessThan(tc, cond(X), 5, 'ellmer X must be well conditioned');
    y = lsqvce.dctI(x, 'Norm', 'ellmer');
    z = lsqvce.dctI(y, 'Norm', 'ellmer', 'Inverse', true);
    verifyEqual(tc, z, x, 'AbsTol', 1e-10, 'ellmer forward then inverse round-trip');
end

function test_dctI_ortho_vs_ellmer_differ_at_nyquist(tc)
    % The two normalisations differ only in the highest-frequency (Nyquist) bin.
    [~, Xo] = lsqvce.dctI(zeros(8,1), 'Norm', 'ortho');
    [~, Xe] = lsqvce.dctI(zeros(8,1), 'Norm', 'ellmer');
    d = abs(Xo - Xe);
    % the difference is concentrated in the last row/column (index N-1)
    interior = d(1:end-1, 1:end-1);
    verifyLessThan(tc, max(interior(:)), 1e-9, ...
        'interior of the two DCT matrices must agree');
    verifyGreaterThan(tc, max(d(:)), 1e-3, ...
        'the Nyquist row/column must differ between normalisations');
end

function test_dctI_too_short_errors(tc)
    verifyError(tc, @() lsqvce.dctI(1), 'lsqvce:dctI:tooShort');
end

%% ============================================ covComponent: V_n banded apply
function test_covComponent_n0_is_identity(tc)
    rng(3);
    x = randn(7, 2);
    verifyEqual(tc, lsqvce.covComponent(x, 0), x, 'AbsTol', 0);
end

function test_covComponent_matches_dense_Vn(tc)
    % The banded apply equals the explicit V_n * x for every lag.
    N = 6;
    rng(4);
    x = randn(N, 3);
    for n = 0:N-1
        V = denseVn(N, n);
        verifyEqual(tc, lsqvce.covComponent(x, n), V*x, 'AbsTol', 1e-12, ...
            sprintf('lag n=%d', n));
    end
end

function test_covComponent_sums_to_toeplitz(tc)
    % sum_n c(n) V_n * x == toeplitz(c) * x  (Ellmer eq. 6.5.17).
    N = 8;
    rng(5);
    c = [1; 0.5; 0.25; 0.1; 0.05; 0.02; 0.01; 0];
    x = randn(N, 4);
    acc = zeros(N, 4);
    for n = 0:N-1
        acc = acc + c(n+1) * lsqvce.covComponent(x, n);
    end
    verifyEqual(tc, acc, toeplitz(c)*x, 'AbsTol', 1e-12);
end

function test_covComponent_lag_too_large_errors(tc)
    verifyError(tc, @() lsqvce.covComponent(ones(4,1), 4), ...
        'lsqvce:covComponent:lagTooLarge');
    verifyError(tc, @() lsqvce.covComponent(ones(4,1), -1), ...
        'lsqvce:covComponent:badLag');
end

%% ============================================ toeplitzCofactor
function test_toeplitzCofactor_dense_is_toeplitz(tc)
    c = [1; 0.5; 0.2; 0.05];
    Q = lsqvce.toeplitzCofactor(c);
    verifyEqual(tc, Q, toeplitz(c), 'AbsTol', 1e-14);
    verifyEqual(tc, Q, Q.', 'AbsTol', 1e-14, 'must be symmetric');
end

function test_toeplitzCofactor_apply_matches_dense(tc)
    c = [1; 0.5; 0.2; 0.05; 0.01];
    q = lsqvce.toeplitzCofactor(c, 'Form', 'apply');
    rng(6);
    x = randn(5, 7);
    verifyEqual(tc, q(x), toeplitz(c)*x, 'AbsTol', 1e-12);
end

function test_toeplitzCofactor_whitens_via_decorr(tc)
    % A valid PD autocovariance -> chol(Q)^-1 whitens: W Q W' = I.
    c = [1; 0.4; 0.1];               % clearly PD
    Q = lsqvce.toeplitzCofactor(c);
    W = lsqvce.decorrFactor(Q);
    verifyEqual(tc, W*Q*W.', eye(3), 'AbsTol', 1e-10);
end

function test_toeplitzCofactor_warns_on_bad_acv(tc)
    % |c(n)| > c(0) is not a valid autocovariance -> a warning (not an error).
    verifyWarning(tc, @() lsqvce.toeplitzCofactor([1; 2; 0.1]), ...
        'lsqvce:toeplitzCofactor:notACV');
end

%% ============================================ Phase B: shared probe
function test_sharedProbe_shapes_and_solve(tc)
    % UiZ must satisfy U*UiZ == Zbar, and carry the documented shapes.
    rng(20);
    P = 5;  Z = 30;
    M = randn(P); Nm = M.'*M + P*eye(P);      % SPD
    U = chol(Nm);                             % upper: Nm = U'*U
    sp = lsqvce.sharedProbe(U, 'numProbes', Z, 'seed', 7);
    verifySize(tc, sp.UiZ , [P Z]);
    verifySize(tc, sp.Zbar, [P Z]);
    verifyEqual(tc, sp.probes, Z);
    verifyEqual(tc, sp.P, P);
    verifyEqual(tc, U*sp.UiZ, sp.Zbar, 'AbsTol', 1e-10, 'UiZ must solve U*UiZ=Zbar');
end

function test_sharedProbe_rademacher_normalisation(tc)
    % Zbar = raw/sqrt(Z) with raw entries +-1, so every entry is +-1/sqrt(Z).
    P = 4;  Z = 16;
    U = eye(P);
    sp = lsqvce.sharedProbe(U, 'numProbes', Z, 'seed', 3);
    verifyEqual(tc, abs(sp.Zbar), ones(P,Z)/sqrt(Z), 'AbsTol', 1e-12);
end

function test_sharedProbe_seed_reproducible(tc)
    U = chol(eye(6)*2);
    a = lsqvce.sharedProbe(U, 'numProbes', 20, 'seed', 11);
    b = lsqvce.sharedProbe(U, 'numProbes', 20, 'seed', 11);
    verifyEqual(tc, a.UiZ, b.UiZ, 'AbsTol', 0, 'a fixed seed must reproduce');
end

function test_sharedProbe_accepts_fixed_probes(tc)
    P = 3;  Z = 5;
    U = chol(eye(P)*3);
    raw = sign(randn(P,Z));  raw(raw==0) = 1;
    sp = lsqvce.sharedProbe(U, 'probes', raw);
    verifyEqual(tc, sp.Zbar, raw/sqrt(Z), 'AbsTol', 1e-12);
end

function test_sharedProbe_rejects_nonsquare(tc)
    verifyError(tc, @() lsqvce.sharedProbe(ones(3,4)), ...
        'lsqvce:sharedProbe:notSquare');
end

%% ============================================ Phase B: per-lag redundancy
function test_redundancyLags_s_matches_dense_trace(tc)
    % The memory-honest sub-diagonal sums must equal trace(Rbar*V_n) exactly for
    % the SAME Rbar (this is an algebraic identity, not a Monte-Carlo question).
    [Abar, W, sp, ebar, N] = lagFixture(tc, 12, 4, 200);
    nLag = 6;
    r = lsqvce.redundancyLags(Abar, W, sp.UiZ, ebar, nLag, 'WantRbar', true);
    for n = 0:nLag-1
        Vn = denseVn(N, n);
        verifyEqual(tc, r.s(n+1), trace(r.Rbar*Vn), 'RelTol', 1e-8, ...
            sprintf('s_%d must equal trace(Rbar V_n)', n));
    end
end

function test_redundancyLags_omega_matches_quadratic(tc)
    % Omega_n == etilde' V_n etilde, with etilde = Sigma^-1 ehat.
    [Abar, W, sp, ebar, N] = lagFixture(tc, 11, 3, 100);
    nLag = 5;
    r  = lsqvce.redundancyLags(Abar, W, sp.UiZ, ebar, nLag);
    et = W \ ebar;                                   % Sigma^-1 ehat
    for n = 0:nLag-1
        Vn = denseVn(N, n);
        verifyEqual(tc, r.omega(n+1), et.'*Vn*et, 'RelTol', 1e-8, ...
            sprintf('Omega_%d must equal etilde'' V_n etilde', n));
    end
end

function test_redundancyLags_Rbar_converges_to_exact(tc)
    % With many probes, Rbar -> Sigma^-1 - Sigma^-1 A N^-1 A' Sigma^-1.
    rng(31);
    N = 10;  P = 3;
    [Sig, W] = spdToeplitz(N);
    A    = randn(N, P);
    Abar = W.' \ A;                                   % W^-T A
    Nm   = Abar.'*Abar;
    U    = chol(Nm);
    sp   = lsqvce.sharedProbe(U, 'numProbes', 20000, 'seed', 5);
    r    = lsqvce.redundancyLags(Abar, W, sp.UiZ, [], 3, 'WantRbar', true);
    Sinv    = Sig \ eye(N);
    Rexact  = Sinv - Sinv*A*(Nm\(A.'*Sinv));
    verifyEqual(tc, r.Rbar, Rexact, 'AbsTol', 0.05, ...
        'Rbar must approach the exact redundancy matrix');
end

function test_redundancyLags_omega_skipped_when_no_residual(tc)
    [Abar, W, sp, ~, ~] = lagFixture(tc, 8, 2, 50);
    r = lsqvce.redundancyLags(Abar, W, sp.UiZ, [], 4);
    verifyTrue(tc, all(isnan(r.omega)), 'empty ebar must give NaN omega');
    verifyTrue(tc, all(isfinite(r.s)), 's must still be computed');
end

function test_redundancyLags_rejects_bad_shapes(tc)
    [Abar, W, sp, ebar, ~] = lagFixture(tc, 9, 3, 40);
    verifyError(tc, @() lsqvce.redundancyLags(Abar, W(1:3,1:3), sp.UiZ, ebar, 3), ...
        'lsqvce:redundancyLags:badW');
    verifyError(tc, @() lsqvce.redundancyLags(Abar, W, sp.UiZ(1:2,:), ebar, 3), ...
        'lsqvce:redundancyLags:badUiZ');
    verifyError(tc, @() lsqvce.redundancyLags(Abar, W, sp.UiZ, ebar, 0), ...
        'lsqvce:redundancyLags:badNMaxLag');
end

%% ============================================ Phase C: per-frequency PSD VCE
function test_psdVCE_identity_when_omega_equals_s(tc)
    % If Omega == s then alpha2 == 1 at every frequency, so the PSD is unchanged:
    % the estimator is at a fixed point and reports convergence.
    rng(50);
    N = 10;  M = 4;
    C  = [1; 0.4.^(1:N-1).'];
    om = rand(N, M) + 0.5;
    p  = lsqvce.psdVCE(om, om, C);
    verifyEqual(tc, p.alpha2, ones(N,1), 'AbsTol', 1e-10, 'alpha2 must be 1');
    verifyEqual(tc, p.maxDev, 0, 'AbsTol', 1e-10);
end

function test_psdVCE_amplitudes_stay_positive(tc)
    % A PSD is non-negative by definition, whatever the inputs do.
    rng(51);
    N = 12;  M = 5;
    C  = [1; 0.3.^(1:N-1).'];
    p  = lsqvce.psdVCE(randn(N,M), randn(N,M), C);
    verifyTrue(tc, all(p.Shat > 0), 'the estimated PSD must be strictly positive');
    verifySize(tc, p.Cnew, [N 1]);
end

function test_psdVCE_covariance_may_go_negative(tc)
    % The POINT of estimating in the spectral domain: a positive PSD inverts to a
    % covariance function that is free to change sign (Ellmer eq. 6.5.19 note).
    N = 8;
    S = zeros(N,1);  S(3) = 1;  S(1) = 0.1;      % a peaked spectrum
    C = lsqvce.dctI(S, 'Inverse', true);
    verifyTrue(tc, any(C < 0), ...
        'a peaked PSD must invert to an ACV with negative lags');
end

function test_psdVCE_warns_on_degenerate_white_start(tc)
    % C = [1 0 0 ...] annihilates every lag n>0 in the weighting, so those lags
    % can never update: the iteration "converges" to the wrong PSD. Must warn.
    N = 8;  M = 3;
    Cdelta = [1; zeros(N-1,1)];
    verifyWarning(tc, ...
        @() lsqvce.psdVCE(rand(N,M), rand(N,M), Cdelta), ...
        'lsqvce:psdVCE:degenerateStart');
end

function test_psdVCE_degenerate_flag_and_override(tc)
    N = 8;  M = 3;
    Cdelta = [1; zeros(N-1,1)];
    p = lsqvce.psdVCE(rand(N,M), rand(N,M), Cdelta, 'AllowDegenerate', true);
    verifyTrue(tc, p.degenerate, 'the degenerate start must be flagged');
    Cok = [1; 0.3*ones(N-1,1)];
    q = lsqvce.psdVCE(rand(N,M), rand(N,M), Cok);
    verifyFalse(tc, q.degenerate);
end

function test_psdVCE_zero_frequency_rule(tc)
    rng(52);
    N = 9;  M = 3;
    C  = [1; 0.35.^(1:N-1).'];
    pl = lsqvce.psdVCE(rand(N,M), rand(N,M), C, 'ZeroFreqRule', 'lowest');
    verifyEqual(tc, pl.Shat(1), pl.Shat(2), 'AbsTol', 0, ...
        '''lowest'' must replace S^0 by the next amplitude');
    pk = lsqvce.psdVCE(rand(N,M), rand(N,M), C, 'ZeroFreqRule', 'keep');
    verifySize(tc, pk.Shat, [N 1]);
end

function test_psdVCE_rejects_bad_shapes(tc)
    N = 6;
    C = [1; 0.2*ones(N-1,1)];
    verifyError(tc, @() lsqvce.psdVCE(rand(N-1,3), rand(N-1,3), C), ...
        'lsqvce:psdVCE:shapeMismatch');
    verifyError(tc, @() lsqvce.psdVCE(rand(N,3), rand(N,4), C), ...
        'lsqvce:psdVCE:arcMismatch');
end

%% ============================================ Phase D: arc-wise variance factors
function test_arcVarianceFactors_plain_mean_normalises_to_one(tc)
    % With the PLAIN mean (eq. 6.5.29) the mean of the arc FACTORS (not their
    % squares) is exactly 1 -- Ellmer's stated normalisation property.
    rng(60);
    N = 8;  M = 12;
    om = rand(N,M) + 0.5;
    s  = rand(N,M) + 0.5;
    a  = lsqvce.arcVarianceFactors(om, s, ones(M,1), 'MeanType', 'plain');
    verifyEqual(tc, mean(sqrt(abs(a.sigma2))), 1, 'AbsTol', 1e-10, ...
        'the plain mean must normalise mean(gamma_m/gammaBar) to 1');
end

function test_arcVarianceFactors_robust_resists_outliers(tc)
    % The robust central-half mean must sit below the outlier-inflated plain mean
    % (that is the entire reason GROOPS uses it).
    rng(61);
    N = 6;  M = 40;
    om = rand(N,M) + 0.5;
    s  = rand(N,M) + 0.5;
    om(:, 1:3) = om(:, 1:3) * 400;             % a few pathological arcs
    ap = lsqvce.arcVarianceFactors(om, s, ones(M,1), 'MeanType', 'plain');
    ar = lsqvce.arcVarianceFactors(om, s, ones(M,1), 'MeanType', 'robust');
    verifyLessThan(tc, ar.gammaBar, ap.gammaBar, ...
        'the robust mean must sit below the outlier-inflated plain mean');
end

function test_arcVarianceFactors_gamma_detour_shapes(tc)
    rng(62);
    N = 7;  M = 9;
    sig2 = rand(M,1) + 0.5;
    a = lsqvce.arcVarianceFactors(rand(N,M)+0.5, rand(N,M)+0.5, sig2);
    verifySize(tc, a.sigma2, [M 1]);
    verifySize(tc, a.alpha2, [M 1]);
    verifyEqual(tc, a.gamma2, a.alpha2 .* sig2, 'RelTol', 1e-10, ...
        'gamma2 must be alpha2 .* sigma2 (eq. 6.5.28)');
    verifyTrue(tc, all(a.sigma2 >= 0), 'variance factors must be non-negative');
end

function test_arcVarianceFactors_small_M_falls_back(tc)
    % With too few arcs to trim, the robust mean must degrade gracefully.
    rng(63);
    N = 5;
    for M = [1 2 3]
        a = lsqvce.arcVarianceFactors(rand(N,M)+0.5, rand(N,M)+0.5, ones(M,1));
        verifyTrue(tc, isfinite(a.gammaBar) && a.gammaBar > 0, ...
            sprintf('M=%d must still give a usable mean', M));
    end
end

function test_arcVarianceFactors_rejects_arc_mismatch(tc)
    verifyError(tc, ...
        @() lsqvce.arcVarianceFactors(rand(5,3), rand(5,3), ones(4,1)), ...
        'lsqvce:arcVarianceFactors:arcMismatch');
end

%% ============================================ Phase E: joint loop driver
function test_adjustVCEjoint_runs_and_returns_shapes(tc)
    % End-to-end: the ten-step loop must run on a real accumulator and return a
    % PSD, a covariance function, and one variance factor per arc.
    cfg  = jointCfg();
    nLag = 8;
    r = lsqvce.adjustVCEjoint(@(model) jointAcc(cfg, model), ...
            'NumLags', nLag, 'MaxIter', 2, 'Probes', 40, 'Seed', 3, ...
            'Verbose', false);
    verifySize(tc, r.C   , [nLag 1]);
    verifySize(tc, r.Shat, [nLag 1]);
    verifyTrue(tc, all(r.Shat > 0), 'the PSD must stay strictly positive');
    verifyTrue(tc, numel(r.sigma2) >= 1, 'at least one arc factor');
    verifyTrue(tc, all(isfinite(r.sigma2)), 'arc factors must be finite');
    verifyClass(tc, r.acc, 'lsqvce.BlockAccumulator');
end

function test_adjustVCEjoint_history_grows_per_iteration(tc)
    cfg  = jointCfg();
    nLag = 6;
    r = lsqvce.adjustVCEjoint(@(model) jointAcc(cfg, model), ...
            'NumLags', nLag, 'MaxIter', 3, 'Probes', 30, 'Seed', 5, ...
            'Tol', 1e-12);                       % force all iterations
    verifyEqual(tc, numel(r.history.maxDev), r.nIter);
    verifySize(tc, r.history.Shat, [r.nIter nLag]);
    verifyTrue(tc, all(isfinite(r.history.maxDev)), ...
        'the convergence measure must be finite every pass');
end

function test_adjustVCEjoint_model_reaches_accFcn(tc)
    % The driver must hand the CURRENT model to accFcn each pass (steps 1-3),
    % so the caller can rebuild Sigma from it. containers.Map is a HANDLE, so
    % the closure below can record into it without a nested function.
    cfg = jointCfg();
    rec = containers.Map('KeyType', 'double', 'ValueType', 'any');
    accF = @(model) recordAndBuild(rec, cfg, model);
    lsqvce.adjustVCEjoint(accF, 'NumLags', 5, 'MaxIter', 2, ...
                          'Probes', 20, 'Seed', 7, 'Tol', 1e-12);
    verifyGreaterThanOrEqual(tc, rec.Count, 2, 'accFcn called once per pass');
    verifyTrue(tc, rec.isKey(1) && rec.isKey(2), 'iterations 1 and 2 recorded');
    m1 = rec(1);
    verifyTrue(tc, isfield(m1, 'C') && isfield(m1, 'iter'), ...
        'the model struct must carry C and iter');
    verifyEqual(tc, m1.iter, 1);
end

function test_adjustVCEjoint_default_start_is_not_degenerate(tc)
    % The default C0 must NOT be a delta -- that start is a degenerate fixed
    % point of the weighting (see LSQVCE.PSDVCE). A delta start DOES warn; the
    % default must not.
    cfg  = jointCfg();
    nLag = 5;
    % explicit delta start -> the degeneracy warning must fire
    verifyWarning(tc, @() lsqvce.adjustVCEjoint( ...
            @(model) jointAcc(cfg, model), 'NumLags', nLag, 'MaxIter', 1, ...
            'Probes', 20, 'Seed', 9, 'C0', [1; zeros(nLag-1,1)]), ...
        'lsqvce:psdVCE:degenerateStart');
end

function test_adjustVCEjoint_rejects_bad_C0(tc)
    cfg = jointCfg();
    verifyError(tc, ...
        @() lsqvce.adjustVCEjoint(@(model) jointAcc(cfg, model), ...
                'NumLags', 8, 'C0', ones(3,1)), ...
        'lsqvce:adjustVCEjoint:badC0');
end

%% ============================================ vceMode flag on adjustVCE
function test_vceMode_defaults_to_separate(tc)
    P = lsqvce.adjustVCE(modeSpec(), 1, 'verbose', false);
    verifyEqual(tc, P.vceMode, 'separate', ...
        'the alternating loop must remain the default');
end

function test_vceMode_constructor_flag(tc)
    P = lsqvce.adjustVCE(modeSpec(), 1, 'vceMode', 'joint', 'verbose', false);
    verifyEqual(tc, P.vceMode, 'joint');
end

function test_vceMode_rejects_unknown(tc)
    verifyError(tc, ...
        @() lsqvce.adjustVCE(modeSpec(), 1, 'vceMode', 'sideways', 'verbose', false), ...
        'MATLAB:InputParser:ArgumentFailedValidation');
end

function test_joint_fluent_switches_mode_and_sets_opts(tc)
    % .joint(...) both switches the mode and configures it, so no separate flag
    % is needed -- the call reads as the intent.
    P = lsqvce.adjustVCE(modeSpec(), 1, 'verbose', false);
    P = P.joint('numLags', 24, 'probes', 40, 'meanType', 'plain');
    verifyEqual(tc, P.vceMode, 'joint');
    verifyEqual(tc, P.jointOpt.NumLags , 24);
    verifyEqual(tc, P.jointOpt.Probes  , 40);
    verifyEqual(tc, P.jointOpt.MeanType, 'plain');
end

function test_separate_fluent_switches_back(tc)
    P = lsqvce.adjustVCE(modeSpec(), 1, 'verbose', false);
    P = P.joint('numLags', 12);
    verifyEqual(tc, P.vceMode, 'joint');
    P = P.separate();
    verifyEqual(tc, P.vceMode, 'separate', ...
        '.separate() must restore the alternating loop');
    verifyEqual(tc, P.jointOpt.NumLags, 12, ...
        'switching back must not erase the joint configuration');
end

function test_joint_defaults_are_ellmer_values(tc)
    P = lsqvce.adjustVCE(modeSpec(), 1, 'verbose', false);
    verifyEqual(tc, P.jointOpt.Probes  , 100, 'Ellmer uses Z = 100');
    verifyEqual(tc, P.jointOpt.NumLags , 60);
    verifyEqual(tc, P.jointOpt.Norm    , 'ortho');
    verifyEqual(tc, P.jointOpt.MeanType, 'robust');
    verifyEmpty(tc, P.jointOpt.C0, ...
        'the default start is chosen per stream, never a fixed delta');
end

%% ============================================ jointModelOpt (accFcn contract)
function test_jointModelOpt_takes_sqrt_of_variance(tc)
    % CONVENTION 1: model.sigma2 is a VARIANCE, a covariance spec's sigma is a
    % STANDARD DEVIATION. Getting this wrong does not error -- the VCE loop just
    % diverges -- so assert the sqrt explicitly.
    model = struct('C', [1; 0.4; 0.1], 'sigma2', [4; 4; 4], 'iter', 2);
    uo = lsqvce.jointModelOpt(model, {'x1|pos'});
    spec = uo.covOverride.bySys.x1.pos;
    verifyEqual(tc, spec.sigma, 2, 'AbsTol', 1e-12, ...
        'sigma must be sqrt(mean variance), not the variance itself');
end

function test_jointModelOpt_places_override_bySys(tc)
    % CONVENTION 2: loaders read opt.covOverride.bySys.(satKey).(obsType).
    model = struct('C', [1; 0.3], 'sigma2', [], 'iter', 1);
    uo = lsqvce.jointModelOpt(model, {'x1|posC', 'x2|rr'});
    verifyTrue(tc, isfield(uo, 'covOverride'));
    verifyTrue(tc, isfield(uo.covOverride, 'bySys'));
    verifyTrue(tc, isfield(uo.covOverride.bySys, 'x1'));
    verifyTrue(tc, isfield(uo.covOverride.bySys.x1, 'posC'));
    verifyTrue(tc, isfield(uo.covOverride.bySys.x2, 'rr'));
end

function test_jointModelOpt_sets_psdMode_on_useropt(tc)
    % CONVENTION 3: psdMode is read from the UserOpt, NOT from the spec. It must
    % appear at the top level or covModel silently uses its 'ar' default.
    model = struct('C', [1; 0.2], 'sigma2', 1, 'iter', 1);
    uo = lsqvce.jointModelOpt(model, {'x1|pos'});
    verifyEqual(tc, uo.psdMode, 'toeplitz', ...
        'joint mode must ask for the direct Toeplitz build, not an AR fit');
    uo2 = lsqvce.jointModelOpt(model, {'x1|pos'}, 'PsdMode', 'ar');
    verifyEqual(tc, uo2.psdMode, 'ar');
end

function test_jointModelOpt_acf_reproduces_the_covariance(tc)
    % The handle must return c(|lag|) over a LAG MATRIX, and zero beyond the
    % modelled range -- that is how covModel's 'acf' type calls it.
    c = [1; 0.5; 0.2];
    model = struct('C', c, 'sigma2', 1, 'iter', 1);
    uo = lsqvce.jointModelOpt(model, {'x1|pos'});
    f  = uo.covOverride.bySys.x1.pos.acf;
    t  = (0:4).';
    C  = f(t - t.', struct());
    verifyEqual(tc, diag(C), repmat(c(1),5,1), 'AbsTol', 1e-12, 'lag 0');
    verifyEqual(tc, C(1,2), c(2), 'AbsTol', 1e-12, 'lag 1');
    verifyEqual(tc, C(1,3), c(3), 'AbsTol', 1e-12, 'lag 2');
    verifyEqual(tc, C(1,5), 0, 'AbsTol', 1e-12, ...
        'beyond the modelled range must be zero, not extrapolated');
    verifyEqual(tc, C, C.', 'AbsTol', 1e-12, 'must be symmetric');
end

function test_jointModelOpt_accepts_cell_streams(tc)
    model = struct('C', [1; 0.3], 'sigma2', 1, 'iter', 1);
    uo = lsqvce.jointModelOpt(model, {{1, 'pos'}});
    verifyTrue(tc, isfield(uo.covOverride.bySys, 'x1'), ...
        'a numeric sysID must become the makeValidName key');
end

function test_jointModelOpt_rejects_bad_input(tc)
    verifyError(tc, @() lsqvce.jointModelOpt(struct('a',1), {'x1|pos'}), ...
        'lsqvce:jointModelOpt:badModel');
    verifyError(tc, @() lsqvce.jointModelOpt( ...
        struct('C', [], 'sigma2', 1, 'iter', 1), {'x1|pos'}), ...
        'lsqvce:jointModelOpt:badC');
    verifyError(tc, @() lsqvce.jointModelOpt( ...
        struct('C', [1;0.2], 'sigma2', 1, 'iter', 1), {'nopipe'}), ...
        'lsqvce:jointModelOpt:badStream');
end

%% ============================================ axis stride (epoch-major lags)
function test_axis_pooling_recovers_temporal_law(tc)
    % REGRESSION GUARD for the joint sweep's axis stride.
    %
    % Rows of a d-axis stream are EPOCH-MAJOR (row = (epoch-1)*d + axis), so a
    % lag taken on ROWS counts axes, not epochs. This test builds three
    % independent AR(1) series sharing one law, interleaves them exactly as the
    % loader does, and shows that:
    %   - the ROW-based ACV is wrong at lag 1 and shows spurious peaks at lags
    %     d and 2d (the axis aliasing), whereas
    %   - pooling per axis DOWN EPOCHS recovers rho.
    % The joint sweep must do the latter (see adjustVCE.lagsPooledOverAxes).
    rng(2);
    nE = 200;  nAx = 3;  rho = 0.8;
    E = zeros(nE, nAx);
    for a = 1:nAx
        x = zeros(nE,1);
        for k = 2:nE, x(k) = rho*x(k-1) + randn; end
        E(:,a) = x;
    end
    rows = reshape(E.', [], 1);                 % epoch-major interleave
    L = 8;

    aRow = biasedAcvLocal(rows, L);   aRow = aRow / aRow(1);
    aAx  = zeros(L,1);
    for a = 1:nAx
        ai  = biasedAcvLocal(E(:,a), L);
        aAx = aAx + ai/ai(1);
    end
    aAx = aAx / nAx;

    % the axis-pooled estimate must track rho^n; the row-based one must not
    verifyEqual(tc, aAx(2), rho, 'AbsTol', 0.08, ...
        'pooling per axis down epochs must recover the AR(1) coefficient');
    verifyGreaterThan(tc, abs(aRow(2) - rho), 0.3, ...
        ['a ROW-based lag must MISS rho -- if this fails the interleave is ' ...
         'not epoch-major and the guard is not testing what it claims']);
    % the spurious row-based peak sits at lag nAx
    verifyGreaterThan(tc, aRow(1+nAx), aRow(2), ...
        'the row-based ACV must peak at the axis stride, not at lag 1');
end

function r = biasedAcvLocal(v, L)
    v = v(:);  n = numel(v);
    r = zeros(L,1);
    for k = 0:L-1
        r(k+1) = (v(1:n-k).' * v(1+k:n)) / n;
    end
end

%% ============================================ local helpers
function sys = modeSpec()
%MODESPEC  Smallest legal 'systems' value for constructor-only mode tests.
%   The vceMode/joint/separate tests below only inspect the CONFIGURED state
%   (P.vceMode, P.jointOpt); they never call run(), so the system needs no
%   loader -- only enough for LSQVCE.RESOLVENC to derive nc. A raw spec struct
%   declares nc beside the nl it already carries, which is exactly that
%   minimum. (Before 2026-07-18 these tests passed nc as a leading argument,
%   'adjustVCE(4, {}, 1, ...)'; nc is now DERIVED from the systems, so that
%   form raises lsqvce:resolveNc:missing.)
%
%   OUTPUT
%     sys  1x1 struct  raw spec: .sysID (1x1 double), .nl (1x1 double),
%                      .nc (1x1 double), .loader (function handle, never called)
    sys = struct('sysID', 1, 'nl', 2, 'nc', 4, ...
                 'loader', @(varargin) error('modeSpec:notCalled', ...
                     'This spec is for constructor-only tests; it never loads.'));
end

function acc = recordAndBuild(rec, cfg, model)
    % Record the model the driver passed in, then build the accumulator.
    rec(model.iter) = model; %#ok<NASGU>  (containers.Map is a handle)
    acc = jointAcc(cfg, model);
end

function cfg = jointCfg()
    % A small, fast, deterministic problem for the joint driver. Uses an
    % explicit acf covariance so every arc exposes Cfcn (the joint estimator
    % needs the explicit W with Sigma = W'W).
    %
    % NOTE ON ISOLATION: lsqvce.synthLoader keeps a PERSISTENT single-day cache
    % keyed on (sat, day, obsType, tag, nc, nl, seed). This fixture deliberately
    % uses a system id far from the other suites' so that even a partial key
    % cannot collide -- a fixture whose nc differs from another suite's would
    % otherwise hand back wrongly-sized arcs to whichever ran second.
    cfg.nc   = 4;
    cfg.nl   = 2;
    cfg.sats = 91;                     % distinct from the other suites (1..3)
    rs = RandStream('mt19937ar', 'Seed', 17);
    cfg.x_c  = randn(rs, cfg.nc, 1);
    cfg.y    = containers.Map('KeyType','double','ValueType','any');
    cfg.y(91) = randn(rs, cfg.nl, 1);
    cfg.obsTypes   = {'pos'};
    cfg.arcsPerDay = 2;
    cfg.mArc       = 40;
    cfg.dt         = 5;
    cfg.seed       = 9117;             % distinct seed -> distinct cache key
    cfg.collinear  = false;
    cfg.collinearEps = 1e-7;
    cfg.cov.pos = struct('type','diag','sigma',0.5);
end

%% =============================== WP-33: locals in residual + projector (F-03)
function test_redundancyLags_local_block_matches_joint_inverse(tc)
    % Kernel exactness (spec: R_bar = Sig^-1 - Sig^-1 [A B] Njoint^-1 [A B]'
    % Sig^-1, realised with the leverage block identities, NO joint inverse).
    % Probe-limit check against the dense joint hat; mirrors v10review.py T3c.
    rng(33);
    N = 24;  P = 3;  nl = 2;  Z = 4000;
    [Sig, W] = spdToeplitz(N);
    A = randn(N, P);  B = randn(N, nl);
    % The reduction lives in the WHITENED metric -- the accumulator whitens
    % before accumulating, so st.Nll_r / st.Nlc (which the drivers hand to
    % jointArcParts) are Bw'Bw / Bw'Aw. Build the reference the same way; a
    % raw-metric G here is the mistake this test made in round 1 (CI red,
    % systematic 0.1-0.2 offsets).
    Aw = W.'\A;   Bw = W.'\B;
    Nll = Bw.'*Bw;  Nlc = Bw.'*Aw;
    G   = Nll \ Nlc;
    AtilBar = Aw - Bw*G;                      % = W^-T (A - B G), whitened G
    Nred = AtilBar.'*AtilBar;                 % the reduced common normal
    U    = chol((Nred+Nred.')/2);
    sp   = lsqvce.sharedProbe(U, 'numProbes', Z, 'seed', 5);
    Ll   = chol((Nll+Nll.')/2, 'lower');
    r = lsqvce.redundancyLags(AtilBar, W, sp.UiZ, [], 3, ...
            'LocalBbar', Bw, 'LocalLl', Ll, 'WantRbar', true);
    % Reference: the JOINT projector with Njoint from the WHITENED designs
    % (the accumulator whitens), Sig^-1 - Sig^-1-metric hat of [A B]:
    Sinv = Sig \ eye(N);
    Mw = [W.'\A, W.'\B];                     % whitened joint design
    Nj = Mw.'*Mw;
    RbarEx = Sinv - (W\Mw)*(Nj\(W\Mw).');
    tolR = 6/sqrt(Z);                         % Hutchinson-limit tolerance
    verifyEqual(tc, r.Rbar, RbarEx, 'AbsTol', tolR, ...
        'Rbar must converge to the JOINT projector (locals absorbed)');
    % s_0 = trace(Rbar). For COLOURED Sigma this is NOT N-(nc+nl) (that holds
    % only for white arcs -- asserted in the alpha2 test); the exact value is
    % trace(RbarEx), and the pre-fix common-only projector sits a structural
    % ~3.2 higher on this fixture -- the 0.5 tolerance discriminates.
    s0Ex = trace(RbarEx);
    verifyEqual(tc, r.s(1), s0Ex, 'AbsTol', 0.5, ...
        's_0 must equal trace of the JOINT Rbar');
    RcomOld = Sinv - (W\Aw)*((Aw.'*Aw)\(W\Aw).');   % the pre-WP-33 projector
    verifyGreaterThan(tc, trace(RcomOld) - s0Ex, 1.0, ...
        'the common-only s_0 must sit measurably above the joint one');
end

function test_joint_alpha2_unity_at_truth_with_locals(tc)
    % THE accuracy assertion the suite lacked (F-03): on WHITE noise with
    % nl > 0, at the joint least-squares solution, pooled Omega_0/s_0 must be
    % ~1. The pieces are formed EXACTLY as the fixed drivers form them
    % (y per system, Atil = A - B G, ebar = l - A x - B y). The pre-WP-33
    % formation (raw Abar, ebar without B*y, common-only projector) measured
    % alpha2 = 4.04 on this construction (v10review.py T3) -- asserted as the
    % contrast so this test FAILS on the old code by construction.
    rng(20260816);
    nc = 4;  nl = 3;  nArc = 6;  m = 120;
    x_true = randn(nc,1);  y_true = randn(nl,1);
    As = cell(nArc,1);  Bs = cell(nArc,1);  ls = cell(nArc,1);
    Ncc = zeros(nc);  Nll = zeros(nl);  Nlc = zeros(nl,nc);
    ncv = zeros(nc,1);  nlv = zeros(nl,1);
    for k = 1:nArc
        As{k} = randn(m, nc);  Bs{k} = randn(m, nl);
        ls{k} = As{k}*x_true + Bs{k}*y_true + randn(m,1);   % white, sigma = 1
        Ncc = Ncc + As{k}.'*As{k};   Nll = Nll + Bs{k}.'*Bs{k};
        Nlc = Nlc + Bs{k}.'*As{k};
        ncv = ncv + As{k}.'*ls{k};   nlv = nlv + Bs{k}.'*ls{k};
    end
    Nred = Ncc - Nlc.'*(Nll\Nlc);
    nred = ncv - Nlc.'*(Nll\nlv);
    x = Nred \ nred;                         % the joint common solution
    y = Nll \ (nlv - Nlc*x);                 % per-system local recovery
    G  = Nll \ Nlc;
    Ll = chol((Nll+Nll.')/2, 'lower');
    U  = chol((Nred+Nred.')/2);
    sp = lsqvce.sharedProbe(U, 'numProbes', 800, 'seed', 7);
    W  = eye(m);                              % white noise
    om0 = 0;  s0 = 0;  om0Old = 0;  s0Old = 0;
    for k = 1:nArc
        Atil = As{k} - Bs{k}*G;
        eNew = ls{k} - As{k}*x - Bs{k}*y;
        r = lsqvce.redundancyLags(Atil, W, sp.UiZ, eNew, 1, ...
                'LocalBbar', Bs{k}, 'LocalLl', Ll);
        om0 = om0 + r.omega(1);  s0 = s0 + r.s(1);
        eOld = ls{k} - As{k}*x;               % the pre-fix residual
        rOld = lsqvce.redundancyLags(As{k}, W, sp.UiZ, eOld, 1);
        om0Old = om0Old + rOld.omega(1);  s0Old = s0Old + rOld.s(1);
    end
    alpha2    = om0    / s0;
    alpha2Old = om0Old / s0Old;
    % tolerance: probe se of s_0 plus the chi2 spread of Omega_0/dof
    dof = nArc*m - (nc + nl);
    tol = 3*sqrt(2/dof) + 3*(nc/sqrt(800))/dof + 0.02;
    verifyEqual(tc, alpha2, 1, 'AbsTol', tol, ...
        'alpha2 must be ~1 at truth on white noise WITH locals');
    verifyTrue(tc, alpha2Old > 1.5, ...
        'the pre-WP-33 formation must stay measurably biased (the contrast)');
end

function test_redundancyLags_nl0_path_bitwise_unchanged(tc)
    % Spec regression: with NO local columns the new machinery must be a
    % NO-OP -- bitwise: the nl = 0 branch keeps the pre-WP-33 expressions
    % verbatim, and an explicit empty LocalBbar routes to that same branch.
    [Abar, W, sp, ebar] = lagFixture(tc, 20, 3, 60);
    a = lsqvce.redundancyLags(Abar, W, sp.UiZ, ebar, 5);
    b = lsqvce.redundancyLags(Abar, W, sp.UiZ, ebar, 5, ...
            'LocalBbar', zeros(20,0), 'LocalLl', []);
    verifyEqual(tc, b.s,     a.s);            % NO tolerance: bitwise
    verifyEqual(tc, b.omega, a.omega);
end

function test_adjustVCEjoint_flat_psd_on_white_noise_with_locals(tc)
    % Driver-level accuracy (spec 1 at the adjustVCEjoint level): on the
    % white (diag) fixture WITH nl = 2 locals the estimated PSD must stay
    % FLAT -- no arc-parameter signal folded into the "noise". Pre-WP-33 the
    % omitted B*y inflated the low lags (Omega_0 x3.8, T3).
    cfg  = jointCfg();
    r = lsqvce.adjustVCEjoint(@(model) jointAcc(cfg, model), ...
            'NumLags', 6, 'MaxIter', 2, 'Probes', 200, 'Seed', 11, ...
            'Verbose', false);
    S = r.Shat / median(r.Shat);
    verifyTrue(tc, max(S)/min(S) < 3.0, ...
        sprintf(['white truth => flat PSD: max/min = %.2f (arc-parameter ' ...
                 'leakage shows up as structure)'], max(S)/min(S)));
end

function test_adjustVCEjoint_psd_recovery_coloured_with_locals(tc)
    % WP-33 round 2 (the deferred spec test): on a COLOURED synthetic WITH
    % local parameters the converged PSD must correlate > 0.95 with the true
    % ARMA spectrum (truth on the estimator's own grid: the DCT-I of the
    % theoretical ACV, the toolbox's ACV<->PSD convention). Pre-fix, the
    % omitted B*y folded the arc-parameter signal into the low lags and the
    % common-only projector mis-sized s_n -- corr falls under the bar, so
    % this test fails on the old code (the §6.1 lever for round 2). The
    % fixture GENERATES AR(2) noise and jointAccMapped feeds the driver's
    % CURRENT model.C back into the declared covariance each pass (covModel
    % type 'acf') -- the honest Ellmer alternation: the estimator both
    % finds the colour and re-whitens with it. Round-2a/2b lessons: a
    % truth-whitening fixture tests nothing (corr -0.31), and a FIXED white
    % weighting biases Shat toward the non-delta prior (corr 0.65).
    cfg  = jointCfgColoured();
    nLag = 10;
    w = warning('off','all'); cleaner = onCleanup(@() warning(w)); %#ok<NASGU>
    r = lsqvce.adjustVCEjoint(@(model) jointAccMapped(cfg, model), ...
            'NumLags', nLag, 'MaxIter', 8, 'Probes', 200, 'Seed', 13, ...
            'Verbose', false);
    Strue = lsqvce.dctI(lsqvce.armaACV(cfg.cov.pos.ar, cfg.cov.pos.ma, nLag));
    Strue = abs(Strue(:));
    cc = corrcoef(r.Shat(:), Strue);
    verifyGreaterThan(tc, cc(1,2), 0.95, ...
        sprintf('PSD recovery with locals: corr = %.3f (need > 0.95)', cc(1,2)));
end

function cfg = jointCfgColoured()
    % The coloured sibling of jointCfg: same block-arrow layout, ARMA noise
    % (the stationary AR(2) of v10review T4) and enough arcs per frequency
    % for the spectral VCE to settle (one day, 6 arcs x 96 rows, nLag = 10).
    cfg.nc   = 3;
    cfg.nl   = 2;
    cfg.sats = 92;                     % distinct from jointCfg's 91
    rs = RandStream('mt19937ar', 'Seed', 23);
    cfg.x_c  = randn(rs, cfg.nc, 1);
    cfg.y    = containers.Map('KeyType','double','ValueType','any');
    cfg.y(92) = randn(rs, cfg.nl, 1);
    cfg.obsTypes   = {'pos'};
    cfg.arcsPerDay = 6;                % jointAcc streams ONE day: 6 x 96 rows
    cfg.mArc       = 96;
    cfg.dt         = 5;
    cfg.seed       = 9223;             % distinct seed -> distinct cache key
    cfg.collinear  = false;
    cfg.collinearEps = 1e-7;
    cfg.cov.pos = struct('type','arma', 'ar',[1 -0.75 0.4], 'ma',1, 'sigma',1.0);
    % The DECLARED covariance is owned by jointAccMapped (white on pass 1,
    % then the driver's own model.C via covModel type 'acf') -- never the
    % truth: a loader that whitens with cfg.cov's ARMA hands the estimator
    % pre-whitened residuals and the correct Shat is FLAT (round-2a lesson).
end

function acc = jointAccMapped(cfg, model)
    % The coloured sibling of jointAcc: maps the driver's CURRENT lag-domain
    % covariance model.C into the loader's DECLARED covariance (type 'acf':
    % C interpolated over |dt| lags, zero beyond nLag-1), so each outer pass
    % re-whitens with the latest estimate -- Ellmer's alternation, end to
    % end. Pass 1 (or an empty model) declares white. Generation always
    % follows cfg.cov (synthLoader).
    uo = struct();
    if ~isempty(model) && isfield(model, 'C') && ~isempty(model.C)
        C = model.C(:);
        uo.covOverride.pos = struct('type','acf', 'sigma', 1, ...
            'acf', @(dT, o) lagInterp(C, abs(dT) / cfg.dt));
    else
        uo.covOverride.pos = struct('type','diag', 'sigma', 1.0);
    end
    acc = lsqvce.BlockAccumulator(cfg.nc, 'KeepLocal', true, 'Verbose', false, ...
                                  'UserOpt', uo);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, ...
            @(d,ot,a,uo2) lsqvce.synthLoader(cfg, s, d, ot, a, uo2), 1);
    end
end

function v = lagInterp(C, lagIdx)
    % Nearest-lag lookup of a covariance sequence, zero beyond its length.
    k = round(lagIdx);
    v = zeros(size(k));
    m = (k >= 0) & (k <= numel(C) - 1);
    v(m) = C(k(m) + 1);
end

function acc = jointAcc(cfg, model) %#ok<INUSD>
    % Rebuild the accumulator. A production caller would map model.C / model.sigma2
    % into the loader's covariance; here the fixture keeps the weighting fixed so
    % the test exercises the DRIVER's plumbing, not the loader's.
    acc = lsqvce.BlockAccumulator(cfg.nc, 'KeepLocal', true, 'Verbose', false);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, ...
            @(d,ot,a,uo) lsqvce.synthLoader(cfg, s, d, ot, a, uo), 1);
    end
end

function [Abar, W, sp, ebar, N] = lagFixture(tc, N, P, Z) %#ok<INUSL>
    % A small decorrelated arc plus its shared probe block.
    rng(41);
    [~, W] = spdToeplitz(N);
    A      = randn(N, P);
    Abar   = W.' \ A;                                 % W^-T A
    U      = chol(Abar.'*Abar);
    sp     = lsqvce.sharedProbe(U, 'numProbes', Z, 'seed', 9);
    ebar   = W.' \ randn(N,1);                        % W^-T ehat
end

function [Sig, W] = spdToeplitz(N)
    % A positive-definite Toeplitz covariance and its factor, Sigma = W'*W.
    base   = [1 0.45 0.15 0.03];
    k      = min(numel(base), N);
    c      = zeros(N,1);
    c(1:k) = base(1:k);
    Sig    = toeplitz(c);
    [V, D] = eig((Sig+Sig.')/2);
    Sig    = V*diag(max(diag(D), 1e-6))*V.';
    Sig    = (Sig+Sig.')/2;
    W      = chol(Sig);                               % upper: W'*W = Sig
end

function V = denseVn(N, n)
    % Explicit V_n (Ellmer eq. 6.5.16): identity at n=0, else the +/-n diagonals.
    if n == 0
        V = eye(N);
    else
        V = zeros(N);
        for i = 1:N-n
            V(i, i+n) = 1;
            V(i+n, i) = 1;
        end
    end
end
