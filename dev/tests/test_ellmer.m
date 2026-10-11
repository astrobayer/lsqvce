function tests = test_ellmer
%TEST_ELLMER  Unit tests for the Ellmer (2018) stochastic-model primitives.
%
%   Phase A: the type-I DCT (LSQVCE.KERNEL.DCTI), the banded lag components
%   (LSQVCE.ENGINE.COVCOMPONENT), and the Toeplitz cofactor builder
%   (LSQVCE.KERNEL.TOEPLITZCOFACTOR). Later phases (per-lag redundancy, per-frequency
%   PSD VCE, arc-wise factors) append their tests here.
%
%   Every check mirrors an identity verified offline against the exact algebra
%   (Ellmer eqs. 6.5.6, 6.5.12-6.5.17): DCT-I orthonormality/round-trip, the
%   V_n banded apply equals the dense V_n, and sum_n c(n) V_n = toeplitz(c).
%
%   Run with:  results = runtests('tests/test_ellmer.m')
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-34: two tests + fixtures for
%   the keyed arc identity in adjustVCEjoint (F-02) -- induced mid-run
%   non-PD skip with keyed alignment/kept-factor/skip-report assertions,
%   and the vanished-enumeration arcMismatch error. With assistance from
%   Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-40 r2 (F-15): the stabilised
%   alternation test (decay + tail averaging clears the floor on the F-15
%   repro fixture) and the bitwise default regression. With assistance
%   from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-37: three tests for the
%   structured Sigma^-1 sums (dense-reference agreement at 1e-12 on the
%   Toeplitz AND the general fallback path; the long-arc speed lever).
%   With assistance from Claude (Anthropic), model Claude Fable 5.
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

function setupOnce(tc)
    % R7-45: the points below seed the global stream with bare rng(n); the
    % state the suite found is put back when it ends, so later suites do
    % not inherit it.
    prev = rng;
    tc.addTeardown(@() rng(prev));
end

%% ============================================ dctI: orthonormal DCT-I
function test_dctI_ortho_is_self_inverse(tc)
% categories: precision
    % The orthonormal type-I DCT is symmetric and self-inverse: X = X.', X*X = I.
    for N = [2 3 5 8 16 33]
        [~, X] = lsqvce.kernel.dctI(zeros(N,1));         % forward matrix
        % tol: orthonormal-transform identity; error ~ N*eps ~ 1e-14 at N<=33, house bounds add headroom (HANDOVER 6.1.2)
        verifyEqual(tc, X, X.', 'AbsTol', 1e-12, sprintf('N=%d symmetric', N));
        verifyEqual(tc, X*X, eye(N), 'AbsTol', 1e-10, sprintf('N=%d involutive', N));
    end
end

function test_dctI_ortho_round_trip(tc)
% categories: precision
    % ortho is self-inverse, so dctI(dctI(x)) == x with no Inverse flag.
    rng(1);
    x = randn(12, 4);
    y = lsqvce.kernel.dctI(x);
    z = lsqvce.kernel.dctI(y);
    % tol: orthonormal-transform identity; error ~ N*eps ~ 1e-14 at N<=33, house bounds add headroom (HANDOVER 6.1.2)
    verifyEqual(tc, z, x, 'AbsTol', 1e-10);
end

function test_dctI_row_in_row_out(tc)
% categories: precision, vectorisation
    x = 1:6;
    y = lsqvce.kernel.dctI(x);
    verifyTrue(tc, isrow(y), 'a row input must give a row output');
    % tol: orthonormal-transform identity; error ~ N*eps ~ 1e-14 at N<=33, house bounds add headroom (HANDOVER 6.1.2)
    verifyEqual(tc, y(:), lsqvce.kernel.dctI(x(:)), 'AbsTol', 1e-12);
end

function test_dctI_ellmer_norm_invertible(tc)
% categories: precision, reference
    % The Ellmer-literal normalisation is NOT orthonormal but is well conditioned
    % (cond ~ 1.5); its inverse round-trips via 'Inverse',true.
    rng(2);
    x = randn(10, 3);
    [~, X] = lsqvce.kernel.dctI(zeros(10,1), 'Norm', 'ellmer');
    verifyLessThan(tc, cond(X), 5, 'ellmer X must be well conditioned');
    y = lsqvce.kernel.dctI(x, 'Norm', 'ellmer');
    z = lsqvce.kernel.dctI(y, 'Norm', 'ellmer', 'Inverse', true);
    % tol: orthonormal-transform identity; error ~ N*eps ~ 1e-14 at N<=33, house bounds add headroom (HANDOVER 6.1.2)
    verifyEqual(tc, z, x, 'AbsTol', 1e-10, 'ellmer forward then inverse round-trip');
end

function test_dctI_ortho_vs_ellmer_differ_at_nyquist(tc)
% categories: reference
    % The two normalisations differ only in the highest-frequency (Nyquist) bin.
    [~, Xo] = lsqvce.kernel.dctI(zeros(8,1), 'Norm', 'ortho');
    [~, Xe] = lsqvce.kernel.dctI(zeros(8,1), 'Norm', 'ellmer');
    d = abs(Xo - Xe);
    % the difference is concentrated in the last row/column (index N-1)
    interior = d(1:end-1, 1:end-1);
    verifyLessThan(tc, max(interior(:)), 1e-9, ...
        'interior of the two DCT matrices must agree');
    verifyGreaterThan(tc, max(d(:)), 1e-3, ...
        'the Nyquist row/column must differ between normalisations');
end

function test_dctI_too_short_errors(tc)
% categories: contract
    verifyError(tc, @() lsqvce.kernel.dctI(1), 'lsqvce:dctI:tooShort');
end

%% ============================================ covComponent: V_n banded apply
function test_covComponent_n0_is_identity(tc)
% categories: precision
    rng(3);
    x = randn(7, 2);
    % tol: 0 -- an EXACT identity claim (same stored value, no arithmetic)
    verifyEqual(tc, lsqvce.engine.covComponent(x, 0), x, 'AbsTol', 0);
end

function test_covComponent_matches_dense_Vn(tc)
% categories: precision
    % The banded apply equals the explicit V_n * x for every lag.
    N = 6;
    rng(4);
    x = randn(N, 3);
    for n = 0:N-1
        V = denseVn(N, n);
        % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
        verifyEqual(tc, lsqvce.engine.covComponent(x, n), V*x, 'AbsTol', 1e-12, ...
            sprintf('lag n=%d', n));
    end
end

function test_covComponent_sums_to_toeplitz(tc)
% categories: precision, reference
    % sum_n c(n) V_n * x == toeplitz(c) * x  (Ellmer eq. 6.5.17).
    N = 8;
    rng(5);
    c = [1; 0.5; 0.25; 0.1; 0.05; 0.02; 0.01; 0];
    x = randn(N, 4);
    acc = zeros(N, 4);
    for n = 0:N-1
        acc = acc + c(n+1) * lsqvce.engine.covComponent(x, n);
    end
    % tol: dense reconstruction of the same algebra; eps-level accumulation, house bounds (HANDOVER 6.1.2)
    verifyEqual(tc, acc, toeplitz(c)*x, 'AbsTol', 1e-12);
end

function test_covComponent_lag_too_large_errors(tc)
% categories: contract
    verifyError(tc, @() lsqvce.engine.covComponent(ones(4,1), 4), ...
        'lsqvce:covComponent:lagTooLarge');
    verifyError(tc, @() lsqvce.engine.covComponent(ones(4,1), -1), ...
        'lsqvce:covComponent:badLag');
end

%% ============================================ toeplitzCofactor
function test_toeplitzCofactor_dense_is_toeplitz(tc)
% categories: precision
    c = [1; 0.5; 0.2; 0.05];
    Q = lsqvce.kernel.toeplitzCofactor(c);
    % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
    verifyEqual(tc, Q, toeplitz(c), 'AbsTol', 1e-14);
    verifyEqual(tc, Q, Q.', 'AbsTol', 1e-14, 'must be symmetric');
end

function test_toeplitzCofactor_apply_matches_dense(tc)
% categories: precision
    c = [1; 0.5; 0.2; 0.05; 0.01];
    q = lsqvce.kernel.toeplitzCofactor(c, 'Form', 'apply');
    rng(6);
    x = randn(5, 7);
    % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
    verifyEqual(tc, q(x), toeplitz(c)*x, 'AbsTol', 1e-12);
end

function test_toeplitzCofactor_whitens_via_decorr(tc)
% categories: precision
    % A valid PD autocovariance -> chol(Q)^-1 whitens: W Q W' = I.
    c = [1; 0.4; 0.1];               % clearly PD
    Q = lsqvce.kernel.toeplitzCofactor(c);
    W = lsqvce.kernel.decorrFactor(Q);
    % tol: dense reconstruction of the same algebra; eps-level accumulation, house bounds (HANDOVER 6.1.2)
    verifyEqual(tc, W*Q*W.', eye(3), 'AbsTol', 1e-10);
end

function test_toeplitzCofactor_warns_on_bad_acv(tc)
% categories: contract
    % |c(n)| > c(0) is not a valid autocovariance -> a warning (not an error).
    verifyWarning(tc, @() lsqvce.kernel.toeplitzCofactor([1; 2; 0.1]), ...
        'lsqvce:toeplitzCofactor:notACV');
end

%% ============================================ Phase B: shared probe
function test_sharedProbe_shapes_and_solve(tc)
% categories: contract, precision
    % UiZ must satisfy U*UiZ == Zbar, and carry the documented shapes.
    rng(20);
    P = 5;  Z = 30;
    M = randn(P); Nm = M.'*M + P*eye(P);      % SPD
    U = chol(Nm);                             % upper: Nm = U'*U
    sp = lsqvce.engine.sharedProbe(U, 'numProbes', Z, 'seed', 7);
    verifySize(tc, sp.UiZ , [P Z]);
    verifySize(tc, sp.Zbar, [P Z]);
    verifyEqual(tc, sp.probes, Z);
    verifyEqual(tc, sp.P, P);
    % tol: probe construction identity (sign/sqrt bookkeeping, one triangular solve); eps-level, house bounds give headroom
    verifyEqual(tc, U*sp.UiZ, sp.Zbar, 'AbsTol', 1e-10, 'UiZ must solve U*UiZ=Zbar');
end

function test_sharedProbe_rademacher_normalisation(tc)
% categories: precision
    % Zbar = raw/sqrt(Z) with raw entries +-1, so every entry is +-1/sqrt(Z).
    P = 4;  Z = 16;
    U = eye(P);
    sp = lsqvce.engine.sharedProbe(U, 'numProbes', Z, 'seed', 3);
    % tol: probe construction identity (sign/sqrt bookkeeping, one triangular solve); eps-level, house bounds give headroom
    verifyEqual(tc, abs(sp.Zbar), ones(P,Z)/sqrt(Z), 'AbsTol', 1e-12);
end

function test_sharedProbe_seed_reproducible(tc)
% categories: precision | tags: provenance
    U = chol(eye(6)*2);
    a = lsqvce.engine.sharedProbe(U, 'numProbes', 20, 'seed', 11);
    b = lsqvce.engine.sharedProbe(U, 'numProbes', 20, 'seed', 11);
    % tol: 0 -- an EXACT identity claim (same stored value, no arithmetic)
    verifyEqual(tc, a.UiZ, b.UiZ, 'AbsTol', 0, 'a fixed seed must reproduce');
end

function test_sharedProbe_accepts_fixed_probes(tc)
% categories: precision
    seedFor(tc, 228);   % R7-45: drew from the global stream unseeded (suite-order dependent)
    P = 3;  Z = 5;
    U = chol(eye(P)*3);
    raw = sign(randn(P,Z));  raw(raw==0) = 1;
    sp = lsqvce.engine.sharedProbe(U, 'probes', raw);
    % tol: probe construction identity (sign/sqrt bookkeeping, one triangular solve); eps-level, house bounds give headroom
    verifyEqual(tc, sp.Zbar, raw/sqrt(Z), 'AbsTol', 1e-12);
end

function test_sharedProbe_rejects_nonsquare(tc)
% categories: contract
    verifyError(tc, @() lsqvce.engine.sharedProbe(ones(3,4)), ...
        'lsqvce:sharedProbe:notSquare');
end

%% ============================================ Phase B: per-lag redundancy
function test_redundancyLags_s_matches_dense_trace(tc)
% categories: precision
    % The memory-honest sub-diagonal sums must equal trace(Rbar*V_n) exactly for
    % the SAME Rbar (this is an algebraic identity, not a Monte-Carlo question).
    [Abar, W, sp, ebar, N] = lagFixture(tc, 12, 4, 200);
    nLag = 6;
    r = lsqvce.kernel.redundancyLags(Abar, W, sp.UiZ, ebar, nLag, 'WantRbar', true);
    for n = 0:nLag-1
        Vn = denseVn(N, n);
        % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
        verifyEqual(tc, r.s(n+1), trace(r.Rbar*Vn), 'RelTol', 1e-8, ...
            sprintf('s_%d must equal trace(Rbar V_n)', n));
    end
end

function test_redundancyLags_omega_matches_quadratic(tc)
% categories: precision
    % Omega_n == etilde' V_n etilde, with etilde = Sigma^-1 ehat.
    [Abar, W, sp, ebar, N] = lagFixture(tc, 11, 3, 100);
    nLag = 5;
    r  = lsqvce.kernel.redundancyLags(Abar, W, sp.UiZ, ebar, nLag);
    et = W \ ebar;                                   % Sigma^-1 ehat
    for n = 0:nLag-1
        Vn = denseVn(N, n);
        % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
        verifyEqual(tc, r.omega(n+1), et.'*Vn*et, 'RelTol', 1e-8, ...
            sprintf('Omega_%d must equal etilde'' V_n etilde', n));
    end
end

function test_redundancyLags_Rbar_converges_to_exact(tc)
% categories: precision
    % With many probes, Rbar -> Sigma^-1 - Sigma^-1 A N^-1 A' Sigma^-1.
    rng(31);
    N = 10;  P = 3;
    [Sig, W] = spdToeplitz(N);
    A    = randn(N, P);
    Abar = W.' \ A;                                   % W^-T A
    Nm   = Abar.'*Abar;
    U    = chol(Nm);
    sp   = lsqvce.engine.sharedProbe(U, 'numProbes', 20000, 'seed', 5);
    r    = lsqvce.kernel.redundancyLags(Abar, W, sp.UiZ, [], 3, 'WantRbar', true);
    Sinv    = Sig \ eye(N);
    Rexact  = Sinv - Sinv*A*(Nm\(A.'*Sinv));
    % tol: Monte-Carlo probe scatter ~ 1/sqrt(Probes) on this fixture; 4/sqrt(Z)-style bar (HANDOVER 6.1.2), not containment
    verifyEqual(tc, r.Rbar, Rexact, 'AbsTol', 0.05, ...
        'Rbar must approach the exact redundancy matrix');
end

function test_redundancyLags_omega_skipped_when_no_residual(tc)
% categories: robustness
    [Abar, W, sp, ~, ~] = lagFixture(tc, 8, 2, 50);
    r = lsqvce.kernel.redundancyLags(Abar, W, sp.UiZ, [], 4);
    verifyTrue(tc, all(isnan(r.omega)), 'empty ebar must give NaN omega');
    verifyTrue(tc, all(isfinite(r.s)), 's must still be computed');
end

function test_redundancyLags_rejects_bad_shapes(tc)
% categories: contract
    [Abar, W, sp, ebar, ~] = lagFixture(tc, 9, 3, 40);
    verifyError(tc, @() lsqvce.kernel.redundancyLags(Abar, W(1:3,1:3), sp.UiZ, ebar, 3), ...
        'lsqvce:redundancyLags:badW');
    verifyError(tc, @() lsqvce.kernel.redundancyLags(Abar, W, sp.UiZ(1:2,:), ebar, 3), ...
        'lsqvce:redundancyLags:badUiZ');
    % WP-39 b1: the scalar constraint moved INTO the typed signature --
    % nMaxLag <= 0 now fails mustBePositive BEFORE the body (the
    % cross-argument shape checks above keep their lsqvce: ids in the
    % body). Contract change release-noted.
    verifyError(tc, @() lsqvce.kernel.redundancyLags(Abar, W, sp.UiZ, ebar, 0), ...
        'MATLAB:validators:mustBePositive');
end

function test_redundancyLags_nugget_matches_the_general_path(tc)
% categories: precision, contract
    % R7 WP-115 part 3b: Sigma = T + diag(c), a stationary law plus the joint
    % estimator's rejected-epoch nuggets. 'Nugget' / 'LawW' take Sigma^-1's
    % sub-diagonal sums from T's Toeplitz path and a Woodbury over the nugget
    % rows; the general O(N^3) fallback on the SAME W is the reference --
    % including a deletion-sized nugget (1e12 T_ii) next to a partial one.
    rs = RandStream('mt19937ar', 'Seed', 115);
    N = 80;  P = 3;  nLag = 7;
    T  = toeplitz(0.8 .^ (0:N-1));
    WT = chol(T);
    c  = zeros(N, 1);  c([5 17 18 40 77]) = [0.3 2 1e12 0.7 1e12];
    W  = chol(T + diag(c));
    Abar = W.' \ randn(rs, N, P);
    sp   = lsqvce.engine.sharedProbe(chol(Abar.'*Abar), 'numProbes', 40, 'seed', 9);
    ebar = W.' \ randn(rs, N, 1);
    rN = lsqvce.kernel.redundancyLags(Abar, W, sp.UiZ, ebar, nLag, 'Nugget', c, 'LawW', WT);
    rG = lsqvce.kernel.redundancyLags(Abar, W, sp.UiZ, ebar, nLag);
    % tol: two exact routes to the same sums (Woodbury + Gohberg-Semencul vs
    % N triangular solves); the mirror agrees to 1e-16, 1e-9 of s_0 leaves
    % eps*cond headroom for the 1e12 nugget (house bounds, HANDOVER 6.1.2)
    verifyEqual(tc, rN.s, rG.s, 'AbsTol', 1e-9 * abs(rG.s(1)), 'the nugget sums equal the general path');
    % tol: exact -- Omega does not touch the sums, the same arithmetic on both calls
    verifyEqual(tc, rN.omega, rG.omega, 'Omega is the same computation');
    verifyError(tc, @() lsqvce.kernel.redundancyLags(Abar, W, sp.UiZ, ebar, nLag, 'Nugget', c(1:end-1), ...
        'LawW', WT), 'lsqvce:redundancyLags:badNugget');
    verifyError(tc, @() lsqvce.kernel.redundancyLags(Abar, W, sp.UiZ, ebar, nLag, 'Nugget', c), ...
        'lsqvce:redundancyLags:badNugget');
end

%% ============================== Phase C: per-bin spectral VCE (positive components)
function test_psdVCE_identity_when_Q_equals_T(tc)
% categories: precision
    % If Q == T then alpha2 == 1 in every bin, so the spectrum is unchanged:
    % the estimator is at a fixed point and reports convergence.
    rs = RandStream('mt19937ar', 'Seed', 50);
    L = 10;  M = 4;
    S  = 0.5 + rand(rs, L, 1);
    q  = rand(rs, L, M) + 0.5;
    p  = lsqvce.engine.psdVCE(q, q, S);
    % tol: fixed-point identity Q == T -> alpha2 = 1; fp division noise only, house exact-identity bound
    verifyEqual(tc, p.alpha2, ones(L,1), 'AbsTol', 1e-10, 'alpha2 must be 1');
    % tol: as above
    verifyEqual(tc, p.Shat, S, 'RelTol', 1e-10, 'the spectrum must not move');
    % tol: as above
    verifyEqual(tc, p.maxDev, 0, 'AbsTol', 1e-10);
end
function test_psdVCE_factors_are_nonnegative_by_construction(tc)
% categories: precision, reference
    % R7 WP-111 part 3 (R7-02): THE property the positive components exist
    % for. Every bin's component toeplitz(K(:, j)) is PSD, so for ANY residual
    % e and ANY PSD redundancy matrix R the bin statistics K' Omega = e' V_j e
    % and K' s = tr(R V_j) are >= 0, and so is alpha2 -- no floor needed. The
    % DCT components were indefinite: the same draw gave negative factors.
    rs = RandStream('mt19937ar', 'Seed', 51);
    m = 48;  L = 9;  M = 5;
    K = lsqvce.kernel.hatACV(L, m);
    Q = zeros(L, M);  T = zeros(L, M);
    for a = 1:M
        e = randn(rs, m, 1);  B = randn(rs, m, m - 4);
        Q(:, a) = K.' * lagSumsLocal(e * e.');
        T(:, a) = K.' * lagSumsLocal(B * B.');
    end
    % tol: 1e-12 relative to the largest statistic -- rounding of a PSD quadratic form
    verifyGreaterThanOrEqual(tc, min([Q(:); T(:)]), -1e-12 * max(abs([Q(:); T(:)])), ...
        'every bin statistic is a PSD quadratic form');
    p = lsqvce.engine.psdVCE(Q, T, ones(L, 1));
    verifyTrue(tc, all(p.alpha2 >= 0) && all(p.Shat > 0), 'factors >= 0, spectrum > 0');
end

function test_hatACV_columns_are_psd_and_exact(tc)
% categories: reference, precision, metamorphic
    % R7 WP-111 part 3: column j of hatACV is the exact ACV of the hat spectrum
    % on bin j (two-sided, unit sample), checked against adaptive quadrature;
    % its Toeplitz matrix is PSD; a flat spectrum is white noise exactly; and
    % the ACV has no lag period (R7-01): a longer arc only APPENDS lags.
    L = 7;  m = 40;
    K = lsqvce.kernel.hatACV(L, m);
    d = 1 / (2*(L-1));
    % tol: quadrature settings -- 1e-14 / 1e-12, two digits below the assertion's bar
    qOpt = {'AbsTol', 1e-14, 'RelTol', 1e-12};
    for j = [1 4 L]
        fj = (j-1) * d;
        hat = @(f) max(0, 1 - abs(f - fj) / d);
        for tau = [0 1 5 23]
            ref = 2 * integral(@(f) hat(f) .* cos(2*pi*f*tau), 0, 0.5, 'Waypoints', fj, qOpt{:});
            % tol: 1e-10 absolute -- the quadrature's own accuracy on a kinked integrand
            verifyEqual(tc, K(tau+1, j), ref, 'AbsTol', 1e-10, sprintf('bin %d, lag %d', j, tau));
        end
        % tol: 1e-12 relative to r(0) -- eigenvalues of an exact PSD Toeplitz, rounding only
        verifyGreaterThanOrEqual(tc, min(eig(toeplitz(K(:, j)))), -1e-12 * K(1, j), ...
            sprintf('bin %d: its component is PSD', j));
    end
    % tol: 1e-12 absolute -- a flat spectrum interpolates exactly: r = [1; 0; ...]
    verifyEqual(tc, K * ones(L, 1), [1; zeros(m-1, 1)], 'AbsTol', 1e-12, 'flat = white');
    K2 = lsqvce.kernel.hatACV(L, 3*m + 1);
    % tol: exact -- the same closed form at the same lags (no resampling, no period)
    verifyEqual(tc, K2(1:m, :), K, 'AbsTol', 0, 'a longer arc only appends lags');
    % the projection form (fold modulo 2(L-1) + one FFT) is K' * X without K
    X = randn(RandStream('mt19937ar', 'Seed', 4), 3*m + 1, 2);
    % tol: 1e-13 of sum |X| -- two summation orders (FFT vs matrix product), rounding only
    verifyEqual(tc, lsqvce.kernel.hatACV(L, 3*m + 1, X), K2.' * X, 'AbsTol', 1e-13 * sum(abs(X(:))));
    verifyError(tc, @() lsqvce.kernel.hatACV(L, m, X), 'lsqvce:hatACV:badX');
end
function test_arcVarianceFactors_clamps_negative_sigma2(tc)
% categories: contract
    % F-05 pin: a constructed negative alpha2 (one arc's Omega column
    % negative) must NOT leave as a negative variance factor: it is clamped
    % at ClampEps, warned about, and flagged in .clamped. Pre-WP-35 the
    % factor left negative and .clamped did not exist -- both assertions
    % fail there (§6.1).
    omegaT = [ones(3,1), ones(3,1), -0.4*ones(3,1)];   % arc 3: alpha2 < 0
    sT     = ones(3, 3);
    sig2   = ones(3, 1);
    muteWarn(tc, 'lsqvce:arcVarianceFactors:clamped');   % provoked on purpose (arc 3)
    a = lsqvce.engine.arcVarianceFactors(omegaT, sT, sig2, 'MeanType', 'plain');
    verifyTrue(tc, a.clamped(3),  'the negative-alpha2 arc must be flagged');
    verifyFalse(tc, any(a.clamped(1:2)), 'healthy arcs must not be flagged');
    verifyEqual(tc, a.sigma2(3), 1e-8, ...
        'the clamped factor must equal the documented default ClampEps');
    verifyTrue(tc, all(a.sigma2 > 0), 'no negative factor may leave');
end

function test_arcVarianceFactors_noop_on_healthy(tc)
% categories: precision, reference
    % F-05 regression: on healthy input the clamp is a NO-OP -- .clamped is
    % all-false and sigma2 matches the closed formula gamma2 / gammaBar^2
    % (eqs. 6.5.28-6.5.30, plain mean) within summation-order ULPs; the
    % clamp branch itself never executes, so the pre-fix numbers are
    % unchanged.
    rng(35);
    M = 6;  nLag = 4;
    alpha2 = 0.5 + rand(M, 1);                       % strictly positive
    sT     = ones(nLag, M);
    omegaT = repmat(alpha2.', nLag, 1);              % per-arc ratio = alpha2
    sig2   = 0.8 + 0.4*rand(M, 1);
    a = lsqvce.engine.arcVarianceFactors(omegaT, sT, sig2, 'MeanType', 'plain');
    verifyFalse(tc, any(a.clamped), 'healthy input must not clamp');
    gamma2   = alpha2 .* sig2;
    gammaBar = mean(sqrt(abs(gamma2)));
    % tol: fixed-point identity om==s -> alpha2=1 by the update law; fp division noise only, house exact-identity bound
    verifyEqual(tc, a.sigma2, gamma2 / gammaBar^2, 'RelTol', 4*eps, ...
        ['converged results on healthy inputs must match the closed ' ...
         'formula (4 eps: summation-order ULPs, CI-measured 2.6e-16)']);
end

function test_psdVCE_freeze_redundancy_holds_weak_bins(tc)
% categories: contract, robustness
    % R7-22: Options.Noise.FreezeRedundancy was documented ("a frequency whose
    % redundancy is below this holds its value") and reached nothing. psdVCE
    % holds a bin whose redundancy Sprior_j * sum_m T(j, m) is below the
    % threshold: every bin's ratio is 2, bin 3's redundancy is 0.5.
    T = [10; 10; 0.5; 10];  Q = 2 * T;  S = ones(4, 1);
    free = lsqvce.engine.psdVCE(Q, T, S);
    held = lsqvce.engine.psdVCE(Q, T, S, 'FreezeRedundancy', 1);
    % tol: 1e-12 absolute -- one division of O(10) numbers
    tc.verifyEqual(free.alpha2, [2; 2; 2; 2], 'AbsTol', 1e-12, 'without a freeze every bin updates');
    % tol: as above
    tc.verifyEqual(held.alpha2, [2; 2; 1; 2], 'AbsTol', 1e-12, 'the weak bin holds');
    tc.verifyEqual(held.frozen, [false; false; true; false]);
    % the redundancy, not the bare denominator: the same bin with a 4x prior is free
    tc.verifyFalse(any(lsqvce.engine.psdVCE(Q, T, [1; 1; 4; 1], 'FreezeRedundancy', 1).frozen));
    % and the Options knob reaches the joint driver (it was forwarded nowhere)
    a = lsqvce.Options(Noise=struct('FreezeRedundancy', 100, 'Smooth', 2)).engineArgs('adjustVCEjoint');
    k = find(strcmp(a.nv(1:2:end), 'FreezeRedundancy'));
    tc.verifyEqual(a.nv{2*k}, 100);
    k = find(strcmp(a.nv(1:2:end), 'Smooth'));
    tc.verifyEqual(a.nv{2*k}, 2);
end
function test_psdVCE_rejects_bad_shapes(tc)
% categories: contract
    S = ones(6, 1);
    verifyError(tc, @() lsqvce.engine.psdVCE(ones(5,3), ones(5,3), S), ...
        'lsqvce:psdVCE:shapeMismatch');
    verifyError(tc, @() lsqvce.engine.psdVCE(ones(6,3), ones(6,4), S), ...
        'lsqvce:psdVCE:arcMismatch');
end
function test_arcVarianceFactors_plain_mean_normalises_to_one(tc)
% categories: precision, reference
    % With the PLAIN mean (eq. 6.5.29) the mean of the arc FACTORS (not their
    % squares) is exactly 1 -- Ellmer's stated normalisation property.
    rng(60);
    N = 8;  M = 12;
    om = rand(N,M) + 0.5;
    s  = rand(N,M) + 0.5;
    a  = lsqvce.engine.arcVarianceFactors(om, s, ones(M,1), 'MeanType', 'plain');
    % tol: exact-path identity (same numbers recombined); house round-trip bound 1e-12/1e-10 (HANDOVER 6.1.2)
    verifyEqual(tc, mean(sqrt(abs(a.sigma2))), 1, 'AbsTol', 1e-10, ...
        'the plain mean must normalise mean(gamma_m/gammaBar) to 1');
end

function test_arcVarianceFactors_robust_resists_outliers(tc)
% categories: contract
    % The robust central-half mean must sit below the outlier-inflated plain mean
    % (that is the entire reason GROOPS uses it).
    rng(61);
    N = 6;  M = 40;
    om = rand(N,M) + 0.5;
    s  = rand(N,M) + 0.5;
    om(:, 1:3) = om(:, 1:3) * 400;             % a few pathological arcs
    ap = lsqvce.engine.arcVarianceFactors(om, s, ones(M,1), 'MeanType', 'plain');
    ar = lsqvce.engine.arcVarianceFactors(om, s, ones(M,1), 'MeanType', 'robust');
    verifyLessThan(tc, ar.gammaBar, ap.gammaBar, ...
        'the robust mean must sit below the outlier-inflated plain mean');
end

function test_arcVarianceFactors_gamma_detour_shapes(tc)
% categories: contract, precision
    rng(62);
    N = 7;  M = 9;
    sig2 = rand(M,1) + 0.5;
    a = lsqvce.engine.arcVarianceFactors(rand(N,M)+0.5, rand(N,M)+0.5, sig2);
    verifySize(tc, a.sigma2, [M 1]);
    verifySize(tc, a.alpha2, [M 1]);
    % tol: fixed-point identity om==s -> alpha2=1 by the update law; fp division noise only, house exact-identity bound
    verifyEqual(tc, a.gamma2, a.alpha2 .* sig2, 'RelTol', 1e-10, ...
        'gamma2 must be alpha2 .* sigma2 (eq. 6.5.28)');
    verifyTrue(tc, all(a.sigma2 >= 0), 'variance factors must be non-negative');
end

function test_arcVarianceFactors_small_M_falls_back(tc)
% categories: robustness
    % With too few arcs to trim, the robust mean must degrade gracefully.
    rng(63);
    N = 5;
    for M = [1 2 3]
        a = lsqvce.engine.arcVarianceFactors(rand(N,M)+0.5, rand(N,M)+0.5, ones(M,1));
        verifyTrue(tc, isfinite(a.gammaBar) && a.gammaBar > 0, ...
            sprintf('M=%d must still give a usable mean', M));
    end
end

function test_arcVarianceFactors_rejects_arc_mismatch(tc)
% categories: contract
    seedFor(tc, 527);   % R7-45: drew from the global stream unseeded (suite-order dependent)
    verifyError(tc, ...
        @() lsqvce.engine.arcVarianceFactors(rand(5,3), rand(5,3), ones(4,1)), ...
        'lsqvce:arcVarianceFactors:arcMismatch');
end

function test_arcVarianceFactors_alpha2_is_blind_to_the_arcs_own_factor(tc)
% categories: metamorphic, precision
    % R6 WP-98: THE INVARIANT THE F-01 FIX RESTS ON, pinned at machine
    % precision. adjustVCEjoint>sweepArcs hands an arc's numerator AND its
    % denominator scaled by the same arc variance, and alpha2_m is their ratio
    % through a LINEAR map (a column sum), so the factor divides out: alpha2 is a
    % function of the data and the law only and can NEVER respond to the factor
    % being estimated. That is why eq. 6.5.28 must multiply by the factor the
    % WEIGHTING used, not by the per-arc one -- multiplying by a number the
    % estimate is blind to made the update a geometric recursion with no fixed
    % point (finding F-01, HANDOVER 12.21 and 12.22).
    % Two related calls, no oracle needed: scale both inputs per column and the
    % raw alpha2 must not move.
    rs = RandStream('mt19937ar', 'Seed', 7);
    n = 61;  M = 6;
    om = rand(rs, n, M);
    sn = rand(rs, n, M) + 1;                       % strictly positive denominators
    s2 = [1 9 1 4 0.25 2].';                       % wildly different current factors
    a0 = lsqvce.engine.arcVarianceFactors(om,          sn,          ones(M, 1));
    a1 = lsqvce.engine.arcVarianceFactors(om .* s2.',  sn .* s2.',  s2);
    % tol: 1e-12 relative -- the cancellation is exact in exact arithmetic, so
    % what is left is the rounding of two column sums
    tc.verifyEqual(a1.alpha2, a0.alpha2, 'RelTol', 1e-12, ...
                   'alpha2 does not see the arc factor it is estimating');
    % and eq. 6.5.28 therefore multiplies a CONSTANT by whatever it is handed
    % tol: 1e-12 relative, the same bound and the same reason as above -- the
    % product is formed from the alpha2 the previous assertion just pinned
    tc.verifyEqual(a1.gamma2, a0.alpha2 .* s2, 'RelTol', 1e-12, ...
                   'gamma2 = alpha2 * whatever factor is passed in');
end

%% ============================================ Phase E: joint loop driver
function test_adjustVCEjoint_runs_and_returns_shapes(tc)
% categories: contract, robustness
    % End-to-end: the ten-step loop must run on a real accumulator and return a
    % PSD, a covariance function, and one variance factor per arc.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % provoked on purpose here (WP-47, gate rule 7)
    cfg  = jointCfg();
    nLag = 8;
    r = lsqvce.engine.adjustVCEjoint(@(model) jointAcc(cfg, model), ...
            'NumLags', nLag, 'MaxIter', 2, 'Probes', 40, 'Seed', 3, ...
            'Verbose', false);
    verifySize(tc, r.C, [nLag 1]);
    verifySize(tc, r.S, [nLag 1]);
    verifyTrue(tc, all(r.S > 0), 'the spectrum must stay strictly positive');
    % tol: 1e-12 -- the law is a shape by construction, r(0) = K(1, :) S = 1
    verifyEqual(tc, r.C(1), 1, 'AbsTol', 1e-12, 'the returned law is normalised');
    % tol: 1e-12 relative -- C is DERIVED from S, the same product evaluated once more
    verifyEqual(tc, r.C, lsqvce.kernel.hatACV(nLag, nLag) * r.S, 'RelTol', 1e-12);
    verifyTrue(tc, numel(r.sigma2) >= 1, 'at least one arc factor');
    verifyTrue(tc, all(isfinite(r.sigma2)), 'arc factors must be finite');
    verifyClass(tc, r.acc, 'lsqvce.engine.BlockAccumulator');
end

function test_adjustVCEjoint_history_grows_per_iteration(tc)
% categories: contract, robustness
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % provoked on purpose here (WP-47, gate rule 7)
    cfg  = jointCfg();
    nLag = 6;
    r = lsqvce.engine.adjustVCEjoint(@(model) jointAcc(cfg, model), ...
            'NumLags', nLag, 'MaxIter', 3, 'Probes', 30, 'Seed', 5, ...
            'Tol', 1e-12);                       % force all iterations
    verifyEqual(tc, numel(r.history.maxDev), r.nIter);
    verifySize(tc, r.history.S, [r.nIter nLag]);
    verifyTrue(tc, all(isfinite(r.history.maxDev)), ...
        'the convergence measure must be finite every pass');
end

function test_adjustVCEjoint_model_reaches_accFcn(tc)
% categories: contract
    % The driver must hand the CURRENT model to accFcn each pass (steps 1-3),
    % so the caller can rebuild Sigma from it. containers.Map is a HANDLE, so
    % the closure below can record into it without a nested function.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % provoked on purpose here (WP-47, gate rule 7)
    cfg = jointCfg();
    rec = containers.Map('KeyType', 'double', 'ValueType', 'any');
    accF = @(model) recordAndBuild(rec, cfg, model);
    lsqvce.engine.adjustVCEjoint(accF, 'NumLags', 5, 'MaxIter', 2, ...
                          'Probes', 20, 'Seed', 7, 'Tol', 1e-12);
    verifyGreaterThanOrEqual(tc, rec.Count, 2, 'accFcn called once per pass');
    verifyTrue(tc, rec.isKey(1) && rec.isKey(2), 'iterations 1 and 2 recorded');
    m1 = rec(1);
    verifyTrue(tc, isfield(m1, 'S') && isfield(m1, 'C') && isfield(m1, 'iter'), ...
        'the model struct must carry S, C and iter');
    verifyEqual(tc, m1.iter, 1);
end

function test_adjustVCEjoint_default_start_is_white(tc)
% categories: contract, precision
    % R7 WP-111 part 3: Ellmer's white start (eq. 6.5.62) is the default and
    % is REGULAR under positive components (the DCT weighting made it a
    % degenerate fixed point, and the default had to be [1; 0.3; ...]). Pass 1
    % receives a flat spectrum and r = [1; 0; ...]; an explicit S0 arrives
    % normalised to r(0) = 1.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % provoked on purpose here (WP-47, gate rule 7)
    cfg = jointCfg();
    rec = containers.Map('KeyType', 'double', 'ValueType', 'any');
    lsqvce.engine.adjustVCEjoint(@(model) recordAndBuild(rec, cfg, model), 'NumLags', 5, ...
                                 'MaxIter', 1, 'Probes', 20, 'Seed', 9);
    m1 = rec(1);
    % tol: 1e-12 -- a flat spectrum interpolates exactly to white noise
    verifyEqual(tc, m1.C, [1; 0; 0; 0; 0], 'AbsTol', 1e-12, 'the default start is white');
    % tol: 1e-12 relative -- one normalisation of equal numbers
    verifyEqual(tc, m1.S / m1.S(1), ones(5, 1), 'RelTol', 1e-12);
    rec2 = containers.Map('KeyType', 'double', 'ValueType', 'any');
    lsqvce.engine.adjustVCEjoint(@(model) recordAndBuild(rec2, cfg, model), 'NumLags', 5, ...
                                 'MaxIter', 1, 'Probes', 20, 'Seed', 9, 'S0', [5; 4; 3; 2; 1]);
    m2 = rec2(1);
    % tol: 1e-12 -- the normalisation r(0) = 1 of the given S0
    verifyEqual(tc, m2.C(1), 1, 'AbsTol', 1e-12, 'S0 arrives as a shape');
    % tol: 1e-12 relative -- the same shape, rescaled
    verifyEqual(tc, m2.S / m2.S(end), [5; 4; 3; 2; 1], 'RelTol', 1e-12);
end
function test_adjustVCEjoint_rejects_bad_S0(tc)
% categories: contract
    cfg = jointCfg();
    for S0 = {ones(3,1), [1; -1; 1; 1; 1; 1; 1; 1], zeros(8,1), [1; NaN; 1; 1; 1; 1; 1; 1]}
        verifyError(tc, ...
            @() lsqvce.engine.adjustVCEjoint(@(model) jointAcc(cfg, model), ...
                    'NumLags', 8, 'S0', S0{1}), ...
            'lsqvce:adjustVCEjoint:badS0');
    end
end
function test_redundancyLags_fft_lag_sums_match_the_loop_form(tc)
% categories: precision, vectorisation
    % R5 M3 3c: all per-lag sub-diagonal sums come from one FFT per column
    % (autocorrelation), replacing the O(N nLag) loop. Pinned against the
    % explicit loop on a dense N x C block with a local block present (both
    % Rt and K paths) and on the residual (Omega), to rounding.
    rs = RandStream('mt19937ar', 'Seed', 33);
    N = 257;  Z = 7;  nl = 3;  nL = 40;
    A = randn(rs, N, 5);  B = randn(rs, N, nl);
    ac = 0.8 .^ (0:N-1);  Sig = toeplitz(ac);  W = chol(Sig);
    Abar = W.' \ A;  Bbar = W.' \ B;
    UiZ  = randn(rs, 5, Z) / sqrt(Z);
    Ll   = chol(Bbar.' * Bbar);
    ebar = randn(rs, N, 1);
    out = lsqvce.kernel.redundancyLags(Abar, W, UiZ, ebar, nL, 'LocalBbar', Bbar, 'LocalLl', Ll, 'WantRbar', true);
    % loop form from the dense Rbar the test path returns (WantRbar)
    sLoop = zeros(nL, 1);  sLoop(1) = trace(out.Rbar);
    for n = 1:nL-1, sLoop(n+1) = 2 * sum(diag(out.Rbar, -n)); end
    et = W \ ebar;  oLoop = zeros(nL, 1);  oLoop(1) = et.' * et;
    for n = 1:nL-1, oLoop(n+1) = 2 * (et(n+1:N).' * et(1:N-n)); end
    % tol: 1e-10 relative -- FFT autocorrelation vs direct sums on O(1)
    % numbers of length 257 (a few hundred ulps of accumulated rounding)
    tc.verifyEqual(out.s, sLoop, 'RelTol', 1e-10, 'AbsTol', 1e-12);
    tc.verifyEqual(out.omega, oLoop, 'RelTol', 1e-10, 'AbsTol', 1e-12);
end

%% ============================================ R5 M3 3d: per-stream policy
function test_joint_policy_none_and_scale_and_pool(tc)
% categories: contract, metamorphic
    % R5 M3 3d (WP-83, R5-D4): the joint driver takes a per-stream POLICY.
    % Fixture: one system (91), one stream 'pos' -> key 'x91|pos'.
    %   'none'  : the stream's arc factors stay 1 and jointModelOpt emits NO
    %             override for it (untouched);
    %   'scale' : factors are estimated, the override carries the level ONLY;
    %   Pool "stream": out.Sgroups holds one law per stream and out.S is it.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');
    cfg = jointCfg();
    key = 'x91|pos';
    E = containers.Map({key}, {'none'});
    rN = lsqvce.engine.adjustVCEjoint(@(model) jointAcc(cfg, model), 'NumLags', 6, ...
             'MaxIter', 2, 'Probes', 80, 'Seed', 3, 'Verbose', false, 'Estimate', E);
    % tol: exact -- no factor is ever written back under 'none'
    tc.verifyEqual(rN.sigma2, ones(size(rN.sigma2)), 'none: factors untouched');
    uoN = lsqvce.engine.jointModelOpt(struct('S', rN.S, 'sigma2', rN.sigma2, 'iter', 2, ...
             'Estimate', E), {key});
    tc.verifyFalse(isfield(uoN.covOverride.bySys, 'x91'), 'none: no override at all');

    E = containers.Map({key}, {'scale'});
    rS = lsqvce.engine.adjustVCEjoint(@(model) jointAcc(cfg, model), 'NumLags', 6, ...
             'MaxIter', 2, 'Probes', 80, 'Seed', 3, 'Verbose', false, 'Estimate', E);
    tc.verifyTrue(any(rS.sigma2 ~= 1), 'scale: factors are estimated');
    uoS = lsqvce.engine.jointModelOpt(struct('S', rS.S, 'sigma2', rS.sigma2, 'iter', 2, ...
             'Estimate', E, 'streamSigma2', containers.Map({key}, {mean(rS.sigma2)})), {key});
    sp = uoS.covOverride.bySys.x91.pos;
    tc.verifyEqual(fieldnames(sp), {'sigmaScale'}, 'scale: the level only, on the declared law');
    % tol: value plumbing identity (sqrt of the mean factor); house round-trip bound 1e-12
    tc.verifyEqual(sp.sigmaScale, sqrt(mean(rS.sigma2)), 'AbsTol', 1e-12);

    rP = lsqvce.engine.adjustVCEjoint(@(model) jointAcc(cfg, model), 'NumLags', 6, ...
             'MaxIter', 2, 'Probes', 80, 'Seed', 3, 'Verbose', false, 'Pool', "stream");
    tc.verifyEqual(keys(rP.Sgroups), {key}, 'one law per stream');
    tc.verifyEqual(keys(rP.Cgroups), {key});
    % tol: exact -- out.S IS the single group's S
    tc.verifyEqual(rP.S, rP.Sgroups(key));
    % parity: with one stream, Pool "stream" and Pool "all" estimate the same law
    rA = lsqvce.engine.adjustVCEjoint(@(model) jointAcc(cfg, model), 'NumLags', 6, ...
             'MaxIter', 2, 'Probes', 80, 'Seed', 3, 'Verbose', false);
    % tol: exact -- same arcs, same probes (seeded), same arithmetic on one group
    tc.verifyEqual(rP.S, rA.S);
end

function test_jointModelOpt_takes_sqrt_of_variance(tc)
% categories: precision
    % CONVENTION 1: model.sigma2 is a VARIANCE, the level a covariance spec
    % carries is a STANDARD DEVIATION. R5 M3 3b: the level is RELATIVE to the
    % stream's declared sigma (sigmaScale, applied by applyCovOverride), so a
    % mean variance factor of 4 arrives as sigmaScale 2 -- never as sigma.
    model = struct('S', [1; 0.4; 0.1], 'sigma2', [4; 4; 4], 'iter', 2);
    uo = lsqvce.engine.jointModelOpt(model, {'x1|pos'});
    spec = uo.covOverride.bySys.x1.pos;
    % tol: value plumbing identity (stored numbers surfaced unchanged); house round-trip bound 1e-12
    verifyEqual(tc, spec.sigmaScale, 2, 'AbsTol', 1e-12, ...
        'the relative level must be sqrt(mean variance), not the variance itself');
    verifyFalse(tc, isfield(spec, 'sigma'), 'the declared sigma of the stream is never overridden');
    % a per-stream factor wins over the mean when the driver supplies it
    m2 = model;  m2.streamSigma2 = containers.Map({'x1|pos'}, {9});
    uo2 = lsqvce.engine.jointModelOpt(m2, {'x1|pos', 'x1|kbr'});
    % tol: value plumbing identity (sqrt of a stored number); house round-trip bound 1e-12
    verifyEqual(tc, uo2.covOverride.bySys.x1.pos.sigmaScale, 3, 'AbsTol', 1e-12);
    % tol: as above
    verifyEqual(tc, uo2.covOverride.bySys.x1.kbr.sigmaScale, 2, 'AbsTol', 1e-12, 'a stream without its own factor gets the mean');
end
function test_jointModelOpt_emits_a_psd_on_normalised_frequency(tc)
% categories: contract, precision
    % R5 M3 3b (H-17) / R7 WP-111 part 3: the joint law leaves as a PSD on
    % NORMALISED frequency (the bins are per SAMPLE; covModel places them on
    % the arc's dt): ordinates 2 S on fNorm = f_j, scaled to r(0) = 1, and
    % its autocovariance IS the hat-basis law K S -- exactly, at every lag
    % and for every dt (the 3b cosine transform of a truncated ACF matched
    % only the shape, to 5e-2).
    S = [3; 2; 0.5; 0.1; 0.02];  L = numel(S);
    uo = lsqvce.engine.jointModelOpt(struct('S', S, 'sigma2', [], 'iter', 1), {'x1|kbr'});
    spec = uo.covOverride.bySys.x1.kbr;
    tc.verifyEqual(spec.type, 'psd');
    tc.verifyEqual(spec.fNorm(1), 0);
    % tol: exact -- the last grid point IS the normalised Nyquist by construction
    tc.verifyEqual(spec.fNorm(end), 0.5);
    m = 30;
    K = lsqvce.kernel.hatACV(L, m);
    rLaw = K * S / (K(1, :) * S);
    for dt = [1 5]
        % on the physical axis f = fNorm/dt the density is S*dt (covModel does this)
        r = lsqvce.kernel.psd2acv(spec.fNorm / dt, spec.S * dt, dt, m);
        % tol: 1e-10 absolute on r(0) = 1 -- psd2acv's small-argument series (x < 1e-2) is
        % accurate to x^4 / 280 relative, the rest is rounding
        tc.verifyEqual(r(:), rLaw, 'AbsTol', 1e-10, sprintf('the law K S at dt = %g', dt));
    end
end
function test_jointModelOpt_places_override_bySys(tc)
% categories: contract
    % CONVENTION 2: loaders read opt.covOverride.bySys.(satKey).(obsType).
    model = struct('S', [1; 0.3], 'sigma2', [], 'iter', 1);
    uo = lsqvce.engine.jointModelOpt(model, {'x1|posC', 'x2|rr'});
    verifyTrue(tc, isfield(uo, 'covOverride'));
    verifyTrue(tc, isfield(uo.covOverride, 'bySys'));
    verifyTrue(tc, isfield(uo.covOverride.bySys, 'x1'));
    verifyTrue(tc, isfield(uo.covOverride.bySys.x1, 'posC'));
    verifyTrue(tc, isfield(uo.covOverride.bySys.x2, 'rr'));
end

function test_jointModelOpt_emits_the_arc_factors_as_a_std(tc)
% categories: contract, precision
    % R7-13 (WP-116 part 2): model.arcFactor (rowWeightKey -> the arc's VARIANCE
    % factor relative to its stream's level) leaves as covOverride.byArc, a
    % STANDARD DEVIATION -- convention 1 again, now per arc. No arcFactor, no
    % byArc: a run without arc factors builds the same UserOpt as before.
    model = struct('S', [1; 0.3], 'sigma2', [], 'iter', 2);
    uo = lsqvce.engine.jointModelOpt(model, {'x1|a'});
    verifyFalse(tc, isfield(uo.covOverride, 'byArc'), 'no factors, no per-arc layer');
    k1 = lsqvce.util.rowWeightKey(1, 3, 'a', 2);
    k2 = lsqvce.util.rowWeightKey(1, 3, 'a', 1);
    model.arcFactor = containers.Map({k1, k2}, {9, 0.25});
    uo = lsqvce.engine.jointModelOpt(model, {'x1|a'});
    % tol: value plumbing identity (sqrt of a stored number); house round-trip bound 1e-12
    verifyEqual(tc, [uo.covOverride.byArc(k1), uo.covOverride.byArc(k2)], [3 0.5], 'AbsTol', 1e-12, ...
        'a variance factor of 9 arrives as a relative std of 3');
end

function test_jointModelOpt_sets_psdMode_on_useropt(tc)
% categories: contract
    % CONVENTION 3: psdMode is read from the UserOpt, NOT from the spec. It must
    % appear at the top level or covModel silently uses its 'ar' default.
    model = struct('S', [1; 0.2], 'sigma2', 1, 'iter', 1);
    uo = lsqvce.engine.jointModelOpt(model, {'x1|pos'});
    verifyEqual(tc, uo.psdMode, 'toeplitz', ...
        'joint mode must ask for the direct Toeplitz build, not an AR fit');
    uo2 = lsqvce.engine.jointModelOpt(model, {'x1|pos'}, 'PsdMode', 'ar');
    verifyEqual(tc, uo2.psdMode, 'ar');
end

function test_jointModelOpt_psd_is_a_shape(tc)
% categories: precision
    % The pooled law is normalised (r(0) -> 1): an S ten times larger gives the
    % SAME PSD; the level rides on sigmaScale (per stream), never on S.
    S = [2; 1; 0.5];
    u1 = lsqvce.engine.jointModelOpt(struct('S', S,     'sigma2', 1, 'iter', 1), {'x1|pos'});
    u2 = lsqvce.engine.jointModelOpt(struct('S', 10*S,  'sigma2', 1, 'iter', 1), {'x1|pos'});
    % tol: 1e-12 relative -- one normalisation of three numbers, evaluated twice
    verifyEqual(tc, u2.covOverride.bySys.x1.pos.S, u1.covOverride.bySys.x1.pos.S, 'RelTol', 1e-12);
    % tol: value plumbing identity (sqrt(1)); house round-trip bound 1e-12
    verifyEqual(tc, u1.covOverride.bySys.x1.pos.sigmaScale, 1, 'AbsTol', 1e-12);
    verifyTrue(tc, all(u1.covOverride.bySys.x1.pos.S >= 0));
end
function test_jointModelOpt_accepts_cell_streams(tc)
% categories: contract
    model = struct('S', [1; 0.3], 'sigma2', 1, 'iter', 1);
    uo = lsqvce.engine.jointModelOpt(model, {{1, 'pos'}});
    verifyTrue(tc, isfield(uo.covOverride.bySys, 'x1'), ...
        'a numeric sysID must become the makeValidName key');
end

function test_jointModelOpt_rejects_bad_input(tc)
% categories: contract
    verifyError(tc, @() lsqvce.engine.jointModelOpt(struct('a',1), {'x1|pos'}), ...
        'lsqvce:jointModelOpt:badModel');
    verifyError(tc, @() lsqvce.engine.jointModelOpt( ...
        struct('S', [], 'sigma2', 1, 'iter', 1), {'x1|pos'}), ...
        'lsqvce:jointModelOpt:badS');
    verifyError(tc, @() lsqvce.engine.jointModelOpt( ...
        struct('S', [1; -0.2], 'sigma2', 1, 'iter', 1), {'x1|pos'}), ...
        'lsqvce:jointModelOpt:badS');
    verifyError(tc, @() lsqvce.engine.jointModelOpt( ...
        struct('S', [1;0.2], 'sigma2', 1, 'iter', 1), {'nopipe'}), ...
        'lsqvce:jointModelOpt:badStream');
end

%% ============================================ axis stride (epoch-major lags)
function test_axis_pooling_recovers_temporal_law(tc)
% categories: precision
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
    % tol: statistical bar on this fixture, sized to the estimator scatter (HANDOVER 6.1.2); not an eps claim
    verifyEqual(tc, aAx(2), rho, 'AbsTol', 0.08, ...
        'pooling per axis down epochs must recover the AR(1) coefficient');
    verifyGreaterThan(tc, abs(aRow(2) - rho), 0.3, ...
        ['a ROW-based lag must MISS rho -- if this fails the interleave is ' ...
         'not epoch-major and the guard is not testing what it claims']);
    % the spurious row-based peak sits at lag nAx
    verifyGreaterThan(tc, aRow(1+nAx), aRow(2), ...
        'the row-based ACV must peak at the axis stride, not at lag 1');
end

function v = lagSumsLocal(M)
    % v(n+1) = sum of M over the +n and -n diagonals (n = 0..m-1), so that
    % v' * r = trace(M * toeplitz(r)) -- the dense reference of redundancyLags.
    m = size(M, 1);  v = zeros(m, 1);
    v(1) = trace(M);
    for n = 1:m-1
        v(n+1) = sum(diag(M, -n)) + sum(diag(M, n));
    end
end

function S = armaSpecLocal(ar, ma, L)
    % The ARMA spectral SHAPE |MA|^2 / |AR|^2 on the L bins f_j = j / (2(L-1)).
    f = (0:L-1).' / (2*(L-1));
    z = @(c) exp(-2i*pi*f*(0:numel(c)-1)) * c(:);
    S = abs(z(ma)).^2 ./ abs(z(ar)).^2;
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
%   loader -- only enough for LSQVCE.ENGINE.RESOLVENC to derive nc. A raw spec struct
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
    % The FIRST model per pass: the final re-accumulation carries the last
    % pass's index and must not overwrite what that pass was weighted with.
    if ~isKey(rec, model.iter), rec(model.iter) = model; end   % containers.Map is a handle
    acc = jointAcc(cfg, model);
end

function cfg = jointCfg()
    % A small, fast, deterministic problem for the joint driver. Uses an
    % explicit acf covariance so every arc exposes Cfcn (the joint estimator
    % needs the explicit W with Sigma = W'W).
    %
    % NOTE ON ISOLATION: synthLoader keeps a PERSISTENT single-day cache
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
% categories: precision
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
    sp   = lsqvce.engine.sharedProbe(U, 'numProbes', Z, 'seed', 5);
    Ll   = chol((Nll+Nll.')/2, 'lower');
    r = lsqvce.kernel.redundancyLags(AtilBar, W, sp.UiZ, [], 3, ...
            'LocalBbar', Bw, 'LocalLl', Ll, 'WantRbar', true);
    % Reference: the JOINT projector with Njoint from the WHITENED designs
    % (the accumulator whitens), Sig^-1 - Sig^-1-metric hat of [A B]:
    Sinv = Sig \ eye(N);
    Mw = [W.'\A, W.'\B];                     % whitened joint design
    Nj = Mw.'*Mw;
    RbarEx = Sinv - (W\Mw)*(Nj\(W\Mw).');
    tolR = 6/sqrt(Z);                         % Hutchinson-limit tolerance
    % tol: 3-sigma (resp. 6/sqrt(Z)) of the estimator's own MC standard error -- self-scaling statistical bar
    verifyEqual(tc, r.Rbar, RbarEx, 'AbsTol', tolR, ...
        'Rbar must converge to the JOINT projector (locals absorbed)');
    % s_0 = trace(Rbar). For COLOURED Sigma this is NOT N-(nc+nl) (that holds
    % only for white arcs -- asserted in the alpha2 test); the exact value is
    % trace(RbarEx), and the pre-fix common-only projector sits a structural
    % ~3.2 higher on this fixture -- the 0.5 tolerance discriminates.
    s0Ex = trace(RbarEx);
    % tol: Monte-Carlo probe scatter ~ 1/sqrt(Probes) on this fixture; 4/sqrt(Z)-style bar (HANDOVER 6.1.2), not containment
    verifyEqual(tc, r.s(1), s0Ex, 'AbsTol', 0.5, ...
        's_0 must equal trace of the JOINT Rbar');
    RcomOld = Sinv - (W\Aw)*((Aw.'*Aw)\(W\Aw).');   % the pre-WP-33 projector
    verifyGreaterThan(tc, trace(RcomOld) - s0Ex, 1.0, ...
        'the common-only s_0 must sit measurably above the joint one');
end

function test_joint_alpha2_unity_at_truth_with_locals(tc)
% categories: precision
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
    sp = lsqvce.engine.sharedProbe(U, 'numProbes', 800, 'seed', 7);
    W  = eye(m);                              % white noise
    om0 = 0;  s0 = 0;  om0Old = 0;  s0Old = 0;
    for k = 1:nArc
        Atil = As{k} - Bs{k}*G;
        eNew = ls{k} - As{k}*x - Bs{k}*y;
        r = lsqvce.kernel.redundancyLags(Atil, W, sp.UiZ, eNew, 1, ...
                'LocalBbar', Bs{k}, 'LocalLl', Ll);
        om0 = om0 + r.omega(1);  s0 = s0 + r.s(1);
        eOld = ls{k} - As{k}*x;               % the pre-fix residual
        rOld = lsqvce.kernel.redundancyLags(As{k}, W, sp.UiZ, eOld, 1);
        om0Old = om0Old + rOld.omega(1);  s0Old = s0Old + rOld.s(1);
    end
    alpha2    = om0    / s0;
    alpha2Old = om0Old / s0Old;
    % tolerance: probe se of s_0 plus the chi2 spread of Omega_0/dof
    dof = nArc*m - (nc + nl);
    tol = 3*sqrt(2/dof) + 3*(nc/sqrt(800))/dof + 0.02;
    % tol: 3-sigma of the estimator's own reported MC standard error (self-scaling bar)
    verifyEqual(tc, alpha2, 1, 'AbsTol', tol, ...
        'alpha2 must be ~1 at truth on white noise WITH locals');
    verifyTrue(tc, alpha2Old > 1.5, ...
        'the pre-WP-33 formation must stay measurably biased (the contrast)');
end

%% ======================= WP-37: structured Sigma^-1 sums, no dense W^-1 (F-09)
function test_joint_signature_typed_arguments(tc)
% categories: contract
    % WP-39 stage (a) (fails on the old code): adjustVCEjoint's typed
    % `arguments` block -- unknown names error, StepSchedule is a strict
    % member set, Probes must be a positive integer. Names were already
    % UpperCamelCase; matching is now EXACT.
    accFcn = @(varargin) struct();               % never reached: parse errors first
    mustErr = { ...
        @() lsqvce.engine.adjustVCEjoint(accFcn, 'probes', 30), ...        % old-style lowercase
        @() lsqvce.engine.adjustVCEjoint(accFcn, 'StepSchedule', 'Decay'), ...
        @() lsqvce.engine.adjustVCEjoint(accFcn, 'Probes', 2.5), ...
        @() lsqvce.engine.adjustVCEjoint(accFcn, 'NoSuchOption', 1) };
    for k = 1:numel(mustErr)
        ok = false;
        try, mustErr{k}(); catch, ok = true; end
        verifyTrue(tc, ok, sprintf('signature case %d must error', k));
    end
end

%% =========================== WP-39 stage (a): typed-signature enforcement
function test_stageB_signatures_psd_and_lags(tc)
% categories: contract
    % WP-39 stage (b) (F-07): psdVCE's typed prior (R7 WP-111 part 3: a
    % spectrum, so non-negative -- the ZeroFreqRule enum it replaced is gone)
    % and redundancyLags' typed positionals.
    ok = false;
    try
        lsqvce.engine.psdVCE(ones(4,2), ones(4,2), [1; -1; 1; 1]);
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeNonnegative');
    end
    verifyTrue(tc, ok, 'a negative prior spectrum must fail mustBeNonnegative');
    ok = false;
    try
        lsqvce.engine.psdVCE(ones(4,2), ones(4,2), ones(4,1), 'ZeroFreqRule', 'keep');
    catch
        ok = true;
    end
    verifyTrue(tc, ok, 'the retired ZeroFreqRule is an unknown option');
    ok = false;
    try
        lsqvce.engine.psdVCE(ones(4,2), ones(4,2), ones(4,1), 'Gamma', 1.5);
    catch e
        ok = startsWith(e.identifier, 'MATLAB:validators:');
    end
    verifyTrue(tc, ok, 'psdVCE Gamma > 1 must fail the typed constraint');
    ok = false;
    try
        lsqvce.kernel.redundancyLags(ones(6,2), eye(6), ones(2,3), [], 0);
    catch e
        ok = startsWith(e.identifier, 'MATLAB:validators:');
    end
    verifyTrue(tc, ok, 'redundancyLags nMaxLag = 0 must fail the signature');
end

function test_joint_signature_types_and_enum_values(tc)
% categories: contract
    % WP-39 (F-07), corrected premise (round a1-2): NV names match
    % case-insensitively even under arguments; the enforced breaks are
    % enum VALUE casing (StepSchedule 'Decay' passed the old strcmpi,
    % fails mustBeMember now -- the fails-on-old lever), typed
    % constraints, and unknown-name rejection.
    dummy = @(s) s;
    ok = false;
    try
        lsqvce.engine.adjustVCEjoint(dummy, 'StepSchedule', 'Decay');
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    verifyTrue(tc, ok, 'StepSchedule ''Decay'' must be rejected (value casing)');
    ok = false;
    try
        lsqvce.engine.adjustVCEjoint(dummy, 'MaxIter', 'three');
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    verifyTrue(tc, ok, 'non-numeric MaxIter must be rejected by the signature');
    ok = false;
    try
        lsqvce.engine.adjustVCEjoint(dummy, 'NoSuchOption', 1);
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    verifyTrue(tc, ok, 'an unknown option name must be rejected');
end

function test_lags_structured_sums_match_dense_toeplitz(tc)
% categories: precision
    % Agreement at 1e-12 (spec): the Levinson + Gohberg-Semencul path must
    % reproduce s_n from the dense reference X = W\eye(N) on the standing
    % TOEPLITZ fixture -- with probes AND locals engaged so every term of
    % s_n is exercised.
    rng(371);
    N = 90;  P = 3;  nl = 2;
    [~, W] = spdToeplitz(N);
    A = randn(N, P);  B = randn(N, nl);
    Aw = W.'\A;  Bw = W.'\B;
    Nll = Bw.'*Bw;  G = Nll\(Bw.'*Aw);  At = Aw - Bw*G;
    U  = chol((At.'*At + (At.'*At).')/2);
    sp = lsqvce.engine.sharedProbe(U, 'numProbes', 40, 'seed', 9);
    Ll = chol((Nll+Nll.')/2, 'lower');
    r  = lsqvce.kernel.redundancyLags(At, W, sp.UiZ, [], 6, ...
             'LocalBbar', Bw, 'LocalLl', Ll);
    % dense reference, built exactly as the pre-WP-37 kernel did:
    X  = W \ eye(N);
    Rt = W \ (At * sp.UiZ);
    K  = W \ (Bw / Ll.');
    sRef = zeros(6,1);
    sRef(1) = sum(X(:).^2) - sum(Rt(:).^2) - sum(K(:).^2);
    for n = 1:5
        idx = 1:(N-n);
        sRef(n+1) = 2*sum( sum(X(idx+n,:).*X(idx,:),2) ...
                         - sum(Rt(idx+n,:).*Rt(idx,:),2) ...
                         - sum(K(idx+n,:).*K(idx,:),2) );
    end
    % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
    verifyEqual(tc, r.s, sRef, 'RelTol', 1e-12, ...
        'GS-path s_n must match the dense reference at 1e-12');
end

function test_lags_structured_sums_match_dense_general(tc)
% categories: precision
    % The NON-Toeplitz fallback (ring-streamed columns of W^-T): same 1e-12
    % agreement on a deliberately non-stationary covariance.
    rng(372);
    N = 70;  P = 3;
    d = linspace(1, 3, N).';                     % varying variance: not Toeplitz
    [Sig0, ~] = spdToeplitz(N);
    Sig = (d*d.') .* Sig0;
    W = chol((Sig+Sig.')/2);
    A = randn(N, P);
    U = chol(A.'*(Sig\A));
    sp = lsqvce.engine.sharedProbe(U, 'numProbes', 30, 'seed', 11);
    Abar = W.'\A;
    r = lsqvce.kernel.redundancyLags(Abar, W, sp.UiZ, [], 5);
    X  = W \ eye(N);
    Rt = W \ (Abar * sp.UiZ);
    sRef = zeros(5,1);
    sRef(1) = sum(X(:).^2) - sum(Rt(:).^2);
    for n = 1:4
        idx = 1:(N-n);
        sRef(n+1) = 2*sum( sum(X(idx+n,:).*X(idx,:),2) ...
                         - sum(Rt(idx+n,:).*Rt(idx,:),2) );
    end
    % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
    verifyEqual(tc, r.s, sRef, 'RelTol', 1e-12, ...
        'fallback-path s_n must match the dense reference at 1e-12');
end

function test_lags_long_arc_speed_no_cubic_term(tc)
% categories: robustness, speed
    % Paragraph-6.1 SPEED lever: on a LONG Toeplitz arc the structured path
    % is O(N^2) -- the pre-WP-37 kernel materialised X = W\eye(N) (O(N^3)
    % work, O(N^2) memory) and blows the budget below by an order of
    % magnitude on the runner (measured there: dense ~12-30 s at N = 4000;
    % structured < 1 s; budget 5 s leaves slow-runner margin without
    % readmitting the cubic term).
    N = 4000;
    base = [1; 0.5*0.8.^(1:N-1).'];
    c = base;  c(1) = 1.5;
    Sig = toeplitz(c);                           % SPD by diagonal dominance margin
    W = chol(Sig);
    A = randn(RandStream('twister','Seed',7), N, 2);
    U = chol(A.'*A);
    sp = lsqvce.engine.sharedProbe(U, 'numProbes', 8, 'seed', 3);
    tStart = tic;
    r = lsqvce.kernel.redundancyLags(W.'\A, W, sp.UiZ, [], 4);
    tElapsed = toc(tStart);
    verifyTrue(tc, isfinite(r.s(1)) && r.s(1) > 0, 'sanity: s_0 computed');
    verifyLessThan(tc, tElapsed, 5.0, ...
        sprintf('structured path took %.2f s at N = %d (cubic would blow this)', ...
                tElapsed, N));
end

function test_redundancyLags_nl0_path_bitwise_unchanged(tc)
% categories: robustness
    % Spec regression: with NO local columns the new machinery must be a
    % NO-OP -- bitwise: the nl = 0 branch keeps the pre-WP-33 expressions
    % verbatim, and an explicit empty LocalBbar routes to that same branch.
    [Abar, W, sp, ebar] = lagFixture(tc, 20, 3, 60);
    a = lsqvce.kernel.redundancyLags(Abar, W, sp.UiZ, ebar, 5);
    b = lsqvce.kernel.redundancyLags(Abar, W, sp.UiZ, ebar, 5, ...
            'LocalBbar', zeros(20,0), 'LocalLl', []);
    verifyEqual(tc, b.s,     a.s);            % NO tolerance: bitwise
    verifyEqual(tc, b.omega, a.omega);
end

function test_adjustVCEjoint_flat_psd_on_white_noise_with_locals(tc)
% categories: contract
    % Driver-level accuracy (spec 1 at the adjustVCEjoint level): on the
    % white (diag) fixture WITH nl = 2 locals the estimated spectrum must stay
    % FLAT -- no arc-parameter signal folded into the "noise", which would
    % lift the LOW bins (pre-WP-33 the omitted B*y inflated the low lags,
    % Omega_0 x3.8, T3).
    % R7 WP-111 part 3: ONE pass. The fixture keeps the declared (white) law
    % fixed, so a second pass multiplies the same bin scatter in again; and
    % the bars are the estimator's own scatter, not a wish: 6 bins on 2 arcs
    % of 40 epochs give each bin a redundancy of only 7-15. Simulated with
    % exact traces (white truth, 6 parameters, 2000 draws): max/min of the
    % bins median 2.10, p99 5.06; mean(bins 1-2) / mean(bins 3-6) p99 2.14,
    % p99.9 2.58. (The DCT-era bar max/min < 3 after two passes sits at
    % p95 13 here -- CI ac38c75.)
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % provoked on purpose here (WP-47, gate rule 7)
    cfg  = jointCfg();
    r = lsqvce.engine.adjustVCEjoint(@(model) jointAcc(cfg, model), ...
            'NumLags', 6, 'MaxIter', 1, 'Probes', 200, 'Seed', 11, ...
            'Verbose', false);
    S = r.S;
    verifyLessThan(tc, mean(S(1:2)) / mean(S(3:end)), 2.6, ...
        sprintf(['white truth => no low-frequency excess: low/high = %.2f ' ...
                 '(arc-parameter leakage lifts the low bins)'], mean(S(1:2)) / mean(S(3:end))));
    verifyLessThan(tc, max(S) / min(S), 6, ...
        sprintf('white truth => flat within the bin scatter: max/min = %.2f', max(S) / min(S)));
end
function test_joint_stabilised_alternation_beats_collapse(tc)
% categories: robustness
    % WP-40 r2 (F-15; fails on the old code: the options are unknown): the
    % closed-loop repro fixture (jointAccMapped, the F-15 vehicle) run WITH
    % 'StepSchedule','decay' + 'AvgTail',4 must stay far above the undamped
    % collapse (which hits 0.18 on this seed). Round-6 correction: the
    % earlier "min ~0.89 across all grids" was measured on a mirror whose
    % DECLARATION path was unfaithful (zero-extended dctI acv); the
    % psd2acv-faithful mirror puts the stabilised seed floor at ~0.59 here
    % and shows the closed loop never beating open-loop pass-1 quality
    % (~0.96) at any tested volume -- full grid in the §9 WP-40 round-6
    % note.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % provoked on purpose here (WP-47, gate rule 7)
    cfg  = jointCfgColoured();
    nLag = 10;
    % R7-30: no blanket warning('off','all') -- it blinded gate rule 7 here
    r = lsqvce.engine.adjustVCEjoint(@(model) jointAccMapped(cfg, model), ...
            'NumLags', nLag, 'MaxIter', 12, 'Probes', 300, 'Seed', 13, ...
            'StepSchedule', 'decay', 'AvgTail', 4, 'Verbose', false);
    Strue = armaSpecLocal(cfg.cov.pos.ar, cfg.cov.pos.ma, nLag);   % the truth on the bins
    cc = corrcoef(r.S(2:end), Strue(2:end));
    % Round-6 bar (FAITHFUL-mirror measured; the 0.8 bar rested on the
    % unfaithful mirror's 0.87-0.91 promise): with the psd2acv declaration
    % the stabilised loop's seed floor is ~0.59 at this volume and the
    % UNDAMPED loop collapses to 0.18 on this very seed (r2f record). The
    % honest property is DAMAGE LIMITATION: the stabilised loop must stay
    % far above the collapse. Production guidance (F-15 §8): open loop or
    % <= 3 damped passes; the closed loop does not beat open-loop quality.
    verifyGreaterThan(tc, cc(1,2), 0.5, ...
        sprintf(['stabilised closed loop: corr = %.3f (damage-limitation ' ...
                 'floor 0.5; undamped collapses to ~0.18 here)'], cc(1,2)));
end

function test_joint_default_schedule_is_plain_step(tc)
% categories: contract
    % Backward compatibility: 'StepSchedule','none' (the default) must be
    % BITWISE the old path -- one pass with and without the explicit default
    % options must agree exactly.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % provoked on purpose here (WP-47, gate rule 7)
    cfg = jointCfg();
    % R7-30: no blanket warning('off','all') -- it blinded gate rule 7 here
    a = lsqvce.engine.adjustVCEjoint(@(model) jointAcc(cfg, model), ...
            'NumLags', 6, 'MaxIter', 1, 'Probes', 100, 'Seed', 5, ...
            'Verbose', false);
    b = lsqvce.engine.adjustVCEjoint(@(model) jointAcc(cfg, model), ...
            'NumLags', 6, 'MaxIter', 1, 'Probes', 100, 'Seed', 5, ...
            'StepSchedule', 'none', 'AvgTail', 0, 'Verbose', false);
    verifyEqual(tc, b.S, a.S);                   % NO tolerance: bitwise
    verifyEqual(tc, b.C,    a.C);
end

function test_adjustVCEjoint_psd_recovery_coloured_with_locals(tc)
% categories: robustness
    % WP-33 round 2, FINAL SHAPE (r2e; full measurement record in HANDOVER
    % §9 WP-33): an ACCURACY FLOOR for the whole coloured chain with local
    % parameters -- corr(Shat, Strue) > 0.85 at MaxIter 1 with the default
    % (non-delta) C0 and fixed-white declared weighting. The Python mirror
    % of the exact pipeline measured, across seeds and sizes: pass-1 corr
    % 0.92-0.96 (min 0.921); more arcs, rows or lags do NOT lift it above
    % ~0.95 (the estimator's truncated-lag DCT ceiling), the no-rewhitening
    % fixed point is corr ~ 0.64 (why r2b/r2c/r2d all landed ~0.6), and a
    % near-delta C0 flattens alpha2 in one pass (the documented degenerate-
    % start throttle). The spec's corr > 0.95 bar is therefore MEASURED
    % unattainable here, and corr does not discriminate F-03 anyway (old-
    % code corr ranged 0.55-0.95 over seeds) -- the §6.1 fail-on-old-code
    % lever for F-03 is round 1's alpha2-at-truth test, NOT this one. This
    % test guards the chain's recovery quality against future regressions.
    % R7 WP-111 part 3: the law is the bin spectrum S (white start) and the
    % truth is the exact ARMA spectrum on the same bins -- the DCT of the
    % truncated ACF it was compared with before is gone with the DCT.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % provoked on purpose here (WP-47, gate rule 7)
    cfg  = jointCfgColoured();
    nLag = 10;
    % R7-30: no blanket warning('off','all') -- it blinded gate rule 7 here
    r = lsqvce.engine.adjustVCEjoint(@(model) jointAcc(cfg, model), ...
            'NumLags', nLag, 'MaxIter', 1, 'Probes', 300, 'Seed', 13, ...
            'Verbose', false);
    Strue = armaSpecLocal(cfg.cov.pos.ar, cfg.cov.pos.ma, nLag);   % the truth on the bins
    cc = corrcoef(r.S(:), Strue);
    verifyGreaterThan(tc, cc(1,2), 0.85, ...
        sprintf('coloured-chain accuracy floor: corr = %.3f (need > 0.85)', ...
                cc(1,2)));
end

%% ========================== WP-34: keyed arc identity in adjustVCEjoint (F-02)
function test_joint_keyed_skip_keeps_factors_aligned(tc)
% categories: precision, robustness
    muteWarn(tc, 'lsqvce:decorrFactor:notPD');   % R7-30: the PD-toggle fixture makes an arc indefinite on purpose
    % Spec test (a)+(b): induce a mid-run skip (one arc's covariance goes
    % non-PD on pass 2) and verify the factors stay ALIGNED BY KEY. Exact
    % assertions: the key sets are pinned and identical across runs; the
    % skip is reported as {2, key, 'nonPD'}; the SKIPPED arc's factor is
    % byte-identical to its own pass-1 value (kept, per the WP contract).
    % The OTHER arcs are compared with tolerance, NOT byte-identically:
    % arcVarianceFactors normalises by gammaBar (eqs. 6.5.29-30), a mean
    % over the PRESENT arcs, so removing one column renormalises every
    % factor by construction -- byte identity across a changed arc set
    % would require decoupling the estimator (deviation pinned in the
    % HANDOVER round note). On the pre-WP-34 code these assertions ERROR
    % (no out.arcKeys/out.skips) and the alignment itself is the F-02
    % shift -- the §6.1 lever.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % provoked on purpose here (WP-47, gate rule 7)
    cfg = jointCfg();
    badArc = 2;                          % skipped on pass 2 in run B
    % R7-30: no blanket warning('off','all') -- it blinded gate rule 7 here
    runA1 = lsqvce.engine.adjustVCEjoint(@(m) jointAccPDToggle(cfg, m, [], badArc), ...
                'NumLags', 6, 'MaxIter', 1, 'Probes', 150, 'Seed', 5, ...
                'Verbose', false);
    runA  = lsqvce.engine.adjustVCEjoint(@(m) jointAccPDToggle(cfg, m, [], badArc), ...
                'NumLags', 6, 'MaxIter', 2, 'Probes', 150, 'Seed', 5, ...
                'Verbose', false);
    % 1 of 2 arcs skipped is a 50 % skip RATE: above the default
    % SkipThreshold the run REFUSES (R5 M1, tested below); this test is
    % about the alignment of the survivors, so the threshold is raised
    % deliberately -- the documented way.
    runB  = lsqvce.engine.adjustVCEjoint(@(m) jointAccPDToggle(cfg, m, 2, badArc), ...
                'NumLags', 6, 'MaxIter', 2, 'Probes', 150, 'Seed', 5, ...
                'SkipThreshold', 1, 'Verbose', false);

    verifyEqual(tc, runB.arcKeys, runA.arcKeys, ...
        'the pinned arc identity must not depend on skips');
    kSkip = find(contains(runA.arcKeys, sprintf('_a%d', badArc)) & ...
                 contains(runA.arcKeys, '_d1_'));
    verifyNumElements(tc, kSkip, 1);
    verifyEqual(tc, size(runB.skips, 1), 1, 'exactly one skip reported');
    verifyEqual(tc, runB.skips{1,1}, 2,        'skip pass');
    verifyEqual(tc, runB.skips{1,2}, runA.arcKeys{kSkip}, 'skip key');
    verifyEqual(tc, runB.skips{1,3}, 'nonPD',  'skip reason');
    verifyEmpty(tc, runA.skips, 'the control run must not skip');
    % the skipped arc KEEPS its pass-1 factor -- byte-identical:
    verifyEqual(tc, runB.sigma2(kSkip), runA1.sigma2(kSkip));
    % the others: same keys, values within the gammaBar renormalisation
    others = setdiff(1:numel(runA.arcKeys), kSkip);
    % tol: gammaBar renormalisation shift when one arc is skipped; 0.15 is
    % the WP-34 approved renormalisation bound (ledger row 34)
    verifyEqual(tc, runB.sigma2(others), runA.sigma2(others), ...
        'RelTol', 0.15, 'non-skipped factors track the no-skip run');
end


function test_joint_skip_rate_above_threshold_refuses(tc)
% categories: contract, robustness
    muteWarn(tc, 'lsqvce:decorrFactor:notPD');   % R7-30: the PD-toggle fixture makes an arc indefinite on purpose
    % R5 M1 (WP-80): a skip RATE above SkipThreshold is a REFUSAL, not a
    % smaller estimate. Same fixture as the alignment test: pass 2 makes one
    % of two arcs non-PD (50 %), which exceeds the default 20 %.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');
    cfg = jointCfg();
    badArc = 2;
    % R7-30: no blanket warning('off','all') -- it blinded gate rule 7 here
    verifyError(tc, @() lsqvce.engine.adjustVCEjoint(@(m) jointAccPDToggle(cfg, m, 2, badArc), ...
                'NumLags', 6, 'MaxIter', 2, 'Probes', 150, 'Seed', 5, 'Verbose', false), ...
                'lsqvce:adjustVCEjoint:skipThreshold');
    % at exactly the rate (0.5) the run is NOT refused: strictly greater refuses
    runB = lsqvce.engine.adjustVCEjoint(@(m) jointAccPDToggle(cfg, m, 2, badArc), ...
                'NumLags', 6, 'MaxIter', 2, 'Probes', 150, 'Seed', 5, ...
                'SkipThreshold', 0.5, 'Verbose', false);
    verifyEqual(tc, size(runB.skips, 1), 1);
end

function test_joint_arc_mismatch_raises(tc)
% categories: contract, reference, robustness
    % Spec test (c): an arc that VANISHES FROM THE ENUMERATION on pass 2 is
    % an identity error, not a skip -- arcMismatch, mirroring vce's
    % groupMismatch contract ("the loader must be deterministic").
    cfg = jointCfg();
    % R7-30: no blanket warning('off','all') -- it blinded gate rule 7 here
    verifyError(tc, @() lsqvce.engine.adjustVCEjoint( ...
            @(m) jointAccDropArc(cfg, m, 2), ...
            'NumLags', 6, 'MaxIter', 2, 'Probes', 150, 'Seed', 5, ...
            'Verbose', false), ...
        'lsqvce:adjustVCEjoint:arcMismatch');
end

function acc = jointAccPDToggle(cfg, model, badIter, badArc)
    % jointAcc with a per-arc covariance override: identity normally, but
    % NON-PD for arc badArc while model.iter == badIter (the WP-34 mid-run
    % skip trigger; the covariance model updating between passes is exactly
    % the realistic F-02 scenario). [] badIter = never (the control run).
    it = 0;
    if ~isempty(model) && isfield(model, 'iter'), it = model.iter; end
    bad = ~isempty(badIter) && it == badIter;
    fcn = @(tb, o) pdToggleCov(numel(tb), o, bad, badArc);
    uo  = struct('covOverride', struct('pos', ...
                 struct('type','custom', 'fcn', fcn, 'sigma', 1)));
    cfgO = cfg;  cfgO.override = uo;
    acc = jointAcc(cfgO, model);
end

function Sig = pdToggleCov(m, o, bad, badArc)
    % INDEFINITE, not negative definite: one -0.5 eigenvalue with max eig
    % 1 > 0 means the ACCUMULATOR survives (decorrFactor floors negative
    % directions with a warning -- the robust production path) while the
    % sweep's strict chol fails and oneArc SKIPS -- exactly the realistic
    % mid-run F-02 trigger (round-1 lesson: -eye kills decorrFactor during
    % accumulation before any sweep runs, notPD hard error).
    Sig = eye(m);
    if bad && isfield(o, 'arc') && isequal(o.arc, badArc)
        Sig(1, 1) = -0.5;                % indefinite: chol fails, eig-max > 0
    end
end

function acc = jointAccDropArc(cfg, model, dropIter)
    % jointAcc whose META drops the last arc from pass dropIter on -- the
    % enumeration itself changes, which must raise arcMismatch. Mirrors
    % jointAcc (addSystem only; the driver streams the data itself).
    it = 0;
    if ~isempty(model) && isfield(model, 'iter'), it = model.iter; end
    drop = it >= dropIter;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'KeepLocal', true, 'Verbose', false);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, ...
            @(d, ot, a, uo2) dropArcLoader(cfg, s, d, ot, a, uo2, drop), 1);
    end
end

function out = dropArcLoader(cfg, s, d, ot, a, uo2, drop)
    % synthLoader pass-through; in META mode (first loader argument is the
    % literal 'meta') the last arc is removed from the enumeration when
    % `drop` is armed -- an identity change, not a skip.
    out = lsqvce.engine.synthLoader(cfg, s, d, ot, a, uo2);
    if drop && (ischar(d) || isstring(d)) && strcmp(char(d), 'meta') ...
            && isstruct(out) && isfield(out, 'arcs')
        out.arcs{1} = out.arcs{1}(1:end-1);
    end
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
    cfg.arcsPerDay = 16;               % ONE day: 16 x 256 rows (the mirror's
    cfg.mArc       = 256;              %  stable regime, F-15 evidence)
    cfg.dt         = 5;
    cfg.seed       = 9223;             % distinct seed -> distinct cache key
    cfg.collinear  = false;
    cfg.collinearEps = 1e-7;
    cfg.cov.pos = struct('type','arma', 'ar',[1 -0.75 0.4], 'ma',1, 'sigma',1.0);
    % r2f: the floor test's premise is FIXED-WHITE DECLARED weighting while
    % GENERATING AR(2). r2e switched the test call to jointAcc but dropped
    % the override, so the loader declared cfg.cov — THE TRUTH — and the
    % r2a trap returned (alpha2 ~ 1, Shat ~ spectrum of the default C0,
    % CI corr 0.182). White must be declared explicitly; jointAcc threads
    % cfg.override into the loader's UserOpt:
    cfg.override = struct('covOverride', struct('pos', ...
                          struct('type','diag','sigma',1.0)));
end

function acc = jointAccMapped(cfg, model) %#ok<DEFNU>
    % NOT called by the r2f floor test: kept as the documented CLOSED-LOOP
    % reproduction fixture for finding F-15 (alternation bistability, §8;
    % queued on WP-40) and the wp33_loop_mirror.py cross-reference.
    % The coloured fixture's CLOSED LOOP: the driver's current law is
    % declared back as the covariance and psdMode 'toeplitz' gives the dense
    % Cfcn oneArc needs. Generation always follows cfg.cov.
    % R7 WP-111 part 3: the driver's law is the bin spectrum model.S; it is
    % declared back exactly as jointModelOpt emits it ('psd' on fNorm, the
    % ordinates 2 S scaled to r(0) = 1), so the loader's covariance IS K S.
    uo = struct();
    if ~isempty(model) && isfield(model, 'S') && ~isempty(model.S)
        J = lsqvce.engine.jointModelOpt(struct('S', model.S), {'x92|pos'});
        sp = J.covOverride.bySys.x92.pos;
        uo.covOverride.pos = struct('type', 'psd', 'fNorm', sp.fNorm, 'S', sp.S, 'sigma', 1);
        uo.psdMode = 'toeplitz';
    else
        uo.covOverride.pos = struct('type','diag', 'sigma', 1.0);
    end
    uo = withArcFactors(uo, model);              % R7-13: the contract since 6.1.6
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'KeepLocal', true, 'Verbose', false, ...
                                  'UserOpt', uo);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, ...
            @(d,ot,a,uo2) lsqvce.engine.synthLoader(cfg, s, d, ot, a, uo2), 1);
    end
end

function acc = jointAcc(cfg, model)
    % Rebuild the accumulator. A production caller would map model.S / model.sigma2
    % into the loader's covariance; here the fixture keeps the LAW fixed so
    % the test exercises the DRIVER's plumbing, not the loader's. A cfg.override
    % reaches the loader as UserOpt (generic hook; the coloured fixture now
    % owns its declaration in jointAccMapped instead). The per-arc factors ARE
    % applied (withArcFactors): since 6.1.6 the driver's update assumes them,
    % and a fixture that drops them compounds them (R7-13; CI 2cc7efb).
    uo = struct();
    if isfield(cfg, 'override'), uo = cfg.override; end
    uo = withArcFactors(uo, model);
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'KeepLocal', true, 'Verbose', false, ...
                                  'UserOpt', uo);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, ...
            @(d,ot,a,uo) lsqvce.engine.synthLoader(cfg, s, d, ot, a, uo), 1);
    end
end

function uo = withArcFactors(uo, model)
    % R7-13 (WP-116 part 2): the per-arc factors of the driver's model, as the
    % ONE translation emits them (jointModelOpt -> covOverride.byArc), added to
    % a fixture's own UserOpt. No stream list: the fixture keeps its own law and
    % level; only the relative per-arc layer is taken. CI 2cc7efb measured why
    % it is needed: without it the factors compounded over the passes and
    % arcFactorsDiverging fired -- the guard doing its job on a fixture that
    % broke the new accFcn contract.
    if isempty(model) || ~isstruct(model) || ~isfield(model, 'arcFactor') || isempty(model.arcFactor) ...
            || model.arcFactor.Count == 0
        return
    end
    J = lsqvce.engine.jointModelOpt(model, {});
    uo.covOverride.byArc = J.covOverride.byArc;
end

function [Abar, W, sp, ebar, N] = lagFixture(tc, N, P, Z) %#ok<INUSL>
    % A small decorrelated arc plus its shared probe block.
    rng(41);
    [~, W] = spdToeplitz(N);
    A      = randn(N, P);
    Abar   = W.' \ A;                                 % W^-T A
    U      = chol(Abar.'*Abar);
    sp     = lsqvce.engine.sharedProbe(U, 'numProbes', Z, 'seed', 9);
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
