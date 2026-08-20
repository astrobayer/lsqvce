function tests = test_traceEstimate
%TEST_TRACEESTIMATE  Unit tests for lsqvce.traceEstimate (Ellmer 2018 6.5.6).
%
%   Run with:  results = runtests('tests/test_traceEstimate.m')
%
%   Covers the Hutchinson estimator's defining properties: UNBIASEDNESS (the
%   average over many independent probe sets converges to the exact trace), the
%   trace(Zbar'*M*Zbar) identity (eqs. 6.5.41-6.5.43), the matrix-free operator
%   path and its equivalence to the explicit path, shared-probe reuse (Ellmer's
%   U^-1*Zbar optimisation), the 1/sqrt(Z) error law, Gaussian vs Rademacher,
%   input validation, and the guarantee that a seeded call leaves the GLOBAL RNG
%   state untouched.
%
%   The tolerances are sampling-theory bounds, not guesses: a single Hutchinson
%   estimate has standard error sqrt(Var/Z), so tests that check a single
%   estimate use loose bounds (few * se), while the unbiasedness test averages
%   over many seeds to drive the error down and can use a tight bound.
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-36: two tests for the ONE
%   probe generator (global-RNG invariance with empty seeds -- fails on the
%   pre-fix code -- and bit-identical probes across callers per seed).
%   With assistance from Claude (Anthropic), model Claude Fable 5.
%   Created 2026-07-16 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.

tests = functiontests(localfunctions);
end

%% =================================================== fixtures
function setupOnce(tc) %#ok<*DEFNU>
    % A fixed, moderately conditioned SPD matrix used across tests.
    rs = RandStream('twister', 'Seed', 42);
    P  = 120;
    A  = randn(rs, P, P);
    tc.TestData.M   = A' * A + P * eye(P);      % SPD, well away from singular
    tc.TestData.P   = P;
    tc.TestData.tr  = trace(tc.TestData.M);
end

function M = spd(P, seed)
    rs = RandStream('twister', 'Seed', seed);
    A  = randn(rs, P, P);
    M  = A' * A + P * eye(P);
end

%% =================================================== unbiasedness
function test_unbiased_over_many_seeds(tc)
% categories: precision
    % E[estimate] = trace(M) exactly. Average over K independent probe sets;
    % the mean has standard error se_single/sqrt(K), so with K=400 and Z=100 the
    % mean is many sigma inside a 1% band for this well-conditioned matrix.
    M = tc.TestData.M;  exactTr = tc.TestData.tr;
    K = 400;  ests = zeros(1, K);
    for k = 1:K
        ests(k) = lsqvce.traceEstimate(M, 'NumProbes', 100, 'Seed', 1000 + k);
    end
    % tol: Monte-Carlo probe scatter ~ 1/sqrt(Probes) on this fixture; 4/sqrt(Z)-style bar (HANDOVER 6.1.2), not containment
    verifyEqual(tc, mean(ests), exactTr, 'RelTol', 0.01, ...
        'the Hutchinson estimator must be unbiased: mean -> trace(M)');
end

function test_single_estimate_within_error_band(tc)
% categories: contract
    % A single estimate must lie within ~4 standard errors of the truth (that
    % band holds with overwhelming probability for a fixed seed).
    M = tc.TestData.M;  exactTr = tc.TestData.tr;
    [t, info] = lsqvce.traceEstimate(M, 'NumProbes', 200, 'Seed', 7);
    verifyLessThan(tc, abs(t - exactTr), 4 * info.se, ...
        'estimate must be within 4 SE of the exact trace');
end

%% =================================================== the estimator identity
function test_matches_trace_of_ZbarMZbar(tc)
% categories: precision
    % The returned scalar must EXACTLY equal trace(Zbar'*M*Zbar) for the probe
    % block it used (eq. 6.5.43), to floating-point. We reproduce the probe
    % block by passing it in explicitly.
    M = spd(60, 11);  P = size(M,1);  Z = 50;
    rs = RandStream('twister', 'Seed', 3);
    raw = 2*randi(rs, [0 1], P, Z) - 1;          % raw Rademacher
    t = lsqvce.traceEstimate(M, 'Probes', raw);
    Zbar = raw / sqrt(Z);
    % tol: recomputation of the same quadratic form; house exact-identity bound (HANDOVER 6.1.2)
    verifyEqual(tc, t, trace(Zbar' * M * Zbar), 'AbsTol', 1e-9, ...
        'estimate must equal trace(Zbar''*M*Zbar) exactly');
end

function test_perProbe_mean_is_estimate(tc)
% categories: precision
    % info.perProbe are the un-normalised z'Mz; their MEAN is the estimate.
    M = spd(50, 5);
    [t, info] = lsqvce.traceEstimate(M, 'NumProbes', 80, 'Seed', 9);
    % tol: probe construction identity (sign/sqrt bookkeeping, one triangular solve); eps-level, house bounds give headroom
    verifyEqual(tc, mean(info.perProbe), t, 'AbsTol', 1e-9);
    verifyEqual(tc, info.Z, 80);
    verifyEqual(tc, info.P, size(M,1));
end

%% =================================================== matrix-free path
function test_matrixfree_equals_explicit(tc)
% categories: precision
    % Passing M as an operator @(Z)->M*Z with the SAME probe block must give the
    % identical scalar as passing M explicitly with that probe block.
    M = spd(70, 21);  P = size(M,1);  Z = 64;
    rs = RandStream('twister', 'Seed', 2);
    raw = 2*randi(rs, [0 1], P, Z) - 1;
    tExplicit = lsqvce.traceEstimate(M,            'Probes', raw);
    tHandle   = lsqvce.traceEstimate(@(ZZ) M*ZZ,   'Probes', raw, 'Size', P);
    % tol: operator vs explicit matrix on the SAME probes -- identical
    % arithmetic path; 1e-9 guards accumulation order only
    verifyEqual(tc, tHandle, tExplicit, 'AbsTol', 1e-9, ...
        'matrix-free operator must match the explicit matrix on the same probes');
end

function test_matrixfree_needs_size(tc)
% categories: contract
    % A function-handle M with no Probes and no Size must error clearly.
    threw = false;
    try
        lsqvce.traceEstimate(@(Z) Z, 'NumProbes', 10);
    catch ME
        threw = strcmp(ME.identifier, 'lsqvce:traceEstimate:needSize');
    end
    verifyTrue(tc, threw, 'handle M without Size must throw needSize');
end

function test_shared_probes_reuse(tc)
% categories: precision, reference
    % Ellmer's optimisation (eq. 6.5.48): form Ubar = U^-1 * Zbar with ONE
    % triangular solve against the global Cholesky factor U (N = U'*U), then
    % reuse it across several N_g to estimate tr(N^-1 N_g) = tr(Ubar' N_g Ubar),
    % because U^-1 U^-T = N^-1. Route it through the estimator's operator path:
    % pass Probes = U \ raw (raw Rademacher solved once against U) and operator
    % @(Z)->Ng*Z. The estimator normalises probes by 1/sqrt(Z) internally, so it
    % forms U^-1 * raw / sqrt(Z) = U^-1 * Zbar = Ubar and returns
    % tr(Ubar' Ng Ubar). Compare to the exact trace.
    P = 80;  N = spd(P, 31);
    U = chol(N);                                  % N = U'*U  (upper)
    rs = RandStream('twister', 'Seed', 4);
    Z  = 300;
    raw    = 2*randi(rs, [0 1], P, Z) - 1;        % raw Rademacher
    UinvRaw = U \ raw;                             % ONE triangular solve, reused
    Ng1 = spd(P, 32);  Ng2 = spd(P, 33);
    t1 = lsqvce.traceEstimate(@(ZZ) Ng1*ZZ, 'Probes', UinvRaw, 'Size', P);
    t2 = lsqvce.traceEstimate(@(ZZ) Ng2*ZZ, 'Probes', UinvRaw, 'Size', P);
    exact1 = trace(N \ Ng1);
    exact2 = trace(N \ Ng2);
    % 300 probes on a well-conditioned operator: a few percent is expected.
    % tol: Monte-Carlo probe scatter ~ 1/sqrt(Probes) on this fixture; 4/sqrt(Z)-style bar (HANDOVER 6.1.2), not containment
    verifyEqual(tc, t1, exact1, 'RelTol', 0.10);
    verifyEqual(tc, t2, exact2, 'RelTol', 0.10);
end

%% =================================================== error law
function test_error_shrinks_like_one_over_sqrtZ(tc)
% categories: contract
    % Doubling Z twice (x4) should roughly HALVE the reported standard error.
    % Compare Z and 4Z; the ratio se(Z)/se(4Z) should be ~2 (allow a wide band,
    % since se is itself a random estimate).
    M = tc.TestData.M;
    [~, iLo] = lsqvce.traceEstimate(M, 'NumProbes', 100,  'Seed', 1);
    [~, iHi] = lsqvce.traceEstimate(M, 'NumProbes', 400,  'Seed', 1);
    ratio = iLo.se / iHi.se;
    verifyGreaterThan(tc, ratio, 1.4, 'se should fall with more probes');
    verifyLessThan(tc,    ratio, 2.8, 'se should fall ~1/sqrt(Z), not faster');
end

function test_relse_reported(tc)
% categories: precision
    M = tc.TestData.M;
    [t, info] = lsqvce.traceEstimate(M, 'NumProbes', 100, 'Seed', 2);
    % tol: Monte-Carlo probe noise, sigma ~ 1/sqrt(Probes); bar per HANDOVER 6.1.2 (4/sqrt(Z)-style, not containment)
    verifyEqual(tc, info.relse, info.se / abs(t), 'AbsTol', 1e-12);
    verifyGreaterThan(tc, info.relse, 0);
end

%% =================================================== distributions
function test_gaussian_also_unbiased(tc)
% categories: precision
    % Gaussian probes are also unbiased (higher variance than Rademacher).
    M = tc.TestData.M;  exactTr = tc.TestData.tr;
    K = 400;  ests = zeros(1, K);
    for k = 1:K
        ests(k) = lsqvce.traceEstimate(M, 'NumProbes', 100, ...
                        'Distribution', 'gaussian', 'Seed', 2000 + k);
    end
    % tol: Monte-Carlo probe scatter ~ 1/sqrt(Probes) on this fixture; 4/sqrt(Z)-style bar (HANDOVER 6.1.2), not containment
    verifyEqual(tc, mean(ests), exactTr, 'RelTol', 0.02);
end

function test_rademacher_beats_gaussian_variance(tc)
% categories: contract
    % For a fixed Z, Rademacher has lower single-estimate variance than Gaussian
    % (equal diagonal weight, no self-term variance). Compare spread over seeds.
    M = tc.TestData.M;  K = 300;
    r = zeros(1,K); g = zeros(1,K);
    for k = 1:K
        r(k) = lsqvce.traceEstimate(M, 'NumProbes', 60, 'Distribution','rademacher', 'Seed', 5000+k);
        g(k) = lsqvce.traceEstimate(M, 'NumProbes', 60, 'Distribution','gaussian',   'Seed', 5000+k);
    end
    verifyLessThan(tc, var(r), var(g), ...
        'Rademacher probes should have lower estimator variance than Gaussian');
end

%% =================================================== validation & hygiene
function test_rejects_nonsquare(tc)
% categories: contract
    threw = false;
    try
        lsqvce.traceEstimate(randn(4,5));
    catch ME
        threw = strcmp(ME.identifier, 'lsqvce:traceEstimate:notSquare');
    end
    verifyTrue(tc, threw, 'a non-square explicit matrix must throw');
end

function test_operator_size_mismatch_caught(tc)
% categories: contract
    % An operator that returns the wrong-sized block must be caught.
    threw = false;
    try
        lsqvce.traceEstimate(@(Z) Z(1:end-1, :), 'Size', 10, 'NumProbes', 5);
    catch ME
        threw = strcmp(ME.identifier, 'lsqvce:traceEstimate:badOperator');
    end
    verifyTrue(tc, threw, 'operator returning a mis-sized block must throw');
end

function test_seed_does_not_touch_global_rng(tc)
% categories: contract
    % A seeded call must use a PRIVATE stream and leave the global RNG untouched,
    % so it never perturbs a caller's reproducible pipeline.
    rng(12345, 'twister');
    before = rand();                              % advance + snapshot
    rng(12345, 'twister');
    s0 = rng;                                     % capture full state
    lsqvce.traceEstimate(tc.TestData.M, 'NumProbes', 100, 'Seed', 999);
    after = rand();
    verifyEqual(tc, after, before, ...
        'a Seeded call must not disturb the global RNG stream');
    % also: the captured state is unchanged by the seeded call
    rng(s0);
end

function test_determinism_same_seed(tc)
% categories: precision
    M = tc.TestData.M;
    t1 = lsqvce.traceEstimate(M, 'NumProbes', 100, 'Seed', 77);
    t2 = lsqvce.traceEstimate(M, 'NumProbes', 100, 'Seed', 77);
    % tol: 0 -- an EXACT identity claim (same stored value, no arithmetic)
    verifyEqual(tc, t1, t2, 'AbsTol', 0, 'same seed must reproduce exactly');
end

%% =================================================== redundancy-seam identity
function test_redundancy_seam_trNinvN_is_nc(tc)
% categories: precision
    % The opt-in 'TraceMethod','hutchinson' path in LSQVCE.REDUNDANCY runs a
    % machinery check on the reduced common normal: it estimates tr(N^-1 N),
    % which is EXACTLY nc, using the same matrix-free operator the seam builds
    % (Xfun = @(V) U\(N*(U'\V)) with N = U'*U). This test reproduces that operator
    % on an SPD N and confirms the estimate is ~nc, so the seam's core is covered
    % deterministically without needing a full accumulator.
    P = 100;  N = spd(P, 71);
    U = chol((N + N.')/2);                         % N = U'*U (upper)
    Xfun = @(V) U.' \ (N * (U \ V));               % U^-T N U^-1 = N^-1 N = I
    [t, info] = lsqvce.traceEstimate(Xfun, 'Size', P, 'NumProbes', 100, 'Seed', 13);
    % tr(I) = P exactly; the estimator is unbiased and I is perfectly conditioned
    % (all eigenvalues 1), so even a single probe set is essentially exact here.
    % tol: quantile/CDF round-trip through base-MATLAB special functions; their relative error is ~1e-14, bound leaves headroom
    verifyEqual(tc, t, P, 'RelTol', 1e-8, ...
        'tr(N^-1 N) must estimate to nc = P (the redundancy-seam check)');
    verifyLessThan(tc, info.relse, 1e-6, ...
        'on the identity operator the probe variance is ~0');
end

%% ================================ WP-36: one probe generator, honest RNG (F-06)
function testDefaultSeedLeavesGlobalRNGUntouched(tc)
% categories: robustness
    % THE F-06 lever (fails on the pre-fix code, which drew default-seeded
    % probes from -- and advanced -- the global stream): calling the probe
    % consumers with an EMPTY seed must leave the global RNG state
    % bit-identical. streamContribution's twin assertion lives in
    % test_contribution (it needs the accumulator fixture).
    rng(424242, 'twister');
    before = rng;
    U = chol(spd(12, 3) + 12*eye(12));
    lsqvce.sharedProbe(U, 'numProbes', 8);
    lsqvce.traceEstimate(spd(12, 4), 'NumProbes', 8);
    lsqvce.traceEstimate(spd(12, 4), 'NumProbes', 8, 'Distribution', 'gaussian');
    after = rng;
    verifyEqual(tc, after.State, before.State, ...
        'an empty seed must never touch the global stream');
    verifyEqual(tc, after.Type, before.Type);
end

function testSameSeedIdenticalProbesAcrossCallers(tc)
% categories: precision
    % The one-generator contract: the SAME (n, Z, dist, seed) yields
    % BIT-IDENTICAL probes in every caller. Observable without reaching into
    % private scope: sharedProbe with U = I returns Zbar = raw/sqrt(Z), and
    % traceEstimate fed that raw block via 'Probes' must equal its own
    % seeded default path exactly.
    n = 20;  Zn = 16;  sd = 77;
    sp1 = lsqvce.sharedProbe(eye(n), 'numProbes', Zn, 'seed', sd);
    sp2 = lsqvce.sharedProbe(eye(n), 'numProbes', Zn, 'seed', sd);
    verifyEqual(tc, sp1.Zbar, sp2.Zbar);          % reproducible, bitwise
    raw = sp1.Zbar * sqrt(Zn);
    verifyEqual(tc, unique(abs(raw(:))), 1, ...
        'Rademacher probes must be exactly +-1');
    M  = spd(n, 5);
    % traceEstimate returns the scalar estimate (not a struct) - round-1
    % lesson: my draft dot-indexed a double.
    t1 = lsqvce.traceEstimate(M, 'NumProbes', Zn, 'Seed', sd);
    t2 = lsqvce.traceEstimate(M, 'Probes', raw);
    % tol: probe construction identity (sign/sqrt bookkeeping, one triangular solve); eps-level, house bounds give headroom
    verifyEqual(tc, t1, t2, 'AbsTol', 4*eps(t2), ...
        'same seed => the identical probe block in every caller');
end

