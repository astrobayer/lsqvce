function out = redundancyLags(Abar, W, UiZ, ebar, nMaxLag, opt)
%LSQVCE.KERNEL.REDUNDANCYLAGS  Per-lag redundancies s_n and residual sums Omega_n.
%
%   out = LSQVCE.KERNEL.REDUNDANCYLAGS(Abar, W, UiZ, ebar, nMaxLag)
%   out = LSQVCE.KERNEL.REDUNDANCYLAGS(..., 'Name', value, ...)
%
%   The per-arc kernel of Ellmer's (2018) spectral variance-component estimation.
%   For one short arc it returns, for every time lag n = 0..nMaxLag-1,
%
%       s_n     = trace( Rbar * V_n )        the redundancy of lag component n
%       Omega_n = ehat' Sigma^-T V_n Sigma^-1 ehat   the matching residual sum
%
%   which are the numerator and denominator of the per-frequency variance factor
%   once transformed to the spectral domain (eqs. 6.5.22-6.5.26).
%
%   THE ALGEBRA (Ellmer eqs. 6.5.45-6.5.56), all verified against exact forms
%     Rtilde = W^-1 * Abar * (U^-1 Zbar)                          (6.5.48)
%     Rbar   = Sigma^-1 - Rtilde*Rtilde'                          (6.5.47)
%     s_0    = trace(Rbar);   s_n = 2 * sum_i Rbar(i+n, i)        (6.5.51-6.5.52)
%     etilde = Sigma^-1 ehat = W^-1 * ebar                        (6.5.54)
%     Omega_0 = etilde'etilde;  Omega_n = 2 * sum_i etilde(i)etilde(i-n)
%                                                                 (6.5.55-6.5.56)
%   Two structural facts make this cheap. (1) V_n is sparse with ones only on the
%   +n/-n diagonals, so the traces and quadratic forms UNROLL into plain sums --
%   no sparse product is ever formed. (2) Rbar is symmetric and V_n = V_n^- +
%   (V_n^-)', so only the lower half is summed and the result doubled.
%
%   MEMORY: the dense N x N Rbar is NEVER formed on the main path. Only the n-th
%   sub-diagonals are needed, and each is a row dot-product,
%       Rbar(i+n,i) = row_{i+n}(X).row_i(X) - row_{i+n}(Rt).row_i(Rt),  X = W^-1,
%   so the work is O(N*Z) against the N x Z Rtilde plus the triangular X the
%   whitener already implies. (FUTURE OPT: for a Toeplitz Sigma the X term has an
%   O(N) structured form; exploited since WP-37: Toeplitz arcs go through
%   Levinson + Gohberg-Semencul, O(N^2) time / O(N) memory, and no dense
%   W^-1 exists outside 'WantRbar'.) Pass 'WantRbar' only for testing.
%
%   WHY UiZ IS AN INPUT (the key optimisation, Ellmer p. 69)
%     U comes from the ACCUMULATED normal equations and is therefore IDENTICAL
%     for every short arc. So U^-1 Zbar is computed ONCE by the caller and passed
%     in, and only the small per-arc triangular solve remains. Passing it in (not
%     recomputing) is the whole point: it turns an O(nArc * P^3) cost into
%     O(P^3) + O(nArc * N*P*Z).
%
%   INPUT
%     Abar     N x P double   the arc's DECORRELATED design matrix, Abar = W^-T A
%                             (whitened exactly as the accumulator whitened it)
%     W        N x N double   the arc's whitening factor with Sigma = W'*W (upper
%                             triangular, as from chol(Sigma)). (WP-41: the old
%                             "struct with fields .W / 'Solver'" sentence
%                             described an option that never existed -- removed.)
%     UiZ      P x Z double   the SHARED U^-1*Zbar block, computed once by the
%                             caller from the factored normal matrix N = U'*U and
%                             a Rademacher probe block scaled by 1/sqrt(Z)
%     ebar     N x 1 double   the arc's decorrelated residual, ebar = W^-T ehat.
%                             Pass [] to skip Omega (returns NaN) when only the
%                             redundancies are wanted.
%     nMaxLag  1x1 int        number of lags to return (1..N). Lag index runs
%                             n = 0..nMaxLag-1.
%
%   OPTIONS (Name/Value)
%   (WP-39 stage (b): canonical UpperCamelCase names below; legacy
%   casings still match -- MATLAB NV semantics -- and enum VALUES are
%   case-sensitive lowercase.)
%     'WantRbar' logical   also return the dense Rbar (N x N). Default false --
%                          for long arcs this is the one big object, and the
%                          unrolled sums never need it materialised beyond the
%                          diagonals actually touched. Use only for testing.
%     'LocalBbar' N x nl double  the arc's DECORRELATED local design W^-T B
%                          (WP-33, fixes F-03). When given (with 'LocalLl'),
%                          the projector gains the EXACT local block: with
%                          Atil = A - B G the joint hat splits (leverage block
%                          identities) as Atil Sxx Atil' + B Nll^-1 B', so
%                              Rbar = Sigma^-1 - Rt Rt' - K K',
%                              K    = W^-1 (W^-T B) Ll^-T,   nl exact columns.
%                          Abar must then be the REDUCED AtilBar = W^-T(A-BG)
%                          (see the private jointArcParts helper) -- only that
%                          pairing probes the operator U factors. Default []:
%                          the pre-WP-33 common-only projector, VERBATIM
%                          (nl = 0 arcs are bitwise unchanged).
%     'LocalLl'  nl x nl double  the SYSTEM's lower Cholesky factor
%                          chol(Nll,'lower') (the factor LSQVCE.ENGINE.BLOCKCOFACTOR inverts);
%                          required with 'LocalBbar'.
%     'Nugget'   N x 1 double  (R7 WP-115 part 3b) Sigma = T + diag(Nugget):
%                          a stationary law T plus an additive per-row
%                          variance (the joint estimator's rejected epochs).
%                          W must then be the factor of the WHOLE Sigma; the
%                          sub-diagonal sums of Sigma^-1 come from T's
%                          Toeplitz fast path and a k x k Woodbury over the
%                          k rows with Nugget > 0, instead of the O(N^3)
%                          general fallback. Requires 'LawW'.
%     'LawW'     N x N double  the upper Cholesky factor of T alone (T = LawW'LawW);
%                          required with 'Nugget'.
%
%   OUTPUT
%     out   1x1 struct
%       .s        nMaxLag x 1  the per-lag redundancies s_n (n = 0 first). NOTE:
%                              s_n for n > 0 may be NEGATIVE -- it is a trace of
%                              an indefinite product, not a variance.
%       .omega    nMaxLag x 1  the per-lag residual sums Omega_n (NaN if ebar
%                              was empty)
%       .N        1x1          the arc length used
%       .probes   1x1          Z, the probe count implied by UiZ
%       .Rbar     N x N | []   dense Rbar if 'WantRbar' was true, else []
%
%   INVARIANTS (tested)
%     s_0     == trace(Rbar)
%     s_n     == trace(Rbar * V_n) for the explicit V_n (LSQVCE.ENGINE.COVCOMPONENT)
%     Omega_n == etilde' * V_n * etilde == ehat' Sigma^-T V_n Sigma^-1 ehat
%     Rbar    -> Sigma^-1 - Sigma^-1 A N^-1 A' Sigma^-1  as Z -> inf  (nl = 0)
%     Rbar    -> Sigma^-1 - Sigma^-1 [A B] Njoint^-1 [A B]' Sigma^-1 as Z ->
%                inf, when 'LocalBbar'/'LocalLl' are given with AtilBar
%                (v10review.py T3c: block realisation == joint inverse)
%
%   EXAMPLE
%     % Sigma = W'W for this arc; A its design; N = U'U the accumulated normals.
%     Zb  = sign(randn(P, 100))/sqrt(100);
%     UiZ = U \ Zb;                       % ONCE for all arcs
%     Abar = W' \ A;   ebar = W' \ ehat;  % per arc (as the accumulator whitened)
%     r = lsqvce.kernel.redundancyLags(Abar, W, UiZ, ebar, 60);
%     r.s(1)          % s_0
%
%   ERRORS
%     lsqvce:redundancyLags:badEbar      the mean-residual argument is malformed
%     lsqvce:redundancyLags:badW         the weight argument is malformed
%     lsqvce:redundancyLags:badUiZ       the probe matrix has the wrong shape
%     MATLAB:validators:mustBePositive / mustBeInteger  nMaxLag not a
%                   positive integer (WP-39: the scalar constraint sits
%                   in the typed signature; the old badNMaxLag body check
%                   is gone -- cross-argument shape checks below keep
%                   their lsqvce: identifiers)
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:redundancyLags:badLocal             LocalBbar has X rows but the arc has X.lsqvce:redundancyLags:badLocal...
%     lsqvce:redundancyLags:badNugget            Nugget / LawW do not match the arc (R7 WP-115 part 3b)
%
%   See also LSQVCE.ENGINE.COVCOMPONENT, LSQVCE.KERNEL.TOEPLITZCOFACTOR, LSQVCE.KERNEL.DCTI,
%            LSQVCE.ENGINE.TRACEESTIMATE, LSQVCE.ENGINE.STREAMCONTRIBUTION
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-131: the Toeplitz fast path
%   is taken only after a structure probe (Sigma y against toeplitz(c1) y);
%   three probe columns let a per-row sigma or a gap through. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-111 part 3: the general
%   (non-Toeplitz) Sigma^-1 lag sums run in BLAS-3 column blocks -- the
%   joint estimator now asks for every lag of an arc. With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-115 part 3b: 'Nugget' /
%   'LawW' -- a stationary law plus a per-row nugget keeps the structured
%   Sigma^-1 sums (Woodbury over the nugget rows). With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 3 (R7-24): the
%   retired leverageBlocks references point at lsqvce.engine.BlockCofactor.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) (F-07): typed
%   `arguments` signature, canonical UpperCamelCase option names (legacy
%   casings keep matching; enum VALUES case-sensitive lowercase). With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-37 (fixes F-09): the dense
%   X = W\eye(N) is gone from the main path -- Sigma^-1 sub-diagonal SUMS
%   come from Levinson + Gohberg-Semencul on Toeplitz (stationary) arcs
%   (O(N^2)/O(N)) and from an nLag-ring column stream of W^-T otherwise
%   (O(N*nLag) memory); 'WantRbar' alone still materialises. With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-18 (CEST, UTC+02:00), WP-33 (fixes F-03, CRITICAL): the
%   projector can now absorb the LOCAL parameter block exactly via
%   'LocalBbar'/'LocalLl' (see OPTIONS); the nl = 0 path is verbatim the old
%   code. With assistance from Claude (Anthropic), model Claude Fable 5.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+kernel (R4 decision R4-D1,
%   WP-64a): the KERNEL tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
arguments
    Abar (:,:) double
    W (:,:) double
    UiZ (:,:) double
    ebar double
    nMaxLag (1,1) double {mustBeInteger, mustBePositive}
    opt.WantRbar (1,1) logical = false
    opt.LocalBbar double = []
    opt.LocalLl double = []
    opt.Nugget double = []
    opt.LawW double = []
end

    % WP-39 stage (b): typed `arguments` signature (names were already
    % UpperCamelCase; NAME matching stays case-insensitive per MATLAB
    % semantics).
    o = opt;

    [N, P] = size(Abar);
    if size(W,1) ~= N || size(W,2) ~= N
        error('lsqvce:redundancyLags:badW', ...
              'W must be %d x %d to match Abar (got %d x %d).', ...
              N, N, size(W,1), size(W,2));
    end
    if size(UiZ,1) ~= P
        error('lsqvce:redundancyLags:badUiZ', ...
             ['UiZ has %d rows but Abar has %d columns: the shared probe block ' ...
              'must be U^-1*Zbar with U the factor of the P x P normal matrix.'], ...
              size(UiZ,1), P);
    end
    % WP-39 b1 r3: the positive-integer check lives in the typed signature
    % ({mustBeInteger, mustBePositive}); the old body check is dead and
    % deleted (rule 8).
    nMaxLag = min(nMaxLag, N);
    Z = size(UiZ, 2);

    % ---- the EXACT local block (WP-33, fixes F-03) ------------------------
    %   With the local design present, the joint hat splits by the leverage
    %   block identities (LSQVCE.ENGINE.BLOCKCOFACTOR) into
    %       [A B] N_joint^-1 [A B]' = Atil Sxx Atil' + B Nll^-1 B',
    %   Atil = A - B G. The FIRST term is what the shared probe covers when
    %   the caller passes AtilBar (see LSQVCE private jointArcParts); the
    %   SECOND is exact and cheap: K = W^-1 W^-T B Ll^-T has only nl columns,
    %   so its sub-diagonal row-dots cost the same O(N*nl) pattern as Rt.
    hasLocal = ~isempty(o.LocalBbar) && size(o.LocalBbar, 2) > 0;
    if hasLocal
        nl = size(o.LocalBbar, 2);
        if size(o.LocalBbar, 1) ~= N
            error('lsqvce:redundancyLags:badLocal', ...
                  'LocalBbar has %d rows but the arc has %d.', ...
                  size(o.LocalBbar,1), N);
        end
        if ~isequal(size(o.LocalLl), [nl nl])
            error('lsqvce:redundancyLags:badLocal', ...
                 ['LocalLl must be the nl x nl lower Cholesky factor of the ' ...
                  'system''s Nll (got %dx%d for nl = %d).'], ...
                  size(o.LocalLl,1), size(o.LocalLl,2), nl);
        end
        K = W \ (o.LocalBbar / o.LocalLl.');    % Sigma^-1-metric local block
    end

    if ~isempty(o.Nugget) && (numel(o.Nugget) ~= N || ~isequal(size(o.LawW), [N N]))
        error('lsqvce:redundancyLags:badNugget', ...
             ['Nugget must hold %d variances and LawW be the %d x %d upper Cholesky factor of ' ...
              'the law without it (got %d values and a %d x %d LawW).'], ...
              N, N, N, numel(o.Nugget), size(o.LawW, 1), size(o.LawW, 2));
    end

    % ---- Rtilde = W^-1 (Abar * UiZ)                              (eq. 6.5.48)
    %   W is upper triangular with Sigma = W'W, so W^-1 X is a backslash solve.
    %   With locals, Abar must be the REDUCED AtilBar = W^-T (A - B G): only
    %   then is the probed operator the one U factors (see the block identity
    %   above).
    Rt = W \ (Abar * UiZ);                       % N x Z

    % ---- s_n WITHOUT forming the dense Rbar                (eqs. 6.5.51-6.5.52)
    %   Rbar = Sigma^-1 - Rt*Rt', and we only ever need its n-th sub-diagonal:
    %       Rbar(i+n,i) = Sinv(i+n,i) - row_{i+n}(Rt) . row_i(Rt)
    %   With X = W^-1 (upper triangular) we have Sigma^-1 = X*X', so likewise
    %       Sinv(i+n,i)  = row_{i+n}(X) . row_i(X)
    %   Both are row dot-products: O(N*Z) and O(N^2_tri) work, never an extra
    %   dense N x N Rbar. This is what keeps a long arc affordable.
    % ---- Sigma^-1 sub-diagonal SUMS without the dense W^-1 (WP-37, F-09) --
    %   Only the SUM over each sub-diagonal of Sigma^-1 enters s_n, so the
    %   N x N inverse never needs to exist. Stationary (Toeplitz) arcs -- the
    %   case every covModel ACF/ARMA/PSD produces -- go through Levinson +
    %   Gohberg-Semencul: O(N^2) time, O(N) memory. A non-Toeplitz W falls
    %   back to blocks of 256 columns of Sigma^-1: O(N^3) time in BLAS-3
    %   (no structure, no shortcut), O(N*256) memory (R7 WP-111 part 3; the
    %   nMaxLag-deep ring buffer before it held N x nLag).
    if ~isempty(o.Nugget) && any(o.Nugget(:) > 0)
        dS = nuggetSubdiagSums(o.LawW, o.Nugget(:), nMaxLag);   % R7 WP-115 part 3b
    else
        dS = sigmaInvSubdiagSums(W, nMaxLag);    % dS(1) = trace(Sigma^-1)
    end
    s = zeros(nMaxLag, 1);
    if hasLocal
        % Rbar = Sinv - Rt*Rt' - K*K': one more row-dot per sub-diagonal.
        % R5 M3 3c: every sub-diagonal sum of Rt Rt' (and K K') is the lag-n
        % AUTOCORRELATION of the columns summed over columns -- ALL lags from
        % one FFT per column, O((Z+nl) N log N) instead of O((Z+nl) N nLag).
        % Exact to rounding (pinned against the loop form in test_ellmer).
        dRtAll = columnLagSums(Rt, nMaxLag);
        dKAll  = columnLagSums(K,  nMaxLag);
        s(1)   = dS(1) - dRtAll(1) - dKAll(1);
        s(2:nMaxLag) = 2 * (dS(2:nMaxLag) - dRtAll(2:nMaxLag) - dKAll(2:nMaxLag));
    else
        dRtAll = columnLagSums(Rt, nMaxLag);
        s(1)   = dS(1) - dRtAll(1);
        s(2:nMaxLag) = 2 * (dS(2:nMaxLag) - dRtAll(2:nMaxLag));
    end

    % ---- Omega_n from the whitened residual              (eqs. 6.5.54-6.5.56)
    omega = nan(nMaxLag, 1);
    if ~isempty(ebar)
        ebar = ebar(:);
        if numel(ebar) ~= N
            error('lsqvce:redundancyLags:badEbar', ...
                  'ebar has %d elements but the arc has %d rows.', numel(ebar), N);
        end
        et = W \ ebar;                           % etilde = Sigma^-1 ehat
        acf = columnLagSums(et, nMaxLag);        % R5 M3 3c: one FFT, all lags
        omega(1) = acf(1);                       % n = 0
        omega(2:nMaxLag) = 2 * acf(2:nMaxLag);
    end

    out = struct('s', s, 'omega', omega, 'N', N, 'probes', Z, 'Rbar', []);
    if logical(o.WantRbar)
        % Testing/diagnostic path ONLY: this is the dense N x N object the main
        % path deliberately avoids (and, since WP-37, never builds otherwise).
        X = W \ eye(N);                          % W^-1, HERE only
        out.Rbar = X * X.' - Rt * Rt.';
        if hasLocal
            out.Rbar = out.Rbar - K * K.';
        end
    end
end

%==========================================================================
function dS = sigmaInvSubdiagSums(W, nL)
%SIGMAINVSUBDIAGSUMS  Sub-diagonal sums of Sigma^-1 = (W'W)^-1, structured.
%   dS(n+1) = sum_i Sigma^-1(i+n, i), n = 0..nL-1 (dS(1) = trace).
%   Toeplitz Sigma (detected from three probe columns): Levinson-Durbin for
%   x = Sigma^-1 e_1, then the Gohberg-Semencul split
%       Sigma^-1 = (1/x_1) * (L(x) L(x)' - L(z) L(z)'),  z = [0; x(N:-1:2)],
%   whose triangular-Toeplitz sub-diagonal sums are weighted correlations:
%       S_n(u) = sum_j u(j+n) u(j) (N-n-j+1),   O(N) per lag.
%   (Validated against the dense inverse to 3e-16, WP-37 round note.)
%   General W: columns of Sigma^-1 in blocks of 256, each column's entries on
%   and below the diagonal added to the lag they sit on.
    N  = size(W, 1);
    dS = zeros(nL, 1);

    % Sigma columns 1, 2 and mid, each O(N)/O(N^2) from the triangular W:
    c1 = W.' * W(:, 1);
    isToep = N >= 3;
    if isToep
        c2  = W.' * W(:, 2);
        jm  = max(3, ceil(N/2));
        cm  = W.' * W(:, jm);
        tol = 1e-10 * max(abs(c1));
        isToep = max(abs(c2(2:N) - c1(1:N-1))) <= tol && ...
                 abs(c2(1) - c1(2)) <= tol && ...
                 max(abs(cm(jm:N) - c1(1:N-jm+1))) <= tol && ...
                 max(abs(cm(1:jm) - c1(jm:-1:1))) <= tol;
    end
    if isToep
        % WP-131: three columns cannot see a change of structure that misses
        % them (a per-row sigma after N/2, a down-weighted epoch, a gap): the
        % fast path then summed the inverse of toeplitz(c1), not of Sigma (+29 %
        % trace; dev/validation/wp131_kernel_mirror.py). One deterministic
        % probe, O(N^2) like the columns: Sigma y against the Toeplitz product.
        y  = sin(sqrt(2) * (1:N).'.^1.5);              % generic, no RNG state
        Sy = W.' * (W * y);
        isToep = norm(Sy - toeplitzTimes(c1, y)) <= 1e-10 * norm(abs(W.') * (abs(W) * abs(y)));
    end

    if isToep
        x  = toeplitzInvFirstCol(c1);            % Sigma^-1 e_1, O(N^2)
        z  = [0; x(N:-1:2)];
        x1 = x(1);
        for n = 0:nL-1
            j = (1:N-n).';
            wgt = (N - n - j + 1);
            dS(n+1) = (sum(x(j+n) .* x(j) .* wgt) - ...
                       sum(z(j+n) .* z(j) .* wgt)) / x1;
        end
        return
    end

    % General fallback (R7 WP-111 part 3): Sigma^-1 in blocks of columns,
    % Sigma^-1 E = W \ (W' \ E), each column's sub-diagonal added to dS --
    % BLAS-3, memory N x 256, cost independent of nL. The former ring buffer
    % of rows of W^-1 held N x nL and took N^2 nL / 2 scalar dot products,
    % i.e. O(N^3) loop iterations once the joint estimator asks for every lag.
    nb = 256;
    for j0 = 1:nb:N
        jj = j0:min(j0 + nb - 1, N);
        E  = zeros(N, numel(jj));
        E(sub2ind(size(E), jj, 1:numel(jj))) = 1;
        Si = W \ (W.' \ E);                     % columns jj of Sigma^-1
        for c = 1:numel(jj)
            n = min(nL, N - jj(c) + 1);
            dS(1:n) = dS(1:n) + Si(jj(c):jj(c) + n - 1, c);
        end
    end
end

function dS = nuggetSubdiagSums(WT, c, nL)
%NUGGETSUBDIAGSUMS  Sub-diagonal sums of (T + diag(c))^-1, T = WT'WT stationary.
%   Woodbury over the k rows z with c > 0 (U their unit columns):
%       Sigma^-1 = T^-1 - H K^-1 H',  H = T^-1 U,  K = diag(1/c) + U' T^-1 U,
%   so the sums are T's (Toeplitz fast path) minus the lag sums of
%   H LK^-T (K = LK LK'), one FFT per column. K holds 1/c, so an epoch with
%   c ~ 1e12 T_ii is deleted to rounding (mirrors/wp115/nugget.py: 1e-16).
    N = size(WT, 1);
    z = find(c > 0);  k = numel(z);
    E = zeros(N, k);
    E(sub2ind([N, k], z(:), (1:k).')) = 1;
    Y = WT.' \ E;                                % WT^-T U
    H = WT \ Y;                                  % T^-1 U
    K = diag(1 ./ c(z)) + Y.' * Y;
    LK = chol((K + K.') / 2, 'lower');
    dS = sigmaInvSubdiagSums(WT, nL) - columnLagSums(H / LK.', nL);
end

function x = toeplitzInvFirstCol(c)
%TOEPLITZINVFIRSTCOL  First column of inv(toeplitz(c)) via Levinson-Durbin.
%   For SPD Sigma = toeplitz(c): Durbin's recursion solves the Yule-Walker
%   system T_{N-1} a = -r(2:N) in O(N^2)/O(N); the normal-equation identity
%   toeplitz(r) * [1; a] = E * e_1 then gives x = [1; a] / (E * c(1)).
    c = c(:);
    N = numel(c);
    if N == 1, x = 1 / c(1); return, end
    r = c / c(1);                                % r(1) = 1
    a = -r(2);
    E = 1 - a^2;
    for n = 2:N-1
        k = -(r(n+1) + a.' * r(n:-1:2)) / E;
        a = [a + k * a(end:-1:1); k];
        E = E * (1 - k^2);
    end
    x = [1; a] / (E * c(1));
end

%==========================================================================
function d = columnLagSums(X, nL)
%COLUMNLAGSUMS  d(n+1) = sum_c sum_i X(i+n, c) X(i, c), n = 0..nL-1, by FFT.
%   The lag-n sub-diagonal sum of X X' equals the sum over columns of each
%   column's lag-n autocorrelation, so all lags come from one real FFT per
%   column of length >= 2N-1 (zero-padded, no circular wrap): the power
%   spectrum summed over columns, inverse-transformed. O(C N log N) for C
%   columns instead of O(C N nL). Exact to rounding; R5 M3 3c (WP-83).
%   INPUT  X N x C double;  nL 1x1 lags wanted (<= N)
%   OUTPUT d nL x 1 double
    [N, C] = size(X);
    d = zeros(nL, 1);
    if N == 0 || C == 0, return, end
    nfft = 2^nextpow2(2*N - 1);
    F = fft(X, nfft, 1);                          % columns, zero-padded
    P = sum(real(F .* conj(F)), 2);               % summed power spectrum
    r = ifft(P, 'symmetric');                     % lag sums, lag 0 first
    d = r(1:nL);
end

function z = toeplitzTimes(c, y)
%TOEPLITZTIMES  toeplitz(c) * y for a symmetric Toeplitz matrix, O(N log N) by
%   the circulant embedding of size 2N (WP-131: the structure probe of
%   sigmaInvSubdiagSums; N x N would cost N^2 memory per arc).
    N = numel(c);
    e = [c(:); 0; c(N:-1:2)];
    z = ifft(fft(e) .* fft([y(:); zeros(N, 1)]), 'symmetric');
    z = z(1:N);
end

