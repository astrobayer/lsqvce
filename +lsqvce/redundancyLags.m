function out = redundancyLags(Abar, W, UiZ, ebar, nMaxLag, varargin)
%LSQVCE.REDUNDANCYLAGS  Per-lag redundancies s_n and residual sums Omega_n.
%
%   out = LSQVCE.REDUNDANCYLAGS(Abar, W, UiZ, ebar, nMaxLag)
%   out = LSQVCE.REDUNDANCYLAGS(..., 'Name', value, ...)
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
%                          chol(Nll,'lower') (as LSQVCE.LEVERAGEBLOCKS.Ll);
%                          required with 'LocalBbar'.
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
%     s_n     == trace(Rbar * V_n) for the explicit V_n (LSQVCE.COVCOMPONENT)
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
%     r = lsqvce.redundancyLags(Abar, W, UiZ, ebar, 60);
%     r.s(1)          % s_0
%
%   ERRORS
%     lsqvce:redundancyLags:badEbar      the mean-residual argument is malformed
%     lsqvce:redundancyLags:badW         the weight argument is malformed
%     lsqvce:redundancyLags:badUiZ       the probe matrix has the wrong shape
%     lsqvce:redundancyLags:badNMaxLag   the lag count is not a positive integer
%
%   See also LSQVCE.COVCOMPONENT, LSQVCE.TOEPLITZCOFACTOR, LSQVCE.DCTI,
%            LSQVCE.TRACEESTIMATE, LSQVCE.STREAMCONTRIBUTION
%
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

    p = inputParser;  p.FunctionName = 'lsqvce.redundancyLags';
    p.addParameter('WantRbar', false, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
    p.addParameter('LocalBbar', [], @(v) isempty(v) || isnumeric(v));
    p.addParameter('LocalLl',   [], @(v) isempty(v) || isnumeric(v));
    p.parse(varargin{:});
    o = p.Results;

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
    if ~(isscalar(nMaxLag) && nMaxLag == floor(nMaxLag) && nMaxLag >= 1)
        error('lsqvce:redundancyLags:badNMaxLag', ...
              'nMaxLag must be a positive integer (got %g).', nMaxLag);
    end
    nMaxLag = min(nMaxLag, N);
    Z = size(UiZ, 2);

    % ---- the EXACT local block (WP-33, fixes F-03) ------------------------
    %   With the local design present, the joint hat splits by the leverage
    %   block identities (LSQVCE.LEVERAGE) into
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
    %   back to streaming the columns of W^-T with an nMaxLag-deep ring
    %   buffer: still O(N^3) time (no structure, no shortcut) but O(N*nLag)
    %   memory -- which makes the memory claim in the help TRUE at last.
    dS = sigmaInvSubdiagSums(W, nMaxLag);        % dS(1) = trace(Sigma^-1)
    s = zeros(nMaxLag, 1);
    if hasLocal
        % Rbar = Sinv - Rt*Rt' - K*K': one more row-dot per sub-diagonal.
        s(1) = dS(1) - sum(Rt(:).^2) - sum(K(:).^2);
        for n = 1:nMaxLag-1
            idx = 1:(N-n);
            dRt   = sum(Rt(idx+n, :) .* Rt(idx, :), 2);
            dK    = sum(K(idx+n, :)  .* K(idx, :),  2);
            s(n+1) = 2 * (dS(n+1) - sum(dRt) - sum(dK));
        end
    else
        s(1) = dS(1) - sum(Rt(:).^2);
        for n = 1:nMaxLag-1
            idx = 1:(N-n);
            dRt   = sum(Rt(idx+n, :) .* Rt(idx, :), 2);  % (Rt Rt') sub-diagonal
            s(n+1) = 2 * (dS(n+1) - sum(dRt));
        end
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
        omega(1) = et.' * et;                    % n = 0
        for n = 1:nMaxLag-1
            omega(n+1) = 2 * (et(n+1:N).' * et(1:N-n));
        end
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
%   General W: stream the columns of W^-T (each is a row of W^-1) with an
%   nL-deep ring buffer and accumulate the offset dot products.
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

    % General fallback: r_j = W' \ e_j is row j of W^-1; keep the last nL.
    ring = zeros(N, nL);
    e    = zeros(N, 1);
    for j = 1:N
        e(j) = 1;
        r = W.' \ e;                             % forward solve, O(N^2)
        e(j) = 0;
        slot = mod(j-1, nL) + 1;
        ring(:, slot) = r;
        dS(1) = dS(1) + r.' * r;
        for n = 1:min(nL-1, j-1)
            dS(n+1) = dS(n+1) + r.' * ring(:, mod(j-1-n, nL) + 1);
        end
    end
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

