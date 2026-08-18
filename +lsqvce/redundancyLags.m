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
%   O(N) structured form; not exploited yet.) Pass 'WantRbar' only for testing.
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
%                             triangular, as from chol(Sigma)). May also be given
%                             as a struct with fields .W (see 'Solver' below).
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
%     Rbar    -> Sigma^-1 - Sigma^-1 A N^-1 A' Sigma^-1  as Z -> inf
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

    p = inputParser;  p.FunctionName = 'lsqvce.redundancyLags';
    p.addParameter('WantRbar', false, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
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

    % ---- Rtilde = W^-1 (Abar * UiZ)                              (eq. 6.5.48)
    %   W is upper triangular with Sigma = W'W, so W^-1 X is a backslash solve.
    Rt = W \ (Abar * UiZ);                       % N x Z

    % ---- s_n WITHOUT forming the dense Rbar                (eqs. 6.5.51-6.5.52)
    %   Rbar = Sigma^-1 - Rt*Rt', and we only ever need its n-th sub-diagonal:
    %       Rbar(i+n,i) = Sinv(i+n,i) - row_{i+n}(Rt) . row_i(Rt)
    %   With X = W^-1 (upper triangular) we have Sigma^-1 = X*X', so likewise
    %       Sinv(i+n,i)  = row_{i+n}(X) . row_i(X)
    %   Both are row dot-products: O(N*Z) and O(N^2_tri) work, never an extra
    %   dense N x N Rbar. This is what keeps a long arc affordable.
    X = W \ eye(N);                              % N x N upper triangular = W^-1
    s = zeros(nMaxLag, 1);
    % n = 0 is trivial: trace(Rbar) = trace(Sinv) - trace(Rt*Rt')
    s(1) = sum(X(:).^2) - sum(Rt(:).^2);         % tr(XX') = sum(X.^2)
    for n = 1:nMaxLag-1
        idx = 1:(N-n);
        dSinv = sum(X(idx+n, :) .* X(idx, :), 2);        % Sinv sub-diagonal
        dRt   = sum(Rt(idx+n, :) .* Rt(idx, :), 2);      % (Rt Rt') sub-diagonal
        s(n+1) = 2 * sum(dSinv - dRt);
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
        % path deliberately avoids. Sigma^-1 = X*X' with X = W^-1.
        out.Rbar = X * X.' - Rt * Rt.';
    end
end
