function [Aw, Bw, lw] = armaWhiten(a, b, s, r, A, B, l, isFirstBatch)
%LSQVCE.KERNEL.ARMAWHITEN  Apply an ARMA whitening filter column-wise (exact).
%
%   [Aw,Bw,lw] = LSQVCE.KERNEL.ARMAWHITEN(a, b, s, r, A, B, l, isFirstBatch)
%
%   Applies the operator F with F'*F = inv(C), where C is the covariance of a
%   process with ARMA(p,q) structure (a, b) and innovation std s. Costs
%   O(m*(p+q)) and never forms the m x m covariance.
%
%   INPUT
%     a             1 x (p+1)  AR polynomial [1 a1 ... ap]
%     b             1 x (q+1)  MA polynomial [1 b1 ... bq]; use 1 for pure AR
%     s             scalar     innovation std, > 0
%     r             >= p x 1   autocovariance r(0..), already scaled by s^2; with
%                              an MA part its length K is the exact startup:
%                              pass it over the filter memory p + N, N the
%                              decay of 1/B(z) (LSQVCE.ENGINE.COVMODEL does)
%     A, B          m x *      design blocks (either may be m x 0 or sparse; sparse in, sparse out)
%     l             m x 1      observation vector
%     isFirstBatch  logical    true for the first batch of an arc
%
%   OUTPUT
%     Aw, Bw, lw    whitened blocks, same sizes as the inputs
%
%   ARGUMENT ORDER -- THE CLASSIC TRAP
%     The filter applied is filter(a, b, x), i.e. NUMERATOR a, DENOMINATOR b.
%     That is the whitening filter A(z)/B(z). Swapping the arguments gives the
%     COLOURING filter B(z)/A(z), which is silently wrong: the result still
%     looks like a residual series but the normal equations are corrupted.
%
%   EXACTNESS
%     This is not an approximation. The operator F satisfies
%         norm(F'*F - inv(C)) / norm(inv(C))  <  1e-13
%     (verified for n = 20..300).
%
%     In fact F IS the Cholesky inverse factor, not merely another square root
%     of inv(C): both are lower triangular, and the triangular square root is
%     unique up to sign, so
%         norm(F - chol(C)^-1) / norm(chol(C)^-1)  ~  1e-14
%     The FIR filter reproduces the dense decorrelation exactly, row for row,
%     at O(m*p) instead of O(m^3) and with no matrix ever formed. The steady-
%     state rows ARE the AR coefficients reversed and divided by s; the leading
%     p rows are the stationary startup block below.
%
%   THE STARTUP BLOCK
%     A zero-state FIR filter mishandles the first p rows: it implicitly
%     assumes x(k) = 0 for k < 1. Those rows are replaced by the exact
%     stationary whitener chol(toeplitz(r(1:p)))^-1, which is the leading
%     p x p block of chol(C)^-1. Only done when isFirstBatch is true.
%     With an MA part the recursion 1/B(z) carries state, so the exact block
%     covers K = numel(r) rows and the recursion continues from their exact
%     innovations u_k = L(k,k) w_k (WP-127); for K >= p + the decay of 1/B(z)
%     to eps the operator is chol(C)^-1 to ~1e-13. Cost O(K^3) once per arc.
%
%   BATCH OVERLAP IS MANDATORY
%     On any batch after the first, the filter restarts from zero state and
%     its first rows are transient garbage. The caller must overlap batches
%     by the filter's memory and discard the overlap: p rows for a pure AR
%     (b = 1); with an MA part, p plus the samples after which the impulse
%     response of 1/B(z) has decayed below eps (WP-134: the help said p, which
%     holds for the AR part only). With that overlap the batched result is
%     identical to filtering the whole arc; with overlap = 0 the relative
%     error in the solution is ~1e-3 instead of ~1e-15.
%     LSQVCE.ENGINE.BLOCKACCUMULATOR does this automatically from
%     arc.filterMemory, which LSQVCE.ENGINE.COVMODEL computes this way.
%
%   EXAMPLE
%     a = [1 -1.6 0.9 -0.25];  s = 0.7;  n = 200;
%     r = s^2 * lsqvce.kernel.armaACV(a, 1, n);
%     X = randn(n, 4);
%     z = zeros(n, 0);
%     [~,~,Y] = lsqvce.kernel.armaWhiten(a, 1, s, r, z, z, X, true);
%     W = lsqvce.kernel.decorrFactor(toeplitz(r));
%     norm(Y.'*Y - (W*X).'*(W*X)) / norm((W*X).'*(W*X))    % ~ 1e-14
%
%   ERRORS
%     lsqvce:armaWhiten:badSigma      the innovation std is not a positive scalar
%     lsqvce:armaWhiten:acvTooShort   the autocovariance is shorter than the
%                                     model order needs for its startup block
%     lsqvce:armaWhiten:sizeMismatch  A and l disagree on their row count
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:armaWhiten:startupNotPD             [warning] The X x X stationary startup block is not positive definite; keeping...
%
%   See also LSQVCE.KERNEL.DECORRFACTOR, LSQVCE.ENGINE.COVMODEL, LSQVCE.KERNEL.LEVINSONDURBIN
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-134: the overlap note names the
%   MA part's memory (p plus the decay of 1/B(z)). With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+kernel (R4 decision R4-D1,
%   WP-64a): the KERNEL tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
%
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-127: with an MA part the exact
%   startup covers numel(r) rows and the recursion continues from their exact
%   innovations (the rows after p were not chol(C)^-1). With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%
%   Sparse blocks: 2026-10-12 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-123.

    narginchk(8,8);
    % WP-123: sparse blocks (makeLoader allows a sparse B) run the dense
    % whitener on their non-zero columns -- N-D reshape and filter refuse
    % sparse input. The dense path below is unchanged.
    if issparse(A) || issparse(B) || issparse(l)
        dense = @(A, B, l) lsqvce.kernel.armaWhiten(a, b, s, r, A, B, l, isFirstBatch);
        [Aw, Bw, lw] = activeColumns(dense, A, B, l);
        return
    end
    a = a(:).';   b = b(:).';   r = r(:);

    p = numel(a) - 1;
    m = size(l, 1);

    if s <= 0 || ~isfinite(s)
        error('lsqvce:armaWhiten:badSigma', 'Innovation std s = %g must be positive.', s);
    end
    if size(A,1) ~= m && ~isempty(A)
        error('lsqvce:armaWhiten:sizeMismatch', ...
              'A has %d rows but l has %d.', size(A,1), m);
    end
    if size(B,1) ~= m && ~isempty(B)
        error('lsqvce:armaWhiten:sizeMismatch', ...
              'B has %d rows but l has %d.', size(B,1), m);
    end

    % --- steady-state filter: NUMERATOR a, DENOMINATOR b -------------------
    Aw = filter(a, b, A, [], 1) / s;
    Bw = filter(a, b, B, [], 1) / s;
    lw = filter(a, b, l, [], 1) / s;

    % --- exact stationary startup ------------------------------------------
    hasMA = numel(b) > 1 && any(b(2:end) ~= 0);
    if ~isFirstBatch, return, end
    if hasMA
        % WP-127: with an MA part the recursion 1/B(z) carries STATE. Replacing
        % the leading p rows alone left that state built from the zero-state
        % transient, so the rows after p were not chol(C)^-1 (F'F off inv(C)
        % by 1.6 % at a = [1 -0.5], b = [1 0.4]; whitened variance 10.7 at row 2
        % for b = [1 0.8]). The leading K = numel(r) rows -- the caller passes
        % the ACV over the filter memory p + the decay of 1/B(z) -- are
        % whitened exactly, and the recursion continues from their EXACT
        % innovations u_k = L(k,k) w_k: beyond K the innovations-algorithm
        % coefficients equal b to rounding. Mirror: dev/validation/
        % wp127_arma_whiten_mirror.py, |F - chol(C)^-1| / |.| <= 3e-13.
        K = min(numel(r), m);
        if K < min(p, m)
            error('lsqvce:armaWhiten:acvTooShort', ...
                ['The autocovariance has %d lags but the startup block needs ' ...
                 'at least %d. Pass r over the filter memory.'], numel(r), min(p, m));
        end
        [Lc, fail] = chol(toeplitz(r(1:K)), 'lower');
        if fail ~= 0
            warnStartup(K);
            return
        end
        Aw = exactMa(a, b, s, Lc, A, Aw);
        Bw = exactMa(a, b, s, Lc, B, Bw);
        lw = exactMa(a, b, s, Lc, l, lw);
        return
    end
    if p > 0
        nS = min(p, m);
        if numel(r) < nS
            error('lsqvce:armaWhiten:acvTooShort', ...
                ['The autocovariance has %d lags but the startup block needs ' ...
                 '%d. Pass r(0..p).'], numel(r), nS);
        end
        G = toeplitz(r(1:nS));
        [Lc, fail] = chol(G, 'lower');
        if fail == 0
            Wp = Lc \ eye(nS);
            if ~isempty(A), Aw(1:nS,:) = Wp * A(1:nS,:); end
            if ~isempty(B), Bw(1:nS,:) = Wp * B(1:nS,:); end
            lw(1:nS,:) = Wp * l(1:nS,:);
        else
            warnStartup(nS);
        end
    end
end

%==========================================================================
function Y = exactMa(a, b, s, Lc, X, Y)
%EXACTMA  Rows 1..K exact (Lc \ X), rows K+1.. the recursion 1/B(z) seeded
%   with the exact innovations: filter(1, b, [filter(b, 1, u); (a * X)(K+1:m)])
%   reproduces u on rows 1..K and continues u_k = (a*X)_k - sum_j b_j u_(k-j).
%   INPUT  Lc K x K lower chol of toeplitz(r(1:K)); X m x n a raw block;
%          Y m x n its zero-state output (overwritten)
    if isempty(X), return, end
    K = size(Lc, 1);
    W = Lc \ X(1:K, :);
    Y(1:K, :) = W;
    if K == size(X, 1), return, end
    u = diag(Lc) .* W;                              % the exact innovations
    v = filter(a, 1, X, [], 1);                     % the AR part: exact for k > p
    z = filter(1, b, [filter(b, 1, u, [], 1); v(K+1:end, :)], [], 1);
    Y(K+1:end, :) = z(K+1:end, :) / s;
end

%==========================================================================
function warnStartup(nS)
%WARNSTARTUP  The stationary startup block is not positive definite.
    lsqvce.util.log.warn('lsqvce:armaWhiten:startupNotPD', ...
        ['The %d x %d stationary startup block is not positive ' ...
         'definite; keeping the filter transient instead. The first ' ...
         '%d rows of this arc are slightly mis-weighted.'], nS, nS, nS);
end
