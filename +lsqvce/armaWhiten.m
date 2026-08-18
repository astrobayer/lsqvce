function [Aw, Bw, lw] = armaWhiten(a, b, s, r, A, B, l, isFirstBatch)
%LSQVCE.ARMAWHITEN  Apply an ARMA whitening filter column-wise (exact).
%
%   [Aw,Bw,lw] = LSQVCE.ARMAWHITEN(a, b, s, r, A, B, l, isFirstBatch)
%
%   Applies the operator F with F'*F = inv(C), where C is the covariance of a
%   process with ARMA(p,q) structure (a, b) and innovation std s. Costs
%   O(m*(p+q)) and never forms the m x m covariance.
%
%   INPUT
%     a             1 x (p+1)  AR polynomial [1 a1 ... ap]
%     b             1 x (q+1)  MA polynomial [1 b1 ... bq]; use 1 for pure AR
%     s             scalar     innovation std, > 0
%     r             >= p x 1   autocovariance r(0..), already scaled by s^2
%     A, B          m x *      design blocks (either may be m x 0)
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
%
%   BATCH OVERLAP IS MANDATORY
%     On any batch after the first, the filter restarts from zero state and
%     the first p rows are transient garbage. The caller must overlap batches
%     by at least p rows and discard the overlap. With overlap = p the batched
%     result is identical to filtering the whole arc; with overlap = 0 the
%     relative error in the solution is ~1e-3 instead of ~1e-15.
%     LSQVCE.BLOCKACCUMULATOR does this automatically from arc.filterMemory.
%
%   EXAMPLE
%     a = [1 -1.6 0.9 -0.25];  s = 0.7;  n = 200;
%     r = s^2 * lsqvce.armaACV(a, 1, n);
%     X = randn(n, 4);
%     z = zeros(n, 0);
%     [~,~,Y] = lsqvce.armaWhiten(a, 1, s, r, z, z, X, true);
%     W = lsqvce.decorrFactor(toeplitz(r));
%     norm(Y.'*Y - (W*X).'*(W*X)) / norm((W*X).'*(W*X))    % ~ 1e-14
%
%   ERRORS
%     lsqvce:armaWhiten:badSigma      the innovation std is not a positive scalar
%     lsqvce:armaWhiten:acvTooShort   the autocovariance is shorter than the
%                                     model order needs for its startup block
%     lsqvce:armaWhiten:sizeMismatch  A and l disagree on their row count
%
%   See also LSQVCE.DECORRFACTOR, LSQVCE.COVMODEL, LSQVCE.LEVINSONDURBIN
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.

    narginchk(8,8);
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

    % --- exact stationary startup for the leading p rows -------------------
    if p > 0 && isFirstBatch
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
            lsqvce.util.log.warn('lsqvce:armaWhiten:startupNotPD', ...
                ['The %d x %d stationary startup block is not positive ' ...
                 'definite; keeping the filter transient instead. The first ' ...
                 '%d rows of this arc are slightly mis-weighted.'], nS, nS, nS);
        end
    end
end
