function G = varACV(A, Sigma, K)
%LSQVCE.VARACV  Cross-autocovariance of a stationary VAR(p) process.
%
%   G = LSQVCE.VARACV(A, Sigma, K)
%
%   For  z_t = sum_{k=1}^p A_k z_{t-k} + e_t,  C{e} = Sigma, returns
%   Gamma(j) = E[z_t z_{t-j}'] for j = 0..K as a d x d x (K+1) array.
%
%   METHOD
%     Companion (state-space) form X_t = F X_{t-1} + G e_t; the stationary
%     state covariance solves the discrete Lyapunov equation
%         Gx = F Gx F' + Q
%     by vectorisation -- a (p*d)^2 linear solve, trivial for the d = 3,
%     p <= 6 cases this exists for. Gamma(0..p-1) are read off Gx; higher
%     lags follow the Yule-Walker recursion Gamma(j) = sum A_k Gamma(j-k).
%
%   VERIFIED (validation/v6verify.py)
%     * Gamma(0) - sum A_k Gamma(k)' = Sigma and the recursion hold to 1e-11.
%     * Matches the empirical ACV of a 4e5-sample simulation to a few percent.
%
%   INPUT
%     A       d x d x p   coefficient matrices A_1..A_p
%     Sigma   d x d       innovation covariance (symmetric positive definite)
%     K       scalar      highest lag
%
%   STATIONARITY IS CHECKED
%     The companion spectral radius must be < 1, or there is no stationary
%     covariance to compute.
%
%   OUTPUT
%     G  d x d x (K+1)  Gamma(0..K), the cross-autocovariance sequence
%
%   ERRORS
%     lsqvce:varACV:badA           the coefficient array is not d x d x p
%     lsqvce:varACV:notStationary  the VAR process is not stationary, so no
%                                  stationary cross-ACV exists
%
%   See also LSQVCE.WHITTLE, LSQVCE.VARWHITEN, LSQVCE.ARMAACV
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.

    narginchk(3,3);
    [d, d2, p] = size(A);
    if d ~= d2
        error('lsqvce:varACV:badA', 'A must be d x d x p (got %d x %d x %d).', d, d2, p);
    end
    validateattributes(Sigma, {'numeric'}, {'size',[d d]}, 'lsqvce.varACV', 'Sigma', 2);
    validateattributes(K, {'numeric'}, {'scalar','integer','nonnegative'}, ...
                       'lsqvce.varACV', 'K', 3);
    Sigma = (Sigma + Sigma.')/2;

    % ---- companion form -----------------------------------------------
    n = p*d;
    F = zeros(n);
    F(1:d, :) = reshape(A, d, n);
    if p > 1
        F(d+1:end, 1:end-d) = eye((p-1)*d);
    end

    rho = max(abs(eig(F)));
    if rho >= 1
        error('lsqvce:varACV:notStationary', ...
              ['The VAR is not stationary: companion spectral radius %.4f ' ...
               '>= 1. No stationary covariance exists.'], rho);
    end

    % ---- discrete Lyapunov by vectorisation ----------------------------
    Q = zeros(n);  Q(1:d, 1:d) = Sigma;
    Gx = reshape((eye(n^2) - kron(F, F)) \ Q(:), n, n);
    Gx = (Gx + Gx.')/2;

    % ---- read off lags 0..p-1, recurse beyond --------------------------
    G = zeros(d, d, K+1);
    for j = 0:min(p-1, K)
        G(:,:,j+1) = Gx(1:d, j*d + (1:d));
    end
    G(:,:,1) = (G(:,:,1) + G(:,:,1).')/2;
    for j = max(p,1):K
        S = zeros(d);
        for k = 1:p
            S = S + A(:,:,k) * G(:,:,j-k+1);
        end
        G(:,:,j+1) = S;
    end
end
