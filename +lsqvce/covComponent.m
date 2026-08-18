function y = covComponent(x, n)
%LSQVCE.COVCOMPONENT  Apply the lag-n covariance component V_n to a vector/block.
%
%   y = LSQVCE.COVCOMPONENT(x, n)
%
%   Returns y = V_n * x WITHOUT forming V_n, where V_n is Ellmer's covariance
%   component matrix (2018, eq. 6.5.16): the N x N identity for n = 0, and for
%   n > 0 the matrix with ones ONLY on the +n and -n diagonals (zeros elsewhere).
%   These are the building blocks of the Toeplitz cofactor matrix,
%
%       Q_ll = sum_{n=0}^{Nmax-1} C^n_xx * V_n            (Ellmer eq. 6.5.17)
%
%   and appear in the per-lag redundancy s_n = tr(R V_n) of the spectral VCE.
%   The apply is BANDED, O(N*K) per component and O(1) extra memory:
%
%       (V_n x)_i = x_{i+n} + x_{i-n}     (edge terms dropped at the boundaries)
%
%   so a length-N arc never materialises an N x N matrix -- essential when N is a
%   full short-arc length (thousands of epochs).
%
%   INPUT
%     x    N x K double   K columns of length N (a row vector is treated as one
%                         column). N >= 1.
%     n    1x1 int        lag index, 0 <= n < N. n = 0 returns x unchanged.
%
%   OUTPUT
%     y    N x K double   V_n * x, same shape as x (row in -> row out)
%
%   INVARIANTS (tested)
%     covComponent(x,0) == x
%     covComponent(x,n) == Vn_dense * x  for the explicit V_n
%     sum_n C(n)*covComponent(x,n) == toeplitz(C) * x
%     trace via probes: E[z .* covComponent(z,n)] -> diag(V_n) (0 for n>0)
%
%   EXAMPLE
%     x = (1:6).';
%     lsqvce.covComponent(x, 2)      % = [x3+0; x4+0; x5+x1; x6+x2; 0+x3; 0+x4]
%
%   ERRORS
%     lsqvce:covComponent:badLag    the lag is not a non-negative integer
%     lsqvce:covComponent:lagTooLarge  the lag exceeds the series length
%
%   See also LSQVCE.TOEPLITZCOFACTOR, LSQVCE.DCTI, LSQVCE.REDUNDANCYLAGS
%
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.

    wasRow = isrow(x);
    if wasRow, x = x(:); end
    N = size(x, 1);

    if ~(isscalar(n) && n == floor(n) && n >= 0)
        error('lsqvce:covComponent:badLag', ...
              'The lag n must be a non-negative integer (got %g).', n);
    end
    if n >= N
        error('lsqvce:covComponent:lagTooLarge', ...
              'Lag n = %d must be < N = %d.', n, N);
    end

    if n == 0
        y = x;
    else
        y = zeros(size(x));
        y(1:N-n, :) = y(1:N-n, :) + x(n+1:N, :);     % +n diagonal
        y(n+1:N, :) = y(n+1:N, :) + x(1:N-n, :);     % -n diagonal
    end

    if wasRow, y = y.'; end
end
