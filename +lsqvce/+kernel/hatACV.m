function K = hatACV(L, m, X)
%LSQVCE.KERNEL.HATACV  The ACV of every bin of a piecewise-linear spectrum: r = K * S.
%
%   K = LSQVCE.KERNEL.HATACV(L, m)
%   Q = LSQVCE.KERNEL.HATACV(L, m, X)       % = K' * X, without forming K
%
%   The spectral law of the joint estimator (R7 WP-111 part 3) is a NON-NEGATIVE
%   spectrum S on L equispaced bins f_j = j / (2(L-1)), j = 0..L-1 (cycles per
%   sample), linear in between. Its autocovariance at sample lags 0..m-1 is
%   r = K * S, where column j of K is the exact ACV of the "hat" spectrum of
%   bin j (1 at f_j, 0 at the neighbouring bins), integrated over the two-sided
%   band [-1/2, 1/2]:
%
%       K(tau, j) = 2 c_j cos(2 pi f_j tau) (1 - cos(2 pi tau d)) / (d (2 pi tau)^2)
%       K(0, j)   = c_j d,       d = 1 / (2(L-1)),   c_j = 2 inside, 1 at both ends
%
%   Why this basis (and not a DCT): every column is the ACV of a non-negative
%   spectrum, so toeplitz(K(:, j)) is positive semi-definite for every size, and
%   so is any combination with S >= 0 -- the Helmert factor of a bin,
%   (K' Omega)_j / (K' s)_j, is then >= 0 by construction and the law never
%   needs a floor (review R7, R7-02: the inverse-DCT components are
%   indefinite). There is no lag period: K is exact at every lag of any arc.
%   S = ones(L, 1) is white noise of unit variance (r = [1; 0; 0; ...]).
%   The same law handed to LSQVCE.ENGINE.COVMODEL as a one-sided 'psd' on
%   fNorm = f_j has ordinates 2 S, and LSQVCE.KERNEL.PSD2ACV reproduces K * S.
%
%   THE PROJECTION FORM. The joint estimator needs only K' * X (X the m x p lag
%   sums of an arc). cos(2 pi f_j tau) has period 2(L-1) in tau for every bin,
%   so X .* env is folded modulo 2(L-1) and one FFT gives all bins: O(m + L log
%   L) time and O(m + L) memory, against m x L for K itself (37 MB for a
%   2160-bin law on a 2160-epoch arc, once per distinct arc length).
%
%   INPUT
%     L   1x1 integer >= 2   number of bins (the Grid)
%     m   1x1 integer >= 1   number of lags (the arc length in samples)
%     X   m x p double       optional: lag sums to project (lag 0 first)
%   OUTPUT
%     K   m x L double       the bin ACVs, lag 0 first; with X: L x p, K' * X
%   ERRORS
%     lsqvce:hatACV:badX     X does not have m rows
%     (the arguments block refuses L < 2 and m < 1 by MATLAB's validators)
%   ACCURACY
%     Closed form; against quadrature of the hat spectra 9e-8 (the quadrature's
%     own error, mirrors/wp111/positive_components.py). The projection form
%     agrees with K' * X to rounding (the FFT's, ~1e-15 of sum |X|).
%   EXAMPLE
%     K = lsqvce.kernel.hatACV(60, 256);
%     r = K * ones(60, 1);          % white: r(1) = 1, r(2:end) ~ 0
%     q = lsqvce.kernel.hatACV(60, 256, randn(256, 1));   % = K' * x
%
%   See also LSQVCE.KERNEL.PSD2ACV, LSQVCE.ENGINE.PSDVCE, LSQVCE.ENGINE.ADJUSTVCEJOINT
%
%   Created 2026-10-11 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- R7 WP-111 part 3 (R7-02).
arguments
    L (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(L, 2)}
    m (1,1) double {mustBeInteger, mustBePositive}
    X double = []
end
    d   = 1 / (2 * (L - 1));
    w   = 2 * pi * (0:m-1).';
    env = (1 - cos(w * d)) ./ (d * w.^2);
    env(1) = d / 2;
    c = 2 * ones(1, L);  c([1 L]) = 1;
    if nargin < 3
        K = 2 * (env .* cos(w * ((0:L-1) * d))) .* c;
        return
    end
    if size(X, 1) ~= m
        error('lsqvce:hatACV:badX', 'X has %d rows but m is %d (one row per lag).', size(X, 1), m);
    end
    P = 2 * (L - 1);                       % the common period of every bin's cosine
    Y = env .* X;
    F = zeros(P, size(X, 2));
    for k0 = 1:P:m                         % fold the lags modulo P
        rows = k0:min(k0 + P - 1, m);
        F(1:numel(rows), :) = F(1:numel(rows), :) + Y(rows, :);
    end
    G = real(fft(F, [], 1));
    K = 2 * c.' .* G(1:L, :);
end
