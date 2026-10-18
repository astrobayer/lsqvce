function d = columnLagSums(X, nL)
%COLUMNLAGSUMS  d(n+1) = sum_c sum_i X(i+n, c) X(i, c), n = 0..nL-1, by FFT.
%   The lag-n sub-diagonal sum of X X' equals the sum over columns of each
%   column's lag-n autocorrelation, so all lags come from one real FFT per
%   column of length >= 2N-1 (zero-padded, no circular wrap): the power
%   spectrum summed over columns, inverse-transformed. O(C N log N) for C
%   columns instead of O(C N nL). Exact to rounding; R5 M3 3c (WP-83).
%   INPUT  X N x C double;  nL 1x1 lags wanted (<= N)
%   OUTPUT d nL x 1 double
%   Shared by LSQVCE.KERNEL.REDUNDANCYLAGS and LSQVCE.KERNEL.TRUNCATIONDEFICIT
%   (moved here from redundancyLags' local functions, WP-156: one
%   implementation, rule 3).
%
%   Updated 2026-10-18 (CEST, UTC+02:00), WP-156: a private kernel function
%   (it was local to redundancyLags). With assistance from Claude (Anthropic),
%   model Claude Opus 5.5.
    [N, C] = size(X);
    d = zeros(nL, 1);
    if N == 0 || C == 0, return, end
    nfft = 2^nextpow2(2*N - 1);
    F = fft(X, nfft, 1);                          % columns, zero-padded
    P = sum(real(F .* conj(F)), 2);               % summed power spectrum
    r = ifft(P, 'symmetric');                     % lag sums, lag 0 first
    d = r(1:nL);
end
