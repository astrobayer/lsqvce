function r = psd2acv(f, S, dt, n)
%LSQVCE.KERNEL.PSD2ACV  One-sided PSD -> autocovariance sequence (Wiener-Khinchin).
%
%   r = LSQVCE.KERNEL.PSD2ACV(f, S, dt, n) returns the autocovariance r(0:n-1) at lag
%   spacing dt, obtained from the one-sided power spectral density S(f) by
%
%       r(tau) = int_0^{fNy} S(f) * cos(2*pi*f*tau) df
%
%   with S the PIECEWISE-LINEAR spectrum through the given ordinates, integrated
%   EXACTLY: per segment [a, b] (centre c, half-width h, x = 2 pi tau h)
%       (S_a + S_b) h cos(2 pi tau c) sin(x)/x - (S_b - S_a) h sin(2 pi tau c) (sin x - x cos x)/x^2
%   (the second factor by its series below x = 1e-2, where it cancels).
%   R7 WP-111 part 2 (finding R7-01): this was trapezoidal quadrature of
%   S cos on the grid, which is periodic in the lag with period 2(nf-1)/(2 fNy dt)
%   samples -- on the uniform grid covModel resampled to, r(m-1) came back as
%   r(1) (even m) or r(0) (odd m), a circulant covariance (measured: Toeplitz
%   cond 3.3e8 for a phi = 0.5 law). The exact integral has no period: r is
%   the ACV of the spectrum S itself, at every lag, and a non-negative S gives a
%   positive semi-definite Toeplitz(r) of every size.
%
%   INPUT
%     f    nf x 1  frequency grid [Hz], strictly increasing, from 0 (WP-153:
%                  the band [0, f(1)] is not in the integral -- an AR(0.99) law
%                  on 512 bins from its first bin lost 34 % of its variance. A
%                  grid whose omitted band, estimated flat as f(1) S(1), exceeds
%                  1e-6 of r(0) is refused; a negligible one -- a habitual start
%                  at 1e-9 Hz -- passes unchanged)
%     S    nf x 1  one-sided PSD [unit^2/Hz], S >= 0
%     dt   scalar  lag spacing [s], > 0 (the sampling interval of the data)
%     n    scalar  number of lags to return, >= 1
%
%   OUTPUT
%     r    n x 1   autocovariance, r(1) = r(0) = variance of the process
%
%   CONVENTION
%     One-sided: var(x) = trapz(f, S) = r(0). If your PSD is two-sided, halve
%     the frequencies and double S, or simply pass S_two_sided*2 over f >= 0.
%     A factor-of-two error here shows up as sigma0 ~ sqrt(2) or 1/sqrt(2),
%     not as an obvious failure -- check r(1) against var(residuals).
%
%   ACCURACY
%     Exact for the piecewise-linear S (1e-11 against a 2e6-point quadrature of
%     the same interpolant, at any lag); r(1) == trapz(f, S). The lags are
%     formed in blocks of at most 4e6 lag-segment products (WP-153): memory
%     O(n + nf), where one n x nf product was 4 GB at 1e5 lags of 512 bins. What remains is
%     how well the linear interpolant represents YOUR spectrum: a grid too
%     coarse near a peak biases r(0) low -- resolve the narrowest feature of S.
%
%   BEST PRACTICE
%     Do not extend f beyond the Nyquist frequency 1/(2*dt) of the data: power
%     above Nyquist is aliased in the observations but not in this integral,
%     which would make r inconsistent with the sampled process.
%
%   EXAMPLE
%     % PSD of an AR(1) with phi = 0.9, innovation sigma = 1, dt = 5 s
%     dt = 5;  a = [1 -0.9];
%     f  = linspace(0, 0.5/dt, 4096).';
%     H  = 1 ./ (exp(-1i*2*pi*f*dt*(0:1)) * a(:));
%     S  = 1^2 * dt * 2 * abs(H).^2;              % one-sided
%     r  = lsqvce.kernel.psd2acv(f, S, dt, 5);
%     r(1)                                        % ~ 1/(1-0.9^2) = 5.263
%
%   ERRORS
%     lsqvce:psd2acv:tooFewPoints         the PSD has too few ordinates to invert
%     lsqvce:psd2acv:noZeroFrequency      the grid starts above 0 and the band below
%                                   it is not negligible: f(1) S(1) > 1e-6 r(0) (WP-153)
%     lsqvce:psd2acv:nonPositiveVariance  the implied variance is not positive,
%                                   so the PSD is not a valid spectrum
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:psd2acv:aboveNyquist                [warning] The PSD extends to X Hz, above the Nyquist frequency X Hz implied by...
%
%   See also LSQVCE.KERNEL.LEVINSONDURBIN, LSQVCE.KERNEL.ARMAACV, LSQVCE.ENGINE.CHECKPSD
%
%   Updated 2026-10-18 (CEST, UTC+02:00), WP-153 (review 12.84, C1, C5): a grid
%   that starts above 0 with a band below it that is not negligible is refused
%   (noZeroFrequency); the lags are
%   formed in blocks, memory O(n + nf). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-111 part 2 (R7-01): the exact
%   integral of the piecewise-linear spectrum replaces trapezoidal quadrature
%   (no lag period, no circulant). With assistance from Claude (Anthropic),
%   model Claude Opus 5.5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+kernel (R4 decision R4-D1,
%   WP-64a): the KERNEL tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    narginchk(4,4);
    f = f(:).';   S = S(:).';

    validateattributes(f,  {'double','single'}, ...
        {'row','finite','nonnegative','increasing'}, mfilename, 'f',  1);
    validateattributes(S,  {'double','single'}, ...
        {'row','finite','nonnegative','numel',numel(f)}, mfilename, 'S', 2);
    validateattributes(dt, {'double','single'}, ...
        {'scalar','finite','positive'}, mfilename, 'dt', 3);
    validateattributes(n,  {'numeric'}, ...
        {'scalar','integer','positive'}, mfilename, 'n',  4);

    if numel(f) < 2
        error('lsqvce:psd2acv:tooFewPoints', ...
              'The frequency grid needs at least 2 points (got %d).', numel(f));
    end

    fNy = 1/(2*dt);
    if f(end) > fNy*(1 + 1e-6)
        lsqvce.util.log.warn('lsqvce:psd2acv:aboveNyquist', ...
            ['The PSD extends to %.4g Hz, above the Nyquist frequency %.4g Hz ' ...
             'implied by dt = %g s. Power above Nyquist is aliased in the data ' ...
             'but not in this integral, so r will not match the sampled ' ...
             'process. Truncate the PSD at Nyquist.'], f(end), fNy, dt);
    end

    h  = (f(2:end) - f(1:end-1)) / 2;           % 1 x nseg half-widths
    c  = (f(2:end) + f(1:end-1)) / 2;           % 1 x nseg centres
    Sa = S(1:end-1);  Sb = S(2:end);
    % WP-153 (C5): the lags in blocks, at most 4e6 products at a time -- one
    % n x nseg block was 4 GB at 1e5 lags of 512 bins, allocated before any
    % guard on the arc's length could refuse it
    r  = zeros(double(n), 1);
    nb = max(1, floor(4e6 / numel(h)));
    for i0 = 1:nb:double(n)
        ii = (i0:min(double(n), i0 + nb - 1)).';
        w  = 2*pi * (ii - 1) * dt;              % block x 1 angular lags
        x  = w .* h;                            % block x nseg
        sx = ones(size(x));
        nz = x ~= 0;  sx(nz) = sin(x(nz)) ./ x(nz);
        g  = x / 3 - x.^3 / 30;                 % (sin x - x cos x)/x^2, small x
        bg = abs(x) >= 1e-2;
        g(bg) = (sin(x(bg)) - x(bg) .* cos(x(bg))) ./ x(bg).^2;
        r(ii) = sum((Sa + Sb) .* h .* cos(w .* c) .* sx - (Sb - Sa) .* h .* sin(w .* c) .* g, 2);
    end

    if r(1) <= 0
        error('lsqvce:psd2acv:nonPositiveVariance', ...
              'r(0) = %g is not positive. The PSD integrates to zero.', r(1));
    end
    % WP-153 (C1): the band [0, f(1)] is not in the integral; estimated flat it
    % carries f(1) S(1) -- refused when that is not negligible against r(0)
    if f(1) > 0 && f(1) * S(1) > 1e-6 * r(1)
        error('lsqvce:psd2acv:noZeroFrequency', ...
              ['The frequency grid starts at %g Hz, and the band [0, %g] Hz below it -- ' ...
               'about %.3g of the variance, estimated flat -- is not in the autocovariance ' ...
               '(an AR(0.99) law on 512 bins from its first bin lost 34 %% this way). Start ' ...
               'the grid at f = 0 with the spectrum''s value there: [0; f] and [S(1); S] ' ...
               'extend it flat.'], f(1), f(1), f(1) * S(1) / r(1));
    end
end
