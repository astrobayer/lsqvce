function r = psd2acv(f, S, dt, n)
%LSQVCE.KERNEL.PSD2ACV  One-sided PSD -> autocovariance sequence (Wiener-Khinchin).
%
%   r = LSQVCE.KERNEL.PSD2ACV(f, S, dt, n) returns the autocovariance r(0:n-1) at lag
%   spacing dt, obtained from the one-sided power spectral density S(f) by
%
%       r(tau) = int_0^{fNy} S(f) * cos(2*pi*f*tau) df
%
%   evaluated with trapezoidal quadrature on the supplied frequency grid.
%
%   INPUT
%     f    nf x 1  frequency grid [Hz], strictly increasing, f >= 0
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
%     The quadrature error is set by the frequency grid, not by n. A grid that
%     is too coarse near a spectral peak biases r(0) low. Use >= 2048 points
%     and make sure the grid resolves the narrowest feature of S. Verify with
%       assert(abs(r(1) - trapz(f,S))/r(1) < 1e-10)   % identity at tau = 0
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
%     lsqvce:psd2acv:nonPositiveVariance  the implied variance is not positive,
%                                   so the PSD is not a valid spectrum
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:psd2acv:aboveNyquist                [warning] The PSD extends to X Hz, above the Nyquist frequency X Hz implied by...
%
%   See also LSQVCE.KERNEL.LEVINSONDURBIN, LSQVCE.KERNEL.ARMAACV, LSQVCE.ENGINE.CHECKPSD
%
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

    tau = (0:double(n)-1).' * dt;               % n x 1
    r   = trapz(f, S .* cos(2*pi*tau*f), 2);    % n x 1

    if r(1) <= 0
        error('lsqvce:psd2acv:nonPositiveVariance', ...
              'r(0) = %g is not positive. The PSD integrates to zero.', r(1));
    end
end
