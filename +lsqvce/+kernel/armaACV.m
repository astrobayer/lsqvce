function r = armaACV(a, b, n)
%LSQVCE.KERNEL.ARMAACV  Autocovariance of an ARMA(p,q) with unit innovation variance.
%
%   r = LSQVCE.KERNEL.ARMAACV(a, b, n) returns r(0:n-1) for the stationary process
%
%       A(z) x(k) = B(z) e(k),      var(e) = 1
%
%   with a = [1 a1 ... ap] and b = [1 b1 ... bq]. For a pure AR pass b = 1.
%   Scale by the innovation variance to get the physical autocovariance:
%   r_phys = sigma^2 * lsqvce.kernel.armaACV(a, b, n).
%
%   INPUT
%     a    1 x (p+1)  AR polynomial, a(1) must be 1 (normalised)
%     b    1 x (q+1)  MA polynomial, b(1) must be 1 (normalised)
%     n    scalar     number of lags, >= 1
%
%   OUTPUT
%     r    n x 1      autocovariance for unit innovation variance
%
%   METHOD
%     Computed from the Wold (impulse-response) expansion via FFT rather than
%     from the recursive Yule-Walker equations. This is deliberate: the
%     recursion loses accuracy for near-unit-root AR polynomials, which is
%     exactly the regime of a steep low-frequency noise floor. The impulse
%     response grows from 512 samples by doubling until the energy of its
%     last max(p,q)+1 samples is below eps^2 of the whole, capped at 2^20
%     samples; its length follows the decay alone, so r is the same law for
%     any n, and one FFT of nh + n - 1 points (no wrap) gives every lag.
%     COST: O((nh + n) log(nh + n)) -- linear-logarithmic in n.
%
%   ACCURACY
%     r(0) to a few eps relative while the response decays within 2^20
%     samples (|root| below about 0.99998); beyond, the truncation is
%     rho^(2^21) of r(0), rho the largest root modulus, and it is reported
%     (lsqvce:armaACV:truncated). Lags beyond the decay length are below
%     eps * r(0) absolute (the FFT's roundoff floor).
%
%   BEST PRACTICE
%     The process must be stationary: all roots of a(z) strictly inside the
%     unit circle. A near-unit root (|root| -> 1) means an extremely long
%     correlation length, which makes the AR whitening filter memory p large
%     and the batch overlap expensive. If lsqvce.engine.checkPSD reports a marginally
%     stable fit, the PSD probably has an unmodelled trend -- detrend the
%     residuals before estimating it.
%
%   EXAMPLE
%     % analytic check against ARMA(1,1): phi = 0.9, theta = 0.5
%     r  = lsqvce.kernel.armaACV([1 -0.9], [1 0.5], 3);
%     r0 = (1 + 2*0.9*0.5 + 0.5^2)/(1 - 0.9^2);   % = 9.7368
%     r1 = 0.9*r0 + 0.5;                          % = 9.2632
%     assert(abs(r(1)-r0)/r0 < 1e-8 && abs(r(2)-r1)/r1 < 1e-8)
%
%   ERRORS
%     lsqvce:armaACV:zeroLeadingAR  the AR polynomial does not start with 1
%     lsqvce:armaACV:nonStationary  the AR part has a root on or inside the unit
%                                   circle, so no stationary ACV exists
%   WARNINGS
%     lsqvce:armaACV:truncated      the response had not decayed at 2^20
%                                   samples; r(0) is low by rho^(2^21) -- once
%                                   per run
%
%   See also LSQVCE.KERNEL.LEVINSONDURBIN, LSQVCE.KERNEL.ARMAWHITEN, LSQVCE.KERNEL.PSD2ACV
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-134: the response length follows
%   its decay alone (from 512, doubling to the 2^20 cap for any n -- the old
%   loop stopped at 65 536 for n <= 128), one FFT of nh + n - 1 points, and
%   a truncation at the cap is reported. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-131: the impulse response stops
%   on the energy of its last max(p,q)+1 samples, not on one sample (a sparse
%   AR polynomial is exactly 0 there long before it has decayed). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+kernel (R4 decision R4-D1,
%   WP-64a): the KERNEL tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    narginchk(3,3);
    a = a(:).';   b = b(:).';

    validateattributes(a, {'double'}, {'row','finite','nonempty'}, mfilename, 'a', 1);
    validateattributes(b, {'double'}, {'row','finite','nonempty'}, mfilename, 'b', 2);
    validateattributes(n, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'n', 3);

    if a(1) == 0
        error('lsqvce:armaACV:zeroLeadingAR', 'a(1) must be nonzero.');
    end
    a = a / a(1);   b = b / b(1);        % normalise

    rts = roots(a);
    if ~isempty(rts) && max(abs(rts)) >= 1
        error('lsqvce:armaACV:nonStationary', ...
            ['The AR polynomial is not stationary (max |root| = %.6f >= 1). ' ...
             'The autocovariance does not exist. Detrend the data or reduce ' ...
             'the AR order.'], max(abs(rts)));
    end

    n = double(n);

    % Impulse response length: grow until it has decayed, capped for safety.
    % WP-134: the length follows the decay alone, not the lag count. Beyond
    % the decay length every r(k) is below eps * r(0), and nfft >= nh + n - 1
    % keeps the FFT correlation linear (no wrap) at every lag returned. The
    % old start max(4n, 512) with at most 7 doublings capped the response at
    % 65 536 samples for n <= 128, short of the 2^20 promised (|root| 0.99999:
    % r(0) 27 % low, and n-dependent), and paid a 4n-tap filter and an 8n-point
    % FFT for a response that had decayed after a few hundred samples.
    nh = 512;
    K  = max(numel(a), numel(b));
    while true
        h = filter(b, a, [1; zeros(nh-1,1)]);
        % WP-131: the TAIL's energy over max(p,q)+1 samples, not one sample --
        % a sparse AR polynomial ([1 0 phi], a seasonal [1 0..0 -phi]) has an
        % h that is exactly 0 at most samples long before it has decayed: the
        % one-sample test stopped at 512 and r(0) of [1 0 0.998] came out 36 %
        % low (dev/validation/wp131_kernel_mirror.py)
        tail = sum(h(end-min(nh, K)+1:end).^2);
        if tail < eps^2 * sum(h.^2) || nh >= 2^20, break, end
        nh = 2*nh;
    end
    rho = max([abs(rts); 0]);
    if rho^(2*nh) > eps          % tol: eps -- the relative truncation of r(0)
        % the cap stopped an undecayed response: the dominant pole's tail
        % beyond nh is rho^(2 nh) of r(0) (exact for an AR(1))
        lsqvce.util.log.warnOnce('lsqvce:armaACV:truncated', ...
            ['The impulse response of an AR root of modulus %.7f has not decayed ' ...
             'after %d samples; the autocovariance is truncated there and r(0) is ' ...
             'low by about %.1e relative. Such a root is a near-unit-root law: ' ...
             'detrend the residuals before estimating it.'], rho, nh, rho^(2*nh));
    end

    nfft = 2^nextpow2(nh + n - 1);
    R    = ifft(abs(fft(h, nfft)).^2, 'symmetric');
    r    = R(1:n);
    r    = r(:);
end
