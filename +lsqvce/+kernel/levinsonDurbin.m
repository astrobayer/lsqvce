function [a, s2, refl] = levinsonDurbin(r, p)
%LSQVCE.KERNEL.LEVINSONDURBIN  Autocovariance -> AR(p) coefficients (Levinson recursion).
%
%   [a, s2] = LSQVCE.KERNEL.LEVINSONDURBIN(r, p) solves the Yule-Walker equations in
%   O(p^2) and returns the AR(p) model whose autocovariance matches r(0:p).
%
%   [a, s2, refl] = LSQVCE.KERNEL.LEVINSONDURBIN(...) also returns the reflection
%   coefficients (partial autocorrelations).
%
%   INPUT
%     r    (>=p+1) x 1  autocovariance r(0..p); r(1) = r(0) = variance > 0
%     p    scalar       AR order, 1 <= p <= numel(r)-1
%
%   OUTPUT
%     a     1 x (p+1)   AR polynomial [1 a1 ... ap]
%     s2    scalar      innovation variance
%     refl  1 x p       reflection coefficients, |refl| < 1 iff r is PD
%
%   WHY THIS MATTERS
%     The AR model turns the dense decorrelation W = chol(C)^-1, which costs
%     O(m^3) and needs the m x m covariance in memory, into an FIR filter of
%     memory p costing O(m*p) and needing no matrix at all. See LSQVCE.KERNEL.ARMAWHITEN.
%
%   CHOOSING p
%     Raise p until the refitted autocovariance matches the input. Use
%     LSQVCE.ENGINE.CHECKPSD, which does exactly this and reports the relative error.
%     Rules of thumb for satellite residuals at 5 s sampling:
%       white / near-white ......... p = 0..2
%       smooth coloured floor ...... p = 20..40
%       floor + orbital harmonics .. p = 60..100, or an ARMA (LSQVCE.KERNEL.AR2ARMA)
%     p scales roughly with the correlation length measured in samples, so
%     halving dt roughly doubles the p you need.
%
%   NUMERICAL BEHAVIOUR
%     The recursion is exact for a positive definite r. If a reflection
%     coefficient reaches |kappa| >= 1 the input autocovariance is not PD --
%     usually because the PSD it came from has zeros, negative values, or a
%     too-coarse frequency grid. The function then truncates at the last valid
%     order and warns rather than returning a non-stationary model.
%
%   EXAMPLE
%     % exact recovery of a known AR(3)
%     a0 = [1 -1.6 0.9 -0.25];   s0 = 0.7;
%     r  = s0^2 * lsqvce.kernel.armaACV(a0, 1, 8);
%     [a, s2] = lsqvce.kernel.levinsonDurbin(r, 3);
%     assert(norm(a - a0) < 1e-10 && abs(sqrt(s2) - s0) < 1e-10)
%
%   ERRORS
%     lsqvce:levinson:notPD          the autocovariance is not positive definite,
%                                    so it is not a valid ACV sequence
%     lsqvce:levinson:orderTooHigh   the requested order exceeds the ACV length
%
%   See also LSQVCE.KERNEL.ARMAACV, LSQVCE.KERNEL.ARMAWHITEN, LSQVCE.ENGINE.CHECKPSD, LSQVCE.KERNEL.PSD2ACV
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+kernel (R4 decision R4-D1,
%   WP-64a): the KERNEL tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    narginchk(2,2);
    r = r(:).';

    validateattributes(r, {'double'}, {'row','finite','nonempty'}, mfilename, 'r', 1);
    validateattributes(p, {'numeric'}, {'scalar','integer','nonnegative'}, mfilename, 'p', 2);

    p = double(p);
    if p > numel(r)-1
        error('lsqvce:levinson:orderTooHigh', ...
            ['AR order p = %d needs %d autocovariance lags, but only %d were ' ...
             'supplied. Pass r(0..p).'], p, p+1, numel(r));
    end
    if r(1) <= 0
        error('lsqvce:levinson:notPD', ...
              'r(0) = %g must be positive (it is the process variance).', r(1));
    end

    a    = 1;
    s2   = r(1);
    refl = zeros(1, p);

    for k = 1:p
        kappa = -( r(k+1) + a(2:end)*r(k:-1:2).' ) / s2;

        if ~isfinite(kappa) || abs(kappa) >= 1
            lsqvce.util.log.warn('lsqvce:levinson:notPD', ...
                ['The autocovariance lost positive definiteness at order %d ' ...
                 '(reflection coefficient %.4g, |kappa| >= 1). Truncating to ' ...
                 'AR(%d). The PSD it came from probably has zeros, negative ' ...
                 'values, or too coarse a frequency grid.'], k, kappa, k-1);
            refl = refl(1:k-1);
            a    = a(:).';
            return
        end

        a       = [a, 0] + kappa*[0, fliplr(a)];   %#ok<AGROW>
        s2      = s2 * (1 - kappa^2);
        refl(k) = kappa;
    end

    a = a(:).';
end
