function out = psdVCE(Q, T, Sprior, opt)
%LSQVCE.ENGINE.PSDVCE  Per-bin spectral amplitudes by variance component estimation.
%
%   out = LSQVCE.ENGINE.PSDVCE(Q, T, Sprior)
%   out = LSQVCE.ENGINE.PSDVCE(Q, T, Sprior, 'Name', value, ...)
%
%   The spectral half of the joint estimator's update (R7 WP-111 part 3, user
%   decision 2026-10-11: POSITIVE COMPONENTS). The law is a non-negative
%   spectrum S on L bins, piecewise linear in frequency, with covariance
%   r = K S (LSQVCE.KERNEL.HATACV). Every bin j is a variance component whose
%   matrix V_j = toeplitz(K(:, j)) is the ACV of a hat spectrum -- positive
%   semi-definite -- so the Helmert factor of a bin,
%
%       alpha2_j = sum_m Q(j, m) / sum_m T(j, m),
%       Q(:, m) = lambda_m K' Omega^m,   T(:, m) = lambda_m K' s^m,
%
%   (Omega^m, s^m the arc's lag sums of the doubly whitened residual and of the
%   redundancy matrix over ALL lags, LSQVCE.KERNEL.REDUNDANCYLAGS; lambda_m the
%   arc's scale) is a ratio of two NON-NEGATIVE quantities, and
%
%       Shat_j = alpha2_j * Sprior_j                                (Helmert)
%
%   stays a valid spectrum without a floor or a zero-frequency rule. Ellmer's
%   (2018) per-frequency factors were taken in the type-I DCT domain, whose
%   components (the columns of the inverse DCT) are INDEFINITE: alpha2_j could
%   be negative, needed an absolute floor, and a floored bin made the declared
%   covariance near-singular -- the level ran away x5e7 with a NEGATIVE lag-1
%   correlation at phi 0.8 / 0.9 (review R7, R7-02, measured on the runner,
%   HANDOVER 12.44). Under the same data the positive-component update
%   converged 18 / 18 with the level 0.91-1.07 x the truth.
%
%   INPUT
%     Q       L x M double   per-arc bin numerators, one COLUMN per arc (>= 0)
%     T       L x M double   per-arc bin denominators (>= 0)
%     Sprior  L x 1 double   the spectrum the arcs were weighted with (>= 0)
%
%   OPTIONS (Name/Value)
%     'Gamma'  1x1 in (0,1], default 1.  Damping of the multiplicative update,
%              Shat = alpha2.^Gamma .* Sprior (WP-40 r2): 1 is the plain step;
%              a closed-loop caller may pass a DECAYING Gamma per pass.
%     'Smooth' 1x1 integer >= 0, default 0.  Half-width of a moving average over
%              neighbouring bins (reflected edges), trading resolution for
%              variance. R7-22: Options.Noise.Smooth.
%     'FreezeRedundancy' 1x1 >= 0, default 0.  A bin whose redundancy
%              Sprior_j * sum_m T(j, m) (its share tr(R V_j) of the total) is
%              below it HOLDS its amplitude (alpha2_j = 1) instead of taking a
%              ratio of two small numbers; a bin whose denominator is zero or
%              not finite always holds. R7-22.
%
%   OUTPUT
%     out   1x1 struct
%       .Shat    L x 1  the updated spectrum (>= 0; bins below 1e-10 x its
%                       maximum are lifted to that, the one numerical floor)
%       .alpha2  L x 1  the per-bin factors (>= 0); the model has settled when
%                       they are all ~ 1
%       .frozen  L x 1  logical, the bins that held their value
%       .Sprior  L x 1  the prior, as given
%       .maxDev  1x1    max |alpha2 - 1| over EVERY bin (there is no
%                       zero-frequency rule left to exclude one: R7-20);
%                       with Smooth > 0 the change of the SHAPE, max |Shat./
%                       Sprior - 1| at r(0) = 1 -- a smoothed update's fixed
%                       point has alpha2 ~= 1 (WP-131)
%       .dev     L x 1  the per-bin deviation .maxDev is the maximum of:
%                       abs(alpha2 - 1), or with Smooth > 0 the bin's
%                       change of shape abs(rho - 1); 0 for a frozen bin
%                       (WP-147: the stop rule reads it bin by bin)
%       .redundancy L x 1  the bin's redundancy Sprior_j * sum_m T(j, m), its
%                       share tr(R V_j) -- the number FreezeRedundancy is
%                       compared with; 0 where the denominator is not usable.
%                       A variance component of redundancy r has a relative
%                       standard error of about sqrt(2 / r) (WP-147)
%
%   EXAMPLE
%     K = lsqvce.kernel.hatACV(L, m);         % per arc length
%     Q(:, a) = lam(a) * (K' * r.omega);  T(:, a) = lam(a) * (K' * r.s);
%     p = lsqvce.engine.psdVCE(Q, T, S);  S = p.Shat;
%
%   ACCURACY
%     The update is exact arithmetic, a ratio of bin sums; the spectrum is a
%     statistical estimate. Q == T leaves alpha2 = 1 and Shat = Sprior to
%     1e-10 (test_ellmer (test_psdVCE_identity_when_Q_equals_T)); a smoothed
%     loop settles to maxDev < 1e-8 (test_ellmer
%     (test_psdVCE_smoothed_fixed_point_settles)). The positive-component law
%     converged 18 / 18 with the level 0.91-1.07 x the truth (HANDOVER 12.44:
%     the candidate before the port, dense exact traces). WP-139.
%
%   ERRORS
%     lsqvce:psdVCE:shapeMismatch   Q or T does not have one row per bin of Sprior
%     lsqvce:psdVCE:arcMismatch     Q and T disagree on the number of arcs
%
%   See also LSQVCE.KERNEL.HATACV, LSQVCE.KERNEL.REDUNDANCYLAGS,
%            LSQVCE.ENGINE.ARCVARIANCEFACTORS, LSQVCE.ENGINE.ADJUSTVCEJOINT
%
%   Updated 2026-10-16 (CEST, UTC+02:00), WP-147: .dev and .redundancy per
%   bin, for the joint loop's stop rule. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, every number traced to the test that pins it. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-131: with Smooth > 0 maxDev
%   measures the change of the spectrum's shape, so a smoothed loop can
%   converge. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-111 part 3 (R7-02, R7-20):
%   REWRITTEN for positive components -- the inputs are per-bin statistics
%   (K' Omega, K' s), the prior is the bin spectrum; the DCT ('Norm'), the
%   absolute floor ('MinAmplitude'), 'ZeroFreqRule' and the degenerate-start
%   check ('AllowDegenerate') are gone with the indefinite components that
%   needed them. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-116 part 1: 'FreezeRedundancy'.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-40 r2 (F-15): optional damped
%   update 'Gamma'. With assistance from Claude (Anthropic), model Claude Fable 5.
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
arguments
    Q (:,:) double
    T (:,:) double
    Sprior (:,1) double {mustBeNonnegative, mustBeFinite}
    opt.Gamma (1,1) double {mustBePositive, mustBeLessThanOrEqual(opt.Gamma, 1)} = 1
    opt.Smooth (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    opt.FreezeRedundancy (1,1) double {mustBeNonnegative} = 0
end
    L = numel(Sprior);
    if size(Q, 1) ~= L || size(T, 1) ~= L
        error('lsqvce:psdVCE:shapeMismatch', ...
             ['Q (%d x %d) and T (%d x %d) must have one ROW per bin to match ' ...
              'Sprior (%d bins).'], size(Q,1), size(Q,2), size(T,1), size(T,2), L);
    end
    if size(Q, 2) ~= size(T, 2)
        error('lsqvce:psdVCE:arcMismatch', 'Q has %d arcs but T has %d.', size(Q,2), size(T,2));
    end

    Qj = sum(Q, 2);  Tj = sum(T, 2);
    ok = isfinite(Qj) & isfinite(Tj) & Tj > 0;
    alpha2 = ones(L, 1);
    alpha2(ok) = max(Qj(ok), 0) ./ Tj(ok);           % both >= 0 in exact arithmetic
    red = zeros(L, 1);
    red(ok) = Sprior(ok) .* Tj(ok);                   % the bin's redundancy (WP-147: returned)
    frozen = ~ok | ~(red >= opt.FreezeRedundancy);
    if opt.FreezeRedundancy > 0
        alpha2(frozen) = 1;
    end
    if opt.Gamma == 1
        Shat = alpha2 .* Sprior;
    else
        Shat = (alpha2 .^ opt.Gamma) .* Sprior;
    end
    if opt.Smooth > 0
        Shat = smoothSpectrum(Shat, opt.Smooth);
    end
    Shat = max(Shat, 1e-10 * max(Shat));              % the ONE numerical floor, relative
    dev = abs(alpha2 - 1);
    if opt.Smooth > 0
        % WP-131: at the fixed point of a SMOOTHED update S = smooth(alpha2 S)
        % alpha2 is not 1 wherever the true spectrum is not smooth, so
        % max|alpha2 - 1| never met Tol (0.35 at every pass from 5 on, 30 passes,
        % "not converged"). The smoothed loop has settled when the SHAPE stops
        % moving: the change of S at r(0) = 1 (c the hatACV weights; the level
        % is the joint loop's levelDev).
        c = 2 * ones(L, 1);  c([1 L]) = 1;
        in = Sprior > 0;
        rho = (Shat(in) ./ Sprior(in)) * ((c.' * Sprior) / (c.' * Shat));
        dev = zeros(L, 1);
        dev(in) = abs(rho - 1);
    end
    out = struct('Shat', Shat, 'alpha2', alpha2, 'frozen', frozen, 'Sprior', Sprior, ...
                 'maxDev', max(dev), 'dev', dev, 'redundancy', red);
end

function S = smoothSpectrum(S, h)
%SMOOTHSPECTRUM  Moving average of half-width h, with reflected edges.
%   Reflection keeps the endpoints (DC and Nyquist) from being pulled toward
%   zero, which a zero-padded average would do.
%   INPUT   S  nx1 double   amplitudes;  h 1x1 int  half-width
%   OUTPUT  S  nx1 double   smoothed amplitudes, same length
    n = numel(S);
    if h < 1 || n < 3, return, end
    h = min(h, floor((n-1)/2));
    idx = [h+1:-1:2, 1:n, n-1:-1:n-h];      % reflected padding
    Sp  = S(idx);
    w   = ones(2*h+1,1) / (2*h+1);
    Sc  = conv(Sp, w, 'same');
    S   = Sc(h+1 : h+n);
end
