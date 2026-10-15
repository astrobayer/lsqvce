function [fac, s, w, scalePath, kReject] = robustWeights(vw, influence, tuning)
%LSQVCE.KERNEL.ROBUSTWEIGHTS  IRLS weight factors from standardised whitened residuals.
%   [fac, s, w, scalePath, kReject] = lsqvce.kernel.robustWeights(vw, influence, tuning)
%
%   THE one influence-function evaluation of the toolbox (R5 M2 round 2b):
%   the block-path driver lsqvce.engine.adjustRobust and the joint
%   estimator's robust pass both call it, so Huber / IGG3 / bisquare and
%   the MAD scale are defined once. Pure kernel: arrays in, arrays out.
%
%   The robust scale is the normalised MAD of the whitened residuals,
%   s = 1.4826 * median(|vw - median(vw)|). When the MAD is degenerate --
%   zero, not finite, or below 1e3*eps of the 75th percentile of
%   |vw - median(vw)| (half the residuals tied, e.g. quantised data or
%   zero-redundancy rows; the percentile is the bulk's spread, which a
%   blunder in the upper quarter cannot move, where the largest deviation
%   IS the worst blunder) -- it falls back to std(vw), then 1, SAYS SO once
%   per run
%   (lsqvce:robustWeights:scaleFallback) and returns the path taken; the
%   std is inflated by the very blunders the MAD exists to expose, so a
%   fallback result is a weaker robust fit. The standardised residual is
%   w = vw / s; the factor multiplies the observation's weight:
%     huber     fac = 1 for |w| <= c, c/|w| beyond                 (c = 1.345)
%     igg3      1 on [0, k0]; (k0/|w|) ((k1-|w|)/(k1-k0))^2 on (k0, k1];
%               0 beyond                                            (k0, k1 = 1.5, 3.0)
%     bisquare  (1 - (w/c)^2)^2 for |w| < c, 0 beyond              (c = 4.685)
%
%   INPUT
%     vw         n x 1 double   whitened (law-standardised) residuals
%     influence  char/string    'huber' | 'igg3' | 'bisquare'
%     tuning     [] | 1xK double  tuning constants; [] = the defaults above;
%                finite and positive, igg3 with k0 < k1
%   OUTPUT
%     fac  n x 1 double  weight factors in [0, 1]
%     s    1x1 double    the robust scale used
%     w    n x 1 double  the standardised residuals vw / s
%     scalePath  1xK char  'mad' | 'std' | 'unit' -- which scale s is
%     kReject    1x1 double  the rejection point in units of s: the |w|
%                beyond which fac = 0 (igg3 k1, bisquare c; Inf for huber,
%                which rejects nothing) -- what the joint estimator's level
%                correction needs (R7 WP-115 part 3b)
%
%   ERRORS
%     lsqvce:robustWeights:badInfluence  unknown influence name
%     lsqvce:robustWeights:badTuning     tuning has the wrong length for the
%                                        influence (huber/bisquare 1, igg3 2),
%                                        a value not finite and positive, or
%                                        igg3 with k0 >= k1
%   WARNINGS
%     lsqvce:robustWeights:scaleFallback  the MAD was degenerate; s is the
%                                        std (or 1) -- once per run
%
%   ACCURACY  exact arithmetic on the inputs; no tolerance of its own.
%
%   EXAMPLE
%     fac = lsqvce.kernel.robustWeights(w, 'huber', []);   % 1 below 1.345 sigma
%
%   See also LSQVCE.ENGINE.ADJUSTROBUST, LSQVCE.ENGINE.ADJUSTVCEJOINT
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-140: the H1 line names the unit
%   in full, lsqvce.kernel.robustWeights, as HANDOVER 6.2.1 asks (it read the
%   bare name). With assistance from Claude (Anthropic), model Claude Opus
%   5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-134 (third review, finding 9):
%   the near-tie test measured the MAD against the LARGEST deviation, the
%   worst blunder -- one residual beyond ~4.5e12 MAD switched the scale to
%   the std it inflates (Huber's estimate x 2e8, 8-sigma rows factor 1);
%   the reference is the 75th percentile now. Tuning values are validated
%   (a c <= 0 zeroed or negated every factor). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-115 part 3b: the rejection
%   point kReject is returned (one table of tuning constants). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34): names
%   its live callers. With assistance from Claude (Anthropic), model Claude
%   Opus 5.5.
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M2, WP-81 round 2b
%   (extracted verbatim from lsqvce.engine.robust so there is ONE definition).
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-115 (finding R7-21): the
%   MAD -> std -> 1 fallback was silent, and a NEAR-tie MAD (3.5e-15) was
%   accepted, giving every ordinary row a factor of ~1e-14 that the collapse
%   guard could not see. Near-ties now count as degenerate, the fallback
%   warns once and returns its path. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.

    vw = vw(:);
    infl = lower(char(influence));
    switch infl
        case 'huber',    tune = 1.345;
        case 'igg3',     tune = [1.5 3.0];
        case 'bisquare', tune = 4.685;
        otherwise
            error('lsqvce:robustWeights:badInfluence', ...
                  'Unknown influence function ''%s'' (huber | igg3 | bisquare).', infl);
    end
    if ~isempty(tuning)
        if numel(tuning) ~= numel(tune)
            error('lsqvce:robustWeights:badTuning', ...
                  '%s takes %d tuning constant(s); %d given.', infl, numel(tune), numel(tuning));
        end
        tune = tuning(:).';
        if ~all(isfinite(tune) & tune > 0) || (strcmp(infl, 'igg3') && tune(1) >= tune(2))
            error('lsqvce:robustWeights:badTuning', ...
                  ['%s tuning constants must be finite and positive (igg3: k0 < k1); ' ...
                   'got [%s].'], infl, num2str(tune));
        end
    end

    switch infl
        case 'huber',    kReject = Inf;
        case 'igg3',     kReject = tune(2);
        case 'bisquare', kReject = tune(1);
    end
    med = median(vw);
    dev = abs(vw - med);
    s   = 1.4826 * median(dev);
    scalePath = 'mad';
    ds  = sort(dev);
    q75 = 0;  if ~isempty(ds), q75 = ds(ceil(0.75 * numel(ds))); end   % the bulk's spread
    % tol: 1e3*eps -- a MAD this far below the bulk's spread is a tie, not a
    % scale (R7-21); the 75th percentile, not max(dev): one blunder beyond
    % 1/(1e3 eps) ~ 4.5e12 MAD moved the max and switched the scale (WP-134)
    if ~isfinite(s) || s <= 1e3 * eps * q75
        s = std(vw);
        scalePath = 'std';
        if s <= 0 || ~isfinite(s), s = 1;  scalePath = 'unit'; end
        lsqvce.util.log.warnOnce('lsqvce:robustWeights:scaleFallback', ...
            ['The MAD of the whitened residuals is degenerate (half or more of them ' ...
             'tie -- quantised data or zero-redundancy rows); the robust scale falls ' ...
             'back to the %s (s = %.3g). That scale is inflated by the blunders the ' ...
             'MAD exists to expose: the robust weights are weaker than intended.'], ...
            scalePath, s);
    end
    w  = vw / s;
    aw = abs(w);
    n  = numel(w);
    switch infl
        case 'huber'
            fac = ones(n, 1);
            m = aw > tune(1);
            fac(m) = tune(1) ./ aw(m);
        case 'igg3'
            k0 = tune(1);  k1 = tune(2);
            fac = ones(n, 1);
            mid = (aw > k0) & (aw <= k1);
            fac(mid) = (k0 ./ aw(mid)) .* ((k1 - aw(mid)) ./ (k1 - k0)).^2;
            fac(aw > k1) = 0;
        case 'bisquare'
            c = tune(1);
            fac = zeros(n, 1);
            m = aw < c;
            fac(m) = (1 - (aw(m)/c).^2).^2;
    end
end
