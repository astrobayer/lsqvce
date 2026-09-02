function [fac, s, w] = robustWeights(vw, influence, tuning)
%ROBUSTWEIGHTS  IRLS weight factors from standardised whitened residuals.
%   [fac, s, w] = lsqvce.kernel.robustWeights(vw, influence, tuning)
%
%   THE one influence-function evaluation of the toolbox (R5 M2 round 2b):
%   the flat driver lsqvce.engine.robust and the block-path driver
%   lsqvce.engine.adjustRobust both call it, so Huber / IGG3 / bisquare and
%   the MAD scale are defined once. Pure kernel: arrays in, arrays out.
%
%   The robust scale is the normalised MAD of the whitened residuals,
%   s = 1.4826 * median(|vw - median(vw)|) (falls back to std(vw), then 1,
%   when the MAD is zero or not finite); the standardised residual is
%   w = vw / s; the factor multiplies the observation's weight:
%     huber     fac = 1 for |w| <= c, c/|w| beyond                 (c = 1.345)
%     igg3      1 on [0, k0]; (k0/|w|) ((k1-|w|)/(k1-k0))^2 on (k0, k1];
%               0 beyond                                            (k0, k1 = 1.5, 3.0)
%     bisquare  (1 - (w/c)^2)^2 for |w| < c, 0 beyond              (c = 4.685)
%
%   INPUT
%     vw         n x 1 double   whitened (law-standardised) residuals
%     influence  char/string    'huber' | 'igg3' | 'bisquare'
%     tuning     [] | 1xK double  tuning constants; [] = the defaults above
%   OUTPUT
%     fac  n x 1 double  weight factors in [0, 1]
%     s    1x1 double    the robust scale used
%     w    n x 1 double  the standardised residuals vw / s
%
%   ERRORS
%     lsqvce:robustWeights:badInfluence  unknown influence name
%     lsqvce:robustWeights:badTuning     tuning has the wrong length for the
%                                        influence (huber/bisquare 1, igg3 2)
%
%   ACCURACY  exact arithmetic on the inputs; no tolerance of its own.
%
%   EXAMPLE
%     fac = lsqvce.kernel.robustWeights(w, 'huber', []);   % 1 below 1.345 sigma
%
%   See also LSQVCE.ENGINE.ADJUSTROBUST, LSQVCE.ENGINE.ROBUST
%
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M2, WP-81 round 2b
%   (extracted verbatim from lsqvce.engine.robust so there is ONE definition).

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
    end

    med = median(vw);
    s   = 1.4826 * median(abs(vw - med));
    if s <= 0 || ~isfinite(s)
        s = std(vw);
        if s <= 0 || ~isfinite(s), s = 1; end
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
