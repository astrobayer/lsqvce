function maxCorr = whitenessCorr(series, maxLag)
%LSQVCE.WHITENESSCORR  Largest absolute autocorrelation over a lag range.
%   maxCorr = lsqvce.whitenessCorr(series, maxLag)
%
%   A cheap whiteness statistic. The series is mean-removed and normalised by
%   its lag-0 energy, so 0 <= maxCorr <= 1. For white noise the sample
%   autocorrelations sit inside roughly +/- k/sqrt(n); compare maxCorr against a
%   band such as whiteFactor/sqrt(numel(series)) to decide whiteness.
%
%   INPUT
%     series   nx1 | 1xn   double    Time series (a residual stream). Required.
%                                     Fewer than 2 elements returns 0.
%     maxLag   1x1         integer   Highest lag to test, >= 1. Required.
%                                     Internally capped at numel(series)-1.
%
%   OUTPUT
%     maxCorr  1x1         double    max over lags 1..maxLag of
%                                     |sum(x_t x_{t+lag})| / sum(x_t^2), in [0,1].
%
%   EXAMPLE
%     w = randn(4000,1);
%     assert(lsqvce.whitenessCorr(w, 5) < 3/sqrt(4000));   % white -> below band
%
%   See also LSQVCE.COLLECTRESIDUALS, LSQVCE.SELECTARACV, LSQVCE.CROSSCORR
%
%   Created 2026-07-15 13:40 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
    narginchk(2, 2);
    validateattributes(series, {'numeric'}, {'vector'}, mfilename, 'series', 1);
    validateattributes(maxLag, {'numeric'}, ...
        {'scalar','integer','positive'}, mfilename, 'maxLag', 2);

    series  = double(series(:));
    series  = series - mean(series);
    nSamp   = numel(series);
    lag0    = series.' * series;
    maxCorr = 0;
    if nSamp < 2 || lag0 == 0
        return
    end
    for lag = 1:min(maxLag, nSamp - 1)
        maxCorr = max(maxCorr, ...
            abs(series(1:nSamp-lag).' * series(1+lag:nSamp)) / lag0);
    end
end
