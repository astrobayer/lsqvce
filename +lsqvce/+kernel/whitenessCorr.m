function maxCorr = whitenessCorr(series, maxLag)
%LSQVCE.KERNEL.WHITENESSCORR  Largest absolute autocorrelation over a lag range.
%   maxCorr = lsqvce.kernel.whitenessCorr(series, maxLag)
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
%   ACCURACY
%     Plain lag products over the lag-0 energy, no approximation of its own;
%     the statistic is a sampling quantity. White noise of n = 4000 stays
%     below 3/sqrt(n), an AR(1) with phi = 0.8 lies above it, both inside
%     [0, 1] (test_streamUtils (testMaxLagCorrWhiteVsCorrelated)), as in the
%     Python mirror v9verify.py (HANDOVER 6.1.4, Python mirrors row). WP-139.
%
%   EXAMPLE
%     w = randn(4000,1);
%     assert(lsqvce.kernel.whitenessCorr(w, 5) < 3/sqrt(4000));   % white -> below band
%
%   See also LSQVCE.ENGINE.COLLECTRESIDUALS, LSQVCE.KERNEL.CROSSCORR
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, every number traced to the test that pins it. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34): the
%   deleted SELECTARACV reference removed. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Created 2026-07-15 13:40 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+kernel (R4 decision R4-D1,
%   WP-64a): the KERNEL tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
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
