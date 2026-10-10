function fired = warnAtBreakdown(fac, driver)
%WARNATBREAKDOWN  Say once per run that a robust pass is at its breakdown point.
%   fired = warnAtBreakdown(fac, driver)
%
%   THE breakdown diagnostic of both robust drivers (user decision
%   2026-10-10, R7 WP-115; HANDOVER 12.34 option 2). The collapse refusal
%   (nAlive < nc + 1) is a model-error detector: it is reached by a bodily
%   residual offset or by near-total contamination, never at 50 %. What a
%   user needs to hear near 50 % is that the M-estimate is losing its
%   footing -- so when the influence ZEROES more than 45 % of the
%   observations (epochs, for the joint loop) this warns, once per run, and
%   changes nothing: no factor, no refusal, no result. Huber never zeroes a
%   factor, so under Huber it cannot fire (its weights only shrink).
%
%   INPUT   fac     n x 1 double   the pass's weight factors in [0, 1]
%           driver  1xK char       'adjustRobust' | 'adjustVCEjoint', for the text
%   OUTPUT  fired   1x1 logical    true when the warning was issued now
%   WARNINGS
%     lsqvce:robust:atBreakdown  more than 45 % of the factors are zero
%
%   EXAMPLE
%     warnAtBreakdown(fac, 'adjustRobust');   % after robustWeights, each pass
%
%   See also LSQVCE.ENGINE.ADJUSTROBUST, LSQVCE.ENGINE.ADJUSTVCEJOINT,
%            LSQVCE.KERNEL.ROBUSTWEIGHTS
%
%   Created 2026-10-10 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- R7 WP-115 part 2.
    fired = false;
    if isempty(fac), return, end
    zeroed = mean(fac(:) == 0);
    if zeroed > 0.45
        fired = lsqvce.util.log.warnOnce('lsqvce:robust:atBreakdown', ...
            ['%s: the influence zeroed %.0f %% of the observations -- at the breakdown ' ...
             'point of the M-estimate (50 %%). The estimate is no longer robust to the ' ...
             'rejected part: check whether the "blunders" are a model error (a bias, a ' ...
             'datum defect) before trusting it. Nothing was changed.'], ...
            driver, 100 * zeroed);
    end
end
