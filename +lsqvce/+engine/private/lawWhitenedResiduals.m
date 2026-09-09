function [vw, keys, sizes] = lawWhitenedResiduals(res, W)
%LAWWHITENEDRESIDUALS  Law-whitened residuals in sweep order, the IRLS factors
%   undone (R6 WP-105: the ONE definition adjustRobust and adjustVCEjoint share).
%   [vw, keys, sizes] = lawWhitenedResiduals(res, W)
%
%   The factors of a robust pass ride INSIDE each arc's whitener
%   (lsqvce.engine.makeLoader, opt.rowWeights), so the residual a sweep
%   returns is w_sweep = sqrt(factor) .* w_law. The influence function wants
%   the LAW-whitened residual -- the standardised one under the current
%   stochastic model -- so the factor of the previous pass is divided out
%   here, per arc, by the key it was written under.
%
%   INPUT
%     res   1x1 struct   lsqvce.engine.residuals(..., 'Whitened', true,
%                        'Collect', true): .sys(k).sysID, .sys(k).obs(j)
%                        with .w (m x 1), .day, .obsType, .arc
%     W     containers.Map   'sysID|day|obsType|arcId' -> m x 1 factors of the
%                        pass the residuals were formed under (may be empty)
%   OUTPUT
%     vw    n x 1 double   law-whitened residuals, arcs in sweep order
%     keys  A x 4 cell     {sysID, obsType, day, arc} per arc, in that order
%     sizes 1 x A double   rows per arc (sum = n), to put factors back
%
%   ERRORS  none of its own.
%
%   See also LSQVCE.ENGINE.ADJUSTROBUST, LSQVCE.ENGINE.ADJUSTVCEJOINT,
%            LSQVCE.KERNEL.ROBUSTWEIGHTS, LSQVCE.UTIL.ROWWEIGHTKEY
%
%   Created 2026-09-09 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R6 WP-105 (lifted verbatim from
%   adjustRobust>collectWhitened, where it lived since R5 M2 2b).
    vw = [];  keys = cell(0, 4);  sizes = [];
    for k = 1:numel(res.sys)
        S = res.sys(k);
        for j = 1:numel(S.obs)
            o = S.obs(j);
            w = o.w(:);
            key = lsqvce.util.rowWeightKey(S.sysID, o.day, o.obsType, o.arc);
            if ~isempty(W) && isKey(W, key)
                w = w ./ sqrt(max(W(key), 1e-12));
            end
            vw = [vw; w]; %#ok<AGROW>
            keys(end+1, :) = {S.sysID, o.obsType, o.day, o.arc}; %#ok<AGROW>
            sizes(end+1) = numel(w); %#ok<AGROW>
        end
    end
end
