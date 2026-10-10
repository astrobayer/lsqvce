function x = collectResiduals(res, satId, obsType, field)
%LSQVCE.ENGINE.COLLECTRESIDUALS  Stack one residual field over the arcs of a stream.
%   x = lsqvce.engine.collectResiduals(res, satId, obsType)
%   x = lsqvce.engine.collectResiduals(res, satId, obsType, field)
%
%   Walks the collected-residual structure returned by
%   LSQVCE.ENGINE.BLOCKACCUMULATOR/residuals and vertically concatenates one field of
%   every arc that matches obsType (and satId, unless satId is empty).
%
%   INPUT
%     res      1x1   struct   Residual structure with res.sys(:).sysID and
%                             res.sys(:).obs(:).{obsType, v, w}. Required.
%     satId    1x1   char|string|numeric|[]   Satellite id to select, or [] to
%                             pool across ALL satellites. Required (may be []).
%     obsType  1xk   char|string   Observation type to select, e.g. 'pos'.
%                             Required.
%     field    1xk   char|string   Which per-arc field to stack:
%                             'w' whitened residuals | 'v' raw residuals.
%                             (default 'w')
%
%   OUTPUT
%     x        Mx1   double   Concatenation of the selected field over all
%                             matching arcs; [] if nothing matches.
%
%   EXAMPLE
%     w = lsqvce.engine.collectResiduals(res, 'C', 'pos');       % GRACE-C, whitened
%     vAll = lsqvce.engine.collectResiduals(res, [], 'rr', 'v'); % raw, both satellites
%
%   ERRORS
%     lsqvce:collectResiduals:badRes  the residual field is absent or the wrong
%                                   shape
%
%   See also LSQVCE.ENGINE.BLOCKACCUMULATOR, LSQVCE.KERNEL.WHITENESSCORR
%
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34): the
%   deleted KERNEL.RESIDUALACV reference removed. With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Created 2026-07-15 13:40 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
    narginchk(3, 4);
    if nargin < 4 || isempty(field), field = 'w'; end
    field = char(field);
    if ~isstruct(res) || ~isfield(res, 'sys')
        error('lsqvce:collectResiduals:badRes', ...
              'res must be a struct with a .sys field.');
    end

    x = [];
    for iSys = 1:numel(res.sys)
        if ~isempty(satId) && ~isequal(res.sys(iSys).sysID, satId), continue, end
        arcObs = res.sys(iSys).obs;
        % WP-43 (T-03a): the help text promises badRes for an ABSENT FIELD,
        % but a missing field used to surface as a raw MATLAB:nonExistentField
        % from the dynamic access below. The documented id is the contract;
        % fire it here and name the offending field.
        if ~isempty(arcObs) && ~isfield(arcObs, field)
            error('lsqvce:collectResiduals:badRes', ...
                  'the arcs carry no field ''%s'' (expected ''w'' or ''v'').', ...
                  field);
        end
        for iArc = 1:numel(arcObs)
            if strcmp(arcObs(iArc).obsType, obsType)
                x = [x; arcObs(iArc).(field)(:)];         %#ok<AGROW>
            end
        end
    end
end
