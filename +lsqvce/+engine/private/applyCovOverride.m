function cm = applyCovOverride(cm, ovr, obsType, opt, day, arcID)
%APPLYCOVOVERRIDE  Merge a runtime covariance override into a base cov spec.
%   cm = applyCovOverride(cm, ovr, obsType, opt)
%   cm = applyCovOverride(cm, ovr, obsType, opt, day, arcID)
%
%   The ONE place the toolbox resolves opt.covOverride. Both loader kinds
%   call it -- LSQVCE.ENGINE.MAKELOADER batches and the raw-closure template
%   LSQVCE.ENGINE.SYNTHLOADER -- so the two shapes LSQVCE.ENGINE.ADJUSTVCEJOINT and
%   LSQVCE.ENGINE.JOINTMODELOPT emit are honoured identically everywhere (R4
%   finding G-06: the template used to read only the flat shape while the
%   VCE engines emit only the per-satellite shape, so a template-derived
%   loader silently ignored every variance component).
%
%   Two shapes, resolved in this order (later wins):
%     (1) ovr.(obsType)                 flat, applies to every satellite
%     (2) ovr.bySys.(satKey).(obsType)  per-satellite, keyed on opt.sysID via
%                                       matlab.lang.makeValidName
%   Both are optional; an absent layer leaves cm untouched. A third layer
%   multiplies the result instead of replacing it:
%     (3) ovr.byArc   containers.Map  LSQVCE.UTIL.ROWWEIGHTKEY(sysID, day,
%                                       obsType, arcID) -> a RELATIVE standard
%                                       deviation of ONE arc (R7-13: the joint
%                                       estimator's per-arc factor, divided by
%                                       its stream's level, which layer (2)
%                                       already carries). Read only when the
%                                       caller passes day and arcID -- a data
%                                       arc; a prior pseudo-arc has no factor.
%
%   INPUT
%     cm       1x1  struct  base covariance spec (see LSQVCE.ENGINE.COVMODEL)
%     ovr      1x1  struct  opt.covOverride as emitted by the engines
%     obsType  1xN  char    observation type this batch belongs to
%     opt      1x1  struct  the loader UserOpt; .sysID selects the bySys key
%     day      1x1          the arc's day     (optional, with arcID)
%     arcID    1x1          the arc's id      (optional, with day)
%
%   OUTPUT
%     cm       1x1  struct  the merged spec, fields of ovr winning over cm
%
%   ERRORS
%     none of its own; LSQVCE.UTIL.MERGESTRUCT errors on non-struct layers
%
%   See also LSQVCE.ENGINE.MAKELOADER, LSQVCE.ENGINE.SYNTHLOADER, LSQVCE.ENGINE.ADJUSTVCEJOINT,
%            LSQVCE.ENGINE.JOINTMODELOPT, LSQVCE.UTIL.MERGESTRUCT, LSQVCE.UTIL.ROWWEIGHTKEY
%
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-116 part 2 (R7-13): layer (3),
%   the per-arc factor, so the joint estimator's arc factors weight the
%   solve. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Created 2026-08-31 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5 -- R4 G-06, Level L1. Moved verbatim
%   from the makeLoader local function so there is ONE implementation.
    if isfield(ovr, obsType)
        cm = lsqvce.util.mergeStruct(cm, ovr.(obsType));
    end
    if isfield(ovr, 'bySys') && isfield(opt, 'sysID')
        satKey = matlab.lang.makeValidName(char(string(opt.sysID)));
        if isfield(ovr.bySys, satKey) && isfield(ovr.bySys.(satKey), obsType)
            cm = lsqvce.util.mergeStruct(cm, ovr.bySys.(satKey).(obsType));
        end
    end
    % R5 M3 3b: an override may carry a RELATIVE level, sigmaScale, instead of
    % an absolute sigma -- the pooled joint law is a SHAPE, and every stream
    % keeps its own declared level times the factor estimated for it. Applied
    % here, once, so both loader kinds agree.
    % R7-13: the per-arc factor scales the stream's level (layer 3).
    if nargin >= 6 && isfield(ovr, 'byArc') && ~isempty(ovr.byArc) && isfield(opt, 'sysID')
        key = lsqvce.util.rowWeightKey(opt.sysID, day, obsType, arcID);
        if isKey(ovr.byArc, key)
            cm.sigmaScale = lsqvce.util.getfielddef(cm, 'sigmaScale', 1) * ovr.byArc(key);
        end
    end
    if isfield(cm, 'sigmaScale') && ~isempty(cm.sigmaScale)
        cm.sigma = lsqvce.util.getfielddef(cm, 'sigma', 1) * cm.sigmaScale;
        cm = rmfield(cm, 'sigmaScale');
    end
end
