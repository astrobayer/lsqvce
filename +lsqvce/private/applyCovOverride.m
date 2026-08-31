function cm = applyCovOverride(cm, ovr, obsType, opt)
%APPLYCOVOVERRIDE  Merge a runtime covariance override into a base cov spec.
%   cm = applyCovOverride(cm, ovr, obsType, opt)
%
%   The ONE place the toolbox resolves opt.covOverride. Both loader kinds
%   call it -- LSQVCE.MAKELOADER batches and the raw-closure template
%   LSQVCE.SYNTHLOADER -- so the two shapes LSQVCE.ADJUSTVCE and
%   LSQVCE.JOINTMODELOPT emit are honoured identically everywhere (R4
%   finding G-06: the template used to read only the flat shape while the
%   VCE engines emit only the per-satellite shape, so a template-derived
%   loader silently ignored every variance component).
%
%   Two shapes, resolved in this order (later wins):
%     (1) ovr.(obsType)                 flat, applies to every satellite
%     (2) ovr.bySys.(satKey).(obsType)  per-satellite, keyed on opt.sysID via
%                                       matlab.lang.makeValidName
%   Both are optional; an absent layer leaves cm untouched.
%
%   INPUT
%     cm       1x1  struct  base covariance spec (see LSQVCE.COVMODEL)
%     ovr      1x1  struct  opt.covOverride as emitted by the engines
%     obsType  1xN  char    observation type this batch belongs to
%     opt      1x1  struct  the loader UserOpt; .sysID selects the bySys key
%
%   OUTPUT
%     cm       1x1  struct  the merged spec, fields of ovr winning over cm
%
%   ERRORS
%     none of its own; LSQVCE.UTIL.MERGESTRUCT errors on non-struct layers
%
%   See also LSQVCE.MAKELOADER, LSQVCE.SYNTHLOADER, LSQVCE.ADJUSTVCE,
%            LSQVCE.JOINTMODELOPT, LSQVCE.UTIL.MERGESTRUCT
%
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
end
