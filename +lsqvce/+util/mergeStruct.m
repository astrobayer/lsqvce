function s = mergeStruct(s, t)
%LSQVCE.UTIL.MERGESTRUCT  Overlay the fields of t onto s.
%
%   s = LSQVCE.UTIL.MERGESTRUCT(s, t) returns s with every field of t copied over,
%   replacing existing values. Fields of s not present in t are kept.
%
%   Used to let a runtime option override a static covariance spec:
%     cm = lsqvce.util.mergeStruct(cfg.cov.pos, opt.covOverride.pos);
%
%   INPUT
%     s  1x1 struct | any  the base; a non-struct is REPLACED by t outright
%     t  1x1 struct | []   the overlay; empty or a non-struct leaves s untouched
%
%   OUTPUT
%     s  1x1 struct        the merged struct (t wins on any shared field)
%
%   NOTE this is a SHALLOW merge: a struct-valued field of t replaces the whole
%   corresponding field of s rather than being merged into it recursively.
%
%   ACCURACY
%     Not applicable: field values are copied, not computed (test_util
%     pins later-wins, added fields and the empty overlay). WP-139.
%
%   EXAMPLE
%     a = struct('type','psd','arOrder',30);
%     b = struct('arOrder',60);
%     c = lsqvce.util.mergeStruct(a, b);      % type='psd', arOrder=60
%
%   See also LSQVCE.UTIL.GETFIELDDEF, FIELDNAMES
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, "Not applicable" with its reason. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 4.8.

    if ~isstruct(t) || isempty(t)
        return
    end
    if ~isstruct(s)
        s = t;
        return
    end
    f = fieldnames(t);
    for k = 1:numel(f)
        s.(f{k}) = t.(f{k});
    end
end
