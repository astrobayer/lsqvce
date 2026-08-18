function v = getfielddef(s, f, d)
%LSQVCE.UTIL.GETFIELDDEF  Field value with a default.
%
%   v = LSQVCE.UTIL.GETFIELDDEF(s, f, d) returns s.(f) if s is a struct that has
%   a nonempty field f, and d otherwise.
%
%   Treats an empty field as absent, so 'BatchOverlap',[] means "use the
%   default" rather than "use zero overlap" -- a distinction that matters a lot
%   here, since zero overlap silently destroys accuracy.
%
%   INPUT
%     s  1x1 struct | any  the struct to read; a non-struct yields the default
%     f  1xk  char         field name
%     d  any                the default, returned when s.(f) is absent or empty
%
%   OUTPUT
%     v  any                s.(f) if present and nonempty, otherwise d
%
%   EXAMPLE
%     p = lsqvce.util.getfielddef(opt, 'arOrder', 30);
%
%   See also LSQVCE.UTIL.MERGESTRUCT, ISFIELD
%
%   Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 4.8.

    if isstruct(s) && isscalar(s) && isfield(s, f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = d;
    end
end
