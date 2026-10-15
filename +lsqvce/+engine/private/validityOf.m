function f = validityOf(tbl, valCols)
%VALIDITYOF  @(row, day) -> is a template row valid on calendar day `day`?
%   f = validityOf(tbl, valCols)
%
%   THE validity rule of the parameter templates, defined ONCE for both
%   tables that carry one: a system's localParams (lsqvce.engine.makeLoader,
%   option ValidityColumns) and the commonParams of a declaration
%   (lsqvce.engine.declaration, data.paramColumns). WP-138: the declaration
%   re-implemented it as mjdOf next to the loader's toMjd -- one rule in two
%   places, each fixed separately in WP-135.
%
%   A validity cell holds [] (always valid), a scalar fractional MJD, or
%   [MJD seconds] (-> MJD + seconds/86400); anything else -- a datetime, a
%   char, a NaN, a longer vector -- is refused, never read as "always".
%   A row is valid on day d when start < d + 1 and end > d (the half-open
%   day overlap, R5-D8 round 3); a start after its end is refused (the row
%   would silently never exist).
%
%   INPUT
%     tbl      NxC  cell     the template table
%     valCols  1x2  double   [cStart cEnd], 1-based columns of tbl; [] = no
%                            validity columns (every row always valid)
%   OUTPUT
%     f        1x1  function_handle   tf = f(row, day), 1x1 logical
%
%   ACCURACY
%     Not applicable: the rule compares the parsed numbers exactly; the one
%     arithmetic step, MJD + seconds/86400, is a single rounding (the cell
%     form's resolution is one ulp of the MJD, ~1e-11 day).
%
%   ERRORS
%     lsqvce:makeLoader:badValidity  valCols is not two integer columns of
%                                    tbl, a cell is not [] / an MJD / [MJD
%                                    seconds], or a row's start is after its
%                                    end. lsqvce.engine.declaration raises
%                                    lsqvce:fit:badValidity for the same rule
%                                    on commonParams (this one is the cause).
%
%   EXAMPLE
%     f = validityOf({[], []; 60000, [60001 43200]}, [1 2]);
%     f(1, 59000)      % true: row 1 is always valid
%     f(2, 60001)      % true: row 2 ends at MJD 60001.5
%
%   See also LSQVCE.ENGINE.MAKELOADER, LSQVCE.ENGINE.DECLARATION
%
%   Created 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135):
%   makeLoader>validityOf and >toMjd moved here verbatim, the row named in
%   every message; declaration>mjdOf retired. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
    if isempty(valCols)
        f = @(row, day) true;
        return
    end
    if numel(valCols) ~= 2 || any(valCols < 1) || any(valCols > size(tbl, 2)) || any(valCols ~= round(valCols))
        error('lsqvce:makeLoader:badValidity', ...
              'ValidityColumns must be [cStart cEnd] within 1..%d (got %s).', size(tbl, 2), mat2str(valCols));
    end
    n = size(tbl, 1);
    t0 = -inf(n, 1);  t1 = inf(n, 1);
    for r = 1:n
        a = toMjd(tbl{r, valCols(1)}, r);  b = toMjd(tbl{r, valCols(2)}, r);
        if ~isempty(a), t0(r) = a; end
        if ~isempty(b), t1(r) = b; end
        if t0(r) > t1(r)
            error('lsqvce:makeLoader:badValidity', 'row %d: validity start %g is after end %g.', r, t0(r), t1(r));
        end
    end
    f = @(row, day) t0(row) < day + 1 && t1(row) > day;
end

function m = toMjd(v, r)
%TOMJD  A validity entry as fractional MJD; [] = always. WP-135: an entry of
%   another type (a datetime, a char) or a NaN was read as [] -- "always
%   valid", silently -- and is refused like a wrong length.
    m = [];
    if isempty(v), return, end
    if isnumeric(v) && ~any(isnan(v(:)))
        if isscalar(v), m = double(v);  return, end
        if numel(v) == 2, m = double(v(1)) + double(v(2)) / 86400;  return, end
    end
    error('lsqvce:makeLoader:badValidity', ...
          'row %d: a validity entry must be [], a fractional MJD, or [MJD seconds].', r);
end
