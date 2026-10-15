function [mixed, repeated] = idDayRules(ids, days)
%IDDAYRULES  The two rules every set of systems meets: one id kind, each day once.
%   [mixed, repeated] = idDayRules(ids, days)
%
%   Defined ONCE for both front doors (WP-138): lsqvce.engine.declaration
%   checks a declaration for lsqvce.fit, lsqvce.engine.BlockAccumulator a
%   direct caller's systems; each re-implemented both rules (WP-135). The
%   rules return the first clash and the DOOR raises -- its own identifier,
%   with the context only it knows (data.systems(2).days, system 'C').
%
%   ONE ID KIND (WP-135, review 3 finding 12): the accumulator concatenates
%   the ids ([store.sysID]), and a character next to a number turned the
%   number into a control character (3 -> char(3)) -- robust factors went
%   under keys the loader never reads, reported as applied.
%   EACH DAY ONCE (finding 11): a repeated day was accumulated twice,
%   silently -- x pulled toward it, sd too small, dof inflated by its rows.
%
%   INPUT
%     ids    1xK  cell     the system ids as values (numeric or char), in
%                          declaration order; {} when only days are judged
%     days   1xD  double   one day list; [] when only ids are judged
%   OUTPUT
%     mixed     0x0 | 1x2  double  [kChar kNum]: the first non-numeric and the
%                                  first numeric id when the kinds mix
%     repeated  0x0 | 1x1  double  the first (smallest) day listed twice
%
%   ACCURACY
%     Not applicable: a classification of the ids and an exact comparison
%     of the days.
%
%   ERRORS
%     none: the door raises (lsqvce:fit:mixedIdTypes / duplicateDay,
%     lsqvce:accum:mixedIdTypes / duplicateDay)
%
%   EXAMPLE
%     [m, r] = idDayRules({'C', 3}, [1 2 2 3])   % m = [1 2], r = 2
%
%   See also LSQVCE.ENGINE.DECLARATION, LSQVCE.ENGINE.BLOCKACCUMULATOR
%
%   Created 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135):
%   declaration>oneIdKind / dayListOf's repeat test and
%   BlockAccumulator.checkSpec's two tests, one definition. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
    mixed = [];
    isNum = cellfun(@isnumeric, ids);
    if any(isNum) && ~all(isNum)
        mixed = [find(~isNum, 1), find(isNum, 1)];
    end
    repeated = [];
    if numel(days) > 1
        [u, ~, j] = unique(days(:));
        twice = u(accumarray(j, 1) > 1);
        if ~isempty(twice), repeated = twice(1); end
    end
end
