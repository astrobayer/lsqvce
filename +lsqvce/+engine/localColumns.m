function loc = localColumns(plan)
%LSQVCE.ENGINE.LOCALCOLUMNS  Which (template row, day, arc) owns each column of B.
%   loc = lsqvce.engine.localColumns(plan)
%
%   The one place the scoped local plan is unrolled into a per-column table:
%   lsqvce.layout, Results.toTable and the prior-group expansion (a 'group'
%   pattern on data.systems(k).priors) all read it, so a column means the
%   same instance everywhere.
%
%   INPUT   plan  1x1 struct   the loader's scoped plan (makeLoader L.plan)
%   OUTPUT  loc   nl x 1 struct  .row (template row), .day ([] for 'system'
%                 rows), .arc ([] unless an 'arc' row) -- in the block order
%                 of B (per day: arc blocks then the day block; then the
%                 system block)
%
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5-D8 round 6.
    loc = struct('row', cell(plan.nl, 1), 'day', cell(plan.nl, 1), 'arc', cell(plan.nl, 1));
    for i = 1:numel(plan.days)
        d = plan.days(i);
        arcRows = plan.rows.arc;  dayRows = plan.rows.day;
        if isfield(plan, 'arcRowsOf'), arcRows = plan.arcRowsOf{i};  dayRows = plan.dayRowsOf{i}; end
        for a = 1:numel(plan.arcId{i})
            for q = 1:numel(arcRows)
                j = plan.offset{i}(a) + q;  loc(j).row = arcRows(q);  loc(j).day = d;  loc(j).arc = plan.arcId{i}(a);
            end
        end
        for q = 1:plan.dayCount(i)
            j = plan.dayOffset(i) + q;  loc(j).row = dayRows(q);  loc(j).day = d;  loc(j).arc = [];
        end
    end
    for q = 1:plan.sysCount
        j = plan.sysOffset + q;  loc(j).row = plan.rows.system(q);  loc(j).day = [];  loc(j).arc = [];
    end
end
