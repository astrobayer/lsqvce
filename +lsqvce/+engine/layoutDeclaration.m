function lay = layoutDeclaration(data)
%LSQVCE.ENGINE.LAYOUTDECLARATION  The column layout of a compiled declaration (engine tier).
%   lay = lsqvce.engine.layoutDeclaration(data)
%
%   ENGINE TIER (v7.0.0, WP-162): the layout of the compiled struct; a user
%   asks the task, task.layout() / lsqvce.layout(task), which calls this
%   function on task.declaration() and adds the names.
%
%   The toolbox derives from your parameter tables (scope tokens, calendar,
%   arcs, validity) which column of A and of every system's B belongs to
%   which instantiated parameter. This function announces that layout BEFORE
%   fit runs, so that data.priors(j).L (columns of A) and
%   data.systems(k).priors(l).L (columns of B) can be written in it -- the
%   only knowledge the user's regularisation and the toolbox share.
%
%   INPUT   data  1x1 struct  a kind-A compiled declaration (the API reference,
%                             chapter The compiled declaration)
%   OUTPUT  lay   1x1 struct
%     .nc            1x1        columns of A
%     .common        nc x 1 struct (.row .day)      the template row and (for 'day' rows) the day of every column of A
%     .systems(k)    1x1 struct per system:
%        .nl         1x1        columns of B
%        .local      nl x 1 struct (.row .day .arc) the template row / day / arc of every column of B
%        .col        @(row, day, arc) -> column     lookup: the column of template row `row` on
%                                                   `day` (and `arc` for 'arc' rows; [] otherwise);
%                                                   0 when that instance does not exist
%     .colCommon     @(row, day) -> column          the same for A ([] day for 'common' rows)
%
%   EXAMPLE  (a Tikhonov prior on local row 7 of system 1, on every day)
%     lay = lsqvce.engine.layoutDeclaration(data);
%     nl  = lay.systems(1).nl;  rows = zeros(0, nl);
%     for d = data.days
%         c = lay.systems(1).col(7, d, []);  r = zeros(1, nl);  r(c) = 1;  rows(end+1, :) = r;
%     end
%     data.systems(1).priors = struct('L', rows, 'x0', 0, 'sigma', 0.01, 'name', 'xaccscale', 'estimate', "scale");
%
%   ACCURACY
%     Not applicable: it counts and indexes columns -- integers from the
%     declaration's plans, the very plans fit builds (declaration); no
%     floating-point result. test_fit pins the columns
%     (testLayoutOfASystemWithoutLocals,
%     testSystemWithoutDaysClaimsOnlyTheGlobalCalendar). WP-139.
%
%   See also LSQVCE.LAYOUT, LSQVCE.ENGINE.FITDECLARATION, LSQVCE.ENGINE.DECLARATION,
%            the API reference
%
%   Updated 2026-11-03 (CET, UTC+01:00), WP-162: the struct layout moved here
%   from lsqvce.layout, which takes a lsqvce.TaskDefinition only. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-11-02 (CET, UTC+01:00), WP-160: a lsqvce.TaskDefinition
%   is laid out with names. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, "Not applicable" with its reason. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (review 3 finding 13): a
%   system without localParams gets nl = 0 and an empty column table (it
%   stopped on the loader's empty plan). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5-D8 round 5.
%
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-124: development history taken
%   out of the help text (it lives in RELEASE_NOTES.md and dev/HANDOVER.md);
%   references point at the API reference. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
    [spec, systems] = lsqvce.engine.declaration(data);
    cp = spec.commonPlan;
    common = struct('row', cell(cp.nc, 1), 'day', cell(cp.nc, 1));
    for j = 1:cp.commonCount, common(j).row = cp.rows.common(j);  common(j).day = []; end
    for i = 1:numel(cp.days)
        if cp.dayCount(i) == 0, continue, end
        rows = cp.rows.day;  if isfield(cp, 'dayRowsOf'), rows = cp.dayRowsOf{i}; end
        for q = 1:cp.dayCount(i)
            j = cp.dayOffset(i) + q;  common(j).row = rows(q);  common(j).day = cp.days(i);
        end
    end
    K = numel(systems);  S = struct('nl', cell(1, K), 'local', cell(1, K), 'col', cell(1, K));
    for k = 1:K
        loc = lsqvce.engine.localColumns(systems{k}.plan);   % [] plan: nl = 0 (WP-135)
        S(k).nl = systems{k}.nl;  S(k).local = loc;
        S(k).col = @(row, day, arc) findCol(loc, row, day, arc);
    end
    lay = struct('nc', cp.nc, 'common', common, 'systems', S, ...
                 'colCommon', @(row, day) findCol(common, row, day, []));
end

function c = findCol(tbl, row, day, arc)
    c = 0;
    for j = 1:numel(tbl)
        e = tbl(j);
        if e.row ~= row, continue, end
        if ~isempty(e.day) && (isempty(day) || e.day ~= day), continue, end
        if isfield(e, 'arc') && ~isempty(e.arc) && (isempty(arc) || e.arc ~= arc), continue, end
        c = j;  return
    end
end
