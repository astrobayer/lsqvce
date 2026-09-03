function lay = layout(data)
%LSQVCE.LAYOUT  The column layout fit will use -- so you can declare priors in it.
%   lay = lsqvce.layout(data)
%
%   The toolbox derives from your parameter tables (scope tokens, calendar,
%   arcs, validity) which column of A and of every system's B belongs to
%   which instantiated parameter. This function announces that layout BEFORE
%   fit runs, so that data.priors(j).L (columns of A) and
%   data.systems(k).priors(l).L (columns of B) can be written in it -- the
%   only knowledge the user's regularisation and the toolbox share (R5-D8).
%
%   INPUT   data  1x1 struct  a kind-A declaration (DATA_CONTRACT.md)
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
%     lay = lsqvce.layout(data);
%     nl  = lay.systems(1).nl;  rows = zeros(0, nl);
%     for d = data.days
%         c = lay.systems(1).col(7, d, []);  r = zeros(1, nl);  r(c) = 1;  rows(end+1, :) = r;
%     end
%     data.systems(1).priors = struct('L', rows, 'x0', 0, 'sigma', 0.01, 'name', 'xaccscale', 'estimate', "scale");
%
%   See also LSQVCE.FIT, LSQVCE.RESULTS/TOTABLE, DATA_CONTRACT.md
%
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5-D8 round 5.
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
        plan = systems{k}.plan;
        loc = lsqvce.engine.localColumns(plan);
        S(k).nl = plan.nl;  S(k).local = loc;
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
