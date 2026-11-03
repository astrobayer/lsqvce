function lay = layout(task)
%LSQVCE.LAYOUT  The column layout fit will use -- so you can declare priors in it.
%   lay = lsqvce.layout(task)
%
%   The toolbox derives from your task (scopes, calendar, arc schedules,
%   validity) which column of A and of every system's B belongs to which
%   instantiated parameter. This function announces that layout BEFORE fit
%   runs, so that a prior operator L (task.addPrior(..., L=...) in the columns
%   of A, sys.addPrior(..., L=...) in the columns of B) can be written in it --
%   the only knowledge your regularisation and the toolbox share. It is
%   task.layout(), the method; both read the same plans fit builds.
%
%   INPUT   task  1x1 lsqvce.TaskDefinition   a streaming task (the API
%                                              reference, chapter The task
%                                              definition)
%   OUTPUT  lay   1x1 struct
%     .nc            1x1        columns of A
%     .common        nc x 1 struct (.name .row .day)   the name, the template
%                               row (its position in addCommon order) and (for
%                               "day" commons) the day of every column of A
%     .systems(k)    1x1 struct per system, in addSystem order:
%        .name       string     the system's name
%        .nl         1x1        columns of B
%        .local      nl x 1 struct (.name .row .day .arc)   the same for every
%                               column of B
%        .col        @(row, day, arc) -> column     lookup: the column of
%                               template row `row` on `day` (and `arc` for
%                               "arc" locals; [] otherwise); 0 when that
%                               instance does not exist
%     .colCommon     @(row, day) -> column          the same for A ([] day
%                               for "common" templates)
%
%   EXAMPLE  (the mean of a "day" common over the calendar is zero)
%     lay = lsqvce.layout(task);
%     c   = find([lay.common.name] == "offset");
%     L   = zeros(1, lay.nc);  L(c) = 1 / numel(c);
%     task.addPrior("offsetMean", L=L, Sigma=1e-3);
%
%   ACCURACY
%     Not applicable: it counts and indexes columns -- integers from the
%     plans fit builds (LSQVCE.ENGINE.LAYOUTDECLARATION on the compiled
%     declaration); no floating-point result. test_taskDefinition pins the
%     names and columns.
%
%   ERRORS
%     lsqvce:layout:structDeclaration   a struct was passed: since 7.0.0 the
%                                       declaration is a lsqvce.TaskDefinition;
%                                       the struct's layout is the engine tier's,
%                                       LSQVCE.ENGINE.LAYOUTDECLARATION
%     MATLAB:validators:mustBeA         neither a struct nor a lsqvce.TaskDefinition
%     lsqvce:TaskDefinition:wrongKind   an in-memory task (fromMatrices,
%                                       fromBatches) has no layout to announce
%     lsqvce:TaskDefinition:incomplete  the task does not compile yet
%
%   See also LSQVCE.TASKDEFINITION, LSQVCE.FIT, LSQVCE.RESULTS,
%            LSQVCE.ENGINE.LAYOUTDECLARATION
%
%   Created 2026-11-03 (CET, UTC+01:00) with assistance from Claude (Anthropic),
%   model Claude Opus 5.5 -- level v7, WP-162: the user tier of layout takes a
%   lsqvce.TaskDefinition only; the layout of a struct is
%   LSQVCE.ENGINE.LAYOUTDECLARATION, where the body of layout (2026-09-03,
%   R5-D8 round 5, to 2026-11-02, WP-160) now lives.
arguments
    task (1,1) {mustBeA(task, ["struct", "lsqvce.TaskDefinition"])}
end
    if isstruct(task)
        error('lsqvce:layout:structDeclaration', ...
              ['lsqvce.layout takes a lsqvce.TaskDefinition since 7.0.0, not a data struct ' ...
               '(RELEASE_NOTES.md, "Migrating a declaration from v6 to v7"); the layout of a ' ...
               'compiled struct is lsqvce.engine.layoutDeclaration(data).']);
    end
    lay = task.layout();
end
