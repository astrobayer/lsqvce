function [L, task] = lsqvce_template(procOpt)
%LSQVCE_TEMPLATE  Starting point for your own least-squares adjustment case.
%
%   L = LSQVCE_TEMPLATE(procOpt)
%   [L, task] = LSQVCE_TEMPLATE(procOpt)
%
%   Copy this file, rename it, and fill in the one stub at the bottom. It
%   declares your case as a LSQVCE.TASKDEFINITION, driven by your procOpt,
%   and returns the cell array of Loaders that the engines consume, compiled
%   from that task:
%
%       [L, task] = myCase(procOpt);
%       lsqvce.util.log.set('fid', procOpt.loggingOverviewFileID, ...
%               'screen', procOpt.loggingLevel > 0);    % one sink for the text
%       R = lsqvce.engine.adjustBlock(L, procOpt.days, ...
%               'Verbose', procOpt.loggingLevel > 0);   % one GLS solve on the Loaders
%
%   ENGINE-TIER TEMPLATE (since 2.0). The user path needs no Loader at all:
%   the task this file declares IS the user-tier declaration, and
%   lsqvce.fit(task, lsqvce.Options(Estimate="noise")) runs it end to end
%   (examples/demo_tutorial.m shows a task from scratch). Keep the Loaders
%   when you want them in hand (procOpt-driven cases, direct engine work):
%   they are compiled from the task by the same declaration step fit uses,
%   [~, L] = lsqvce.engine.declaration(task.declaration()).
%
%   THIS FILE HAS TWO REGISTERS -- READ THIS BEFORE RUNNING IT
%     * declareTask (below) is REAL, RUNNABLE code. It declares one system per
%       entry of procOpt.observationGroups, its local parameters and one
%       stream per observation type. You edit it only if your topology differs
%       from one-system-per-group (see MERGING SYSTEMS below).
%     * myBuild is a STUB. It ERRORS until you write it. That is deliberate:
%       only you know how to read your data, and a plausible fake would be
%       worse than an honest stop.
%   So: running this as shipped gets you Loaders (and a task) that error on
%   the first day read. That is the expected state until the stub is filled.
%
%   WHAT YOU MUST WRITE                                       (the stub)
%     myBuild      one stream's day: design rows + observations   REQUIRED
%                  + their error information (Sigma= of each arc: a sigma per
%                  epoch, or a d x d covariance per epoch for multi-axis streams;
%                  the toolbox composes it with the stream's temporal law)
%
%   WHAT YOU MUST DECLARE                                     (in your procOpt)
%     .observationGroups   one entry per system: its streams and identity
%     .days                the day list
%     .estimation          your common parametrisation
%
%   WHAT THE TOOLBOX DERIVES -- DO NOT STATE IT TWICE
%     nc  from the common parameters   |   nl  from the local parameters
%   The counts are stated nowhere: each parameter is declared ONCE, by name,
%   and Build reads the counts as ctx.NumCommon / ctx.NumLocal and the columns
%   of a named parameter as ctx.cols(names) -- never an offset by hand.
%
%   THE PARAMETRISATION TABLES. procOpt.estimation.commonParams (Rc x C) and
%   g.localParams (Rl x C, per group) hold one parameter TEMPLATE per row, in
%   the row convention of the scoped tables: column 1 the scope token
%   ('common' | 'day' for commons, 'arc' | 'day' | 'system' for locals),
%   column 2 the NAME -- unique, as every parameter of a task needs one (WP-161:
%   before, only the row counts were read) -- and columns 3..C yours, which
%   travel as the parameter's Info
%   (variables col3 .. colC; ctx.Common.Info / ctx.Local.Info in Build). Each
%   row becomes one addCommon / addLocal call, in table order (the columns of
%   A follow the commons' order). An 'arc' local needs an arc schedule: give
%   the system Arcs=... in declareTask (none is declared as shipped -- the
%   arcs are the ones myBuild adds).
%
%   INPUT
%     procOpt  1x1  struct  your processing configuration (see PROCOPT below).
%                           Each system's Settings: myBuild reads it as
%                           ctx.Settings, and every Loader keeps it as
%                           L{k}.procOpt. Only LSQVCE.FIT surfaces it to an
%                           OnPass hook (as its context's processingOptions);
%                           at engine tier hand it over yourself --
%                           'OnPassContext', struct('processingOptions',
%                           {{procOpt}}) -- and ADJUSTBLOCK has no OnPass at all.
%
%   OUTPUT
%     L        1xS  cell    one Loader per system, ready for LSQVCE.ENGINE.ADJUSTBLOCK
%                           or LSQVCE.ENGINE.ADJUSTROBUST.
%     task     1x1  lsqvce.TaskDefinition   the declaration the Loaders are
%                           compiled from -- R = lsqvce.fit(task, opts)
%
%   PROCOPT -- IT IS IMMUTABLE INPUT, NOT SHARED STATE
%     procOpt is captured when the task is declared and piped through
%     unchanged for the whole run. So it must NOT carry live cursors: a field
%     like 'iteridx' would freeze at its value on the first pass while the
%     loop iterates around it, and your builder would silently read the wrong
%     pass forever. The engine tracks these itself -- an OnPass hook is handed
%     the live pass number as ctx.pass and the driver's name as ctx.driver,
%     and the system identity rides on the Loader as L{k}.sysID. It must not
%     have a field 'prepared' either: that name is the toolbox's.
%
%     Group your settings into nested structs; only the fields your builder
%     needs have to be there:
%
%       procOpt.days            1xD   day list             (the task's Days)
%       procOpt.loggingLevel    1x1   0 silent | 1 progress | 2 detail
%                                     NOTE 3 is REJECTED (merged into 2), so a
%                                     pipeline that can emit 3 must clamp first.
%       procOpt.loggingOverviewFileID 1x1  fid
%                                     NOT an engine argument: route it to the
%                                     toolbox-wide logger before the solve,
%                                     lsqvce.util.log.set('fid', ...) (user
%                                     tier: lsqvce.util.log's file sink).
%       procOpt.showAndSaveGraphics   1x1  gates hook attachment (an IF in your
%                                     driver, not a toolbox option)
%       procOpt.observationGroups 1xS cell of structs      (drives the loop below)
%                                     .observationTypes  the stream names
%                                     .systemId          numeric or ONE character
%                                     .name              optional, unique label
%                                     .localParams       Rl x C cell (may be {})
%       procOpt.estimation      1x1   struct  .commonParams, Rc x C cell
%       procOpt.data            1x1   struct  paths, maximumDegree, ...
%       procOpt.epochs          1x1   struct  mEpoch, dt, arcsPerDay, ...
%       procOpt.forces .constants .ephemeris .earthRotation .metadata
%                                     your own sub-structs, opaque to the toolbox
%
%     THE FIVE SHARED VALUES. days, loggingLevel, loggingOverviewFileID,
%     showAndSaveGraphics and observationGroups are needed by your pipeline
%     BEFORE any toolbox object exists, so they live here -- and the toolbox is
%     FED from them (see the call at the top of this help). Never type the same
%     value into both: that is exactly the drift this interface was built to
%     remove.
%
%   ACCURACY
%     Not applicable: a declaration and a stub, no arithmetic; the Loaders are
%     the ones lsqvce.fit compiles from the same task.
%
%   ERRORS
%     lsqvce_template:notImplemented   a stub was called before you wrote it
%     lsqvce_template:badProcOpt       procOpt is missing a required field
%     lsqvce:TaskDefinition:duplicateName   two parameters (or two systems) share a
%                                      name; the other declaration mistakes are
%                                      lsqvce:TaskDefinition:* too, raised by the
%                                      add-call (LSQVCE.TASKDEFINITION)
%     lsqvce:DayContext:*              a mistake of myBuild while a day is read
%                                      (a name ctx.cols does not know, ...)
%
%   EXAMPLE
%     procOpt = myPipelineConfig();
%     [L, task] = lsqvce_template(procOpt);
%     R = lsqvce.engine.adjustBlock(L, procOpt.days);   % the engine tier
%     R = lsqvce.fit(task);                              % or the user tier
%     R.report();
%
%   See also LSQVCE.TASKDEFINITION, LSQVCE.SYSTEMDEFINITION, LSQVCE.DAYCONTEXT,
%            LSQVCE.DAYDATA, LSQVCE.NOISE, LSQVCE.ENGINE.DECLARATION,
%            LSQVCE.ENGINE.ADJUSTBLOCK, LSQVCE.FIT, DEMO_TUTORIAL
%
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: the case is declared as a
%   lsqvce.TaskDefinition from procOpt (one addCommon / addLocal per table row,
%   one system per group, one stream per type); the stub is myBuild(day, ctx);
%   the Loaders come from lsqvce.engine.declaration(task.declaration()), and the
%   task is the second output. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-102 (v6): myBuildDay returns
%   one cell per arc (the cells are the arcs; arcIdx is gone). With
%   assistance from Claude (Anthropic), model Claude Fable 5.1.
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-95: one stub, not two; the
%   overview fid goes to lsqvce.util.log.set (there is no 'logFid'); procOpt
%   reaches a hook only through 'OnPassContext' (ctx.processingOptions is a
%   LSQVCE.FIT service) and there is no ctx.sysID; the parametrisation tables
%   are Rc x C / Rl x C, the column count your own convention.
%   With assistance from Claude (Anthropic), model Claude Opus 5.
%   Updated 2026-09-06 (CEST, UTC+02:00), R6 WP-90: the error information
%   comes back from myBuildDay with the day (D.Cep) -- the myEpochCov stub,
%   the epochCovariance hook and the retired plot-hook / adjustVCE text are
%   gone. With assistance from Claude (Anthropic), model Claude Fable 5.1.
%   Created 2026-07-19 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
%
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-124: references to the data
%   contract point at the API reference. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.

    narginchk(1, 1);
    requireFields(procOpt, {'observationGroups', 'days'});
    task = declareTask(procOpt);
    [~, L] = lsqvce.engine.declaration(task.declaration());   % fit's own declaration step
end

%% ======================================================= RUNNABLE: the task
function task = declareTask(procOpt)
%DECLARETASK  The task: commons, one system per observation group. REAL CODE -- edit to taste.
%
%   Each entry of procOpt.observationGroups describes ONE system: which
%   observation types it carries, how it is identified and its local
%   parameters. This loop turns each into a system of the task.
%
%   MERGING SYSTEMS -- the one topology decision you cannot avoid
%     Two satellites must share ONE system whenever an observation's design
%     couples both of their LOCAL parameters (an inter-satellite range or
%     range-rate does exactly this). The block-arrow model has one local block
%     per system, so a coupling observation cannot span two of them. In that
%     case give both satellites' streams the SAME group entry and stack their
%     locals: [y_A ; y_B]. Streams that are merely correlated do NOT need
%     merging -- VCE handles that through the variance components.
%
%   INPUT
%     procOpt  1x1 struct  the processing configuration (see the file header)
%   OUTPUT
%     task     1x1 lsqvce.TaskDefinition
    task = lsqvce.TaskDefinition(Days=procOpt.days);

    % --- the common parameters: one addCommon per row, in table order ------
    C = procOpt.estimation.commonParams;
    for r = 1:size(C, 1)
        task.addCommon(string(C{r, 2}), Scope=string(C{r, 1}), Info=rowInfo(C(r, :)));
    end

    groups = procOpt.observationGroups;
    for k = 1:numel(groups)
        g = groups{k};

        % --- the system --------------------------------------------------------
        % Id MUST be numeric or a SINGLE char: the accumulator indexes
        % systems by [store.sysID], which flattens a multi-char string into
        % separate elements. Carry a readable label in the name instead (it
        % must be unique among the systems). procOpt is the system's Settings:
        % myBuild reads its settings FROM ctx.Settings -- not from a variable
        % captured in this workspace. Two copies of one value is the bug this
        % interface exists to prevent. No Arcs: the arcs are the ones myBuild
        % adds; declare a schedule (Arcs=...) only for 'arc'-scoped locals
        % (nl must be known before any read).
        sys = task.addSystem(groupName(g, k), Settings=procOpt, Id=g.systemId);

        % --- its local parameters: one addLocal per row, in table order --------
        Lp = {};
        if isfield(g, 'localParams'), Lp = g.localParams; end
        for r = 1:size(Lp, 1)
            sys.addLocal(string(Lp{r, 2}), Scope=string(Lp{r, 1}), Info=rowInfo(Lp(r, :)));
        end

        % --- what this group observes: one stream per observation type --------
        % obsTypes are the STREAM names. They are the unit of variance-component
        % estimation, so name them after what they physically are, not after the
        % file they came from.
        obsTypes = cellstr(g.observationTypes);
        for j = 1:numel(obsTypes)
            sys.addStream(obsTypes{j}, Build=@myBuild, Noise=aprioriNoise(obsTypes{j}, g));
        end
    end
end

%% ============================================================== STUB: your data
function out = myBuild(day, ctx)  %#ok<STOUT> -- a stub: it raises until you write it
%MYBUILD  Read ONE day of ONE observation type. YOU WRITE THIS.
%
%   It must be DETERMINISTIC. The accumulator reads every day a SECOND time for
%   the residual pass and must get byte-identical data back; a random draw or a
%   changing file would corrupt the residuals silently, not loudly.
%
%   INPUT
%     day      1x1   double               the day identifier from procOpt.days
%     ctx      1x1   lsqvce.DayContext    what the toolbox knows about this call:
%                    ctx.Stream     which stream to build (1x1 string)
%                    ctx.Settings   your procOpt, piped unchanged
%                    ctx.NumCommon  nc, the columns of A (addCommon order;
%                                   a 'day' common at ctx.cols(name))
%                    ctx.NumLocal   nl, the columns of B of this system
%                    ctx.Common / ctx.Local   the parameter tables (Name,
%                                   Scope, OnDay, Info = your columns 3..C)
%                    ctx.cols(names[, Arc=id])   the columns of named parameters
%
%   OUTPUT
%     out  1x1 lsqvce.DayData -- out = ctx.newDay(), then ONE addArc PER ARC of
%          the day (without an Arcs schedule the arcs are the ones you add, in
%          order; a day with one arc adds one):
%            out.addArc(A, l, t, Local={names, P}, Sigma=s)
%       A        m_a x nc    common design rows of the arc
%       l        m_a x 1     observed minus computed
%       t        m_a x 1     epoch timestamps [s]
%       Local    {names, P}  the local partials BY NAME, P m_a x numel(names)
%                            (several rows {names1, P1; names2, P2}); or B=B, an
%                            m_a x nl matrix (may be sparse) whose columns come
%                            from ctx.cols(names) -- never an offset by hand
%       Sigma    THE ERROR INFORMATION of the arc: nEpoch_a x 1 (the per-epoch
%                standard deviation of a scalar stream), nEpoch_a x d (per-axis
%                standard deviations) or nEpoch_a x d x d (the full per-epoch
%                covariance, e.g. from a PPP solution). The toolbox composes it
%                with the temporal law: C = D (R_t kron I_d) D. The axis count
%                is derived from it. (optional; omit when every observation is
%                equally accurate)
%       Excluded m_a x 1 logical: rows of weight 0, kept for their index
%     A multi-axis stream WITHOUT Sigma states its axes per epoch once, on the
%     stream: addStream(..., Axes=d).
%
%   ROW ORDER FOR MULTI-AXIS STREAMS: rows are EPOCH-MAJOR, i.e.
%   [x1 y1 z1 | x2 y2 z2 | ...]. The temporal correlation model strides by
%   the axis count, so a stream stored axis-major comes out with its time law
%   stretched by the axis count -- and it will not announce itself.
    error('lsqvce_template:notImplemented', ...
        ['myBuild is a stub -- write it for your data.\n' ...
         'Asked for day %s, stream ''%s''. Return ctx.newDay() with one addArc per ' ...
         'arc (Sigma= for the error information); see `help lsqvce_template`.'], ...
         num2str(day), char(ctx.Stream));
end

%% ========================================================== small local helpers
function n = aprioriNoise(obsType, ~)  %#ok<INUSD>
%APRIORINOISE  The a-priori temporal law of one stream, a LSQVCE.NOISE.
%   White and unit-level by default: the VCE loop discovers the truth, so a
%   deliberately plain start makes the discovery visible. Replace with your
%   best prior if you have one (lsqvce.Noise.arma, .psd, ...). The per-epoch
%   accuracy of the observations is not declared here: myBuild returns it with
%   every arc (Sigma=), for scalar and multi-axis streams alike.
%
%   INPUT
%     obsType   1xk char     the stream name
%     g         1x1 struct   the observation-group entry (unused here)
%   OUTPUT
%     n         1x1 lsqvce.Noise
    n = lsqvce.Noise.white(1);
end

function T = rowInfo(row)
%ROWINFO  Columns 3..C of one parametrisation row as its Info (variables col3 ..).
%   INPUT   row  1xC cell   {scope, name, yours...}
%   OUTPUT  T    1x(C-2) table, one cell variable per column (table() for C <= 2)
    T = table();
    for c = 3:numel(row)
        T.(sprintf('col%d', c)) = row(c);
    end
end

function s = groupName(g, k)
%GROUPNAME  A human-readable label for a group (the system's unique name).
%   INPUT   g 1x1 struct; k 1x1 double  the group index
%   OUTPUT  s 1x1 string
    if isfield(g, 'name') && ~isempty(g.name)
        s = string(g.name);
    else
        s = string(sprintf('system-%d', k));
    end
end

function requireFields(s, names)
%REQUIREFIELDS  Error clearly when procOpt is missing something required.
%   Fails at the TOP with the field name, rather than deep inside the
%   declaration with a message about an argument the user never typed.
%
%   INPUT   s 1x1 struct; names 1xN cellstr
%   OUTPUT  none (throws lsqvce_template:badProcOpt)
    missing = names(~isfield(s, names));
    if ~isempty(missing)
        error('lsqvce_template:badProcOpt', ...
            'procOpt is missing required field(s): %s.', strjoin(missing, ', '));
    end
end
