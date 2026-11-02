classdef SystemDefinition < handle
%LSQVCE.SYSTEMDEFINITION  One observing system of a task: local parameters, streams, priors.
%
%   sys = task.addSystem(name, Name=Value)        made by LSQVCE.TASKDEFINITION only
%   sys.addLocal(names, Name=Value)
%   sys.addStream(name, Build=@build, Noise=noise, Name=Value)
%   sys.addPrior(name, Name=Value)
%
%   A system is one satellite, one receiver, or a pair whose observation
%   depends on the local parameters of both (a range-rate between two
%   spacecraft is ONE system). Systems share only the common parameters;
%   each has its own locals, streams, calendar and priors.
%
%   WHERE THE DAILY WORK GOES. The system's Prepare (given to addSystem) runs
%   once per day, before the first stream of the day; each stream's Build
%   (given to addStream) runs once per day and returns the stream's arcs:
%       P   = prepare(day, ctx)       % its result reaches every Build as ctx.Prepared
%       out = build(day, ctx)         % out = ctx.newDay(); out.addArc(...)
%   (LSQVCE.DAYCONTEXT, LSQVCE.DAYDATA). With Options.CacheDays = 1 a day is
%   built twice per pass (accumulation, residual sweep); CacheDays = Inf builds
%   it once per run. Both functions must return the same result for the same day.
%
%   METHODS
%     addLocal(sys, names, Name=Value)      local parameters, eliminated arc by arc
%       Scope  "arc" (default) | "day" | "system"   constant over one arc, one day of
%                                           the system, or the whole run; "arc" needs
%                                           an Arcs schedule on the system or stream
%       Valid  [start end] fractional MJD | [] (default: always)   when the row exists
%       Info   table, one row per name      your own columns, readable in ctx.Local
%     addStream(sys, name, Name=Value)      an observation stream (one observation type)
%       Build      @(day, ctx) -> LSQVCE.DAYDATA   required; once per stream and day
%       Noise      LSQVCE.NOISE              required; the temporal law and its level
%       Estimate   "law" | "scale" | "none"  what Options.Estimate="noise" may change: the
%                                           shape and level (default "law"), the level
%                                           only, or nothing
%       ArcFactors logical                   one variance factor per arc (default:
%                                           Options.ArcFactors)
%       Axes       integer >= 1              axes per epoch, rows epoch-major (default 1)
%       Arcs       as addSystem's Arcs       this stream's own schedule (default: the
%                                           system's)
%     addPrior(sys, name, Name=Value)       a prior group on this system's locals
%       On, L, Value, Sigma, Estimate, From  as LSQVCE.TASKDEFINITION/ADDPRIOR, with L
%                                           q x nl in the columns task.layout() announces
%
%   INPUT   (addSystem's, stored read-only)
%     Name      1x1 string     the system's name, unique in the task
%     Prepare   function_handle | []   @(day, ctx)
%     Settings  1x1 struct     your values; ctx.Settings in Prepare and Build
%     Days      1 x D double   the system's calendar ([] = the task's)
%     Arcs      [] | double | function_handle   the arc schedule
%     Id        [] | number | char   the system's identifier
%   OUTPUT
%     sys       1x1 lsqvce.SystemDefinition   a handle: add-calls change the task
%
%   ACCURACY
%     Not applicable: bookkeeping (names, scopes, handles); no floating-point
%     result. Its compiled system is pinned bitwise against the hand-written
%     one in test_taskDefinition.
%
%   ERRORS
%     lsqvce:TaskDefinition:duplicateName    a local, stream or prior name used twice
%     lsqvce:TaskDefinition:badName          a stream name that is not a MATLAB identifier
%     lsqvce:TaskDefinition:badValid         Valid is not [] or [start end], start <= end
%     lsqvce:TaskDefinition:infoRows         Info has not one row per name
%     lsqvce:TaskDefinition:badHandle        Build is not a function handle of (day, ctx)
%     lsqvce:TaskDefinition:missingArgument  addStream without Build or Noise; a prior
%                                            without Sigma
%     lsqvce:TaskDefinition:priorTarget      a prior gives both On and L, or neither
%     lsqvce:TaskDefinition:badArcs          (at declaration) some streams have a
%                                            schedule and others none
%     lsqvce:TaskDefinition:incomplete       (at declaration) a system without a stream,
%                                            or without a calendar
%     MATLAB:validators:*                    Scope, Estimate, Axes of the wrong kind
%
%   EXAMPLE
%     task = lsqvce.TaskDefinition(Days=59396:59426);
%     task.addCommon(compose("C%d_%d", [2 2 2].', [0 1 2].'));
%     grace = task.addSystem("GRACE-FO", Prepare=@prepareGrace, Arcs=4);
%     grace.addLocal(["x0" "y0" "z0" "vx0" "vy0" "vz0"], Scope="arc");
%     grace.addLocal("accScale", Scope="day");
%     grace.addPrior("accScale", On="^accScale", Value=1, Sigma=0.02, Estimate="scale");
%     grace.addStream("kbr", Build=@buildKbr, Noise=lsqvce.Noise.white(1e-7));
%
%   See also LSQVCE.TASKDEFINITION, LSQVCE.NOISE, LSQVCE.DAYCONTEXT, LSQVCE.DAYDATA
%
%   Created 2026-11-02 (CET, UTC+01:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- level v7, WP-160.

    properties (SetAccess = private)
        Name (1,1) string = ""
        Prepare = []
        Settings (1,1) struct = struct()
        Days double = []
        Arcs = []
        Id = []
    end

    properties (Access = private)
        Locals = struct('Name', {}, 'Scope', {}, 'Valid', {}, 'Info', {})
        Streams = struct('Name', {}, 'Build', {}, 'Noise', {}, 'Estimate', {}, ...
                         'ArcFactors', {}, 'Axes', {}, 'Arcs', {})
        Priors = struct('Name', {}, 'On', {}, 'L', {}, 'Value', {}, 'Sigma', {}, ...
                        'Estimate', {}, 'From', {})
    end

    methods (Access = {?lsqvce.TaskDefinition})
        function obj = SystemDefinition(name, opts)
            obj.Name = name;
            obj.Prepare = opts.Prepare;
            obj.Settings = opts.Settings;
            obj.Days = [];
            if ~isempty(opts.Days), obj.Days = reshape(opts.Days, 1, []); end
            obj.Arcs = opts.Arcs;
            obj.Id = opts.Id;
        end
    end

    methods
        function addLocal(obj, names, opts)
            arguments
                obj (1,1) lsqvce.SystemDefinition
                names {mustBeText, mustBeNonzeroLengthText}
                opts.Scope (1,1) string {mustBeMember(opts.Scope, ["arc", "day", "system"])} = "arc"
                opts.Valid double = []
                opts.Info table = table()
            end
            names = lsqvce.TaskDefinition.nameList(names);
            lsqvce.TaskDefinition.checkNew(names, obj.localNames(), 'local parameter');
            valid = lsqvce.TaskDefinition.checkValid(opts.Valid);
            info = lsqvce.TaskDefinition.infoRows(opts.Info, numel(names));
            for k = 1:numel(names)
                obj.Locals(end + 1) = struct('Name', names(k), 'Scope', opts.Scope, ...
                                             'Valid', valid, 'Info', info{k});
            end
        end

        function addStream(obj, name, opts)
            arguments
                obj (1,1) lsqvce.SystemDefinition
                name (1,1) string {mustBeNonzeroLengthText}
                opts.Build = []
                opts.Noise lsqvce.Noise {mustBeScalarOrEmpty} = lsqvce.Noise.empty
                opts.Estimate string {mustBeScalarOrEmpty, ...
                                      mustBeMember(opts.Estimate, ["law", "scale", "none"])} = string.empty
                opts.ArcFactors logical {mustBeScalarOrEmpty} = logical.empty
                opts.Axes (1,1) double {mustBeInteger, mustBePositive} = 1
                opts.Arcs = []
            end
            if ~isvarname(name)
                error('lsqvce:TaskDefinition:badName', ...
                      ['Stream name "%s" must be a valid MATLAB identifier (letters, digits, ' ...
                       'underscores, starting with a letter).'], name);
            end
            lsqvce.TaskDefinition.checkNew(name, string({obj.Streams.Name}), 'stream');
            if isempty(opts.Build) || isempty(opts.Noise)
                error('lsqvce:TaskDefinition:missingArgument', ...
                      'addStream("%s", ...) needs Build=@build and Noise=lsqvce.Noise.<law>(...).', name);
            end
            lsqvce.TaskDefinition.checkHandle(opts.Build, 'Build');
            obj.Streams(end + 1) = struct('Name', name, 'Build', opts.Build, 'Noise', opts.Noise, ...
                                          'Estimate', opts.Estimate, 'ArcFactors', opts.ArcFactors, ...
                                          'Axes', opts.Axes, 'Arcs', {opts.Arcs});
        end

        function addPrior(obj, name, opts)
            arguments
                obj (1,1) lsqvce.SystemDefinition
                name (1,1) string {mustBeNonzeroLengthText}
                opts.On (1,1) string = ""
                opts.L double = []
                opts.Value double = 0
                opts.Sigma double {mustBeScalarOrEmpty, mustBePositive, mustBeFinite} = []
                opts.Estimate (1,1) string {mustBeMember(opts.Estimate, ["none", "scale"])} = "none"
                opts.From (1,1) double {mustBeInteger, mustBePositive} = 1
            end
            obj.Priors(end + 1) = lsqvce.TaskDefinition.priorEntry(name, opts, string({obj.Priors.Name}));
        end
    end

    methods (Access = {?lsqvce.TaskDefinition})
        function n = localNames(obj)
            n = string({obj.Locals.Name});
        end

        function tf = usesValidity(obj)
            tf = ~all(cellfun(@isempty, {obj.Locals.Valid}));
        end

        function s = summary(obj)
            sc = string({obj.Locals.Scope});
            prep = 'no Prepare';
            if ~isempty(obj.Prepare), prep = ['Prepare ' func2str(obj.Prepare)]; end
            s = sprintf('system "%s": locals %d (arc %d, day %d, system %d); streams %s; %s; %d prior(s)', ...
                        obj.Name, numel(sc), sum(sc == "arc"), sum(sc == "day"), sum(sc == "system"), ...
                        strjoin(string({obj.Streams.Name}), ', '), prep, numel(obj.Priors));
        end

        function s = compile(obj, commons, useValid, taskHasDays)
            % COMPILE  The engines' system struct; Build and Prepare behind generated handles.
            if isempty(obj.Streams)
                error('lsqvce:TaskDefinition:incomplete', ...
                      'System "%s" has no stream: add one with addStream.', obj.Name);
            end
            if isempty(obj.Days) && ~taskHasDays
                error('lsqvce:TaskDefinition:incomplete', ...
                      'System "%s" has no calendar: give the task Days=..., or the system Days=....', obj.Name);
            end
            own = ~cellfun(@isempty, {obj.Streams.Arcs});
            sched = own | ~isempty(obj.Arcs);
            if any(sched) && ~all(sched)
                error('lsqvce:TaskDefinition:badArcs', ...
                      ['System "%s": give every stream an arc schedule (Arcs on the system or the ' ...
                       'stream), or none.'], obj.Name);
            end
            arcs = obj.Arcs;
            if any(own)
                arcs = struct();
                for j = 1:numel(obj.Streams)
                    a = obj.Streams(j).Arcs;
                    if isempty(a), a = obj.Arcs; end
                    arcs.(char(obj.Streams(j).Name)) = a;
                end
            end
            cov = struct();
            spec = struct('System', obj.Name, 'HasPrepare', ~isempty(obj.Prepare), ...
                          'CommonNames', string({commons.Name}), 'CommonScopes', string({commons.Scope}), ...
                          'CommonInfo', {{commons.Info}}, 'LocalNames', obj.localNames(), ...
                          'LocalScopes', string({obj.Locals.Scope}), 'LocalInfo', {{obj.Locals.Info}}, ...
                          'Prepare', {obj.Prepare}, 'Streams', struct('Name', {}, 'Build', {}, ...
                          'Axes', {}, 'Arcs', {}));
            for j = 1:numel(obj.Streams)
                st = obj.Streams(j);
                c = st.Noise.spec();
                if ~isempty(st.Estimate), c.estimate = st.Estimate; end
                if ~isempty(st.ArcFactors), c.arcFactors = st.ArcFactors; end
                cov.(char(st.Name)) = c;
                a = st.Arcs;
                if isempty(a), a = obj.Arcs; end
                spec.Streams(j) = struct('Name', st.Name, 'Build', st.Build, 'Axes', st.Axes, 'Arcs', {a});
            end
            prep = [];
            if ~isempty(obj.Prepare)
                prep = @(day, common, local, procOpt, info) ...
                    lsqvce.DayContext.runPrepare(spec, day, procOpt, info);
            end
            loc = [];
            if ~isempty(obj.Locals), loc = lsqvce.TaskDefinition.templateTable(obj.Locals, useValid); end
            s = struct('buildDay', @(day, obsType, common, local, procOpt, info) ...
                           lsqvce.DayContext.runBuild(spec, day, obsType, procOpt, info), ...
                       'covariance', cov, 'localParams', {loc}, 'arcs', {arcs}, ...
                       'days', obj.Days, 'id', obj.Id, 'name', char(obj.Name), ...
                       'processingOptions', obj.Settings, ...
                       'priors', lsqvce.TaskDefinition.compilePriors(obj.Priors), ...
                       'prepareDay', {prep});
        end
    end
end
