classdef TaskDefinition < handle
%LSQVCE.TASKDEFINITION  Declare a least-squares task: parameters, systems, streams, priors.
%
%   task = lsqvce.TaskDefinition(Days=days)
%   task = lsqvce.TaskDefinition(Days=days, Name=label)
%   task = lsqvce.TaskDefinition.fromMatrices(A, l)
%   task = lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p)
%   task = lsqvce.TaskDefinition.fromBatches(@batchFcn, NumBatches=K, NumUnknowns=u)
%   task = lsqvce.TaskDefinition.fromBatches(@batchFcn, NumBatches=K, X0=x0)
%
%   A TaskDefinition is a handle object you fill with flat add-calls. Every
%   call checks its arguments at once and names the one at fault; nothing is
%   read until LSQVCE.FIT(task, ...) runs. disp(task) prints what was
%   declared: parameter templates per scope, systems, streams and priors.
%
%   WHERE THE DAILY WORK GOES. A system's Prepare function (addSystem) runs
%   once per system and day, before its streams; each stream's Build function
%   (addStream) runs once per stream and day and returns that day's arcs.
%   Both are called as  f(day, ctx)  with a LSQVCE.DAYCONTEXT.
%
%   The task compiles to the declaration struct of the engines
%   (task.declaration()); the estimate, the noise model and every diagnostic
%   are the engines' own -- this class does no arithmetic.
%
%   METHODS
%     addCommon(task, names, Name=Value)    common parameters (every system sees them)
%       Scope  "common" (default) | "day"   constant over the run, or one per calendar day
%       Valid  [start end] fractional MJD | [] (default: always)   when a "day" row exists
%       Info   table, one row per name      your own columns, readable in ctx.Common
%     sys = addSystem(task, name, Name=Value)   an observing system (LSQVCE.SYSTEMDEFINITION)
%       Prepare  @(day, ctx) -> anything    once per system and day; result = ctx.Prepared
%       Settings struct                     paths, models, thresholds; ctx.Settings
%       Days     integer row vector         the system's calendar (default: the task's)
%       Arcs     count | id vector | @(day, stream) -> ids   the arc schedule; needed by
%                                           Scope="arc" locals (default: none, the arcs
%                                           are what Build returns)
%       Id       number | single character  an identifier (give every system one, or none)
%     addPrior(task, name, Name=Value)      a prior group on the common parameters
%       On       regular expression on parameter names (give On or L)
%       L        q x nc operator in the columns task.layout() announces
%       Value    a-priori value, scalar or q x 1 (default 0)
%       Sigma    one standard deviation for the group (required)
%       Estimate "none" (default) | "scale"   estimate the group's sigma
%       From     pass at which the group becomes active (default 1)
%     lay = layout(task)                    the columns of A and B, with names
%     data = declaration(task)              the engine-tier declaration struct
%   STATIC CONSTRUCTORS (in-memory tasks; both take addPrior with L only)
%     fromMatrices(A, l, Weights=p)         A m x u, l m x 1, p m x 1 weights (1./sigma.^2)
%     fromBatches(@batchFcn, NumBatches=K, NumUnknowns=u)   [A, l, p] = batchFcn(k)
%     fromBatches(@batchFcn, NumBatches=K, X0=x0)          [J, w, p] = batchFcn(k, x)
%
%   INPUT
%     Days  1 x D double   integer day identifiers (e.g. MJD), each once
%     Name  1x1 string     a label for the run (default "")
%   OUTPUT
%     task  1x1 lsqvce.TaskDefinition   the declaration, empty until filled
%
%   ACCURACY
%     Not applicable: bookkeeping -- names, scopes and handles turned into the
%     engines' declaration; no floating-point result. A fit of a
%     TaskDefinition equals the fit of the hand-written declaration BITWISE
%     (test_taskDefinition, testCompiledTaskFitsBitwiseLikeTheStruct).
%
%   ERRORS  (raised by the call that makes the mistake, before any data is read)
%     lsqvce:TaskDefinition:duplicateName    a name is used twice: commons, the
%                                            locals of a system, a local named like
%                                            a common, systems, streams, priors
%     lsqvce:TaskDefinition:badName          a stream name that is not a valid
%                                            MATLAB identifier (it becomes a field)
%     lsqvce:TaskDefinition:badValid         Valid is not [] or [start end] with
%                                            start <= end
%     lsqvce:TaskDefinition:infoRows         Info has not one row per name
%     lsqvce:TaskDefinition:badHandle        Prepare or Build is not a function
%                                            handle taking (day, ctx)
%     lsqvce:TaskDefinition:missingArgument  a required Name-Value argument is
%                                            missing (Build, Noise, Sigma, NumBatches)
%     lsqvce:TaskDefinition:priorTarget      a prior gives both On and L, or neither
%     lsqvce:TaskDefinition:wrongKind        addCommon / addSystem / a prior with On
%                                            on an in-memory task
%     lsqvce:TaskDefinition:badArcs          a system mixes scheduled and
%                                            unscheduled streams
%     lsqvce:TaskDefinition:incomplete       declaration() of a task without a
%                                            common parameter, a system, a stream
%                                            per system, or a calendar
%     MATLAB:validators:*                    a Name-Value value of the wrong kind
%                                            (Scope, Estimate, Days, ...)
%
%   EXAMPLE
%     task = lsqvce.TaskDefinition(Days=1:3);
%     task.addCommon(["a" "b"]);
%     sys = task.addSystem("line");
%     sys.addStream("y", Build=@buildLine, Noise=lsqvce.Noise.white(0.1));
%     R = lsqvce.fit(task);
%     % with  function out = buildLine(day, ctx)
%     %           t = (0:99).' + 100 * day;   out = ctx.newDay();
%     %           out.addArc([ones(100, 1), t / 1e3], 2 + 3e-3 * t + 0.1 * randn(100, 1), t);
%     %       end
%
%   See also LSQVCE.SYSTEMDEFINITION, LSQVCE.NOISE, LSQVCE.DAYCONTEXT, LSQVCE.DAYDATA,
%            LSQVCE.FIT, LSQVCE.LAYOUT, LSQVCE.OPTIONS
%
%   Created 2026-11-02 (CET, UTC+01:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- level v7, WP-160.

    properties (SetAccess = private)
        Kind (1,1) string = "streaming"   % "streaming" | "matrices" | "batched"
        Days double = []                  % the default calendar (1 x D)
        Name (1,1) string = ""            % the run label
    end

    properties (Access = private)
        Commons = struct('Name', {}, 'Scope', {}, 'Valid', {}, 'Info', {})
        Systems = {}                      % 1 x M cell of lsqvce.SystemDefinition
        Priors = struct('Name', {}, 'On', {}, 'L', {}, 'Value', {}, 'Sigma', {}, ...
                        'Estimate', {}, 'From', {})
        InMemory = struct()               % the matrices / batched declaration
    end

    methods
        function obj = TaskDefinition(opts)
            arguments
                opts.Days double {mustBeInteger, mustBeFinite} = []
                opts.Name (1,1) string = ""
            end
            if ~isempty(opts.Days), obj.Days = reshape(opts.Days, 1, []); end
            obj.Name = opts.Name;
        end

        function addCommon(obj, names, opts)
            arguments
                obj (1,1) lsqvce.TaskDefinition
                names {mustBeText, mustBeNonzeroLengthText}
                opts.Scope (1,1) string {mustBeMember(opts.Scope, ["common", "day"])} = "common"
                opts.Valid double = []
                opts.Info table = table()
            end
            obj.requireKind("streaming", 'addCommon');
            names = lsqvce.TaskDefinition.nameList(names);
            lsqvce.TaskDefinition.checkNew(names, string({obj.Commons.Name}), 'common parameter');
            valid = lsqvce.TaskDefinition.checkValid(opts.Valid);
            info = lsqvce.TaskDefinition.infoRows(opts.Info, numel(names));
            for k = 1:numel(names)
                obj.Commons(end + 1) = struct('Name', names(k), 'Scope', opts.Scope, ...
                                              'Valid', valid, 'Info', info{k});
            end
        end

        function sys = addSystem(obj, name, opts)
            arguments
                obj (1,1) lsqvce.TaskDefinition
                name (1,1) string {mustBeNonzeroLengthText}
                opts.Prepare = []
                opts.Settings (1,1) struct = struct()
                opts.Days double {mustBeInteger, mustBeFinite} = []
                opts.Arcs = []
                opts.Id = []
            end
            obj.requireKind("streaming", 'addSystem');
            have = string(cellfun(@(s) s.Name, obj.Systems, 'UniformOutput', false));
            lsqvce.TaskDefinition.checkNew(name, have, 'system');
            if ~isempty(opts.Prepare)
                lsqvce.TaskDefinition.checkHandle(opts.Prepare, 'Prepare');
            end
            sys = lsqvce.SystemDefinition(name, opts);
            obj.Systems{end + 1} = sys;
        end

        function addPrior(obj, name, opts)
            arguments
                obj (1,1) lsqvce.TaskDefinition
                name (1,1) string {mustBeNonzeroLengthText}
                opts.On (1,1) string = ""
                opts.L double = []
                opts.Value double = 0
                opts.Sigma double {mustBeScalarOrEmpty, mustBePositive, mustBeFinite} = []
                opts.Estimate (1,1) string {mustBeMember(opts.Estimate, ["none", "scale"])} = "none"
                opts.From (1,1) double {mustBeInteger, mustBePositive} = 1
            end
            if opts.On ~= "" && obj.Kind ~= "streaming"
                error('lsqvce:TaskDefinition:wrongKind', ...
                      ['addPrior(..., On=...) matches parameter names; an in-memory task has ' ...
                       'none -- give the operator L instead.']);
            end
            obj.Priors(end + 1) = lsqvce.TaskDefinition.priorEntry(name, opts, string({obj.Priors.Name}));
        end

        function lay = layout(obj)
            % LAYOUT  The columns of A and of every system's B, with names (see LSQVCE.LAYOUT).
            obj.requireKind("streaming", 'layout');
            lay = lsqvce.engine.layoutDeclaration(obj.declaration());
            names = string({obj.Commons.Name});
            for j = 1:numel(lay.common)
                lay.common(j).name = names(lay.common(j).row);
            end
            for k = 1:numel(lay.systems)
                lay.systems(k).name = obj.Systems{k}.Name;
                lnames = obj.Systems{k}.localNames();
                for c = 1:numel(lay.systems(k).local)
                    lay.systems(k).local(c).name = lnames(lay.systems(k).local(c).row);
                end
            end
        end

        function data = declaration(obj)
            % DECLARATION  The engine-tier declaration struct this task compiles to.
            switch obj.Kind
                case "matrices"
                    data = obj.InMemory;
                case "batched"
                    data = obj.InMemory;
                otherwise
                    data = obj.streamingDeclaration();
            end
            if ~isempty(obj.Priors)
                data.priors = lsqvce.TaskDefinition.compilePriors(obj.Priors);
            end
        end

        function disp(obj)
            if obj.Kind ~= "streaming"
                fprintf('  lsqvce.TaskDefinition %s (%s), %d prior group(s)\n', ...
                        lsqvce.TaskDefinition.quoted(obj.Name), obj.Kind, numel(obj.Priors));
                return
            end
            sc = string({obj.Commons.Scope});
            fprintf('  lsqvce.TaskDefinition %s, %d calendar day(s)\n', ...
                    lsqvce.TaskDefinition.quoted(obj.Name), numel(obj.Days));
            fprintf('    common parameters: %d template(s) (common %d, day %d)\n', ...
                    numel(sc), sum(sc == "common"), sum(sc == "day"));
            for k = 1:numel(obj.Systems)
                fprintf('    %s\n', obj.Systems{k}.summary());
            end
            for j = 1:numel(obj.Priors)
                fprintf('    prior %s\n', lsqvce.TaskDefinition.priorSummary(obj.Priors(j)));
            end
        end
    end

    methods (Static)
        function task = fromMatrices(A, l, opts)
            arguments
                A (:,:) {mustBeNumeric}
                l (:,1) {mustBeNumeric}
                opts.Weights (:,1) double = []
            end
            task = lsqvce.TaskDefinition();
            task.Kind = "matrices";
            task.InMemory = struct('A', A, 'l', l);
            if ~isempty(opts.Weights), task.InMemory.p = opts.Weights; end
        end

        function task = fromBatches(batchFcn, opts)
            arguments
                batchFcn (1,1) function_handle
                opts.NumBatches double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []
                opts.NumUnknowns double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []
                opts.X0 (:,1) double = []
            end
            if isempty(opts.NumBatches)
                error('lsqvce:TaskDefinition:missingArgument', ...
                      'fromBatches needs NumBatches=K, the number of batches batchFcn serves.');
            end
            task = lsqvce.TaskDefinition();
            task.Kind = "batched";
            task.InMemory = struct('batchFcn', batchFcn, 'numBatches', opts.NumBatches);
            if ~isempty(opts.NumUnknowns), task.InMemory.numUnknowns = opts.NumUnknowns; end
            if ~isempty(opts.X0), task.InMemory.x0 = opts.X0; end
        end
    end

    methods (Static, Hidden)
        function names = nameList(names)
            % NAMELIST  Text input -> 1 x n string row.
            names = reshape(string(names), 1, []);
        end

        function checkNew(names, have, what)
            % CHECKNEW  Refuse a name used twice (within NAMES or against HAVE).
            names = reshape(string(names), 1, []);
            [u, ~, j] = unique(names);
            dup = u(accumarray(j(:), 1) > 1);
            if isempty(dup), dup = intersect(names, have); end
            if ~isempty(dup)
                error('lsqvce:TaskDefinition:duplicateName', ...
                      'The %s name "%s" is used twice -- every name is unique.', what, dup(1));
            end
        end

        function v = checkValid(v)
            % CHECKVALID  [] or [start end] with start <= end (fractional MJD).
            if isempty(v)
                v = [];
                return
            end
            if numel(v) ~= 2 || v(1) > v(2)
                error('lsqvce:TaskDefinition:badValid', ...
                      'Valid must be [] (always) or [start end] in fractional MJD with start <= end.');
            end
            v = reshape(v, 1, 2);
        end

        function rows = infoRows(T, n)
            % INFOROWS  One struct per name from the Info table (empty structs without one).
            rows = repmat({struct()}, 1, n);
            if isempty(T) && width(T) == 0, return, end
            if height(T) ~= n
                error('lsqvce:TaskDefinition:infoRows', ...
                      'Info must have one row per name (%d), not %d.', n, height(T));
            end
            for k = 1:n, rows{k} = table2struct(T(k, :)); end
        end

        function checkHandle(h, what)
            % CHECKHANDLE  A function handle that can take (day, ctx).
            if ~isa(h, 'function_handle') || ~isscalar(h)
                error('lsqvce:TaskDefinition:badHandle', ...
                      '%s must be a function handle @(day, ctx), not a %s.', what, class(h));
            end
            try
                n = nargin(h);
            catch ME %#ok<NASGU> -- PROCEED: a handle not on the path yet is checked by fit,
                return                     % which refuses it under lsqvce:makeLoader:buildDay
            end
            if n >= 0 && n < 2
                error('lsqvce:TaskDefinition:badHandle', ...
                      '%s must take two arguments, (day, ctx); %s takes %d.', what, func2str(h), n);
            end
        end

        function e = priorEntry(name, opts, have)
            % PRIORENTRY  One prior group, validated (shared by the task and its systems).
            lsqvce.TaskDefinition.checkNew(name, have, 'prior');
            if (opts.On ~= "") == ~isempty(opts.L)
                error('lsqvce:TaskDefinition:priorTarget', ...
                      'Prior "%s": give either On (a name pattern) or L (an operator), not both or neither.', ...
                      name);
            end
            if isempty(opts.Sigma)
                error('lsqvce:TaskDefinition:missingArgument', ...
                      'Prior "%s" needs Sigma=s, one standard deviation for the group.', name);
            end
            e = struct('Name', name, 'On', opts.On, 'L', opts.L, 'Value', opts.Value, ...
                       'Sigma', opts.Sigma, 'Estimate', opts.Estimate, 'From', opts.From);
        end

        function P = compilePriors(G)
            % COMPILEPRIORS  Prior entries -> the engines' prior struct array.
            P = struct('group', {}, 'L', {}, 'x0', {}, 'sigma', {}, 'from', {}, ...
                       'name', {}, 'estimate', {});
            for j = 1:numel(G)
                grp = [];
                if G(j).On ~= "", grp = char(G(j).On); end
                P(j) = struct('group', grp, 'L', G(j).L, 'x0', G(j).Value, 'sigma', G(j).Sigma, ...
                              'from', G(j).From, 'name', char(G(j).Name), ...
                              'estimate', G(j).Estimate);
            end
        end

        function s = priorSummary(g)
            % PRIORSUMMARY  One line per prior group for disp.
            if g.On ~= ""
                what = sprintf('On %s', g.On);
            else
                what = sprintf('L %d x %d', size(g.L, 1), size(g.L, 2));
            end
            s = sprintf('%s: %s, sigma %g, estimate %s, from pass %d', ...
                        g.Name, what, g.Sigma, g.Estimate, g.From);
        end

        function s = quoted(name)
            % QUOTED  "name" or (unnamed).
            if name == "", s = '(unnamed)'; else, s = sprintf('"%s"', name); end
        end
    end

    methods (Access = private)
        function requireKind(obj, kind, what)
            if obj.Kind ~= kind
                error('lsqvce:TaskDefinition:wrongKind', ...
                      '%s is for a streaming task; this task was made by from%s.', what, ...
                      lsqvce.TaskDefinition.kindWord(obj.Kind));
            end
        end

        function data = streamingDeclaration(obj)
            % STREAMINGDECLARATION  The kind-A struct: commonParams, days, systems, paramColumns.
            if isempty(obj.Commons) || isempty(obj.Systems)
                error('lsqvce:TaskDefinition:incomplete', ...
                      'A task needs at least one addCommon and one addSystem before it can be fitted.');
            end
            cnames = string({obj.Commons.Name});
            useValid = ~all(cellfun(@isempty, {obj.Commons.Valid}));
            for k = 1:numel(obj.Systems)
                useValid = useValid || obj.Systems{k}.usesValidity();
            end
            data = struct();
            data.commonParams = lsqvce.TaskDefinition.templateTable(obj.Commons, useValid);
            if ~isempty(obj.Days), data.days = obj.Days; end
            S = cell(1, numel(obj.Systems));
            for k = 1:numel(obj.Systems)
                S{k} = obj.Systems{k}.compile(obj.Commons, useValid, ~isempty(obj.Days));
            end
            data.systems = [S{:}];
            if useValid, data.paramColumns = struct('name', 2, 'start', 3, 'end', 4); end
            % a local named like a common would make ctx.cols and the priors ambiguous
            for k = 1:numel(obj.Systems)
                clash = intersect(obj.Systems{k}.localNames(), cnames);
                if ~isempty(clash)
                    error('lsqvce:TaskDefinition:duplicateName', ...
                          'System "%s" has a local parameter "%s", which is also a common parameter.', ...
                          obj.Systems{k}.Name, clash(1));
                end
            end
        end
    end

    methods (Static, Hidden)
        function T = templateTable(entries, useValid)
            % TEMPLATETABLE  {scope, name[, start, end]} rows of the engines' template table.
            n = numel(entries);
            T = cell(n, 2 + 2 * useValid);
            for r = 1:n
                T{r, 1} = char(entries(r).Scope);
                T{r, 2} = char(entries(r).Name);
                if useValid
                    v = entries(r).Valid;
                    if isempty(v), T{r, 3} = []; T{r, 4} = []; else, T{r, 3} = v(1); T{r, 4} = v(2); end
                end
            end
        end

        function w = kindWord(kind)
            if kind == "matrices", w = 'Matrices'; else, w = 'Batches'; end
        end
    end
end
