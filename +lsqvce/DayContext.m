classdef DayContext
%LSQVCE.DAYCONTEXT  What Prepare and Build know about the day they serve (read-only).
%
%   P   = prepare(day, ctx)                 a system's Prepare (addSystem)
%   out = build(day, ctx)                   a stream's Build (addStream)
%   out  = ctx.newDay()                     an empty LSQVCE.DAYDATA for this stream
%   cols = ctx.cols(names)                  columns of named parameters on this day
%   cols = ctx.cols(names, Arc=id)          the same for Scope="arc" locals of arc id
%
%   The toolbox makes one DayContext per call of your Prepare and Build. It
%   carries the day, the names of the system and stream, your Settings, the
%   result of Prepare, the parameter templates by name, and the arc schedule
%   of the day. Build fills ctx.newDay() with arcs (LSQVCE.DAYDATA/ADDARC) and
%   returns it.
%
%   COLUMNS BY NAME. Local partials are normally given BY NAME
%   (out.addArc(A, l, t, Local={names, P})): the toolbox places them. When the
%   partials already form one matrix B (m x ctx.NumLocal), take every column
%   from ctx.cols and never compute an offset:
%       B = sparse(m, ctx.NumLocal);
%       B(:, ctx.cols(["x0" "y0" "z0"], Arc=id)) = Pxyz;
%       B(:, ctx.cols("accScale")) = Pacc;
%       out.addArc(A, l, t, B=B);
%   ctx.cols of common parameters gives columns of A (a Scope="day" common
%   moves with the day); of local parameters, columns of B. One call names
%   commons or locals, not both.
%
%   PROPERTIES (read-only)
%     Day        1x1 double     the day being built (an entry of the calendar)
%     System     1x1 string     the system's name
%     Stream     1x1 string     the stream's name ("" in Prepare)
%     Settings   1x1 struct     the system's Settings
%     Prepared   any            what Prepare returned for this day ([] in Prepare,
%                               or without Prepare)
%     NumCommon  1x1 double     nc, the columns of A
%     NumLocal   1x1 double     nl, the columns of B of this system (0 without locals)
%     ArcIds     1 x A double   the stream's scheduled arc ids of this day; in
%                               Prepare, the system's scheduled arcs of the day ([]
%                               when the system declares no schedule: the arcs are
%                               what Build adds)
%     Common     table          Name, Scope, OnDay, Info of every common template
%     Local      table          Name, Scope, OnDay, Info of every local template
%   METHODS
%     newDay(ctx)                  -> LSQVCE.DAYDATA
%     cols(ctx, names, Arc=id)     names 1 x n string -> 1 x n double columns;
%                                  Arc (a scheduled arc id) is needed by Scope="arc"
%                                  locals and ignored by the others
%
%   INPUT   (made by the toolbox; you never construct one)
%   OUTPUT  ctx  1x1 lsqvce.DayContext
%
%   ACCURACY
%     Not applicable: integer bookkeeping over the engines' column plans; the
%     columns equal the engines' layout (test_taskDefinition).
%
%   ERRORS  (raised while the day is read; fit passes them unchanged)
%     lsqvce:DayContext:unknownName   a name that is not a parameter of this task
%                                     and system, or not instantiated on this day
%     lsqvce:DayContext:arcRequired   a Scope="arc" local asked for without Arc=id,
%                                     or with an id the day does not schedule
%     lsqvce:DayContext:mixedTables   one cols call names commons and locals
%     lsqvce:DayContext:badReturn     Build returned something other than the
%                                     LSQVCE.DAYDATA of ctx.newDay()
%
%   EXAMPLE
%     function out = buildKbr(day, ctx)
%         P = ctx.Prepared;                         % what Prepare returned today
%         out = ctx.newDay();
%         for a = 1:numel(P.arc)
%             out.addArc(P.arc(a).A, P.arc(a).omc, P.arc(a).t, ...
%                        Local={"kbrBias", ones(numel(P.arc(a).t), 1)});
%         end
%     end
%
%   See also LSQVCE.DAYDATA, LSQVCE.SYSTEMDEFINITION, LSQVCE.TASKDEFINITION
%
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: ArcIds in Prepare are the
%   system's scheduled arcs of the day (they were empty). With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Created 2026-11-02 (CET, UTC+01:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- level v7, WP-160.

    properties (SetAccess = private)
        Day double = []
        System (1,1) string = ""
        Stream (1,1) string = ""
        Settings (1,1) struct = struct()
        Prepared = []
        NumCommon double = 0
        NumLocal double = 0
        ArcIds double = []
    end

    properties (Dependent)
        Common
        Local
    end

    properties (Access = private)
        Spec = struct()        % the compiled system (names, scopes, info)
        Info = []              % the engines' column slice of this day
        StreamIndex = 0
    end

    methods (Access = private)
        function obj = DayContext(spec, k, day, procOpt, info)
            obj.Spec = spec;
            obj.StreamIndex = k;
            obj.Day = day;
            obj.System = spec.System;
            if spec.HasPrepare && isfield(procOpt, 'prepared')
                obj.Prepared = procOpt.prepared;
                procOpt = rmfield(procOpt, 'prepared');
            end
            obj.Settings = procOpt;
            obj.Info = info;
            if isstruct(info) && isfield(info, 'nc'), obj.NumCommon = info.nc; end
            if isstruct(info) && isfield(info, 'nl'), obj.NumLocal = info.nl; end
            if k > 0
                st = spec.Streams(k);
                obj.Stream = st.Name;
                obj.ArcIds = lsqvce.DayContext.scheduleOf(st.Arcs, day, char(st.Name));
            elseif isstruct(info) && isfield(info, 'arcId')
                obj.ArcIds = reshape(double(info.arcId), 1, []);   % Prepare: the system's arcs of the day
            else
                obj.ArcIds = lsqvce.DayContext.scheduleOf(spec.SystemArcs, day, '');
            end
        end
    end

    methods
        function out = newDay(obj)
            out = lsqvce.DayData(obj);
        end

        function c = cols(obj, names, opts)
            arguments
                obj (1,1) lsqvce.DayContext
                names {mustBeText, mustBeNonzeroLengthText}
                opts.Arc double {mustBeScalarOrEmpty} = []
            end
            names = reshape(string(names), 1, []);
            isC = ismember(names, obj.Spec.CommonNames);
            if any(isC) && ~all(isC)
                error('lsqvce:DayContext:mixedTables', ...
                      ['ctx.cols(...) names common and local parameters in one call: commons are ' ...
                       'columns of A, locals of B -- ask for each separately.']);
            end
            if all(isC)
                c = obj.commonColumns(names);
            else
                c = obj.localColumns(names, opts.Arc);
            end
        end

        function T = get.Common(obj)
            s = obj.Spec;
            on = true(numel(s.CommonNames), 1);
            for r = 1:numel(s.CommonNames)
                if s.CommonScopes(r) == "day", on(r) = ~isempty(obj.commonDayColumn(r)); end
            end
            T = table(s.CommonNames(:), s.CommonScopes(:), on, s.CommonInfo(:), ...
                      'VariableNames', {'Name', 'Scope', 'OnDay', 'Info'});
        end

        function T = get.Local(obj)
            s = obj.Spec;
            on = false(numel(s.LocalNames), 1);
            for r = 1:numel(s.LocalNames)
                on(r) = obj.localPresent(r);
            end
            T = table(s.LocalNames(:), s.LocalScopes(:), on, s.LocalInfo(:), ...
                      'VariableNames', {'Name', 'Scope', 'OnDay', 'Info'});
        end
    end

    methods (Static, Hidden)
        function D = runBuild(spec, day, obsType, procOpt, info)
            % RUNBUILD  The generated buildDay: the stream's Build, its arcs as the engines' day.
            k = find([spec.Streams.Name] == string(obsType), 1);
            ctx = lsqvce.DayContext(spec, k, day, procOpt, info);
            out = spec.Streams(k).Build(day, ctx);
            if ~(isa(out, 'lsqvce.DayData') && isscalar(out) && out.belongsTo(ctx))
                error('lsqvce:DayContext:badReturn', ...
                      ['Build of stream "%s" (system "%s") on day %s returned a %s: return the ' ...
                       'DayData of ctx.newDay(), filled with addArc.'], ...
                      ctx.Stream, ctx.System, string(day), class(out));
            end
            D = out.toEngineDay();
        end

        function P = runPrepare(spec, day, procOpt, info)
            % RUNPREPARE  The generated prepareDay: the system's Prepare with a context.
            ctx = lsqvce.DayContext(spec, 0, day, procOpt, info);
            P = spec.Prepare(day, ctx);
        end

        function ids = scheduleOf(arcs, day, obsType)
            % SCHEDULEOF  The scheduled arc ids of a day ([] = no schedule), as the engines read it.
            if isempty(arcs)
                ids = [];
            elseif isa(arcs, 'function_handle')
                if nargin(arcs) == 1, ids = arcs(day); else, ids = arcs(day, obsType); end
            elseif isscalar(arcs)
                ids = 1:arcs;
            else
                ids = arcs;
            end
            ids = reshape(double(ids), 1, []);
        end
    end

    methods (Hidden)
        function a = streamAxes(obj)
            % STREAMAXES  The stream's declared axes per epoch (addStream Axes).
            a = obj.Spec.Streams(obj.StreamIndex).Axes;
        end

        function c = localColumns(obj, names, arcId)
            % LOCALCOLUMNS  Columns of B of the named locals (arcId for Scope="arc" ones).
            s = obj.Spec;
            c = zeros(1, numel(names));
            for k = 1:numel(names)
                r = find(s.LocalNames == names(k), 1);
                if isempty(r)
                    error('lsqvce:DayContext:unknownName', ...
                          '"%s" is not a local parameter of system "%s".', names(k), s.System);
                end
                c(k) = obj.localColumn(r, arcId, names(k));
            end
        end
    end

    methods (Access = private)
        function c = commonColumns(obj, names)
            s = obj.Spec;
            c = zeros(1, numel(names));
            for k = 1:numel(names)
                r = find(s.CommonNames == names(k), 1);
                if s.CommonScopes(r) == "common"
                    c(k) = sum(s.CommonScopes(1:r) == "common");
                else
                    cc = obj.commonDayColumn(r);
                    if isempty(cc)
                        error('lsqvce:DayContext:unknownName', ...
                              'The day common "%s" is not instantiated on day %s.', names(k), string(obj.Day));
                    end
                    c(k) = cc;
                end
            end
        end

        function c = commonDayColumn(obj, r)
            c = [];
            I = obj.Info;
            if ~isstruct(I) || ~isfield(I, 'commonDayOffset') || isempty(I.commonDayOffset), return, end
            q = find(I.commonDayRows == r, 1);
            if ~isempty(q), c = I.commonDayOffset + q; end
        end

        function tf = localPresent(obj, r)
            I = obj.Info;
            tf = false;
            if ~isstruct(I) || ~isfield(I, 'sysRows'), return, end
            switch obj.Spec.LocalScopes(r)
                case "arc",    tf = ~isempty(I.arcId) && any(I.arcRows == r);
                case "day",    tf = ~isempty(I.dayOffset) && any(I.dayRows == r);
                otherwise,     tf = any(I.sysRows == r);
            end
        end

        function c = localColumn(obj, r, arcId, name)
            I = obj.Info;
            if ~obj.localPresent(r)
                error('lsqvce:DayContext:unknownName', ...
                      'The local "%s" of system "%s" is not instantiated on day %s.', ...
                      name, obj.Spec.System, string(obj.Day));
            end
            switch obj.Spec.LocalScopes(r)
                case "arc"
                    if isempty(arcId) || ~any(I.arcId == arcId)
                        error('lsqvce:DayContext:arcRequired', ...
                              ['The arc-scoped local "%s" needs the id of one of this day''s scheduled ' ...
                               'arcs (%s): ctx.cols(names, Arc=id).'], name, mat2str(I.arcId));
                    end
                    c = I.offsetOf(arcId) + find(I.arcRows == r, 1);
                case "day"
                    c = I.dayOffset + find(I.dayRows == r, 1);
                otherwise
                    c = I.sysOffset + find(I.sysRows == r, 1);
            end
        end
    end
end
