classdef DayData < handle
%LSQVCE.DAYDATA  The arcs one Build returns for one stream and one day.
%
%   out = ctx.newDay()                              made by LSQVCE.DAYCONTEXT only
%   out.addArc(A, l, t)                             one arc, no local parameters
%   out.addArc(A, l, t, Local={names, P; ...})      local partials by name
%   out.addArc(A, l, t, B=B)                        local partials as one m x nl matrix
%   out.addArc(..., Sigma=s, Excluded=e, Arc=id)
%
%   Build adds the arcs of its day one by one and returns OUT. Every row is
%   RAW: observed minus computed, never whitened or scaled by a weight -- the
%   toolbox whitens each arc from the stream's LSQVCE.NOISE and its per-epoch
%   Sigma, and eliminates the local parameters arc by arc. An arc may be
%   empty (0 rows): a data gap.
%
%   LOCAL PARTIALS, ONE OF TWO WAYS PER ARC
%     Local={names, P}   P m x numel(names); several rows {names1, P1; names2, P2}.
%                        The scope of each name decides its column: this arc's
%                        block (Scope="arc"), this day's (Scope="day") or the
%                        system's (Scope="system"). The usual way.
%     B=B                m x ctx.NumLocal (may be sparse), every column taken from
%                        ctx.cols(names, Arc=id) -- for partials that already form
%                        one matrix. Never computes an offset by hand.
%   Neither: the arc sees no local parameter (B is all zero).
%
%   ARC SCHEDULE. Without a declared schedule (no Arcs on the system or
%   stream) the arcs are those added, in order. With one, the arcs are
%   ctx.ArcIds: addArc fills them in order, Arc=id places an arc out of
%   order, and a scheduled id never added is a data gap.
%
%   INPUT (addArc)
%     A         m x nc double      design of the common parameters, nc = ctx.NumCommon,
%                                  columns in addCommon order (Scope="day" commons at
%                                  ctx.cols(names))
%     l         m x 1 double       observed minus computed
%     t         m x 1 double       epoch of each row in SECONDS
%   NAME-VALUE (addArc)
%     Local     n x 2 cell         {names, partials} rows (default: none)
%     B         m x nl double      the local design as one matrix (default: none)
%     Sigma     e x 1 | e x d | e x d x d double   the per-epoch accuracy: standard
%                                  deviations per epoch, per axis, or covariances
%                                  (e = m / d epochs); relative to the law's sigma
%                                  (default: none, D = I)
%     Excluded  m x 1 logical      rows of weight 0, zero in A, B and l, kept for their
%                                  index (default: none)
%     Arc       scheduled arc id   place this arc (default: the next scheduled one)
%   OUTPUT
%     out       1x1 lsqvce.DayData   a handle: addArc adds to it
%
%   ACCURACY
%     Not applicable: it places given numbers into the engines' arc cells
%     without arithmetic. Local= and the equivalent B= give the bitwise same
%     fit (test_taskDefinition).
%
%   ERRORS  (raised while the day is read; fit passes them unchanged)
%     lsqvce:DayContext:localAndB     one arc gives both Local and B
%     lsqvce:DayContext:badLocal      Local is not an n x 2 cell of {names, P} with P
%                                     m x numel(names), or names a column twice
%     lsqvce:DayContext:badArc        Arc without a schedule, an id the day does not
%                                     schedule, an arc added twice, or more arcs than
%                                     scheduled
%     lsqvce:DayContext:unknownName   a Local name that is not a local parameter of the
%                                     system on this day
%     lsqvce:DayContext:arcRequired   an arc-scoped Local without a scheduled arc
%     lsqvce:makeLoader:*             the engines' checks of the day (row counts,
%                                     Sigma shapes), with day and stream named
%
%   EXAMPLE
%     out = ctx.newDay();
%     out.addArc(A1, l1, t1, Local={["x0" "y0" "z0"], P1; "kbrBias", ones(numel(t1), 1)});
%     out.addArc(A2, l2, t2, Local={["x0" "y0" "z0"], P2; "kbrBias", ones(numel(t2), 1)}, ...
%                Excluded=flagged2);
%
%   See also LSQVCE.DAYCONTEXT, LSQVCE.SYSTEMDEFINITION, LSQVCE.NOISE
%
%   Created 2026-11-02 (CET, UTC+01:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- level v7, WP-160.

    properties (SetAccess = private)
        Day double = []
        Stream (1,1) string = ""
        System (1,1) string = ""
    end

    properties (Access = private)
        Ctx                       % the lsqvce.DayContext of this call
        Arcs = struct('Pos', {}, 'A', {}, 'B', {}, 'l', {}, 't', {}, 'Sigma', {}, 'Excluded', {})
        Next = 1                  % the next scheduled position
    end

    methods (Access = {?lsqvce.DayContext})
        function obj = DayData(ctx)
            obj.Ctx = ctx;
            obj.Day = ctx.Day;
            obj.Stream = ctx.Stream;
            obj.System = ctx.System;
        end
    end

    methods
        function addArc(obj, A, l, t, opts)
            arguments
                obj (1,1) lsqvce.DayData
                A (:,:) {mustBeNumeric}
                l (:,1) {mustBeNumeric}
                t (:,1) {mustBeNumeric}
                opts.Local cell = {}
                opts.B {mustBeNumeric} = []
                opts.Sigma {mustBeNumeric} = []
                opts.Excluded {mustBeNumericOrLogical} = []
                opts.Arc double {mustBeScalarOrEmpty} = []
            end
            if ~isempty(opts.Local) && ~isempty(opts.B)
                error('lsqvce:DayContext:localAndB', ...
                      ['Stream "%s", day %s: an arc gives its local partials either by name ' ...
                       '(Local=) or as a matrix (B=), not both.'], obj.Stream, string(obj.Day));
            end
            pos = obj.position(opts.Arc);
            ids = obj.Ctx.ArcIds;
            arcId = [];
            if ~isempty(ids), arcId = ids(pos); end
            m = size(A, 1);
            if ~isempty(opts.B)
                B = opts.B;
            elseif ~isempty(opts.Local)
                B = obj.placeLocal(opts.Local, m, arcId);
            elseif obj.Ctx.NumLocal == 0
                B = zeros(m, 0);
            else
                B = sparse(m, obj.Ctx.NumLocal);
            end
            obj.Arcs(end + 1) = struct('Pos', pos, 'A', A, 'B', B, 'l', l, 't', t, ...
                                       'Sigma', opts.Sigma, 'Excluded', opts.Excluded);
        end
    end

    methods (Hidden)
        function tf = belongsTo(obj, ctx)
            % BELONGSTO  True for the DayData of this context's day and stream.
            tf = isequal(obj.Day, ctx.Day) && obj.Stream == ctx.Stream && obj.System == ctx.System;
        end

        function D = toEngineDay(obj)
            % TOENGINEDAY  The engines' day: 1 x nArc cells A, B, l, t (+ Cep, excluded, nAxis).
            ids = obj.Ctx.ArcIds;
            n = numel(obj.Arcs);
            if ~isempty(ids), n = numel(ids); end
            [A, B, l, t, C, E] = deal(cell(1, n));
            for k = 1:numel(obj.Arcs)
                a = obj.Arcs(k);
                A{a.Pos} = a.A;  B{a.Pos} = a.B;  l{a.Pos} = a.l;  t{a.Pos} = a.t;
                C{a.Pos} = a.Sigma;  E{a.Pos} = a.Excluded;
            end
            D = struct('A', {A}, 'B', {B}, 'l', {l}, 't', {t});
            if ~all(cellfun(@isempty, C)), D.Cep = C; end
            if ~all(cellfun(@isempty, E)), D.excluded = E; end
            ax = obj.Ctx.streamAxes();
            if ax > 1, D.nAxis = ax; end
        end
    end

    methods (Access = private)
        function pos = position(obj, arc)
            ids = obj.Ctx.ArcIds;
            used = [obj.Arcs.Pos];
            if isempty(ids)
                if ~isempty(arc)
                    error('lsqvce:DayContext:badArc', ...
                          ['Stream "%s", day %s: Arc=%g needs an arc schedule (Arcs on the system or ' ...
                           'stream); without one the arcs are those added, in order.'], ...
                          obj.Stream, string(obj.Day), arc);
                end
                pos = numel(obj.Arcs) + 1;
                return
            end
            if isempty(arc)
                pos = obj.Next;
            else
                pos = find(ids == arc, 1);
                if isempty(pos)
                    error('lsqvce:DayContext:badArc', ...
                          'Stream "%s", day %s: arc %g is not scheduled (scheduled: %s).', ...
                          obj.Stream, string(obj.Day), arc, mat2str(ids));
                end
            end
            if pos > numel(ids) || any(used == pos)
                error('lsqvce:DayContext:badArc', ...
                      ['Stream "%s", day %s: the day schedules the arcs %s; an arc was added twice ' ...
                       'or after the last one.'], obj.Stream, string(obj.Day), mat2str(ids));
            end
            obj.Next = pos + 1;
        end

        function B = placeLocal(obj, spec, m, arcId)
            if size(spec, 2) ~= 2
                error('lsqvce:DayContext:badLocal', ...
                      'Local must be an n x 2 cell {names, partials; ...} (stream "%s").', obj.Stream);
            end
            B = sparse(m, obj.Ctx.NumLocal);
            seen = zeros(1, 0);
            for r = 1:size(spec, 1)
                names = spec{r, 1};
                P = spec{r, 2};
                if ~(ischar(names) || isstring(names) || iscellstr(names)) || ~isnumeric(P)
                    error('lsqvce:DayContext:badLocal', ...
                          'Local row %d must be {names, partials} with text names and numeric partials.', r);
                end
                names = reshape(string(names), 1, []);
                if size(P, 1) ~= m || size(P, 2) ~= numel(names)
                    error('lsqvce:DayContext:badLocal', ...
                          'Local row %d: the partials must be %d x %d (rows x names), not %s.', ...
                          r, m, numel(names), mat2str(size(P)));
                end
                c = obj.Ctx.localColumns(names, arcId);
                if any(ismember(c, seen)) || numel(unique(c)) < numel(c)
                    error('lsqvce:DayContext:badLocal', ...
                          'Local names one parameter twice in one arc (stream "%s").', obj.Stream);
                end
                seen = [seen, c]; %#ok<AGROW> a handful of rows
                B(:, c) = P;
            end
        end
    end
end
