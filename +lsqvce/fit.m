function R = fit(data, opts)
%LSQVCE.FIT  Fit a least-squares adjustment (optionally with VCE) to declared data.
%
%   R = lsqvce.fit(data)
%   R = lsqvce.fit(data, opts)
%
%   The one verb (R4, HANDOVER §11.4). You describe the problem with a plain
%   struct `data` (DATA_CONTRACT.md: which unknowns are common, which days,
%   which observing systems and the ONE function buildDay each of them points
%   to), choose what to do with an LSQVCE.OPTIONS object, and get an
%   LSQVCE.RESULTS back:
%
%       opts = lsqvce.Options(Model="vce", GroupBy="satObs");
%       R    = lsqvce.fit(data, opts);
%       R.report
%
%   fit does three things and nothing else: (1) validates `data` against the
%   contract BEFORE any data is read (fail fast, one error identifier per
%   rule -- ERRORS below), (2) builds one Loader per system through
%   LSQVCE.ENGINE.MAKELOADER (the field names are makeLoader's canonical
%   option names in lowerCamelCase; the mapping table lives here and nowhere
%   else), (3) dispatches on opts.Model to the engine that does the work:
%
%       Model      engine                    what
%       "gls"      LSQVCE.ENGINE.ADJUSTBLOCK        one GLS solve with the declared
%                                            covariance
%       "vce"      LSQVCE.ENGINE.ADJUSTVCE          alternate stochastic model <->
%                                            variance components (separate)
%       "pooled"   LSQVCE.ENGINE.ADJUSTVCE.joint    pooled spectrum, no arc factors
%       "ellmer"   LSQVCE.ENGINE.ADJUSTVCEJOINT     the ten-step joint co-estimation
%                                            with keyed arc factors (Ellmer
%                                            2018 §6.5; R4 WP-56c). fit builds
%                                            the accFcn the driver wants -- a
%                                            BlockAccumulator over the same
%                                            loaders whose UserOpt carries the
%                                            current model through
%                                            LSQVCE.ENGINE.JOINTMODELOPT -- and
%                                            bridges the driver's out struct
%                                            into a Results: x / sigma0 /
%                                            diagnostics via diagnoseAccum +
%                                            diagnose, the converged noise
%                                            model on R.noise (.C .Shat .sigma2
%                                            .arcKeys), the pass record on
%                                            R.history (.kind = 'ellmer',
%                                            .maxDev, .Shat, .nIter, .skips),
%                                            R.converged. summary()/
%                                            convergence() are adjustVCE's
%                                            readers and refuse (noHistory).
%
%   THE v4 CONTRACT (Level L11, WP-75) STATES NOTHING THE DATA ALREADY
%   STATES, and the parameter tables are SCOPED: column 1 of BOTH tables is
%   a scope token, the rows are TEMPLATES, and the toolbox expands them over
%   the calendar -- nc and nl are DERIVED from the expansion, the observation
%   types = fieldnames(covariance), u = numel(x0) on the nonlinear batched
%   kind. commonParams tokens: 'common' (one column, always) | 'day' (one
%   column per calendar day that carries at least one scheduled arc in at
%   least one system). localParams tokens: 'arc' (one block per (day, arc)
%   of the system) | 'day' (one block per day with arcs) | 'system' (once).
%   Binding column layout -- A: ['common' rows, table order | day blocks,
%   calendar order]; B: per day ['arc' blocks, sorted ids | 'day' block],
%   then the 'system' block. buildDay's 6th argument carries every offset
%   (below), so no offset is ever computed in user code.
%
%   TWO MORE KINDS OF DATA (R4 WP-56b) -- told apart by the fields present,
%   for problems WITHOUT local unknowns and with diagonal weights:
%
%   (B) IN-MEMORY MATRICES
%         data.A   m x u double     design matrix
%         data.l   m x 1 double     observed minus computed
%         data.p   m x 1 double     weights (diagonal of P); optional (unit)
%   (C) BATCHED CALLBACK
%         data.batchFcn     fh      [A, l, p] = batchFcn(k), k = 1..nBatch
%                                   (the LSQVCE.ENGINE.ADJUST contract)
%         data.numBatches   1x1     nBatch
%         data.numUnknowns  1x1     u -- LINEAR form only. With
%         data.x0           u x 1   -> NONLINEAR: [A, w, p] = batchFcn(k, x)
%                                   (the LSQVCE.ENGINE.SOLVENONLINEAR
%                                   contract); u = numel(x0), so numUnknowns
%                                   must then be ABSENT (stated once).
%   Dispatch for both: opts.Model must be "gls" (variance components need
%   the groups of a streaming system); opts.Robust = "irls" selects
%   LSQVCE.ENGINE.ROBUST, data.x0 selects LSQVCE.ENGINE.SOLVENONLINEAR, else LSQVCE.ENGINE.ADJUST.
%   The Results carries sol / diag / source (the driver's third output:
%   verdict V, ROB, or the nonlinear IT record). Its sigma0 is the
%   A-POSTERIORI std of unit weight and sol.Cxx = s0sq * D.Qxx, i.e. the
%   streaming Results contract (the batch drivers' own sol.sigma0 is the
%   a-priori value; it is kept as sol.sigma0Prior). R.acc is empty, so
%   residuals() and contribution() refuse with lsqvce:Results:noAccumulator
%   -- the per-observation residuals are already in R.diag.v. The options
%   forwarded to LSQVCE.ENGINE.DIAGNOSE come from the Options object explicitly
%   (no KeepUnmatched on this path).
%
%   R.systems (the existing Results field -- no second name) carries the
%   normalised system specs the engine ran: .sysID .nl .dataFcn (the loader
%   closure fit built) .dayList .name; every lower layer stays reachable on
%   R (R.acc, R.diag).
%
%   INPUT
%     data  1x1  struct  the data declaration (binding rules: DATA_CONTRACT.md
%                        §2). Summary of the streaming kind:
%              .commonParams  RxC cell    the COMMON templates, one ROW per
%                                         template; column 1 = 'common'|'day',
%                                         columns 2..C your convention. nc is
%                                         DERIVED by the expansion. Required.
%                                         Passed to buildDay as its 3rd arg.
%              .days     1xD double   default day list (numeric integers);
%                                     required unless every system has .days
%              .systems  1xM struct   one entry per observing system, fields:
%                REQUIRED
%                  .buildDay          fh        your day builder
%                  .covariance        struct    .(obsType) = covariance spec;
%                                               the observation types ARE its
%                                               fieldnames (their order is the
%                                               processing order)
%                OPTIONAL (default)
%                  .arcs (1)          1x1 | 1xA | struct | @(day,obsType)
%                                     ('fromDay' is ENGINE-tier only)
%                  .localParams  ({}) RxC cell of LOCAL templates; column 1 =
%                                     'arc'|'day'|'system', columns 2..C your
%                                     convention; nl DERIVED by the expansion;
%                                     buildDay's 4th arg. Not declared ->
%                                     nl = 0 (B is m x 0).
%                  .days ([]) .name ('') .id (the system's index)
%                  .processingOptions (struct())   buildDay's 5th arg
%                  .epochCovariance ([])  @(day,obsType) per-epoch D hook:
%                                     the exact separable-descale path of the
%                                     VCE Shape step (DATA_CONTRACT §7).
%              buildDay's 6th argument (info) carries the expansion: .nl,
%              .arcId / .offsetOf(arcId) / .count, .dayOffset / .dayCount,
%              .sysOffset / .sysCount for B; .nc, .commonCount,
%              .commonDayOffset / .commonDayCount for THIS day's columns in A.
%     opts  1x1  lsqvce.Options   what to do (default: every engine default,
%                        i.e. Model "gls")
%
%   OUTPUT
%     R     1x1  lsqvce.Results  solution, diagnostics, provenance; for Model
%                        "vce"/"pooled" also R.history / R.noise / R.converged
%
%   ACCURACY
%     No numerics of its own. The claim is EXACTNESS of the dispatch: for the
%     same declaration, fit(data, opts) returns the results of the direct
%     driver call it stands for -- adjustBlock with
%     opts.engineArgs("adjustBlock"), and adjustVCE with the configurators
%     applied from opts.engineArgs("adjustVCE") -- to solver precision:
%     x, sigma0 and Cxx agree within RelTol 1e-12 on the CI fixtures
%     (dev/tests/test_fit.m, reference). It is the SAME code path; the 1e-12
%     rather than isequal only guards the multithreaded-BLAS reduction-order
%     nondeterminism recorded in HANDOVER §10.8 (WP-47 r1e).
%
%   ERRORS
%     Your declaration is wrong (raised BEFORE any data is read):
%       lsqvce:fit:unknownField    a field name the contract does not know
%                                  (names the field; typos are errors on
%                                  purpose -- and so are the retired fields:
%                                  common, numCommon, numLocal,
%                                  observationTypes, systemId, cache, loader,
%                                  plotPostSolve, plotPostModel, and since
%                                  v4 localsPerArc: declare 'arc' rows in
%                                  localParams instead)
%       lsqvce:fit:noSystems       data.systems missing, not a struct, or empty
%       lsqvce:fit:noCommon        data.commonParams missing, or not a 2-D
%                                  cell of template rows
%       lsqvce:fit:badScope        a scope-token column 1 does not carry:
%                                  commonParams needs 'common'|'day' (raised
%                                  here), localParams needs
%                                  'arc'|'day'|'system' (raised by makeLoader
%                                  and relabelled)
%       lsqvce:fit:missingField    a required system field is absent, or
%                                  covariance carries no observation type
%       lsqvce:fit:noDays          no day list anywhere
%       lsqvce:fit:badId           id is not a numeric scalar or single char
%       lsqvce:fit:duplicateSystem two systems share an id (explicit or
%                                  index default)
%       lsqvce:fit:builderKind     buildDay is not a function handle
%       lsqvce:fit:badArcs         arcs = 'fromDay' (ENGINE-tier only: it
%                                  reads data before the run; declare a
%                                  resolver @(day,obsType) instead)
%       lsqvce:fit:<makeLoader id> a makeLoader construction error, relabelled
%                                  (e.g. paramType, badArcs, localPlan);
%                                  the original is attached as the cause
%     Your choice is not available here:
%       lsqvce:fit:notYet          Robust ~= "none" together with Model ~= "gls" --
%                                  (R5 M2 2b) the one-loop composition of M4;
%                                  before 2b this refused Robust on the streaming
%                                  kind altogether --
%                                  the streaming kind screens outliers
%                                  through its diagnostics (report())
%       lsqvce:fit:modelKind       Model other than "gls" on kinds B / C
%       lsqvce:batchSystem:weightMatrix  kind B/C returned a full weight matrix
%                                  (R5 M2: the block path takes diagonal weights;
%                                  a full P is a covariance law of the streaming
%                                  kind) -- raised by lsqvce.engine.batchSystem
%       lsqvce:fit:badEstimate     (R5 M3 3d) covariance.<type>.estimate is not
%                                  "none" | "scale" | "law"
%       lsqvce:fit:badMatrix       kind B shapes disagree (A m x u, l m x 1,
%                                  p m x 1)
%       lsqvce:fit:missingField    kind C lacks batchFcn / numBatches, or the
%                                  linear form lacks numUnknowns
%       lsqvce:fit:uCount          kind C carries BOTH x0 and numUnknowns --
%                                  u is stated once (numel(x0))
%     Your buildDay is wrong (raised WHILE READING, by the engine, ids
%     unchanged): lsqvce:makeLoader:badDay / dayShape / badArc / badArcId /
%     buildDay, lsqvce:accum:arcBadNc / arcBadNl
%
%   EXAMPLE (the GRACE parameterization that motivated v4: arc-wise ICs,
%   day-wise ACC scale/bias/drift, KBR empiricals per arc)
%     data.commonParams = [repmat({'common'}, 6, 1), num2cell((1:6).')];
%     data.days = 1:3;
%     s.buildDay = @myBuildDay;                  % uses the info argument
%     s.covariance.pos = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',0.7);
%     s.covariance.kbr = struct('type','diag','sigma',0.2);
%     s.localParams = [ ...
%         [repmat({'arc'}, 6, 1),  num2cell((1:6).')];  ...  % ICs per arc
%         {'arc', 'kbrBias'; 'arc', 'kbrCos1cpr'};      ...  % KBR per arc
%         [repmat({'day'}, 3, 1),  {'accScale'; 'accBias'; 'accDrift'}]];
%     s.arcs = 2;
%     data.systems = s;          % nl = 3 days * (2 arcs * 8 + 3) = 57, derived
%     R = lsqvce.fit(data, lsqvce.Options(Model="vce"));
%     R.report
%
%   See also LSQVCE.OPTIONS, LSQVCE.RESULTS, LSQVCE.ENGINE.MAKELOADER,
%            LSQVCE.ENGINE.ADJUSTBLOCK, LSQVCE.ENGINE.ADJUSTVCE, LSQVCE.ENGINE.ADJUSTVCEJOINT,
%            LSQVCE.ENGINE.ADJUST, LSQVCE.ENGINE.ROBUST, LSQVCE.ENGINE.SOLVENONLINEAR
%
%   Updated 2026-08-31 (CEST, UTC+02:00), v4 SCOPED tables (Level L11,
%   WP-75, four pinned user decisions): scope column on BOTH tables
%   (commonParams 'common'|'day', localParams 'arc'|'day'|'system'), rows
%   are templates expanded over the calendar (nc/nl derived; fit builds the
%   COMMON plan over the global day calendar and passes it to makeLoader);
%   localsPerArc left this surface (unknownField -- declare 'arc' rows);
%   localCount retired, badScope new. With assistance from Claude
%   (Anthropic), model Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-72): every derivable count derived (nc / nl / obsTypes / u-with-x0),
%   commonParams required and RxC, localParams | localsPerArc | neither
%   (nl = 0), systemId -> id with the index default, the raw-loader escape
%   hatch and the v2 count/behaviour fields removed from this surface
%   (cache / plot hooks / 'fromDay' stay ENGINE-tier on makeLoader).
%   With assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-09-01 (CEST, UTC+02:00), R4 WP-56c (Level L3): Model "ellmer"
%   wired -- accFcn from the loaders + jointModelOpt, out -> Results bridge
%   (diagnoseAccum + diagnose; noise/history/converged set). With assistance
%   from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), R4 WP-56b (Level L4): kinds B
%   (in-memory matrices) and C (batched callback, linear or nonlinear with
%   x0) through adjust / robust / solveNonlinear, returned as Results without
%   an accumulator. With assistance from Claude (Anthropic), model Claude
%   Fable 5.
%   Created 2026-08-31 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5 -- R4 Level L2, WP-56 (decision R4-D4).
arguments
    data (1,1) struct
    opts (1,1) lsqvce.Options = lsqvce.Options()
end

    opts.validate();

    % ---- kinds B / C (R4 WP-56b): no systems, no loaders -------------------
    f = fieldnames(data);
    if ~any(strcmp(f, 'systems')) && (any(strcmp(f, 'A')) || any(strcmp(f, 'l')))
        R = fitMatrix(data, opts);
        return
    elseif ~any(strcmp(f, 'systems')) && any(strcmp(f, 'batchFcn'))
        R = fitBatch(data, opts);
        return
    end

    % ---- kind A: streaming --------------------------------------------------
    lsqvce.util.log.inventory('reset');          % R5 M1: this run's warning inventory
    spec = validateData(data);
    systems = buildSystems(spec);
    days = spec.days;

    % R5 M2 2b: robust IRLS is a POLICY on the block path -- on the streaming
    % kind it is Model "gls" with reweighting (H-02: the same accumulator);
    % combining it with a VCE model is a v5 loop feature, refused for now.
    if ~strcmp(opts.Robust, "none")
        if ~strcmp(opts.Model, "gls")
            error('lsqvce:fit:notYet', ...
                ['Robust = "%s" together with Model = "%s" is the one-loop policy ' ...
                 'composition of R5 M4; run the robust pass with Model "gls".'], ...
                opts.Robust, opts.Model);
        end
        a = opts.engineArgs('adjustRobust');
        prov = opts.toStruct();  prov.kind = 'stream';
        R = lsqvce.engine.adjustRobust(systems, days, a.nv{:}, 'Provenance', prov);
        return
    end

    switch char(opts.Model)
        case 'gls'
            a = opts.engineArgs('adjustBlock');
            R = lsqvce.engine.adjustBlock(systems, days, a.nv{:});
        case 'ellmer'
            R = fitEllmer(systems, days, opts);
        otherwise                                   % 'vce' | 'pooled'
            a = opts.engineArgs('adjustVCE');
            % constructor-only configuration since 2.0 (R4-D2): the groups
            % ride as Name-Value structs; VceMode in a.nv selects the mode
            P = lsqvce.engine.adjustVCE(systems, days, a.nv{:}, ...
                    'Shape', a.shape, 'Scale', a.scale, 'Outer', a.outer, ...
                    'Joint', a.joint, 'Axes', a.axes);
            R = P.run();
    end
end

%==========================================================================
function spec = validateData(data)
%VALIDATEDATA  Enforce DATA_CONTRACT.md §2 (v3) on the raw struct; fail fast.
%   OUTPUT spec 1x1 struct  .commonParams (RxC cell), .days (double),
%                           .systems (1xM), .ids (1xM cellstr, resolved)
    topAllowed = {'commonParams', 'days', 'systems'};
    rejectUnknown(fieldnames(data), topAllowed, 'data');

    if ~isfield(data, 'systems') || ~isstruct(data.systems) || isempty(data.systems)
        error('lsqvce:fit:noSystems', ...
            'data.systems must be a non-empty 1xM struct array (one entry per observing system).');
    end
    sys = data.systems(:).';

    % ---- the common unknowns: the SCOPED table is the declaration (v4).
    % Column 1 = scope token 'common' | 'day'; rows are TEMPLATES ('day' rows
    % expand over the global day calendar); nc is DERIVED by the common plan.
    if ~isfield(data, 'commonParams') || isempty(data.commonParams) || ...
            ~iscell(data.commonParams) || ~ismatrix(data.commonParams)
        error('lsqvce:fit:noCommon', ...
            ['data.commonParams is required: an R x C cell, one ROW per common ' ...
             'template, column 1 the scope token ''common'' | ''day''. ' ...
             'Minimal table: [repmat({''common''}, nc, 1), num2cell((1:nc).'')].']);
    end
    for r = 1:size(data.commonParams, 1)
        v = data.commonParams{r, 1};
        if ~(ischar(v) || (isstring(v) && isscalar(v))) || ...
                ~any(strcmp(char(v), {'common', 'day'}))
            error('lsqvce:fit:badScope', ...
                ['data.commonParams row %d: column 1 must be the scope token ' ...
                 '''common'' | ''day'' (''arc'' is system-specific and lives in ' ...
                 'localParams).'], r);
        end
    end

    % ---- days
    days = [];
    if isfield(data, 'days') && ~isempty(data.days)
        days = double(data.days(:).');
        if ~all(isfinite(days) & days == round(days))
            error('lsqvce:fit:noDays', 'data.days must be numeric integers (the loader cache is keyed on them).');
        end
    end

    % ---- systems
    sysRequired = {'buildDay', 'covariance'};
    sysAllowed  = [sysRequired, {'localParams', 'arcs', 'days', ...
                   'id', 'name', 'processingOptions', 'epochCovariance'}];
    rejectUnknown(fieldnames(sys), sysAllowed, 'data.systems');
    ids = cell(1, numel(sys));
    for k = 1:numel(sys)
        s = sys(k);
        for f = sysRequired
            if ~isfield(s, f{1}) || isempty(s.(f{1}))
                error('lsqvce:fit:missingField', ...
                    'data.systems(%d) needs the field ''%s''.', k, f{1});
            end
        end
        if ~isa(s.buildDay, 'function_handle')
            error('lsqvce:fit:builderKind', ...
                'data.systems(%d): buildDay must be a function handle.', k);
        end
        if ~isstruct(s.covariance) || ~isscalar(s.covariance) || isempty(fieldnames(s.covariance))
            error('lsqvce:fit:missingField', ...
                ['data.systems(%d).covariance must be a scalar struct with one field per ' ...
                 'observation type (the fieldnames ARE the observation types).'], k);
        end
        if isfield(s, 'arcs') && (ischar(s.arcs) || isstring(s.arcs))
            error('lsqvce:fit:badArcs', ...
                ['data.systems(%d).arcs = ''%s'': the string form (''fromDay'') is ' ...
                 'ENGINE-tier only -- it reads data before the run. Declare a count, ' ...
                 'an id vector, or a resolver @(day,obsType).'], k, char(s.arcs));
        end
        if isfield(s, 'id') && ~isempty(s.id)
            id = s.id;
            if ~((isnumeric(id) && isscalar(id)) || (ischar(id) && isscalar(id)) || ...
                 (isstring(id) && isscalar(id) && strlength(id) == 1))
                error('lsqvce:fit:badId', ...
                    'data.systems(%d).id must be a numeric scalar or a single character.', k);
            end
            ids{k} = char(string(id));
        else
            ids{k} = char(string(k));               % the index default
        end
        if isempty(days) && ~(isfield(s, 'days') && ~isempty(s.days))
            error('lsqvce:fit:noDays', ...
                'no day list: set data.days or data.systems(%d).days.', k);
        end
    end
    if numel(unique(ids)) < numel(ids)
        error('lsqvce:fit:duplicateSystem', ...
            'system ids must be unique across data.systems (explicit ids and index defaults share one namespace).');
    end
    spec = struct('commonParams', {data.commonParams}, 'days', days, ...
                  'systems', sys, 'ids', {ids});
    spec.commonPlan = buildCommonPlan(spec);
end

function plan = buildCommonPlan(spec)
%BUILDCOMMONPLAN  Expand the SCOPED commonParams over the GLOBAL day calendar
%   (L11 WP-75, v4). Layout, binding: the 'common' rows first (table order),
%   then one 'day' block per calendar day (sorted union of every declared day
%   list) whose SCHEDULE carries at least one arc in at least one system --
%   schedules are resolvable without reading data by contract, so the plan is
%   built before any read. nc is DERIVED.
%   OUTPUT plan 1x1 struct  .nc .commonCount .days (calendar) .dayOffset
%          .dayCount (1xD) .rows (.common / .day template row indices)
    cp  = spec.commonParams;
    tok = cellfun(@(v) char(string(v)), cp(:, 1), 'UniformOutput', false);
    iCommon = find(strcmp(tok, 'common')).';
    iDay    = find(strcmp(tok, 'day')).';
    nCom = numel(iCommon);  nDayRows = numel(iDay);
    cal = spec.days;
    for k = 1:numel(spec.systems)
        s = spec.systems(k);
        if isfield(s, 'days') && ~isempty(s.days)
            cal = union(cal, double(s.days(:).'));
        end
    end
    cal = sort(cal(:).');
    D = numel(cal);
    dayOffset = zeros(1, D);  dayCount = zeros(1, D);
    off = nCom;
    if nDayRows > 0
        for i = 1:D
            if calendarDayHasArcs(spec.systems, cal(i))
                dayOffset(i) = off;  dayCount(i) = nDayRows;
                off = off + nDayRows;
            end
        end
    end
    plan = struct('nc', off, 'commonCount', nCom, 'days', cal, ...
                  'dayOffset', dayOffset, 'dayCount', dayCount, ...
                  'rows', struct('common', iCommon, 'day', iDay));
end

function tf = calendarDayHasArcs(sys, day)
%CALENDARDAYHASARCS  Does at least one system's SCHEDULE carry an arc on DAY?
%   Mirrors the engine's arc resolution on the declared forms (default 1;
%   count; id vector; per-type struct; @(day,obsType) resolver) -- resolvable
%   pre-run by contract.
    tf = false;
    for k = 1:numel(sys)
        s = sys(k);
        if isfield(s, 'days') && ~isempty(s.days) && ~any(double(s.days) == day)
            continue                                   % this system skips the day
        end
        if ~isfield(s, 'arcs') || isempty(s.arcs)
            tf = true;  return                         % default: one arc per day
        end
        a = s.arcs;
        types = fieldnames(s.covariance);
        for j = 1:numel(types)
            aj = a;
            if isstruct(a)
                if ~isfield(a, types{j}), continue, end
                aj = a.(types{j});
            end
            if isa(aj, 'function_handle')
                if ~isempty(aj(day, types{j})), tf = true;  return, end
            elseif isnumeric(aj)
                if (isscalar(aj) && aj >= 1) || (~isscalar(aj) && ~isempty(aj))
                    tf = true;  return
                end
            end
        end
    end
end

function rejectUnknown(names, allowed, where)
%REJECTUNKNOWN  Unknown field names are errors (typo protection, rule 5).
    bad = setdiff(names, allowed);
    if ~isempty(bad)
        error('lsqvce:fit:unknownField', ...
            '%s has unknown field ''%s''. Allowed: %s.', where, bad{1}, strjoin(allowed, ', '));
    end
end

%==========================================================================
function systems = buildSystems(spec)
%BUILDSYSTEMS  One Loader per system through LSQVCE.ENGINE.MAKELOADER.
%   Derivations (v4): ObservationTypes = fieldnames(covariance); nc rides in
%   as the COMMON PLAN (fit expands the scoped commonParams over the global
%   calendar); nl is derived by makeLoader's scope plan from the scoped
%   LocalParams table, or NumLocal = 0 when none is declared.
    sys  = spec.systems;
    days = spec.days;
    % ---- the ONE mapping table: lowerCamelCase field -> makeLoader NV name
    map = {'covariance', 'Covariance'; 'arcs', 'Arcs'; ...
           'localParams', 'LocalParams'; ...
           'days', 'Days'; 'name', 'Name'; ...
           'processingOptions', 'ProcessingOptions'; ...
           'epochCovariance', 'EpochCovariance'};
    systems = cell(1, numel(sys));
    for k = 1:numel(sys)
        s  = sys(k);
        nv = {};
        for r = 1:size(map, 1)
            f = map{r, 1};
            if isfield(s, f) && ~isempty(s.(f))
                nv = [nv, {map{r, 2}, s.(f)}];                       %#ok<AGROW>
            end
        end
        nv = [nv, {'ObservationTypes', fieldnames(s.covariance).'}]; %#ok<AGROW>
        nv = [nv, {'SystemId', idValue(spec.ids{k})}];               %#ok<AGROW>
        nv = [nv, {'CommonParams', spec.commonParams}];              %#ok<AGROW>
        nv = [nv, {'Scoped', true, 'CommonPlan', spec.commonPlan}];  %#ok<AGROW>
        if ~(isfield(s, 'days') && ~isempty(s.days)) && ~isempty(days)
            nv = [nv, {'Days', days}];                               %#ok<AGROW>
        end
        if ~(isfield(s, 'localParams') && ~isempty(s.localParams))
            nv = [nv, {'NumLocal', 0}];                              %#ok<AGROW>
        end
        try
            systems{k} = lsqvce.engine.makeLoader(s.buildDay, nv{:});
        catch e
            if startsWith(e.identifier, 'lsqvce:makeLoader:')
                e2 = MException(strrep(e.identifier, 'lsqvce:makeLoader:', 'lsqvce:fit:'), ...
                                'data.systems(%d): %s', k, e.message);
                throw(addCause(e2, e));
            end
            rethrow(e);
        end
    end
end

function v = idValue(idChar)
%IDVALUE  A resolved id back to the type makeLoader stores: numeric where the
%   char is a number ('1' -> 1, index defaults and explicit numerics), the
%   single char otherwise ('C' -> 'C').
    n = str2double(idChar);
    if isscalar(idChar) && isnan(n)
        v = idChar;
    else
        v = n;
    end
end

%========================================================== kinds B / C
function R = fitMatrix(data, opts)
%FITMATRIX  Kind B: in-memory A, l (, p) -> one-batch driver -> Results.
    rejectUnknown(fieldnames(data), {'A', 'l', 'p'}, 'data');
    for f = {'A', 'l'}
        if ~isfield(data, f{1}) || isempty(data.(f{1}))
            error('lsqvce:fit:missingField', 'data lacks the required field %s.', f{1});
        end
    end
    A = data.A;  l = data.l(:);
    p = [];
    if isfield(data, 'p') && ~isempty(data.p), p = data.p(:); end
    if ~isnumeric(A) || ~ismatrix(A) || size(A, 1) ~= numel(l) || ...
            (~isempty(p) && numel(p) ~= numel(l))
        error('lsqvce:fit:badMatrix', ...
            'data.A must be m x u, data.l m x 1 and data.p (if given) m x 1 (m = %d).', numel(l));
    end
    batchFcn = @(k) deal(A, l, p);
    R = runBatch(size(A, 2), batchFcn, 1, [], opts, 'matrix');
end

function R = fitBatch(data, opts)
%FITBATCH  Kind C: batched callback (linear, or nonlinear with x0) -> Results.
%   v3: u is stated ONCE -- numUnknowns on the linear form, numel(x0) on the
%   nonlinear form (numUnknowns must then be absent, lsqvce:fit:uCount).
    rejectUnknown(fieldnames(data), {'batchFcn', 'numUnknowns', 'numBatches', 'x0'}, 'data');
    for f = {'batchFcn', 'numBatches'}
        if ~isfield(data, f{1}) || isempty(data.(f{1}))
            error('lsqvce:fit:missingField', 'data lacks the required field %s.', f{1});
        end
    end
    x0 = [];
    if isfield(data, 'x0') && ~isempty(data.x0), x0 = data.x0(:); end
    if isempty(x0)
        if ~isfield(data, 'numUnknowns') || isempty(data.numUnknowns)
            error('lsqvce:fit:missingField', ...
                'the linear batched kind lacks numUnknowns (u; the nonlinear form derives it from x0).');
        end
        u = data.numUnknowns;
    else
        if isfield(data, 'numUnknowns') && ~isempty(data.numUnknowns)
            error('lsqvce:fit:uCount', ...
                'u is stated once: the nonlinear form derives u = numel(x0), so omit numUnknowns.');
        end
        u = numel(x0);
    end
    R = runBatch(u, data.batchFcn, data.numBatches, x0, opts, 'batch');
end

function R = runBatch(u, batchFcn, nBatch, x0, opts, kind)
%RUNBATCH  The ONE dispatch for kinds B / C.
%   R5 M2 (WP-81, round 2a): the LINEAR, non-robust case rides the block
%   path -- the batches become the DAYS of one local-free system (finding
%   H-02: a batch problem is a streaming problem with nl = 0), so one
%   accumulator serves every kind and the streaming Results (residuals(),
%   contribution(), the exact per-observation sweep) is available on kinds
%   B / C. Robust IRLS and the nonlinear form still run on the flat drivers
%   until rounds 2b / 2c give them a policy on the block path.
    if ~strcmp(opts.Model, "gls")
        error('lsqvce:fit:modelKind', ...
            ['Model "%s" needs the streaming kind (variance components are per ' ...
             'group of a system); the %s kind is Model "gls".'], opts.Model, kind);
    end
    % R5 M2 2c: DatumB (minimum-constraint datum), Sigma0 (a-priori sigma0)
    % and Sparse are FLAT-path options with no block-path home yet (the
    % block accumulator has no datum machinery; M5's priors give it one).
    % Round 2a routed every linear problem to the block path and would have
    % IGNORED them -- the silent outcome the doctrine forbids. Until M5 a
    % problem that sets any of them runs on the flat drivers, as before.
    flatOnly = ~isempty(opts.DatumB) || opts.Sigma0 ~= 1 || opts.Sparse;
    if isempty(x0) && ~flatOnly
        R = runBatchOnBlockPath(u, batchFcn, nBatch, opts, kind);   % gls, and IRLS (2b)
        return
    end
    if ~isempty(x0) && ~flatOnly && strcmp(opts.Robust, "none")
        % R5 M2 2c': Gauss-Newton as a policy on the block path
        lsqvce.util.log.inventory('reset');
        a = opts.engineArgs('adjustNonlinear');
        prov = opts.toStruct();  prov.kind = kind;
        R = lsqvce.engine.adjustNonlinear(u, batchFcn, nBatch, x0, a.nv{:}, 'Provenance', prov);
        return
    end
    if ~isempty(x0)
        a = opts.engineArgs("solveNonlinear");
        [sol, D, ~, IT] = lsqvce.engine.solveNonlinear(u, batchFcn, nBatch, x0, a.nv{:});
        source = IT;
    elseif strcmp(opts.Robust, "irls")
        a = opts.engineArgs("robust");
        [sol, D, ROB] = lsqvce.engine.robust(u, batchFcn, nBatch, a.nv{:});
        source = ROB;
    else
        a = opts.engineArgs("adjust");
        [sol, D, V] = lsqvce.engine.adjust(u, batchFcn, nBatch, a.nv{:});
        source = V;
    end
    % The Results contract is the STREAMING one (BlockAccumulator.solve):
    % sol.sigma0 = a-posteriori std of unit weight, sol.Cxx = covariance
    % scaled by it. The batch drivers return Accumulator.solve's struct,
    % where sigma0 is the A-PRIORI Sigma0 and there is no Cxx (R, QxxCorr,
    % s0sq instead) -- so the bridge derives both exactly as lsqvce.engine.diagnose
    % does (D.Qxx is the datum-corrected cofactor). When diagnose skipped
    % Qxx (very large u) Cxx stays absent and coefficients() refuses with
    % lsqvce:Results:noCovariance -- honest, never a plausible substitute.
    sol.sigma0Prior = sol.sigma0;
    sol.sigma0      = sqrt(sol.s0sq);
    if isfield(D, 'Qxx') && ~isempty(D.Qxx)
        sol.Cxx = sol.s0sq * D.Qxx;
    end
    resolved = opts.toStruct();
    resolved.kind = kind;
    R = lsqvce.Results(sol, D, [], source, batchFcn, nBatch, [], lsqvce.engine.runMeta(struct(), resolved));
end

%======================================================= Model "ellmer"
function R = fitEllmer(systems, days, opts)
%FITELLMER  Model "ellmer": LSQVCE.ENGINE.ADJUSTVCEJOINT over the built systems.
%   Two bridges (R4 WP-56c): (1) the accFcn the driver wants -- a fresh
%   BlockAccumulator per outer pass whose UserOpt carries the current model
%   through LSQVCE.ENGINE.JOINTMODELOPT for every (system, obsType) stream; (2) the
%   driver's out struct -> a Results with the streaming contract.
    nc    = lsqvce.engine.resolveNc(systems);
    specs = accSpecs(systems, days);
    streams = streamKeys(systems, specs);
    a = opts.engineArgs("adjustVCEjoint");
    psdMode = char(opts.PsdMode);
    sepStreams = separableStreams(systems, specs);        % R5 M3 (H-15)
    [estMap, afMap] = streamPolicy(systems, specs);       % R5 M3 3d (R5-D4)
    accFcn = @(model) jointAccumulator(nc, specs, opts, ...
                 lsqvce.engine.jointModelOpt(model, streams, 'PsdMode', psdMode, ...
                                             'Separable', sepStreams));
    out = lsqvce.engine.adjustVCEjoint(accFcn, a.nv{:}, 'Estimate', estMap, 'ArcFactors', afMap);

    [dsol, batchFcn, nB] = lsqvce.engine.diagnoseAccum(out.acc, out.sol);
    D = lsqvce.engine.diagnose(dsol, batchFcn, 'nBatch', nB, opts.engineArgs("diagnose").nv{:});
    resolved = opts.toStruct();
    resolved.kind = 'stream';
    seeds = struct('jointSeed', opts.Seed);
    R = lsqvce.Results(out.sol, D, out.acc, dsol, batchFcn, nB, specs, lsqvce.engine.runMeta(seeds, resolved));
    hist = out.history;
    hist.kind  = 'ellmer';
    hist.nIter = out.nIter;
    hist.skips = out.skips;
    % arcKeys is a CELL: braces keep struct() from expanding it into a struct
    % ARRAY (the DATA_CONTRACT §2.4 pitfall -- it bit this very line in L3 r1).
    noise = struct('C', out.C, 'Shat', out.Shat, 'sigma2', out.sigma2, 'arcKeys', {out.arcKeys}, ...
                   'Cgroups', out.Cgroups);            % R5 M3 3d: one C per pool group
    R.setIterative(hist, noise, out.converged);
end

function specs = accSpecs(systems, days)
%ACCSPECS  BlockAccumulator.addSystems specs from the built Loaders.
    specs = struct('sysID', {}, 'nl', {}, 'dataFcn', {}, 'dayList', {}, 'name', {});
    for k = 1:numel(systems)
        L = systems{k};
        dl = days;
        if isfield(L, 'days') && ~isempty(L.days), dl = L.days; end
        nm = '';
        if isfield(L, 'name'), nm = char(L.name); end
        specs(k) = struct('sysID', L.sysID, 'nl', L.nl, 'dataFcn', L.closure, ...
                          'dayList', dl, 'name', nm);
    end
end

function keys = streamKeys(systems, specs)
%STREAMKEYS  'satKey|obsType' for every (system, observation type) stream --
%   the names LSQVCE.ENGINE.ADJUSTVCE / LSQVCE.ENGINE.JOINTMODELOPT use.
    keys = {};
    for k = 1:numel(specs)
        ot = cellstr(systems{k}.obsTypes);
        satKey = matlab.lang.makeValidName(char(string(specs(k).sysID)));
        for j = 1:numel(ot)
            keys{end+1} = sprintf('%s|%s', satKey, ot{j});  %#ok<AGROW>
        end
    end
end

function acc = jointAccumulator(nc, specs, opts, uo)
%JOINTACCUMULATOR  One solve-ready accumulator under the current joint model.
    acc = lsqvce.engine.BlockAccumulator(nc, 'Backend', char(opts.Backend), 'KeepLocal', true, ...
              'Reg', char(opts.Reg), 'RegLambda', opts.RegLambda, 'RegMatrix', opts.RegMatrix, ...
              'Parallel', char(opts.Parallel), 'DistributeDays', opts.DistributeDays, ...
              'ParallelThreshold', opts.ParallelThreshold, 'Verbose', false, 'UserOpt', uo);
    acc.addSystems(specs);
end

%======================================================= kinds B / C on the block path
function R = runBatchOnBlockPath(u, batchFcn, nBatch, opts, kind)
%RUNBATCHONBLOCKPATH  Kinds B / C (linear, Model "gls") as ONE local-free
%   streaming system through LSQVCE.ENGINE.BATCHSYSTEM (R5 M2 WP-81 2a/2c').
    L = lsqvce.engine.batchSystem(u, batchFcn, nBatch, 'Name', kind);
    lsqvce.util.log.inventory('reset');
    prov = opts.toStruct();  prov.kind = kind;   % provenance = the user's resolved Options
    if strcmp(opts.Robust, "none")
        a = opts.engineArgs('adjustBlock');
        R = lsqvce.engine.adjustBlock({L}, 1:nBatch, a.nv{:}, 'Provenance', prov);
    else
        a = opts.engineArgs('adjustRobust');
        R = lsqvce.engine.adjustRobust({L}, 1:nBatch, a.nv{:}, 'Provenance', prov);
    end
end

function keys = separableStreams(systems, specs)
%SEPARABLESTREAMS  'satKey|obsType' of every stream whose base law is separable.
%   R5 M3 (finding H-15): those streams keep their type / axis count / Cep
%   under the joint override; the estimated law reaches them as an AR(p)
%   temporal model (LSQVCE.ENGINE.JOINTMODELOPT 'Separable').
    keys = {};
    for k = 1:numel(systems)
        L = systems{k};
        if ~isfield(L, 'cov') || ~isstruct(L.cov), continue, end
        satKey = matlab.lang.makeValidName(char(string(specs(k).sysID)));
        ots = fieldnames(L.cov);
        for j = 1:numel(ots)
            cm = L.cov.(ots{j});
            if isfield(cm, 'type') && strcmpi(char(cm.type), 'separable')
                keys{end+1} = sprintf('%s|%s', satKey, ots{j}); %#ok<AGROW>
            end
        end
    end
end

function [estMap, afMap] = streamPolicy(systems, specs)
%STREAMPOLICY  Per-stream 'estimate' / 'arcFactors' from the covariance specs.
%   R5 M3 3d (WP-83, R5-D4): a stream spec may carry
%     .estimate    "none" | "scale" | "law"   default "law" under Model "ellmer"
%     .arcFactors  logical                    default true
%   Returned as containers.Map keyed 'satKey|obsType' for the joint driver.
    estMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    afMap  = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for k = 1:numel(systems)
        L = systems{k};
        if ~isfield(L, 'cov') || ~isstruct(L.cov), continue, end
        satKey = matlab.lang.makeValidName(char(string(specs(k).sysID)));
        ots = fieldnames(L.cov);
        for j = 1:numel(ots)
            cm = L.cov.(ots{j});
            key = sprintf('%s|%s', satKey, ots{j});
            if isfield(cm, 'estimate') && ~isempty(cm.estimate)
                e = lower(char(string(cm.estimate)));
                if ~any(strcmp(e, {'none', 'scale', 'law'}))
                    error('lsqvce:fit:badEstimate', ...
                          'covariance.%s.estimate must be "none", "scale" or "law" (got "%s").', ots{j}, e);
                end
                estMap(key) = e;
            end
            if isfield(cm, 'arcFactors') && ~isempty(cm.arcFactors)
                afMap(key) = logical(cm.arcFactors);
            end
        end
    end
end
