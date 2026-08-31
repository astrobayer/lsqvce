function R = fit(data, opts)
%LSQVCE.FIT  Fit a least-squares adjustment (optionally with VCE) to declared data.
%
%   R = lsqvce.fit(data)
%   R = lsqvce.fit(data, opts)
%
%   THE ONE VERB (R4, HANDOVER §11.4). You describe the problem with a plain
%   struct `data` (DATA_CONTRACT.md: what is common, which days, which
%   observing systems and the ONE function buildDay each of them points to),
%   choose what to do with an LSQVCE.OPTIONS object, and get an
%   LSQVCE.RESULTS back:
%
%       opts = lsqvce.Options(Model="vce", GroupBy="satObs");
%       R    = lsqvce.fit(data, opts);
%       R.report
%
%   fit does three things and nothing else: (1) validates `data` against the
%   contract BEFORE any data is read (fail fast, one error identifier per
%   rule -- ERRORS below), (2) builds one Loader per system through LSQVCE.MAKELOADER (the
%   field names are makeLoader's canonical option names in lowerCamelCase;
%   the mapping table lives here and nowhere else), (3) dispatches on
%   opts.Model to the engine that does the work:
%
%       Model      engine                    what
%       "gls"      LSQVCE.ADJUSTBLOCK        one GLS solve with the declared
%                                            covariance
%       "vce"      LSQVCE.ADJUSTVCE          alternate stochastic model <->
%                                            variance components (separate)
%       "pooled"   LSQVCE.ADJUSTVCE.joint    pooled spectrum, no arc factors
%       "ellmer"   LSQVCE.ADJUSTVCEJOINT     NOT WIRED YET (WP-56c): fit
%                                            refuses with lsqvce:fit:notYet
%
%   TWO MORE KINDS OF DATA (R4 WP-56b) -- told apart by the fields present,
%   for problems WITHOUT local unknowns and with diagonal weights:
%
%   (B) IN-MEMORY MATRICES
%         data.A   m x u double     design matrix
%         data.l   m x 1 double     observed minus computed
%         data.p   m x 1 double     weights (diagonal of P); optional (unit)
%   (C) BATCHED CALLBACK
%         data.numUnknowns  1x1     u
%         data.numBatches   1x1     nBatch
%         data.batchFcn     fh      [A, l, p] = batchFcn(k), k = 1..nBatch
%                                   (the LSQVCE.ADJUST contract), or, with
%         data.x0           u x 1   -> NONLINEAR: [A, w, p] = batchFcn(k, x)
%                                   (the LSQVCE.SOLVENONLINEAR contract)
%   Dispatch for both: opts.Model must be "gls" (variance components need
%   the groups of a streaming system); opts.Robust = "irls" selects
%   LSQVCE.ROBUST, data.x0 selects LSQVCE.SOLVENONLINEAR, else LSQVCE.ADJUST.
%   The Results carries sol / diag / source (the driver's third output:
%   verdict V, ROB, or the nonlinear IT record). Its sigma0 is the
%   A-POSTERIORI std of unit weight and sol.Cxx = s0sq * D.Qxx, i.e. the
%   streaming Results contract (the batch drivers' own sol.sigma0 is the
%   a-priori value; it is kept as sol.sigma0Prior). R.acc is empty, so
%   residuals() and contribution() refuse with lsqvce:Results:noAccumulator
%   -- the per-observation residuals are already in R.diag.v. The options
%   forwarded to LSQVCE.DIAGNOSE come from the Options object explicitly
%   (no KeepUnmatched on this path).
%
%   R.systems (the existing Results field -- no second name) carries the
%   normalised system specs the engine ran: .sysID .nl .dataFcn (the loader
%   closure fit built) .dayList .name; every lower layer stays reachable on
%   R (R.acc, R.diag).
%
%   INPUT
%     data  1x1  struct  the data declaration (binding rules: DATA_CONTRACT.md
%                        §2). Summary:
%              .common   1x1 struct  exactly ONE of
%                          .commonParams  Rx14 cell  (nc = R) | .numCommon 1x1
%              .days     1xD double   default day list (numeric integers);
%                                     required unless every system has .days
%              .systems  1xM struct   one entry per observing system, fields:
%                REQUIRED
%                  .buildDay          fh        your day builder (or .loader)
%                  .observationTypes  1xT cellstr/string
%                  .covariance        struct    .(obsType) = covariance spec
%                  .systemId          1x1 numeric | 1 char, unique
%                  ONE of .localsPerArc | .numLocal | .localParams  (nl; 0 valid)
%                OPTIONAL (default)
%                  .arcs (1) .days ([]) .name ('') .cache (true)
%                  .epochCovariance ([]) .processingOptions (struct())
%                  .plotPostSolve ([]) .plotPostModel ([])
%                ENGINE ESCAPE HATCH
%                  .loader   fh   a raw 4-arg closure instead of .buildDay
%                                 (all systems must then use .loader)
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
%     (tests/test_fit.m, reference). It is the SAME code path; the 1e-12
%     rather than isequal only guards the multithreaded-BLAS reduction-order
%     nondeterminism recorded in HANDOVER §10.8 (WP-47 r1e).
%
%   ERRORS
%     Your declaration is wrong (raised BEFORE any data is read):
%       lsqvce:fit:unknownField    a field name the contract does not know
%                                  (names the field; typos are errors on purpose)
%       lsqvce:fit:noSystems       data.systems missing, not a struct, or empty
%       lsqvce:fit:missingField    a required system field is absent
%       lsqvce:fit:localCount      not exactly one of localsPerArc|numLocal|
%                                  localParams
%       lsqvce:fit:commonCount     data.common has not exactly one of
%                                  commonParams|numCommon
%       lsqvce:fit:noDays          no day list anywhere
%       lsqvce:fit:duplicateSystem two systems share a systemId
%       lsqvce:fit:builderKind     buildDay/loader is not a function handle,
%                                  both given, or systems mix the two kinds
%       lsqvce:fit:<makeLoader id> a makeLoader construction error, relabelled
%                                  (e.g. noCov, badArcs, localPlan); the
%                                  original is attached as the cause
%     Your choice is not available here:
%       lsqvce:fit:notYet          Model "ellmer" (WP-56c) or Robust ~= "none"
%                                  on the STREAMING kind -- use the engine
%                                  directly
%       lsqvce:fit:modelKind       Model other than "gls" on kinds B / C
%       lsqvce:fit:badMatrix       kind B shapes disagree (A m x u, l m x 1,
%                                  p m x 1)
%       lsqvce:fit:missingField    kind C lacks numUnknowns / numBatches
%     Your buildDay is wrong (raised WHILE READING, by the engine, ids
%     unchanged): lsqvce:makeLoader:badDay / dayShape / badArc / badArcId /
%     buildDay, lsqvce:accum:arcBadNc / arcBadNl
%
%   EXAMPLE
%     data.common.numCommon = 6;   data.days = 1:3;
%     s.buildDay = @(d, ot) myBuildDay(d, ot);
%     s.observationTypes = {'pos','kbr'};
%     s.covariance.pos = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',0.7);
%     s.covariance.kbr = struct('type','diag','sigma',0.2);
%     s.systemId = 1;   s.numLocal = 4;   s.arcs = 2;
%     data.systems = s;
%     R = lsqvce.fit(data, lsqvce.Options(Model="vce"));
%     R.report
%
%   See also LSQVCE.OPTIONS, LSQVCE.RESULTS, LSQVCE.MAKELOADER,
%            LSQVCE.ADJUSTBLOCK, LSQVCE.ADJUSTVCE, LSQVCE.ADJUSTVCEJOINT,
%            LSQVCE.ADJUST, LSQVCE.ROBUST, LSQVCE.SOLVENONLINEAR
%
%   Updated 2026-08-31 (CEST, UTC+02:00), R4 WP-56b (Level L4): kinds B
%   (in-memory matrices) and C (batched callback, linear or nonlinear with
%   x0) through adjust / robust / solveNonlinear, returned as Results without
%   an accumulator. Model "ellmer" stays refused (WP-56c). With assistance
%   from Claude (Anthropic), model Claude Fable 5.
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
    if ~strcmp(opts.Model, "gls") && ~strcmp(opts.Model, "vce") && ~strcmp(opts.Model, "pooled")
        error('lsqvce:fit:notYet', ...
            ['Model "%s" is not wired through fit yet (WP-56c); call ' ...
             'lsqvce.adjustVCEjoint directly for the Ellmer joint loop.'], opts.Model);
    end
    if ~strcmp(opts.Robust, "none")
        error('lsqvce:fit:notYet', ...
            ['Robust = "%s" is the in-memory / batched kinds'' path (lsqvce.robust); ' ...
             'the streaming kind screens outliers through its diagnostics (report()).'], opts.Robust);
    end

    spec = validateData(data);
    [systems, days] = buildSystems(spec);

    switch char(opts.Model)
        case 'gls'
            a = opts.engineArgs('adjustBlock');
            R = lsqvce.adjustBlock(systems, days, a.nv{:});
        otherwise                                   % 'vce' | 'pooled'
            a = opts.engineArgs('adjustVCE');
            P  = lsqvce.adjustVCE(systems, days, a.nv{:});
            nv = struct2nv(a.shape);   P.stochasticModel(nv{:});
            nv = struct2nv(a.scale);   P.varianceComponents(nv{:});
            nv = struct2nv(a.outer);   P.iteration(nv{:});
            if strcmp(a.mode, 'joint')
                nv = struct2nv(a.joint);   P.joint(nv{:});
            end
            P.componentsPerEpoch(a.axes);
            R = P.run();
    end
end

%==========================================================================
function spec = validateData(data)
%VALIDATEDATA  Enforce DATA_CONTRACT.md §2 on the raw struct; fail fast.
%   OUTPUT spec 1x1 struct  .common (struct), .days (double), .systems (1xM)
    topAllowed = {'common', 'days', 'systems'};
    rejectUnknown(fieldnames(data), topAllowed, 'data');

    if ~isfield(data, 'systems') || ~isstruct(data.systems) || isempty(data.systems)
        error('lsqvce:fit:noSystems', ...
            'data.systems must be a non-empty 1xM struct array (one entry per observing system).');
    end
    sys = data.systems(:).';

    % ---- common block: nc stated exactly once
    common = struct('commonParams', {{}}, 'numCommon', []);
    if isfield(data, 'common')
        c = data.common;
        if ~isstruct(c) || ~isscalar(c)
            error('lsqvce:fit:commonCount', 'data.common must be a scalar struct.');
        end
        rejectUnknown(fieldnames(c), {'commonParams', 'numCommon'}, 'data.common');
        hasP = isfield(c, 'commonParams') && ~isempty(c.commonParams);
        hasN = isfield(c, 'numCommon')    && ~isempty(c.numCommon);
        if hasP == hasN
            error('lsqvce:fit:commonCount', ...
                'data.common must carry exactly ONE of commonParams (Rx14 cell) or numCommon (1x1).');
        end
        if hasP, common.commonParams = c.commonParams; end
        if hasN, common.numCommon    = c.numCommon;    end
    else
        error('lsqvce:fit:commonCount', ...
            'data.common is missing: state the common unknowns ONCE as commonParams or numCommon.');
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
    sysRequired = {'observationTypes', 'covariance', 'systemId'};
    sysAllowed  = [sysRequired, {'buildDay', 'loader', 'localsPerArc', 'numLocal', ...
                   'localParams', 'arcs', 'days', 'name', 'cache', 'epochCovariance', ...
                   'processingOptions', 'plotPostSolve', 'plotPostModel'}];
    rejectUnknown(fieldnames(sys), sysAllowed, 'data.systems');
    ids  = cell(1, numel(sys));
    kind = zeros(1, numel(sys));                    % 1 buildDay | 2 loader
    for k = 1:numel(sys)
        s = sys(k);
        for f = sysRequired
            if ~isfield(s, f{1}) || isempty(s.(f{1}))
                error('lsqvce:fit:missingField', ...
                    'data.systems(%d) needs the field ''%s''.', k, f{1});
            end
        end
        hasB = isfield(s, 'buildDay') && ~isempty(s.buildDay);
        hasL = isfield(s, 'loader')   && ~isempty(s.loader);
        if hasB == hasL
            error('lsqvce:fit:builderKind', ...
                'data.systems(%d) must carry exactly ONE of buildDay or loader.', k);
        end
        if (hasB && ~isa(s.buildDay, 'function_handle')) || (hasL && ~isa(s.loader, 'function_handle'))
            error('lsqvce:fit:builderKind', ...
                'data.systems(%d): buildDay / loader must be a function handle.', k);
        end
        kind(k) = 1 + hasL;
        nLocalWays = (isfield(s, 'localsPerArc') && ~isempty(s.localsPerArc)) + ...
                     (isfield(s, 'numLocal')     && ~isempty(s.numLocal)) + ...
                     (isfield(s, 'localParams')  && ~isempty(s.localParams));
        if nLocalWays ~= 1
            error('lsqvce:fit:localCount', ...
                ['data.systems(%d) must declare its local unknowns in exactly ONE way: ' ...
                 'localsPerArc, numLocal (0 is valid), or localParams.'], k);
        end
        id = s.systemId;
        if ~((isnumeric(id) && isscalar(id)) || (ischar(id) && isscalar(id)) || ...
             (isstring(id) && isscalar(id) && strlength(id) == 1))
            error('lsqvce:fit:missingField', ...
                'data.systems(%d).systemId must be a numeric scalar or a single character.', k);
        end
        ids{k} = char(string(id));
        if isempty(days) && ~(isfield(s, 'days') && ~isempty(s.days))
            error('lsqvce:fit:noDays', ...
                'no day list: set data.days or data.systems(%d).days.', k);
        end
    end
    if numel(unique(ids)) < numel(ids)
        error('lsqvce:fit:duplicateSystem', 'systemId values must be unique across data.systems.');
    end
    if numel(unique(kind)) > 1
        error('lsqvce:fit:builderKind', ...
            'all systems must use buildDay, or all must use loader -- not a mixture.');
    end
    spec = struct('common', common, 'days', days, 'systems', sys, 'kind', kind(1));
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
function [systems, days] = buildSystems(spec)
%BUILDSYSTEMS  One Loader per system through makeLoader (kind 1), or the
%   struct-array raw-closure form adjustBlock/adjustVCE accept (kind 2).
    sys  = spec.systems;
    days = spec.days;
    nc   = [];
    if ~isempty(spec.common.numCommon), nc = spec.common.numCommon; end
    if spec.kind == 2
        systems = struct('sysID', {}, 'nl', {}, 'nc', {}, 'loader', {}, 'days', {});
        for k = 1:numel(sys)
            s = sys(k);
            if isempty(nc)
                nc = size(spec.common.commonParams, 1);
            end
            nl = [];
            if isfield(s, 'numLocal') && ~isempty(s.numLocal), nl = s.numLocal; end
            if isempty(nl)
                error('lsqvce:fit:localCount', ...
                    'data.systems(%d): a raw loader needs numLocal (the plan/table forms belong to buildDay).', k);
            end
            d = days;
            if isfield(s, 'days') && ~isempty(s.days), d = double(s.days(:).'); end
            systems(k) = struct('sysID', s.systemId, 'nl', nl, 'nc', nc, ...
                                'loader', s.loader, 'days', d);
        end
        return
    end
    % ---- the ONE mapping table: lowerCamelCase field -> makeLoader NV name
    map = {'observationTypes', 'ObservationTypes'; 'covariance', 'Covariance'; ...
           'arcs', 'Arcs'; 'localsPerArc', 'LocalsPerArc'; 'systemId', 'SystemId'; ...
           'numLocal', 'NumLocal'; 'localParams', 'LocalParams'; 'days', 'Days'; ...
           'name', 'Name'; 'cache', 'Cache'; 'epochCovariance', 'EpochCovariance'; ...
           'processingOptions', 'ProcessingOptions'; 'plotPostSolve', 'PlotPostSolve'; ...
           'plotPostModel', 'PlotPostModel'};
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
        if ~(isfield(s, 'days') && ~isempty(s.days)) && ~isempty(days)
            nv = [nv, {'Days', days}];                               %#ok<AGROW>
        end
        if ~isempty(spec.common.commonParams)
            nv = [nv, {'CommonParams', spec.common.commonParams}];   %#ok<AGROW>
        else
            nv = [nv, {'NumCommon', nc}];                            %#ok<AGROW>
        end
        try
            systems{k} = lsqvce.makeLoader(s.buildDay, nv{:});
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

function nv = struct2nv(s)
%STRUCT2NV  Struct -> Name-Value cell, field order preserved.
    f  = fieldnames(s);
    nv = cell(1, 2*numel(f));
    for k = 1:numel(f)
        nv{2*k-1} = f{k};
        nv{2*k}   = s.(f{k});
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
    rejectUnknown(fieldnames(data), {'batchFcn', 'numUnknowns', 'numBatches', 'x0'}, 'data');
    for f = {'batchFcn', 'numUnknowns', 'numBatches'}
        if ~isfield(data, f{1}) || isempty(data.(f{1}))
            error('lsqvce:fit:missingField', 'data lacks the required field %s.', f{1});
        end
    end
    x0 = [];
    if isfield(data, 'x0') && ~isempty(data.x0), x0 = data.x0(:); end
    R = runBatch(data.numUnknowns, data.batchFcn, data.numBatches, x0, opts, 'batch');
end

function R = runBatch(u, batchFcn, nBatch, x0, opts, kind)
%RUNBATCH  The ONE dispatch for kinds B / C: adjust | robust | solveNonlinear.
    if ~strcmp(opts.Model, "gls")
        error('lsqvce:fit:modelKind', ...
            ['Model "%s" needs the streaming kind (variance components are per ' ...
             'group of a system); the %s kind is Model "gls".'], opts.Model, kind);
    end
    if ~isempty(x0)
        a = opts.engineArgs("solveNonlinear");
        [sol, D, ~, IT] = lsqvce.solveNonlinear(u, batchFcn, nBatch, x0, a.nv{:});
        source = IT;
    elseif strcmp(opts.Robust, "irls")
        a = opts.engineArgs("robust");
        [sol, D, ROB] = lsqvce.robust(u, batchFcn, nBatch, a.nv{:});
        source = ROB;
    else
        a = opts.engineArgs("adjust");
        [sol, D, V] = lsqvce.adjust(u, batchFcn, nBatch, a.nv{:});
        source = V;
    end
    % The Results contract is the STREAMING one (BlockAccumulator.solve):
    % sol.sigma0 = a-posteriori std of unit weight, sol.Cxx = covariance
    % scaled by it. The batch drivers return Accumulator.solve's struct,
    % where sigma0 is the A-PRIORI Sigma0 and there is no Cxx (R, QxxCorr,
    % s0sq instead) -- so the bridge derives both exactly as lsqvce.diagnose
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
    R = lsqvce.Results(sol, D, [], source, batchFcn, nBatch, [], runMeta(struct(), resolved));
end
