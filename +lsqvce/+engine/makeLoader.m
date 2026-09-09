function L = makeLoader(buildDay, opt)
%LSQVCE.ENGINE.MAKELOADER  Turn a day-builder into an accumulator-ready Loader.
%   L = lsqvce.engine.makeLoader(buildDay, Name, Value, ...)
%
%   You write the ONE part that is genuinely yours -- reading your data for a
%   day and returning its design matrices, observations and their error
%   information -- and makeLoader supplies the protocol the accumulator
%   expects: the 'meta' discovery call, the per-arc slicing, the day cache,
%   and the covariance / whitener of every arc through LSQVCE.ENGINE.COVMODEL.
%   The result L plugs straight into LSQVCE.ENGINE.ADJUSTBLOCK,
%   LSQVCE.ENGINE.ADJUSTROBUST or LSQVCE.ENGINE.BLOCKACCUMULATOR/addSystem;
%   LSQVCE.FIT builds one per declared system (LSQVCE.ENGINE.DECLARATION).
%
%   THE PART YOU WRITE
%     buildDay(day, obsType) -> struct                            (classic)
%     buildDay(day, obsType, common, local, procOpt) -> struct    (parameter-driven)
%     buildDay(day, obsType, common, local, procOpt, info) -> struct   (with the layout)
%       THE CELLS ARE THE ARCS (v6, R6 WP-102): every field is a 1 x nArc cell,
%       one cell per arc of the day, cell a = arc a. A day with one arc returns
%       {A}; a day without data returns empty cells (nArc = 0).
%         .A       {1 x nArc}  each m_a x nc     common design matrix of arc a
%         .B       {1 x nArc}  each m_a x nl     local design matrix (may be sparse)
%         .l       {1 x nArc}  each m_a x 1      observations (observed minus computed)
%         .t       {1 x nArc}  each m_a x 1      epoch of every row [s]
%         .Cep     {1 x nArc}  OPTIONAL: THE ERROR INFORMATION of every arc (R6
%                            WP-90) -- per arc nEpoch_a x 1 (the per-epoch standard
%                            deviation of a scalar stream), nEpoch_a x d (per-axis
%                            standard deviations) or nEpoch_a x d x d (the full
%                            per-epoch covariance). nEpoch_a * d must be m_a; the
%                            axis count is DERIVED from it and is one for the day.
%                            The toolbox composes it with the declared temporal law
%                            per arc: C = D (R_t kron I_d) D.
%         .nAxis   1x1       the axis count of a multi-axis stream WITHOUT Cep
%                            (rows epoch-major, m_a = nEpoch_a * nAxis); derived
%                            from Cep when that is present and must then agree
%       An EMPTY cell (0 rows in A, B, l, t) is a DATA GAP: that arc of the
%       schedule is skipped by every engine. Any other field is an error --
%       in particular the v5 fields arcIdx / arcId, whose job the cells do.
%     It must be DETERMINISTIC: the accumulator reads each day a second time
%     for the residuals and must get identical data back.
%
%     buildDay is called with as many of the SIX arguments as it DECLARES, so
%     a @(day,obsType) builder keeps working unchanged; a builder that also
%     takes (common, local, procOpt) receives the parametrisation tables and
%     processing options passed to makeLoader; the 6th argument is the day's
%     expansion slice (offsets of every block of A and B, the rows present on
%     the day; the Manual, chapter 3 section 6.3). A varargin builder receives
%     all six. common / local have ONE ROW per template and C columns of YOUR
%     convention; nc and nl are DERIVED (see 'CommonParams' / 'LocalParams').
%
%   INPUT
%     buildDay   1x1 function_handle   the day builder above. Required.
%     Name, Value options ([] or omitted -> defaults; canonical UpperCamelCase,
%     lowerCamelCase casings keep matching -- MATLAB NV semantics):
%       'ObservationTypes' 1xT cellstr/string  the streams, e.g. {'pos','kbr'}. Required.
%       'Covariance'   1x1 struct   .(obsType) = the TEMPORAL law of that stream, a
%                                   spec for LSQVCE.ENGINE.COVMODEL (diag | arma | psd |
%                                   acf | var | custom; sigma = its LEVEL). One field per
%                                   obsType. Required. The per-epoch error information is
%                                   NOT declared here -- the day returns it (.Cep).
%       'Arcs'         var          the arc SCHEDULE -- which arc ids each day carries.
%                                   OPTIONAL since v6 (R6 WP-102): omitted ([]), the arcs
%                                   of a day are the cells buildDay returns, ids 1:nArc,
%                                   learned by reading the day when the engine's 'meta'
%                                   probe asks for it (the day is read right after, so
%                                   the cache makes this free). Declared, it is resolved
%                                   WITHOUT reading data and cell a is arc ids(a):
%                                     1x1 integer   -> arcs 1:Arcs for every obsType
%                                     1xA vector    -> those arc ids for every obsType
%                                     1x1 struct    -> .(obsType) = any form on this list
%                                     function_handle -> @(day) or @(day,obsType), that
%                                                      day's arc ids (a ragged schedule)
%                                   and the day must return exactly as many cells as the
%                                   schedule names (lsqvce:makeLoader:arcCount); an empty
%                                   cell is the data gap of that id. REQUIRED by 'arc'
%                                   template rows and by 'LocalsPerArc' (nl must be known
%                                   before any day is read).                (default [])
%       'Scoped'       1x1 logical  v4 contract (L11 WP-75): LocalParams is SCOPED --
%                                   column 1 the token 'arc' | 'day' | 'system', rows are
%                                   templates expanded over the calendar (nl DERIVED;
%                                   layout per day [arc blocks, sorted ids | 'day'
%                                   block], trailing 'system' block; a day whose schedule
%                                   has no arcs gets no blocks). Default false: free tables.
%       'CommonPlan'   1x1 struct   v4: the COMMON scope plan LSQVCE.ENGINE.DECLARATION
%                                   builds over the GLOBAL calendar (.nc .commonCount .days
%                                   .dayOffset .dayCount .rows); sets nc and the common
%                                   slice of buildDay's 6th argument.
%       'ValidityColumns' [c1 c2]   R5-D8 round 3: the columns of LocalParams holding a
%                                   row's validity ([MJD sec] or fractional MJD; [] =
%                                   always); a 'day' / 'arc' row is instantiated only on
%                                   the days it overlaps.                (default [])
%       'Priors'       1xG struct   R5-D8 rounds 5/6: priors on this system's LOCAL
%                                   parameters, one struct per group: .group (a regexp on
%                                   the name column of LocalParams; L built by the loader)
%                                   OR .L (q x nl in the columns of B); .x0 .sigma .from
%                                   .name .estimate. Each group rides as its own stream
%                                   'prior_<name>' (pseudo-observations, unit law).
%       'NameColumn'   1x1          the name column of LocalParams        (default 2)
%       'LocalsPerArc' var          engine tier (WP-52): a per-arc local count, 1x1 or
%                                   @(day,arcId), when arcs carry DIFFERENT counts; nl
%                                   derived, offsets handed to buildDay. Needs 'Days' and
%                                   a declared 'Arcs' schedule.               (default [])
%       'SystemId'     1x1          system identifier carried on L (default 1). NUMERIC or a
%                                   SINGLE char -- the accumulator indexes systems by
%                                   [store.sysID]; a longer label goes in 'Name'.
%       'NumLocal'     1x1          local count when no LocalParams is given (0 is valid:
%                                   no local block, B is m x 0). Given with a table, it
%                                   must agree.                          (default [])
%       'NumCommon'    1x1          common count when no CommonParams / CommonPlan is
%                                   given; nc is SHARED by every system. (default [])
%       'CommonParams' RxC cell     the common templates, piped to buildDay (3rd arg);
%                                   nc derived from the rows (v3) or the plan (v4)
%       'LocalParams'  RxC cell     the local templates, piped to buildDay (4th arg)
%       'ProcessingOptions' struct  free-form options piped to buildDay (5th arg)
%       'Days'         1xD          default day list carried on L        (default [])
%       'Name'         char/string  display name                         (default '')
%       'Cache'        1x1 logical  keep built days resident             (default true)
%       'CacheDays'    1x1          resident built DAYS (R5 M4 4d): 1 = the current day
%                                   with every one of its streams (R6 WP-102: the slot
%                                   count is CacheDays x streams, the key is day|stream),
%                                   Inf = the calendar (one build per run) (default 1)
%
%   RUNTIME COVARIANCE OVERRIDE (opt.covOverride)
%     Every batch merges a runtime override into the base spec, so an outer
%     loop can reshape the law per pass WITHOUT rebuilding the loader. Two
%     shapes, later wins: opt.covOverride.(obsType) (flat) and
%     opt.covOverride.bySys.(satKey).(obsType) (per system, satKey =
%     matlab.lang.makeValidName of the sysID) -- the shape
%     LSQVCE.ENGINE.JOINTMODELOPT emits. The override carries the temporal law
%     and a level; the day's Cep composes with it unchanged.
%
%   OUTPUT
%     L  1x1 struct  the Loader:
%        .closure  @(day,obsType,arcID,opt)  the 4-arg callback the accumulator drives
%                                            (answers the 'meta' probe too)
%        .meta     @([day])  -> .obsTypes, .arcs (the arc ids per stream on that day)
%        .batch    @(day,obsType,arcID,opt)  one arc's batch (.A .B .l .t .Cfcn .Wfcn
%                                            .filterMemory; .Cep and .nAxis when the
%                                            day carries error information); [] for a
%                                            data gap (an empty cell)
%        .day      @(day,obsType)  the cached, validated whole day
%        .isLoader .sysID .nc .nl .obsTypes .cov .name .days .plan .priors
%        .common .local .procOpt   the tables and options piped to buildDay
%
%   THE LOADER CONTRACT (what the engines require)
%     The ONLY hard requirement is the 4-arg data callback that also answers
%     the 'meta' / 'day' sentinels -- .closure. Every other field is
%     OPTIONAL and read defensively, so a hand-written raw closure (e.g.
%     LSQVCE.ENGINE.SYNTHLOADER) is a first-class loader: pass it to
%     acc.addSystem or to adjustBlock's struct form without wrapping it here.
%
%   EXAMPLE
%     cov.pos = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',0.7);
%     L = lsqvce.engine.makeLoader(@(d,ot) myBuildDay(d,ot), ...
%             'ObservationTypes',{'pos'}, 'Covariance',cov, 'Arcs',2, 'SystemId','C', 'NumLocal',4);
%     R = lsqvce.engine.adjustBlock({L}, 1:3);
%
%   ACCURACY
%     No arithmetic of its own beyond slicing: the whitening it hands out is
%     LSQVCE.ENGINE.COVMODEL's. The composition of an arc's Cep with the
%     temporal law is pinned through THIS unit against the dense reference
%     D (R_t kron I_d) D at 1e-8 relative on the whitened normals
%     (test_streamUtils, R6 WP-90); the cell contract is pinned by parity
%     with the v5 row-index form to 1e-12 on x / sigma0 / Cxx (test_fit,
%     R6 WP-102), and an arc's Cep is its own cell (finding E-01 cannot recur).
%
%   ERRORS
%     Raised while BUILDING the loader (your declaration is wrong):
%       lsqvce:makeLoader:noObsTypes   'ObservationTypes' missing or empty
%       lsqvce:makeLoader:noCov        'Covariance' missing, or no entry for a stream
%       lsqvce:makeLoader:noArcs       a per-stream 'Arcs' struct has no field for a stream
%       lsqvce:makeLoader:badArcs      'Arcs' is a string (the v5 'fromDay' -- omit 'Arcs'
%                                      instead), or a resolver returned something that is
%                                      not unique positive integers (WP-52)
%       lsqvce:makeLoader:localPlan    'LocalsPerArc' / scoped rows without 'Days', 'arc'
%                                      rows or 'LocalsPerArc' without a declared 'Arcs'
%                                      schedule, or a bad count
%       lsqvce:makeLoader:paramType    CommonParams / LocalParams is not a 2-D cell
%       lsqvce:makeLoader:badScope     a scoped LocalParams row whose column 1 is not
%                                      'arc' | 'day' | 'system'
%       lsqvce:makeLoader:ncConflict   NumCommon disagrees with the table / plan
%       lsqvce:makeLoader:nlConflict   NumLocal disagrees with the table / plan
%       lsqvce:makeLoader:badValidity  ValidityColumns out of range or a bad entry
%       lsqvce:makeLoader:badPrior     a Priors group is malformed (fields, width,
%                                      pattern, a row coupling two days)
%       lsqvce:makeLoader:priorType    a stream name collides with a prior group's
%
%     Raised while READING a day (your buildDay is wrong):
%       lsqvce:makeLoader:buildDay     your buildDay itself threw; the original error
%                                      is attached as the cause
%       lsqvce:makeLoader:badDay       a required field (.A .B .l .t) is missing or is not
%                                      a 1 x nArc cell, or the day carries a field that is
%                                      not .A .B .l .t .Cep .nAxis (the v5 arcIdx / arcId
%                                      among them)
%       lsqvce:makeLoader:dayShape     the cells do not agree in number, or an arc's A, B,
%                                      l, t do not share a row count
%       lsqvce:makeLoader:arcCount     a declared 'Arcs' schedule names n arcs for the day
%                                      and buildDay returned a different number of cells
%       lsqvce:makeLoader:emptyDay     without a declared schedule, a calendar day that
%                                      owns a 'day' block returned no arc in any stream
%                                      (take the day off the calendar, or declare 'Arcs')
%       lsqvce:makeLoader:badArc       an arc of a multi-axis day does not hold whole epochs
%       lsqvce:makeLoader:badDayCep    .Cep is not a 1 x nArc cell, or an arc's entry is not
%                                      nEpoch x 1 / nEpoch x d / nEpoch x d x d, or
%                                      nEpoch * d is not the arc's row count (R6 WP-90)
%       lsqvce:makeLoader:axisConflict .nAxis disagrees with .Cep or with a 'd' on the
%                                      spec, or does not divide the rows (R6 WP-90)
%       lsqvce:makeLoader:cepTwice     the spec carries a Cep AND the day returns one
%                                      (R6 WP-90: one fact, one authority)
%       lsqvce:makeLoader:noCep        [warning, once per run] a multi-axis day without
%                                      Cep: composed with D = I (R6 WP-90)
%       lsqvce:makeLoader:rowWeightSize  (R5 M2) the robust weight vector of an arc does
%                                      not match the arc's row count
%       lsqvce:covModel:*              the temporal law is malformed (see covModel)
%
%   See also LSQVCE.ENGINE.ADJUSTBLOCK, LSQVCE.ENGINE.COVMODEL, LSQVCE.ENGINE.DECLARATION,
%            LSQVCE.ENGINE.BLOCKACCUMULATOR, LSQVCE.ENGINE.ADJUSTVCEJOINT, LSQVCE.ENGINE.SYNTHLOADER
%
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-102 (v6): THE CELLS ARE THE
%   ARCS -- A, B, l, t (and Cep) come back as 1 x nArc cells, arcIdx / arcId
%   are gone, an empty cell is a data gap, 'Arcs' is optional (the cells of a
%   day are its schedule; a declared schedule must match their number), the
%   string form 'fromDay' is retired, and the day cache holds CacheDays x
%   streams slots. With assistance from Claude (Anthropic), model Claude
%   Fable 5.1.
%   Updated 2026-09-06 (CEST, UTC+02:00), R6 WP-90: the error information
%   travels with the day (.Cep for every stream, nAxis derived, validated at
%   the day; the Cep slice by arc POSITION -- finding E-01); the
%   'EpochCovariance' hook and the dead plot-hook helpers are gone; the help
%   rewritten on the v5 surface. With assistance from Claude (Anthropic),
%   model Claude Fable 5.1.
%   Updated 2026-08-31 (CEST, UTC+02:00), v4 scoped tables (Level L11,
%   WP-75): 'Scoped' + 'CommonPlan' options; buildScopePlan generalizes the
%   WP-52 per-arc plan to 'arc'/'day'/'system' template rows (nl derived).
%   With assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-21 (CEST, UTC+02:00), WP-52: ragged schedules are
%   NATIVE (per-day arc resolution, arcs addressed by .arcId, 'localsPerArc'
%   owning the per-arc bookkeeping). Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (c) c1 (F-07): typed
%   `arguments` signature, canonical UpperCamelCase option names. Claude Fable 5.
%   Created 2026-07-16 12:20 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
arguments
    buildDay
    opt.ObservationTypes = {}
    opt.Covariance (1,1) struct = struct()
    opt.Arcs = []
    opt.LocalsPerArc = []
    opt.SystemId = 1
    opt.NumLocal double {mustBeScalarOrEmpty} = []
    opt.NumCommon double {mustBeScalarOrEmpty} = []
    opt.Days double = []
    opt.Name (1,1) string = ""
    opt.Cache (1,1) logical = true
    opt.ValidityColumns double = []   % R5-D8 round 3: [cStart cEnd] columns of localParams holding the
                                      % validity ([MJD sec] pairs or fractional MJD; [] = always): a
                                      % 'day'/'arc' row is instantiated only on days its validity overlaps
    opt.Priors = []               % R5-D8 round 5/6: struct array, the priors on this system's LOCAL
                                  % parameters: per GROUP either .group (a regexp matched against the
                                  % NAME column of localParams -> L built by the loader, one row per
                                  % matched template row per instance) or .L (q x nl in the columns
                                  % of B, used as is); .x0, .sigma, .from, .name, .estimate; each
                                  % group rides as its own observation type 'prior_<name>'
    opt.NameColumn (1,1) double {mustBeInteger, mustBePositive} = 2   % the name column of localParams
    opt.CacheDays (1,1) double {mustBePositive} = 1   % R5 M4 4d: built days kept resident;
                                                     % Inf = every day of the run (the
                                                     % user's memory decision, see fit)
    opt.CommonParams cell = {}
    opt.LocalParams cell = {}
    opt.ProcessingOptions (1,1) struct = struct()
    opt.Scoped (1,1) logical = false
    opt.CommonPlan = []
end

    narginchk(1, inf);
    if ~isa(buildDay, 'function_handle')
        error('lsqvce:makeLoader:buildDay', 'buildDay must be a function handle.');
    end
    % WP-39 stage (c) c1: typed `arguments` signature (canonical
    % UpperCamelCase; legacy casings keep matching -- MATLAB NV
    % semantics; enum VALUES case-sensitive lowercase).
    o = struct('observationTypes', {opt.ObservationTypes}, 'covariance', opt.Covariance, 'arcs', opt.Arcs, 'localsPerArc', opt.LocalsPerArc, 'systemId', opt.SystemId, 'numLocal', opt.NumLocal, 'numCommon', opt.NumCommon, 'days', opt.Days, 'name', char(opt.Name), 'cache', opt.Cache, 'cacheDays', opt.CacheDays, 'priors', opt.Priors, 'nameColumn', opt.NameColumn, 'validityColumns', opt.ValidityColumns, 'commonParams', {opt.CommonParams}, 'localParams', {opt.LocalParams}, 'processingOptions', opt.ProcessingOptions, 'scoped', opt.Scoped, 'commonPlan', opt.CommonPlan);

    % ---- map long option names onto the short internal fields -------------
    % Canonical long names only. The resolved value is written onto the SHORT
    % field name so the rest of makeLoader (which reads o.obsTypes/o.cov/
    % o.sysID/... and the local vars) stays unchanged.
    o.obsTypes = o.observationTypes;
    o.cov      = o.covariance;
    o.sysID    = o.systemId;
    o.nl       = o.numLocal;
    o.nc       = o.numCommon;
    o.common   = o.commonParams;
    o.local    = o.localParams;
    o.procOpt  = o.processingOptions;

    obsTypes = cellstr(o.obsTypes);
    if isempty(obsTypes)
        error('lsqvce:makeLoader:noObsTypes', 'Provide ''observationTypes''.');
    end
    for k = 1:numel(obsTypes)
        if ~isfield(o.cov, obsTypes{k})
            error('lsqvce:makeLoader:noCov', ...
                'covariance has no field for observationType ''%s''.', obsTypes{k});
        end
    end
    % R6 WP-102 (v6): the schedule is OPTIONAL -- omitted, the cells a day
    % returns ARE its arcs (that is what the v5 string form 'fromDay' meant,
    % so the string form is retired under the no-aliases rule).
    if ischar(o.arcs) || isstring(o.arcs)
        error('lsqvce:makeLoader:badArcs', ...
              ['''Arcs'' = ''%s'': the string form is retired (v6). Omit ''Arcs'' -- the ' ...
               'cells buildDay returns are the arcs of the day -- or declare a count, an ' ...
               'id vector, a per-stream struct or a resolver @(day,obsType).'], char(o.arcs));
    end

    % ---- parametrisation arrays: validate, derive nc/nl -------------------
    % Each is a cell array with ONE ROW per unknown and 14 columns of meaning
    % (the meaning is the caller's convention; makeLoader only checks shape and
    % pipes the arrays through to buildDay). Empty ({}) means "not used" -- the
    % builder then defines the design itself and nc/nl come from 'nl'/adjustBlock.
    checkParamArray(o.common, 'common');
    checkParamArray(o.local,  'local');

    % nc DERIVED from common rows (the common field is SHARED across systems, so
    % every loader must pass the SAME 'common'); nl DERIVED from local rows.
    ncDerived = size(o.common, 1);        % 0 if common == {}
    nlDerived = size(o.local,  1);        % 0 if local  == {}

    ncFinal = o.nc;
    if ~isempty(o.common) && ~(o.scoped && ~isempty(o.commonPlan))
        % commonParams given -> it SETS nc (v3). Under a SCOPED common plan
        % (v4) the rows are TEMPLATES: nc comes from the plan, below.
        if ~isempty(o.nc) && o.nc ~= ncDerived
            error('lsqvce:makeLoader:ncConflict', ...
                ['''numCommon'' = %d but the ''commonParams'' parametrisation ' ...
                 'has %d rows. Omit ''numCommon'' (it is derived from ' ...
                 '''commonParams'') or make them agree.'], ...
                 o.nc, ncDerived);
        end
        ncFinal = ncDerived;
    end

    nlFinal = o.nl;
    if ~isempty(o.local) && ~o.scoped
        % local given -> it SETS nl (v3). Under 'Scoped' the rows are
        % TEMPLATES: nl comes from the scope plan, below.
        if ~isempty(o.nl) && o.nl ~= nlDerived
            error('lsqvce:makeLoader:nlConflict', ...
                ['''numLocal'' = %d but the ''localParams'' parametrisation has %d rows. ' ...
                 'Omit ''numLocal'' (it is derived from ''localParams'') or make them agree.'], ...
                 o.nl, nlDerived);
        end
        nlFinal = nlDerived;
    end

    % WP-52: the per-arc local plan. Declaring 'localsPerArc' makes makeLoader
    % OWN the bookkeeping the user used to hand-roll: it walks the day list
    % once, assigns each (day, arc) its global ordinal and column offset,
    % DERIVES nl by summing, and hands each day its slice as buildDay's 6th
    % argument. nl can then never disagree with the schedule, because it is no
    % longer stated twice.
    plan = [];
    if o.scoped && ~isempty(o.local)
        % L11 WP-75 (v4): the SCOPED local table -- column 1 is the scope
        % token ('arc' | 'day' | 'system'), rows are templates expanded over
        % the calendar. Generalizes the WP-52 per-arc plan; nl is DERIVED.
        if ~isempty(o.localsPerArc)
            error('lsqvce:makeLoader:localPlan', ...
                '''Scoped'' local tables and ''localsPerArc'' are one mechanism twice: declare the ''arc'' rows in the table.');
        end
        plan = buildScopePlan(o.local, o.arcs, obsTypes, o.days, o.validityColumns);
        if ~isempty(nlFinal) && nlFinal ~= plan.nl
            error('lsqvce:makeLoader:nlConflict', ...
                '''numLocal'' = %d but the scoped plan derives %d. Omit ''numLocal''.', nlFinal, plan.nl);
        end
        nlFinal = plan.nl;
    end
    if ~isempty(o.localsPerArc)
        plan = buildLocalPlan(o.localsPerArc, o.arcs, obsTypes, o.days);
        if ~isempty(nlFinal) && nlFinal ~= plan.nl
            error('lsqvce:makeLoader:nlConflict', ...
                ['''numLocal'' = %d but the ''localsPerArc'' plan needs %d ' ...
                 '(%d arc(s) over %d day(s)). Omit ''numLocal'' -- it is ' ...
                 'DERIVED from the plan.'], nlFinal, plan.nl, ...
                 plan.nArcTot, numel(o.days));
        end
        nlFinal = plan.nl;
    end

    cache = containers.Map('KeyType','char', 'ValueType','any');   % handle: mutable
    % R6 WP-102: the cache is keyed day|stream, so ONE resident day is one
    % slot per data stream (the meaning 'CacheDays' has always claimed; with
    % one slot the discovery probe of a 3-stream day evicted what the engine
    % was about to read). The prior streams added below never enter it.
    capacity = o.cache * o.cacheDays * numel(obsTypes);

    % R5-D8 round 2: priors on local template rows are pseudo-observations of
    % THIS system -- the observation type 'prior' (one arc per day, unit-whitened
    % rows, A-part zero, B-part the constrained local columns of that day). The
    % accumulator, the elimination, the sweeps and the redundancy count treat
    % them exactly as data (Koch: a prior is a group with a redundancy).
    L.priors = validatePriors(o.priors, plan, o.sysID, o.local, o.nameColumn);
    for g = 1:numel(L.priors)
        ot = L.priors(g).obsType;
        if any(strcmp(obsTypes, ot))
            error('lsqvce:makeLoader:priorType', ...
                  '''%s'' is reserved for prior group %d; rename your stream or the group.', ot, g);
        end
        obsTypes{end+1} = ot;                                     %#ok<AGROW>
        o.cov.(ot) = struct('type', 'diag', 'sigma', 1, 'estimate', L.priors(g).estimate, 'arcFactors', false);
    end
    o.obsTypes = obsTypes;
    L.isLoader = true;
    L.sysID    = o.sysID;
    L.nl       = nlFinal;
    L.obsTypes = obsTypes;
    L.cov      = o.cov;
    L.name     = char(o.name);
    % Parametrisation piped to buildDay + the DERIVED counts (nc is 0 when
    % 'common' is not used; it is the SHARED common count and every loader must
    % agree on it -- pass nc = L.nc to adjustBlock, or the explicit scalar).
    L.common   = o.common;
    L.local    = o.local;
    L.procOpt  = o.procOpt;
    if o.scoped && ~isempty(o.commonPlan)
        % L11 WP-75 (v4): fit builds the COMMON plan over the GLOBAL day
        % calendar (it spans systems) and passes it in; nc is the plan's.
        if ~isempty(ncFinal) && ncFinal ~= o.commonPlan.nc
            error('lsqvce:makeLoader:ncConflict', ...
                '''numCommon'' = %d but the common plan derives %d. Omit ''numCommon''.', ncFinal, o.commonPlan.nc);
        end
        ncFinal = o.commonPlan.nc;
    end
    L.nc       = ncFinal;
    L.plan     = plan;                                 % R5-D8 round 2: priorArc reads the layout
    L.commonPlan = o.commonPlan;                       % R6 WP-102: the emptyDay rule reads it
    if ~isempty(o.days), L.days = o.days; end

    % WP-52: L.meta() answers for a DAY. Called without one it uses the first
    % declared day (or 1), which is what the old day-blind accessor meant.
    % The day accessor is built FIRST and captured by value: L.meta needs it
    % without a declared schedule (R6 WP-102), and an anonymous function
    % cannot capture a struct field that does not exist yet.
    dayAccessor = @(day, obsType) dayOf(buildDay, cache, capacity, day, obsType, ...
                                        o.common, o.local, o.procOpt, plan, o.commonPlan);
    L.meta    = @(varargin) metaOf(obsTypes, o.arcs, ...
                                   metaDay(varargin, o.days), dayAccessor, L);
    L.localPlan = plan;               % WP-52: [] unless 'localsPerArc' declared
    L.day     = dayAccessor;
    L.batch   = @(day, obsType, arcID, opt) batchOf(L.day, o.arcs, o.cov, day, obsType, arcID, opt);
    L.closure = @(day, obsType, arcID, opt) dispatch(L, o.arcs, day, obsType, arcID, opt);
end

%==========================================================================
function checkParamArray(C, name)
%CHECKPARAMARRAY  Validate a parametrisation cell array: empty, or a 2-D cell
%   with one ROW per unknown (the column count and meaning are the caller's
%   convention -- v3 contract, L10 WP-72; the fixed 14 columns are retired).
%   Empty is allowed (means the builder defines the design itself).
    if isempty(C), return, end
    if ~iscell(C) || ~ismatrix(C)
        error('lsqvce:makeLoader:paramType', ...
            '''%s'' must be a 2-D cell array (one row per unknown).', name);
    end
end

function out = dispatch(L, arcs, day, obsType, arcID, opt)
%DISPATCH  Route a loader call: return metaOf(...) for the 'meta' probe, else
%   delegate to L.batch(day,obsType,arcID,opt). OUTPUT out = meta struct | batch struct.
    if nargin < 6 || isempty(opt), opt = struct(); end
    if (ischar(day) || isstring(day)) && strcmp(day, 'meta')
        % WP-52: the probe carries the DAY (every engine passes it), so a
        % ragged schedule is answered natively -- see metaOf/arcVec.
        out = metaOf(L.obsTypes, arcs, obsType, L.day, L);
        return
    end
    if startsWith(obsType, 'prior_')                  % R5-D8 round 5: a local prior group's pseudo-arc
        out = priorArc(L, day, obsType, arcID, opt, L.cov);
        return
    end
    out = L.batch(day, obsType, arcID, opt);
end

function m = metaOf(obsTypes, arcs, day, dayFcn, L)
%METAOF  The (obsTypes, arcs-per-type) probe answer for ONE day.
%   With a declared schedule the answer needs no data. Without one (R6
%   WP-102) every data stream of the day is read -- the engine reads them
%   right after, and the cache holds a whole day, so nothing is read twice --
%   and a day the plan gave a 'day' block that turns out empty in EVERY
%   stream is refused here, before any solve could go rank-deficient on it.
    m.obsTypes = obsTypes;
    m.arcs     = cell(1, numel(obsTypes));
    nData = 0;
    for k = 1:numel(obsTypes)
        if startsWith(obsTypes{k}, 'prior_')
            m.arcs{k} = priorArcsOfDay(L, day, obsTypes{k});   % R5-D8 round 5: one arc, or none
        else
            m.arcs{k} = arcVec(arcs, obsTypes{k}, day, dayFcn);
            nData = nData + numel(m.arcs{k});
        end
    end
    if isempty(arcs) && nData == 0 && dayOwnsBlock(L, day)
        error('lsqvce:makeLoader:emptyDay', ...
              ['day %s carries no arc in any stream of system %s, yet it owns a ''day'' ' ...
               'block (its ''day'' template rows or the common ''day'' rows are instantiated ' ...
               'on every calendar day when no schedule is declared). Take the day off this ' ...
               'system''s calendar (''days''), or declare ''Arcs'' with an empty list for it.'], ...
              string(day), char(string(L.sysID)));
    end
end

function tf = dayOwnsBlock(L, day)
%DAYOWNSBLOCK  Does the local plan or the common plan hold a 'day' block for DAY?
    tf = false;
    if ~isempty(L.plan) && isfield(L.plan, 'dayCount')
        i = find(L.plan.days == day, 1);
        if ~isempty(i) && L.plan.dayCount(i) > 0, tf = true;  return, end
    end
    if ~isempty(L.commonPlan)
        j = find(L.commonPlan.days == day, 1);
        if ~isempty(j) && L.commonPlan.dayCount(j) > 0, tf = true; end
    end
end

function v = priorArcsOfDay(L, day, obsType)
%PRIORARCSOFDAY  [1] when prior group obsType has rows on this day, else [].
    v = zeros(1, 0);
    if isempty(L.priors) || isempty(L.plan), return, end
    i = find(L.plan.days == day, 1);
    g = find(strcmp({L.priors.obsType}, obsType), 1);
    if isempty(i) || isempty(g) || ~any(L.priors(g).dayOf == i), return, end
    v = 1;
end

function P = validatePriors(P, plan, sysID, localParams, nameCol)
%VALIDATEPRIORS  The local prior groups: fields, widths, one day per row (R5-D8 round 5/6).
%   Each element is one GROUP, declared EITHER by .group -- a regular
%   expression matched against the name column of localParams; the loader
%   builds L with one row per matched template row per instance (day / arc /
%   once), from the layout it owns -- OR by .L (q x nl in the columns of B,
%   used as is). Never both. Then .x0 (q x 1, default 0), .sigma (> 0,
%   default 1: the ONE sigma the rows share), .from (>= 1, default 1), .name
%   (default the pattern or 'prior<g>'), .estimate "none" | "scale" (default
%   "none": sigma fixed; "scale": sigma estimated every pass by the Helmert /
%   Koch-Kusche rule, exactly as for data.priors -- the group is a stream of
%   its own). Every row must lie inside ONE day's local block (or the system
%   block); a row coupling two days is refused -- couple days through common
%   parameters.
    if isempty(P), P = struct('L', {}, 'x0', {}, 'sigma', {}, 'from', {}, 'name', {}, 'estimate', {}, 'obsType', {}, 'dayOf', {});  return, end
    if ~isstruct(P)
        error('lsqvce:makeLoader:badPrior', 'Priors must be a struct array (L, x0, sigma, from, name, estimate).');
    end
    if isempty(plan)
        error('lsqvce:makeLoader:badPrior', 'Priors need SCOPED localParams (the v4 tables).');
    end
    out = struct('L', {}, 'x0', {}, 'sigma', {}, 'from', {}, 'name', {}, 'estimate', {}, 'obsType', {}, 'dayOf', {});
    for g = 1:numel(P)
        e = P(g);
        bad = setdiff(fieldnames(e), {'L', 'x0', 'sigma', 'from', 'name', 'estimate', 'group'});
        if ~isempty(bad)
            error('lsqvce:makeLoader:badPrior', 'system %s, Priors(%d): unknown field ''%s''.', char(string(sysID)), g, bad{1});
        end
        hasL = isfield(e, 'L') && ~isempty(e.L);
        hasG = isfield(e, 'group') && ~isempty(e.group);
        if hasL && hasG
            error('lsqvce:makeLoader:badPrior', 'system %s, Priors(%d): give EITHER group (a name pattern) OR L, not both.', char(string(sysID)), g);
        end
        if hasG
            e.L = rowsOfPattern(e.group, localParams, nameCol, plan);
            if isempty(e.L)
                error('lsqvce:makeLoader:badPrior', ...
                      'system %s, Priors(%d): the pattern ''%s'' matches no name in column %d of localParams.', ...
                      char(string(sysID)), g, char(string(e.group)), nameCol);
            end
            if ~isfield(e, 'name') || isempty(e.name), e.name = char(string(e.group)); end
        end
        if ~isfield(e, 'L') || isempty(e.L) || ~isnumeric(e.L) || size(e.L, 2) ~= plan.nl
            error('lsqvce:makeLoader:badPrior', ...
                  'system %s, Priors(%d): L must be q x %d (the columns of this system''s B; see lsqvce.layout); got %s.', ...
                  char(string(sysID)), g, plan.nl, mat2str(size(e.L)));
        end
        Lg = full(double(e.L));  q = size(Lg, 1);
        s = 1;  if isfield(e, 'sigma') && ~isempty(e.sigma), s = double(e.sigma); end
        if ~isscalar(s) || ~isfinite(s) || s <= 0
            error('lsqvce:makeLoader:badPrior', 'system %s, Priors(%d): sigma must be a positive scalar.', char(string(sysID)), g);
        end
        x0 = zeros(q, 1);
        if isfield(e, 'x0') && ~isempty(e.x0)
            x0 = double(e.x0(:));
            if isscalar(x0), x0 = repmat(x0, q, 1); end   % R6 E-14: a scalar x0 is every row's
            if numel(x0) ~= q                              % (the data.priors rule, fit>priorsOf)
                error('lsqvce:makeLoader:badPrior', 'system %s, Priors(%d): x0 must have %d entries (or be a scalar).', char(string(sysID)), g, q);
            end
        end
        fr = 1;  if isfield(e, 'from') && ~isempty(e.from), fr = double(e.from); end
        nm = sprintf('prior%d', g);
        if isfield(e, 'name') && ~isempty(e.name), nm = char(string(e.name)); end
        est = 'none';
        if isfield(e, 'estimate') && ~isempty(e.estimate)
            est = lower(char(string(e.estimate)));
            if ~any(strcmp(est, {'none', 'scale'}))
                error('lsqvce:makeLoader:badPrior', 'system %s, Priors(%d): estimate must be "none" or "scale".', char(string(sysID)), g);
            end
        end
        % the day of every row: the block its nonzero columns lie in
        dayOf = zeros(q, 1);
        for r = 1:q
            cols = find(Lg(r, :));
            if isempty(cols)
                error('lsqvce:makeLoader:badPrior', 'system %s, Priors(%d): row %d of L is all zero.', char(string(sysID)), g, r);
            end
            d = arrayfun(@(c) dayOfColumn(plan, c), cols);
            if numel(unique(d)) ~= 1
                error('lsqvce:makeLoader:badPrior', ...
                      ['system %s, Priors(%d): row %d couples the local blocks of different days; ' ...
                       'a local prior row lives in one day -- couple days through common parameters.'], char(string(sysID)), g, r);
            end
            dayOf(r) = d(1);
        end
        out(end+1) = struct('L', Lg, 'x0', x0, 'sigma', s, 'from', fr, 'name', nm, 'estimate', est, ...
                            'obsType', ['prior_' matlab.lang.makeValidName(nm)], 'dayOf', dayOf); %#ok<AGROW>
    end
    P = out;
end

function Lg = rowsOfPattern(pattern, localParams, nameCol, plan)
%ROWSOFPATTERN  L for a 'group' pattern: one unit row per column of B whose
%   template row's name (column nameCol of localParams) matches the regexp.
    if nameCol > size(localParams, 2)
        error('lsqvce:makeLoader:badPrior', 'the name column %d does not exist in localParams (%d columns).', nameCol, size(localParams, 2));
    end
    names = localParams(:, nameCol);
    hit = false(size(names));
    for r = 1:numel(names)
        v = names{r};
        if ischar(v) || (isstring(v) && isscalar(v))
            hit(r) = ~isempty(regexp(char(v), char(string(pattern)), 'once'));
        end
    end
    loc = lsqvce.engine.localColumns(plan);
    cols = find(arrayfun(@(c) ~isempty(c.row) && hit(c.row), loc));
    Lg = zeros(numel(cols), plan.nl);
    for k = 1:numel(cols), Lg(k, cols(k)) = 1; end
end

function i = dayOfColumn(plan, c)
%DAYOFCOLUMN  Calendar index of the day whose local block holds column c
%   (the system block counts as day 1).
    if c > plan.sysOffset, i = 1;  return, end
    for i = 1:numel(plan.days)
        lo = Inf;  hi = -Inf;
        if ~isempty(plan.offset{i}), lo = min(lo, plan.offset{i}(1) + 1);  hi = max(hi, plan.offset{i}(end) + plan.count{i}(end)); end
        if plan.dayCount(i) > 0, lo = min(lo, plan.dayOffset(i) + 1);  hi = max(hi, plan.dayOffset(i) + plan.dayCount(i)); end
        if c >= lo && c <= hi, return, end
    end
    error('lsqvce:makeLoader:badPrior', 'column %d of L belongs to no day''s local block.', c);
end

function out = priorArc(L, day, obsType, arcID, opt, cov)
%PRIORARC  The pseudo-arc of ONE local prior group on one day (R5-D8 round 5):
%   its rows whose columns lie in this day's block -- A = 0, B = L/sigma,
%   l = x0/sigma, unit covariance. A group with from > opt.pass is not in.
    out = [];
    if arcID ~= 1 || isempty(L.priors) || isempty(L.plan), return, end
    i = find(L.plan.days == day, 1);
    g = find(strcmp({L.priors.obsType}, obsType), 1);
    if isempty(i) || isempty(g), return, end
    P = L.priors(g);
    pass = 1;
    if isstruct(opt) && isfield(opt, 'pass') && ~isempty(opt.pass), pass = opt.pass; end
    if P.from > pass, return, end
    rows = find(P.dayOf == i);
    if isempty(rows), return, end
    out.B = P.L(rows, :) / P.sigma;
    out.l = P.x0(rows) / P.sigma;
    m = numel(out.l);
    out.A = zeros(m, L.nc);
    out.t = (1:m).';
    % the group's law is its declared unit spec WITH the run's override (the
    % estimated level, Estimate="noise" with estimate "scale") -- the same
    % path every data arc takes, so the level whitens the next pass
    cm = cov.(obsType);
    if isstruct(opt) && isfield(opt, 'covOverride')
        cm = applyCovOverride(cm, opt.covOverride, obsType, opt);
    end
    [out.Cfcn, out.Wfcn, info2] = lsqvce.engine.covModel(cm, out.t, struct());
    out.filterMemory = info2.filterMemory;
end
function v = arcVec(arcs, obsType, day, dayFcn)
%ARCVEC  Resolve the arc-ID vector for one obsType and ONE DAY.
%   Static forms (count / vector / per-obsType struct) ignore the day; the
%   function_handle form @(day[,obsType]) is a ragged schedule (WP-52). No
%   data is read for a declared schedule. With NO schedule (R6 WP-102) the
%   day is read and its cells are its arcs, ids 1:nArc.
%   OUTPUT v 1xA double  the arc IDs present on THIS day.
    if isempty(arcs)
        if isempty(dayFcn)
            error('lsqvce:makeLoader:localPlan', ...
                  ['the arc set of day %s (stream %s) is needed before any day is read, ' ...
                   'and no ''Arcs'' schedule is declared.'], string(day), obsType);
        end
        D = dayFcn(day, obsType);
        v = 1:numel(D.A);
    elseif isstruct(arcs)
        if ~isfield(arcs, obsType)
            error('lsqvce:makeLoader:noArcs', 'arcs has no field for ''%s''.', obsType);
        end
        v = resolveArcSpec(arcs.(obsType), obsType, day);
    else
        v = resolveArcSpec(arcs, obsType, day);
    end
end

function v = resolveArcSpec(spec, obsType, day)
%RESOLVEARCSPEC  One declared arc spec -> the arc-ID row vector for one day.
    if isa(spec, 'function_handle')
        if nargin(spec) == 1
            v = spec(day);
        else
            v = spec(day, obsType);
        end
        v = checkArcIds(v, obsType, day, 'the ''arcs'' function');
    elseif isscalar(spec) && isnumeric(spec)
        v = 1:spec;
    else
        v = spec(:).';
    end
end

function ids = dayArcIds(arcs, obsType, day, D)
%DAYARCIDS  The arc ids of a BUILT day, in cell order (R6 WP-102).
%   No schedule: 1:nArc. A declared schedule: its ids, which must be as many
%   as the day's cells (cell a is arc ids(a); an empty cell is that id's gap).
    nArc = numel(D.A);
    if isempty(arcs)
        ids = 1:nArc;
        return
    end
    ids = arcVec(arcs, obsType, day, []);
    if numel(ids) ~= nArc
        error('lsqvce:makeLoader:arcCount', ...
              ['buildDay(%s,''%s'') returned %d arc cell(s) but the declared schedule names ' ...
               '%d arc(s) (%s): the cells are the scheduled arcs in order -- return an EMPTY ' ...
               'cell for an arc the day does not carry, or fix the schedule.'], ...
              string(day), obsType, nArc, numel(ids), mat2str(ids));
    end
end

function plan = buildLocalPlan(spec, arcs, obsTypes, days)
%BUILDLOCALPLAN  Assign every (day, arc) its local-block ordinal and column
%   offset, and DERIVE nl (WP-52).
%   The arc set of a day is the UNION over the loader's obsTypes: an arc seen
%   by two observation types is ONE arc state and therefore ONE local block.
%   Requires a DECLARED schedule -- static arcs or an @(day[,obsType])
%   resolver -- because nl must be known at construction, without a read.
%   INPUT
%     spec      1x1 double | function_handle  locals per arc, or @(day,arc)
%     arcs      the arcs spec; obsTypes 1xT cell; days 1xD double
%   OUTPUT
%     plan  1x1 struct  .days (1xD) .arcId {1xD cell} .offset {1xD cell}
%                       .count {1xD cell} .ordinal {1xD cell}
%                       .nl (1x1) .nArcTot (1x1)
    if isempty(days)
        error('lsqvce:makeLoader:localPlan', ...
            ['''localsPerArc'' needs the day list: pass ''days'' so the ' ...
             'local blocks can be numbered over the whole adjustment.']);
    end
    if isempty(arcs)
        error('lsqvce:makeLoader:localPlan', ...
            ['''localsPerArc'' needs a declared ''Arcs'' schedule: nl must be known ' ...
             'WITHOUT reading data. Declare a count, an id vector or an @(day,obsType) resolver.']);
    end
    days = days(:).';
    D = numel(days);
    plan = struct('days', days, 'arcId', {cell(1,D)}, 'offset', {cell(1,D)}, ...
                  'count', {cell(1,D)}, 'ordinal', {cell(1,D)}, ...
                  'nl', 0, 'nArcTot', 0);
    off = 0;  ord = 0;
    for i = 1:D
        d = days(i);
        ids = zeros(1,0);
        for k = 1:numel(obsTypes)
            ids = union(ids, arcVec(arcs, obsTypes{k}, d, []));
        end
        ids = sort(ids(:).');
        cnt = zeros(1, numel(ids));
        offs = zeros(1, numel(ids));
        ords = zeros(1, numel(ids));
        for a = 1:numel(ids)
            if isa(spec, 'function_handle')
                cnt(a) = spec(d, ids(a));
            else
                cnt(a) = spec;
            end
            if ~isscalar(cnt(a)) || cnt(a) < 0 || cnt(a) ~= round(cnt(a))
                error('lsqvce:makeLoader:localPlan', ...
                    ['''localsPerArc'' must give a non-negative integer ' ...
                     'count (day %s, arc %d gave %g).'], ...
                     string(d), ids(a), cnt(a));
            end
            offs(a) = off;   off = off + cnt(a);
            ord = ord + 1;   ords(a) = ord;
        end
        plan.arcId{i} = ids;  plan.offset{i} = offs;
        plan.count{i} = cnt;  plan.ordinal{i} = ords;
    end
    plan.nl = off;  plan.nArcTot = ord;
end

function d = metaDay(args, days)
%METADAY  Which day L.meta() answers for: the one given, else the first
%   declared day, else 1.
    if ~isempty(args) && ~isempty(args{1})
        d = args{1};
    elseif ~isempty(days)
        d = days(1);
    else
        d = 1;
    end
end

function info = arcInfoOf(plan, day)
%ARCINFOOF  The local-block slice buildDay gets for ONE day (WP-52; scoped
%   fields since L11 WP-75).
%   OUTPUT info 1x1 struct | []  .arcId .offset .count .ordinal (1xA each),
%          .nl (total over the adjustment), .offsetOf @(arcId)->column offset;
%          for a SCOPED plan also .dayOffset/.dayCount (this day's 'day'
%          block; [] / 0 when the day has no arcs) and .sysOffset/.sysCount
%          (the trailing 'system' block).
    info = [];
    if isempty(plan), return, end
    i = find(plan.days == day, 1);
    if isempty(i) && ~(isfield(plan, 'sysCount') && plan.sysCount > 0), return, end
    if isempty(i)
        arcIds = zeros(1,0); offs = zeros(1,0); cnts = zeros(1,0); ords = zeros(1,0);
    else
        arcIds = plan.arcId{i}; offs = plan.offset{i};
        cnts = plan.count{i};   ords = plan.ordinal{i};
    end
    info = struct('arcId', arcIds, 'offset', offs, ...
                  'count', cnts, 'ordinal', ords, ...
                  'nl', plan.nl, 'offsetOf', [], ...
                  'dayIndex', i, 'nDays', numel(plan.days));   % R5-D8 round 5: which day of the calendar
    info.offsetOf = @(a) offs(find(arcIds == a, 1));
    if isfield(plan, 'sysCount')                    % scoped plan (v4)
        dOff = [];  dCnt = 0;
        if ~isempty(i) && plan.dayCount(i) > 0
            dOff = plan.dayOffset(i);  dCnt = plan.dayCount(i);
        end
        info.dayOffset = dOff;   info.dayCount = dCnt;
        info.sysOffset = plan.sysOffset;  info.sysCount = plan.sysCount;
        % R5-D8 round 3: WHICH template rows are present on this day (their
        % order inside every arc block / the day block); without validity
        % columns these are plan.rows.arc / plan.rows.day on every day
        if ~isempty(i) && isfield(plan, 'arcRowsOf')
            info.arcRows = plan.arcRowsOf{i};  info.dayRows = plan.dayRowsOf{i};
        else
            info.arcRows = plan.rows.arc;      info.dayRows = plan.rows.day;
        end
        info.sysRows = plan.rows.system;
    end
end

function info = addCommonInfo(info, commonPlan, day)
%ADDCOMMONINFO  Merge the COMMON-plan slice for ONE day into buildDay's 6th
%   argument (L11 WP-75, v4): .nc, .commonCount (the leading 'common' block),
%   .commonDayOffset/.commonDayCount (THIS day's block; [] / 0 when the
%   calendar carries no block for it).
    if isempty(commonPlan), return, end
    if isempty(info), info = struct(); end
    info.nc          = commonPlan.nc;
    info.commonCount = commonPlan.commonCount;
    j = find(commonPlan.days == day, 1);
    if isempty(j) || commonPlan.dayCount(j) == 0
        info.commonDayOffset = [];  info.commonDayCount = 0;
    else
        info.commonDayOffset = commonPlan.dayOffset(j);
        info.commonDayCount  = commonPlan.dayCount(j);
    end
    info.commonDayRows = zeros(1, 0);                               % R5-D8 round 3
    if isfield(commonPlan, 'rows'), info.commonDayRows = commonPlan.rows.day; end   % (engine-tier plans may omit rows)
    if ~isempty(j) && isfield(commonPlan, 'dayRowsOf'), info.commonDayRows = commonPlan.dayRowsOf{j}; end
end

function plan = buildScopePlan(tbl, arcs, obsTypes, days, valCols)
    if nargin < 5, valCols = []; end
%BUILDSCOPEPLAN  Expand a SCOPED local table over the calendar (L11 WP-75).
%   Column 1 of TBL is the scope token: 'arc' | 'day' | 'system'; rows are
%   TEMPLATES. Layout, binding: per day (declared order) [arc blocks in
%   sorted-id order, each = the 'arc' rows | the 'day' block = the 'day'
%   rows], then ONE trailing 'system' block. A day whose DECLARED schedule
%   carries no arcs gets NO blocks at all; without a schedule (R6 WP-102)
%   every calendar day owns its 'day' block -- the calendar is the
%   declaration, and a day that then returns no data is refused at its
%   'meta' probe (emptyDay). 'arc' rows need the schedule (nl must be known
%   before any day is read). nl is DERIVED. Generalizes the WP-52 plan
%   (that plan = a table of n 'arc' rows).
%   INPUT  tbl RxC cell (C >= 1); arcs / obsTypes / days as buildLocalPlan
%   OUTPUT plan 1x1 struct  the WP-52 fields (.days .arcId .offset .count
%          .ordinal .nl .nArcTot) + .dayOffset .dayCount (1xD) + .sysOffset
%          .sysCount + .rows (struct of row indices per scope, so callers can
%          map block columns back to template rows)
    tok = cell(size(tbl, 1), 1);
    for r = 1:size(tbl, 1)
        v = tbl{r, 1};
        if ~(ischar(v) || (isstring(v) && isscalar(v))) || ...
                ~any(strcmp(char(v), {'arc', 'day', 'system'}))
            error('lsqvce:makeLoader:badScope', ...
                ['localParams row %d: column 1 must be the scope token ' ...
                 '''arc'' | ''day'' | ''system'' (got %s).'], r, class(v));
        end
        tok{r} = char(v);
    end
    iArc = find(strcmp(tok, 'arc')).';  iDay = find(strcmp(tok, 'day')).';
    iSys = find(strcmp(tok, 'system')).';
    nArcRows = numel(iArc);  nDayRows = numel(iDay);  nSysRows = numel(iSys);
    if (nArcRows + nDayRows) > 0 && isempty(days)
        error('lsqvce:makeLoader:localPlan', ...
            ['scoped ''arc''/''day'' rows need the day list: pass ''days'' so ' ...
             'the blocks can be numbered over the whole adjustment.']);
    end
    declared = ~isempty(arcs);
    if nArcRows > 0 && ~declared
        error('lsqvce:makeLoader:localPlan', ...
            ['scoped ''arc'' rows need a declared ''Arcs'' schedule: nl counts one block per ' ...
             '(day, arc) and must be known WITHOUT reading data. Declare a count, an id ' ...
             'vector or an @(day,obsType) resolver.']);
    end
    days = days(:).';
    D = numel(days);
    plan = struct('days', days, 'arcId', {cell(1,D)}, 'offset', {cell(1,D)}, ...
                  'count', {cell(1,D)}, 'ordinal', {cell(1,D)}, ...
                  'dayOffset', zeros(1,D), 'dayCount', zeros(1,D), ...
                  'sysOffset', 0, 'sysCount', nSysRows, ...
                  'nl', 0, 'nArcTot', 0, ...
                  'rows', struct('arc', iArc, 'day', iDay, 'system', iSys), ...
                  'arcRowsOf', {cell(1,D)}, 'dayRowsOf', {cell(1,D)});
    % R5-D8 round 3: validity-aware expansion -- a 'day'/'arc' template row is
    % instantiated only on the days its validity overlaps (valCols name the
    % columns; [] = every row is valid always, the v4 behaviour)
    valid = validityOf(tbl, valCols);              % @(row, day) -> logical
    off = 0;  ord = 0;
    for i = 1:D
        d = days(i);
        ids = zeros(1,0);
        if declared
            for k = 1:numel(obsTypes)
                ids = union(ids, arcVec(arcs, obsTypes{k}, d, []));
            end
        end
        ids = sort(ids(:).');
        A = numel(ids);
        arcRows = iArc(arrayfun(@(r) valid(r, d), iArc));
        dayRows = iDay(arrayfun(@(r) valid(r, d), iDay));
        nA = numel(arcRows);  nD = numel(dayRows);
        offs = zeros(1, A);  cnts = repmat(nA, 1, A);  ords = zeros(1, A);
        for a = 1:A
            offs(a) = off;  off = off + nA;
            ord = ord + 1;  ords(a) = ord;
        end
        plan.arcId{i} = ids;  plan.offset{i} = offs;
        plan.count{i} = cnts; plan.ordinal{i} = ords;
        plan.arcRowsOf{i} = arcRows;  plan.dayRowsOf{i} = dayRows;
        if (A > 0 || ~declared) && nD > 0           % a day WITH arcs gets its block
            plan.dayOffset(i) = off;  plan.dayCount(i) = nD;
            off = off + nD;
        end
    end
    plan.sysOffset = off;
    off = off + nSysRows;
    plan.nl = off;  plan.nArcTot = ord;
end

function v = checkArcIds(v, obsType, day, who)
%CHECKARCIDS  A resolved arc list must be a row of unique positive integers.
    if isempty(v), v = zeros(1,0); return, end
    v = v(:).';
    ok = isnumeric(v) && all(isfinite(v)) && all(v >= 1) && all(v == round(v));
    if ~ok || numel(unique(v)) ~= numel(v)
        error('lsqvce:makeLoader:badArcs', ...
            ['%s returned an invalid arc list for day %s, obsType %s: ' ...
             'expected unique positive integers.'], who, string(day), obsType);
    end
end

function D = dayOf(buildDay, cache, capacity, day, obsType, common, local, procOpt, plan, commonPlan)
%DAYOF  Return the built+validated day struct for (day,obsType), using a
%   cache of CAPACITY resident days (0 = none, 1 = the single slot of v4,
%   Inf = every day of the run -- R5 M4 4d: a two-stage pass builds every day
%   twice, accumulate + sweep, and a VCE run 2n+1 times; with the calendar
%   resident it is ONCE per run). OUTPUT D 1x1 struct  the day data.
    if nargin < 10, commonPlan = []; end
    key = sprintf('%d|%s', day, obsType);
    useCache = capacity > 0;
    if useCache && isKey(cache, key)
        D = cache(key);  return
    end
    arcInfo = arcInfoOf(plan, day);           % WP-52: [] unless a plan exists
    arcInfo = addCommonInfo(arcInfo, commonPlan, day);   % L11: v4 common slice
    D = callBuildDay(buildDay, day, obsType, common, local, procOpt, arcInfo);
    D = validateDay(D, day, obsType);
    if useCache
        while cache.Count >= capacity          % evict the oldest (insertion order)
            ks = cache.keys;  remove(cache, ks{1});
        end
        cache(key) = D;              %#ok<NASGU>
    end
end

function D = callBuildDay(buildDay, day, obsType, common, local, procOpt, arcInfo)
%CALLBUILDDAY  Call buildDay with as many of the 6 args as it declares, so an
%   OLD 2-arg builder @(day,obsType) keeps working while a NEW builder can take
%   @(day,obsType,common,local,procOpt,arcInfo). Anonymous handles report a
%   fixed nargin; a named function that uses varargin reports a negative
%   nargin, in which case we pass all six. WP-52 added the 6th (arcInfo, the
%   per-day local-block plan; empty unless 'localsPerArc' was declared).
    n = nargin(buildDay);
    args = {day, obsType, common, local, procOpt, arcInfo};
    if n < 0 || n >= 6
        D = buildDay(args{:});                 % varargin, or takes all six
    else
        D = buildDay(args{1:max(2,n)});        % pass exactly what it declares (>=2)
    end
end

function D = validateDay(D, day, obsType)
%VALIDATEDAY  Enforce the v6 day contract (R6 WP-102): A, B, l, t are 1 x nArc
%   cells -- THE CELLS ARE THE ARCS -- agreeing in number, each arc's four
%   entries sharing a row count (0 rows = a data gap); Cep, when given, a
%   matching cell of per-arc error information; nAxis derived from Cep or
%   declared, and every arc holding whole epochs. No other field is allowed
%   (arcIdx / arcId are the v5 way of saying what the cells say). Checked
%   HERE, at the day, with day and stream in the message, before any shape
%   could surface per arc in a kernel id. OUTPUT D the checked day.
    allowed = {'A', 'B', 'l', 't', 'Cep', 'nAxis'};
    bad = setdiff(fieldnames(D), allowed);
    if ~isempty(bad)
        hint = '';
        if any(strcmp(bad{1}, {'arcIdx', 'arcId'}))
            hint = ' -- retired in v6: return one cell per arc (an empty cell for a gap)';
        end
        error('lsqvce:makeLoader:badDay', ...
              'buildDay(%s,''%s'') returned the unknown field .%s%s. Allowed: %s.', ...
              string(day), obsType, bad{1}, hint, strjoin(allowed, ', '));
    end
    for f = {'A', 'B', 'l', 't'}
        if ~isfield(D, f{1})
            error('lsqvce:makeLoader:badDay', ...
                'buildDay(%s,''%s'') must return field .%s.', string(day), obsType, f{1});
        end
        if ~iscell(D.(f{1})) || (~isvector(D.(f{1})) && ~isempty(D.(f{1})))
            error('lsqvce:makeLoader:badDay', ...
                ['buildDay(%s,''%s''): .%s must be a 1 x nArc cell, one cell per arc of the ' ...
                 'day ({%s} for a single arc; v6, the cells are the arcs); got a %s.'], ...
                string(day), obsType, f{1}, f{1}, class(D.(f{1})));
        end
        D.(f{1}) = reshape(D.(f{1}), 1, []);
    end
    nArc = numel(D.A);
    if numel(D.B) ~= nArc || numel(D.l) ~= nArc || numel(D.t) ~= nArc
        error('lsqvce:makeLoader:dayShape', ...
              'buildDay(%s,''%s''): A, B, l, t must hold the same number of arc cells (A has %d, B %d, l %d, t %d).', ...
              string(day), obsType, nArc, numel(D.B), numel(D.l), numel(D.t));
    end
    mArc = zeros(1, nArc);
    for a = 1:nArc
        m = size(D.A{a}, 1);
        if numel(D.l{a}) ~= m || numel(D.t{a}) ~= m || size(D.B{a}, 1) ~= m
            error('lsqvce:makeLoader:dayShape', ...
                'buildDay(%s,''%s''): arc %d: A, B, l, t must share the same number of rows.', ...
                string(day), obsType, a);
        end
        D.l{a} = D.l{a}(:);  D.t{a} = D.t{a}(:);
        mArc(a) = m;
    end
    % R6 WP-90: the ERROR INFORMATION of the day -- per arc since v6. An
    % entry is nEpoch x 1 (a scalar stream's per-epoch sigma), nEpoch x d
    % (per-axis sigmas, a diagonal D) or nEpoch x d x d (the full per-epoch
    % covariance). The axis count is DERIVED from it and is ONE for the day;
    % an explicit D.nAxis must agree (one fact, one authority).
    hasCep = isfield(D, 'Cep') && ~isempty(D.Cep);
    hasAx  = isfield(D, 'nAxis') && ~isempty(D.nAxis);
    if hasAx && ~(isnumeric(D.nAxis) && isscalar(D.nAxis) && D.nAxis >= 1 && D.nAxis == round(D.nAxis))
        error('lsqvce:makeLoader:axisConflict', ...
              'buildDay(%s,''%s''): D.nAxis must be a positive integer scalar.', string(day), obsType);
    end
    if hasCep
        if ~iscell(D.Cep) || numel(D.Cep) ~= nArc
            error('lsqvce:makeLoader:badDayCep', ...
                  ['buildDay(%s,''%s''): D.Cep must be a 1 x nArc cell like A, one entry per arc ' ...
                   '(the day has %d arcs); got a %s of %d.'], ...
                   string(day), obsType, nArc, class(D.Cep), numel(D.Cep));
        end
        D.Cep = reshape(D.Cep, 1, []);
        dCep = [];
        for a = 1:nArc
            C = D.Cep{a};
            if mArc(a) == 0
                if ~isempty(C)
                    error('lsqvce:makeLoader:badDayCep', ...
                          'buildDay(%s,''%s''): arc %d is a gap (0 rows) but Cep{%d} is not empty.', ...
                          string(day), obsType, a, a);
                end
                continue
            end
            sz = size(C);
            ok = isnumeric(C) && ~isempty(C) && (numel(sz) == 2 || (numel(sz) == 3 && sz(2) == sz(3)));
            if ok && numel(sz) == 2 && sz(1) == 1 && sz(2) > 1 && mArc(a) > 1
                ok = false;                              % a ROW vector: nEpoch x 1 is a column
            end
            if ~ok
                error('lsqvce:makeLoader:badDayCep', ...
                      ['buildDay(%s,''%s''): Cep{%d} must be nEpoch x 1 (a per-epoch sigma), nEpoch x d ' ...
                       '(per-axis sigmas) or nEpoch x d x d (per-epoch covariances); got %s.'], ...
                       string(day), obsType, a, mat2str(sz));
            end
            if sz(1) * sz(2) ~= mArc(a)
                error('lsqvce:makeLoader:badDayCep', ...
                      ['buildDay(%s,''%s''): Cep{%d} has %d epochs of %d axes = %d rows, but arc %d has ' ...
                       '%d rows (rows are epoch-major: nEpoch * nAxis).'], ...
                       string(day), obsType, a, sz(1), sz(2), sz(1) * sz(2), a, mArc(a));
            end
            if isempty(dCep)
                dCep = sz(2);
            elseif sz(2) ~= dCep
                error('lsqvce:makeLoader:axisConflict', ...
                      'buildDay(%s,''%s''): Cep{%d} has %d axes, Cep{1..} %d; one axis count per day.', ...
                      string(day), obsType, a, sz(2), dCep);
            end
        end
        if ~isempty(dCep)
            if hasAx && D.nAxis ~= dCep
                error('lsqvce:makeLoader:axisConflict', ...
                      ['buildDay(%s,''%s''): D.nAxis = %d but D.Cep has %d axes; the axis count is ' ...
                       'derived from Cep -- state it once.'], string(day), obsType, D.nAxis, dCep);
            end
            D.nAxis = dCep;
        end
    end
    if isfield(D, 'nAxis') && ~isempty(D.nAxis) && D.nAxis > 1
        if mod(sum(mArc), D.nAxis) ~= 0
            error('lsqvce:makeLoader:axisConflict', ...
                  'buildDay(%s,''%s''): %d rows are not a multiple of nAxis = %d (rows must be epoch-major).', ...
                  string(day), obsType, sum(mArc), D.nAxis);
        end
        for a = 1:nArc
            if mod(mArc(a), D.nAxis) ~= 0
                error('lsqvce:makeLoader:badArc', ...
                      'buildDay(%s,''%s''): arc %d has %d rows, not whole epochs of %d axes (rows must be epoch-major).', ...
                      string(day), obsType, a, mArc(a), D.nAxis);
            end
        end
    end
end

function out = batchOf(dayFcn, arcs, cov, day, obsType, arcID, opt)
%BATCHOF  One arc (arcID) of a day as a batch struct (A,B,l,...), with the
%   covariance/whitening hook cov applied. OUTPUT out 1x1 struct the batch,
%   or [] for a data gap (an empty cell -- every engine skips it).
    D = dayFcn(day, obsType);
    % R6 WP-102: the cells are the arcs. Without a schedule arc a is cell a;
    % with one, cell a is the a-th scheduled id (dayArcIds checks the count).
    % An id the schedule does not name is a mis-declared STATIC schedule and
    % stays an error -- a different mistake from a ragged one.
    ids = dayArcIds(arcs, obsType, day, D);
    pos = find(ids == arcID, 1);
    if isempty(pos)
        error('lsqvce:makeLoader:badArc', ...
            'Arc %d does not exist (day %s, obsType %s has arcs %s).', ...
            arcID, string(day), obsType, mat2str(ids));
    end
    if isempty(D.l{pos})
        out = [];                                    % the data gap of this arc
        return
    end
    out.A = D.A{pos};
    out.B = D.B{pos};
    out.l = D.l{pos};
    out.t = D.t{pos};

    cm = cov.(obsType);
    if isfield(opt,'covOverride')
        cm = applyCovOverride(cm, opt.covOverride, obsType, opt);
    end

    % R6 WP-90: the ERROR INFORMATION travels with the day -- per arc since
    % v6 (WP-102): D.Cep{pos} is THIS arc's, validated at the day, and is
    % handed to covModel as cm.Cep, where ONE composition joins it with the
    % declared temporal law. A Cep on the spec AND on the day is two
    % authorities for one fact: refused.
    d = lsqvce.util.getfielddef(D, 'nAxis', 1);
    specCep = isfield(cm, 'Cep') && ~isempty(cm.Cep);
    isVar = strcmpi(lsqvce.util.getfielddef(cm, 'type', ''), 'var');
    if isfield(D, 'Cep') && ~isempty(D.Cep)
        if specCep
            error('lsqvce:makeLoader:cepTwice', ...
                  ['covariance.%s carries a Cep and buildDay(%s,''%s'') returns one: state the ' ...
                   'per-epoch error information ONCE, on the day.'], obsType, string(day), obsType);
        end
        cm.Cep = D.Cep{pos};
    elseif d > 1 && ~specCep && ~isVar
        % a multi-axis stream without per-epoch scaling is composed with
        % D = I -- said once per run, so the omission is distinguishable from
        % a choice (R5: never overrule, never assume silently)
        lsqvce.util.log.warnOnce('lsqvce:makeLoader:noCep', ...
            ['buildDay(%s,''%s'') returns a %d-axis day without Cep: the axes are weighted by the ' ...
             'temporal law alone (D = I). Return D.Cep (per arc nEpoch x %d sigmas or nEpoch x %d x %d) ' ...
             'when the epochs are not equally accurate.'], string(day), obsType, d, d, d, d);
        cm.Cep = ones(numel(out.l) / d, d);
    end
    dDecl = lsqvce.util.getfielddef(cm, 'd', []);
    if ~isempty(dDecl) && ~isVar && dDecl ~= d
        error('lsqvce:makeLoader:axisConflict', ...
              ['covariance.%s declares d = %d but the day has %d axes (nAxis / Cep); the axis ' ...
               'count is stated once, on the day.'], obsType, dDecl, d);
    end

    [out.Cfcn, out.Wfcn, info] = lsqvce.engine.covModel(cm, out.t, opt);
    out.filterMemory = info.filterMemory;
    if isfield(info, 'Cep') && ~isempty(info.Cep)
        out.Cep = info.Cep;                          % the per-epoch descale travels with the arc
    end

    % R5 M2 (WP-81 round 2b): per-observation ROBUST weights ride INSIDE the
    % whitener, so every consumer of Wfcn (accumulation, residuals, the
    % diagnostics sweep) sees the same reweighted, law-whitened rows. The
    % driver lsqvce.engine.adjustRobust supplies opt.rowWeights, a
    % containers.Map keyed 'sysID|day|obsType|arcId' -> m_arc x 1 factors in
    % [0, 1]; an absent key means unit weights (pass 1). Zero factors are
    % floored at 1e-12 so a row keeps a (negligible) presence and its
    % residual can still be recovered as w / sqrt(factor).
    if isfield(opt, 'rowWeights') && ~isempty(opt.rowWeights)
        key = lsqvce.util.rowWeightKey(opt.sysID, day, obsType, arcID);
        if isKey(opt.rowWeights, key)
            fac = max(opt.rowWeights(key), 1e-12);
            if numel(fac) ~= numel(out.l)
                error('lsqvce:makeLoader:rowWeightSize', ...
                      ['rowWeights(%s) has %d entries; the arc has %d rows -- the ' ...
                       'builder must return the same arc on every call.'], ...
                      key, numel(fac), numel(out.l));
            end
            W0 = out.Wfcn;  sq = sqrt(fac(:));
            out.Wfcn = @(idx, A, B, l, o) weightedWhiten(W0, sq, idx, A, B, l, o);
        end
    end

    % Multi-axis (epoch-major) streams must tell the accumulator the axis count
    % so it batches on whole epochs and the whitener treats axes correctly:
    % the day's (derived) axis count, else the cov model's own d (a 'var' law).
    if d > 1
        out.nAxis = d;
    elseif isfield(info, 'nAxis') && info.nAxis > 1
        out.nAxis = info.nAxis;
    end
end

%==========================================================================
function [A, B, l] = weightedWhiten(W0, sq, idx, A, B, l, o)
%WEIGHTEDWHITEN  Law-whiten a batch, then scale its rows by sqrt(factor).
%   idx are the arc rows of this batch (overlap included), so sq(idx) is
%   aligned with W0's output rows.
    [A, B, l] = W0(idx, A, B, l, o);
    sc = sq(idx(:));
    A = sc .* A;  B = sc .* B;  l = sc .* l;
end

%==========================================================================
function f = validityOf(tbl, valCols)
%VALIDITYOF  @(row, day) -> is template row valid on calendar day `day`?
%   valCols = [cStart cEnd] (1-based columns of the table); a cell holds
%   [MJD sec] (-> MJD + sec/86400), a scalar fractional MJD, or [] (always).
%   Overlap rule: start < day + 1 && end > day (half-open day). R5-D8 round 3.
    if isempty(valCols)
        f = @(row, day) true;
        return
    end
    if numel(valCols) ~= 2 || any(valCols < 1) || any(valCols > size(tbl, 2)) || any(valCols ~= round(valCols))
        error('lsqvce:makeLoader:badValidity', ...
              'ValidityColumns must be [cStart cEnd] within 1..%d (got %s).', size(tbl, 2), mat2str(valCols));
    end
    n = size(tbl, 1);
    t0 = -inf(n, 1);  t1 = inf(n, 1);
    for r = 1:n
        a = toMjd(tbl{r, valCols(1)});  b = toMjd(tbl{r, valCols(2)});
        if ~isempty(a), t0(r) = a; end
        if ~isempty(b), t1(r) = b; end
        if t0(r) > t1(r)
            error('lsqvce:makeLoader:badValidity', 'row %d: validity start %g is after end %g.', r, t0(r), t1(r));
        end
    end
    f = @(row, day) t0(row) < day + 1 && t1(row) > day;
end

function m = toMjd(v)
    m = [];
    if isempty(v) || ~isnumeric(v), return, end
    if isscalar(v), m = double(v);  return, end
    if numel(v) == 2, m = double(v(1)) + double(v(2)) / 86400;  return, end
    error('lsqvce:makeLoader:badValidity', 'a validity entry must be [], a fractional MJD, or [MJD seconds].');
end
