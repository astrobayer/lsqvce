function L = makeLoader(buildDay, opt)
%LSQVCE.MAKELOADER  Turn a day-builder into an accumulator-ready Loader.
%   L = lsqvce.makeLoader(buildDay, Name, Value, ...)
%
%   You write the ONE part that is genuinely yours -- reading your data for a
%   day and returning its design matrix and observations -- and makeLoader
%   supplies the protocol the accumulator expects (the 'meta' discovery call,
%   the per-arc slicing, a single-day cache, and the covariance/whitener via
%   LSQVCE.COVMODEL). The result L plugs straight into LSQVCE.ADJUSTBLOCK or
%   LSQVCE.BLOCKACCUMULATOR/addSystem.
%
%   THE PART YOU WRITE
%     buildDay(day, obsType) -> struct with              (the classic 2-arg form)
%     buildDay(day, obsType, common, local, procOpt) -> struct   (parameter-driven)
%         .A       m x nc    common design matrix for the whole day
%         .B       m x nl    local design matrix (may be sparse)
%         .l       m x 1     observations
%         .t       m x 1     epoch timestamps [s]
%         .arcIdx  1 x nArc  cell of row-index vectors, one per arc (optional; if
%                            omitted, the whole day is one arc). Arc LENGTHS may
%                            differ freely -- they are just row counts, and the
%                            whitener is built per arc.
%         .arcId   1 x nArc  OPTIONAL arc ids naming the arcs THIS day carries
%                            (WP-52). With it, arcs are addressed by id rather
%                            than by position, so a day may omit arcs in the
%                            middle of the schedule; an id this day does not
%                            carry is treated as a data gap.
%     It must be DETERMINISTIC: the accumulator reads each day a second time for
%     residuals, and must get identical data back.
%
%     buildDay(day, obsType, common, local, procOpt, arcInfo) -> struct
%         also receives the per-day local-block plan when 'localsPerArc' is
%         declared (WP-52; empty otherwise).
%
%     buildDay is called with as many of the SIX arguments as it DECLARES, so an
%     existing @(day,obsType) builder keeps working unchanged; a builder that
%     also takes (common, local, procOpt) receives the parametrisation arrays and
%     processing options passed to makeLoader (below). A builder that uses
%     varargin receives all six. common/local have ONE ROW per unknown and 14
%     columns; nc and nl are DERIVED from their row counts (see 'common'/'local').
%
%   INPUT
%     buildDay  1x1  function_handle   @(day,obsType) day builder (above).
%                                      Required.
%     Name/Value options ([] or omitted -> defaults). Each user option has a
%     CANONICAL long-form name; the terse alias in (parentheses) is DEPRECATED
%     but still accepted for one version (a one-time warning nudges migration):
%       'observationTypes' 1xT cellstr/string  observation types, e.g.
%                                       {'pos','kbr'}. Required.  (alias: 'obsTypes')
%       'covariance' 1x1 struct         .(obsType) = a covariance spec for
%                                       LSQVCE.COVMODEL (e.g.
%                                       struct('type','arma','ar',..,'sigma',..)).
%                                       One field per obsType. Required. (alias: 'cov')
%       'arcs'     var  which arcs each day carries. STATIC forms (unchanged):
%                        1x1 integer    -> arcs 1:arcs for every obsType
%                        1xA vector     -> those arc ids for every obsType
%                        1x1 struct     -> .(obsType) = any form on this list
%                                       (default 1 -> a single arc per day)
%                       RAGGED forms (WP-52; the accumulator probes 'meta'
%                       PER DAY, so the schedule may differ from day to day):
%                        function_handle -> @(day) or @(day,obsType) returning
%                                       that day's arc ids. No data read.
%                        'fromDay'      -> read the day and take its .arcId,
%                                       or 1:numel(.arcIdx). Convenient, but it
%                                       reads data during the probe (the
%                                       single-day cache makes the following
%                                       arc reads free); prefer the resolver
%                                       when the schedule is known.
%                       With ragged arcs the day may declare .arcId (see THE
%                       PART YOU WRITE) so arcs are addressed BY ID and absent
%                       ids are skipped as data gaps -- no padding needed.
%       'localsPerArc' var  let makeLoader OWN the per-arc local bookkeeping
%                                       (WP-52). 1x1 integer, or @(day,arcId)
%                                       returning the count for that arc. It
%                                       numbers every (day, arc) over the whole
%                                       adjustment, DERIVES nl (so 'numLocal'
%                                       must be omitted or agree), and hands
%                                       each day its slice as buildDay's 6th
%                                       argument arcInfo:
%                                         .arcId .offset .count .ordinal 1xA
%                                         .nl        1x1  total local count
%                                         .offsetOf  @(arcId) -> column offset
%                                       Fill B(rows, off + (1:k)) and never
%                                       compute an offset. Needs 'days' and a
%                                       schedule resolvable without reading
%                                       data (not 'fromDay').  (default [])
%       'systemId' 1x1  system identifier carried on L    (default 1) (alias: 'sysID')
%                       MUST be NUMERIC or a SINGLE char: the accumulator stores
%                       systems in an array indexed by [store.sysID], which
%                       flattens a multi-char string into separate elements. Use
%                       a number (e.g. 1) or a single letter (e.g. 'C'); carry a
%                       longer human-readable label in 'name' instead.
%       'numLocal' 1x1  local-parameter count carried on L    (default []: set it
%                                       or supply 'localParams'; passing both
%                                       must agree.
%       'numCommon' 1x1 COMMON-parameter count carried on L   (default []: set
%                                       it or supply 'commonParams'; passing both
%                                       must agree).
%                                       nc is SHARED by every system, so all
%                                       loaders in one adjustment must agree.
%                                       LSQVCE.ADJUSTBLOCK and LSQVCE.ADJUSTVCE
%                                       DERIVE it from here -- there is no nc
%                                       argument on them, so it is stated once.
%       'commonParams' RxC cell  parametrisation of the COMMON unknowns, one ROW
%                                       per unknown, C = 14 columns (meaning is
%                                       your convention). Piped to buildDay as its
%                                       3rd arg. nc is DERIVED as size(.,1) and
%                                       carried as L.nc. SHARED across systems, so
%                                       every loader must pass the SAME value.
%                                       (default {})  (alias: 'common')
%       'localParams'  RxC cell  parametrisation of the LOCAL unknowns, one ROW per
%                                       unknown, C = 14 columns. Piped to buildDay
%                                       as its 4th arg. nl is DERIVED as size(.,1).
%                                       (default {})  (alias: 'local')
%       'processingOptions' 1x1 struct  free-form processing options, piped to
%                                       buildDay as its 5th arg (carried as
%                                       L.procOpt). (default struct()) (alias: 'procOpt')
%       'days'     1xD  default day list carried on L                (default [])
%       'name'     1xk  display name                                 (default '')
%       'cache'    1x1  logical, keep a one-day cache               (default true)
%       'epochCovariance' 1x1 function_handle  @(day,obsType) -> a struct exposing
%                                       the FIXED per-epoch scaling D directly, in
%                                       the form LSQVCE.RESIDUALACV wants:
%                                       .arc(i).rows/.epochs and .sigma (scalar
%                                       stream) or .Cep (nEpoch x d x d). Optional;
%                                       lets LSQVCE.ADJUSTVCE re-estimate the
%                                       temporal law by exact separable de-scaling
%                                       instead of the whitened-residual path.
%                                                                    (default [])
%
%   RUNTIME COVARIANCE OVERRIDE (opt.covOverride)
%     The batch merges a runtime override into the base cov spec, so an outer
%     loop can reshape the model per pass WITHOUT rebuilding the loader. Two
%     shapes, later wins:
%       opt.covOverride.(obsType)                flat, all satellites
%       opt.covOverride.bySys.(satKey).(obsType) per-satellite (keyed on the
%                                                injected opt.sysID via
%                                                matlab.lang.makeValidName). This
%                                                is how LSQVCE.ADJUSTVCE
%                                                threads a distinct sigma_g / AR
%                                                law into each stream.
%
%   OUTPUT
%     L    1x1  struct  A Loader object:
%            .closure  1x1 function_handle  @(day,obsType,arcID,opt) the 4-arg
%                                           callback the accumulator drives
%                                           (handles the 'meta' sentinel too)
%            .meta     1x1 function_handle  @() -> .obsTypes, .arcs
%            .batch    1x1 function_handle  @(day,obsType,arcID,opt) -> one arc's
%                                           batch (.A .B .l .t .Cfcn .Wfcn
%                                           .filterMemory)
%            .day      1x1 function_handle  @(day,obsType) -> the cached whole day
%            .sysID .nl .obsTypes .cov .name .days   the carried settings
%            .nc       1x1  common-parameter count DERIVED from 'common' rows
%                           (0 if 'common' not given); nl likewise from 'local'
%            .common .local .procOpt   the parametrisation arrays and processing
%                                      options piped to buildDay (3rd/4th/5th args)
%            .epochCov 1x1 function_handle | []   optional per-epoch cov hook
%            .plotPostSolve .plotPostModel   containers.Map obsType->@(ctx)->fig,
%                                            or empty Maps if no hooks given
%            .isLoader 1x1 logical  true (marks it as a Loader for adjustBlock)
%
%   THE LOADER CONTRACT (what the engines actually require)
%     A Loader is consumed by LSQVCE.BLOCKACCUMULATOR, LSQVCE.ADJUSTBLOCK and
%     LSQVCE.ADJUSTVCE. The ONLY hard requirement is the 4-arg data
%     callback @(day,obsType,arcID,opt) that also answers the 'meta'/'day'
%     sentinels -- .closure above. Every OTHER field is OPTIONAL: the engines
%     read them defensively (isfield guards on .epochCov, .plotPostSolve,
%     .plotPostModel, .days, .common/.local/.procOpt, per-field .cov), so a
%     Loader missing any of them simply does not use that feature -- it never
%     errors. This is deliberate so that:
%       (a) makeLoader can add new optional fields without breaking old loaders;
%       (b) a HAND-WRITTEN raw 4-arg closure (e.g. GRACELOADER bound over a cfg,
%           or LSQVCE.SYNTHLOADER) is a FIRST-CLASS loader -- you can pass it
%           straight to acc.addSystem or to adjustBlock's struct-array form
%           (.sysID/.nl/.loader) without wrapping it here. Such a closure carries
%           none of the L.* struct fields, so it fires no plot hooks and exposes
%           no parametrisation pipe; wrap it with makeLoader ONLY when you want
%           those conveniences. Both paths are supported and tested.
%
%   EXAMPLE
%     cov.pos = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',0.7);
%     L = lsqvce.makeLoader(@(d,ot) myBuildDay(d,ot), ...
%             'observationTypes',{'pos'}, 'covariance',cov, 'arcs',2, 'systemId','C', 'numLocal',4);
%     R = lsqvce.adjustBlock({L}, 1:3);
%
%   ERRORS
%     Raised while BUILDING the loader (your declaration is wrong):
%       lsqvce:makeLoader:noObsTypes   'observationTypes' missing or empty
%       lsqvce:makeLoader:noCov        'covariance' missing, or no entry for an
%                                      observation type
%       lsqvce:makeLoader:noArcs       'arcs' missing and not derivable
%       lsqvce:makeLoader:badArcs      an arc resolver returned something that
%                                      is not unique positive integers (WP-52)
%       lsqvce:makeLoader:badArcId     a day's .arcId does not name every
%                                      .arcIdx entry exactly once (WP-52)
%       lsqvce:makeLoader:localPlan    'localsPerArc' without 'days', combined
%                                      with 'fromDay', or a bad count (WP-52)
%       lsqvce:makeLoader:paramType    commonParams/localParams is not a 2-D cell
%       lsqvce:makeLoader:paramCols    ... or does not have the required columns
%       lsqvce:makeLoader:ncConflict   numCommon disagrees with the commonParams
%                                      row count -- state it ONCE, in the table
%       lsqvce:makeLoader:nlConflict   the same for numLocal/localParams
%       lsqvce:makeLoader:badPlotHandle  a plot hook is neither a function handle
%                                      nor a containers.Map keyed by obsType
%
%     Raised while READING a day (your buildDay is wrong):
%       lsqvce:makeLoader:buildDay     your buildDay itself threw; the original
%                                      error is attached as the cause
%       lsqvce:makeLoader:badDay       a required field (.A .B .l .t) is missing
%       lsqvce:makeLoader:dayShape     A, B, l, t do not share a row count
%       lsqvce:makeLoader:badArc       an .arcIdx entry is not a valid row-index
%                                      vector for this day
%       lsqvce:makeLoader:noD          a multi-axis stream has no per-epoch D:
%                                      supply an 'epochCovariance' hook
%       lsqvce:makeLoader:badDayCep    the hook returned a Cep that is neither
%                                      nEpoch x d nor nEpoch x d x d
%
%   See also LSQVCE.ADJUSTBLOCK, LSQVCE.SYNTHLOADER, LSQVCE.COVMODEL,
%            LSQVCE.BLOCKACCUMULATOR, LSQVCE.ADJUSTVCE
%
%   Updated 2026-08-21 (CEST, UTC+02:00), WP-52: ragged schedules are
%   NATIVE -- per-day arc resolution ('arcs' as @(day,obsType) or
%   'fromDay'), arcs addressed by .arcId with absent ids as data gaps,
%   and 'localsPerArc' owning the per-arc local bookkeeping (nl derived,
%   offsets handed to buildDay). All additive. With assistance from
%   Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (c) c1 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   Created 2026-07-16 12:20 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
%   Updated 2026-07-17 (CEST, UTC+02:00): documented the Loader contract
%   (optional fields; raw 4-arg closures are first-class). Claude Opus 4.8.
%   Updated 2026-07-17 (CEST, UTC+02:00): expanded terse option names to
%   canonical long forms (observationTypes/covariance/systemId/numLocal/
%   commonParams/localParams/processingOptions/epochCovariance); the old terse
%   names still work for one version (warn once). Claude Opus 4.8.
arguments
    buildDay
    opt.ObservationTypes = {}
    opt.Covariance (1,1) struct = struct()
    opt.Arcs = 1
    opt.LocalsPerArc = []
    opt.SystemId = 1
    opt.NumLocal double {mustBeScalarOrEmpty} = []
    opt.NumCommon double {mustBeScalarOrEmpty} = []
    opt.Days double = []
    opt.Name (1,1) string = ""
    opt.Cache (1,1) logical = true
    opt.EpochCovariance = []
    opt.CommonParams cell = {}
    opt.LocalParams cell = {}
    opt.ProcessingOptions (1,1) struct = struct()
    opt.PlotPostSolve = []
    opt.PlotPostModel = []
end

    narginchk(1, inf);
    if ~isa(buildDay, 'function_handle')
        error('lsqvce:makeLoader:buildDay', 'buildDay must be a function handle.');
    end
    % WP-39 stage (c) c1: typed `arguments` signature (canonical
    % UpperCamelCase; legacy casings keep matching -- MATLAB NV
    % semantics; enum VALUES case-sensitive lowercase).
    o = struct('observationTypes', {opt.ObservationTypes}, 'covariance', opt.Covariance, 'arcs', opt.Arcs, 'localsPerArc', opt.LocalsPerArc, 'systemId', opt.SystemId, 'numLocal', opt.NumLocal, 'numCommon', opt.NumCommon, 'days', opt.Days, 'name', char(opt.Name), 'cache', opt.Cache, 'epochCovariance', opt.EpochCovariance, 'commonParams', {opt.CommonParams}, 'localParams', {opt.LocalParams}, 'processingOptions', opt.ProcessingOptions, 'plotPostSolve', opt.PlotPostSolve, 'plotPostModel', opt.PlotPostModel);

    % ---- map long option names onto the short internal fields -------------
    % Canonical long names only. The resolved value is written onto the SHORT
    % field name so the rest of makeLoader (which reads o.obsTypes/o.cov/
    % o.sysID/... and the local vars) stays unchanged.
    o.obsTypes = o.observationTypes;
    o.cov      = o.covariance;
    o.sysID    = o.systemId;
    o.nl       = o.numLocal;
    o.nc       = o.numCommon;
    o.epochCov = o.epochCovariance;
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

    % ---- parametrisation arrays: validate, derive nc/nl -------------------
    % Each is a cell array with ONE ROW per unknown and 14 columns of meaning
    % (the meaning is the caller's convention; makeLoader only checks shape and
    % pipes the arrays through to buildDay). Empty ({}) means "not used" -- the
    % builder then defines the design itself and nc/nl come from 'nl'/adjustBlock.
    NCOLPARAM = 14;
    checkParamArray(o.common, 'common', NCOLPARAM);
    checkParamArray(o.local,  'local',  NCOLPARAM);

    % nc DERIVED from common rows (the common field is SHARED across systems, so
    % every loader must pass the SAME 'common'); nl DERIVED from local rows.
    ncDerived = size(o.common, 1);        % 0 if common == {}
    nlDerived = size(o.local,  1);        % 0 if local  == {}

    ncFinal = o.nc;
    if ~isempty(o.common)                 % commonParams given -> it SETS nc
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
    if ~isempty(o.local)                  % local given -> it SETS nl
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

    L.isLoader = true;
    L.sysID    = o.sysID;
    L.nl       = nlFinal;
    L.obsTypes = obsTypes;
    L.cov      = o.cov;
    L.name     = char(o.name);
    L.epochCov = o.epochCov;          % optional @(day,obsType)->epoch-cov struct
    % Parametrisation piped to buildDay + the DERIVED counts (nc is 0 when
    % 'common' is not used; it is the SHARED common count and every loader must
    % agree on it -- pass nc = L.nc to adjustBlock, or the explicit scalar).
    L.common   = o.common;
    L.local    = o.local;
    L.procOpt  = o.procOpt;
    L.nc       = ncFinal;
    % Plot hooks normalised to a containers.Map obsType -> handle (a single
    % handle is broadcast to every obsType of this loader). Empty -> empty Map.
    L.plotPostSolve = normalisePlotSpec(o.plotPostSolve, obsTypes);
    L.plotPostModel = normalisePlotSpec(o.plotPostModel, obsTypes);
    if ~isempty(o.days), L.days = o.days; end

    % WP-52: L.meta() answers for a DAY. Called without one it uses the first
    % declared day (or 1), which is what the old day-blind accessor meant.
    L.meta    = @(varargin) metaOf(obsTypes, o.arcs, ...
                                   metaDay(varargin, o.days), L.day);
    L.localPlan = plan;               % WP-52: [] unless 'localsPerArc' declared
    L.day     = @(day, obsType) dayOf(buildDay, cache, o.cache, day, obsType, ...
                                      o.common, o.local, o.procOpt, plan);
    L.batch   = @(day, obsType, arcID, opt) batchOf(L.day, o.cov, day, obsType, arcID, opt);
    L.closure = @(day, obsType, arcID, opt) dispatch(L, o.arcs, day, obsType, arcID, opt);
end

%==========================================================================
function tf = isPlotSpec(v)
%ISPLOTSPEC  A valid plot-hook spec: empty, a function handle (broadcast to all
%   obsTypes), or a containers.Map / struct keyed obsType -> handle.
    tf = isempty(v) || isa(v, 'function_handle') || ...
         isa(v, 'containers.Map') || isstruct(v);
end

function M = normalisePlotSpec(spec, obsTypes)
%NORMALISEPLOTSPEC  Turn a plot-hook spec into a containers.Map obsType->handle.
%   Empty -> empty Map. A single handle -> that handle for every obsType. A Map
%   or struct -> validated pass-through (a handle value is required).
    M = containers.Map('KeyType','char', 'ValueType','any');
    if isempty(spec), return, end
    if isa(spec, 'function_handle')
        for k = 1:numel(obsTypes), M(obsTypes{k}) = spec; end
        return
    end
    if isa(spec, 'containers.Map')
        ks = spec.keys;
        for i = 1:numel(ks)
            requireHandle(spec(ks{i}), ks{i});
            M(char(ks{i})) = spec(ks{i});
        end
        return
    end
    if isstruct(spec)
        f = fieldnames(spec);
        for i = 1:numel(f)
            requireHandle(spec.(f{i}), f{i});
            M(f{i}) = spec.(f{i});
        end
    end
end

function requireHandle(v, key)
%REQUIREHANDLE  Assert v is a function handle (a plot hook for KEY); errors otherwise.
    if ~isa(v, 'function_handle')
        error('lsqvce:makeLoader:badPlotHandle', ...
            'Plot hook for ''%s'' must be a function handle @(ctx)->figHandle.', ...
            char(string(key)));
    end
end

function checkParamArray(C, name, ncol)
%CHECKPARAMARRAY  Validate a parametrisation cell array: empty, or a 2-D cell
%   with exactly NCOL columns (one row per unknown). Empty is allowed (means the
%   builder defines the design itself).
    if isempty(C), return, end
    if ~iscell(C) || ~ismatrix(C)
        error('lsqvce:makeLoader:paramType', ...
            '''%s'' must be a 2-D cell array (one row per unknown).', name);
    end
    if size(C, 2) ~= ncol
        error('lsqvce:makeLoader:paramCols', ...
            '''%s'' must have %d columns (got %d).', name, ncol, size(C, 2));
    end
end

function out = dispatch(L, arcs, day, obsType, arcID, opt)
%DISPATCH  Route a loader call: return metaOf(...) for the 'meta' probe, else
%   delegate to L.batch(day,obsType,arcID,opt). OUTPUT out = meta struct | batch struct.
    if nargin < 6 || isempty(opt), opt = struct(); end
    if (ischar(day) || isstring(day)) && strcmp(day, 'meta')
        % WP-52: the probe carries the DAY (every engine passes it), so a
        % ragged schedule is answered natively -- see metaOf/arcVec.
        out = metaOf(L.obsTypes, arcs, obsType, L.day);
        return
    end
    out = L.batch(day, obsType, arcID, opt);
end

function m = metaOf(obsTypes, arcs, day, dayFcn)
%METAOF  Build the meta struct {obsTypes, arcs{:}} the accumulator expects.
%   INPUT  obsTypes 1xK cell; arcs the spec; day the day the probe asks about
%          (the 'meta' sentinel carries it in the obsType slot); dayFcn the
%          loader's day accessor, used only by the 'fromDay' resolver.
%   OUTPUT m 1x1 struct  .obsTypes, .arcs{k}
    m.obsTypes = obsTypes;
    m.arcs     = cell(1, numel(obsTypes));
    for k = 1:numel(obsTypes)
        m.arcs{k} = arcVec(arcs, obsTypes{k}, day, dayFcn);
    end
end

function v = arcVec(arcs, obsType, day, dayFcn)
%ARCVEC  Resolve the arc-ID vector for one obsType and ONE DAY.
%   Static forms (count / vector / per-obsType struct) ignore the day, as
%   before. WP-52 adds two ragged forms:
%     function_handle  @(day,obsType) -> arc ids     (no data read)
%     'fromDay'        read the day and take its .arcId, or 1:numel(.arcIdx)
%   A per-obsType struct may hold any of these per field.
%   OUTPUT v 1xA double  the arc IDs present on THIS day.
    if isstruct(arcs)
        if ~isfield(arcs, obsType)
            error('lsqvce:makeLoader:noArcs', 'arcs has no field for ''%s''.', obsType);
        end
        v = resolveArcSpec(arcs.(obsType), obsType, day, dayFcn);
    else
        v = resolveArcSpec(arcs, obsType, day, dayFcn);
    end
end

function v = resolveArcSpec(spec, obsType, day, dayFcn)
%RESOLVEARCSPEC  One arc spec -> the arc-ID row vector for one day.
    if isa(spec, 'function_handle')
        if nargin(spec) == 1
            v = spec(day);
        else
            v = spec(day, obsType);
        end
        v = checkArcIds(v, obsType, day, 'the ''arcs'' function');
    elseif (ischar(spec) || isstring(spec)) && strcmpi(char(spec), 'fromday')
        if isempty(dayFcn)
            error('lsqvce:makeLoader:noArcs', ...
                  '''arcs'',''fromDay'' needs the loader''s day accessor.');
        end
        D = dayFcn(day, obsType);
        if isfield(D, 'arcId') && ~isempty(D.arcId)
            v = D.arcId(:).';
        else
            v = 1:numel(D.arcIdx);
        end
        v = checkArcIds(v, obsType, day, '''arcs'',''fromDay''');
    elseif isscalar(spec) && isnumeric(spec)
        v = 1:spec;
    else
        v = spec(:).';
    end
end

function plan = buildLocalPlan(spec, arcs, obsTypes, days)
%BUILDLOCALPLAN  Assign every (day, arc) its local-block ordinal and column
%   offset, and DERIVE nl (WP-52).
%   The arc set of a day is the UNION over the loader's obsTypes: an arc seen
%   by two observation types is ONE arc state and therefore ONE local block.
%   Requires a schedule resolvable WITHOUT reading data -- static arcs or an
%   @(day[,obsType]) resolver -- because nl must be known at construction.
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
    if (ischar(arcs) || isstring(arcs)) && strcmpi(char(arcs), 'fromday')
        error('lsqvce:makeLoader:localPlan', ...
            ['''localsPerArc'' cannot be combined with ''arcs'',''fromDay'': ' ...
             'nl must be known WITHOUT reading data. Declare the schedule ' ...
             'with an @(day,obsType) resolver instead.']);
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
%ARCINFOOF  The local-block slice buildDay gets for ONE day (WP-52).
%   OUTPUT info 1x1 struct | []  .arcId .offset .count .ordinal (1xA each),
%          .nl (total over the adjustment), .offsetOf @(arcId)->column offset.
    info = [];
    if isempty(plan), return, end
    i = find(plan.days == day, 1);
    if isempty(i), return, end
    info = struct('arcId', plan.arcId{i}, 'offset', plan.offset{i}, ...
                  'count', plan.count{i}, 'ordinal', plan.ordinal{i}, ...
                  'nl', plan.nl, 'offsetOf', []);
    ids = plan.arcId{i};  offs = plan.offset{i};
    info.offsetOf = @(a) offs(find(ids == a, 1));
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

function D = dayOf(buildDay, cache, useCache, day, obsType, common, local, procOpt, plan)
%DAYOF  Return the built+validated day struct for (day,obsType), using a
%   single-slot cache when useCache. OUTPUT D 1x1 struct  the day data.
    key = sprintf('%d|%s', day, obsType);
    if useCache && isKey(cache, key)
        D = cache(key);  return
    end
    arcInfo = arcInfoOf(plan, day);           % WP-52: [] unless 'localsPerArc'
    D = callBuildDay(buildDay, day, obsType, common, local, procOpt, arcInfo);
    D = validateDay(D, day, obsType);
    if useCache
        remove(cache, cache.keys);   % keep at most ONE resident day
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
%VALIDATEDAY  Assert a built day has fields A,B,l,t with matching row counts and
%   fill a default single-arc arcIdx; errors otherwise. OUTPUT D the checked day.
    for f = {'A','B','l','t'}
        if ~isfield(D, f{1})
            error('lsqvce:makeLoader:badDay', ...
                'buildDay(%d,''%s'') must return field .%s.', day, obsType, f{1});
        end
    end
    m = size(D.A, 1);
    if numel(D.l) ~= m || numel(D.t) ~= m || size(D.B,1) ~= m
        error('lsqvce:makeLoader:dayShape', ...
            'buildDay(%d,''%s''): A,B,l,t must share the same number of rows.', day, obsType);
    end
    if ~isfield(D, 'arcIdx') || isempty(D.arcIdx)
        D.arcIdx = {1:m};                       % whole day = one arc
    end
    % WP-52: the optional .arcId names the arcs this day carries, so a ragged
    % schedule needs no padding. It must name EVERY entry of .arcIdx exactly
    % once -- a mismatch would silently mis-address whole arcs.
    if isfield(D, 'arcId') && ~isempty(D.arcId)
        ids = D.arcId(:).';
        if numel(ids) ~= numel(D.arcIdx) || ...
                numel(unique(ids)) ~= numel(ids) || ...
                ~all(isfinite(ids) & ids >= 1 & ids == round(ids))
            error('lsqvce:makeLoader:badArcId', ...
                ['buildDay(%s,''%s'') returned %d .arcId entr(ies) for %d ' ...
                 'arcIdx entr(ies); .arcId must name every arc exactly once ' ...
                 'with unique positive integers.'], ...
                 string(day), obsType, numel(ids), numel(D.arcIdx));
        end
    end
end

function out = batchOf(dayFcn, cov, day, obsType, arcID, opt)
%BATCHOF  Slice one arc (arcID) out of a day into a batch struct (A,B,l,...),
%   applying the covariance/whitening hook cov. OUTPUT out 1x1 struct the batch.
    D = dayFcn(day, obsType);
    % WP-52: a day that declares .arcId is addressed BY ID -- an id this day
    % does not carry is a DATA GAP (return empty; every engine skips it), not
    % an error. Without .arcId the addressing stays positional and an
    % out-of-range id remains an error: that guard catches a mis-declared
    % STATIC schedule, which is a different mistake from a ragged one.
    if isfield(D, 'arcId') && ~isempty(D.arcId)
        pos = find(D.arcId(:).' == arcID, 1);
        if isempty(pos)
            out = [];
            return
        end
    else
        if arcID < 1 || arcID > numel(D.arcIdx)
            error('lsqvce:makeLoader:badArc', ...
                'Arc %d does not exist (day %d, obsType %s has %d arcs).', ...
                arcID, day, obsType, numel(D.arcIdx));
        end
        pos = arcID;
    end
    r = D.arcIdx{pos};
    out.A = D.A(r,:);
    out.B = D.B(r,:);
    out.l = D.l(r);
    out.t = D.t(r);

    cm = cov.(obsType);
    if isfield(opt,'covOverride')
        cm = applyCovOverride(cm, opt.covOverride, obsType, opt);
    end

    % Separable model needs the FIXED per-epoch scaling D for THIS arc. Source,
    % in priority order: an already-supplied cm.Cep; a per-day D.Cep on the day
    % struct (sliced to the arc's epochs); or the loader's epochCov hook. The
    % day-struct route keeps everything in one deterministic buildDay.
    if strcmpi(lsqvce.util.getfielddef(cm, 'type', ''), 'separable') ...
            && ~isfield(cm, 'Cep')
        d = lsqvce.util.getfielddef(cm, 'd', []);
        if isfield(D, 'Cep') && ~isempty(D.Cep)
            cm.Cep = sliceEpochCov(D.Cep, D.arcIdx{arcID}, d, r);
        elseif isfield(D, 'nAxis') && ~isempty(d)
            % epochs of this arc are its rows folded by d (epoch-major)
            cm.Cep = [];   % leave empty -> covModel errors clearly
        end
    end

    [out.Cfcn, out.Wfcn, info] = lsqvce.covModel(cm, out.t, opt);
    out.filterMemory = info.filterMemory;

    % Multi-axis (epoch-major) streams must tell the accumulator the axis count
    % so it batches on whole epochs and the whitener treats axes correctly.
    % Priority: an explicit D.nAxis from buildDay, else the cov model's own d.
    if isfield(D, 'nAxis') && ~isempty(D.nAxis)
        out.nAxis = D.nAxis;
    elseif isfield(info, 'nAxis') && info.nAxis > 1
        out.nAxis = info.nAxis;
    end
end

%==========================================================================
function Cep = sliceEpochCov(dayCep, arcRows, d, ~)
%SLICEEPOCHCOV  Select the day's per-epoch scaling for one arc's epochs.
%   dayCep is per-DAY-epoch: nEpochDay x d  or  nEpochDay x d x d. arcRows are
%   the arc's day rows in epoch-major order; its day-epochs are the unique
%   ceil(arcRows/d). Returns Cep for exactly those epochs, in order.
    if isempty(d) || d < 1
        error('lsqvce:makeLoader:noD', ...
              'A separable model needs cm.d to slice the per-epoch covariance.');
    end
    ep = unique(ceil(arcRows(:)/d), 'stable');
    if ismatrix(dayCep) && size(dayCep,2) == d
        Cep = dayCep(ep, :);
    elseif ndims(dayCep) == 3 && size(dayCep,2) == d && size(dayCep,3) == d
        Cep = dayCep(ep, :, :);
    else
        error('lsqvce:makeLoader:badDayCep', ...
              'D.Cep must be nEpochDay x %d or nEpochDay x %d x %d.', d, d, d);
    end
end

%==========================================================================
function cm = applyCovOverride(cm, ovr, obsType, opt)
%APPLYCOVOVERRIDE  Merge a runtime covariance override into the base spec.
%   Supports two shapes, resolved in this order (later wins):
%     (1) ovr.(obsType)               -- flat, applies to every satellite
%     (2) ovr.bySys.(satKey).(obsType)-- per-satellite, keyed by opt.sysID
%   The per-satellite layer is how LSQVCE.ADJUSTVCE threads a distinct
%   sigma_g / AR law into each (satellite,obsType) stream while all loaders
%   share one UserOpt. Both layers are optional and fully backward compatible.
    if isfield(ovr, obsType)
        cm = lsqvce.util.mergeStruct(cm, ovr.(obsType));
    end
    if isfield(ovr, 'bySys') && isfield(opt, 'sysID')
        satKey = matlab.lang.makeValidName(char(string(opt.sysID)));
        if isfield(ovr.bySys, satKey) && isfield(ovr.bySys.(satKey), obsType)
            cm = lsqvce.util.mergeStruct(cm, ovr.bySys.(satKey).(obsType));
        end
    end
end
