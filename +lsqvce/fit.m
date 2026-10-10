function R = fit(data, opts, nv)
%LSQVCE.FIT  Fit a least-squares adjustment (optionally with VCE) to declared data.
%
%   R = lsqvce.fit(data)                        every default (one GLS solve)
%   R = lsqvce.fit(data, opts)                  a lsqvce.Options object
%   R = lsqvce.fit(data, Name, Value, ...)      the Options properties, straight on fit
%   R = lsqvce.fit(data, opts, Name, Value)     both: the pairs override opts
%
%   THE ONE VERB. You describe the problem with a plain struct `data`
%   (manual, chapter 3), say what to do with lsqvce.Options (chapter 4), and
%   get a lsqvce.Results back (chapter 5):
%
%       R = lsqvce.fit(data, Estimate="noise", Pool="obsType", OnPass=@myInspect);
%       R.report;  R.plot("psd");  [C, L] = R.toTable(data);
%
%   THREE DATA KINDS, told apart by the fields of `data`, all on ONE accumulator:
%     (A) streaming   data.commonParams, data.days, data.systems(k) with
%                     buildDay / covariance / localParams / priors /
%                     processingOptions and the optional arc schedule `arcs`
%                     (v6: buildDay returns one CELL per arc -- A, B, l, t,
%                     Cep as 1 x nArc cells -- and the cells are the schedule
%                     unless `arcs` declares one); optional data.priors,
%                     data.paramColumns
%     (B) in memory   data.A (m x u), data.l (m x 1), optional data.p (m x 1
%                     weights); optional data.priors. A FULL m x m SPD weight
%                     matrix is kind C's: return it from batchFcn (R6 WP-101)
%     (C) batched     data.batchFcn @(k) -> [A, l, p] (or @(k, x) -> [J, w, p] with
%                     data.x0 for Gauss-Newton), data.numUnknowns, data.numBatches;
%                     optional data.priors
%   WHAT RUNS (opts.Estimate / opts.Robust / data.x0):
%     Estimate "none"    one generalised solve under the declared laws (adjustBlock)
%     Estimate "noise"   the joint variance-component estimator (adjustVCEjoint):
%                        one spectral shape per Pool group, one level per stream,
%                        per-arc factors when ArcFactors; prior groups with
%                        estimate "scale" get their sigma by Helmert / Koch-Kusche
%     Robust "huber" | "igg3" | "bisquare"   IRLS on the same loop: alone it is
%                        adjustRobust; with Estimate "noise" the robust pass runs
%                        INSIDE the joint loop (R6 WP-105), R.noise.robust carries
%                        the final factors
%     data.x0            Gauss-Newton with damping (adjustNonlinear)
%   REGULARISATION is declared, never an option: data.priors(j) on the common
%   parameters and data.systems(k).priors(l) on the local ones -- one struct per
%   group, by `group` (a name pattern) or `L` (in the columns lsqvce.layout
%   announces), x0, sigma, from, name, estimate. A free network's datum is the
%   prior on its null space.
%   WATCHING A RUN: Options.OnPass fires after every pass with the pass state,
%   the live accumulator, `data` and every system's processingOptions.
%
%   INPUT
%     data   1x1 struct    the declaration (kind A, B or C as above)
%     opts   1x1 lsqvce.Options   default lsqvce.Options()
%     Name, Value          any property of lsqvce.Options
%   OUTPUT
%     R      1x1 lsqvce.Results   x, sigma0, sol (Cxx, dof, vtv), diag, acc, history,
%                                 noise (C, Cgroups, stream, model, priors), stopReason,
%                                 meta (options as resolved, warnings, seeds, version);
%                                 methods report, summary, convergence, residuals,
%                                 contribution, compare, plot, toTable, coefficients
%
%   ERRORS
%     Your declaration is wrong (raised BEFORE any data is read):
%       lsqvce:fit:unknownField    a field name the contract does not know
%                                  (names the field; typos are errors on
%                                  purpose -- and so are the retired fields:
%                                  common, numCommon, numLocal,
%                                  observationTypes, systemId, cache, loader,
%                                  plotPostSolve, plotPostModel, since v4
%                                  localsPerArc (declare 'arc' rows in
%                                  localParams instead), and since v5.3
%                                  epochCovariance (the day returns Cep))
%       lsqvce:fit:noSystems       data.systems missing, not a struct, or empty
%       lsqvce:fit:noCommon        data.commonParams missing, or not a 2-D
%                                  cell of template rows
%       lsqvce:fit:badScope        a scope-token column 1 does not carry:
%                                  commonParams needs 'common'|'day' (raised
%                                  here), localParams needs
%                                  'arc'|'day'|'system' (raised by makeLoader
%                                  and relabelled)
%       lsqvce:fit:badParamColumns data.paramColumns is malformed: not a
%                                  struct, an unknown field, only one of
%                                  start / end, or a column the table lacks
%       lsqvce:fit:missingField    a required system field is absent, or
%                                  covariance carries no observation type
%       lsqvce:fit:noDays          no day list anywhere
%       lsqvce:fit:badId           id is not a numeric scalar or single char
%       lsqvce:fit:duplicateSystem two systems share an id (explicit or
%                                  index default)
%       lsqvce:fit:builderKind     buildDay is not a function handle
%       lsqvce:fit:badArcs         arcs is a string (the v5 'fromDay' is retired:
%                                  omit arcs -- the cells are the arcs), or a
%                                  resolver returned something that is not
%                                  unique positive integers
%       lsqvce:fit:<makeLoader id> a makeLoader construction error, relabelled
%                                  (e.g. paramType, badArcs, localPlan);
%                                  the original is attached as the cause
%     Your choice is not available here:
%       lsqvce:fit:notYet          Robust ~= "none" on a NONLINEAR batch problem
%                                  (x0): IRLS inside Gauss-Newton is not written;
%                                  linearise first. (Robust with Estimate "noise"
%                                  RUNS since R6 WP-105: the robust pass inside
%                                  the joint loop)
%       lsqvce:fit:modelKind       Estimate = "noise" on kinds B / C, which have
%                                  no streams to estimate a law for
%       lsqvce:batchSystem:weightMatrix  a batch returned a weight argument that
%                                  is neither m x 1 nor an m x m SPD matrix. A
%                                  full SPD P IS accepted on kind C and whitened
%                                  by its Cholesky factor (R5 M5 5b') -- raised
%                                  by lsqvce.engine.batchSystem
%       lsqvce:fit:priorNotActive  (warning, R5-D8) a prior with from > 1 in a
%                                  single-pass run is not applied
%       lsqvce:fit:badPrior        (R5 M5 5b) data.priors is malformed (unknown
%                                  field, L not p x nc, sigma not a positive scalar,
%                                  x0 of the wrong length)
%       lsqvce:fit:badEstimate     (R5 M3 3d) covariance.<type>.estimate is not
%                                  "none" | "scale" | "law"
%       lsqvce:fit:lawOnFixedStructure  (R5 M4 4b) estimate "law" on a stream whose
%                                  sigma is per-row (a fixed structure: level only)
%       lsqvce:fit:varUnderNoise   (R7-07) a 'var' law under Estimate = "noise": the
%                                  joint sweep has no vector-AR path
%       lsqvce:fit:badMatrix       kind B shapes disagree (A m x u, l m x 1,
%                                  p m x 1)
%       lsqvce:fit:missingField    kind C lacks batchFcn / numBatches, or the
%                                  linear form lacks numUnknowns
%       lsqvce:fit:uCount          kind C carries BOTH x0 and numUnknowns --
%                                  u is stated once (numel(x0))
%     Your buildDay is wrong (raised WHILE READING, by the engine, ids
%     unchanged): lsqvce:makeLoader:badDay / dayShape / badArc / badArcId /
%     buildDay / badDayCep / axisConflict / cepTwice (v5.3: the day's Cep and
%     axis count), lsqvce:covModel:cepOnVar, lsqvce:accum:arcBadNc / arcBadNl;
%     lsqvce:makeLoader:noCep is a once-per-run warning (a multi-axis day
%     without Cep is composed with D = I)
%
%   See also LSQVCE.OPTIONS, LSQVCE.RESULTS, LSQVCE.LAYOUT, LSQVCE.ENGINE.DECLARATION,
%            LSQVCE.ENGINE.MAKELOADER, LSQVCE.ENGINE.ADJUSTBLOCK, LSQVCE.ENGINE.ADJUSTVCEJOINT,
%            LSQVCE.ENGINE.ADJUSTROBUST, LSQVCE.ENGINE.ADJUSTNONLINEAR
%
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34, R7-37):
%   ERRORS lists badParamColumns; the joint seed recorded is the one used.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-116 part 1: a 'var' law under Estimate = "noise" is refused by
%   name (lsqvce:fit:varUnderNoise, R7-07). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%
%   Created 2026-08-31 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Fable 5 -- R4 Level L2, WP-56 (decision R4-D4).
%   Help rewritten 2026-09-04 for v5.2 (R5 M5 / R5-D8), model Claude Fable 5.1:
%   the v5 surface (Estimate, priors, paramColumns, OnPass, layout, toTable).
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-105: Robust together with
%   Estimate "noise" runs (the robust pass inside the joint loop; the notYet
%   refusal stays for x0 problems only). With assistance from Claude
%   (Anthropic), model Claude Fable 5.1.
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-102 (v6): the day contract is
%   ONE CELL PER ARC and the arc schedule is optional. With assistance from
%   Claude (Anthropic), model Claude Fable 5.1.
%   Updated 2026-09-08 (CEST, UTC+02:00), R6 WP-99: R.noise.model{k}.C is the
%   law of the stream's pool group and EMPTY where no law was estimated, instead
%   of falling back to another stream's fitted spectrum. With assistance from
%   Claude (Anthropic), model Claude Opus 5.
%   Updated 2026-09-06 (CEST, UTC+02:00), R6 WP-90: the error information
%   travels with the day (D.Cep); epochCovariance retired; the See also on
%   the units that exist. With assistance from Claude (Anthropic), model
%   Claude Fable 5.1.
arguments
    data (1,1) struct
    opts (1,1) lsqvce.Options = lsqvce.Options()
    nv.?lsqvce.Options
end
    f = fieldnames(nv);
    for k = 1:numel(f), opts.(f{k}) = nv.(f{k}); end

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
    [spec, systems] = lsqvce.engine.declaration(data, opts.CacheDays);   % R5-D8 4: the one declaration step
    days = spec.days;
    [pri, priorGroups] = priorsOf(data, lsqvce.engine.resolveNc(systems), spec);   % R5 M5 5b / D8-6: data.priors
    if strcmp(opts.Estimate, "none") && ~isempty(priorGroups) && any([priorGroups.from] > 1)
        lsqvce.util.log.warnOnce('lsqvce:fit:priorNotActive', ...
            ['%d prior group(s) declare from > 1 and a single solve (Estimate="none") has one pass: ' ...
             'they are not applied. Use Estimate="noise" for a pass schedule.'], sum([priorGroups.from] > 1));
        pri = stackPriors(priorGroups, 1);
    end

    % R5 M2 2b: robust IRLS is a POLICY on the block path -- on the streaming
    % kind it is Estimate "none" with reweighting (H-02: the same accumulator).
    % R6 WP-105: with Estimate "noise" the robust pass runs INSIDE the joint
    % loop (adjustVCEjoint 'Influence'), so it falls through to fitEllmer.
    if ~strcmp(opts.Robust, "none") && strcmp(opts.Estimate, "none")
        a = opts.engineArgs('adjustRobust');
        prov = opts.toStruct();  prov.kind = 'stream';
        R = lsqvce.engine.adjustRobust(systems, days, a.nv{:}, 'Provenance', prov, 'Prior', pri, ...
                                       'OnPassContext', onPassContext(data, systems));
        return
    end

    switch char(opts.Estimate)
        case 'none'
            a = opts.engineArgs('adjustBlock');
            prov = opts.toStruct();  prov.kind = 'stream';
            R = lsqvce.engine.adjustBlock(systems, days, a.nv{:}, 'Provenance', prov, 'Prior', pri);
        otherwise
            % R5 M4 4b (decision R5-D4): ONE estimator. The three iterative
            % models are POLICIES of the joint (Ellmer) estimator:
            %   "ellmer"  per-arc factors on, laws pooled as Options.Pool says
            %   "vce"     one law per stream (Pool "stream"), no arc factors --
            %             the successor of the separate mode, whose AR law was
            %             fitted to the BIASED raw residual ACV (finding H-03)
            %   "pooled"  one law for all streams (Pool "all"), no arc factors
            % (the separate-mode engine adjustVCE was deleted in M4 4c)
            R = fitEllmer(systems, days, opts, pri, onPassContext(data, systems), priorGroups);
    end
end

%==========================================================================
function R = fitMatrix(data, opts)
%FITMATRIX  Kind B: in-memory A, l (, p) -> one-batch driver -> Results.
    rejectUnknown(fieldnames(data), {'A', 'l', 'p', 'priors'}, 'data');
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
    R = runBatch(size(A, 2), batchFcn, 1, [], opts, 'matrix', data);
end

function R = fitBatch(data, opts)
%FITBATCH  Kind C: batched callback (linear, or nonlinear with x0) -> Results.
%   v3: u is stated ONCE -- numUnknowns on the linear form, numel(x0) on the
%   nonlinear form (numUnknowns must then be absent, lsqvce:fit:uCount).
    rejectUnknown(fieldnames(data), {'batchFcn', 'numUnknowns', 'numBatches', 'x0', 'priors'}, 'data');
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
    R = runBatch(u, data.batchFcn, data.numBatches, x0, opts, 'batch', data);
end

function R = runBatch(u, batchFcn, nBatch, x0, opts, kind, data)
%RUNBATCH  The ONE dispatch for kinds B / C -- every case on the block path
%   (R5 M5 5b): linear (gls / IRLS) through the batch adapter, nonlinear
%   through adjustNonlinear; a datum or any regularisation is a declared
%   prior (data.priors), so the flat drivers have no user-surface caller left.
    if ~strcmp(opts.Estimate, "none")
        error('lsqvce:fit:modelKind', ...
            ['Estimate = "noise" needs the streaming kind (variance components are ' ...
             'per stream of a system); the %s kind is Estimate "none".'], kind);
    end
    [pri, priorGroups] = priorsOf(data, u);
    if ~isempty(priorGroups) && any([priorGroups.from] > 1)
        lsqvce.util.log.warnOnce('lsqvce:fit:priorNotActive', ...
            ['%d prior group(s) declare from > 1 and this kind runs one pass (or Gauss-Newton ' ...
             'passes of one linearisation): they are not applied.'], sum([priorGroups.from] > 1));
        pri = stackPriors(priorGroups, 1);
    end
    if isempty(x0)
        R = runBatchOnBlockPath(u, batchFcn, nBatch, opts, kind, pri, data);
        return
    end
    if ~strcmp(opts.Robust, "none")
        error('lsqvce:fit:notYet', ...
              ['Robust = "%s" on a NONLINEAR batch problem (x0) is the one-loop ' ...
               'policy composition of a later level; run the robust pass on the ' ...
               'linearised problem.'], opts.Robust);
    end
    lsqvce.util.log.inventory('reset');
    a = opts.engineArgs('adjustNonlinear');
    prov = opts.toStruct();  prov.kind = kind;
    R = lsqvce.engine.adjustNonlinear(u, batchFcn, nBatch, x0, a.nv{:}, 'Provenance', prov, 'Prior', pri, ...
                                      'OnPassContext', struct('data', data, 'processingOptions', {{}}));
end
function R = fitEllmer(systems, days, opts, pri, opc, priorGroups)
    if nargin < 4, pri = []; end
    if nargin < 5, opc = struct(); end
    if nargin < 6, priorGroups = []; end
%FITELLMER  Model "ellmer": LSQVCE.ENGINE.ADJUSTVCEJOINT over the built systems.
%   Two bridges (R4 WP-56c): (1) the accFcn the driver wants -- a fresh
%   BlockAccumulator per outer pass whose UserOpt carries the current model
%   through LSQVCE.ENGINE.JOINTMODELOPT for every (system, obsType) stream; (2) the
%   driver's out struct -> a Results with the streaming contract.
    nc    = lsqvce.engine.resolveNc(systems);
    specs = accSpecs(systems, days);
    streams = streamKeys(systems, specs);
    a = opts.engineArgs("adjustVCEjoint");
    psdMode = char(opts.Noise.Whitener);
    pool = opts.Pool;  afDefault = opts.ArcFactors;      % R5 M5: explicit policy knobs
    [estMap, afMap] = streamPolicy(systems, specs, streams, afDefault);   % R5 M3 3d (R5-D4)
    % R5 M7 7b: prior groups with estimate "scale" get one Helmert step per
    % pass (from the previous pass's solution) and the prior is re-stacked
    priState = containers.Map('KeyType', 'char', 'ValueType', 'any');
    priState('groups') = priorGroups;  priState('pri') = pri;  priState('rec') = {};
    accFcn = @(model) jointAccumulator(nc, specs, opts, ...
                 withPass(lsqvce.engine.jointModelOpt(model, streams, 'PsdMode', psdMode), model), ...
                 priorFor(priState, model));
    nvj = a.nv;
    out = lsqvce.engine.adjustVCEjoint(accFcn, nvj{:}, 'Estimate', estMap, 'ArcFactors', afMap, ...
                                       'OnPassContext', opc);

    [dsol, batchFcn, nB] = lsqvce.engine.diagnoseAccum(out.acc, out.sol);
    D = lsqvce.engine.diagnose(dsol, batchFcn, 'nBatch', nB, opts.engineArgs("diagnose").nv{:});
    resolved = opts.toStruct();
    resolved.kind = 'stream';
    seeds = struct('jointSeed', out.seed);       % R7-37: the seed USED (pinned when opts.Seed is empty)
    R = lsqvce.Results(out.sol, D, out.acc, dsol, batchFcn, nB, specs, lsqvce.engine.runMeta(seeds, resolved));
    hist = out.history;
    hist.kind  = 'ellmer';
    hist.nIter = out.nIter;
    hist.skips = out.skips;
    % arcKeys is a CELL: braces keep struct() from expanding it into a struct
    % ARRAY (the Manual ch. 3 §2.4 pitfall -- it bit this very line in L3 r1).
    noise = struct('C', out.C, 'Shat', out.Shat, 'sigma2', out.sigma2, 'arcKeys', {out.arcKeys}, ...
                   'Cgroups', out.Cgroups);            % R5 M3 3d: one C per pool group
    % R5 M4 4b: the per-stream summary every model shares -- .stream (keys),
    % .model{k}.sigma = the stream's LEVEL (declared sigma x its estimated
    % factor) and its law (the pool group's C); what the separate mode's
    % Results carried and the gravity demo reads
    [noise.stream, noise.model] = streamModels(systems, specs, out, char(pool), estMap);
    noise.priors = priState('groups');               % R5 M7 7b: the prior groups with their final sigmas
    noise.robust = out.robust;                       % R6 WP-105: influence, scale, factor, keys ([] under "none")
    noise.priorHistory = priState('rec');            % one record per pass that updated them
    reason = "";  if isfield(out, 'userStop') && out.userStop, reason = "userStop"; end
    R.setIterative(hist, noise, out.converged, reason);
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
%   the names LSQVCE.ENGINE.JOINTMODELOPT / LSQVCE.ENGINE.ADJUSTVCEJOINT use.
    keys = {};
    for k = 1:numel(specs)
        ot = cellstr(systems{k}.obsTypes);
        satKey = matlab.lang.makeValidName(char(string(specs(k).sysID)));
        for j = 1:numel(ot)
            keys{end+1} = sprintf('%s|%s', satKey, ot{j});  %#ok<AGROW>
        end
    end
end

function acc = jointAccumulator(nc, specs, opts, uo, pri)
%JOINTACCUMULATOR  One solve-ready accumulator under the current joint model.
    acc = lsqvce.engine.BlockAccumulator(nc, 'Backend', char(opts.Backend), 'KeepLocal', true, ...
              'Parallel', char(opts.Parallel.Mode), 'DistributeDays', opts.Parallel.DistributeDays, ...
              'ParallelThreshold', opts.Parallel.Threshold, 'Verbose', false, 'UserOpt', uo);
    acc.setPrior(pri);
    acc.addSystems(specs);
end

function [pri, groups] = priorsOf(data, nc, spec)
%PRIORSOF  data.priors -> ONE lsqvce.engine.CommonPrior (R5 M5 5b) + the groups.
%   data.priors(k), one struct per GROUP, declared EITHER by .group -- a
%   regular expression matched against the NAME column of commonParams
%   (data.paramColumns.name, default 2); fit builds L with one row per
%   matched template row per instance from the plan it owns -- OR by .L
%   (p x nc in the columns of A, used as is); never both. Then .x0 (p x 1,
%   default 0), .sigma (scalar > 0, default 1: the ONE sigma the rows share --
%   relative weights inside a group go into L), .from (activation pass,
%   default 1), .name (label; default the pattern or 'prior<k>'), .estimate
%   "none" | "scale" ("scale": sigma updated every pass of Estimate="noise" by
%   the Helmert / Koch-Kusche rule, see priorUpdate). Every group is a
%   pseudo-observation set L x ~ x0 with weight 1/sigma^2; the groups are
%   stacked with their rows divided by sigma into one prior of sigma 1
%   (exact). A DATUM (minimum constraints B'x = 0 of a free network) is the
%   group L = B', x0 = 0: exact for any weight when B spans the null space.
%   OUTPUT  pri     [] | CommonPrior     the stacked prior
%           groups  1xG struct (.L .x0 .sigma .name .estimate .from)
    if nargin < 3, spec = []; end
    pri = [];  groups = struct('L', {}, 'x0', {}, 'sigma', {}, 'name', {}, 'estimate', {}, 'from', {});
    if ~isfield(data, 'priors') || isempty(data.priors), return, end
    P = data.priors;
    if ~isstruct(P)
        error('lsqvce:fit:badPrior', 'data.priors must be a struct array (group | L, x0, sigma, from, name, estimate).');
    end
    for k = 1:numel(P)
        f = fieldnames(P(k));
        bad = setdiff(f, {'L', 'x0', 'sigma', 'name', 'estimate', 'from', 'group'});
        if ~isempty(bad)
            error('lsqvce:fit:badPrior', 'data.priors(%d) has an unknown field ''%s''.', k, bad{1});
        end
        hasL = isfield(P(k), 'L') && ~isempty(P(k).L);
        hasG = isfield(P(k), 'group') && ~isempty(P(k).group);
        if hasL && hasG
            error('lsqvce:fit:badPrior', 'data.priors(%d): give EITHER group (a name pattern) OR L, not both.', k);
        end
        nm = sprintf('prior%d', k);
        if hasG
            if isempty(spec)
                error('lsqvce:fit:badPrior', 'data.priors(%d): a group pattern needs the streaming kind (a commonParams table).', k);
            end
            L = commonRowsOfPattern(P(k).group, spec);
            if isempty(L)
                error('lsqvce:fit:badPrior', 'data.priors(%d): the pattern ''%s'' matches no name in commonParams.', k, char(string(P(k).group)));
            end
            nm = char(string(P(k).group));
        elseif hasL
            if ~isnumeric(P(k).L) || size(P(k).L, 2) ~= nc
                error('lsqvce:fit:badPrior', ...
                      'data.priors(%d).L must be a p x %d numeric matrix (nc = %d).', k, nc, nc);
            end
            L = double(P(k).L);
        else
            error('lsqvce:fit:badPrior', 'data.priors(%d) needs group (a name pattern) or L.', k);
        end
        s = 1;
        if isfield(P(k), 'sigma') && ~isempty(P(k).sigma), s = double(P(k).sigma); end
        if ~isscalar(s) || ~isfinite(s) || s <= 0
            error('lsqvce:fit:badPrior', 'data.priors(%d).sigma must be a positive scalar.', k);
        end
        x0 = zeros(size(L, 1), 1);
        if isfield(P(k), 'x0') && ~isempty(P(k).x0)
            x0 = double(P(k).x0(:));
            if isscalar(x0), x0 = repmat(x0, size(L, 1), 1); end
            if numel(x0) ~= size(L, 1)
                error('lsqvce:fit:badPrior', 'data.priors(%d).x0 must have %d entries (or be a scalar).', k, size(L, 1));
            end
        end
        e = 'none';
        if isfield(P(k), 'estimate') && ~isempty(P(k).estimate)
            e = lower(char(string(P(k).estimate)));
            if ~any(strcmp(e, {'none', 'scale'}))
                error('lsqvce:fit:badPrior', 'data.priors(%d).estimate must be "none" or "scale" (got "%s").', k, e);
            end
        end
        if isfield(P(k), 'name') && ~isempty(P(k).name), nm = char(string(P(k).name)); end
        fr = 1;
        if isfield(P(k), 'from') && ~isempty(P(k).from)
            fr = double(P(k).from);
            if ~isscalar(fr) || ~isfinite(fr) || fr < 1 || fr ~= round(fr)
                error('lsqvce:fit:badPrior', 'data.priors(%d).from must be a pass number >= 1.', k);
            end
        end
        groups(end+1) = struct('L', L, 'x0', x0, 'sigma', s, 'name', nm, 'estimate', e, 'from', fr); %#ok<AGROW>
    end
    pri = stackPriors(groups);
end

function L = commonRowsOfPattern(pattern, spec)
%COMMONROWSOFPATTERN  L for a 'group' pattern on commonParams: one unit row
%   per column of A whose template row's name matches the regexp.
    cp = spec.commonPlan;  T = spec.commonParams;
    nameCol = 2;
    if isfield(spec, 'nameColumn') && ~isempty(spec.nameColumn), nameCol = spec.nameColumn; end
    if nameCol > size(T, 2)
        error('lsqvce:fit:badPrior', 'the name column %d does not exist in commonParams (%d columns).', nameCol, size(T, 2));
    end
    hit = false(size(T, 1), 1);
    for r = 1:size(T, 1)
        v = T{r, nameCol};
        if ischar(v) || (isstring(v) && isscalar(v))
            hit(r) = ~isempty(regexp(char(v), char(string(pattern)), 'once'));
        end
    end
    rowOf = zeros(cp.nc, 1);
    rowOf(1:cp.commonCount) = cp.rows.common(:);
    for i = 1:numel(cp.days)
        if cp.dayCount(i) == 0, continue, end
        rows = cp.rows.day;  if isfield(cp, 'dayRowsOf'), rows = cp.dayRowsOf{i}; end
        rowOf(cp.dayOffset(i) + (1:cp.dayCount(i))) = rows(1:cp.dayCount(i));
    end
    cols = find(rowOf > 0 & hit(max(rowOf, 1)));
    L = zeros(numel(cols), cp.nc);
    for k = 1:numel(cols), L(k, cols(k)) = 1; end
end
function pri = stackPriors(groups, pass)
%STACKPRIORS  The groups ACTIVE at this pass (from <= pass; default: all)
%   with their CURRENT sigmas as one CommonPrior (sigma 1).
    if nargin < 2, pass = Inf; end
    if isempty(groups), pri = [];  return, end
    act = groups([groups.from] <= pass);
    if isempty(act), pri = [];  return, end
    Ks = arrayfun(@(g) g.L / g.sigma, act, 'UniformOutput', false);
    xs = arrayfun(@(g) g.x0 / g.sigma, act, 'UniformOutput', false);
    pri = lsqvce.engine.CommonPrior(vertcat(Ks{:}), vertcat(xs{:}), 1);
end

function [groups, rec] = priorUpdate(groups, lastSol)
%PRIORUPDATE  One Helmert / Koch-Kusche step for every prior group with
%   estimate "scale" (R5 M7 7b): with the previous pass's x, Cxx, s0sq,
%     Omega_k = || (L_k x - x0_k) / sigma_k ||^2                 (the group's residual sum)
%     r_k     = p_k - tr( L_k N^-1 L_k' ) / sigma_k^2,   N^-1 = Cxx / s0sq   (its redundancy)
%     sigma_k^2 <- sigma_k^2 * Omega_k / r_k                    (the factor on the current sigma)
%   The fixed point Omega_k = r_k is the classic condition "the prior's
%   weighted residual sum equals its redundancy": lambda = s0sq / sigma_k^2 is
%   then the regularisation parameter the data support. A group whose
%   redundancy is not positive (it absorbs nothing) is left unchanged and
%   flagged in rec.
%   OUTPUT  rec 1xG struct (.name .sigma .omega .r .updated)
    rec = struct('name', {}, 'sigma', {}, 'omega', {}, 'r', {}, 'updated', {});
    Ninv = lastSol.Cxx / max(lastSol.s0sq, realmin);
    for k = 1:numel(groups)
        g = groups(k);
        v = (g.L * lastSol.x - g.x0) / g.sigma;
        om = v.' * v;
        r  = size(g.L, 1) - trace(g.L * Ninv * g.L.') / g.sigma^2;
        upd = strcmp(g.estimate, 'scale') && isfinite(r) && r > 0 && isfinite(om) && om > 0;
        if upd, groups(k).sigma = g.sigma * sqrt(om / r); end
        rec(end+1) = struct('name', g.name, 'sigma', groups(k).sigma, 'omega', om, 'r', r, 'updated', upd); %#ok<AGROW>
    end
end
function R = runBatchOnBlockPath(u, batchFcn, nBatch, opts, kind, pri, data)
%RUNBATCHONBLOCKPATH  Kinds B / C (linear) as ONE local-free streaming system
%   through LSQVCE.ENGINE.BATCHSYSTEM (R5 M2 WP-81 2a/2c'; priors M5 5b).
    if nargin < 6, pri = []; end
    if nargin < 7, data = struct(); end
    L = lsqvce.engine.batchSystem(u, batchFcn, nBatch, 'Name', kind);
    lsqvce.util.log.inventory('reset');
    prov = opts.toStruct();  prov.kind = kind;   % provenance = the user's resolved Options
    if strcmp(opts.Robust, "none")
        a = opts.engineArgs('adjustBlock');
        R = lsqvce.engine.adjustBlock({L}, 1:nBatch, a.nv{:}, 'Provenance', prov, 'Prior', pri);
    else
        a = opts.engineArgs('adjustRobust');
        R = lsqvce.engine.adjustRobust({L}, 1:nBatch, a.nv{:}, 'Provenance', prov, 'Prior', pri, ...
                                       'OnPassContext', struct('data', data, 'processingOptions', {{}}));
    end
end

function [estMap, afMap] = streamPolicy(systems, specs, streams, afDefault)
%STREAMPOLICY  Per-stream 'estimate' / 'arcFactors' from the covariance specs.
%   R5 M3 3d (WP-83, R5-D4): a stream spec may carry
%     .estimate    "none" | "scale" | "law"   default "law"
%     .arcFactors  logical                    default afDefault (the Model's:
%                                            true for "ellmer", false for "vce"
%                                            and "pooled" -- R5 M4 4b)
%   Returned as containers.Map keyed 'satKey|obsType' for the joint driver.
    estMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    afMap  = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for k = 1:numel(streams)
        % R6 WP-99: seed the estimate map too, so it is TOTAL over the streams.
        % The default "law" is then stated ONCE, here, instead of being restated
        % by every reader of the map (R1: one authority per fact).
        estMap(char(string(streams{k}))) = 'law';
        afMap(char(string(streams{k}))) = logical(afDefault);
        if contains(char(string(streams{k})), '|prior_')     % R5-D8 round 5: a local prior group is a
            afMap(char(string(streams{k}))) = false;          % stream; its estimate policy comes from
        end                                                   % its cov spec (the group's "estimate")
    end
    for k = 1:numel(systems)
        L = systems{k};
        if ~isfield(L, 'cov') || ~isstruct(L.cov), continue, end
        satKey = matlab.lang.makeValidName(char(string(specs(k).sysID)));
        ots = fieldnames(L.cov);
        for j = 1:numel(ots)
            cm = L.cov.(ots{j});
            key = sprintf('%s|%s', satKey, ots{j});
            % R7-07: a 'var' law carries its own cross-axis structure, which the
            % joint sweep cannot use -- it would treat the d-interleaved block-
            % Toeplitz Cfcn as a per-axis temporal law over epochs (wrong Sigma,
            % silently) and apply an estimated level as a std to a covariance.
            % Refused by name until the sweep supports it.
            if strcmpi(char(string(lsqvce.util.getfielddef(cm, 'type', ''))), 'var')
                error('lsqvce:fit:varUnderNoise', ...
                     ['covariance.%s of system %s is a ''var'' law: Estimate = "noise" ' ...
                      'cannot estimate it (the joint sweep has no vector-AR path). Run it ' ...
                      'with Estimate = "none", or declare the axes as a scalar temporal ' ...
                      'law plus the per-epoch Cep.'], ots{j}, satKey);
            end
            % a law with a NON-SCALAR sigma (per-row weights) is a fixed
            % structure: it can carry a LEVEL, not a spectral shape -- its
            % default policy is "scale"; asking for "law" is refused (R5 M4 4b)
            sig = lsqvce.util.getfielddef(cm, 'sigma', 1);
            fixedStructure = isnumeric(sig) && numel(sig) > 1;
            if fixedStructure, estMap(key) = 'scale'; end
            if isfield(cm, 'estimate') && ~isempty(cm.estimate)
                e = lower(char(string(cm.estimate)));
                if ~any(strcmp(e, {'none', 'scale', 'law'}))
                    error('lsqvce:fit:badEstimate', ...
                          'covariance.%s.estimate must be "none", "scale" or "law" (got "%s").', ots{j}, e);
                end
                if fixedStructure && strcmp(e, 'law')
                    error('lsqvce:fit:lawOnFixedStructure', ...
                         ['covariance.%s declares a per-row sigma (%d values): a fixed ' ...
                          'structure can take a LEVEL ("scale") but not a spectral shape ' ...
                          '("law"). Declare one scalar sigma if the temporal law is to be estimated.'], ...
                          ots{j}, numel(sig));
                end
                estMap(key) = e;
            end
            if isfield(cm, 'arcFactors') && ~isempty(cm.arcFactors)
                afMap(key) = logical(cm.arcFactors);
            end
        end
    end
end

function [names, models] = streamModels(systems, specs, out, pool, estMap)
%STREAMMODELS  Per-stream noise summary of a joint run (R5 M4 4b).
%   names{k}  'satKey|obsType';  models{k} struct('type','psd','sigma',level,
%   'C', the stream's law, 'sigma2', its arc factors) -- the LEVEL is the
%   declared sigma times sqrt(mean arc factor) of the stream.
%
%   R6 WP-99: .C is the law of the stream's POOL GROUP -- and it is EMPTY for a
%   stream whose law was not estimated ("none" or "scale"), because there is no
%   such law: what that stream used is the law YOU declared, which the toolbox
%   did not touch. Until this round the field fell back to the global out.C, so
%   a "none" stream reported a law fitted from somebody else's arcs; measured
%   on a four-stream run with ONE "law" stream, all four .C came back identical
%   (lag-1 1.313 for every one of them). A plausible wrong value returned in
%   silence is the one outcome the doctrine forbids (R4).
    names = {};  models = {};
    sMean = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:numel(out.arcKeys)
        tok = regexp(out.arcKeys{i}, '^s(.+?)_(.+)_d[^_]+_a', 'tokens', 'once');   % R5-D8 6: obsType may carry underscores
        if isempty(tok), continue, end
        sk = sprintf('%s|%s', matlab.lang.makeValidName(tok{1}), tok{2});
        if ~isKey(sMean, sk), sMean(sk) = []; end
        sMean(sk) = [sMean(sk); out.sigma2(i)];
    end
    for k = 1:numel(systems)
        L = systems{k};
        if ~isfield(L, 'cov') || ~isstruct(L.cov), continue, end
        satKey = matlab.lang.makeValidName(char(string(specs(k).sysID)));
        ots = fieldnames(L.cov);
        for j = 1:numel(ots)
            sk = sprintf('%s|%s', satKey, ots{j});
            s0 = lsqvce.util.getfielddef(L.cov.(ots{j}), 'sigma', 1);
            if ~isscalar(s0), s0 = 1; end
            f = 1;  s2 = [];
            if isKey(sMean, sk) && ~isempty(sMean(sk))
                s2 = sMean(sk);  f = sqrt(max(mean(s2), eps));
            end
            switch lower(pool)
                case 'obstype', gk = ots{j};
                case 'stream',  gk = sk;
                otherwise,      gk = 'all';
            end
            % R6 WP-99: only a stream whose law was ESTIMATED has one to report.
            % estMap is total over the streams (streamPolicy seeds it), so this
            % reads the policy rather than restating its default.
            C = [];
            if isKey(estMap, sk) && strcmp(estMap(sk), 'law') && isKey(out.Cgroups, gk)
                C = out.Cgroups(gk);
            end
            names{end+1}  = sk;                                                %#ok<AGROW>
            models{end+1} = struct('type', 'psd', 'sigma', s0 * f, 'C', C, 'sigma2', s2); %#ok<AGROW>
        end
    end
end

function opc = onPassContext(data, systems)
%ONPASSCONTEXT  What lsqvce.fit adds to every OnPass ctx (R5 M5 5c'): the
%   declaration as given and every system's processingOptions, in the order
%   of data.systems (a cell, one struct per system).
    po = cell(1, numel(systems));
    for k = 1:numel(systems)
        if isfield(systems{k}, 'procOpt'), po{k} = systems{k}.procOpt; else, po{k} = struct(); end
    end
    opc = struct('data', data, 'processingOptions', {po});
end

function pri = priorFor(priState, model)
%PRIORFOR  The stacked prior for THIS pass: the groups updated from the
%   previous pass's solution when they carry estimate "scale" (R5 M7 7b).
    groups = priState('groups');
    if isempty(groups), pri = [];  return, end
    pass = 1;
    if isfield(model, 'iter'), pass = model.iter; end
    if isfield(model, 'lastSol') && ~isempty(model.lastSol) && any(strcmp({groups.estimate}, 'scale'))
        active = [groups.from] < pass;                 % only groups that were IN the last solve
        [groups(active), rec] = priorUpdate(groups(active), model.lastSol);
        priState('groups') = groups;
        recs = priState('rec');  recs{end+1} = rec;  priState('rec') = recs;
    end
    pri = stackPriors(groups, pass);                   % R5-D8: activation by pass (from)
    priState('pri') = pri;
end

function uo = withPass(uo, model)
%WITHPASS  The pass number rides on the loader's UserOpt (R5-D8 round 2: the
%   local priors' 'from' is honoured by lsqvce.engine.makeLoader>priorArc),
%   and so do the robust factors of the previous pass (R6 WP-105:
%   model.rowWeights -> UserOpt.epochWeights, folded into every arc's
%   per-epoch error information BEFORE the temporal law by
%   lsqvce.engine.makeLoader>batchOf -- the Koch convention).
    uo.pass = 1;
    if isfield(model, 'iter'), uo.pass = model.iter; end
    if isfield(model, 'rowWeights') && ~isempty(model.rowWeights) && model.rowWeights.Count > 0
        uo.epochWeights = model.rowWeights;     % per epoch, BEFORE the temporal law (Koch; makeLoader>batchOf)
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
