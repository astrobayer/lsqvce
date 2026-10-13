function red = redundancy(acc, sol, opt)
%LSQVCE.ENGINE.REDUNDANCY  Per-group and per-observation redundancy, exactly.
%
%   red = LSQVCE.ENGINE.REDUNDANCY(acc, sol, 'Name', value, ...)
%
%   A second pass over the data that computes, for every observation,
%
%       h_i = a_i' inv(N) a_i        leverage
%       r_i = 1 - h_i                partial redundancy
%
%   and aggregates r_i into groups. This is the machinery behind LSQVCE.ENGINE.VCE and
%   LSQVCE.ENGINE.OUTLIERS -- both are just different uses of the same numbers, so they
%   share one data pass.
%
%   INPUT
%     acc      1x1 lsqvce.engine.BlockAccumulator   after solve
%     sol      1x1 struct                            its solution (x, Cxx, s0sq)
%     Name, Value                                    the OPTIONS below
%
%   OPTIONS
%   (WP-39 stage (b): canonical UpperCamelCase names below; legacy
%   casings still match -- MATLAB NV semantics -- and enum VALUES are
%   case-sensitive lowercase.)
%     'GroupBy'    'arc' (default) | 'obsType' | 'system' | 'day' | function handle
%                  @(sysID, obsType, arc, day) -> group key (numeric or string)
%     'MinRedundancy'  warn if any group has r_g below this (default 50)
%     'Collect'    also return the per-observation h and r (default true)
%     'Verbose'    default true
%     'Parallel'   'serial' | 'parfor' | 'spmd'. Default: the accumulator's
%                  'Parallel' option. The pass is parallel OVER ARCS: the
%                  block cofactor (one commons factor, the local blocks) is
%                  built once and broadcast,
%                  then each arc-task runs independently and the partial sums
%                  are reduced in task order (proven exact in
%                  validation/v4verify.py). This is THE pass to parallelise:
%                  at global SH d/o 90 it is ~4 h serial per Helmert
%                  iteration; with 16 workers plausibly 30-40 min. 'spmd'
%                  reduces on the workers (multi-node) and forces
%                  'Collect', false.
%     'TraceMethod' 'exact' (default) | 'hutchinson'. 'exact' computes the
%                  redundancies from the per-arc leverage below (exact,
%                  O(nc^2)/arc). 'hutchinson' ALSO runs the Monte-Carlo trace
%                  estimator of Ellmer (2018, 6.5.6) via LSQVCE.ENGINE.TRACEESTIMATE and
%                  returns it in red.traceEstimate. EXPERIMENTAL and OPT-IN: the
%                  exact per-group r_g remain authoritative (the estimator gives
%                  a per-group scalar only, not the per-observation r_i the
%                  snooping and the sum(r)=m-u identity need). Intended as an
%                  accelerator/cross-check at very large nc; read the accuracy
%                  caveat in LSQVCE.ENGINE.BLOCKCOFACTOR.
%     'Probes'     1x1  number of Hutchinson probes Z (default 100).
%     'TraceSeed'  1x1  RNG seed for reproducible probes (default []).
%
%   OUTPUT
%     red.group(k).key       group identifier
%     red.group(k).m         observations in the group
%     red.group(k).r         redundancy r_g = sum(r_i)
%     red.group(k).vtv       weighted residual square sum
%     red.group(k).tcc       common trace tr(inv(N) N_g^cc): the group's share
%                            of the common solution (sums to nc over the data,
%                            less what a common prior absorbs, tr(N^-1 K'K))
%     red.group(k).s0sq      vtv / r    -- the variance component of this group
%     red.group(k).sigma0    sqrt(s0sq)
%     red.obs(j).h, .r       per-observation, if Collect
%     red.rTotal             sum over all groups
%     red.priorRedundancy    r_K of a common prior's nk rows (0 without one)
%     red.check              |sum(r_g) + r_K - (m + nk - u)|, must be ~0
%     red.traceEstimate      present only with 'TraceMethod','hutchinson': a
%                            struct (.value/.expected/.relErr/.se/.relSE/.probes/
%                            .note) with the Hutchinson estimate of tr(N^-1 N) on
%                            the reduced common normal (a known = nc machinery
%                            check on THIS adjustment's conditioning).
%
%   THE GLOBAL CHECK
%     sum(r_i) over ALL observations must equal m - u exactly. This is the
%     strongest available test that the block-arrow inverse algebra is right,
%     and it is verified on every call. A nonzero red.check means something is
%     structurally wrong -- do not trust the variance components.
%
%   COST -- KNOW THIS BEFORE YOU RUN IT
%     One pass costs O(nc^3/3) + O(m*nc^2/2), i.e. about the SAME as one
%     accumulation pass. Measured, at 25 GFLOP/s:
%
%       regional mascon    nc=300    m=2e6      ~4 s
%       regional (large)   nc=2000   m=2e6      ~3 min
%       global SH d/o 60   nc=3721   m=1e7      ~45 min
%       global SH d/o 90   nc=8281   m=1e7      ~4 h
%       global SH d/o 120  nc=14641  m=1e7      ~12 h
%
%     So VCE roughly DOUBLES the cost per Helmert iteration, and you need 2-4
%     iterations. That is the honest price of exact per-arc redundancy. It is
%     not reducible by a stochastic estimator -- that was tried and gives
%     20-50% error on realistically conditioned normals (see LSQVCE.ENGINE.BLOCKCOFACTOR).
%
%   CHOOSING THE GROUP
%     'arc' is the default and the physically meaningful unit: an arc is a
%     contiguous data span whose noise level genuinely can differ from its
%     neighbours. Do NOT group by batch -- a batch is an artefact of BatchSize,
%     chosen for computational reasons, and grouping by it makes your variance
%     components depend on a tuning parameter.
%
%     If arcs are short, r_g gets small and the variance component gets noisy:
%     its relative precision is sqrt(2/r_g). At r_g = 50 that is 20%. Pool to
%     'obsType' or 'day' in that case; the warning fires automatically.
%
%   EXAMPLE
%     sol = acc.solve();
%     red = lsqvce.engine.redundancy(acc, sol, 'GroupBy', 'arc');
%     disp(red.check)                             % must be ~1e-9
%
%     for k = 1:numel(red.group)
%         g = red.group(k);
%         lsqvce.util.log.msg('%-16s r=%8.1f  sigma=%.4f\n', string(g.key), g.r, g.sigma0);
%     end
%
%   ERRORS
%     lsqvce:makeKeyFcn:badGroupBy  the grouping key is not a known name or
%                   handle (WP-38: the one shared key resolver)
%     lsqvce:redundancy:noStore     nothing was retained; re-run with 'store','all'
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:redundancy:hutchFailed              [warning] Hutchinson trace check failed (X); exact r_g are unaffected.valueexpe...
%     lsqvce:redundancy:hutchinsonExperimental   [warning] 'TraceMethod','hutchinson' is EXPERIMENTAL
%     lsqvce:redundancy:identityViolated         [warning] sum(r_i) + r_prior = X but m + nk - u = X (relative error X)
%     lsqvce:redundancy:smallGroup               [warning] X group(s) have redundancy below X. The variance component of a group...
%     lsqvce:redundancy:spmdNoCollect            [warning] 'spmd' reduces on the workers, so per-observation h/r/v are not shipp...
%
%   See also LSQVCE.ENGINE.BLOCKCOFACTOR, LSQVCE.ENGINE.VCE, LSQVCE.ENGINE.OUTLIERS, LSQVCE.ENGINE.RESIDUALS
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-130: the identity's u counts a
%   rank-deficient local block by its rank (BlockAccumulator.localRankOf). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-125: group(k).tcc, the common trace
%   tr(inv(N) N_g^cc) of the group. With assistance from Claude (Anthropic),
%   model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-114 (R7-15): the identity counts a common
%   prior (red.priorRedundancy). With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 3 (R7-24): the
%   leverage comes from lsqvce.engine.BlockCofactor -- one commons factor
%   broadcast, not one per system (leverage / leverageBlocks retired). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) (F-07): typed
%   `arguments` signature, canonical UpperCamelCase option names (legacy
%   casings keep matching; enum VALUES case-sensitive lowercase). With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-38 (fixes F-08): the key
%   presets and the batch-and-whiten skeleton live in the ONE private
%   makeKeyFcn / arcBatchWhiten (byte-identical numerics; the GroupBy
%   error id is now lsqvce:makeKeyFcn:badGroupBy). With assistance from
%   Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-37: the arc-task enumeration
%   caches the meta calls and PREALLOCATES the task array (the old
%   tasks(end+1) growth copied O(nTask^2)). With assistance from Claude
%   (Anthropic), model Claude Fable 5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
arguments
    acc
    sol
    opt.GroupBy = 'arc'
    opt.MinRedundancy (1,1) double {mustBePositive} = 50
    opt.Collect (1,1) logical = true
    opt.Verbose (1,1) logical = true
    opt.Parallel (1,1) string {mustBeMember(opt.Parallel, ["","serial","parfor","spmd"])} = ""
    opt.TraceMethod (1,1) string {mustBeMember(opt.TraceMethod, ["exact","hutchinson"])} = "exact"
    opt.Probes (1,1) double {mustBeInteger, mustBePositive} = 100
    opt.TraceSeed double {mustBeScalarOrEmpty} = []
end

    narginchk(2,inf);

    % WP-39 stage (b): typed `arguments` signature (canonical UpperCamel;
    % NAME matching case-insensitive per MATLAB semantics; enum VALUES
    % case-sensitive lowercase). Lowercase mapping keeps the body untouched.
    o = struct('groupBy', opt.GroupBy, 'minRedundancy', opt.MinRedundancy, ...
               'collect', opt.Collect, 'verbose', opt.Verbose, ...
               'parallel', char(opt.Parallel), ...
               'traceMethod', char(opt.TraceMethod), 'probes', opt.Probes, ...
               'traceSeed', opt.TraceSeed);

    useHutch = strcmpi(o.traceMethod, 'hutchinson');
    if useHutch
        % The exact per-arc path yields per-OBSERVATION r_i (needed for snooping,
        % the sum(r_i)=m-u identity, and the (w,r) diagnostics). The Hutchinson
        % estimator yields only a per-GROUP scalar trace, so it cannot replace the
        % per-arc path wholesale. In this release the estimator is wired as an
        % OPT-IN GLOBAL cross-check (a stochastic estimate of the total common
        % redundancy) computed alongside the exact per-group r_g, which remain
        % authoritative. The estimate is returned in red.traceEstimate for the
        % caller to compare. Set nc large and compare cost/accuracy before relying
        % on it. Full per-group estimator substitution is future work (it needs a
        % group-level common-normal assembly; see HISTORY_v1.md).
        lsqvce.util.log.warn('lsqvce:redundancy:hutchinsonExperimental', ...
            ['''TraceMethod'',''hutchinson'' is EXPERIMENTAL: the exact per-arc ' ...
             'redundancies are still computed and remain authoritative; the ' ...
             'Hutchinson estimator (Ellmer 2018 6.5.6) is returned in ' ...
             'red.traceEstimate as a stochastic estimate of the TOTAL common ' ...
             'redundancy for comparison. See LSQVCE.ENGINE.TRACEESTIMATE and the ' ...
             'accuracy note in LSQVCE.ENGINE.BLOCKCOFACTOR.']);
    end

    % the execution mode: explicit option > accumulator option > serial
    mode = lower(o.parallel);
    if isempty(mode)
        mode = lower(lsqvce.util.getfielddef(acc.opt, 'Parallel', 'serial'));
    end
    thr = lsqvce.util.getfielddef(acc.opt, 'ParallelThreshold', 8);

    if strcmp(mode, 'spmd') && o.collect
        lsqvce.util.log.warn('lsqvce:redundancy:spmdNoCollect', ...
            ['''spmd'' reduces on the workers, so per-observation h/r/v are ' ...
             'not shipped back. Forcing ''Collect'', false. Use ''parfor'' ' ...
             'if you need the per-observation series.']);
        o.collect = false;
    end

    keyFcn = makeKeyFcn(o.groupBy);

    x   = sol.x(:);
    ids = acc.systemIDs();
    if isempty(ids)
        error('lsqvce:redundancy:noStore', ...
              'Construct the accumulator with ''KeepLocal'', true.');
    end

    % ---- phase 1: the broadcast state (serial, cheap) ---------------------
    % The cofactor is the ONLY O(nc^3)/O(nl^3) work: ONE commons factor and
    % each system's local blocks (lsqvce.engine.BlockCofactor, WP-86 part 3 --
    % before, every system carried its own copy of the nc x nc factor and the
    % nc x nl products, nSys x 0.7 GB at a GRAVFIRE month). Broadcast to every
    % task; the accumulator itself is NEVER captured by a task: only Q, these
    % plain arrays and the loader.
    Q       = lsqvce.engine.BlockCofactor.fromAccumulator(acc);
    nSys    = numel(ids);
    stores  = cell(nSys, 1);
    ys      = cell(nSys, 1);
    nParTot = numel(x);
    for k = 1:nSys
        st        = acc.getStore(ids(k));
        ys{k}     = lsqvce.engine.BlockAccumulator.localSolveOf(st, x);   % WP-128
        stores{k} = struct('sysID', st.sysID, 'dataFcn', st.dataFcn, ...
                           'uo', st.uo, 'dayList', st.dayList);
        nParTot   = nParTot + lsqvce.engine.BlockAccumulator.localRankOf(st);   % WP-130: the rank
    end

    % ---- phase 2: enumerate the independent arc-tasks (meta calls only) ---
    % WP-37 (hot-loop toll): the meta calls are made ONCE and cached, the
    % task count is summed, and the struct array is PREALLOCATED -- the
    % previous tasks(end+1) growth re-copied the whole array per arc
    % (O(nTask^2) copies on big campaigns).
    metas = cell(nSys, 1);
    nTask = 0;
    for k = 1:nSys
        st = stores{k};
        days = st.dayList(:).';
        mk = cell(numel(days), 1);
        for di = 1:numel(days)
            mk{di} = st.dataFcn('meta', days(di), [], st.uo);
            for ot = 1:numel(mk{di}.obsTypes)
                nTask = nTask + numel(mk{di}.arcs{ot});
            end
        end
        metas{k} = mk;
    end
    tasks = repmat(struct('k',0,'day',0,'obsType','','arc',0), nTask, 1);
    t = 0;
    for k = 1:nSys
        st = stores{k};
        days = st.dayList(:).';
        for di = 1:numel(days)
            meta = metas{k}{di};
            for ot = 1:numel(meta.obsTypes)
                for a = meta.arcs{ot}(:).'
                    t = t + 1;
                    tasks(t) = struct('k', k, 'day', days(di), ...
                        'obsType', meta.obsTypes{ot}, 'arc', a);
                end
            end
        end
    end

    % ---- phase 3: the data pass, one task per arc --------------------------
    % Each task loads its arc, whitens it, and returns the per-arc sums
    % (m, sum r, v'v, sum h) plus, if Collect, the per-observation series.
    % validation/v4verify.py proves partial sums over any partition of these
    % tasks reconstruct the serial result -- the schedule is free.
    taskFcn = @(j) arcTask(tasks(j), stores{tasks(j).k}, Q, tasks(j).k, ...
                           ys{tasks(j).k}, x, keyFcn, o.collect);

    switch mode
        case 'spmd'
            [P, info] = lsqvce.engine.mapreduce(nTask, taskFcn, 'Mode', 'spmd', ...
                            'Threshold', thr, 'Verbose', o.verbose, ...
                            'Reduce',  @foldTask, 'Init', emptyPartial(), ...
                            'Combine', @combinePartials);
            if iscell(P)                     % fell back to a collecting mode
                Q = emptyPartial();
                for j = 1:numel(P), Q = foldTask(Q, P{j}); end
                P = Q;
            end
        otherwise
            [raw, info] = lsqvce.engine.mapreduce(nTask, taskFcn, 'Mode', mode, ...
                              'Threshold', thr, 'Verbose', o.verbose);
            P = emptyPartial();
            for j = 1:nTask                  % ordered: identical to serial
                P = foldTask(P, raw{j});
            end
    end

    if lsqvce.util.verbosity.progress(o.verbose) && ~strcmp(info.mode, 'serial')
        lsqvce.util.log.msg('%s%d arc-tasks on %d workers (%s)\n', lsqvce.util.verbosity.banner('Redundancy'), ...
                nTask, info.nWorkers, info.mode);
    end

    % ---- assemble ---------------------------------------------------------
    mTot    = P.mTot;
    obsList = struct([]);
    if o.collect && ~isempty(P.obs)
        obsList = [P.obs{:}];
    end

    G = struct([]);
    for k = 1:numel(P.keys)
        e = struct('key', P.keyVals{k}, 'm', P.m(k), 'r', P.r(k), ...
                   'vtv', P.vtv(k), 'tcc', P.tcc(k));
        e.s0sq   = e.vtv / max(e.r, realmin);
        e.sigma0 = sqrt(e.s0sq);
        if isempty(G), G = e; else, G(end+1) = e; end            %#ok<AGROW>
    end

    red = struct('group', G, 'obs', obsList, ...
                 'rTotal', sum([G.r]), 'm', mTot, 'u', nParTot);

    % ---- optional Hutchinson trace estimator (Ellmer 2018 6.5.6) ----------
    % Opt-in, EXPERIMENTAL. Demonstrates the estimator on the reduced common
    % normal against a KNOWN answer: tr(Ncc_red^-1 * Ncc_red) = nc exactly, so
    % the estimate must return ~nc. This validates the estimator machinery end
    % to end on the actual normal matrix of THIS adjustment (probe count,
    % conditioning) without yet replacing the exact per-group r_g. A future
    % release will estimate the per-group tr(N^-1 N_g^cc) directly (needs a
    % group-level common-normal assembly). The exact r_g above are unchanged.
    if useHutch
        Ncc_red = acc.reducedNormalMatrix();      % TOTAL, summed over all systems
        try
            U = chol((Ncc_red + Ncc_red.')/2);        % Ncc_red = U'*U
            nc = size(Ncc_red, 1);
            % tr(N^-1 N) = tr(U^-1 (U'U) U^-T) = tr(I) = nc. With N = U'U the
            % identity operator is U^-T N U^-1 (NOT U^-1 N U^-T, which is a
            % similarity that is not the identity and biases the trace): apply
            % the upper solve on the LEFT and the transpose solve on the RIGHT.
            Xfun = @(V) U.' \ (Ncc_red * (U \ V));     % U^-T N U^-1 = I, matrix-free
            [tEst, tInfo] = lsqvce.engine.traceEstimate(Xfun, 'Size', nc, ...
                                'NumProbes', o.probes, 'Seed', o.traceSeed);
            red.traceEstimate = struct('value', tEst, 'expected', nc, ...
                'relErr', abs(tEst - nc)/max(nc,1), 'se', tInfo.se, ...
                'relSE', tInfo.relse, 'probes', o.probes, ...
                'note', ['machinery check on Ncc_red: tr(N^-1 N)=nc; not a ' ...
                         'per-group redundancy (see help)']);
            if lsqvce.util.verbosity.progress(o.verbose)
                lsqvce.util.log.msg(['%sHutchinson check: tr(N^-1 N) est %.2f ' ...
                    'vs nc = %d  (relErr %.2e, relSE %.2e, Z = %d)\n'], lsqvce.util.verbosity.banner('Redundancy'), ...
                    tEst, nc, red.traceEstimate.relErr, tInfo.relse, o.probes);
            end
        catch ME
            lsqvce.util.log.warn('lsqvce:redundancy:hutchFailed', ...
                'Hutchinson trace check failed (%s); exact r_g are unaffected.', ...
                ME.identifier);
            red.traceEstimate = struct('value', NaN, 'expected', size(Ncc_red,1), ...
                'relErr', NaN, 'se', NaN, 'relSE', NaN, 'probes', o.probes, ...
                'note', 'estimator failed; see warning');
        end
    end

    % ---- the global identity: sum(r_i) + r_K = m + nk - u (R7 WP-114, R7-15:
    % a common prior's rows and their redundancy r_K, on the solve's factor) --
    nk = 0;  red.priorRedundancy = 0;
    if ~isempty(acc.prior), nk = acc.prior.nk;  red.priorRedundancy = acc.prior.redundancy(acc.commonFactor().'); end
    rExp      = mTot + nk - nParTot;
    red.check = abs(red.rTotal + red.priorRedundancy - rExp);
    relChk    = red.check / max(rExp, 1);

    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg('%sgroups = %d   sum(r) = %.4f   r_prior = %.4f   m+nk-u = %d   check = %.2e\n', ...
                lsqvce.util.verbosity.banner('Redundancy'), numel(G), red.rTotal, red.priorRedundancy, rExp, red.check); end

    if relChk > 1e-6
        lsqvce.util.log.warn('lsqvce:redundancy:identityViolated', ...
            ['sum(r_i) + r_prior = %.6f but m + nk - u = %d (relative error %.2e). The redundancy identity ' ...
             'MUST hold exactly. Do not trust the variance components. Likely causes: the loader returned ' ...
             'different data on this pass, or the stored Schur blocks are stale.'], ...
             red.rTotal + red.priorRedundancy, rExp, relChk);
    end

    % ---- small-group warning ---------------------------------------------
    small = find([G.r] < o.minRedundancy);
    if ~isempty(small) && o.verbose
        lsqvce.util.log.warn('lsqvce:redundancy:smallGroup', ...
            ['%d group(s) have redundancy below %g. The variance component of ' ...
             'a group has relative precision sqrt(2/r_g), so at r_g = %.0f that ' ...
             'is %.0f%%. Consider pooling: ''GroupBy'', ''obsType'' or ''day''.'], ...
             numel(small), o.minRedundancy, min([G(small).r]), ...
             100*sqrt(2/max(min([G(small).r]),1)));
    end
end

%==========================================================================
function t = arcTask(tk, st, Q, s, y, x, keyFcn, collect)
%ARCTASK  One independent unit of the redundancy pass: load, whiten, sum.
    uo = st.uo;
    uo.obsType = tk.obsType;
    uo.day     = tk.day;
    uo.arc     = tk.arc;

    arc = st.dataFcn(tk.day, tk.obsType, tk.arc, uo);
    if isempty(arc) || ~isfield(arc,'l') || isempty(arc.l)
        t = [];                                    % data gap
        return
    end

    [h, v, hc] = arcLeverageAndResidual(arc, x, y, uo, Q, s);
    r = 1 - h;

    key = keyFcn(st.sysID, tk.obsType, tk.arc, tk.day);
    t = struct('ckey', char(string(key)), 'key', key, 'm', numel(v), ...
               'r', sum(r), 'vtv', v.'*v, 'h', sum(h), 'tcc', sum(hc), 'rec', []);
    if collect
        t.rec = struct('sysID', st.sysID, 'obsType', tk.obsType, ...
                       'day', tk.day, 'arc', tk.arc, 'h', h, 'r', r, ...
                       'v', v, 'key', key);
    end
end

%==========================================================================
function P = emptyPartial()
%EMPTYPARTIAL  The identity element for the group-redundancy fold: an empty
%   partial-accumulation struct that foldTask/reducePartials build up.
%   OUTPUT P  1x1 struct with fields
%            .keys    1x0 cell    group keys seen so far
%            .keyVals 1x0 cell    aligned group key values
%            .m       0x1 double  per-group observation counts
%            .r       0x1 double  per-group summed redundancy
%            .vtv     0x1 double  per-group weighted residual sum of squares
%            .tcc     0x1 double  per-group common trace tr(inv(N) N_g^cc)
%            .hTot    1x1 double  total leverage accumulated (0)
%            .mTot    1x1 double  total observations accumulated (0)
%            .obs     1x0 cell    per-observation records (when stored)
    P = struct('keys', {{}}, 'keyVals', {{}}, 'm', [], 'r', [], 'vtv', [], 'tcc', [], ...
               'hTot', 0, 'mTot', 0, 'obs', {{}});
end

function P = foldTask(P, t)
%FOLDTASK  Fold one task result into a partial. Applied in TASK ORDER for
%   serial and parfor, and per worker (ascending) for spmd.
    if isempty(t), return, end
    [tf, k] = ismember(t.ckey, P.keys);
    if ~tf
        P.keys{end+1}    = t.ckey;
        P.keyVals{end+1} = t.key;
        P.m(end+1) = 0;  P.r(end+1) = 0;  P.vtv(end+1) = 0;  P.tcc(end+1) = 0;
        k = numel(P.keys);
    end
    P.m(k)   = P.m(k)   + t.m;
    P.r(k)   = P.r(k)   + t.r;
    P.vtv(k) = P.vtv(k) + t.vtv;
    P.tcc(k) = P.tcc(k) + t.tcc;
    P.hTot   = P.hTot + t.h;
    P.mTot   = P.mTot + t.m;
    if ~isempty(t.rec)
        P.obs{end+1} = t.rec;
    end
end

function P = combinePartials(P, Q)
%COMBINEPARTIALS  Associative merge of two partials (the spmd combiner).
    for k = 1:numel(Q.keys)
        [tf, j] = ismember(Q.keys{k}, P.keys);
        if ~tf
            P.keys{end+1}    = Q.keys{k};
            P.keyVals{end+1} = Q.keyVals{k};
            P.m(end+1) = 0;  P.r(end+1) = 0;  P.vtv(end+1) = 0;  P.tcc(end+1) = 0;
            j = numel(P.keys);
        end
        P.m(j)   = P.m(j)   + Q.m(k);
        P.r(j)   = P.r(j)   + Q.r(k);
        P.vtv(j) = P.vtv(j) + Q.vtv(k);
        P.tcc(j) = P.tcc(j) + Q.tcc(k);
    end
    P.hTot = P.hTot + Q.hTot;
    P.mTot = P.mTot + Q.mTot;
    P.obs  = [P.obs, Q.obs];
end

%==========================================================================

%==========================================================================
function [h, v, hc] = arcLeverageAndResidual(arc, x, y, uo, Q, s)
%ARCLEVERAGEANDRESIDUAL  Whiten one arc, then get h and v, batched.
%   WP-38 (fixes F-08): the batch/overlap/whitening skeleton lives in the
%   ONE private arcBatchWhiten; this callback keeps the numeric payload
%   and the exact accumulation order (byte-identical outputs).

    m = size(arc.A, 1);
    h = zeros(m, 1);
    v = zeros(m, 1);
    hc = zeros(m, 1);
    arcBatchWhiten(arc, uo, @perBatch);

    function perBatch(rows, keep, Ab, Bb, lb)
        % the block cofactor's leverage of system s's rows; h in [0,1] up to
        % rounding (the sum(r) = m - u identity below is the integrity check)
        [hb, hcb] = Q.leverage(s, Ab, Bb);
        hb = min(max(hb, 0), 1);
        vb = lb - Ab*x - Bb*y;
        h(rows) = hb(keep);
        hc(rows) = hcb(keep);
        v(rows) = vb(keep);
    end
end
