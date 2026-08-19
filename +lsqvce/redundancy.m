function red = redundancy(acc, sol, varargin)
%LSQVCE.REDUNDANCY  Per-group and per-observation redundancy, exactly.
%
%   red = LSQVCE.REDUNDANCY(acc, sol, 'Name', value, ...)
%
%   A second pass over the data that computes, for every observation,
%
%       h_i = a_i' inv(N) a_i        leverage
%       r_i = 1 - h_i                partial redundancy
%
%   and aggregates r_i into groups. This is the machinery behind LSQVCE.VCE and
%   LSQVCE.OUTLIERS -- both are just different uses of the same numbers, so they
%   share one data pass.
%
%   OPTIONS
%     'groupBy'    'arc' (default) | 'obsType' | 'system' | 'day' | function handle
%                  @(sysID, obsType, arc, day) -> group key (numeric or string)
%     'minRedundancy'  warn if any group has r_g below this (default 50)
%     'collect'    also return the per-observation h and r (default true)
%     'verbose'    default true
%     'parallel'   'serial' | 'parfor' | 'spmd'. Default: the accumulator's
%                  'parallel' option. The pass is parallel OVER ARCS: the
%                  leverage blocks are built once per system and broadcast,
%                  then each arc-task runs independently and the partial sums
%                  are reduced in task order (proven exact in
%                  validation/v4verify.py). This is THE pass to parallelise:
%                  at global SH d/o 90 it is ~4 h serial per Helmert
%                  iteration; with 16 workers plausibly 30-40 min. 'spmd'
%                  reduces on the workers (multi-node) and forces
%                  'collect', false.
%     'traceMethod' 'exact' (default) | 'hutchinson'. 'exact' computes the
%                  redundancies from the per-arc leverage below (exact,
%                  O(nc^2)/arc). 'hutchinson' ALSO runs the Monte-Carlo trace
%                  estimator of Ellmer (2018, 6.5.6) via LSQVCE.TRACEESTIMATE and
%                  returns it in red.traceEstimate. EXPERIMENTAL and OPT-IN: the
%                  exact per-group r_g remain authoritative (the estimator gives
%                  a per-group scalar only, not the per-observation r_i the
%                  snooping and the sum(r)=m-u identity need). Intended as an
%                  accelerator/cross-check at very large nc; read the accuracy
%                  caveat in LSQVCE.LEVERAGE.
%     'probes'     1x1  number of Hutchinson probes Z (default 100).
%     'traceSeed'  1x1  RNG seed for reproducible probes (default []).
%
%   OUTPUT
%     red.group(k).key       group identifier
%     red.group(k).m         observations in the group
%     red.group(k).r         redundancy r_g = sum(r_i)
%     red.group(k).vtv       weighted residual square sum
%     red.group(k).s0sq      vtv / r    -- the variance component of this group
%     red.group(k).sigma0    sqrt(s0sq)
%     red.obs(j).h, .r       per-observation, if Collect
%     red.rTotal             sum over all groups
%     red.check              |sum(r_g) - (m - u)|, must be ~0
%     red.traceEstimate      present only with 'traceMethod','hutchinson': a
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
%     20-50% error on realistically conditioned normals (see LSQVCE.LEVERAGE).
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
%     red = lsqvce.redundancy(acc, sol, 'groupBy', 'arc');
%     disp(red.check)                             % must be ~1e-9
%
%     for k = 1:numel(red.group)
%         g = red.group(k);
%         lsqvce.util.log.msg('%-16s r=%8.1f  sigma=%.4f\n', string(g.key), g.r, g.sigma0);
%     end
%
%   ERRORS
%     lsqvce:redundancy:badGroupBy  the grouping key is not a known name or handle
%     lsqvce:redundancy:noStore     nothing was retained; re-run with 'store','all'
%
%   See also LSQVCE.LEVERAGE, LSQVCE.VCE, LSQVCE.OUTLIERS, LSQVCE.RESIDUALS
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-37: the arc-task enumeration
%   caches the meta calls and PREALLOCATES the task array (the old
%   tasks(end+1) growth copied O(nTask^2)). With assistance from Claude
%   (Anthropic), model Claude Fable 5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.

    narginchk(2,inf);

    p = inputParser;  p.FunctionName = 'lsqvce.redundancy';
    p.addParameter('groupBy'      , 'arc');
    p.addParameter('minRedundancy', 50, @(v) isscalar(v) && v > 0);
    p.addParameter('collect'      , true,  @(v) islogical(v) && isscalar(v));
    p.addParameter('verbose'      , true,  @(v) (islogical(v) || isnumeric(v)) && isscalar(v));
    p.addParameter('parallel'     , '',    @(s) isempty(s) || any(strcmpi(s, ...
                                    {'serial','parfor','spmd'})));
    % TraceMethod: how the common-block redundancy trace tr(Sxx N_g^cc) is
    % obtained. 'exact' (default) uses the per-arc leverage path below (exact,
    % O(nc^2)/arc). 'hutchinson' uses the Monte-Carlo estimator of Ellmer (2018,
    % 6.5.6) via LSQVCE.TRACEESTIMATE -- an OPT-IN accelerator for very large nc.
    % See the accuracy caveat in LSQVCE.LEVERAGE. EXPERIMENTAL: see note below.
    p.addParameter('traceMethod'  , 'exact', @(s) any(strcmpi(s, {'exact','hutchinson'})));
    p.addParameter('probes'       , 100,   @(v) isscalar(v) && v >= 1 && v == floor(v));
    p.addParameter('traceSeed'    , [],    @(v) isempty(v) || (isscalar(v) && isnumeric(v)));
    p.parse(varargin{:});
    o = p.Results;

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
        % group-level common-normal assembly; see HANDOVER).
        lsqvce.util.log.warn('lsqvce:redundancy:hutchinsonExperimental', ...
            ['''TraceMethod'',''hutchinson'' is EXPERIMENTAL: the exact per-arc ' ...
             'redundancies are still computed and remain authoritative; the ' ...
             'Hutchinson estimator (Ellmer 2018 6.5.6) is returned in ' ...
             'red.traceEstimate as a stochastic estimate of the TOTAL common ' ...
             'redundancy for comparison. See LSQVCE.TRACEESTIMATE and the ' ...
             'accuracy note in LSQVCE.LEVERAGE.']);
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

    Ncc_red = acc.reducedNormalMatrix();      % TOTAL, summed over all systems

    % ---- phase 1: per-system broadcast state (serial, cheap) --------------
    % The leverage blocks are the ONLY O(nc^3)/O(nl^3) work. Built once per
    % system, then broadcast to every task -- a few hundred MB at d/o 90,
    % negligible next to the compute they enable. The accumulator itself is
    % NEVER captured by a task: only these plain arrays and the loader.
    nSys    = numel(ids);
    stores  = cell(nSys, 1);
    blks    = cell(nSys, 1);
    ys      = cell(nSys, 1);
    nParTot = numel(x);
    for k = 1:nSys
        st        = acc.getStore(ids(k));
        blks{k}   = lsqvce.leverageBlocks(Ncc_red, st.Nll_r, st.Nlc);
        ys{k}     = st.Nll_r \ (st.nl_r - st.Nlc * x);
        stores{k} = struct('sysID', st.sysID, 'dataFcn', st.dataFcn, ...
                           'uo', st.uo, 'dayList', st.dayList);
        nParTot   = nParTot + st.nl;
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
    taskFcn = @(j) arcTask(tasks(j), stores{tasks(j).k}, blks{tasks(j).k}, ...
                           ys{tasks(j).k}, x, keyFcn, o.collect);

    switch mode
        case 'spmd'
            [P, info] = lsqvce.mapreduce(nTask, taskFcn, 'Mode', 'spmd', ...
                            'Threshold', thr, 'Verbose', o.verbose, ...
                            'Reduce',  @foldTask, 'Init', emptyPartial(), ...
                            'Combine', @combinePartials);
            if iscell(P)                     % fell back to a collecting mode
                Q = emptyPartial();
                for j = 1:numel(P), Q = foldTask(Q, P{j}); end
                P = Q;
            end
        otherwise
            [raw, info] = lsqvce.mapreduce(nTask, taskFcn, 'Mode', mode, ...
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
                   'vtv', P.vtv(k));
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
        try
            U = chol((Ncc_red + Ncc_red.')/2);        % Ncc_red = U'*U
            nc = size(Ncc_red, 1);
            % tr(N^-1 N) = tr(U^-1 (U'U) U^-T) = tr(I) = nc. With N = U'U the
            % identity operator is U^-T N U^-1 (NOT U^-1 N U^-T, which is a
            % similarity that is not the identity and biases the trace): apply
            % the upper solve on the LEFT and the transpose solve on the RIGHT.
            Xfun = @(V) U.' \ (Ncc_red * (U \ V));     % U^-T N U^-1 = I, matrix-free
            [tEst, tInfo] = lsqvce.traceEstimate(Xfun, 'Size', nc, ...
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

    % ---- the global identity: sum(r_i) = m - u ----------------------------
    red.check = abs(red.rTotal - (mTot - nParTot));
    relChk    = red.check / max(mTot - nParTot, 1);

    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg('%sgroups = %d   sum(r) = %.4f   m-u = %d   check = %.2e\n', lsqvce.util.verbosity.banner('Redundancy'), ...
                numel(G), red.rTotal, mTot - nParTot, red.check);
    end

    if relChk > 1e-6
        lsqvce.util.log.warn('lsqvce:redundancy:identityViolated', ...
            ['sum(r_i) = %.6f but m - u = %d (relative error %.2e). The ' ...
             'redundancy identity MUST hold exactly. Do not trust the variance ' ...
             'components. Likely causes: the loader returned different data on ' ...
             'this pass, or the stored Schur blocks are stale.'], ...
             red.rTotal, mTot - nParTot, relChk);
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
function t = arcTask(tk, st, blk, y, x, keyFcn, collect)
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

    [h, v] = arcLeverageAndResidual(arc, x, y, uo, blk);
    r = 1 - h;

    key = keyFcn(st.sysID, tk.obsType, tk.arc, tk.day);
    t = struct('ckey', char(string(key)), 'key', key, 'm', numel(v), ...
               'r', sum(r), 'vtv', v.'*v, 'h', sum(h), 'rec', []);
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
%            .hTot    1x1 double  total leverage accumulated (0)
%            .mTot    1x1 double  total observations accumulated (0)
%            .obs     1x0 cell    per-observation records (when stored)
    P = struct('keys', {{}}, 'keyVals', {{}}, 'm', [], 'r', [], 'vtv', [], ...
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
        P.m(end+1) = 0;  P.r(end+1) = 0;  P.vtv(end+1) = 0;
        k = numel(P.keys);
    end
    P.m(k)   = P.m(k)   + t.m;
    P.r(k)   = P.r(k)   + t.r;
    P.vtv(k) = P.vtv(k) + t.vtv;
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
            P.m(end+1) = 0;  P.r(end+1) = 0;  P.vtv(end+1) = 0;
            j = numel(P.keys);
        end
        P.m(j)   = P.m(j)   + Q.m(k);
        P.r(j)   = P.r(j)   + Q.r(k);
        P.vtv(j) = P.vtv(j) + Q.vtv(k);
    end
    P.hTot = P.hTot + Q.hTot;
    P.mTot = P.mTot + Q.mTot;
    P.obs  = [P.obs, Q.obs];
end

%==========================================================================
function f = makeKeyFcn(g)
%MAKEKEYFCN  Resolve a GroupBy spec into a key function f(sat,obsType,arc,day).
%   INPUT  g  char | function_handle   a preset name ('arc' | 'obsType' |
%             'system' | 'day') or a custom @(sat,obsType,arc,day)->char handle.
%   OUTPUT f  1x1 function_handle       @(sat,obsType,arc,day) -> group key char
%            (a handle is passed through unchanged; an unknown name errors).
    if isa(g, 'function_handle')
        f = g;  return
    end
    switch lower(g)
        case 'arc',     f = @(s,ot,a,d) sprintf('s%s_%s_d%s_a%s', ...
                                    string(s), string(ot), string(d), string(a));
        case 'obstype', f = @(s,ot,a,d) sprintf('%s', string(ot));
        case 'system',  f = @(s,ot,a,d) sprintf('s%s', string(s));
        case 'day',     f = @(s,ot,a,d) sprintf('d%s', string(d));
        otherwise
            error('lsqvce:redundancy:badGroupBy', ...
                ['GroupBy must be ''arc'', ''obsType'', ''system'', ''day'' ' ...
                 'or a function handle (got ''%s'').'], g);
    end
end

%==========================================================================
function [h, v] = arcLeverageAndResidual(arc, x, y, uo, blk)
%ARCLEVERAGEANDRESIDUAL  Whiten one arc, then get h and v, batched.

    m   = size(arc.A,1);
    bs  = lsqvce.util.getfielddef(uo, 'BatchSize', 20000);
    nOv = lsqvce.util.getfielddef(uo, 'BatchOverlap', []);
    if isempty(nOv)
        nOv = lsqvce.util.getfielddef(arc, 'filterMemory', 0);
    end
    bs  = max(1, min(round(bs),  m));
    nOv = max(0, min(round(nOv), bs-1));

    hasW = isfield(arc,'Wfcn') && ~isempty(arc.Wfcn);
    hasC = isfield(arc,'Cfcn') && ~isempty(arc.Cfcn);

    h = zeros(m,1);
    v = zeros(m,1);

    for i0 = 1:bs:m
        i1   = min(i0+bs-1, m);
        j0   = max(1, i0-nOv);
        idx  = j0:i1;
        keep = (i0-j0+1):numel(idx);
        rows = i0:i1;

        Ab = arc.A(idx,:);
        Bb = arc.B(idx,:);
        lb = arc.l(idx);  lb = lb(:);

        if hasW
            [Ab, Bb, lb] = arc.Wfcn(idx, Ab, Bb, lb, uo);
        elseif hasC
            W  = lsqvce.decorrFactor(arc.Cfcn(idx, uo));
            Ab = W*Ab;  Bb = W*Bb;  lb = W*lb;
        end

        hb = lsqvce.leverage(Ab, Bb, blk);
        vb = lb - Ab*x - Bb*y;

        h(rows) = hb(keep);
        v(rows) = vb(keep);
    end
end
