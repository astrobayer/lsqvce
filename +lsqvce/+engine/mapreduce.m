function [out, info] = mapreduce(nTask, taskFcn, opt)
%LSQVCE.ENGINE.MAPREDUCE  Execution-strategy layer: serial | parfor | spmd.
%
%   [out, info] = LSQVCE.ENGINE.MAPREDUCE(nTask, taskFcn, 'Name', value, ...)
%
%   Runs taskFcn(j) for j = 1:nTask under one of three schedules and combines
%   the results. This is the ONE place in the toolbox where parallelism lives:
%   the algorithms (LSQVCE.ENGINE.BLOCKACCUMULATOR/addSystems, LSQVCE.ENGINE.REDUNDANCY) only describe
%   independent tasks and a reduction, and this layer decides how to run them.
%   Do not sprinkle parfor through the algorithms -- it couples parallelism to
%   the math and makes both untestable.
%
%   WHY THIS IS EXACT
%     Every parallelised quantity in this toolbox is a plain SUM (or ordered
%     collection) of independent per-task contributions -- proven in
%     dev/validation/v4verify.py: partial sums over any partition reconstruct the
%     serial result. Changing the schedule changes WHICH core computes a term,
%     never the formula of a term (its last bits can still move -- see
%     DETERMINISM).
%
%   INPUT
%     nTask    number of independent tasks
%     taskFcn  @(j) -> result of task j. Must be pure: no shared mutable
%              state, no handle-object mutation. Handle objects reachable
%              from the closure are SERIALISED to the workers, so keep the
%              captured state small (see LSQVCE.ENGINE.REDUNDANCY: broadcast the
%              leverage blocks, never the accumulator).
%
%   OPTIONS
%     'Mode'      'serial' (default) | 'parfor' | 'spmd'
%     'Threshold' minimum nTask for a parallel schedule (default 8). Below it
%                 pool startup + serialisation exceed the saving; falls back
%                 to serial.
%     'Reduce'    @(acc, taskResult) -> acc, applied IN TASK ORDER.
%                 [] (default): out is the {nTask x 1} cell of raw results.
%     'Init'      initial accumulator value for 'Reduce'
%     'Combine'   @(partial1, partial2) -> partial. REQUIRED for 'spmd': each
%                 worker reduces its own tasks with 'Reduce' in ascending task
%                 order, then partials are combined across workers. Must be
%                 associative.
%     'AutoPool'  start a pool if none is open (default true). If starting
%                 fails, falls back to serial with a warning.
%     'Verbose'   default false
%
%   OUTPUT
%     out   cell array of task results, or the reduced value
%     info  .mode      the schedule actually used
%           .fallback  '' or the reason a parallel request ran serially
%           .nWorkers  pool size used (1 for serial)
%
%   DETERMINISM -- READ THIS ONCE
%     'serial' and 'parfor' share the IDENTICAL reduction order: parfor only
%     scrambles which worker computes a task; the per-task results are then
%     reduced on the client in task order. Differences against serial can
%     therefore come only from the arithmetic inside a task, not from the
%     schedule -- and they do come: a worker runs single-threaded, the client
%     multithreaded, so a task's matrix products and large sums can differ in
%     their last bits. Hence parfor
%     matches serial to 1e-14 relative, NOT bitwise (test_parallel
%     (testMapreduceParforMatchesSerial); WP-140). 'spmd' reduces on the
%     workers (that is its point -- per-task results are never shipped), so
%     its summation order differs and results match serial only to ~1e-10
%     relative. Regression tests against 'parfor' or 'spmd' need a tolerance,
%     never isequal.
%
%   EXAMPLE
%     % sum of expensive per-task scalars, exactly like the redundancy pass
%     f = @(j) struct('r', someCostlyThing(j));
%     [tot, info] = lsqvce.engine.mapreduce(nArc, f, 'Mode', 'parfor', ...
%                       'Reduce', @(a, t) a + t.r, 'Init', 0);
%
%   ACCURACY
%     No arithmetic of its own: it schedules taskFcn and applies Reduce in
%     task order. The serial schedule equals the hand-written loop exactly
%     (test_parallel (testMapreduceSerialOrderedReduce)); parfor equals
%     serial to 1e-14 relative and spmd, which reduces on the workers, to
%     1e-10 (test_parallel (testMapreduceParforMatchesSerial,
%     testMapreduceSpmdMatchesSerialToTolerance); PCT only, HANDOVER 6.1.3).
%     WP-139.
%
%   ERRORS
%     lsqvce:mapreduce:badTask      a task descriptor is malformed
%     lsqvce:mapreduce:serialFallback  [warning, once per run] a parallel mode was
%                                   requested but no pool could be used (R5 M1)
%
%   See also LSQVCE.ENGINE.BLOCKACCUMULATOR, LSQVCE.ENGINE.REDUNDANCY, LSQVCE.ENGINE.VCE
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-140: the determinism paragraph
%   says why parfor is not bitwise (single-threaded workers, a multithreaded
%   client: 1e-14 in test_parallel); v4verify's path is dev/validation; the
%   2026-08-19 stamp moved after See also. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, every number traced to the test that pins it. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 2 (R7-41): the log
%   queue is handed to every parfor task and spmd block, and their warning
%   ids are merged into the run inventory in task order. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (c) c1 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

arguments
    nTask
    taskFcn
    opt.Mode (1,1) string {mustBeMember(opt.Mode, ["serial","parfor","spmd"])} = "serial"
    opt.Threshold (1,1) double {mustBeNonnegative} = 8
    opt.Reduce = []
    opt.Init = []
    opt.Combine = []
    opt.AutoPool (1,1) logical = true
    opt.Verbose (1,1) logical = false
end

    narginchk(2, inf);
    validateattributes(nTask, {'numeric'}, {'scalar','integer','nonnegative'}, ...
                       mfilename, 'nTask', 1);
    if ~isa(taskFcn, 'function_handle')
        error('lsqvce:mapreduce:badTask', 'taskFcn must be a function handle.');
    end

    % WP-39 stage (c) c1: typed `arguments` signature (canonical
    % UpperCamelCase; legacy casings keep matching -- MATLAB NV
    % semantics; enum VALUES case-sensitive lowercase).
    o = struct('mode', char(opt.Mode), 'threshold', opt.Threshold, 'reduce', opt.Reduce, 'init', opt.Init, 'combine', opt.Combine, 'autoPool', opt.AutoPool, 'verbose', opt.Verbose);

    mode = lower(o.mode);
    info = struct('mode', mode, 'fallback', '', 'nWorkers', 1);

    % ---- gate the parallel request -----------------------------------------
    mode0 = mode;
    if ~strcmp(mode, 'serial')
        [ok, why, nW] = poolReady(o.autoPool);
        if nTask < o.threshold
            mode = 'serial';
            info.fallback = sprintf('nTask = %d < Threshold = %d', nTask, o.threshold);
        elseif ~ok
            mode = 'serial';
            info.fallback = why;
            % R5 M1: a parallel request that runs serially is REPORTED once per
            % run, never only recorded in a struct nobody prints.
            lsqvce.util.log.warnOnce('lsqvce:mapreduce:serialFallback', ...
                'Parallel mode ''%s'' requested but %s; running serially.', mode0, why);
        elseif strcmp(mode, 'spmd') && isempty(o.combine)
            mode = 'parfor';
            info.fallback = '''spmd'' requires a ''Combine'' function; using parfor';
        else
            info.nWorkers = nW;
        end
        if ~isempty(info.fallback)
            if lsqvce.util.verbosity.progress(o.verbose)
                lsqvce.util.log.msg('%s%s -> running %s\n', lsqvce.util.verbosity.banner('Parallel execution'), info.fallback, mode);
            end
        end
    end
    info.mode = mode;
    doReduce  = ~isempty(o.reduce);

    % ---- execute -----------------------------------------------------------
    switch mode
        case 'serial'
            if doReduce
                out = o.init;
                for j = 1:nTask
                    out = o.reduce(out, taskFcn(j));
                end
            else
                out = cell(nTask, 1);
                for j = 1:nTask
                    out{j} = taskFcn(j);
                end
            end

        case 'parfor'
            raw = cell(nTask, 1);
            % Worker-aware logging: install a DataQueue sink so any
            % lsqvce.util.log.msg/.warn on a worker lands in the CLIENT's log
            % file, in arrival order (safe no-op without the PCT). Removed after.
            logQ = lsqvce.util.log.startParallel();
            cleanupQ = onCleanup(@() lsqvce.util.log.stopParallel());  %#ok<NASGU>
            % feval(taskFcn, j), not taskFcn(j): R2026a's parfor otherwise
            % classifies the referenced handle as a SLICED variable and rejects
            % it ("sliced variable must not be a function handle"). feval keeps
            % taskFcn a broadcast variable, which is what we want.
            % R7-41: the queue is HANDED to every task (a worker cannot see the
            % client's persistent state) and every task returns the warning ids
            % it fired; the client merges them in task order, the order a
            % serial run fires them in. Before, logQ was never referenced in the
            % body: the sink was inert and worker warnings never reached
            % R.meta.warnings.
            wids = cell(nTask, 1);
            parfor j = 1:nTask
                [raw{j}, wids{j}] = workerTask(taskFcn, j, logQ);  %#ok<PFBNS>
            end
            for j = 1:nTask
                lsqvce.util.log.mergeInventory(wids{j});
            end
            clear cleanupQ                                       % drain & remove sink
            if doReduce
                out = o.init;
                for j = 1:nTask                    % ordered: same as serial
                    out = o.reduce(out, raw{j});
                end
            else
                out = raw;
            end

        case 'spmd'
            red  = o.reduce;  ini = o.init;  cmb = o.combine;
            logQ = lsqvce.util.log.startParallel();                  % R7-41, as parfor
            cleanupQ = onCleanup(@() lsqvce.util.log.stopParallel());  %#ok<NASGU>
            spmd
                onW = lsqvce.util.log.attachWorker(logQ);
                loc = ini;
                for j = spmdIndex:spmdSize:nTask   % worker-local, ascending
                    loc = red(loc, feval(taskFcn, j));
                end
                loc = spmdReduce(cmb, loc);        % associative tree combine
                wid = {};
                if onW, wid = lsqvce.util.log.inventory(); end
            end
            out = loc{1};
            for w = 1:numel(wid)
                lsqvce.util.log.mergeInventory(wid{w});
            end
            clear cleanupQ
    end
end

%==========================================================================
function [ok, why, nW] = poolReady(autoPool)
%POOLREADY  Is a parallel pool usable (and open it if allowed)?
    ok = false;  why = '';  nW = 1;
    if exist('gcp', 'file') ~= 2 || ...
            ~license('test', 'Distrib_Computing_Toolbox')
        why = 'Parallel Computing Toolbox not available';
        return
    end
    try
        pool = gcp('nocreate');
        if isempty(pool)
            if ~autoPool
                why = 'no pool open and AutoPool = false';
                return
            end
            pool = parpool;                        % default profile
        end
        nW = pool.NumWorkers;
        ok = true;
    catch me
        why = sprintf('could not open a pool (%s)', me.message);
    end
end

%==========================================================================
function [out, ids] = workerTask(taskFcn, j, q)
%WORKERTASK  One parfor task with the log sink attached (R7-41): the worker
%   forwards its log text to the client's queue q and returns the warning
%   identifiers this task fired, for lsqvce.util.log.mergeInventory. On the
%   client (a parfor without a pool) it is a plain call and returns {}.
    onW = lsqvce.util.log.attachWorker(q);
    if onW   % detach afterwards: a later log call on this worker must not send to a stale queue
        cu = onCleanup(@() lsqvce.util.log.attachWorker([])); %#ok<NASGU>
    end
    out = feval(taskFcn, j);
    ids = {};
    if onW, ids = lsqvce.util.log.inventory(); end
end

