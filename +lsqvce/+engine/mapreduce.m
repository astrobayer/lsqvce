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
%     validation/v4verify.py: partial sums over any partition reconstruct the
%     serial result. Changing the schedule changes WHICH core computes a term,
%     never the terms themselves.
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
%     therefore come only from multithreaded BLAS inside a task, not from the
%     schedule. 'spmd' reduces on the workers (that is its point -- per-task
%     results are never shipped), so its summation order differs and results
%     match serial only to ~1e-10 relative, NOT bit-identical. Regression
%     tests against 'spmd' need a tolerance, never isequal.
%
%   EXAMPLE
%     % sum of expensive per-task scalars, exactly like the redundancy pass
%     f = @(j) struct('r', someCostlyThing(j));
%     [tot, info] = lsqvce.engine.mapreduce(nArc, f, 'Mode', 'parfor', ...
%                       'Reduce', @(a, t) a + t.r, 'Init', 0);
%
%   ERRORS
%     lsqvce:mapreduce:badTask      a task descriptor is malformed
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (c) c1 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   See also LSQVCE.ENGINE.BLOCKACCUMULATOR, LSQVCE.ENGINE.REDUNDANCY, LSQVCE.ENGINE.VCE
%
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
    if ~strcmp(mode, 'serial')
        [ok, why, nW] = poolReady(o.autoPool);
        if nTask < o.threshold
            mode = 'serial';
            info.fallback = sprintf('nTask = %d < Threshold = %d', nTask, o.threshold);
        elseif ~ok
            mode = 'serial';
            info.fallback = why;
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
            parfor j = 1:nTask
                raw{j} = feval(taskFcn, j);                      %#ok<PFBNS>
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
            spmd
                loc = ini;
                for j = spmdIndex:spmdSize:nTask   % worker-local, ascending
                    loc = red(loc, feval(taskFcn, j));
                end
                loc = spmdReduce(cmb, loc);        % associative tree combine
            end
            out = loc{1};
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
