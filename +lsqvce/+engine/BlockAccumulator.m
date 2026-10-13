classdef BlockAccumulator < handle
%LSQVCE.ENGINE.BLOCKACCUMULATOR  Batched least squares with per-system pre-elimination.
%
%   Solves the multi-system adjustment
%
%       l_s = A_s*x_c + B_s*y_s + v_s ,      C{v_s} = sigma0^2 * Q_s
%
%   where x_c are COMMON parameters shared by every system (gravity field
%   coefficients, ...) and y_s are LOCAL parameters specific to system s
%   (accelerometer bias and scale, ambiguities, ...). The local parameters are
%   pre-eliminated by Schur complement after each system, so the global normal
%   matrix stays nc x nc no matter how many systems or local parameters exist.
%
%   Data arrive through a lazy callback. Nothing is precomputed and nothing is
%   retained: peak memory is one arc of (A,B,l) plus one batch covariance.
%
%   CONSTRUCTOR
%     acc = LSQVCE.ENGINE.BLOCKACCUMULATOR(nc, 'Name', value, ...)
%
%   OPTIONS
%     'Backend'      'qr' (default) | 'normal'
%                    'qr' keeps cond(A); 'normal' squares it. Use 'qr' unless
%                    you have measured that the conditioning is benign.
%     'BatchSize'    max observation rows per batch (default 20000)
%                    With a filter covariance the cost is O(m*p), so this can
%                    be the whole arc. With a dense covariance it is O(b^3):
%                    b = 5000 is already a 200 MB matrix and ~40 GFLOP.
%     'BatchOverlap' rows of overlap between batches (default [] = auto)
%                    Leave empty. The accumulator then reads arc.filterMemory
%                    from the loader, which LSQVCE.ENGINE.COVMODEL sets. Forcing
%                    it to 0 with a filter covariance costs ~12 orders of
%                    magnitude of accuracy. See BATCH OVERLAP below.
%     'Reg'          'none' (default) | 'tikhonov' | 'tsvd' | 'lm' | 'auto'
%                    Regularisation of the LOCAL block, applied before the
%                    Schur complement. See LSQVCE.ENGINE.APPLYREGULARISATION.
%     'RegLambda'    Tikhonov/LM weight; 0 with 'tikhonov' triggers GCV.
%     'RegMatrix'    L in lambda*||L*y||^2 (default identity).
%     'TruncTol'     rank tolerance relative to lambda_max (default 1e-10)
%     'CondWarn'     cond(Nll) above this -> 'ill-conditioned' (default 1e8)
%     'CondAbort'    cond(Nll) above this -> 'rank-deficient'  (default 1e14)
%     'LocalBlockSize'  [] | k: the locals form independent blocks of k
%                    columns (arc-wise locals), so Nll is pre-eliminated
%                    block by block -- O(nBlk k^3) instead of O(nl^3); used
%                    only when k divides nl (default [] = dense)
%     'KeepLocal'    keep the Schur blocks so the local parameters can be
%                    recovered and residuals computed (default true)
%     'Verbose'      print one line per system (default true)
%     'Parallel'     'serial' (default) | 'parfor' | 'spmd'
%                    Execution strategy for addSystems and for the
%                    redundancy/VCE data pass (LSQVCE.ENGINE.REDUNDANCY, LSQVCE.ENGINE.VCE).
%                    'serial' is the tested reference. See PARALLELISM below.
%     'ParallelThreshold'  minimum task count before a parallel schedule is
%                    used (default 8); below it, pool overhead loses.
%     'DistributeDays'  false (default). With true and 'Parallel' active,
%                    addSystem distributes the accumulation over DAYS: one raw
%                    partial system per day, merged in day order. This is the
%                    right axis for few-systems / many-days campaigns
%                    (GRACE-FO: 2 satellites, hundreds of days), where the
%                    addSystems ceiling (= nSys) is useless. The QR backend
%                    merges compressed triangular factors directly -- no
%                    normal-form detour, the conditioning advantage survives
%                    (validation/v5verify.py). addSystems disables it inside
%                    its workers: system-level OR day-level, never nested.
%     'UserOpt'      struct forwarded verbatim to every loader call
%
%   METHODS
%     addSystem(sysID, nl, dataFcn, dayList)   accumulate and pre-eliminate
%     addSystems(specs)                        several systems, parallelisable
%     sol = solve()                            solve for the common parameters
%     y   = backSubstitute(sysID, x)           recover one system's locals
%     Cyy = localCovariance(sysID, Cxx, s0sq)  their covariance (sol.Cxx, sol.s0sq)
%     res = residuals(sol, ...)                second pass over the data
%     tbl = report()                           per-system diagnostics table
%
%   THE LOADER CONTRACT
%     dataFcn is called with exactly two signatures:
%
%       meta = dataFcn('meta', day, [], opt)
%         .obsTypes  {1 x nOT} cellstr, e.g. {'pos','kbr'}
%         .arcs      {1 x nOT} cell, arcs{k} = arc IDs for obs type k
%
%       arc = dataFcn(day, obsType, arcID, opt)
%         .A             m x nc   common design
%         .B             m x nl   local design
%         .l             m x 1    observed minus computed
%         .t             m x 1    epochs (optional)
%         .Wfcn          @(idx,A,B,l,opt) -> whitened A,B,l   (preferred)
%         .Cfcn          @(idx,opt) -> dense batch covariance (fallback)
%         .filterMemory  required BatchOverlap (set by LSQVCE.ENGINE.COVMODEL)
%
%     The accumulator visits day -> obsType -> arc and releases each arc
%     immediately. Return an empty arc (l = []) for a data gap; it is skipped.
%
%     opt carries UserOpt plus the injected context sysID, nc, nl, day,
%     obsType, arc, BatchSize. Put STATIC loader state (paths, sensor
%     parameters, covariance specs) in a cfg struct captured by the closure,
%     not in opt.
%
%   WHY THE LOADER IS A CALLBACK
%     A day of GRACE-FO data across several observation types does not fit in
%     memory alongside its covariance, and a whole campaign certainly does not.
%     The callback lets the accumulator pull one arc at a time. If even a
%     single day of A is too large, make the loader arc-aware and build only
%     the rows of the requested arc.
%
%   BATCH OVERLAP -- THE ONE THING TO GET RIGHT
%     A whitening filter restarts its state at every batch, so the first p rows
%     of each batch are transient. Overlapping by the filter memory p and
%     discarding the overlap makes the batched result identical to filtering
%     the whole arc. Measured (AR(3) noise, BatchSize = 100):
%
%         BatchOverlap = p (auto) ->  1.6e-15  relative error
%         BatchOverlap = 0        ->  1.3e-03  relative error
%
%     A factor of 1e12, and it fails silently. Leave 'BatchOverlap' unset.
%
%   PARALLELISM -- WHERE IT PAYS AND WHERE IT DOES NOT
%     The block-arrow structure makes the key quantities SUMS of independent
%     per-arc / per-system contributions -- a clean reduction, no
%     approximation (proven in validation/v4verify.py). Two places pay:
%
%       * The redundancy / VCE data pass (LSQVCE.ENGINE.REDUNDANCY), parallel over
%         arcs. This is where the hours go at global SH scale: ~4 h/iter at
%         d/o 90 serial; with 16 workers plausibly 30-40 min. Build your
%         workflow around parallelising THIS.
%       * addSystems, parallel over systems. Ceiling = number of systems
%         (2 for GRACE-FO) -- a rounding error next to the redundancy pass.
%       * addSystem with 'DistributeDays': parallel over days WITHIN a
%         system, one raw partial per day merged in day order (v4.1). The
%         axis that actually scales for a 2-satellite mission over months.
%         Day-partial grouping changes the summation order vs the fused
%         loop: agreement to ~1e-13 relative, not bit-identical
%         (measured, v5verify.py tests 2 and 6).
%
%     Deliberately NOT parallelised:
%       * batches within an arc: the filter overlap couples them; sequential
%         by construction.
%       * the dense chol/qr/\ of the reduced normals: already multithreaded
%         by BLAS/LAPACK. A parfor there would fight the BLAS threads.
%
%     'serial' and 'parfor' reduce in the SAME order (per-task partials,
%     combined on the client in task order), so they differ only through
%     multithreaded-BLAS rounding inside a task. 'spmd' reduces on the
%     workers and matches serial to ~1e-10 relative, not bit-identical --
%     regression tests against it need a tolerance, never isequal.
%
%   REGULARISATION BIASES THE COMMON BLOCK
%     Pre-elimination forms Ncc - Nlc'*inv(Nll)*Nlc, so any bias introduced
%     into inv(Nll) propagates into x_c. Damping the local block is NOT a free
%     way to dispose of a nuisance parameter. If a local parameter is
%     inestimable by construction, remove it. LSQVCE.ENGINE.LOCALDIAGNOSTICS returns
%     .nullDir, the combination the data does not constrain.
%
%   EXAMPLE
%     nc = 12;  nl = 6;  days = 1:30;
%
%     opt = struct('psdMode','ar', 'arOrder',40);
%     acc = lsqvce.engine.BlockAccumulator(nc, ...
%             'Backend','qr', 'BatchSize',20000, 'Reg','auto', 'UserOpt',opt);
%
%     for s = [1 2]
%         loader = @(d,ot,a,uo) myLoader(cfg, s, d, ot, a, uo);
%         acc.addSystem(s, nl, loader, days);
%     end
%
%     sol = acc.solve();                     % common parameters
%     y   = acc.backSubstitute(1, sol.x);    % satellite 1 bias / scale
%     res = acc.residuals(sol);              % second pass
%     disp(acc.report())
%
%     lsqvce.util.log.msg('sigma0 = %.4f (want ~1)\n', sol.sigma0);
%
%   ERRORS
%     The system declaration is wrong (raised at construction)
%       lsqvce:accum:cancellation     (warning) v'Pv on the prior-folded
%                                     solve destroyed by cancellation
%                                     (negative or > 8 digits lost) -> NaN,
%                                     superseded by lsqvce.engine.residuals (WP-49)
%       lsqvce:accum:noSystems        no system specs at all
%       lsqvce:accum:badSpecs         a spec is not a struct with the required
%                                     fields
%       lsqvce:accum:duplicateSystem  two systems share a sysID; ids must be
%                                     unique and numeric or a SINGLE char
%       lsqvce:accum:badLoader        a spec's loader is not a function handle
%       lsqvce:accum:noArcs           a system declares no arcs
%       lsqvce:accum:noDays           no day list and none derivable
%       lsqvce:accum:badMeta          the loader's meta() returned the wrong shape
%       lsqvce:accum:metaMismatch     meta() disagrees with the spec (obsTypes
%                                     or arc count)
%       lsqvce:accum:badThresholds    a conditioning threshold is not a positive
%                                     scalar
%
%     A day returned by YOUR buildDay does not fit the declaration
%       lsqvce:accum:badArc           an arc index is outside 1..nArc
%       lsqvce:accum:arcBadNc         the arc's A has the wrong column count
%       lsqvce:accum:arcBadNl         the arc's B has the wrong column count
%       lsqvce:accum:arcRowMismatch   A, B and l disagree on their row count
%       lsqvce:accum:arcNotFinite     the arc carries NaN or Inf
%       lsqvce:accum:notWholeEpochs   a multi-axis arc does not contain whole
%                                     epochs: rows MUST be epoch-major and the
%                                     arc must set nAxis = d
%
%     The regularisation prior is malformed
%       lsqvce:accum:badPrior         the prior is not a struct of the expected
%                                     shape
%       lsqvce:accum:priorSize        its dimension disagrees with nc
%       lsqvce:accum:priorNotPD       the prior covariance is not positive
%                                     definite
%
%     Asked for something the accumulator does not have
%       lsqvce:accum:noStore          residuals or diagnostics were requested but
%                                     nothing was retained ('store','none')
%       lsqvce:accum:unknownSystem    no system carries the requested sysID
%       lsqvce:accum:badX             the estimate handed back has the wrong length
%       lsqvce:accum:badCxx           likewise for the covariance
%       lsqvce:accum:localNotRedundant  an arc's local block has no redundancy,
%                                     so its variance factor is undefined
%
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-112 (R7-03, R7-16): on the
%   QR backend the whole path stays square-root -- locals eliminated on the
%   factor (QRSystem.reducedRows), contributions merged as rows, the prior
%   folded in as rows; the rank decision is QRSystem.RankTol on the pivots.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-08-31 (CEST, UTC+02:00), R4 G-07 (Level L2): nl = 0 is a
%   valid system (user decision) -- addSystem accepts it, preeliminate
%   returns the unreduced common block with consistent empty local shapes.
%   With assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (c) c2 (F-07):
%   typed `arguments` signatures, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:accum:blockSchurFailed              [warning] A local block of system X is singular; falling back to the dense Schu...
%     lsqvce:accum:localSummarySkipped           [warning] locals[X]
%     lsqvce:accum:pinvFallback                  [warning] System X
%
%   See also LSQVCE.ENGINE.COVMODEL, LSQVCE.ENGINE.RESIDUALS, LSQVCE.ENGINE.CHECKPSD, LSQVCE.ENGINE.QRSYSTEM
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-128: a local block eliminated on
%   the QR factor keeps that factor (store fields Rll, Rlc, rl); its cofactor,
%   solution and covariance come from it, never from the squared Nll. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-127: datumRows -- a datum may
%   span several prior groups (blind groups need orthogonality only to the
%   groups the data sees); Dm with pinv; arcNotFinite also for A and B. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-12 (CEST, UTC+02:00), R7 WP-119: datumPart tests each
%   declared prior GROUP (CommonPrior.group) -- a datum stacked with a
%   regularising group gets the minimum-constraint covariance too. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 4: systemStore's
%   help names lsqvce.engine.verdict. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 3 (R7-24):
%   reducedNormalMatrix names lsqvce.engine.BlockCofactor (leverageBlocks
%   retired). With assistance from Claude (Anthropic), model Claude Opus
%   5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-86 part 2 (R7-24):
%   commonFactor, the reduced commons factor of lsqvce.engine.BlockCofactor
%   (the QR backend's own square root, the prior folded in by QR). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34, R7-35,
%   R7-40): the deleted COVMODELFACTORY reference replaced, LocalBlockSize
%   documented, localCovariance takes the solve's s0sq. With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    properties (SetAccess = private)
        nc                          % number of common parameters
        backend = 'qr'              % 'normal' | 'qr'
        sys                         % global (common-only) accumulator
        opt                         % option struct
        diag  = struct([])          % per-system diagnostics
        store = struct([])          % per-system Schur blocks (if KeepLocal)
        prior = []                  % lsqvce.engine.CommonPrior on the common parameters
        windowMax = 0               % longest whitening window seen (rows), MEASURED
                                    % by accumArc (R5 M1, H-06): the covariance grid
                                    % is sized from this, never from a probe
    end

    methods
        %------------------------------------------------------------------
        function obj = BlockAccumulator(nc, opt)
        %BLOCKACCUMULATOR  Construct a streaming block-arrow accumulator over nc
        %   common parameters. See the classdef header for the full option list
        %   ('Backend','KeepLocal','Verbose',...). INPUT nc 1x1 common count;
        %   Name/Value options. OUTPUT obj the accumulator handle.
            % WP-39 c2: typed `arguments` signature (canonical UpperCamel --
            % unchanged names; legacy casings keep matching, MATLAB NV
            % semantics; enum VALUES case-sensitive lowercase). The r remap
            % keeps the ctor body untouched; Verbose keeps its historical
            % logical-or-level contract (untyped here, resolved below).
        arguments
            nc (1,1) double {mustBeInteger, mustBePositive}
            opt.Backend (1,1) string {mustBeMember(opt.Backend, ["normal","qr"])} = "qr"
            opt.BatchSize (1,1) double {mustBePositive} = 20000
            opt.BatchOverlap double {mustBeScalarOrEmpty} = []
            opt.Reg (1,1) string {mustBeMember(opt.Reg, ["none","tikhonov","tsvd","lm","auto"])} = "none"
            opt.RegLambda (1,1) double {mustBeNonnegative} = 0
            opt.RegMatrix double = []
            opt.TruncTol (1,1) double {mustBePositive, mustBeLessThan(opt.TruncTol,1)} = 1e-10
            opt.CondWarn (1,1) double {mustBePositive} = 1e8
            opt.CondAbort (1,1) double {mustBePositive} = 1e14
            opt.LocalBlockSize double {mustBeScalarOrEmpty} = []
            opt.KeepLocal (1,1) logical = true
            opt.Verbose = true
            opt.Parallel (1,1) string {mustBeMember(opt.Parallel, ["serial","parfor","spmd"])} = "serial"
            opt.ParallelThreshold (1,1) double {mustBeNonnegative} = 8
            opt.DistributeDays (1,1) logical = false
            opt.UserOpt (1,1) struct = struct()
        end
            r = struct('nc', nc, 'Backend', char(opt.Backend), ...
                       'BatchSize', opt.BatchSize, 'BatchOverlap', opt.BatchOverlap, ...
                       'Reg', char(opt.Reg), 'RegLambda', opt.RegLambda, ...
                       'RegMatrix', opt.RegMatrix, 'TruncTol', opt.TruncTol, ...
                       'CondWarn', opt.CondWarn, 'CondAbort', opt.CondAbort, ...
                       'LocalBlockSize', opt.LocalBlockSize, 'KeepLocal', opt.KeepLocal, ...
                       'Verbose', opt.Verbose, 'Parallel', char(opt.Parallel), ...
                       'ParallelThreshold', opt.ParallelThreshold, ...
                       'DistributeDays', opt.DistributeDays, 'UserOpt', opt.UserOpt);
            if r.CondWarn > r.CondAbort
                error('lsqvce:accum:badThresholds', ...
                      'CondWarn (%.1e) must not exceed CondAbort (%.1e).', ...
                      r.CondWarn, r.CondAbort);
            end

            % Resolve the verbosity ONCE, so obj.opt.Verbose is always a
            % numeric level and every downstream gate is a level comparison
            % rather than a logical test (a logical test would silence the
            % stage detail at the higher levels).
            r.Verbose   = lsqvce.util.verbosity.resolve(r.Verbose);

            obj.nc      = double(r.nc);
            obj.opt     = r;
            obj.backend = lower(r.Backend);
            obj.sys     = obj.newSys(obj.nc);
        end

        %------------------------------------------------------------------
        function addSystem(obj, sysID, nl, dataFcn, dayList)
        %ADDSYSTEM  Accumulate one system, then pre-eliminate its local block.
        %
        %   addSystem(sysID, nl, dataFcn, dayList)
        %
        %     sysID    numeric identifier (satellite number, ...)
        %     nl       number of local parameters of this system (0 allowed:
        %              the local block is then empty and pre-elimination is
        %              the identity -- R4 G-07)
        %     dataFcn  loader callback (see the class help)
        %     dayList  vector of day identifiers, passed straight to the loader
        %
        %   Systems are independent: the order in which they are added changes
        %   only the floating-point rounding, not the result.

            narginchk(5,5);
            % R4 G-07 (user decision 2026-08-31): nl = 0 is VALID -- a system
            % with no local unknowns skips the Schur step (preeliminate).
            validateattributes(nl, {'numeric'}, {'scalar','integer','nonnegative'}, ...
                               mfilename, 'nl', 3);
            if ~isa(dataFcn, 'function_handle')
                error('lsqvce:accum:badLoader', 'dataFcn must be a function handle.');
            end
            if isempty(dayList)
                error('lsqvce:accum:noDays', 'dayList is empty for system %s.', ...
                      string(sysID));
            end
            if ~isempty(obj.diag) && any([obj.diag.sysID] == sysID)
                error('lsqvce:accum:duplicateSystem', ...
                      'System %s has already been added.', string(sysID));
            end

            out = lsqvce.engine.BlockAccumulator.processSystem(obj.nc, obj.backend, obj.opt, ...
                                                sysID, nl, dataFcn, dayList);
            obj.commit(out);
        end

        %------------------------------------------------------------------
        function addSystems(obj, specs)
        %ADDSYSTEMS  Add several systems at once, optionally in parallel.
        %
        %   addSystems(specs)
        %
        %     specs   struct array, one element per system, with fields
        %               .sysID     numeric identifier
        %               .nl        number of local parameters
        %               .dataFcn   loader callback (see the class help)
        %               .dayList   vector of day identifiers
        %
        %   Semantically identical to calling addSystem once per element, in
        %   order. With 'Parallel', 'parfor' the per-system accumulation and
        %   pre-elimination run concurrently; the reduced contributions are
        %   then MERGED IN INPUT ORDER on the client, so the merge sequence --
        %   and with it the rounding -- is the same as the serial call. Order
        %   independence of the merge is a tested invariant anyway; the
        %   ordered merge just makes regression comparisons tight.
        %
        %   THE LOADER MUST BE WORKER-SAFE
        %     Each worker gets its own copy of the loader closure. A
        %     'persistent' cache inside the loader (as in LSQVCE.ENGINE.SYNTHLOADER)
        %     becomes worker-local automatically and is safe. A hand-rolled
        %     SHARED cache (a file, a global) would race -- do not use one.
        %
        %   CEILING
        %     The speed-up is bounded by the number of systems (2 for
        %     GRACE-FO; a constellation is where this pays). The big win at
        %     global SH scale is the redundancy/VCE pass -- see LSQVCE.ENGINE.REDUNDANCY.
        %
        %   'spmd' is not meaningful across a handful of systems and is run
        %   as 'parfor'.

            narginchk(2,2);
            need = {'sysID','nl','dataFcn','dayList'};
            if ~isstruct(specs) || ~all(isfield(specs, need))
                error('lsqvce:accum:badSpecs', ...
                    'specs must be a struct array with fields %s.', ...
                    strjoin(need, ', '));
            end
            newIDs = [specs.sysID];
            oldIDs = [];
            if ~isempty(obj.diag), oldIDs = [obj.diag.sysID]; end
            if numel(unique(newIDs)) < numel(newIDs) || ...
                    any(ismember(newIDs, oldIDs))
                error('lsqvce:accum:duplicateSystem', ...
                      'Duplicate system IDs in specs or already added.');
            end

            mode = lower(obj.opt.Parallel);
            if strcmp(mode, 'spmd'), mode = 'parfor'; end

            % The task must not touch obj: capture only plain values, so the
            % workers never serialise the accumulated state.
            ncL   = obj.nc;
            beL   = obj.backend;
            optL  = obj.opt;
            taskSilenced = ~strcmp(mode, 'serial');
            if taskSilenced
                % Under parfor/spmd the workers' lines interleave unpredictably,
                % so the task is silenced and an ORDERED summary is printed below
                % instead. In SERIAL there is nothing to interleave, and muting
                % here unconditionally is what made the level-2 per-day detail
                % unreachable: it is produced inside processSystem, which only
                % ever ran with Verbose forced off. Caught by the messaging
                % wiring test, not by any pass/fail count.
                optL.Verbose = false;
            end
            optL.DistributeDays = false;       % no nested pools: system-level
                                               % parfor OR day-level, not both
            task  = @(s) lsqvce.engine.BlockAccumulator.processSystem(ncL, beL, optL, ...
                             specs(s).sysID, specs(s).nl, ...
                             specs(s).dataFcn, specs(s).dayList, ...
                             lsqvce.util.getfielddef(specs(s), 'name', ''));

            [out, info] = lsqvce.engine.mapreduce(numel(specs), task, 'Mode', mode, ...
                              'Threshold', obj.opt.ParallelThreshold, ...
                              ... % mapreduce takes a LOGICAL; its one line
                              ... % reports the execution mode, which belongs
                              ... % with the other progress output.
                              'Verbose', lsqvce.util.verbosity.progress(obj.opt.Verbose));

            for s = 1:numel(specs)             % deterministic: input order
                obj.commit(out{s});
                % This line SUBSTITUTES for the one processSystem would have
                % printed itself, so emit it only when that one was suppressed.
                % Printing both would double every system in serial mode.
                if taskSilenced && lsqvce.util.verbosity.detail(obj.opt.Verbose)
                    d = out{s}.dg;
                    lsqvce.util.log.msg('%s', lsqvce.util.verbosity.fieldBlock(sprintf('System %s', string(d.sysID)), ...
                        {'observations',   sprintf('%d', out{s}.keep.mObs); ...
                         'local unknowns', sprintf('%d', d.nl); ...
                         'cond(Nll)',      sprintf('%.3e', d.condNll); ...
                         'regularisation', char(d.regUsed)}));
                end
            end
            if lsqvce.util.verbosity.progress(obj.opt.Verbose) && ~strcmp(info.mode, 'serial')
                lsqvce.util.log.msg('%s%d systems on %d workers (%s)\n', ...
                        lsqvce.util.verbosity.banner('Systems'), numel(specs), info.nWorkers, info.mode);
            end
        end

        %------------------------------------------------------------------
        function sol = solve(obj)
        %SOLVE  Solve the reduced system for the common parameters.
        %
        %   sol.x        nc x 1   common parameter estimate
        %   sol.Cxx      nc x nc  covariance, scaled by the a-posteriori sigma0^2
        %   sol.sigma0   scalar   a-posteriori std deviation of unit weight
        %   sol.s0sq     scalar   sigma0^2
        %   sol.dof      scalar   redundancy
        %   sol.vtv      scalar   weighted residual square sum
        %   sol.diag     struct   per-system diagnostics
        %
        %   SANITY CHECK
        %     sigma0 should be ~1. If not, the covariance model is wrong by
        %     that factor: sigma0 = 2 means the noise is twice what you told
        %     the accumulator. A sigma0 well below 1 usually means the same
        %     observation was accumulated twice, or the weights are too loose.

            if isempty(obj.diag)
                error('lsqvce:accum:noSystems', 'No systems were added.');
            end

            if isempty(obj.prior)
                [x, Cxx, s0sq, dof, vtv] = obj.sys.solve();
                priorAbsorbed = 0;  rK = NaN;
                % rcond of the reduced common normal (no common-side prior here).
                if lsqvce.util.verbosity.detail(obj.opt.Verbose)
                    Ncc0 = obj.sys.commonBlocks();
                    rc0  = rcond(full((Ncc0 + Ncc0.')/2));
                    rc1  = rc0;                       % no regularisation applied
                end
            elseif strcmp(obj.backend, 'qr')
                % R7 WP-112 (R7-03): the prior's pseudo-observations K x ~ x0
                % (weight 1/sigma) are stacked onto the factor as ROWS and the
                % augmented system re-triangularised -- x, vtv (the residual
                % row) and Cxx from the square root; nothing is squared.
                [x, Cxx, s0sq, dof, vtv, rK, priorAbsorbed] = obj.solveQRWithPrior();
                if lsqvce.util.verbosity.detail(obj.opt.Verbose)
                    [Nk, ~] = obj.prior.normals();
                    Ncc = obj.sys.commonBlocks();
                    rc0 = rcond(full((Ncc + Ncc.')/2));  rc1 = rcond(full((Ncc + Nk + (Ncc + Nk).')/2));
                end
            else
                % Fold the prior in as pseudo-observations, then solve.
                [Nk, nk] = obj.prior.normals();
                [Ncc, nc_, ltl, m] = obj.sys.commonBlocks();

                N = Ncc + Nk;
                n = nc_ + nk;
                N = (N + N.')/2;

                if lsqvce.util.verbosity.detail(obj.opt.Verbose)
                    rc0 = rcond(full((Ncc + Ncc.')/2));   % BEFORE the prior
                    rc1 = rcond(N);                       % AFTER the prior
                end

                [Lc, fail] = chol(N, 'lower');
                if fail > 0
                    error('lsqvce:accum:priorNotPD', ...
                        ['Even with the prior, the common normal matrix is not ' ...
                         'positive definite. The prior does not span the null ' ...
                         'space of the data.']);
                end

                x   = Lc.' \ (Lc \ n);
                vtv = ltl + obj.prior.x0.'*obj.prior.x0/obj.prior.sigma^2 - n.'*x;
                % WP-49: same rule and same guard as the backends -- one
                % failure mode, one identifier per namespace. Negative OR
                % > 8 digits lost against the total quadratic mass -> NaN;
                % lsqvce.engine.residuals() supersedes it from a second pass.
                lPlTot = ltl + obj.prior.x0.'*obj.prior.x0/obj.prior.sigma^2;
                poisoned = vtv < 0 || (lPlTot > 0 && ...
                    log10(max(lPlTot / max(vtv, realmin), 1)) > 8);
                if poisoned
                    lsqvce.util.log.warn('lsqvce:accum:cancellation', ...
                        ['v''Pv = %.3e on the prior-folded solve was ' ...
                         'destroyed by cancellation (total quadratic mass ' ...
                         '%.6e). RETURNING NaN; lsqvce.engine.residuals() ' ...
                         'recomputes it honestly and SUPERSEDES it.'], ...
                         vtv, lPlTot);
                    vtv = NaN;
                end

                [rK, priorAbsorbed] = obj.prior.redundancy(Lc);
                dof  = m + obj.prior.nk - obj.nc;
                s0sq = vtv / max(dof, 1);

                Ri  = Lc \ eye(obj.nc);
                Cxx = s0sq * (Ri.' * Ri);
                Cxx = (Cxx + Cxx.')/2;
            end

            % R7 WP-114 (R7-05): a DATUM prior (the data do not see its
            % directions) fixes x for any sigma, and so must the covariance --
            % the minimum-constraint cofactor, not s0^2 inv(N + K'K/sigma^2),
            % which carries sigma^2 K'(KK')^-2 K on the datum directions.
            dp = [];  if ~isempty(obj.prior), dp = obj.datumPart(); end
            if ~isempty(dp)
                Rmi = dp.Rm \ eye(obj.nc);
                Cxx = s0sq * (Rmi * Rmi.' - dp.Dm * dp.Dm.');
                Cxx = (Cxx + Cxx.')/2;
            end
            sol = struct('x', x, 'Cxx', Cxx, 's0sq', s0sq, 'sigma0', sqrt(s0sq), ...
                         'dof', dof, 'vtv', vtv, 'diag', obj.diag, ...
                         'priorAbsorbed', priorAbsorbed, 'priorRedundancy', rK);

            if lsqvce.util.verbosity.detail(obj.opt.Verbose)
                % Conditioning: rcond is MATLAB's cheap 1-norm reciprocal
                % condition estimate (0 = singular, 1 = perfectly conditioned).
                % Show before and, if a common-side prior regularised the system,
                % after -- so the effect of the regularisation is visible.
                if ~isempty(obj.prior)
                    lsqvce.util.log.msg(['%srcond(N) = %.3e  ->  common prior ' ...
                            '(regularisation)  ->  rcond(N_reg) = %.3e\n'], lsqvce.util.verbosity.banner('Conditioning'), rc0, rc1);
                else
                    lsqvce.util.log.msg('%srcond(N) = %.3e   (no regularisation applied)\n', ...
                        lsqvce.util.verbosity.banner('Conditioning'), rc0);
                end
                lsqvce.util.log.msg(['           nc = %d   dof = %d (redundancy = obs - params)   ' ...
                        'sigma0 = %.6g\n'], obj.nc, dof, sqrt(s0sq));
                lsqvce.util.log.msg(['           sigma0 is the a-posteriori variance factor; ' ...
                        'it should be ~1 (>1 => noise larger than the model says).\n']);
                if ~isempty(obj.prior)
                    lsqvce.util.log.msg(['           prior absorbed %.1f of %d DoF ' ...
                             '(r_K = %.1f)\n'], priorAbsorbed, obj.prior.nk, rK);
                    if priorAbsorbed < 0.01*obj.prior.nk
                        lsqvce.util.log.msg(['           the prior is doing almost nothing ' ...
                                 '-- you may not need it\n']);
                    elseif priorAbsorbed > 0.99*obj.prior.nk
                        lsqvce.util.log.msg(['           the prior is doing almost EVERYTHING ' ...
                                 '-- the data is barely contributing\n']);
                    end
                end
                if sqrt(s0sq) < 0.5 || sqrt(s0sq) > 2
                    lsqvce.util.log.msg(['           sigma0 is far from 1 -- the covariance ' ...
                             'model is off by roughly this factor.\n']);
                end
                % Compact per-system LOCAL parameter summary (estimates + sds),
                % written to the log after the common solve. Full vectors are
                % available via R.backSubstitute / acc.backSubstitute(sysID, x).
                obj.logLocalSummary(x, Cxx, s0sq);
            end
        end

        %------------------------------------------------------------------
        function logLocalSummary(obj, x, Cxx, s0sq)
        %LOGLOCALSUMMARY  Write a compact per-system summary of the LOCAL
        %   parameter estimates and their standard deviations to the log.
        %   Summary only (count, value range, mean sd) -- systems can carry
        %   thousands of locals; the full vectors come from backSubstitute.
        %   INPUT  x   nc x 1  common estimate;  Cxx nc x nc common covariance;
        %          s0sq 1x1 the solve's variance factor (R7-40)
        %   OUTPUT none (writes via lsqvce.util.log.msg)
            if ~obj.opt.KeepLocal, return, end
            ids = obj.systemIDs();
            for i = 1:numel(ids)
                try
                    y = obj.backSubstitute(ids(i), x);
                catch ME
                    % REPORTED, never skipped in silence (R5 M1, H-07): a
                    % local block that cannot be recovered is a finding.
                    lsqvce.util.log.warn('lsqvce:accum:localSummarySkipped', ...
                        ['locals[%s]: back-substitution failed (%s: %s); ' ...
                         'the summary line is omitted.'], string(ids(i)), ...
                        ME.identifier, ME.message);
                    continue
                end
                if isempty(y), continue, end
                sdText = '';
                if ~isempty(Cxx)                      % the cofactor is a CONDITION, not a catch
                    Cyy = obj.localCovariance(ids(i), Cxx, s0sq);
                    sd  = sqrt(max(real(diag(Cyy)), 0));
                    sdText = sprintf('   mean sd = %.3e', mean(sd));
                end
                lsqvce.util.log.msg(['           locals[%s]: nl = %d   ' ...
                        'est in [%.3e, %.3e]%s\n'], string(ids(i)), numel(y), ...
                        min(y), max(y), sdText);
            end
        end

        %------------------------------------------------------------------
        function y = backSubstitute(obj, sysID, x)
        %BACKSUBSTITUTE  Recover the local parameters of one system.
        %
        %   y = backSubstitute(sysID, x)   solves   Nll * y = nl - Nlc * x
        %
        %   Uses the REGULARISED Nll that was actually applied during
        %   pre-elimination, so y is consistent with the reduced normal
        %   equations that produced x. If the block was regularised, y carries
        %   that bias -- check acc.diag(k).regUsed.

            narginchk(3,3);
            s = obj.store(obj.findStore(sysID));
            if numel(x) ~= obj.nc
                error('lsqvce:accum:badX', 'x has %d elements, expected %d.', ...
                      numel(x), obj.nc);
            end
            if ~isempty(s.Rll)                          % WP-128: on the QR factor
                y = s.Rll \ (s.rl - s.Rlc * x(:));
            else
                y = s.Nll_r \ (s.nl_r - s.Nlc * x(:));
            end
        end

        %------------------------------------------------------------------
        function Cyy = localCovariance(obj, sysID, Cxx, s0sq)
        %LOCALCOVARIANCE  Covariance of one system's local parameters.
        %
        %       Cyy = localCovariance(sysID, Cxx, s0sq)   (sol.Cxx, sol.s0sq)
        %       Cyy = s0^2*inv(Nll) + inv(Nll)*Nlc * Cxx * Nlc'*inv(Nll)
        %
        %   s0sq is the SOLVE's variance factor, the one Cxx carries. R7-40: it
        %   used to be the store's s0sq_hint = l'l / dof after eliminating the
        %   locals only -- the commons unsolved, so it carried the whole common
        %   signal and the verbose locals line read ~300x too large (the
        %   user-facing residuals path was right; both now call this).
        %
        %   The second term propagates the uncertainty of the common block into
        %   the locals. Omitting it (as a naive inv(Nll) would) understates the
        %   local uncertainty, sometimes badly, because locals and commons are
        %   correlated through Nlc.

            narginchk(4,4);
            s = obj.store(obj.findStore(sysID));
            if ~isequal(size(Cxx), [obj.nc obj.nc])
                error('lsqvce:accum:badCxx', 'Cxx must be %d x %d.', obj.nc, obj.nc);
            end

            if ~isempty(s.Rll)                          % WP-128: on the QR factor,
                Ri = s.Rll \ eye(size(s.Rll, 1));       % never the squared Nll
                G = s.Rll \ s.Rlc;  Qll = Ri * Ri.';
            else
                G   = s.Nll_r \ s.Nlc;                     % nl x nc
                Qll = s.Nll_r \ eye(size(s.Nll_r,1));      % cofactor matrix
            end
            Cyy = s0sq * Qll + G * Cxx * G.';
            Cyy = (Cyy + Cyy.')/2;
        end

        %------------------------------------------------------------------
        function res = residuals(obj, sol, varargin)
        %RESIDUALS  Second pass over the data to compute residuals.
        %   See LSQVCE.ENGINE.RESIDUALS for the option list.
        %
        %   VERBOSITY IS INHERITED. LSQVCE.ENGINE.RESIDUALS defaults its own 'verbose'
        %   to TRUE, which is right for someone calling it directly but wrong for
        %   a call made THROUGH an accumulator: the caller already stated a level
        %   when it built this object, and a second pass over the data is not a
        %   separate conversation. Without this, an adjustVCE run with verbose
        %   off still emitted a per-day residual line for every system and day --
        %   544 lines in the 2026-07-18 suite log, the only output a "silent"
        %   run produced. An explicit 'Verbose' in varargin still wins.
            if ~any(strcmpi(varargin(1:2:end), 'Verbose'))
                varargin = [varargin, {'Verbose', obj.opt.Verbose}];
            end
            res = lsqvce.engine.residuals(obj, sol, varargin{:});
        end

        %------------------------------------------------------------------
        function red = redundancy(obj, sol, varargin)
        %REDUNDANCY  Per-group / per-observation redundancy.
        %   See LSQVCE.ENGINE.REDUNDANCY.
            red = lsqvce.engine.redundancy(obj, sol, varargin{:});
        end

        %------------------------------------------------------------------
        function out = outliers(obj, sol, varargin)
        %OUTLIERS  Data snooping. See LSQVCE.ENGINE.OUTLIERS.
            out = lsqvce.engine.outliers(obj, sol, varargin{:});
        end

        %------------------------------------------------------------------
        function N = reducedNormalMatrix(obj)
        %REDUCEDNORMALMATRIX  The TOTAL reduced common normal matrix.
        %
        %   Summed over ALL systems, and INCLUDING the prior if one was set.
        %   This is what LSQVCE.ENGINE.BLOCKCOFACTOR factors: the leverage of an
        %   observation depends on the JOINT inverse, so a single system's own
        %   Ncc is the wrong thing to use.
            N = obj.sys.commonBlocks();
            if ~isempty(obj.prior)
                N = N + obj.prior.normals();
            end
            N = (N + N.')/2;
        end

        %------------------------------------------------------------------
        function Rc = commonFactor(obj)
        %COMMONFACTOR  Upper-triangular Rc with Rc'*Rc = the TOTAL reduced
        %   common normal matrix, the prior included (WP-86 part 2: the commons
        %   factor of lsqvce.engine.BlockCofactor). On the QR backend it is the
        %   backend's own square root with the prior rows folded in by QR --
        %   never chol of the squared normals; on the normal backend chol of
        %   reducedNormalMatrix, the factor solve() uses. Positive diagonal.
            if strcmp(obj.backend, 'qr')
                Rc = obj.sys.squareRoot();
                if ~isempty(obj.prior)
                    Rc = triu(qr(full([Rc; obj.prior.K / obj.prior.sigma]), 0));
                    Rc = Rc(1:obj.nc, :);
                end
            else
                Rc = chol(obj.reducedNormalMatrix());   % the matrix solve() factored
            end
            sg = sign(diag(Rc));  sg(sg == 0) = 1;
            Rc = sg .* Rc;
        end

        %------------------------------------------------------------------
        function setPrior(obj, pri)
        %SETPRIOR  Attach a regularisation prior on the COMMON parameters.
        %
        %   acc.setPrior(lsqvce.engine.CommonPrior(K, x0, sigma))
        %
        %   The prior is applied in solve(); it does not change accumulation.
        %   Setting it AFTER accumulation is therefore fine, and lets a VCE
        %   iteration vary sigma without re-reading the data.
            if ~isempty(pri) && ~isa(pri, 'lsqvce.engine.CommonPrior')
                error('lsqvce:accum:badPrior', ...
                      'The prior must be an lsqvce.engine.CommonPrior (or [] to clear).');
            end
            if ~isempty(pri) && pri.nc ~= obj.nc
                error('lsqvce:accum:priorSize', ...
                      'The prior has %d columns but nc = %d.', pri.nc, obj.nc);
            end
            obj.prior = pri;
        end

        %------------------------------------------------------------------
        function info = systemStore(obj, sysID)
        %SYSTEMSTORE  Public read-only accessor to one system's replay handles.
        %
        %   info = acc.systemStore(sysID) returns struct(sysID, nl, dataFcn,
        %   dayList, uo): everything an external diagnostics adapter needs to
        %   RE-ENUMERATE and re-whiten that system's arcs (e.g. to build the
        %   full-parameter batches lsqvce.engine.verdict.diagnose consumes). Requires
        %   'KeepLocal', true. Returns the loader handle and options only --
        %   never the accumulated normals -- so it cannot corrupt the solve.
            s = obj.store(obj.findStore(sysID));
            info = struct('sysID', s.sysID, 'nl', s.nl, ...
                          'dataFcn', s.dataFcn, 'dayList', s.dayList, ...
                          'uo', s.uo);
        end

        %------------------------------------------------------------------
        function ids = systems(obj)
        %SYSTEMS  Public list of system IDs, in insertion order.
            if isempty(obj.store), ids = [];
            else,                  ids = [obj.store.sysID]; end
        end

        %------------------------------------------------------------------
        function tbl = report(obj)
        %REPORT  Per-system diagnostics as a table.
        %
        %   Triage: any row with status ~= 'ok' had a poorly conditioned local
        %   block, and regUsed says what was done about it.
            if isempty(obj.diag)
                tbl = table();  return
            end
            d   = obj.diag;
            tbl = table([d.sysID].', [d.nl].', [d.nArc].', [d.condNll].', ...
                        [d.rankDef].', {d.status}.', {d.regUsed}.', ...
                        [d.lambdaUsed].', ...
                        'VariableNames', {'sysID','nLocal','nArc','condNll', ...
                                          'rankDef','status','regUsed','lambda'});
        end
    end

    %======================================================================
    methods (Access = public, Hidden)   % callable by the lsqvce.engine.residuals
                                        % function; Hidden keeps them off the
                                        % public API surface (?function is not a
                                        % valid access grant -- a function is
                                        % not a class)
        function s = getStore(obj, sysID)
        %GETSTORE  The kept local-block store for system sysID (needs KeepLocal).
            s = obj.store(obj.findStore(sysID));
        end
        function noteWindow(obj, w)
        %NOTEWINDOW  Record the longest whitening window (rows) seen so far.
        %   Called once per accumulated day; read by the joint estimator (the covariance
        %   grid is sized from this MEASURED length instead
        %   of probing the loader (R5 M1, finding H-06).
        %   INPUT  w 1x1 double  rows of the window just whitened
            obj.windowMax = max(obj.windowMax, w);
        end

        function ids = systemIDs(obj)
        %SYSTEMIDS  Row vector of the sysIDs added so far ([] if none).
            if isempty(obj.store), ids = []; else, ids = [obj.store.sysID]; end
        end

        function dp = datumPart(obj)
        %DATUMPART  The minimum-constraint cofactor's pieces when the common
        %   prior holds a DATUM, else [] (R7 WP-114, R7-05; per group since
        %   WP-119). The unit is the declared GROUP (CommonPrior.group, one per
        %   data.priors(k)): group g is a datum (DATUMROWS) when the data AND
        %   every group the data can see are blind to its rows,
        %       |N Kg'| <= 1e-10 |N| |Kg|    (on the QR backend |R Kg'| <= 1e-8 |R| |Kg|
        %                                     on its factor, N never formed)
        %       |Ko Kg'| <= 1e-10 |Ko| |Kg|  (Ko: the non-blind groups' rows; a
        %                                     datum may span several groups, WP-127)
        %   -- with every term PSD, the same as Kg in the null space of N plus
        %   the other groups. Kd = the datum groups' rows, Ki = the rest,
        %   Ni = N + Ki'Ki: span(Kd') is the null space of Ni, so x is Ni's
        %   minimum-constraint solution whatever the datum's sigma, and its
        %   cofactor is pinv(Ni). Computed stably with the datum rows
        %   re-weighted to the data's scale, c = sqrt(|Ni|_F) / |Kd|_F
        %   (|R|_F / |Kd|_F on the QR backend):
        %       pinv(Ni) = inv(Rm'Rm) - Dm*Dm',  Rm'Rm = Ni + c^2 Kd'Kd,
        %       Dm = Kd' pinv(Kd Kd') / c
        %   (both terms O(1/|N|); subtracting sigma^2 Kd'(Kd Kd')^-2 Kd from the
        %   solve's cofactor instead cancels: 5e-5 off at sigma 1e6 in the
        %   mirror, mirrors/wp114/datum_cofactor.py, against 3e-16 here;
        %   mirrors/wp119/datum_split.py: 1e-15 with a regularising group).
        %   A group only PARTLY in the null space is not a datum -- its x moves
        %   with its sigma, so its covariance stays the solve's (a ridge on
        %   every parameter is a Bayesian prior on the null part, not a datum).
        %   OUTPUT  dp  [] | struct .Rm (nc x nc upper) .Dm (nc x nd)
            dp = [];
            if isempty(obj.prior), return, end
            [isDatum, S, nS] = obj.datumRows();
            if ~any(isDatum), return, end
            K = full(obj.prior.K) / obj.prior.sigma;
            isQR = strcmp(obj.backend, 'qr');
            Kd = K(isDatum, :);  Ki = K(~isDatum, :);
            if isQR
                c = nS / norm(Kd, 'fro');
                Rm = triu(qr([S; Ki; c * Kd], 0));
                Rm = Rm(1:obj.nc, :);
            else
                Ni = S + Ki.' * Ki;
                c = sqrt(norm(Ni, 'fro')) / norm(Kd, 'fro');
                [Rm, fail] = chol((Ni + Ni.')/2 + c^2 * (Kd.' * Kd));
                if fail > 0, return, end
            end
            dp = struct('Rm', Rm, 'Dm', (Kd.' * pinv(Kd * Kd.')) / c);   % pinv: WP-127
        end

        function [isDatum, S, nS] = datumRows(obj)
        %DATUMROWS  Which rows of the common prior are a DATUM (WP-127).
        %   A group is BLIND when the data cannot see its rows (|S Kg'| <=
        %   tolS |S| |Kg|, S the factor on QR, N on the normal backend); a blind
        %   group is a datum when it is also orthogonal to the rows of every
        %   NON-blind group. WP-119 asked orthogonality to ALL other groups, so
        %   a datum split over two groups (translation and rotation of a
        %   free network: Kt Kr' = -sum y_i, nonzero unless the net is centred)
        %   failed for both halves and the covariance fell back to
        %   inv(N + K'K/sigma^2), growing with sigma (0.9 at sigma 1, 9e5 at
        %   1e3; dev/validation/wp127_datum_split_mirror.py).
        %   OUTPUT  isDatum  nk x 1 logical over obj.prior's rows (all false
        %                    without a prior); S, nS the data block and its norm
            isDatum = false(0, 1);  S = [];  nS = 0;
            if isempty(obj.prior), return, end
            K = full(obj.prior.K) / obj.prior.sigma;
            grp = obj.prior.group(:);
            if strcmp(obj.backend, 'qr')                 % on the factor
                S = obj.sys.squareRoot();  tolS = 1e-8;
            else
                S = obj.sys.commonBlocks();  S = (S + S.')/2;  tolS = 1e-10;
            end
            nS = norm(S, 'fro');
            gs = unique(grp).';
            blind = false(size(grp));
            for g = gs
                in = grp == g;  Kg = K(in, :);
                blind(in) = norm(S * Kg.', 'fro') <= tolS * nS * norm(Kg, 'fro');
            end
            isDatum = false(size(grp));
            Ko = K(~blind, :);                           % the groups the data CAN see
            for g = gs
                in = grp == g;
                if ~blind(find(in, 1)), continue, end
                Kg = K(in, :);
                if ~isempty(Ko) && norm(Ko * Kg.', 'fro') > 1e-10 * norm(Ko, 'fro') * norm(Kg, 'fro')
                    continue
                end
                isDatum(in) = true;
            end
        end
    end

    %======================================================================
    %======================================================================
    methods (Static, Hidden)   % public so the mapreduce worker context can
                               % invoke it; hidden because it is not user API

        %------------------------------------------------------------------
        function out = processSystem(nc, backend, opt, sysID, nl, dataFcn, dayList, sysName)
        %PROCESSSYSTEM  The complete per-system pipeline as a PURE function.
        %
        %   Accumulate one system's arcs and pre-eliminate its local block.
        %   Static and side-effect free by design: it never touches the
        %   accumulator object, so LSQVCE.ENGINE.MAPREDUCE can run it on a worker
        %   without serialising any accumulated state. addSystem and
        %   addSystems both go through here -- one code path, one test
        %   surface.

            nl  = double(nl);

            uo           = opt.UserOpt;
            uo.sysID     = sysID;
            uo.nc        = nc;
            uo.nl        = nl;
            uo.BatchSize = opt.BatchSize;
            if ~isempty(opt.BatchOverlap)
                uo.BatchOverlap = opt.BatchOverlap;
            end

            dayList = dayList(:).';
            % The readable label now ARRIVES with the spec (2026-07-19) instead
            % of being sniffed out of UserOpt: the Loader is where the name is
            % set, so carrying it through the spec is the one path that cannot
            % go stale. Fall back to the sysID when no name was given.
            if nargin < 8 || isempty(sysName)
                sysName = char(string(sysID));
            else
                sysName = char(sysName);
            end
            % Also publish it on UserOpt so the SECOND pass over the data
            % (lsqvce.engine.residuals, which sees only the stored uo) can label its
            % lines the same way. One name, one source.
            uo.systemName = sysName;

            distrib = opt.DistributeDays && ~strcmpi(opt.Parallel, 'serial') ...
                      && numel(dayList) >= opt.ParallelThreshold;
            window = 0;                          % longest whitening window (H-06)

            if distrib
                % -- day-distributed: one RAW partial system per day. Proven:
                %    validation/v5verify.py (merge == fused), v6verify.py
                %    test 8 (tree-grouping independence -> spmd is safe).
                task = @(j) lsqvce.engine.BlockAccumulator.processDay(nc, backend, ...
                                sysID, nl, dataFcn, dayList(j), uo);
                if strcmpi(opt.Parallel, 'spmd')
                    % reduce on the workers: mergeRaw is associative in every
                    % derived quantity, so the tree shape does not matter
                    ini = struct('loc', ...
                              lsqvce.engine.BlockAccumulator.newSysFor(backend, nc + nl), ...
                              'nArc', 0, 'window', 0);
                    P = lsqvce.engine.mapreduce(numel(dayList), task, ...
                            'Mode', 'spmd', ...
                            'Threshold', opt.ParallelThreshold, ...
                            'Reduce',  @lsqvce.engine.BlockAccumulator.foldDayPart, ...
                            'Init',    ini, ...
                            'Combine', @lsqvce.engine.BlockAccumulator.foldDayPart, ...
                            'Verbose', false);
                    loc = P.loc;  nArc = P.nArc;
                    window = P.window;
                else
                    % parfor/serial: collect, then merge on the client in
                    % day order (deterministic rounding)
                    out  = lsqvce.engine.mapreduce(numel(dayList), task, ...
                               'Mode', opt.Parallel, ...
                               'Threshold', opt.ParallelThreshold, ...
                               'Verbose', false);
                    loc  = lsqvce.engine.BlockAccumulator.newSysFor(backend, nc + nl);
                    nArc = 0;
                    for j = 1:numel(dayList)
                        loc.mergeRaw(out{j}.loc);
                        nArc = nArc + out{j}.nArc;
                        window = max(window, out{j}.window);
                    end
                end
            else
                loc  = lsqvce.engine.BlockAccumulator.newSysFor(backend, nc + nl);
                nArc = 0;
                nDay = numel(dayList);
                stageDetail = lsqvce.util.verbosity.detail(opt.Verbose);
                % The {system, type} index needs a SYSTEM ordinal. A numeric
                % sysID is already one; anything else falls back to 1 rather
                % than inventing an ordering the caller did not give.
                if isnumeric(sysID) && isscalar(sysID)
                    sysIdx = double(sysID);
                else
                    sysIdx = 1;
                end
                if stageDetail && nDay > 0
                    lsqvce.util.log.msg('%s', ...
                        lsqvce.util.verbosity.heading('Accumulate the normal matrix'));
                end
                for d = dayList
                    tDay = tic;
                    nA0  = nArc;
                    [nAdded, perType] = lsqvce.engine.BlockAccumulator.accumDay(loc, ...
                                            nc, nl, dataFcn, d, sysID, uo);
                    nArc = nArc + nAdded;
                    window = max(window, lsqvce.engine.BlockAccumulator.windowOf(perType, uo));
                    if stageDetail
                        % ONE LINE PER TYPE, under a day header. There used to be
                        % a numbered "Accumulate observations of type X" step
                        % above each line, which only restated the type the line
                        % already names.
                        %
                        % The index is {system, type}. The FIRST element was
                        % hardcoded to 1, so every system's types printed as
                        % {1,k} -- system 2's lines claimed to belong to system 1.
                        lsqvce.util.log.msg('%s', lsqvce.util.verbosity.dayLine(d));
                        for t = 1:numel(perType)
                            lsqvce.util.log.msg('%s', lsqvce.util.verbosity.typeLine( ...
                                sysIdx, t, perType(t).obsType, sysName, '', ...
                                perType(t).seconds));
                        end
                    elseif lsqvce.util.verbosity.progress(opt.Verbose)
                        % Level 1 keeps the compact one-line-per-day form: at
                        % progress level the point is the shape of the run, not
                        % which observation type took how long.
                        lsqvce.util.log.msg(['%ssystem %s, day %s : ' ...
                            '%d arc(s) this day (%d total)   [%.2f s]\n'], ...
                            lsqvce.util.verbosity.banner('Normal equations'), ...
                            string(sysID), string(d), nArc - nA0, nArc, toc(tDay));
                    end
                end
                if lsqvce.util.verbosity.progress(opt.Verbose) && nDay > 0
                    lsqvce.util.log.msg(['%ssystem %s: %d day(s), %d arc(s) ' ...
                        'accumulated\n'], lsqvce.util.verbosity.banner('Normal equations'), ...
                        string(sysID), nDay, nArc);
                end
            end

            if nArc == 0
                error('lsqvce:accum:noArcs', ...
                    ['System %s contributed no observations: every arc was ' ...
                     'empty, or the loader returned no obsTypes.'], string(sysID));
            end

            [red, dg, keep] = lsqvce.engine.BlockAccumulator.preeliminate(loc, nl, sysID, nc, opt, sysName);

            dg.sysID = sysID;
            dg.nl    = nl;
            dg.nArc  = nArc;

            keep.sysID   = sysID;
            keep.nl      = nl;
            keep.dataFcn = dataFcn;
            keep.dayList = dayList;
            keep.uo      = uo;

            out = struct('red', red, 'dg', dg, 'keep', keep, 'window', window);
        end

        %------------------------------------------------------------------
        function [nArc, perType] = accumDay(loc, nc, nl, dataFcn, d, sysID, uo)
        %ACCUMDAY  Accumulate one day of one system into loc.
        %   The single code surface for the fused-serial AND the
        %   day-distributed paths.
        %
        %   OUTPUT
        %     nArc     1x1 double   arcs accumulated this day (data gaps skipped)
        %     perType  1xT struct   one entry per observation type actually seen,
        %                           for the level-2 stage detail:
        %                           .obsType 1xN char   the type
        %                           .nArc    1x1 double arcs accumulated for it
        %                           .seconds 1x1 double wall time spent on it
        %                           Collected unconditionally -- it is one tic
        %                           per obsType, not per arc, so it costs nothing
        %                           measurable and the caller decides whether to
        %                           print it.
            nArc = 0;
            meta = dataFcn('meta', d, [], uo);
            lsqvce.engine.BlockAccumulator.validateMeta(meta, d, sysID);

            perType = struct('obsType', {}, 'nArc', {}, 'seconds', {}, ...
                             'rows', {}, 'filterMemory', {});
            for ot = 1:numel(meta.obsTypes)
                uo.obsType = meta.obsTypes{ot};
                arcs = meta.arcs{ot};
                tType = tic;
                nBefore = nArc;
                rowsType = 0;  memType = 0;            % MEASURED (R5 M1, H-06)

                for a = arcs(:).'
                    uo.day = d;  uo.arc = a;
                    arc = dataFcn(d, uo.obsType, a, uo);
                    if isempty(arc) || ~isfield(arc,'l') || isempty(arc.l)
                        continue                             % data gap
                    end
                    lsqvce.engine.BlockAccumulator.validateArc(arc, nc, nl, d, ...
                                                uo.obsType, a, sysID);
                    lsqvce.engine.BlockAccumulator.accumArc(loc, arc, uo);
                    nArc = nArc + 1;
                    rowsType = rowsType + size(arc.l, 1);
                    memType  = max(memType, lsqvce.util.getfielddef(arc, 'filterMemory', 0));
                end                                          % arc released here

                perType(end+1) = struct('obsType', char(uo.obsType), ...
                                        'nArc', nArc - nBefore, ...
                                        'seconds', toc(tType), ...
                                        'rows', rowsType, ...
                                        'filterMemory', memType);  %#ok<AGROW>
            end
        end

        %------------------------------------------------------------------
        function out = processDay(nc, backend, sysID, nl, dataFcn, d, uo)
        %PROCESSDAY  One day as a pure task: a fresh RAW partial system.
        %   Returned whole (a QRSystem carries (dim+1)^2 doubles, a
        %   NormalEqSystem about the same) and merged on the client in day
        %   order via mergeRaw. Arcs are never split, so the filter state
        %   never crosses a task boundary.
            loc  = lsqvce.engine.BlockAccumulator.newSysFor(backend, nc + nl);
            [nArc, perType] = lsqvce.engine.BlockAccumulator.accumDay(loc, nc, nl, dataFcn, d, sysID, uo);
            out  = struct('loc', loc, 'nArc', nArc, ...
                          'window', lsqvce.engine.BlockAccumulator.windowOf(perType, uo));
        end

        %------------------------------------------------------------------
        function w = windowOf(perType, uo)
        %WINDOWOF  Longest whitening window (rows) a day's streams imply.
        %   covModel precomputes its autocovariance out to
        %   min(BatchSize + BatchOverlap, rows of the day's epoch vector), so
        %   the covariance grid must cover that length (R5 M1, finding H-06:
        %   measured from what was accumulated, never probed from the loader).
        %   INPUT  perType 1xT struct (.rows .filterMemory); uo the user opt
        %   OUTPUT w 1x1 double  rows; 0 when nothing was accumulated
            w = 0;
            bs = lsqvce.util.getfielddef(uo, 'BatchSize', 20000);
            for t = 1:numel(perType)
                nOv = lsqvce.util.getfielddef(uo, 'BatchOverlap', []);
                if isempty(nOv), nOv = perType(t).filterMemory; end
                w = max(w, min(bs + nOv, perType(t).rows));
            end
        end

        %------------------------------------------------------------------
        function P = foldDayPart(P, Q)
        %FOLDDAYPART  Associative combine of two day partials (spmd reduce).
        %   Grouping independence proven: validation/v6verify.py test 8.
            P.loc.mergeRaw(Q.loc);
            P.nArc = P.nArc + Q.nArc;
            P.window = max(lsqvce.util.getfielddef(P, 'window', 0), ...
                           lsqvce.util.getfielddef(Q, 'window', 0));
        end
    end

    %======================================================================
    methods (Access = private)
        %------------------------------------------------------------------
        function [x, Cxx, s0sq, dof, vtv, rK, absorbed] = solveQRWithPrior(obj)
        %SOLVEQRWITHPRIOR  The prior-folded solve on the QR backend, square-root
        %   throughout (R7 WP-112): [R; K/sigma, x0/sigma] re-triangularised,
        %   refused by QRSystem.RankTol on its pivots (priorNotPD), vtv the
        %   squared residual row (NaN when the backend's own vtv is poisoned).
        %   OUTPUT  as solve(), plus rK / absorbed (CommonPrior.redundancy).
            nc = obj.nc;
            w  = 1 / obj.prior.sigma;
            G  = triu(qr([obj.sys.R; full(obj.prior.K) * w, obj.prior.x0 * w], 0));
            G  = G(1:min(size(G, 1), nc + 1), :);
            dR = abs(diag(G(1:min(nc, size(G, 1)), 1:nc)));
            if numel(dR) < nc || min(dR) <= lsqvce.engine.QRSystem.RankTol * max(dR)
                error('lsqvce:accum:priorNotPD', ...
                    ['Even with the prior, the common normal matrix is not ' ...
                     'positive definite. The prior does not span the null ' ...
                     'space of the data.']);
            end
            R1  = G(1:nc, 1:nc);
            x   = R1 \ G(1:nc, end);
            vtv = 0;
            if size(G, 1) > nc, vtv = G(nc + 1, end)^2; end
            if obj.sys.vtvPoisoned, vtv = NaN; end           % warned at merge()
            [rK, absorbed] = obj.prior.redundancy(R1.');     % lower: R1' R1 = N
            dof  = obj.sys.m + obj.prior.nk - nc;
            s0sq = vtv / max(dof, 1);
            Ri   = R1 \ eye(nc);
            Cxx  = s0sq * (Ri * Ri.');
            Cxx  = (Cxx + Cxx.')/2;
        end


        %------------------------------------------------------------------
        function commit(obj, out)
        %COMMIT  Merge one processed system into the accumulator (serial).
            obj.sys.merge(out.red);
            obj.noteWindow(lsqvce.util.getfielddef(out, 'window', 0));
            if isempty(obj.diag), obj.diag = out.dg;
            else,                 obj.diag(end+1) = out.dg; end
            if obj.opt.KeepLocal
                if isempty(obj.store), obj.store = out.keep;
                else,                  obj.store(end+1) = out.keep; end
            end
        end

        %------------------------------------------------------------------
        function k = findStore(obj, sysID)
        %FINDSTORE  Index of system sysID in obj.store; errors if absent/not kept.
            if isempty(obj.store)
                error('lsqvce:accum:noStore', ...
                    ['The local blocks were not kept. Construct the accumulator ' ...
                     'with ''KeepLocal'', true to enable backSubstitute, ' ...
                     'localCovariance and residuals.']);
            end
            k = find([obj.store.sysID] == sysID, 1);
            if isempty(k)
                error('lsqvce:accum:unknownSystem', ...
                      'System %s was never added.', string(sysID));
            end
        end

        %------------------------------------------------------------------
        function s = newSys(obj, n)
        %NEWSYS  A fresh per-system solver of this accumulators backend, size n.
            s = lsqvce.engine.BlockAccumulator.newSysFor(obj.backend, n);
        end
    end

    %======================================================================
    methods (Static, Access = private)

        %------------------------------------------------------------------
        function s = newSysFor(backend, n)
        %NEWSYSFOR  Construct an n x n solver object for the named backend (normal|qr).
            switch backend
                case 'normal', s = lsqvce.engine.NormalEqSystem(n);
                case 'qr',     s = lsqvce.engine.QRSystem(n);
            end
        end

        %------------------------------------------------------------------
        function validateMeta(meta, d, sysID)
        %VALIDATEMETA  Assert a days meta struct has matching obsTypes/arcs; errors otherwise.
            if ~isstruct(meta) || ~isfield(meta,'obsTypes') || ~isfield(meta,'arcs')
                error('lsqvce:accum:badMeta', ...
                    ['System %s, day %s: the meta call must return a struct with ' ...
                     'fields ''obsTypes'' and ''arcs''.'], string(sysID), string(d));
            end
            if numel(meta.obsTypes) ~= numel(meta.arcs)
                error('lsqvce:accum:metaMismatch', ...
                    ['System %s, day %s: obsTypes has %d entries but arcs has %d.'], ...
                     string(sysID), string(d), numel(meta.obsTypes), numel(meta.arcs));
            end
        end

        %------------------------------------------------------------------
        function validateArc(arc, nc, nl, d, ot, a, sysID)
        %VALIDATEARC  Assert an arc struct has consistent A/B/l sizes for (nc,nl,d); errors otherwise.
            ctx = sprintf('system %s, day %s, obsType %s, arc %s', ...
                          string(sysID), string(d), string(ot), string(a));

            if ~isfield(arc,'A') || ~isfield(arc,'B')
                error('lsqvce:accum:badArc', '%s: the arc needs fields A, B and l.', ctx);
            end

            m = numel(arc.l);
            if size(arc.A,1) ~= m || size(arc.B,1) ~= m
                error('lsqvce:accum:arcRowMismatch', ...
                    ['%s: A has %d rows, B has %d rows, l has %d elements -- all ' ...
                     'three must agree.'], ctx, size(arc.A,1), size(arc.B,1), m);
            end
            if size(arc.A,2) ~= nc
                error('lsqvce:accum:arcBadNc', ...
                      '%s: A has %d columns but nc = %d.', ctx, size(arc.A,2), nc);
            end
            if size(arc.B,2) ~= nl
                error('lsqvce:accum:arcBadNl', ...
                      '%s: B has %d columns but nl = %d.', ctx, size(arc.B,2), nl);
            end
            if ~all(isfinite(arc.l))
                error('lsqvce:accum:arcNotFinite', ...
                    ['%s: l contains NaN or Inf. Remove the bad epochs in the ' ...
                     'loader -- the accumulator cannot know which are valid.'], ctx);
            end
            % WP-127: A and B too -- a NaN in A of a system without locals went
            % through the QR into every pivot and came out as x = NaN, blamed
            % on cancellation (the only finiteness check sat in the local block)
            if ~all(isfinite(arc.A(:))) || ~all(isfinite(arc.B(:)))
                error('lsqvce:accum:arcNotFinite', ...
                    ['%s: the design (A or B) contains NaN or Inf. Remove the bad ' ...
                     'epochs in the loader.'], ctx);
            end
        end

        %------------------------------------------------------------------
        function accumArc(acc, arc, uo)
        %ACCUMARC  Split one arc into row batches and accumulate.
        %
        %   The whitening filter restarts its state at every batch, so batches
        %   overlap by the filter memory and the overlap rows are discarded.
        %   With overlap >= p this is exactly equivalent to filtering the whole
        %   arc; with overlap = 0 the transient leaks into the normal equations.

            % R6 WP-91 (E-09): the batch bounds, the overlap, the whole-epoch
            % rounding and the whitening branches live in the ONE iterator
            % arcBatchWhiten (engine-private); this callback keeps the
            % accumulation order byte-identical.
            arcBatchWhiten(arc, uo, ...
                @(rows, keep, Ab, Bb, lb) acc.add([Ab(keep,:), Bb(keep,:)], lb(keep)));
        end

        %------------------------------------------------------------------
        function [red, dg, keep] = preeliminate(loc, nl, sysID, nc, opt, sysName)
            % sysName is the readable label for the log ONLY; everything
            % numerical keys off sysID. Optional so the older 5-argument
            % call still works and simply prints the id.
            if nargin < 6 || isempty(sysName), sysName = char(string(sysID)); end
        %PREELIMINATE  Schur complement of the local block.

            ic = 1:nc;                                       %#ok<*PROPLC>
            il = nc + (1:nl);

            [Nll, Nlc, Ncc, nl_rhs, nc_rhs, ltl, mObs] = loc.blocks(ic, il);

            if nl == 0
                % R4 G-07: no local block. The Schur complement is the
                % identity (Ncc_red = Ncc, dof_red = mObs); the stored blocks
                % are the consistent EMPTY shapes so backSubstitute returns a
                % 0x1 y, localCovariance a 0x0 Cyy, and residuals form
                % l - A*x - B*y with B m x 0. localDiagnostics /
                % applyRegularisation are bypassed (eig of an empty matrix is
                % not a diagnostic); dg carries their full field set so every
                % consumer (report, commit, verbose block) reads it unchanged.
                if mObs == 0
                    error('lsqvce:accum:localNotRedundant', ...
                        'System %s has no observations.', string(sysID));
                end
                dg = struct('eig', zeros(0,1), 'V', zeros(0,0), 'condNll', 1, ...
                            'rankDef', 0, 'spd', true, 'nullDir', zeros(0,1), ...
                            'diagDom', 1, 'status', 'ok', 'regUsed', 'none', ...
                            'lambdaUsed', 0);
                Ncc  = (Ncc + Ncc.')/2;
                if isa(loc, 'lsqvce.engine.QRSystem')     % R7 WP-112: the factor itself
                    red = lsqvce.engine.ReducedContribution([], [], [], mObs, loc.reducedRows(0));
                else
                    red = lsqvce.engine.ReducedContribution(Ncc, nc_rhs, ltl, mObs);
                end
                keep = struct('Nll_r', zeros(0,0), 'Nlc', zeros(0,nc), ...
                              'nl_r', zeros(0,1), 'ltl', ltl, 'mObs', mObs, ...
                              's0sq_hint', max(ltl,0)/max(mObs,1), ...
                              'Rll', zeros(0), 'Rlc', zeros(0, nc), 'rl', zeros(0, 1));
                if lsqvce.util.verbosity.detail(opt.Verbose)
                    lsqvce.util.log.msg('%s', lsqvce.util.verbosity.fieldBlock(sprintf('System %s "%s"', string(sysID), sysName), ...
                        {'observations',   sprintf('%d', mObs); ...
                         'local unknowns', '0 (no pre-elimination)'}));
                end
                return
            end

            if mObs <= nl
                error('lsqvce:accum:localNotRedundant', ...
                    ['System %s has %d observations but %d local parameters. The ' ...
                     'local block cannot be pre-eliminated.'], string(sysID), mObs, nl);
            end

            dg = lsqvce.engine.localDiagnostics(Nll, opt);
            [Nll_r, nl_r, dg] = lsqvce.engine.applyRegularisation(Nll, nl_rhs, dg, opt);

            % R7 WP-112 (R7-03): on the QR backend an UNREGULARISED local block
            % is eliminated on the factor (QRSystem.reducedRows) -- the reduced
            % contribution leaves as square-root rows and merge never squares
            % it. The normal-form Schur below remains for a regularised block
            % (the regularisation is defined on Nll) and for a local factor
            % whose pivots fail QRSystem.RankTol (the pinv fallback's case).
            if isa(loc, 'lsqvce.engine.QRSystem') && strcmp(dg.regUsed, 'none')
                [rows, R11, R12, r1] = loc.reducedRows(nl);
                dR = abs(diag(R11));
                if numel(dR) == nl && min(dR) > lsqvce.engine.QRSystem.RankTol * max(dR)
                    red  = lsqvce.engine.ReducedContribution([], [], [], mObs - nl, rows);
                    keep = struct('Nll_r', Nll_r, 'Nlc', Nlc, 'nl_r', nl_r, ...
                                  'ltl', ltl, 'mObs', mObs, ...
                                  's0sq_hint', max(red.ltl, 0)/max(mObs - nl, 1), ...
                                  'Rll', R11, 'Rlc', R12, 'rl', r1);   % WP-128: the exact factor
                    if lsqvce.util.verbosity.detail(opt.Verbose)
                        hdr = sprintf('System %s "%s"', string(sysID), sysName);
                        lsqvce.util.log.msg('%s', lsqvce.util.verbosity.fieldBlock(hdr, ...
                            {'observations',    sprintf('%d', mObs); ...
                             'local unknowns',  sprintf('%d (eliminated on the QR factor)', nl); ...
                             'cond(Nll)',       sprintf('%.3e', dg.condNll)}));
                    end
                    return
                end
            end

            % Arc-wise local parameters make Nll BLOCK-DIAGONAL: arc a only
            % touches its own columns of B. Exploiting that turns an O(nl^3)
            % dense inversion into nBlk independent O(blk^3) solves.
            % Measured at 120 arcs x 6 params (nl = 720): 21x faster, 120x less
            % memory, agreement 7e-16. At full mission length the dense path is
            % simply not viable while the block path stays flat.
            blkSize = opt.LocalBlockSize;
            useBlk  = ~isempty(blkSize) && blkSize > 0 && mod(nl, blkSize) == 0 ...
                      && nl > blkSize;

            % Pre-elimination must never abort. Ordinary ill conditioning is the
            % QR backend's job and the regularisation above; only a genuinely
            % singular block reaches the fallback.
            RHS = [Nlc, nl_r];

            if useBlk
                nBlk = nl / blkSize;
                Z    = zeros(nl, nc+1);
                ok   = true;
                for b = 1:nBlk
                    ib = (b-1)*blkSize + (1:blkSize);
                    Nb = Nll_r(ib, ib);
                    if rcond(Nb) < eps
                        ok = false;  break
                    end
                    Z(ib,:) = Nb \ RHS(ib,:);
                end
                if ok
                    Ncc_red = Ncc    - Nlc.'*Z(:,1:nc);
                    ncc_red = nc_rhs - Nlc.'*Z(:,end);
                    ltl_red = ltl    - nl_r.'*Z(:,end);
                    dof_red = mObs   - nl;
                    Ncc_red = (Ncc_red + Ncc_red.')/2;
                    red  = lsqvce.engine.ReducedContribution(Ncc_red, ncc_red, ltl_red, dof_red);
                    keep = struct('Nll_r', Nll_r, 'Nlc', Nlc, 'nl_r', nl_r, ...
                                  'ltl', ltl, 'mObs', mObs, ...
                                  's0sq_hint', max(ltl_red,0)/max(dof_red,1), ...
                                  'Rll', zeros(0), 'Rlc', zeros(0, nc), 'rl', zeros(0, 1));
                    if lsqvce.util.verbosity.detail(opt.Verbose)
                        lsqvce.util.log.msg('%s', lsqvce.util.verbosity.fieldBlock(sprintf('System %s "%s"', string(sysID), sysName), ...
                            {'observations',    sprintf('%d', mObs); ...
                             'local unknowns',  sprintf('%d (%d blocks of %d)', nl, nBlk, blkSize); ...
                             'cond(Nll)',       sprintf('%.3e', dg.condNll); ...
                             'regularisation',  char(dg.regUsed)}));
                    end
                    return
                end
                lsqvce.util.log.warn('lsqvce:accum:blockSchurFailed', ...
                    ['A local block of system %s is singular; falling back to the ' ...
                     'dense Schur complement.'], string(sysID));
            end

            if ~all(isfinite(Nll_r(:))) || rcond(Nll_r) < eps
                Z = pinv(Nll_r, opt.TruncTol*norm(Nll_r)) * RHS;
                dg.regUsed    = 'pinv-fallback';
                dg.lambdaUsed = opt.TruncTol;
                lsqvce.util.log.warn('lsqvce:accum:pinvFallback', ...
                    ['System %s: the local block is still singular after ''%s'' ' ...
                     'regularisation (cond %.3e). Falling back to the pseudo-' ...
                     'inverse. The local parameters are not separately estimable; ' ...
                     'the weakest direction is\n    [%s]\nConsider removing a ' ...
                     'parameter rather than damping it.'], ...
                     string(sysID), opt.Reg, dg.condNll, ...
                     strjoin(compose('%+.3f', dg.nullDir(:).'), ' '));
            else
                Z = Nll_r \ RHS;
            end

            Ncc_red = Ncc    - Nlc.'*Z(:,1:nc);
            ncc_red = nc_rhs - Nlc.'*Z(:,end);
            ltl_red = ltl    - nl_r.'*Z(:,end);
            dof_red = mObs   - nl;

            Ncc_red = (Ncc_red + Ncc_red.')/2;

            red = lsqvce.engine.ReducedContribution(Ncc_red, ncc_red, ltl_red, dof_red);

            keep = struct('Nll_r', Nll_r, 'Nlc', Nlc, 'nl_r', nl_r, ...
                          'ltl', ltl, 'mObs', mObs, ...
                          's0sq_hint', max(ltl_red, 0)/max(dof_red, 1), ...
                          'Rll', zeros(0), 'Rlc', zeros(0, nc), 'rl', zeros(0, 1));

            if lsqvce.util.verbosity.detail(opt.Verbose)
                lsqvce.util.log.msg('%s', lsqvce.util.verbosity.fieldBlock(sprintf('System %s "%s"', string(sysID), sysName), ...
                    {'observations',    sprintf('%d', mObs); ...
                     'local unknowns',  sprintf('%d', nl); ...
                     'cond(Nll)',       sprintf('%.3e', dg.condNll); ...
                     'rank defect',     sprintf('%d', dg.rankDef); ...
                     'regularisation',  sprintf('%s (lambda = %.3e)', ...
                                                char(dg.regUsed), dg.lambdaUsed)}));
            end
        end
    end
end
