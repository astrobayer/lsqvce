function R = adjustBlock(systems, days, opt)
%LSQVCE.ENGINE.ADJUSTBLOCK  One-call block-arrow adjustment: accumulate, solve, judge.
%   R = lsqvce.engine.adjustBlock(systems, days)
%   R = lsqvce.engine.adjustBlock(systems, days, Name, Value, ...)
%
%   The single front door over the streaming pipeline. It builds a
%   LSQVCE.ENGINE.BLOCKACCUMULATOR, adds every system, solves the reduced common
%   system, bridges the result into the diagnostics (LSQVCE.ENGINE.VERDICT.SOURCE /
%   LSQVCE.ENGINE.VERDICT.DIAGNOSE), and returns a LSQVCE.RESULTS object that carries the
%   solution, the per-observation diagnostics, and the accumulator itself. The
%   six-step dance (accumulator -> addSystem -> solve -> diagnoseAccum -> diagnose
%   -> report) becomes one line; the lower layers stay reachable on R for
%   experts (R.acc, R.diag).
%
%   This performs ONE GLS solve with whatever covariance the loaders supply.
%   The iterative separable refinement (temporal whitening + VCE reweighting)
%   is layered on top by LSQVCE.ENGINE.ADJUSTVCEJOINT; adjustBlock is its base.
%
%   INPUT
%     systems  1xM   struct|cell  The systems to combine. Two accepted forms:
%              (a) struct array with fields
%                    .sysID   1x1        system identifier (numeric/char/string)
%                    .nl      1x1        number of LOCAL parameters
%                    .nc      1x1        number of COMMON parameters (shared by
%                                        every system; declared here beside .nl)
%                    .loader  1x1        a Loader from LSQVCE.ENGINE.MAKELOADER, or a raw
%                                        4-arg loader closure
%                                        @(day,obsType,arcID,opt)
%                    .days    1xD  (opt) per-system day list; else the DAYS arg
%              (b) cell array of Loader objects (LSQVCE.ENGINE.MAKELOADER), each of which
%                  already carries its own .sysID, .nl, .nc and (optionally) .days.
%                              Required.
%     days     1xD   numeric   Default day list for systems that do not carry
%                              their own. Required (may be [] if every system
%                              specifies .days).
%
%   NOTE ON nc. There is no nc argument. The common-parameter count is DERIVED
%   from the systems (LSQVCE.ENGINE.RESOLVENC): from L.nc for a Loader, which is itself
%   the number of rows in its 'commonParams' table, or from the .nc field of a
%   raw spec. Stating it separately would let the same number be written twice
%   and disagree, which is exactly the bug this removes.
%
%     Name/Value options (WP-39: canonical UpperCamelCase names -- legacy
%     casings still match, MATLAB NV semantics; enum VALUES are
%     case-sensitive lowercase):
%       'Backend'    string 1x1  "qr" | "normal"                (default "qr")
%       'Reg'        string 1x1  "none"|"tikhonov"|"tsvd"|"lm"|"auto"
%                                                               (default "none")
%       'RegLambda'  double 1x1 >= 0  regularisation weight     (default 0)
%       'RegMatrix'  double p x nc    regularisation matrix     (default [])
%       'Parallel'   string 1x1  "serial"|"parfor"|"spmd"       (default "serial")
%       'DistributeDays' logical 1x1  parallelise over the DAYS of each
%                                system (LSQVCE.ENGINE.BLOCKACCUMULATOR 'DistributeDays':
%                                one raw partial per day, merged in day order;
%                                the axis that scales for a 2-satellite mission).
%                                Needs Parallel ~= "serial" and at least
%                                ParallelThreshold days.       (default false)
%       'ParallelThreshold' double 1x1  tasks below which no pool is used
%                                                              (default 8)
%       'Verbose'    logical 1x1  accumulator chatter           (default false)
%       'Prior'      [] | lsqvce.engine.CommonPrior  the prior on the common
%                                block (lsqvce.fit stacks data.priors)   (default [])
%
%   PLOTTING (R5 M5 5c): the post-solve plot hooks of v4 are retired; plot
%     from the returned Results -- R.plot("residuals" | "contribution" | ...).
%
%   OUTPUT
%     R        1x1   lsqvce.Results   Solution + diagnostics + accumulator.
%                              R.x, R.sigma0, R.report(), R.redundancy, R.acc, ...
%                              (see LSQVCE.RESULTS).
%
%   EXAMPLE
%     cov.pos = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',0.7);
%     cov.kbr = struct('type','diag','sigma',0.2);
%     Lc = lsqvce.engine.makeLoader(@(d,ot) buildDayC(d,ot), 'sysID','C', 'nl',4, ...
%              'obsTypes',{'pos','kbr'}, 'cov',cov, 'arcs',2);
%     R  = lsqvce.engine.adjustBlock({Lc}, 1:3);
%     R.report();          % levels 0-4
%     xhat = R.x;          % common-parameter estimate
%
%   ERRORS
%     lsqvce:adjustBlock:badSystems  systems is neither a struct array nor a
%                                    cell of loaders
%     lsqvce:adjustBlock:badCell     a cell entry is not a loader struct
%     lsqvce:adjustBlock:badStruct   a system struct is missing a required field
%     lsqvce:adjustBlock:badLoader   a loader's .loader is not a function handle
%     lsqvce:adjustBlock:noDays      no day list, and none derivable from the
%                                    loaders (set 'days' on either)
%
%   'Provenance'  struct  (R5 M2) the resolved user option struct recorded in
%                 Results.meta.options (fit passes opts.toStruct() + .kind);
%                 struct() = record the driver's own resolved options
%
%   See also LSQVCE.RESULTS, LSQVCE.ENGINE.MAKELOADER, LSQVCE.ENGINE.ADJUSTVCEJOINT,
%            LSQVCE.ENGINE.BLOCKACCUMULATOR, LSQVCE.ENGINE.VERDICT.SOURCE, LSQVCE.ENGINE.VERDICT.DIAGNOSE,
%            LSQVCE.REPORT, LSQVCE.RESULTS
%
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 4: the verdict comes
%   from lsqvce.engine.verdict.fromAccumulator (diagnoseAccum + diagnose
%   absorbed). With assistance from Claude (Anthropic), model Claude Opus
%   5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-35): Prior
%   documented; ShowSeconds / ExportDir, documented but never accepted,
%   removed. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-08-31 (CEST, UTC+02:00), R4 WP-62 (Level L4): 'DistributeDays'
%   and 'ParallelThreshold' exposed (pass-through to the accumulator; the
%   engine's own defaults). With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (a) (F-07): typed
%   `arguments` signature; option names are UpperCamelCase (Backend, Reg,
%   RegLambda, RegMatrix, Parallel, Verbose, ShowSeconds, ExportDir) with
%   case-insensitive NAME matching kept (MATLAB NV semantics, measured
%   round a1-2), and enum VALUES are case-sensitive lowercase
%   (mustBeMember; strcmpi used to accept any casing) -- release note in
%   HANDOVER §9 WP-39. With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Created 2026-07-16 12:20 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
%   Updated 2026-07-17 (CEST, UTC+02:00): added the postSolve plot-hook seam
%   ('showSeconds'/'exportDir'), entry-point parity with adjustVCE.
%   Renamed 2026-07-17 from RECOVER -> ADJUSTBLOCK (block-arrow adjustment;
%   pairs with LSQVCE.ENGINE.ADJUST and LSQVCE.ENGINE.ADJUSTVCEJOINT). The old name was REMOVED on
%   2026-07-18 along with every other alias (pre-release; one name per thing).
%   Claude (Anthropic), model Claude Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
arguments
    systems
    days
    opt.Backend  (1,1) string {mustBeMember(opt.Backend, ["normal","qr"])} = "qr"
    opt.Reg      (1,1) string {mustBeMember(opt.Reg, ["none","tikhonov","tsvd","lm","auto"])} = "none"
    opt.RegLambda (1,1) double {mustBeNonnegative} = 0
    opt.RegMatrix double = []
    opt.Parallel (1,1) string {mustBeMember(opt.Parallel, ["serial","parfor","spmd"])} = "serial"
    opt.DistributeDays (1,1) logical = false
    opt.ParallelThreshold (1,1) double {mustBeInteger, mustBeNonnegative} = 8
    opt.Verbose  (1,1) logical = false
    opt.Prior = []                       % R5 M5 5b: lsqvce.engine.CommonPrior on the common
                                         % block (data.priors, stacked by lsqvce.fit), or []
    opt.Diagnostics cell = {}          % R5 M5 5b': Name-Value pairs for lsqvce.engine.verdict.diagnose
    opt.Provenance struct = struct()     % R5 M2: the RESOLVED user options to record
                                         % in Results.meta instead of this driver's own
                                         % option struct (lsqvce.fit passes opts.toStruct()
                                         % with .kind); empty = record o (engine-tier use)
end
    narginchk(2, inf);
    % WP-39 stage (a): typed `arguments` signature; names UpperCamelCase
    % (was backend/reg/regLambda/regMatrix/parallel/verbose/showSeconds/
    % exportDir -- release note in HANDOVER §9 WP-39; NV NAME matching
    % stays case-insensitive per MATLAB semantics). Enum values are now
    % case-sensitive lowercase (mustBeMember), where strcmpi used to accept
    % any casing.
    o = struct('backend', char(opt.Backend), 'reg', char(opt.Reg), ...
               'regLambda', opt.RegLambda, 'regMatrix', opt.RegMatrix, ...
               'parallel', char(opt.Parallel), 'verbose', opt.Verbose, ...
               'distributeDays', opt.DistributeDays, ...
               'parallelThreshold', opt.ParallelThreshold);

    [specs, ~] = normalizeSystems(systems, days);

    % ORDER MATTERS: normalizeSystems validates the SHAPE of `systems` first,
    % so a malformed entry reports the specific lsqvce:adjustBlock:badCell /
    % :badStruct that names what is wrong with it. Deriving nc first made every
    % such mistake surface as a confusing lsqvce:resolveNc:missing ("no system
    % declares a common-parameter count") when the real problem was that the
    % entry was not a system at all. nc is not needed until the accumulator is
    % constructed below, so nothing is lost by deriving it here.
    nc = lsqvce.engine.resolveNc(systems);

    acc = lsqvce.engine.BlockAccumulator(nc, 'Backend',o.backend, 'KeepLocal',true, ...
              'Verbose',o.verbose, 'Reg',o.reg, 'RegLambda',o.regLambda, ...
              'RegMatrix',o.regMatrix, 'Parallel',o.parallel, ...
              'DistributeDays',o.distributeDays, ...
              'ParallelThreshold',o.parallelThreshold);
    acc.setPrior(opt.Prior);
    acc.addSystems(specs);
    sol = acc.solve();

    [D, dsol, batchFcn, nB] = lsqvce.engine.verdict.fromAccumulator(acc, sol, opt.Diagnostics{:});

    % WP-44: run provenance through the ONE composer -- adjustBlock resolves
    % no RNG seed (empty struct), and o IS its resolved option struct.
    prov = o;                                   % engine-tier call: its own options
    if isfield(opt, 'Provenance') && ~isempty(fieldnames(opt.Provenance))
        prov = opt.Provenance;                  % fit: the user's resolved Options
    end
    R = lsqvce.Results(sol, D, acc, dsol, batchFcn, nB, specs, ...
                       lsqvce.engine.runMeta(struct(), prov));
    % R5 M5 5c: the post-solve plot hooks are retired -- plotting is
    % lsqvce.Results.plot(kind) on the returned Results.
end

%==========================================================================