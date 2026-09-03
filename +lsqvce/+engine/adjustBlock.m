function R = adjustBlock(systems, days, opt)
%LSQVCE.ENGINE.ADJUSTBLOCK  One-call block-arrow adjustment: accumulate, solve, judge.
%   R = lsqvce.engine.adjustBlock(systems, days)
%   R = lsqvce.engine.adjustBlock(systems, days, Name, Value, ...)
%
%   The single front door over the streaming pipeline. It builds a
%   LSQVCE.ENGINE.BLOCKACCUMULATOR, adds every system, solves the reduced common
%   system, bridges the result into the diagnostics (LSQVCE.ENGINE.DIAGNOSEACCUM /
%   LSQVCE.ENGINE.DIAGNOSE), and returns a LSQVCE.RESULTS object that carries the
%   solution, the per-observation diagnostics, and the accumulator itself. The
%   six-step dance (accumulator -> addSystem -> solve -> diagnoseAccum -> diagnose
%   -> report) becomes one line; the lower layers stay reachable on R for
%   experts (R.acc, R.diag).
%
%   This performs ONE GLS solve with whatever covariance the loaders supply.
%   The iterative separable refinement (temporal whitening + VCE reweighting)
%   is layered on top by LSQVCE.ENGINE.ADJUSTVCE; adjustBlock is its base.
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
%       'ShowSeconds' double 1x1  seconds to display each plot-hook figure
%                            (default 0 = export only; < 0 = show and
%                            LEAVE OPEN -- the v1 contract, unchanged).
%                            Engine-side knob, mirrors LSQVCE.ENGINE.ADJUSTVCE.
%       'ExportDir'  string 1x1  folder for plot-hook PNGs ("" = do not export).
%
%   PLOT HOOKS (single-solve, postSolve only)
%     If a loader was built with LSQVCE.ENGINE.MAKELOADER's 'plotPostSolve', adjustBlock
%     fires that hook once per (system,obsType) stream AFTER the solve, exactly
%     as the LSQVCE.ENGINE.ADJUSTVCE postSolve seam does, keyed 'satKey|obsType'. The
%     hook receives ONE ctx struct (.stage='postSolve', .stream, .satKey,
%     .obsType, .sysID, .pass=1, .R (the Results), .x, .exportDir, .showSeconds)
%     and returns a figure; the engine owns the timed show + PNG export via
%     LSQVCE.UTIL.SHOWEXPORT. A throwing or figure-less hook is logged
%     ('lsqvce:adjustBlock:plotHookFailed') and skipped -- plots NEVER break the
%     solve. There is NO 'plotPostModel' seam here: adjustBlock does ONE GLS solve
%     with the loaders' covariance and estimates no stochastic model, so
%     postModel has nothing to fire on; use LSQVCE.ENGINE.ADJUSTVCE for that.
%     Raw 4-arg closures (not Loader objects) carry no hooks and are skipped.
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
%   See also LSQVCE.RESULTS, LSQVCE.ENGINE.MAKELOADER, LSQVCE.ENGINE.ADJUSTVCE,
%            LSQVCE.ENGINE.BLOCKACCUMULATOR, LSQVCE.ENGINE.DIAGNOSEACCUM, LSQVCE.ENGINE.DIAGNOSE,
%            LSQVCE.REPORT, LSQVCE.UTIL.SHOWEXPORT
%
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
%   pairs with LSQVCE.ENGINE.ADJUST and LSQVCE.ENGINE.ADJUSTVCE). The old name was REMOVED on
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
    opt.ShowSeconds (1,1) double = 0
    opt.ExportDir (1,1) string = ""
    opt.Prior = []                       % R5 M5 5b: lsqvce.engine.CommonPrior on the common
                                         % block (data.priors, stacked by lsqvce.fit), or []
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
               'parallelThreshold', opt.ParallelThreshold, ...
               'showSeconds', opt.ShowSeconds, 'exportDir', char(opt.ExportDir));

    [specs, loaders] = normalizeSystems(systems, days);

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

    [dsol, batchFcn, nB] = lsqvce.engine.diagnoseAccum(acc, sol);
    D = lsqvce.engine.diagnose(dsol, batchFcn, 'nBatch', nB);

    % WP-44: run provenance through the ONE composer -- adjustBlock resolves
    % no RNG seed (empty struct), and o IS its resolved option struct.
    prov = o;                                   % engine-tier call: its own options
    if isfield(opt, 'Provenance') && ~isempty(fieldnames(opt.Provenance))
        prov = opt.Provenance;                  % fit: the user's resolved Options
    end
    R = lsqvce.Results(sol, D, acc, dsol, batchFcn, nB, specs, ...
                       lsqvce.engine.runMeta(struct(), prov));

    % ---- postSolve plot hooks (single-solve; no postModel seam) -----------
    firePostSolvePlots(loaders, R, double(o.showSeconds), char(o.exportDir));
end

%==========================================================================
function firePostSolvePlots(loaders, R, showSeconds, exportDir)
%FIREPOSTSOLVEPLOTS  Call each Loader's plotPostSolve hook once per
%   (system,obsType) stream, then show/export via lsqvce.util.showExport.
%   Mirrors lsqvce.engine.adjustVCE.firePlots(postSolve): keyed 'satKey|obsType',
%   ONE ctx struct in, a figure out, and every call wrapped so a throwing or
%   figure-less hook can never break the solve. Raw closures (no .plotPostSolve)
%   and loaders without hooks are simply skipped.
    for i = 1:numel(loaders)
        L = loaders{i};
        if ~isstruct(L) || ~isfield(L, 'plotPostSolve') || isempty(L.plotPostSolve)
            continue
        end
        hookMap = L.plotPostSolve;
        if ~isa(hookMap, 'containers.Map') || hookMap.Count == 0, continue, end
        satKey  = lsqvce.engine.satKeyOf(L.sysID);
        obsTypes = L.obsTypes;
        for k = 1:numel(obsTypes)
            ot = obsTypes{k};
            if ~isKey(hookMap, ot), continue, end
            key = sprintf('%s|%s', satKey, ot);
            ctx = struct('stage','postSolve', 'stream',key, 'satKey',satKey, ...
                         'obsType',ot, 'sysID',L.sysID, 'pass',1, ...
                         'R',R, 'x',R.x, 'exportDir',exportDir, ...
                         'showSeconds',showSeconds);
            hookFcn = hookMap(ot);
            try
                figH = hookFcn(ctx);
            catch ME
                lsqvce.util.log.warn('lsqvce:adjustBlock:plotHookFailed', ...
                    'postSolve plot hook for %s threw (%s); skipping.', ...
                    key, ME.identifier);
                continue
            end
            name = sprintf('postSolve_%s_pass1', matlab.lang.makeValidName(key));
            lsqvce.util.showExport(figH, name, showSeconds, exportDir);
        end
    end
end

%==========================================================================