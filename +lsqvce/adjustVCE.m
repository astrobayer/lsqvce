classdef adjustVCE < handle
%LSQVCE.ADJUSTVCE  Iterative separable-GLS engine (stochastic model + variance components) over adjustBlock.
%
%   P = lsqvce.adjustVCE(systems, days)
%   P = lsqvce.adjustVCE(systems, days, Name, Value, ...)
%
%   The general separable generalised-least-squares engine. It layers the
%   two-unknown covariance loop on top of LSQVCE.ADJUSTBLOCK: a single GLS solve is
%   adjustBlock's job; discovering the covariance that solve should use is this
%   class's job. Build it fluently, then call RUN.
%
%       R = lsqvce.adjustVCE({Lc, Ld}, 1:30) ...
%              .groupComponentsBy('satObs') ...
%              .stochasticModel('MaxOrder',60, 'AcvTol',1e-3) ...
%              .varianceComponents('MaxIter',8, 'Tol',1e-3) ...
%              .run();
%       R.report();                 % the usual level-0..4 verdict
%       R.history                   % the sigma / AR-order trajectory
%
%   THE MODEL (per group g)
%   -----------------------
%       C_g = sigma_g^2 * D_g (R_{t,g} kron I_d) D_g
%
%     D_g        FIXED per-epoch scaling from the loader (e.g. the PPP epoch
%                covariance). Never estimated here.
%     R_{t,g}    TEMPORAL correlation, the SAME on every one of the d axes.
%                ESTIMATED -- this is the stochastic model.
%     sigma_g^2  one scalar variance component per group. ESTIMATED by VCE --
%                this is the variance components.
%     kron I_d   no cross-axis correlation (what makes the whitener separable
%                and cheap; see LSQVCE.SEPARABLEWHITEN).
%
%   THE TWO-LAYER LOOP
%   ------------------
%       initial weights (sigma_g^2 = 1, R_{t,g} = white)
%              |
%         (*)  v  SOLVE ......... lsqvce.adjustBlock (block-arrow GLS)
%              |
%              |  STOCHASTIC MODEL (within group): re-estimate R_{t,g}
%              |    collectResiduals -> residualAcv -> selectArAcv
%              |    -> new AR polynomial per (sat,obsType) stream
%              |
%              |  VARIANCE COMPONENTS (between groups): VCE (Helmert), one sigma_g^2 per group
%              |    lsqvce.vce: sigma_g^2 = ||w_g||^2 / r_g  (redundancy pass)
%              |
%              |  REWEIGHT: push sigma_g and R_{t,g} back into the loaders as a
%              |    per-obsType covariance override -> next pass whitens with them
%              |
%              +-- repeat from (*) until CONVERGED
%
%   CONVERGENCE (both, per group)
%     WHITE     every stream's whitened residuals are serially uncorrelated
%               (selectArAcv reports order 0 / 'white', or the achieved AR order
%               has stopped changing).
%     BALANCED  every sigma_g^2 ~= 1 -- VCE no longer wants to move weight.
%
%   WIRING, NOT NEW NUMERICS
%     Every numerical step already exists and is tested: LSQVCE.ADJUSTBLOCK,
%     LSQVCE.SEPARABLEWHITEN, LSQVCE.SELECTARACV, LSQVCE.RESIDUALACV,
%     LSQVCE.COLLECTRESIDUALS, LSQVCE.REDUNDANCY, LSQVCE.VCE, LSQVCE.DIAGNOSEACCUM.
%     This class only alternates and connects them, and threads the estimated
%     covariance back through the loader's opt.covOverride hook.
%
%   INPUT
%     systems  1xM   cell|struct  Loader objects from LSQVCE.MAKELOADER (cell),
%                              or the struct-array spec that LSQVCE.ADJUSTBLOCK
%                              accepts. Required. Each obsType stream becomes a
%                              candidate VCE group (see groupComponentsBy).
%                              The number of COMMON parameters is DERIVED from
%                              these (LSQVCE.RESOLVENC) -- from a Loader's
%                              'commonParams' rows, or a raw spec's .nc field --
%                              so it is stated exactly once, where the system is
%                              defined, and cannot disagree with itself.
%     days     1xD   numeric   Default day list for systems without their own.
%                              Required (may be [] if every loader carries days).
%
%     Name/Value options ([] or omitted -> defaults):
%       'GroupBy'   char|fh  grouping granularity for the variance-component step (same name
%                            and meaning as in lsqvce.vce/redundancy/outliers):
%                            'satObs' (default) one component per
%                            (satellite,obsType); 'obsType'; 'system'; 'arc';
%                            or @(sysID,obsType,arc,day)->key. Also settable via
%                            groupComponentsBy().
%       'MaxOuter'  1x1  outer stochastic model <-> variance components iterations             (default 6)
%       'Backend'   char 'qr' | 'normal'                          (default 'qr')
%       'Reg'       char passed to adjustBlock                        (default 'none')
%       'Parallel'  char 'serial'|'parfor'|'spmd'                 (default 'serial')
%       'verbose'   1x1  logical | 0..2 message level             (default 1)
%                        0 SILENT     nothing
%                        1 PROGRESS   outer passes, stage start/done + timing,
%                                     the convergence verdict
%                        2 DETAIL     everything else: per-day accumulation and
%                                     residuals, per-stream AR order and
%                                     whiteness, per-group shares and sigmas,
%                                     system conditioning, inner VCE iterations
%                        A LOGICAL still works and keeps its old meaning
%                        (false -> 0, true -> 1). See LSQVCE.UTIL.VERBOSITY.
%       'LogFid'    1x1  file identifier for ALL toolbox output ("standard
%                        output"). Default: leaves the global log state as-is
%                        (initially 1 = command window). Set to an fopen fid to
%                        capture status text AND warnings to a file. Sets the
%                        toolbox-wide state; see LSQVCE.UTIL.LOG.
%       'LogFlag'   1x1  logical: also show output on screen (default: leave
%                        as-is, initially true). false + a file LogFid => quiet
%                        run, everything captured to the file only.
%       'ShowSeconds' 1x1  seconds to display each hook figure (default 0 =
%                        export only, no pause). < 0 = show and LEAVE OPEN.
%       'ExportDir' char  folder for the hook PNGs ('' = do not export).
%
%   PLOT HOOKS (per system + obsType)
%     Per-(system,obsType) plots are driven by function handles you attach to
%     the LOADERS, via LSQVCE.MAKELOADER's 'plotPostSolve' and 'plotPostModel'
%     (a single @(ctx)->figHandle for all the loader's obsTypes, or a
%     containers.Map/struct keyed by obsType). The engine fires them:
%       * postSolve -- after each solve, BEFORE the residual pass;
%       * postModel -- ONCE after the stochastic model is estimated.
%     ctx is a FIXED, DOCUMENTED struct (ctx.schema, currently 1); fields are
%     only ADDED across versions, never renamed, so a hook keeps working. It
%     carries, on BOTH stages: ctx.x (common ESTIMATES), ctx.sd (their STANDARD
%     DEVIATIONS), ctx.common (common parametrisation labels), ctx.procOpt (the
%     loader's processingOptions), ctx.model (the fitted noise model), ctx.sigmaG
%     (the stream's VCE weight), plus ctx.R/stream/satKey/obsType/sysID/pass. On
%     postModel it ALSO carries the estimated ctx.psd (+ctx.psdFreq), ctx.acv and
%     ctx.arOrder. See PLOTCONTEXT for the full field list. Local parameters are
%     NOT in ctx (they are summarised to the log after each solve instead).
%     Your handle only BUILDS and RETURNS a figure (or returns nothing to just
%     record); the engine shows it for 'ShowSeconds' and exports
%     <ExportDir>/<stage>_<sys>_<obsType>_pass<k>.png.
%     A hook that errors is logged and skipped -- plots never break
%     the solve. Set the two engine knobs via the constructor or P.plots(...).
%
%   FLUENT CONFIGURATORS (each returns P, so they chain)
%     P.groupComponentsBy(key)          which observations share a variance component
%     P.stochasticModel('MaxOrder',p,'AcvTol',t)   the temporal noise model
%     P.varianceComponents('MaxIter',n,'Tol',t)    the Helmert variance components
%     P.iteration('MaxOuter',k,'Tol',t) the outer convergence gate
%     P.plots('ShowSeconds',s,'ExportDir',d)  configure the plot-hook show/export
%
%   RUN
%     R = P.run()   ->  1x1 lsqvce.Results, with the iteration trajectory
%                       attached (R.history, R.converged). R is the SAME object
%                       LSQVCE.ADJUSTBLOCK returns, so R.report / R.residuals /
%                       R.acc all work unchanged.
%
%   OUTPUT of run() -- R.history (struct)
%     .sigma      nOuter x nGroup   the variance-component trajectory
%     .groups     1xnGroup cellstr  the VCE group keys, aligned with .sigma
%     .contribution nOuter x nGroup Level-1 relative REDUNDANCY share per group
%                 per pass, share_g = t_g/sum_g t_g with t_g = tr(N^-1 N_g^cc)
%                 (the cc-attribution trace; NOT the redundancy r_g -- WP-32,
%                 F-01: r_g = n_g - local_g - t_g. The Level-2 diagonal in
%                 .contribDiag is the ABSORBED-information view instead), rows
%                 sum to 1. WEIGHT-DEPENDENT: it moves across passes, the final
%                 row is authoritative. NaN entries mark a group with no finite
%                 redundancy. Plot it over passes to see the share settle.
%     .arOrder    nOuter x nStream  selected AR order per stream per outer pass
%     .streams    1xnStream cell    the (sat,obsType) stream keys
%     .white      nOuter x nStream  logical: stream declared white this pass
%     .whitePval  nOuter x nStream  Ljung-Box whiteness p-value per stream per
%                 pass (NaN where trivially white); p > WhiteAlpha => white.
%     .converged  1x1 logical       both criteria met
%     .nOuter     1x1               outer passes actually run
%
%   ENTRY POINT / WHAT YOU ADAPT  --  see the block at the bottom of this help,
%   under "USER-SPECIFIC ADAPTATION".
%
%   EXAMPLE (synthetic, end to end)
%     [cfg, L] = demo_grace_dataLoader();          % GRACE-FO synthetic (examples/)
%     R = lsqvce.adjustVCE(L, cfg.days) ...
%            .groupComponentsBy('obsType').run();
%     R.report();
%     disp(R.history.sigma(end,:));   % ~ 1 per group at truth
%     % or just run the driver:  demo_grace
%
%   ERRORS
%     lsqvce:adjustVCE:badSystems     systems is not a loader cell/struct array
%     lsqvce:adjustVCE:badComponents  groupComponentsBy got a key that is
%                                     neither a known name nor a function handle
%     lsqvce:adjustVCE:validateBadResults  validate() was given something that
%                                     is not a Results from this problem
%
%     Errors raised DEEPER are not re-wrapped: a malformed loader surfaces as
%     lsqvce:makeLoader:*, and a bad covariance spec as lsqvce:covModel:*. That
%     is deliberate -- the identifier names the layer that can actually explain
%     the problem.
%
%   See also LSQVCE.ADJUSTBLOCK, LSQVCE.VCE, LSQVCE.REDUNDANCY,
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (c) c2 (F-07):
%   typed `arguments` signatures, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-41 (docs only, F-12 claims
%   half + F-14): the 'joint' mode claim corrected (pooled spectrum +
%   stream scales; NO per-arc factors -- that is ADJUSTVCEJOINT); the
%   Level-1 history field documented as the cc-trace share it is (not a
%   redundancy); v1-history prose retargeted to HISTORY_v1.md. With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-40 (fixes F-12, gate half):
%   'joint' convergence additionally requires the LAST joint pass to meet
%   the outer Tol (hist.jointMaxDev(end) <= Tol) -- shapeStable alone
%   degenerated from pass 2. With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   Updated 2026-08-18 (CEST, UTC+02:00), WP-33 (fixes F-03, CRITICAL):
%   sweepStreamLags subtracts the per-system local fit B*y and absorbs the
%   local block exactly into the per-lag projector (shared jointArcParts +
%   redundancyLags 'LocalBbar'/'LocalLl'; axis pooling threads the sliced
%   local rows symmetrically). nl = 0 paths are verbatim unchanged. With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%            LSQVCE.SELECTARACV, LSQVCE.RESIDUALACV,
%            LSQVCE.SEPARABLEWHITEN, LSQVCE.MAKELOADER, LSQVCE.RESULTS
%
% ==========================================================================
%   USER-SPECIFIC ADAPTATION (read before running on your own data)
% ==========================================================================
%   1. LOADER: write buildDay(day,obsType) -> .A .B .l .t (.arcIdx) and wrap it
%      with LSQVCE.MAKELOADER, one Loader per satellite. Give each obsType an
%      initial 'cov' spec (type 'diag' for positions with a per-epoch sigma,
%      type 'arma'/'diag' for range-rate). The loop refines .ar and .sigma; it
%      does NOT invent D -- your buildDay/cov must supply the fixed per-epoch
%      scaling.
%   2. nAxes PER OBSTYPE: 3-D position streams have nAxes=3 (x,y,z share one
%      R_t); range-rate is scalar, nAxes=1. Declare it via nAxesOf() below
%      (default: 'pos'->3, everything else->1). Override with P.componentsPerEpoch(map).
%   3. VCE GRANULARITY: default one sigma per (satellite,obsType). If a group
%      has too little redundancy (VCE warns), coarsen with
%      P.groupComponentsBy('obsType').
%   4. COUPLED OBSERVATIONS (e.g. inter-satellite range-rate). The block-arrow
%      model has ONE local block per system, so an observation that depends on
%      TWO satellites' arc states (KBR/range-rate between C and D) cannot be a
%      third system. Merge the two satellites into ONE system with a stacked
%      local block [y_C | y_D] (nl = nl_C + nl_D) carrying all obsTypes; each
%      obsType's buildDay fills the local sub-columns it touches (posC -> first
%      nl_C, posD -> next nl_D, rr -> both). Give the three streams DISTINCT
%      obsType names ('posC','posD','rr') and group VCE by 'obsType'. Give the
%      merged system a NUMERIC sysID (or a single char) -- a multi-char string
%      like 'CD' breaks the accumulator's [store.sysID] indexing. See
%      demo_grace_dataLoader (in examples/) for a worked template. (If instead
%      the range-rate depends only on the gravity field, add it as a third
%      system with nl = 0.)
%   5. PRIORS: a decent starting sigma and AR order shorten the loop but are not
%      required; sigma_g^2 = 1 and white start both converge.
% ==========================================================================
%
%   Created 2026-07-16 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
%   Renamed 2026-07-17 from SeparableProblem -> adjustVCE (pairs with the
%   in-memory shortcut LSQVCE.ADJUST: adjust is "just solve it", adjustVCE is
%   "solve AND discover the stochastic model"). The old name was REMOVED on
%   2026-07-18 along with every other alias: this is pre-release, and two ways
%   to say one thing is the ambiguity that produced the nc duplication.
%   Claude (Anthropic), model Claude Opus 4.8.

    properties (SetAccess = immutable)
        nc              % 1x1  common parameter count
        systems         % 1xM  cell/struct of loaders (as given)
        days            % 1xD  default day list
    end
    properties (SetAccess = private)
        vceGroups   = 'satObs';   % variance-component grouping key
        vceMode     = 'separate'; % 'separate' (alternate stochastic model / variance components) | 'joint'
        jointOpt              % struct: NumLags, Probes, Seed, Norm, MeanType, C0
        shapeOpt              % struct: MaxOrder, AcvTol, MaxLag
        scaleOpt              % struct: MaxIter, Tol, MinRedundancy
        outerOpt              % struct: MaxOuter, Tol
        axesMap               % containers.Map obsType -> nAxes
        blockOpt            % struct: Backend, Reg, Parallel
        plotOpt               % struct: ShowSeconds, ExportDir (engine-side)
        verbose     = 1     % 0 silent | 1 progress | 2 detail
                            % (a logical is accepted and resolved: false->0,
                            %  true->1; see LSQVCE.UTIL.VERBOSITY)
    end

    methods
        % ---- construction ------------------------------------------------
        function obj = adjustVCE(systems, days, opt)
        arguments
            systems
            days
            opt.GroupBy = 'satObs'
            opt.VceMode (1,1) string {mustBeMember(opt.VceMode, ["separate","joint"])} = "separate"
            opt.MaxOuter (1,1) double {mustBeInteger, mustBePositive} = 6
            opt.Backend (1,1) string {mustBeMember(opt.Backend, ["normal","qr"])} = "qr"
            opt.Reg (1,1) string {mustBeMember(opt.Reg, ["none","tikhonov","tsvd","lm","auto"])} = "none"
            opt.Parallel (1,1) string {mustBeMember(opt.Parallel, ["serial","parfor","spmd"])} = "serial"
            opt.Verbose = true             % historical logical OR level 0..3
            opt.LogFid double {mustBeScalarOrEmpty} = []
            opt.LogFlag logical = logical([])
            opt.ShowSeconds (1,1) double = 0
            opt.ExportDir (1,1) string = ""
        end
            % WP-39 c2: the old narginchk(2, inf) is dead -- the arguments
            % block enforces the arity first (missing positionals surface
            % MATLAB:minrhs) -- and is deleted (rule 8).
            obj.nc      = lsqvce.resolveNc(systems);
            obj.systems = systems;
            obj.days    = days;

            % WP-39 c2: typed `arguments` signature (canonical UpperCamel --
            % was lowercase; legacy casings keep matching, MATLAB NV semantics;
            % enum VALUES case-sensitive lowercase). The o remap keeps the ctor
            % body untouched; Verbose keeps the logical-or-level contract.
            o = struct('groupBy', opt.GroupBy, 'vceMode', char(opt.VceMode), ...
                       'maxOuter', opt.MaxOuter, 'backend', char(opt.Backend), ...
                       'reg', char(opt.Reg), 'parallel', char(opt.Parallel), ...
                       'verbose', opt.Verbose, 'logFid', opt.LogFid, ...
                       'logFlag', opt.LogFlag, 'showSeconds', opt.ShowSeconds, ...
                       'exportDir', char(opt.ExportDir));
            % Route all toolbox output (status + warnings) to logfid, on screen
            % only if logflag. These set the toolbox-wide log state (see
            % lsqvce.util.log); omit them to leave the current state untouched.
            if ~isempty(o.logFid) || ~isempty(o.logFlag)
                args = {};
                if ~isempty(o.logFid),  args = [args, {'fid',    o.logFid}];  end
                if ~isempty(o.logFlag), args = [args, {'screen', o.logFlag}]; end
                lsqvce.util.log.set(args{:});
            end
            % (o built above from the typed opt -- WP-39 c2)

            obj.vceGroups  = o.groupBy;
            obj.vceMode    = lower(char(o.vceMode));
            obj.jointOpt   = struct('NumLags',60, 'Probes',100, 'Seed',[], ...
                                    'Norm','ortho', 'MeanType','robust', 'C0',[], ...
                                    'PsdMode','toeplitz', 'Smooth',0);
            % Resolve ONCE here, so obj.verbose is always a numeric level
            % and every downstream check is a level comparison, never a
            % logical test that would silence detail at higher levels.
            obj.verbose    = lsqvce.util.verbosity.resolve(o.verbose);
            obj.blockOpt = struct('Backend',o.backend, 'Reg',o.reg, ...
                                    'Parallel',o.parallel);
            obj.shapeOpt   = struct('MaxOrder',60, 'AcvTol',1e-3, 'MaxLag',[], ...
                                    'WhiteAlpha',0.01);
            obj.scaleOpt   = struct('MaxIter',10, 'Tol',5e-3, 'MinRedundancy',50);
            obj.outerOpt   = struct('MaxOuter',o.maxOuter, 'Tol',5e-3);
            obj.plotOpt    = struct('ShowSeconds', double(o.showSeconds), ...
                                    'ExportDir', char(o.exportDir));
            obj.axesMap    = containers.Map('KeyType','char','ValueType','double');
            obj.axesMap('pos') = 3;      % 3-D positions share one temporal law
        end

        % ---- fluent configurators ---------------------------------------
        function obj = groupComponentsBy(obj, key)
        %GROUPCOMPONENTSBY  Which observations share one variance component.
        %   'satObs' (default) | 'obsType' | 'system' | 'arc' | a function handle
        %   @(sysID,obsType,arc,day) -> group key. This is the grouping of the
        %   Helmert variance-component estimation: every observation mapping to
        %   the same key gets one common variance factor.
            obj.vceGroups = key;
        end
        function obj = stochasticModel(obj, varargin)
        %STOCHASTICMODEL  Configure the TEMPORAL noise model (the covariance shape).
        %   'MaxOrder' (default 60)  highest AR order considered
        %   'AcvTol'   (1e-3)        autocovariance fit tolerance
        %   'MaxLag'   ([] = auto)   lags used to estimate the ACV
        %   'WhiteAlpha' (0.01)      Ljung-Box significance for the whiteness gate
        %   This sets the shape of the noise (how it is correlated in time), as
        %   opposed to VARIANCECOMPONENTS, which sets the relative WEIGHT of each
        %   observation group.
            obj.shapeOpt = mergeNV(obj.shapeOpt, varargin);
        end
        function obj = varianceComponents(obj, varargin)
        %VARIANCECOMPONENTS  Configure the Helmert variance-component estimation.
        %   'MaxIter'       (default 10)  inner VCE iterations
        %   'Tol'           (5e-3)        convergence on max|sigma-1|
        %   'MinRedundancy' (50)          below this a group cannot be estimated
        %   This sets the relative WEIGHT between observation groups, as opposed
        %   to STOCHASTICMODEL, which sets the temporal correlation shape.
            obj.scaleOpt = mergeNV(obj.scaleOpt, varargin);
        end
        function obj = iteration(obj, varargin)
        %ITERATION  Configure the outer loop over the stochastic model and the
        %   variance components.
        %   'MaxOuter' (default 6)  outer passes
        %   'Tol'      (5e-3)       convergence on max|sigma-1| and AR-order drift
            obj.outerOpt = mergeNV(obj.outerOpt, varargin);
        end
        function obj = joint(obj, varargin)
        %JOINT  Switch to the JOINT stochastic model and configure it.
        %   Estimates the temporal temporal law and the variance factors TOGETHER in one Helmert
        %   VCE per outer pass (Ellmer 2018 §6.5.7), instead of alternating them.
        %   Calling this sets vceMode to 'joint'; no separate flag is needed.
        %
        %   'numLags'  ('auto')       number of frequencies in the noise model.
        %                             'auto' sizes it from the data as
        %                             floor(mMax/2)+1 for the longest arc mMax,
        %                             which is the MINIMUM for the PSD-derived
        %                             covariance to be full rank: a PSD on nf
        %                             frequencies spans only 2*(nf-1) lags, so a
        %                             smaller nf gives a RANK DEFICIENT (not
        %                             merely ill-conditioned) covariance. Pass a
        %                             number to fix it, e.g. 60 for Ellmer's Nmax
        %                             -- but then keep batches under 2*(nf-1) rows.
        %   'probes'   (100)          Hutchinson probes Z
        %   'seed'     ([])           RNG seed for reproducible probes
        %   'norm'     ('ortho')      DCT normalisation, 'ortho' | 'ellmer'
        %   'meanType' ('robust')     arc-factor mean, 'robust' | 'plain'
        %   'smooth'   (0)            half-width of a moving average over the
        %                             estimated PSD amplitudes; trades spectral
        %                             resolution for variance (see LSQVCE.PSDVCE).
        %                             0 = off. Try 1-3 if the spectrum is noisy.
        %   'psdMode'  ('toeplitz')   how the estimated PSD becomes a covariance:
        %                             'toeplitz' builds the Toeplitz Sigma from
        %                             the PSD DIRECTLY (Ellmer eqs. 6.5.6/6.5.17)
        %                             -- no AR fit, so no stationarity condition
        %                             to violate. Cost O(b^3) per batch.
        %                             'ar'/'arma' fit a filter to the spectrum:
        %                             cheaper per batch (O(m*p)) but an estimated
        %                             PSD can yield a NON-STATIONARY AR whose
        %                             autocovariance does not exist, which is
        %                             rejected (see 'jointPsdRejected'). Prefer
        %                             'toeplitz' unless arcs are long enough that
        %                             the dense factorisation hurts.
        %   'c0'       ([])           starting covariance function; the default
        %                             is deliberately NOT a delta (that start is
        %                             degenerate -- see LSQVCE.PSDVCE).
        %
        %   WHEN TO USE WHICH
        %     'separate' (default) alternates stochastic model then variance components, each with its own
        %     gate, and stops on a Ljung-Box whiteness test. Robust, interpretable,
        %     and the right choice for general work.
        %     'joint' pools the per-lag sums over ALL streams and updates the
        %     spectrum and the per-stream scales together each outer pass --
        %     it does NOT estimate per-arc factors (WP-41 claims fix, F-12:
        %     the old sentence promised "every variance component at once").
        %     The FULL ten-step co-estimation with keyed per-arc factors is
        %     LSQVCE.ADJUSTVCEJOINT. The two can converge to DIFFERENT
        %     models, because the components are coupled -- scales and
        %     frequency amplitudes multiply the same Sigma. Use
        %     LSQVCE.ADJUSTVCEJOINT for Ellmer/ITSG fidelity.
            obj.vceMode  = 'joint';
            obj.jointOpt = mergeNV(obj.jointOpt, varargin);
        end
        function obj = separate(obj)
        %SEPARATE  Switch back to the default alternating stochastic model <-> variance components loop.
            obj.vceMode = 'separate';
        end
        function obj = componentsPerEpoch(obj, map)
        %COMPONENTSPEREPOCH  Declare how many components each observation type
        %   contributes per epoch: 3 for a kinematic position (x,y,z share one
        %   temporal law), 1 for a scalar such as range-rate. Accepts a struct
        %   obsType -> count, or a containers.Map.
            if isstruct(map)
                f = fieldnames(map);
                for k = 1:numel(f), obj.axesMap(f{k}) = map.(f{k}); end
            elseif isa(map,'containers.Map')
                obj.axesMap = map;
            else
                error('lsqvce:adjustVCE:badComponents', ...
                      ['componentsPerEpoch expects a struct or containers.Map ' ...
                       'mapping obsType to a component count.']);
            end
        end
        function obj = plots(obj, opt)
        %PLOTS  Configure the engine-side plot behaviour (the PER-ITEM plot
        %   HANDLES live on the loaders, via makeLoader's 'plotPostSolve' /
        %   'plotPostModel'; this only sets how they are shown and exported).
        %     'ShowSeconds' seconds to display each figure (default 0 = export
        %                   only, no pause). < 0 = show and LEAVE OPEN.
        %     'ExportDir'   folder for the PNGs ('' = do not export).
            % WP-39 c2: typed `arguments` signature; the defaults are the
            % CURRENT object state (arguments defaults may reference the
            % earlier-declared obj). The "< 0 = show and LEAVE OPEN"
            % contract is unchanged.
        arguments
            obj
            opt.ShowSeconds (1,1) double = obj.plotOpt.ShowSeconds
            opt.ExportDir (1,1) string = obj.plotOpt.ExportDir
        end
            obj.plotOpt.ShowSeconds = double(opt.ShowSeconds);
            obj.plotOpt.ExportDir   = char(opt.ExportDir);
        end

        % ---- the loop ----------------------------------------------------
        function R = run(obj)
        %RUN  Execute the stochastic model <-> variance components loop and return an lsqvce.Results.
            v = obj.verbose;

            % Discover the streams (sat,obsType) and their initial covariance.
            [streams, ~] = obj.discoverStreams();
            nStream = numel(streams);

            % The estimated covariance override, one field per (sat,obsType),
            % pushed into every loader via opt.covOverride on each pass.
            ovr = struct();        % ovr.(satKey).(obsType) = struct('ar',..,'sigma',..)

            hist.streams = streams;
            hist.arOrder = [];
            hist.white   = [];
            hist.whitePval = [];
            hist.sigma   = [];
            hist.groups  = {};
            hist.contribution = [];   % nOuter x nGroup relative trace share (Level-1)
            hist.jointUpdated = [];   % nOuter x 1: streams updated per joint pass
            hist.jointMaxDev  = [];   % nOuter x 1: max|alpha2-1| per joint pass

            R = [];  converged = false;  prevOrder = -ones(1,nStream);
            nOuter = obj.outerOpt.MaxOuter;
            vcFinal = struct('sigma',[], 'groups',{{}});   % last VCE result

            for it = 1:nOuter
                if lsqvce.util.verbosity.progress(v)
                    lsqvce.util.log.msg(['\n===== Iteration %d of %d ' ...
                        '=====================================\n'], it, nOuter);
                    lsqvce.util.log.msg(['      variance components -> adjustment -> ' ...
                        'stochastic model -> convergence\n']);
                end

                % ---------- VARIANCE COMPONENTS: VCE (Helmert) between groups ------------
                % Pin the group order so the sigma vector aligns 1:1 with the
                % keys covOverrideOpt maps against (no reliance on VCE's own
                % discovery order).
                if lsqvce.util.verbosity.progress(v)
                    lsqvce.util.log.msg('%sestimating between groups ...\n', ...
                        lsqvce.util.verbosity.banner('Variance components'));
                end
                tScale = tic;
                groupKeys = obj.vceGroupKeyList(streams, obj.vceGroupKeyFcn());
                groupKeys = unique(groupKeys, 'stable');
                % The accumulator VCE rebuilds on EVERY inner iteration, so its
                % per-day lines repeat ~4x per outer pass (~30x per run). Only a
                % diagnostic run wants that; below level 3 it is a flood that
                % buries the stage structure it sits inside.
                % THE ACCUMULATOR VCE REBUILDS STAYS SILENT AT EVERY LEVEL.
                % It is rebuilt on each inner Helmert iteration (~4 per outer
                % pass), and its per-day block is IDENTICAL each time apart from
                % the weights -- so printing it at level 3 repeated the same
                % wall four times and buried the thing level 3 actually adds.
                % What differs between iterations is the SIGMAS, and those are
                % reported by lsqvce.vce itself (it receives diagnostic(v) just
                % below). The one accumulation worth seeing is the SOLVE's, in
                % blockSolveWith, which runs once per pass and prints at level 2.
                accFcn = @(sig) obj.buildAccumulator(sig, ovr, 0);
                vc = lsqvce.vce(accFcn, 'GroupBy', obj.vceGroupKey(), ...
                          'Groups', groupKeys, ...
                          'MaxIter', obj.scaleOpt.MaxIter, ...
                          'Tol',     obj.scaleOpt.Tol, ...
                          'MinRedundancy', obj.scaleOpt.MinRedundancy, ...
                          ... % vce takes a LOGICAL and prints its INNER Helmert
                          ... % iterations -- the per-iteration sigmas. THIS is
                          ... % what level 3 adds over level 2; at levels 1-2 the
                          ... % 'Variance components' start/done pair is the
                          ... % summary of the same thing.
                          'Verbose', lsqvce.util.verbosity.detail(v));
                % fold the per-group sigma into the override (variance-factor part)
                ovr = obj.applyScale(ovr, streams, vc);
                vcFinal = vc;                        % remember for post-model plots
                hist.sigma(it,:) = vc.sigma(:).';    %#ok<AGROW>
                hist.groups      = vc.groups;

                % ---------- CONTRIB (Level 1): relative cc-trace share ------
                % share_g = r_g / sum_g r_g = tr(N^-1 N_g^cc) / (m-u), which
                % sums to 1 over groups (Sigma_g share_g = 1). FREE: vc.r is
                % already returned by lsqvce.vce (per-group redundancy, aligned
                % with vc.groups). It is WEIGHT-DEPENDENT, so it moves across
                % passes -- the FINAL pass is the authoritative share; a share
                % still drifting late is itself a signal the loop has not settled.
                share = contributionShare(vc.r);
                hist.contribution(it,:) = share(:).';   %#ok<AGROW>
                if lsqvce.util.verbosity.progress(v)
                    lsqvce.util.log.msg(['%sdone: %d group(s), max|sigma-1| = %.3e' ...
                        '   %.2f s\n'], lsqvce.util.verbosity.banner('Variance components'), numel(vc.sigma), ...
                        max(abs(vc.sigma(:)-1)), toc(tScale));
                end

                % ---------- SOLVE: block-arrow adjustment -----------------
                if lsqvce.util.verbosity.progress(v)
                    lsqvce.util.log.msg('%sblock-arrow, current weights ...\n', ...
                        lsqvce.util.verbosity.banner('Adjustment'));
                end
                tSolve = tic;
                % Build the Results for THIS weighting (variance factor applied). The
                % accumulator logs its own rcond / sigma0 / locals detail.
                R = obj.blockSolveWith(ovr);
                if lsqvce.util.verbosity.progress(v)
                    lsqvce.util.log.msg('%sdone   %.2f s\n', ...
                        lsqvce.util.verbosity.banner('Adjustment'), toc(tSolve));
                end

                % ---------- PLOT HOOKS: post-solve (before residuals) -----
                % Fire per (system,obsType) plot handles on the current estimate,
                % BEFORE the residual pass, exactly as requested. Never lets a
                % plotting error break the solve.
                obj.firePlots('postSolve', streams, R, ovr, vc, it);

                % ---------- STOCHASTIC MODEL: re-estimate R_{t,g} per stream ----------
                % In 'joint' mode the temporal shape is NOT fitted separately:
                % it was co-estimated with the scale in the step above, so the
                % per-stream AR fit is replaced by the spectral update. The
                % history columns still get filled (with the lag count and the
                % joint convergence flag) so downstream plots keep working.
                if strcmpi(obj.vceMode, 'joint')
                    if lsqvce.util.verbosity.progress(v)
                        lsqvce.util.log.msg('%sspectrum and variances together ...\n', ...
                            lsqvce.util.verbosity.banner('Joint estimation'));
                    end
                    tShape = tic;
                    [ovr, jinfo] = obj.applyJoint(ovr, streams, R, v);
                    ordr   = repmat(jinfo.nLag, 1, numel(streams));
                    white  = jinfo.converged & true(1, numel(streams));
                    pwhite = nan(1, numel(streams));
                    hist.arOrder(it,:)   = ordr(:).';    %#ok<AGROW>
                    hist.white(it,:)     = white(:).';   %#ok<AGROW>
                    hist.whitePval(it,:) = pwhite(:).';  %#ok<AGROW>
                    % Joint-mode diagnostics: how many streams actually received
                    % an updated noise model this pass, and the spectral
                    % convergence measure. A stream count BELOW numel(streams)
                    % means some stream contributed no usable arcs (warned above).
                    hist.jointUpdated(it,1) = jinfo.nUpdated;   %#ok<AGROW>
                    hist.jointMaxDev(it,1)  = jinfo.maxDev;     %#ok<AGROW>
                    if lsqvce.util.verbosity.progress(v)
                        lsqvce.util.log.msg(['%sdone: %d lag(s), ' ...
                            'max|alpha2-1| = %.3e   [%.2f s]\n'], lsqvce.util.verbosity.banner('Joint estimation'), ...
                            jinfo.nLag, jinfo.maxDev, toc(tShape));
                    end
                else
                if lsqvce.util.verbosity.progress(v)
                    lsqvce.util.log.msg('%sre-estimating temporal correlation per stream ...\n', ...
                        lsqvce.util.verbosity.banner('Stochastic model'));
                end
                tShape = tic;
                res = R.residuals('Whitened', true);
                [ovr, ordr, white, pwhite] = obj.applyShape(ovr, streams, res, v);
                hist.arOrder(it,:)   = ordr(:).';    %#ok<AGROW>
                hist.white(it,:)     = white(:).';   %#ok<AGROW>
                hist.whitePval(it,:) = pwhite(:).';  %#ok<AGROW>
                if lsqvce.util.verbosity.progress(v)
                    lsqvce.util.log.msg(['%sdone: %d/%d stream(s) white, ' ...
                        'AR order in [%d, %d]   [%.2f s]\n'], lsqvce.util.verbosity.banner('Stochastic model'), ...
                        sum(white), numel(white), min(ordr), max(ordr), toc(tShape));
                end
                end

                % ---------- outer convergence gate ------------------------
                dsig   = max(abs(vc.sigma(:) - 1));
                balanced = dsig < obj.outerOpt.Tol;
                shapeStable = all(ordr(:).' == prevOrder) || all(white);
                % WP-40 (fixes F-12, gate half): in JOINT mode shapeStable
                % degenerates from pass 2 (the nLag values it compares are
                % constants), so the spectral convergence jinfo.maxDev was
                % never enforced -- a run could declare convergence with the
                % PSD still moving. The gate now additionally requires the
                % LAST joint pass to satisfy the outer Tol.
                jointSettled = true;
                if strcmp(obj.vceMode, 'joint') && ~isempty(hist.jointMaxDev)
                    jointSettled = hist.jointMaxDev(end) <= obj.outerOpt.Tol;
                end
                prevOrder = ordr(:).';

                if lsqvce.util.verbosity.progress(v)
                    % THE CONVERGENCE REPORT IS PER STREAM, not a lone scalar.
                    % max|sigma-1| says how far the WORST stream is from unity
                    % but not WHICH stream, so a loop stalled on one stubborn
                    % stream reads identically to one drifting everywhere. Each
                    % row carries BOTH halves of the stopping rule -- the
                    % variance factor (driven to 1 by the VCE) and the whiteness
                    % of the residuals after whitening -- because the gate needs
                    % both for every stream, so a row failing one but not the
                    % other says immediately which half is stuck. Level 1
                    % deliberately: this is the report you read to judge whether
                    % the run is healthy, not a diagnostic extra.
                    lsqvce.util.log.msg('%s', lsqvce.util.verbosity.heading('Convergence'));
                    for k = 1:numel(streams)
                        [sk, otk] = obj.splitStream(streams{k});
                        lsqvce.util.log.msg('%s', convergenceRow( ...
                            obj.sysLabelOf(sk), otk, ...
                            obj.currentSigma(ovr, sk, otk), ...
                            white(k), pwhite(k)));
                    end
                    % NB: `converged` is not assigned until the gate below, so
                    % recompute the predicate rather than read a stale value.
                    if balanced && jointSettled && (all(white) || shapeStable)
                        tail = 'converged';
                    else
                        tail = 'another pass';
                    end
                    lsqvce.util.log.msg('  max|sigma-1| = %.3e (tol %.1e) -> %s\n', ...
                        dsig, obj.outerOpt.Tol, tail);
                end
                if lsqvce.util.verbosity.detail(v)
                    % Level-1 contribution: current-pass relative info share per
                    % group. Weight-dependent, so it MOVES across passes; the last
                    % pass is authoritative and a still-drifting share signals the
                    % loop has not settled. Free -- reuses vc.r already in hand.
                    lsqvce.util.log.msg('%s', lsqvce.util.verbosity.heading('Contribution'));
                    keys = groupKeyStrings(vc.groups);
                    for g = 1:numel(share)
                        if g <= numel(keys), key = keys{g};
                        else,                key = sprintf('g%d', g);
                        end
                        lsqvce.util.log.msg('%s', ...
                            lsqvce.util.verbosity.labelledRow(key, ...
                                sprintf('%5.1f%%', 100*share(g))));
                    end
                end
                if lsqvce.util.verbosity.detail(v)
                    % The per-STREAM table is already printed at level 1 above.
                    % What level 3 adds is the per-GROUP view: with a coarse
                    % groupComponentsBy several streams SHARE one component, and
                    % the group is what the VCE actually solved for, so the two
                    % views differ -- and the difference is the thing to inspect.
                    lsqvce.util.log.msg('%sper-group sigma: %s\n', ...
                        lsqvce.util.verbosity.banner('Variance components'), ...
                        sigmaLine(vc.groups, vc.sigma));
                end
                if balanced && jointSettled && (all(white) || shapeStable)
                    converged = true;
                    if lsqvce.util.verbosity.progress(v)
                        lsqvce.util.log.msg('%sCONVERGED after %d pass(es)\n', lsqvce.util.verbosity.banner('Convergence'), it);
                    end
                    break
                elseif lsqvce.util.verbosity.progress(v)
                    lsqvce.util.log.msg('%snot converged -> another iteration\n', ...
                        lsqvce.util.verbosity.banner('Convergence'));
                end
            end

            if ~converged
                lsqvce.util.log.warn('lsqvce:adjustVCE:notConverged', ...
                    ['The loop did not converge in %d iteration(s): the ' ...
                     'variance factors and the residual whiteness were not ' ...
                     'both satisfied. Inspect R.history: a stream that never ' ...
                     'whitens needs a higher stochasticModel MaxOrder or ARMA ' ...
                     '(LSQVCE.AR2ARMA); a variance factor that will not settle ' ...
                     'usually means too little group redundancy -- coarsen ' ...
                     'groupComponentsBy.'], nOuter);
            end

            hist.converged = converged;
            hist.nOuter    = size(hist.sigma,1);
            % The whiteness significance travels WITH the history. R.convergence
            % needs it to draw the threshold line, and a Results must be
            % self-describing: requiring the builder object to interpret a
            % finished result is exactly the coupling that moved convergence()
            % onto Results in the first place (HISTORY_v1.md 1.1a decision 6).
            hist.whiteAlpha = obj.shapeOpt.WhiteAlpha;
            noise = obj.noiseSummary(ovr, streams);
            R.setIterative(hist, noise, converged);

            % ---------- PLOT HOOKS: post-model (once, after convergence) ------
            % The stochastic model (AR/PSD per stream) and the VCE weights are now
            % final; fire the per (system,obsType) post-model handles (PSDs, VCE
            % weight bars, ...). vcFinal is the last VCE result.
            obj.firePlots('postModel', streams, R, ovr, vcFinal, hist.nOuter);
        end

        % ---- empirical stochastic modelling (named front door over the loop) ----
        function R = empiricalStochastic(obj, opt)
        %EMPIRICALSTOCHASTIC  Empirical stochastic modelling of the observation noise, as one call.
        %
        %   R = P.empiricalStochastic()
        %   R = P.empiricalStochastic('Name', value, ...)
        %
        %   THE METHOD THIS PACKAGES
        %     Empirical stochastic modelling (Ellmer 2018, TU Graz, §6.5;
        %     operationalised in Kvas et al. 2019, JGR, ITSG-Grace2018) does NOT
        %     assume the observation noise covariance -- it DISCOVERS it from the
        %     postfit residuals and iterates to self-consistency, then VALIDATES
        %     the model by the agreement of formal and empirical parameter errors.
        %     This method runs that exact sequence as a single, named entry point
        %     over the adjustVCE loop. It adds NO new numerics: it is the
        %     canonical ordering of steps the class already implements, plus the
        %     acceptance test (see VALIDATE).
        %
        %   THE STEPS (each mapped to the toolbox call that performs it)
        %     0  HOMOGENISE (preflight, positions only). Kinematic-orbit positions
        %        must already carry their per-epoch PPP covariance as the FIXED D_g
        %        (rotated into the orbit frame: along/across/radial). D is NOT
        %        re-estimated. This method only CHECKS that every multi-axis
        %        stream supplies an epochCov hook / Cep, and warns if one does not
        %        (you would then be whitening 3-D positions with a scalar sigma).
        %        Supply D via the loader's buildDay/cov and epochCov -- see
        %        LSQVCE.MAKELOADER and the class help item (1).
        %     1  START WHITE. sigma_g^2 = 1 and R_{t,g} = I unless a prior override
        %        is already present (a decent a-priori shortens the loop; white is
        %        a legitimate start -- Kvas et al. 2019).
        %     2  SOLVE. One block-arrow GLS pass under the current covariance
        %        (LSQVCE.ADJUSTBLOCK, via blockSolveWith).
        %     3  POSTFIT RESIDUALS. Descaled by the fixed D so their ACV reflects
        %        the temporal law only (LSQVCE.COLLECTRESIDUALS,
        %        LSQVCE.RESIDUALACV).
        %     4  temporal law (colour). Re-fit the per-stream temporal law from the
        %        residual ACV; whiteness decided by a Ljung-Box portmanteau test
        %        (LSQVCE.SELECTARACV, applyShape). PSD<->ACV is LSQVCE.PSD2ACV.
        %     5  variance factor (level). One Helmert variance component per group
        %        (LSQVCE.VCE), folded back into the covariance override.
        %     6  ITERATE 2..5 until BOTH criteria hold: every stream white AND every
        %        sigma_g^2 ~= 1. Typically ~3 outer passes (Kvas et al. 2019).
        %     7  VALIDATE. Compare the FORMAL parameter errors (sqrt(diag(Cxx)))
        %        against an EMPIRICAL error scale, and report per-stream whiteness
        %        p-values and per-group sigma. Good formal/empirical agreement is
        %        the indicator that the stochastic model -- hence the whole
        %        adjustment -- is trustworthy, and removes the need for external
        %        error calibration. See VALIDATE for the returned struct.
        %
        %   RELATION TO RUN
        %     empiricalStochastic() delegates steps 1..6 to RUN (identical loop,
        %     identical Results), then runs step 0 as a preflight and step 7 as a
        %     postflight. Use RUN when you just want the solution; use this method
        %     when you want the procedure named, the homogenisation precondition
        %     checked, and the validation attached. Everything RUN attaches
        %     (R.history, R.noise, R.converged) is present unchanged.
        %
        %   OPTIONS (Name/Value; all optional)
        %     'Homogenise'  logical  run the step-0 preflight check (default true).
        %                            false to skip it (e.g. a positions-free,
        %                            scalar-only problem where D is a plain sigma).
        %     'validate'    logical  run the step-7 validation (default true).
        %     'alpha'       1x1      significance for the whiteness p-values in the
        %                            validation summary (default 0.01, matching the
        %                            loop's WhiteAlpha).
        %     'report'      logical  print the one-screen validation summary at the
        %                            end (default = obj.verbose).
        %
        %   OUTPUT
        %     R   1x1 lsqvce.Results   exactly what RUN returns, PLUS an
        %         R.history.empiricalStochastic sub-struct (only when 'validate' is
        %         true) with the validation fields documented in VALIDATE.
        %         R.report(), R.residuals(), R.acc, R.x, R.sigma0 all work unchanged.
        %
        %   EXAMPLE (synthetic, end to end)
        %     [cfg, L] = demo_grace_dataLoader();
        %     R = lsqvce.adjustVCE(L, cfg.days) ...
        %            .groupComponentsBy('obsType') ...
        %            .empiricalStochastic();    % homogenise -> loop -> validate
        %     R.report();                       % the level-0..4 verdict
        %     R.history.empiricalStochastic     % the validation summary
        %
        %   REFERENCES
        %     Ellmer, M. (2018). Contributions to GRACE Gravity Field Recovery.
        %       Doctoral thesis, TU Graz. doi:10.3217/978-3-85125-646-8. (§6.5)
        %     Kvas, A. et al. (2019). ITSG-Grace2018. J. Geophys. Res. Solid
        %       Earth, 124(8). doi:10.1029/2019JB017415.
        %
        %   See also RUN, VALIDATE, LSQVCE.ITERATENOISEMODEL, LSQVCE.VCE,
        %            LSQVCE.SELECTARACV, LSQVCE.RESIDUALACV,
        %            LSQVCE.SEPARABLEWHITEN, LSQVCE.MAKELOADER
        %
        %   Created 2026-07-16 (CEST, UTC+02:00) with assistance from
        %   Claude (Anthropic), model Claude Opus 4.8.
            % WP-39 c2: typed `arguments` signature; Alpha's and Report's
            % defaults are the CURRENT object state ('Report' is a plain
            % on/off flag whose default follows whether the object prints
            % at all, level >= 1, not the level itself).
        arguments
            obj
            opt.Homogenise (1,1) logical = true
            opt.Validate (1,1) logical = true
            opt.Alpha (1,1) double {mustBePositive, mustBeLessThan(opt.Alpha,1)} = obj.shapeOpt.WhiteAlpha
            opt.Report (1,1) logical = lsqvce.util.verbosity.progress(obj.verbose)
        end
            o = struct('homogenise', opt.Homogenise, 'validate', opt.Validate, ...
                       'alpha', opt.Alpha, 'report', opt.Report);
            v = obj.verbose;

            % ---- step 0: HOMOGENISE (preflight) ----------------------------
            if o.homogenise
                obj.checkHomogenised(v);
            end

            % ---- steps 1..6: the SOLVE<->stochastic model <-> variance components loop ----------------
            if lsqvce.util.verbosity.progress(v)
                lsqvce.util.log.msg('\n===== empirical stochastic modelling =====\n');
            end
            R = obj.run();

            % ---- step 7: VALIDATE (postflight) -----------------------------
            if o.validate
                E = obj.validate(R, 'Alpha', o.alpha, 'Report', o.report);
                h = R.history;
                h.empiricalStochastic = E;
                R.setIterative(h, R.noise, R.converged);
            end
        end


        % ---- the validation (formal vs empirical error) ------------------
        function E = validate(obj, R, opt)
        %VALIDATE  The acceptance test: formal vs empirical parameter errors.
        %
        %   E = P.validate(R)
        %   E = P.validate(R, 'Name', value, ...)
        %
        %   WHY THIS IS THE ACCEPTANCE TEST
        %     A self-consistent covariance model (sigma0 -> 1, residuals white) is
        %     necessary but not sufficient: it can be self-consistent yet wrong in
        %     shape. The external check is that the FORMAL errors of the
        %     estimate -- which depend on both the geometry AND the assumed noise
        %     -- match the EMPIRICAL error scale. Agreement means the noise model
        %     fed into the normals was realistic, so the formal covariance can be
        %     trusted directly and no separate error calibration is needed
        %     (Kvas et al. 2019). This method gathers the evidence for that call;
        %     it does NOT re-solve.
        %
        %   WHAT IT COMPUTES (no new solve; reads R)
        %     * per-group sigma_g^2 and its distance from 1 (the variance-factor evidence;
        %       from R.history.sigma, the last outer pass);
        %     * per-stream whiteness p-value and pass/fail at Alpha (the whiteness
        %       evidence; a Ljung-Box portmanteau on the whitened residual ACV);
        %     * the formal parameter std, formalStd = sqrt(diag(R.sol.Cxx)), and
        %       its ratio to sigma0 (formal errors already carry sigma0, so at a
        %       correct model the empirical/formal ratio ~ 1);
        %     * a single scalar verdict: MODEL ACCEPTED when sigma0 ~ 1, every
        %       group |sigma_g^2 - 1| < SigmaTol, and every stream white at Alpha.
        %
        %   OPTIONS (Name/Value)
        %     'alpha'     1x1  whiteness significance (default obj.shapeOpt.WhiteAlpha).
        %     'SigmaTol'  1x1  tolerance on |sigma_g^2 - 1| for the verdict
        %                      (default = obj.outerOpt.Tol).
        %     'Sigma0Tol' 1x1  tolerance on |sigma0 - 1| for the verdict (default 0.05).
        %     'report'    logical  print the one-screen summary (default obj.verbose).
        %
        %   OUTPUT
        %     E  1x1 struct with fields
        %        .sigma0        1x1      global variance factor of the final solve
        %        .groups        1xG cell  VCE group keys
        %        .sigma         1xG      final per-group sigma_g^2 (last outer pass)
        %        .sigmaOffset   1xG      |sigma_g^2 - 1|
        %        .streams       1xS cell  (satKey|obsType) stream keys
        %        .pWhite        1xS      Ljung-Box p-value per stream (NaN if the
        %                                stream had too few samples to test)
        %        .white         1xS logical  stream white at Alpha
        %        .formalStd     nc x 1   formal std of the common parameters
        %        .empOverForm   1x1      empirical/formal error ratio proxy
        %                                (= sigma0; ~1 when the model is realistic,
        %                                since Cxx already includes sigma0^2)
        %        .accepted      1x1 logical  the overall verdict
        %        .alpha, .sigmaTol, .sigma0Tol  the thresholds used
        %
        %   EXAMPLE
        %     R = P.groupComponentsBy('obsType').run();
        %     E = P.validate(R);
        %     if E.accepted, disp('stochastic model accepted -- formal errors trustworthy'); end
        %
        %   See also EMPIRICALSTOCHASTIC, RUN, LSQVCE.REDUNDANCY, LSQVCE.SYSTEMATIC,
        %            LSQVCE.CHECKPSD
        %
        %   Created 2026-07-16 (CEST, UTC+02:00) with assistance from
        %   Claude (Anthropic), model Claude Opus 4.8.
            % WP-39 c2: typed `arguments` signature (defaults from the
            % CURRENT object state); the arguments block replaces
            % narginchk (a missing R now surfaces MATLAB:minrhs), while
            % the isa contract check on R stays in the body with its
            % documented lsqvce: identifier (the b1-r3 rule).
        arguments
            obj
            R
            opt.Alpha (1,1) double {mustBePositive, mustBeLessThan(opt.Alpha,1)} = obj.shapeOpt.WhiteAlpha
            opt.SigmaTol (1,1) double {mustBePositive} = obj.outerOpt.Tol
            opt.Sigma0Tol (1,1) double {mustBePositive} = 0.05
            opt.Report (1,1) logical = lsqvce.util.verbosity.progress(obj.verbose)
        end
            if ~isa(R, 'lsqvce.Results')
                error('lsqvce:adjustVCE:validateBadResults', ...
                      'validate expects an lsqvce.Results (from run/empiricalStochastic).');
            end
            o = struct('alpha', opt.Alpha, 'sigmaTol', opt.SigmaTol, ...
                       'sigma0Tol', opt.Sigma0Tol, 'report', opt.Report);

            E = struct();
            E.alpha     = o.alpha;
            E.sigmaTol  = o.sigmaTol;
            E.sigma0Tol = o.sigma0Tol;
            E.sigma0    = R.sigma0;

            % ---- variance-factor evidence: per-group sigma from the last outer pass ---
            h = R.history;
            if isstruct(h) && ~isempty(fieldnames(h)) && isfield(h, 'sigma') ...
                    && ~isempty(h.sigma)
                E.groups      = h.groups;
                E.sigma       = h.sigma(end, :);
                E.sigmaOffset = abs(E.sigma - 1);
            else
                E.groups = {};  E.sigma = [];  E.sigmaOffset = [];
            end

            % ---- whiteness evidence: per-stream whiteness on whitened residuals -
            [streams, ~] = obj.discoverStreams();
            E.streams = streams;
            nS = numel(streams);
            E.pWhite = nan(1, nS);
            E.white  = false(1, nS);
            res = R.residuals('Whitened', true);
            for k = 1:nS
                [satKey, ot] = obj.splitStream(streams{k});
                w = pluckResidual(res, ot, obj.sysIdOfKey(satKey));
                if isempty(w) || numel(w) < 8, continue, end
                nAx = obj.nAxesOf(ot);
                [acv, nSamp] = biasedAcvAxes(w, min(60, floor(numel(w)/nAx/4)), nAx);
                [tf, pval]   = ljungBoxWhiteP(acv, nSamp, o.alpha);
                E.pWhite(k)  = pval;
                E.white(k)   = tf;
            end

            % ---- formal parameter errors (the external check) ---------------
            if isfield(R.sol, 'Cxx') && ~isempty(R.sol.Cxx)
                E.formalStd = sqrt(max(diag(R.sol.Cxx), 0));
            else
                E.formalStd = nan(numel(R.x), 1);
            end
            % Cxx already carries sigma0^2, so sigma0 IS the empirical/formal
            % error-scale proxy: sigma0 ~ 1 <=> formal errors match the residual
            % scatter. (A full empirical/formal comparison per parameter needs an
            % external truth ensemble; sigma0 is the single-solution surrogate.)
            E.empOverForm = R.sigma0;

            % ---- the verdict ------------------------------------------------
            okLevel = abs(E.sigma0 - 1) < o.sigma0Tol;
            okScale = isempty(E.sigmaOffset) || all(E.sigmaOffset < o.sigmaTol);
            testable = ~isnan(E.pWhite);
            okShape = ~any(testable) || all(E.white(testable));
            E.accepted = okLevel && okScale && okShape;

            if o.report
                obj.reportEmpirical(E);
            end
        end
    end

    % ====================================================================
    %  PRIVATE MACHINERY
    % ====================================================================
    methods (Access = private)


        function checkHomogenised(obj, verbose)
        %CHECKHOMOGENISED  Step 0: every multi-axis (position) stream must supply
        %   the fixed per-epoch D via an epochCov hook. Warn (do not error) if one
        %   does not -- the loop still runs, but would whiten 3-D positions with a
        %   scalar sigma, which is not the separable model.
            [streams, ~] = obj.discoverStreams();
            hooks = obj.epochCovHooks();          % keys 'satKey|obsType' for multi-axis
            missing = {};
            for k = 1:numel(streams)
                [~, ot] = obj.splitStream(streams{k});
                if obj.nAxesOf(ot) > 1 && ~isKey(hooks, streams{k})
                    missing{end+1} = streams{k}; %#ok<AGROW>
                end
            end
            if isempty(missing)
                if verbose
                    lsqvce.util.log.msg(['[esm   0] homogenisation OK: every ' ...
                        'multi-axis stream carries a fixed per-epoch D.\n']);
                end
            else
                lsqvce.util.log.warn('lsqvce:adjustVCE:notHomogenised', ...
                    ['Homogenisation (step 0): %d multi-axis stream(s) have no per-epoch ' ...
                     'covariance (epochCov/Cep): %s. Positions should carry their ' ...
                     'PPP epoch covariance as the FIXED D (rotated to ' ...
                     'along/across/radial); without it the loop whitens them with ' ...
                     'a scalar sigma. Supply D via the loader (see ' ...
                     'LSQVCE.MAKELOADER) or set ''Homogenise'',false to silence ' ...
                     'this if D is deliberately scalar.'], ...
                     numel(missing), strjoin(missing, ', '));
            end
        end

        function reportEmpirical(~, E)
        %REPORTEMPIRICAL  One-screen empirical-stochastic-model validation summary.
            L = @lsqvce.util.log.msg;
            L('\n----- empirical stochastic-model validation -----\n');
            L('  global sigma0            : %.4f   (want ~1; |.-1| < %.3g)\n', ...
                    E.sigma0, E.sigma0Tol);
            if ~isempty(E.sigma)
                L('  Variance components  per-group sigma^2 (last pass):\n');
                for g = 1:numel(E.sigma)
                    flag = ''; if E.sigmaOffset(g) >= E.sigmaTol, flag = '  <-- off'; end
                    L('    %-20s %8.4f%s\n', E.groups{g}, E.sigma(g), flag);
                end
            end
            L('  Stochastic model  per-stream whiteness (Ljung-Box, alpha = %.3g):\n', E.alpha);
            for s = 1:numel(E.streams)
                if isnan(E.pWhite(s))
                    L('    %-20s   p =   n/a   (too few samples)\n', E.streams{s});
                else
                    tag = 'WHITE'; if ~E.white(s), tag = 'coloured'; end
                    L('    %-20s   p = %6.3f   %s\n', E.streams{s}, E.pWhite(s), tag);
                end
            end
            L('  formal param std         : min %.3g  median %.3g  max %.3g\n', ...
                    min(E.formalStd), median(E.formalStd), max(E.formalStd));
            L('  empirical/formal proxy   : %.4f   (= sigma0; ~1 when realistic)\n', ...
                    E.empOverForm);
            if E.accepted
                L('  VERDICT                  : MODEL ACCEPTED -- formal errors trustworthy\n');
            else
                L('  VERDICT                  : NOT ACCEPTED -- inspect the off/coloured rows above\n');
            end
            L('----------------------------------------------\n');
        end

        function [streams, initCov] = discoverStreams(obj)
        %DISCOVERSTREAMS  List (satKey,obsType) pairs and their starting cov spec.
            streams = {};  initCov = {};
            L = obj.loaderList();
            for i = 1:numel(L)
                satKey = obj.satKeyOf(L{i}.sysID);
                for k = 1:numel(L{i}.obsTypes)
                    ot = L{i}.obsTypes{k};
                    streams{end+1} = sprintf('%s|%s', satKey, ot);  %#ok<AGROW>
                    if isfield(L{i}.cov, ot)
                        initCov{end+1} = L{i}.cov.(ot);              %#ok<AGROW>
                    else
                        initCov{end+1} = struct('type','diag','sigma',1); %#ok<AGROW>
                    end
                end
            end
        end

        function L = loaderList(obj)
        %LOADERLIST  Normalise obj.systems to a cell of Loader structs.
            if iscell(obj.systems)
                L = obj.systems;
            elseif isstruct(obj.systems)
                L = cell(1, numel(obj.systems));
                for i = 1:numel(obj.systems), L{i} = obj.systems(i).loader; end
            else
                error('lsqvce:adjustVCE:badSystems', ...
                    'systems must be a cell of Loaders or an adjustBlock struct array.');
            end
        end

        function acc = buildAccumulator(obj, sigma, ovr, lvl)
        %BUILDACCUMULATOR  The accFcn VCE calls: fresh accumulator, current cov.
        %   INPUT (beyond sigma/ovr)
        %     lvl  1x1 double  verbosity LEVEL for this accumulator, 0..3. Passed
        %                      in rather than read from obj because the two call
        %                      sites differ: BLOCKSOLVEWITH accumulates ONCE per
        %                      outer pass, so its per-day detail is exactly what
        %                      level 2 is for; the accFcn handed to VCE is rebuilt
        %                      on EVERY inner iteration (~4 per pass, ~30 per run),
        %                      so the same detail there is a flood and is held back
        %                      to level 3.
        %
        %   sigma is VCE's per-group scale; ovr carries the current AR shape and
        %   the previously folded scale. VCE multiplies sigma INTO the group, so
        %   here we thread sigma through opt.covOverride as an extra scale.
            if nargin < 4 || isempty(lvl)
                lvl = 0;      % default SILENT: a nested accumulator speaks only
            end                % when a call site explicitly asks it to.
            uo = obj.covOverrideOpt(ovr, sigma);
            acc = lsqvce.BlockAccumulator(obj.nc, ...
                      'Backend', obj.blockOpt.Backend, 'KeepLocal', true, ...
                      'Reg', obj.blockOpt.Reg, 'Parallel', obj.blockOpt.Parallel, ...
                      ... % The level is passed IN by the caller, because the two
                      ... % call sites need different ones -- see below. It used to
                      ... % be hardcoded 'Verbose', false, which muted the
                      ... % accumulator on every pass regardless of what the user
                      ... % asked for, so the level-2 per-day detail could never
                      ... % appear through adjustVCE at all.
                      'Verbose', lvl, 'UserOpt', uo);
            specs = obj.blockSpecs();
            acc.addSystems(specs);
        end

        function R = blockSolveWith(obj, ovr)
        %BLOCKSOLVEWITH  One full solve with the current covariance override,
        %   bridged into an lsqvce.Results exactly as LSQVCE.ADJUSTBLOCK does. We
        %   cannot use adjustBlock() directly because it does not expose UserOpt (the
        %   channel the covariance override travels on), so we replicate its
        %   six-step dance here with the accumulator built by buildAccumulator.
            % THE solve for this pass: one accumulation, so its per-day detail
            % is what the stage level exists to show. Silent below that.
            lvl = 0;
            if lsqvce.util.verbosity.detail(obj.verbose), lvl = obj.verbose; end
            acc = obj.buildAccumulator([], ovr, lvl);     % [] -> no extra scale
            sol = acc.solve();
            [dsol, batchFcn, nB] = lsqvce.diagnoseAccum(acc, sol);
            D = lsqvce.diagnose(dsol, batchFcn, 'nBatch', nB);
            R = lsqvce.Results(sol, D, acc, dsol, batchFcn, nB, obj.blockSpecs());
        end

        function uo = covOverrideOpt(obj, ovr, sigma)
        %COVOVERRIDEOPT  Pack ovr (+ optional VCE sigma) into a UserOpt struct.
        %   Emits opt.covOverride.bySys.(satKey).(obsType) so each stream gets its
        %   OWN AR law (from ovr, the temporal law) and its own scale. When sigma is
        %   given (VCE's accFcn call) it is applied as an extra multiplicative
        %   std on top of ovr, mapped to each stream through the VCE group key so
        %   the alignment is exact regardless of group ordering.
            [streams,~] = obj.discoverStreams();
            bySys = struct();

            % base: the accumulated temporal law (+ any folded scale) from ovr
            for k = 1:numel(streams)
                [satKey, ot] = obj.splitStream(streams{k});
                if isfield(ovr, satKey) && isfield(ovr.(satKey), ot)
                    if ~isfield(bySys, satKey), bySys.(satKey) = struct(); end
                    bySys.(satKey).(ot) = ovr.(satKey).(ot);
                end
            end

            % extra VCE scale, mapped stream -> group key -> sigma index
            if ~isempty(sigma)
                keyFcn = obj.vceGroupKeyFcn();
                gkeys  = obj.vceGroupKeyList(streams, keyFcn);   % 1xnStream keys
                uk     = unique(gkeys, 'stable');                % group order
                for k = 1:numel(streams)
                    [satKey, ot] = obj.splitStream(streams{k});
                    gi = find(strcmp(uk, gkeys{k}), 1);
                    if isempty(gi) || gi > numel(sigma), continue, end
                    sExtra = sqrt(max(sigma(gi), eps));
                    if ~isfield(bySys, satKey) || ~isfield(bySys.(satKey), ot)
                        bySys = setField(bySys, satKey, ot, 'sigma', 1);
                    end
                    base = 1;
                    if isfield(bySys.(satKey).(ot),'sigma')
                        base = bySys.(satKey).(ot).sigma;
                    end
                    bySys.(satKey).(ot).sigma = base * sExtra;
                end
            end

            uo = struct('covOverride', struct('bySys', bySys));
            % In JOINT mode the estimated noise model is a PSD. covModel reads
            % psdMode from the UserOpt (not from the covariance spec), and its
            % DEFAULT is 'ar' -- which fits an AR to the estimated spectrum and
            % can produce a NON-STATIONARY polynomial whose autocovariance does
            % not exist. Ellmer builds the Toeplitz covariance from the PSD
            % DIRECTLY (eqs. 6.5.6, 6.5.17), with no AR step and hence no
            % stationarity requirement, so the joint mode says so explicitly.
            if strcmpi(obj.vceMode, 'joint')
                uo.psdMode = obj.jointOpt.PsdMode;
            end
        end

        function fh = vceGroupKeyFcn(obj)
        %VCEGROUPKEYFCN  The group-key function actually used, as a handle over
        %   (sysID,obsType,arc,day). Mirrors vceGroupKey but always a handle.
            k = obj.vceGroupKey();
            if isa(k,'function_handle')
                fh = k;
            else
                switch lower(char(k))
                    case 'obstype', fh = @(s,ot,a,d) char(string(ot));
                    case 'system',  fh = @(s,ot,a,d) char(string(s));
                    otherwise,      fh = @(s,ot,a,d) sprintf('%s|%s', string(s), string(ot));
                end
            end
        end

        function keys = vceGroupKeyList(obj, streams, keyFcn)
        %VCEGROUPKEYLIST  Group key per stream (arc/day-independent groupings).
            L = obj.loaderList();
            sysOf = containers.Map('KeyType','char','ValueType','any');
            for i = 1:numel(L)
                sysOf(obj.satKeyOf(L{i}.sysID)) = L{i}.sysID;
            end
            keys = cell(1, numel(streams));
            for k = 1:numel(streams)
                [satKey, ot] = obj.splitStream(streams{k});
                sysID = satKey;
                if isKey(sysOf, satKey), sysID = sysOf(satKey); end
                keys{k} = char(string(keyFcn(sysID, ot, 1, 1)));
            end
        end

        function specs = blockSpecs(obj)
        %BLOCKSPECS  Build the spec array the accumulator's addSystems wants.
            L = obj.loaderList();
            % Same field set as the struct assigned below, and as
            % ADJUSTBLOCK's oneSpec: a mismatch raises 'Subscripted
            % assignment between dissimilar structures' at RUN time, not
            % at parse time, so nothing catches it until the loop runs.
            tmpl = struct('sysID',[], 'nl',[], 'dataFcn',[], ...
                          'dayList',[], 'name','');
            specs = repmat(tmpl, 1, numel(L));
            for i = 1:numel(L)
                dl = obj.days;
                if isfield(L{i}, 'days') && ~isempty(L{i}.days), dl = L{i}.days; end
                nm = '';
                if isfield(L{i}, 'name') && ~isempty(L{i}.name)
                    nm = char(L{i}.name);      % readable label for the log
                end
                specs(i) = struct('sysID', L{i}.sysID, 'nl', double(L{i}.nl), ...
                                  'dataFcn', L{i}.closure, 'dayList', dl, ...
                                  'name', nm);
            end
        end

        function ovr = applyScale(obj, ovr, streams, vc)
        %APPLYSCALE  Fold VCE's per-group sigma into the covariance override.
        %   Maps each stream to its VCE group via the same key function used to
        %   pin the group order, so sigma indexing is exact for any grouping.
            gkeys = obj.vceGroupKeyList(streams, obj.vceGroupKeyFcn());
            vcKeys = cellstr(string(vc.groups));
            for k = 1:numel(streams)
                [satKey, ot] = obj.splitStream(streams{k});
                g = find(strcmp(vcKeys, gkeys{k}), 1);
                if isempty(g) || g > numel(vc.sigma), continue, end
                sg = sqrt(max(vc.sigma(g), eps));
                if ~isfield(ovr, satKey), ovr.(satKey) = struct(); end
                if ~isfield(ovr.(satKey), ot), ovr.(satKey).(ot) = struct(); end
                base = 1;
                if isfield(ovr.(satKey).(ot),'sigma'), base = ovr.(satKey).(ot).sigma; end
                ovr.(satKey).(ot).sigma = base * sg;
            end
        end

        function [ovr, info] = applyJoint(obj, ovr, streams, R, v)
        %APPLYJOINT  Co-estimate the temporal stochastic model and the variance components in one step.
        %   The 'joint' counterpart of APPLYSHAPE (Ellmer 2018 §6.5.7 steps 6-9).
        %   Instead of fitting an AR law per stream and scaling separately, this
        %   runs the REAL spectral estimator: one shared Monte-Carlo probe
        %   U^-1*Zbar (LSQVCE.SHAREDPROBE) reused across every arc, the per-lag
        %   redundancies s_n = tr(Rbar V_n) and residual sums Omega_n per arc
        %   (LSQVCE.REDUNDANCYLAGS), and the per-frequency amplitude VCE
        %   (LSQVCE.PSDVCE). The amplitudes and the scale therefore move
        %   TOGETHER -- they are coupled variance components, and updating them
        %   alternately can settle on a different model.
        %
        %   The estimated PSD is written into ovr as an explicit 'psd' covariance
        %   override, the same channel APPLYSHAPE uses for its AR law.
        %
        %   Returns INFO with .nLag, .maxDev (max|alpha2-1| over streams) and
        %   .converged, used for the history columns and the verbose banner.
            nLag   = obj.jointOpt.NumLags;
            maxDev = 0;
            updated = {};
            acc    = R.acc;
            sol    = R.sol;

            % ---- step 6: the shared probe, ONCE for every arc --------------
            N = acc.reducedNormalMatrix();
            [U, pd] = chol((N + N.')/2);
            if pd > 0
                lsqvce.util.log.warn('lsqvce:adjustVCE:jointNotPD', ...
                    ['The reduced normal matrix is not positive definite, so the ' ...
                     'joint spectral update was skipped this pass. The previous ' ...
                     'noise model is kept.']);
                info = struct('nLag', nLag, 'maxDev', NaN, ...
                              'updated', {{}}, 'nUpdated', 0, 'converged', false);
                return
            end
            sp = lsqvce.sharedProbe(U, 'numProbes', obj.jointOpt.Probes, ...
                                    'seed', obj.jointOpt.Seed);

            % ---- steps 7-8: per-stream arc sweep --------------------------
            for k = 1:numel(streams)
                [satKey, ot] = obj.splitStream(streams{k});
                sysID = obj.sysIdOfKey(satKey);
                [omegaT, sT, nArc, skip] = obj.sweepStreamLags(acc, sol, sp.UiZ, ...
                                                         sysID, ot, nLag);
                if nArc == 0
                    % A stream contributing NO arcs is silently unmodelled --
                    % its noise law simply never updates while every other
                    % stream moves, which quietly biases the whole solution.
                    % Say so, with the reason, rather than dropping it.
                    lsqvce.util.log.warn('lsqvce:adjustVCE:jointStreamSkipped', ...
                       ['Joint update: stream %s contributed NO usable arcs. ' ...
                        'Of %d arc(s): %d empty, %d with neither Cfcn nor Wfcn, ' ...
                        '%d Cfcn threw, %d non-PD covariance, %d Wfcn threw, ' ...
                        '%d Wfcn wrong shape, %d Wfcn singular, %d axis-pooling ' ...
                        'empty.%s Its noise model is left UNCHANGED this ' ...
                        'pass.'], streams{k}, ...
                        skip.seen, skip.empty, skip.noCov, skip.cfcnErr, ...
                        skip.notPD, skip.wfcnErr, skip.wfcnShape, ...
                        skip.wfcnSing, skip.axisPool, firstErrText(skip.lastErr));
                    continue
                end
                updated{end+1} = streams{k};   %#ok<AGROW>

                % starting covariance for this stream: the user's C0, else the
                % measured lag-0-normalised Omega (never a delta -- that start
                % is a degenerate fixed point, see LSQVCE.PSDVCE).
                C0 = obj.jointOpt.C0;
                L  = size(omegaT, 1);
                if isempty(C0)
                    m0 = mean(omegaT, 2);
                    C0 = m0 / max(abs(m0(1)), eps);
                else
                    C0 = C0(:);
                    L  = min(L, numel(C0));
                    C0 = C0(1:L);  omegaT = omegaT(1:L,:);  sT = sT(1:L,:);
                end

                % ---- step 9: per-frequency amplitude VCE ------------------
                pv = lsqvce.psdVCE(C0 .* omegaT, C0 .* sT, C0, ...
                                   'Norm', obj.jointOpt.Norm, ...
                                   'Smooth', obj.jointOpt.Smooth);
                maxDev = max(maxDev, pv.maxDev);

                % Write the estimated PSD as this stream's covariance override.
                %
                % TWO DIFFERENT NUMBERS, do not conflate them:
                %   nLag  -- how many lags the ESTIMATOR resolves. Small is
                %            correct: each lag is estimated from finite data, so
                %            more lags means a noisier alpha2_j.
                %   nf    -- the frequency grid the COVARIANCE is built on. This
                %            must satisfy nf >= m/2 + 1 for the Toeplitz matrix
                %            to be full rank on an m-row window, otherwise it is
                %            rank deficient by construction.
                % The fix is NOT to estimate more lags but to INTERPOLATE the
                % coarse estimate onto a fine enough grid: verified that a
                % 12-lag spectrum resampled to 451 frequencies gives min eig
                % +5.5e-02 at m = 900, where the coarse grid gives -2.9e-13.
                sig  = obj.currentSigma(ovr, satKey, ot);
                % Frequency axis in PHYSICAL units. Nyquist is 1/(2*dt), NOT
                % 0.5: building the grid on normalised frequency claims power
                % above Nyquist that the sampled process cannot carry, and
                % psd2acv rightly warns that the resulting r will not match the
                % data. dt comes from the arcs actually swept.
                dt = skip.dt;
                if ~isfinite(dt) || dt <= 0, dt = 1; end
                fNyq = 1/(2*dt);
                fEst = (0:L-1).' / (2*(L-1)) * (2*fNyq);   % 0 .. fNyq
                nf   = max(L, obj.covGridSize(R.acc));     % grid to build on
                fCov = linspace(0, fNyq, nf).';
                Scov = interp1(fEst, pv.Shat(:), fCov, 'linear', 'extrap');
                Scov = max(Scov, eps * max(Scov));       % keep it a valid PSD
                cm   = struct('type','psd', 'f', fCov, 'S', Scov, 'sigma', sig, ...
                              'psdMode', obj.jointOpt.PsdMode);
                % A rough estimated PSD can fit a NON-STATIONARY AR, whose
                % autocovariance does not exist -- covModel would then throw on
                % every later use (lsqvce:armaACV:nonStationary) and abort the
                % whole run. Probe the model here and KEEP THE PREVIOUS one for
                % this stream if the new spectrum is not usable.
                if ~obj.psdIsUsable(cm)
                    lsqvce.util.log.warn('lsqvce:adjustVCE:jointPsdRejected', ...
                       ['Joint update: the PSD estimated for stream %s does not ' ...
                        'yield a usable (stationary) noise model, so the ' ...
                        'previous model is kept for this stream. Raise ' ...
                        '''numLags'' or start from a smoother ''c0'' if this ' ...
                        'persists.'], streams{k});
                    continue
                end
                if ~isfield(ovr, satKey), ovr.(satKey) = struct(); end
                ovr.(satKey).(ot) = cm;
                if lsqvce.util.verbosity.detail(v)
                    lsqvce.util.log.msg(['%s%s: %d arc(s), %d ' ...
                        'frequencies, max|alpha2-1| = %.3e\n'], lsqvce.util.verbosity.banner('Joint estimation'), ...
                        streams{k}, nArc, L, pv.maxDev);
                end
            end
            info = struct('nLag', nLag, 'maxDev', maxDev, ...
                          'updated', {updated}, 'nUpdated', numel(updated), ...
                          'converged', maxDev <= obj.outerOpt.Tol);
        end

        function nf = covGridSize(obj, acc)
        %COVGRIDSIZE  Frequency-grid size needed so the covariance is full rank.
        %   A PSD sampled at nf frequencies spans exactly 2*(nf-1) independent
        %   lags. Any window longer than that yields a RANK DEFICIENT covariance
        %   (see LSQVCE.COVMODEL's resolution check), so
        %       nf >= m/2 + 1
        %   with m the longest window the covariance is ever built on. That is
        %   the BATCH length, not the arc length: LSQVCE.COVMODEL uses
        %   min(BatchSize + BatchOverlap, numel(t)), and a batch may span several
        %   arcs.
        %
        %   This sizes the GRID THE COVARIANCE IS BUILT ON, not the number of
        %   lags the estimator resolves. The estimate stays coarse (jointOpt
        %   NumLags) and is INTERPOLATED onto this finer grid; estimating
        %   hundreds of lags instead would only make alpha2_j noisier and cost a
        %   dense factorisation per probe. Falls back to 128 when no length can
        %   be determined.
            mMax = 0;
            % (a) the arcs themselves
            try
                ids = acc.systemIDs();
                for k = 1:numel(ids)
                    st = acc.getStore(ids(k));
                    uo = st.uo;
                    d  = st.dayList(1);
                    meta = st.dataFcn('meta', d, [], uo);
                    for ot = 1:numel(meta.obsTypes)
                        a = meta.arcs{ot};
                        if isempty(a), continue, end
                        uo.obsType = meta.obsTypes{ot};
                        uo.day = d;  uo.arc = a(1);
                        arc = st.dataFcn(d, meta.obsTypes{ot}, a(1), uo);
                        if ~isempty(arc) && isfield(arc,'A')
                            mMax = max(mMax, size(arc.A,1));
                        end
                    end
                end
            catch
                mMax = 0;
            end
            % (b) ARCS REMAIN INDEPENDENT -- nothing here models correlation
            % ACROSS arcs, and a batch never spans two arcs (see
            % BlockAccumulator.accumArc, "split ONE arc into row batches").
            % But covModel PRECOMPUTES its autocovariance vector r once per
            % stream, sized min(BatchSize + overlap, numel(t)) where t is the
            % whole per-day epoch vector, so that one r can serve every batch.
            % psd2acv must therefore be able to generate r out to that many
            % lags, which is longer than a single arc. Scale by the arcs per day
            % to size the GRID for that precompute; it does not widen the
            % covariance model itself.
            try
                nArcDay = 1;
                ids = acc.systemIDs();
                if ~isempty(ids)
                    st = acc.getStore(ids(1));
                    md = st.dataFcn('meta', st.dayList(1), [], st.uo);
                    for ot = 1:numel(md.obsTypes)
                        nArcDay = max(nArcDay, numel(md.arcs{ot}));
                    end
                end
                mMax = mMax * nArcDay;
            catch
                % keep the single-arc estimate
            end
            if mMax < 2
                nf = 128;
            else
                nf = floor(mMax/2) + 1;
            end
        end

        function [om, sn] = lagsPooledOverAxes(obj, Abar, W, UiZ, ebar, nAx, nE, nLag, Bbar, LlS) %#ok<INUSL>
        %LAGSPOOLEDOVERAXES  Per-lag Omega and s down EPOCHS, pooled over axes.
        %   Rows are epoch-major (row = (epoch-1)*nAx + axis). Under the
        %   separable model the per-epoch de-scale leaves exactly R_t kron I_d,
        %   so the axes are decorrelated and share ONE temporal law: estimate
        %   each axis independently down epochs and average, mirroring
        %   BIASEDACVAXES on the stochastic-model side.
        %
        %   Slicing matters, and NOT the obvious way: W(rows,rows) stays upper
        %   triangular but is NOT a valid factor of the sliced covariance --
        %   Sigma = W'W does not survive taking a subset of rows/columns
        %   (verified: max |Wa'Wa - Sigma(rows,rows)| = 2.8 on a 18x18 test).
        %   The correct per-axis whitener comes from reconstructing
        %   Sigma = W'*W, slicing THAT, and re-factoring it.
        %
        %   INPUT   Abar  m x P    decorrelated design
        %           W     m x m    factor with Sigma = W'W
        %           UiZ   P x Z    the shared probe block
        %           ebar  m x 1    decorrelated residual
        %           nAx   1x1      axes per epoch
        %           nE    1x1      epochs in this arc
        %           nLag  1x1      lags wanted
        %   OUTPUT  om, sn  nLagUsed x 1   pooled Omega_n and s_n, or [] if no
        %                                  axis produced a usable estimate
            maxLag = min(nLag, nE);
            accOm  = zeros(maxLag, 1);
            accS   = zeros(maxLag, 1);
            nOk    = 0;
            Sig    = W.' * W;                    % the arc covariance
            for ax = 1:nAx
                rows = ax : nAx : (nE*nAx);      % this axis, down epochs
                if numel(rows) < 2, continue, end
                Sa = Sig(rows, rows);
                Sa = (Sa + Sa.')/2;
                [Wa, pdv] = chol(Sa);            % re-factor, do NOT slice W
                if pdv > 0, continue, end
                try
                    % WP-33: the local rows are sliced exactly like the common
                    % ones (symmetric treatment under the separable model).
                    if nargin >= 10 && ~isempty(Bbar) && size(Bbar,2) > 0
                        ra = lsqvce.redundancyLags(Abar(rows,:), Wa, UiZ, ...
                                 ebar(rows), maxLag, ...
                                 'LocalBbar', Bbar(rows,:), 'LocalLl', LlS);
                    else
                        ra = lsqvce.redundancyLags(Abar(rows,:), Wa, UiZ, ...
                                                   ebar(rows), maxLag);
                    end
                catch
                    continue
                end
                o = padLags(ra.omega, maxLag);
                s = padLags(ra.s    , maxLag);
                accOm = accOm + o;
                accS  = accS  + s;
                nOk   = nOk + 1;
            end
            if nOk == 0
                om = [];  sn = [];
            else
                om = accOm / nOk;
                sn = accS  / nOk;
            end
        end

        function ok = psdIsUsable(obj, cm)
        %PSDISUSABLE  Can this covariance spec actually build a covariance?
        %   With psdMode 'toeplitz' the PSD becomes a Toeplitz Sigma directly and
        %   essentially always works. With 'ar'/'arma' an estimated PSD can fit a
        %   NON-STATIONARY AR (max |root| >= 1) whose autocovariance does not
        %   exist; LSQVCE.COVMODEL then throws every time the model is used, which
        %   would abort the accumulation rather than merely degrade it. Building
        %   the model once here, on a short probe window AND WITH THE SAME
        %   psdMode the run will use, converts that fatal case into a rejected
        %   update.
        %
        %   INPUT   cm  1x1 struct   a covariance spec (type 'psd', .f, .S, .sigma)
        %   OUTPUT  ok  1x1 logical  true when a finite covariance can be built
            ok = false;
            probeOpt = struct('psdMode', obj.jointOpt.PsdMode);
            try
                % CHEAP probe. Factorising a full-size dense Toeplitz matrix here
                % costs O(m^3) PER STREAM PER PASS and dominated the runtime
                % (74 s vs 17 s when this probed 256 rows). The covariance is now
                % built on a grid sized by covGridSize, which guarantees full
                % rank by construction, so the probe only has to catch a spec
                % that cannot build a covariance AT ALL (non-finite amplitudes,
                % a non-stationary AR under psdMode 'ar'). A short window
                % suffices for that.
                nProbe = 32;
                t = (0:nProbe-1).' ;
                [Cfcn, Wfcn] = lsqvce.covModel(cm, t, probeOpt);
                if ~isempty(Cfcn)
                    C = Cfcn(1:nProbe, probeOpt);
                    if isvector(C), C = diag(C(:)); end
                    ok = all(isfinite(C(:)));
                elseif ~isempty(Wfcn)
                    Aw = Wfcn(1:nProbe, eye(nProbe), zeros(nProbe,0), ...
                              zeros(nProbe,1), probeOpt);
                    ok = all(isfinite(Aw(:)));
                end
            catch
                ok = false;        % any failure here means "do not adopt it"
            end
        end

        function [omegaT, sT, nArc, skip] = sweepStreamLags(obj, acc, sol, UiZ, sysID, ot, nLag) %#ok<INUSL>
        %SWEEPSTREAMLAGS  Per-lag Omega and s for every arc of ONE stream.
        %   Re-streams the arcs of (sysID, obsType) exactly as LSQVCE.REDUNDANCY
        %   does, decorrelates each with its own W (Sigma = W'W), forms the
        %   residual, and calls LSQVCE.REDUNDANCYLAGS. Returns one COLUMN per
        %   arc, zero-padded to nLag (eqs. 6.5.22-6.5.23 set the contribution to
        %   zero for lags beyond the arc length).
        %
        %   SKIP is a struct counting arcs dropped for each reason (empty, no
        %   Cfcn, non-PD covariance). A stream whose arcs are ALL dropped would
        %   otherwise vanish from the estimate without trace -- APPLYJOINT warns
        %   using these counts rather than silently leaving the stream unmodelled.
            omegaT = zeros(nLag, 0);
            sT     = zeros(nLag, 0);
            nArc   = 0;
            skip   = struct('empty',0, 'noCov',0, 'notPD',0, 'seen',0, ...
                            'cfcnErr',0, 'wfcnErr',0, 'wfcnShape',0, ...
                            'wfcnSing',0, 'axisPool',0, 'lastErr','', 'dt',NaN);
            x      = sol.x(:);
            st     = acc.getStore(sysID);
            uo     = st.uo;
            % WP-33: this system's local solution and Schur pieces, ONCE per
            % sweep -- exactly as LSQVCE.REDUNDANCY recovers them. The arc
            % residual becomes ebar = W^-T(l - A x - B y) and the per-lag
            % projector absorbs the local block exactly (jointArcParts +
            % redundancyLags 'LocalBbar'/'LocalLl'); nl = 0 stays verbatim.
            nlS = size(st.Nlc, 1);
            if nlS > 0
                ySys = st.Nll_r \ (st.nl_r - st.Nlc * x(1:size(st.Nlc,2)));
                Gsys = st.Nll_r \ st.Nlc;
                [LlS, plS] = chol((st.Nll_r + st.Nll_r.')/2, 'lower');
                if plS > 0
                    error('lsqvce:adjustVCE:localNotPD', ...
                         ['System %s: the local normal matrix is not ' ...
                          'positive definite -- the local block cannot be ' ...
                          'absorbed into the projector.'], string(sysID));
                end
            else
                ySys = [];  Gsys = [];  LlS = [];
            end
            for d = st.dayList(:).'
                meta = st.dataFcn('meta', d, [], uo);
                oti  = find(strcmp(meta.obsTypes, ot), 1);
                if isempty(oti), continue, end
                for a = meta.arcs{oti}(:).'
                    skip.seen = skip.seen + 1;
                    uo.obsType = ot;  uo.day = d;  uo.arc = a;
                    arc = st.dataFcn(d, ot, a, uo);
                    if isempty(arc) || ~isfield(arc,'l') || isempty(arc.l)
                        skip.empty = skip.empty + 1;  continue
                    end
                    m   = size(arc.A, 1);
                    % The explicit triangular W with Sigma = W'W is what
                    % REDUNDANCYLAGS needs (eqs. 6.5.48, 6.5.54). Two routes:
                    %   (a) Cfcn -> the dense arc covariance -> chol.
                    %   (b) Wfcn only. This is the NORMAL case for multi-axis
                    %       position streams: a 'separable' covariance
                    %       D*(R_t kron I_d)*D is applied as an OPERATOR and
                    %       never formed densely, so covModel returns Wfcn and
                    %       leaves Cfcn empty. Recover W by whitening an
                    %       identity block -- W*I = W -- which costs one extra
                    %       m x m apply per arc and works for ANY whitener.
                    W = [];
                    Sig = [];
                    if isfield(arc,'Cfcn') && ~isempty(arc.Cfcn)
                        % Cfcn can THROW: it may build the covariance from an AR
                        % fitted to the current (estimated) PSD, and a rough PSD
                        % gives a non-stationary AR whose autocovariance does not
                        % exist (lsqvce:armaACV:nonStationary). That must not
                        % abort the sweep -- fall through to Wfcn instead.
                        try
                            Sig = arc.Cfcn(1:m, uo);
                        catch ME
                            skip.cfcnErr = skip.cfcnErr + 1;
                            if isempty(skip.lastErr)
                                skip.lastErr = sprintf('%s: %s', ...
                                                       ME.identifier, ME.message);
                            end
                            Sig = [];
                        end
                    end
                    if ~isempty(Sig)
                        if isvector(Sig), Sig = diag(Sig(:)); end
                        [W, pdv] = chol((Sig + Sig.')/2);
                        if pdv > 0
                            skip.notPD = skip.notPD + 1;  W = [];
                        end
                    end
                    if isempty(W) && isfield(arc,'Wfcn') && ~isempty(arc.Wfcn)
                        % Wfcn APPLIES the whitener: Wop*A is already Abar, and
                        % Wop'*Wop = Sigma^-1. REDUNDANCYLAGS wants the factor W
                        % with Sigma = W'W, and Wop = inv(W'), so W = inv(Wop)'.
                        % Recover Wop by whitening an identity block.
                        Wop = [];
                        try
                            % Pass correctly SHAPED dummies: a whitener filters
                            % A, B and l together and may index all three, so []
                            % can throw. B is m x 0 (no local columns) and l is
                            % m x 1; only the A slot carries the identity we
                            % actually want whitened.
                            Wop = arc.Wfcn(1:m, eye(m), zeros(m,0), zeros(m,1), uo);
                        catch ME
                            skip.wfcnErr = skip.wfcnErr + 1;
                            if isempty(skip.lastErr)
                                skip.lastErr = sprintf('%s: %s', ...
                                                       ME.identifier, ME.message);
                            end
                        end
                        if isempty(Wop)
                            W = [];
                        elseif ~isequal(size(Wop), [m m])
                            skip.wfcnShape = skip.wfcnShape + 1;
                            if isempty(skip.lastErr)
                                skip.lastErr = sprintf(...
                                    'Wfcn returned %dx%d, expected %dx%d', ...
                                    size(Wop,1), size(Wop,2), m, m);
                            end
                            W = [];
                        elseif ~all(isfinite(Wop(:))) || rcond(Wop) < eps
                            skip.wfcnSing = skip.wfcnSing + 1;  W = [];
                        else
                            W = inv(Wop).';                  %#ok<MINV>
                        end
                    end
                    if isempty(W)
                        % Nothing above produced a usable factor. Count it as a
                        % missing covariance only if neither route even existed;
                        % the specific failure counters were already bumped.
                        if (~isfield(arc,'Cfcn') || isempty(arc.Cfcn)) && ...
                           (~isfield(arc,'Wfcn') || isempty(arc.Wfcn))
                            skip.noCov = skip.noCov + 1;
                        end
                        continue
                    end
                    % WP-33: full-vector residual + reduced common design.
                    Barc = zeros(m, 0);
                    if isfield(arc,'B') && ~isempty(arc.B), Barc = arc.B; end
                    [Abar, Bbar, ebar] = jointArcParts(arc.A, Barc, ...
                        arc.l(:), W, x(1:size(arc.A,2)), ySys, Gsys);

                    % AXIS STRIDE. Rows are EPOCH-MAJOR:
                    %     row = (epoch-1)*nAx + axis
                    % so a lag computed on ROWS counts axes, not epochs -- row
                    % lag 1 is the next AXIS of the same epoch. Estimating that
                    % way stretches the temporal law by a factor nAx and folds
                    % axis structure into what must be a purely temporal ACV.
                    % The separable model gives every axis the SAME R_t (after
                    % the per-epoch de-scale the covariance is exactly
                    % R_t kron I_d), so the correct estimate takes each axis
                    % down epochs and POOLS them -- the same thing the stochastic model
                    % path does via biasedAcvAxes. For nAx == 1 this reduces to
                    % the previous behaviour exactly.
                    nAx = 1;
                    if isfield(arc,'nAxis') && ~isempty(arc.nAxis)
                        nAx = max(1, round(arc.nAxis));
                    end
                    nE = floor(m / nAx);
                    if nAx > 1 && nE >= 2
                        [om, sn] = obj.lagsPooledOverAxes(Abar, W, UiZ, ebar, ...
                                                          nAx, nE, nLag, ...
                                                          Bbar, LlS);
                    else
                        rr = lsqvce.redundancyLags(Abar, W, UiZ, ebar, ...
                                                   min(nLag, m), ...
                                                   'LocalBbar', Bbar, ...
                                                   'LocalLl', LlS);
                        om = rr.omega;  sn = rr.s;
                    end
                    if isempty(om)
                        % Every guard above passed but the axis pooling produced
                        % nothing. Without its own counter this shows up as an
                        % all-zero skip tally, which is exactly how a regression
                        % here hid behind a green test once already.
                        skip.axisPool = skip.axisPool + 1;
                        continue
                    end
                    if isnan(skip.dt) && isfield(arc,'t') && numel(arc.t) > 1
                        % Sampling interval in EPOCHS, not rows: with nAx axes
                        % per epoch the row-to-row time step is zero within an
                        % epoch, so take the step between epoch starts.
                        tt = arc.t(:);
                        if nAx > 1 && numel(tt) >= 2*nAx
                            skip.dt = median(diff(tt(1:nAx:end)));
                        else
                            skip.dt = median(diff(tt));
                        end
                    end
                    nArc = nArc + 1;
                    omegaT(:, nArc) = padLags(om, nLag);   %#ok<AGROW>
                    sT(:, nArc)     = padLags(sn, nLag);   %#ok<AGROW>
                end
            end
        end

        function [ovr, ordr, white, pwhite] = applyShape(obj, ovr, streams, res, v)
        %APPLYSHAPE  Re-estimate the AR temporal law per stream from residuals.
        %   Two paths for the descaled ACV:
        %     (a) a loader that exposes an epochCov hook (full per-epoch D) ->
        %         LSQVCE.RESIDUALACV (exact separable de-scale);
        %     (b) the scalar-sigma case (type 'diag') -> pool the RAW residuals
        %         with LSQVCE.COLLECTRESIDUALS, descale by the current scalar
        %         sigma, and take their biased ACV. Both paths estimate the TOTAL
        %         temporal law each pass, so the new AR REPLACES the old (no
        %         cross-pass filter composition).
        %   Also returns PWHITE(k): the Ljung-Box whiteness p-value per stream
        %   (NaN where the stream is trivially white / has no usable ACV), for
        %   the convergence diagnostic. A stream is declared white when p exceeds
        %   the WhiteAlpha threshold.
            n = numel(streams);
            ordr = zeros(1,n);  white = false(1,n);  pwhite = nan(1,n);
            dayFcnMap = obj.epochCovHooks();
            for k = 1:n
                [satKey, ot] = obj.splitStream(streams{k});
                nAx = obj.nAxesOf(ot);
                maxLag = obj.shapeOpt.MaxLag;
                if isempty(maxLag), maxLag = obj.shapeOpt.MaxOrder + 20; end
                acv = [];  nSamp = [];
                if isKey(dayFcnMap, streams{k})            % path (a)
                    try
                        [acv, nSamp] = lsqvce.residualAcv(res, ...
                                  obj.sysIdOfKey(satKey), ot, nAx, maxLag, ...
                                  dayFcnMap(streams{k}));
                    catch ME
                        if lsqvce.util.verbosity.detail(v)
                            lsqvce.util.log.msg('%s%s: epochCov acv failed (%s)\n', ...
                                lsqvce.util.verbosity.banner('Stochastic model'), ...
                                streams{k}, ME.identifier);
                        end
                    end
                end
                if isempty(acv)                            % path (b): scalar D
                    % Estimate the TOTAL temporal law from RAW residuals descaled
                    % by the current scalar sigma, so the new AR REPLACES the old
                    % one (no cross-pass filter composition). For multi-axis
                    % streams the residuals are epoch-major; pool the biased ACV
                    % PER AXIS down epochs (the axes share one temporal law).
                    vraw = lsqvce.collectResiduals(res, obj.sysIdOfKey(satKey), ot, 'v');
                    sig  = obj.currentSigma(ovr, satKey, ot);
                    [acv, nSamp] = biasedAcvAxes(vraw / sig, ...
                               min(maxLag, max(1, floor(numel(vraw)/max(nAx,1))-1)), nAx);
                end
                if isempty(acv) || numel(acv) < 2 || all(acv(2:end) == 0)
                    ordr(k) = 0;  white(k) = true;
                    ovr = setDiag(ovr, satKey, ot, obj.currentSigma(ovr, satKey, ot));
                    continue
                end
                if isempty(nSamp) || nSamp < 2, nSamp = numel(acv); end
                so = struct('pMax', obj.shapeOpt.MaxOrder, 'AcvTol', obj.shapeOpt.AcvTol);
                [a, s2, p, sel] = lsqvce.selectArAcv(acv, nSamp, so);
                ordr(k)  = p;

                % Whiteness gate: a LJUNG-BOX portmanteau test on the residual
                % ACF. The single-lag 2/sqrt(n) band is wrong here -- with ~100
                % lags the max sample autocorrelation of TRUE white noise exceeds
                % 2 sigma routinely (multiple comparisons). Ljung-Box tests all
                % lags jointly: p > alpha means "cannot reject whiteness".
                [lbWhite, lbP] = ljungBoxWhiteP(acv, nSamp, obj.shapeOpt.WhiteAlpha);
                pwhite(k) = lbP;                          % record for the diagnostic
                isWhite = (p == 0) || strcmp(sel.verdict, 'white') || lbWhite;
                white(k) = isWhite;
                if isWhite, p = 0; ordr(k) = 0; end        % treat as white downstream

                if ~isfield(ovr, satKey), ovr.(satKey) = struct(); end
                if ~isfield(ovr.(satKey), ot), ovr.(satKey).(ot) = struct(); end
                sigCur = obj.currentSigma(ovr, satKey, ot);

                if nAx > 1
                    % 3-D (or d-D) stream: SEPARABLE model C = D (R_t kron I_d) D.
                    % The per-epoch D is supplied by the loader (D.Cep / epochCov);
                    % here we set ONLY the shared scalar temporal AR shape and
                    % preserve the current scale -- VCE owns the scale, so the stochastic model
                    % must be scale-neutral or the two steps fight and neither
                    % converges. (LSQVCE.SEPARABLEWHITEN applies the per-axis
                    % filter AFTER de-scaling, so axes never mix.)
                    ovr.(satKey).(ot).type  = 'separable';
                    ovr.(satKey).(ot).d     = nAx;
                    if p == 0
                        ovr.(satKey).(ot).ar = 1;             % white in time
                    else
                        ovr.(satKey).(ot).ar = a;
                    end
                    ovr.(satKey).(ot).sigma = sigCur;         % scale unchanged
                    if isfield(ovr.(satKey).(ot), 'ma')       % drop stale field
                        ovr.(satKey).(ot) = rmfield(ovr.(satKey).(ot), 'ma');
                    end
                elseif p == 0
                    ovr = setDiag(ovr, satKey, ot, sigCur);   % scalar, white
                else
                    ovr.(satKey).(ot).type  = 'arma';         % scalar, coloured
                    ovr.(satKey).(ot).ar    = a;
                    ovr.(satKey).(ot).ma    = 1;
                    ovr.(satKey).(ot).sigma = sigCur;         % scale unchanged
                end
                if lsqvce.util.verbosity.detail(v)
                    if k == 1
                        % Title once, then one row per stream -- the same shape
                        % as Contribution and Convergence, so the three blocks
                        % of a pass read as one document.
                        lsqvce.util.log.msg('%s', lsqvce.util.verbosity.heading('Stochastic model'));
                    end
                    lsqvce.util.log.msg('%s', lsqvce.util.verbosity.streamLine(ot, ...
                        obj.sysLabelOf(satKey), ...
                        sprintf('AR(%d)   %s%s', p, sel.verdict, ...
                                ternary(nAx>1, '   (separable)', ''))));
                end
            end
        end

        function noise = noiseSummary(obj, ovr, streams)
        %NOISESUMMARY  Compact per-stream final noise model, for R.noise.
            noise = struct('stream', {streams}, 'model', {cell(1,numel(streams))});
            for k = 1:numel(streams)
                [satKey, ot] = obj.splitStream(streams{k});
                if isfield(ovr, satKey) && isfield(ovr.(satKey), ot)
                    noise.model{k} = ovr.(satKey).(ot);
                else
                    noise.model{k} = struct('type','diag','sigma',1);
                end
            end
        end

        % ---- small helpers ----------------------------------------------
        function key = vceGroupKey(obj)
        %VCEGROUPKEY  Translate 'satObs' to a redundancy GroupBy handle.
            if ischar(obj.vceGroups) || isstring(obj.vceGroups)
                switch lower(char(obj.vceGroups))
                    case 'satobs'
                        key = @(sysID,obsType,arc,day) ...
                              sprintf('%s|%s', string(sysID), string(obsType));
                    otherwise
                        key = char(obj.vceGroups);   % 'obsType'|'system'|'arc'|'day'
                end
            else
                key = obj.vceGroups;                 % function handle
            end
        end

        function s = currentSigma(~, ovr, satKey, ot)
        %CURRENTSIGMA  The scale currently in the override for a stream (else 1).
            s = 1;
            if isfield(ovr, satKey) && isfield(ovr.(satKey), ot) && ...
               isfield(ovr.(satKey).(ot), 'sigma')
                s = ovr.(satKey).(ot).sigma;
            end
            if ~(isscalar(s) && isfinite(s) && s > 0), s = 1; end
        end

        function n = nAxesOf(obj, ot)
        %NAXESOF  Axes per epoch for an obsType. Explicit axes() map wins; else
        %   any obsType whose name starts with 'pos' is treated as 3-D (x,y,z
        %   share one temporal law), everything else scalar. This makes merged
        %   pair-systems with 'posC'/'posD' work without an explicit axes() call.
            if isKey(obj.axesMap, ot)
                n = obj.axesMap(ot);
            elseif numel(ot) >= 3 && strcmpi(ot(1:3), 'pos')
                n = 3;
            else
                n = 1;
            end
        end

        function m = plotHooks(obj, stage)
        %PLOTHOOKS  Map stream 'satKey|obsType' -> plot handle for a STAGE
        %   ('postSolve' | 'postModel'), collected from every loader's
        %   plotPostSolve / plotPostModel obsType->handle map. Mirrors
        %   epochCovHooks. Absent hooks are simply omitted.
            m = containers.Map('KeyType','char','ValueType','any');
            field = ['plot', upper(stage(1)), stage(2:end)];   % plotPostSolve/Model
            L = obj.loaderList();
            for i = 1:numel(L)
                if ~isfield(L{i}, field) || isempty(L{i}.(field)), continue, end
                hookMap = L{i}.(field);
                if ~isa(hookMap, 'containers.Map') || hookMap.Count == 0, continue, end
                satKey = obj.satKeyOf(L{i}.sysID);
                for k = 1:numel(L{i}.obsTypes)
                    ot = L{i}.obsTypes{k};
                    if isKey(hookMap, ot)
                        m(sprintf('%s|%s', satKey, ot)) = hookMap(ot);
                    end
                end
            end
        end

        function firePlots(obj, stage, streams, R, ovr, vc, pass)
        %FIREPLOTS  Call each stream's plot handle for a stage, then show/export
        %   the returned figure. A plot hook must NEVER break the solve, so every
        %   call is wrapped: a throwing or figure-less hook is logged and skipped.
            hooks = obj.plotHooks(stage);
            if hooks.Count == 0, return, end
            ss  = obj.plotOpt.ShowSeconds;
            dir = obj.plotOpt.ExportDir;
            for k = 1:numel(streams)
                key = streams{k};
                if ~isKey(hooks, key), continue, end
                [satKey, ot] = obj.splitStream(key);
                ctx = obj.plotContext(stage, satKey, ot, key, R, ovr, vc, pass);
                hookFcn = hooks(key);                 % handle from the map
                try
                    % A hook may BUILD and RETURN a figure, or it may only do a
                    % side effect (record/log) and return nothing. Request one
                    % output but tolerate a no-output hook: MATLAB raises
                    % TooManyOutputs/maxlhs for a handle that assigns no outputs,
                    % in which case we call it for its side effect and treat the
                    % figure as empty (showExport is then a no-op).
                    try
                        figH = hookFcn(ctx);          % user builds the figure
                    catch MEout
                        if any(strcmp(MEout.identifier, {'MATLAB:TooManyOutputs', ...
                                'MATLAB:maxlhs', 'MATLAB:unassignedOutputs'}))
                            hookFcn(ctx);             % side-effect-only hook
                            figH = [];
                        else
                            rethrow(MEout);
                        end
                    end
                catch ME
                    lsqvce.util.log.warn('lsqvce:adjustVCE:plotHookFailed', ...
                        'Plot hook %s for %s threw (%s); skipping.', ...
                        stage, key, ME.identifier);
                    continue
                end
                name = sprintf('%s_%s_pass%d', stage, matlab.lang.makeValidName(key), pass);
                lsqvce.util.showExport(figH, name, ss, dir);
            end
        end

        function ctx = plotContext(obj, stage, satKey, ot, key, R, ovr, vc, pass)
        %PLOTCONTEXT  The single struct handed to a plot hook. FIXED, DOCUMENTED
        %   SCHEMA (ctx.schema): fields are only ADDED across versions, never
        %   renamed or removed, so a hook written against schema 1 keeps working.
        %
        %   ALWAYS PRESENT (schema 1)
        %     .schema      1x1     schema version (currently 1)
        %     .stage       char    'postSolve' | 'postModel'
        %     .stream      char    'satKey|obsType' (this stream's key)
        %     .satKey      char    the satellite/system key
        %     .obsType     char    the observation type
        %     .sysID       1x1     numeric/char system id
        %     .nAxes       1x1     axes per epoch (3 for positions, else 1)
        %     .pass        1x1     outer-pass index
        %     .R           1x1     lsqvce.Results (x, sol.Cxx, acc, residuals...)
        %     .x           nc x 1  common-parameter ESTIMATES
        %     .sd          nc x 1  common-parameter STANDARD DEVIATIONS
        %                          (sqrt(diag(sol.Cxx)); [] if unavailable)
        %     .common      cell    common parametrisation labels from the loader
        %                          ({} if the loader supplied none)
        %     .procOpt     struct  the loader's processingOptions
        %     .model       struct  the fitted noise model for this stream
        %                          (.type/.ar/.sigma/.d ...)
        %     .sigmaG      1x1     the VCE sigma of this stream's group (NaN if n/a)
        %     .exportDir   char    where showExport writes PNGs ('' = no export)
        %     .showSeconds 1x1     how long showExport displays a figure
        %
        %   POSTMODEL ONLY (the stochastic model is final here)
        %     .psd         Kx1     estimated one-sided PSD of the fitted model
        %                          (S(w)=sigma^2|B|^2/|A|^2 over .psdFreq)
        %     .psdFreq     Kx1     normalised frequencies [0, pi] for .psd
        %     .acv         Lx1     model autocovariance (sigma^2 * armaACV)
        %     .arOrder     1x1     fitted AR order (numel(a)-1)
            ctx = struct();
            ctx.schema   = 1;
            ctx.stage    = stage;            % 'postSolve' | 'postModel'
            ctx.stream   = key;              % 'satKey|obsType'
            ctx.satKey   = satKey;
            ctx.obsType  = ot;
            ctx.sysID    = obj.sysIdOfKey(satKey);
            ctx.nAxes    = obj.nAxesOf(ot);
            ctx.pass     = pass;             % outer-pass index
            ctx.R        = R;                % lsqvce.Results (x, sol, acc, ...)
            ctx.x        = R.x;              % common-parameter estimate

            % standard deviations of the common parameters (sqrt of the diagonal
            % of the a-posteriori covariance). Guarded: a Results without a full
            % Cxx leaves this empty rather than erroring.
            ctx.sd = [];
            try
                if ~isempty(R.sol) && isfield(R.sol, 'Cxx') && ~isempty(R.sol.Cxx)
                    ctx.sd = sqrt(max(real(diag(R.sol.Cxx)), 0));
                end
            catch
            end

            % common parametrisation labels + processingOptions from the loader
            % that owns this system (piped straight through, untouched).
            [ctx.common, ctx.procOpt] = obj.loaderParametrisation(ctx.sysID);

            ctx.exportDir   = obj.plotOpt.ExportDir;
            ctx.showSeconds = obj.plotOpt.ShowSeconds;
            % the fitted noise model for this stream (AR/PSD + sigma), if any
            if isfield(ovr, satKey) && isfield(ovr.(satKey), ot)
                ctx.model = ovr.(satKey).(ot);
            else
                ctx.model = struct('type','diag','sigma',1);
            end
            % the VCE weight of this stream's group (post-model is where it matters)
            ctx.sigmaG = obj.streamSigmaG(key, vc);

            % --- postModel: the estimated PSD / ACV / AR order are now final ---
            if strcmp(stage, 'postModel')
                [ctx.psd, ctx.psdFreq, ctx.acv, ctx.arOrder] = obj.modelSpectrum(ctx.model);
            end
        end

        function [common, procOpt] = loaderParametrisation(obj, sysID)
        %LOADERPARAMETRISATION  The common labels and processingOptions of the
        %   loader owning system SYSID. Empdefaults if not found. Piped to hooks.
        %   OUTPUT common cell (labels or {});  procOpt struct (or struct())
            common = {};  procOpt = struct();
            L = obj.systems;
            if ~iscell(L), L = {L}; end
            for i = 1:numel(L)
                Li = L{i};
                if isstruct(Li) && isfield(Li, 'systemId') && isequaln(Li.systemId, sysID)
                    if isfield(Li, 'common'),  common  = Li.common;  end
                    if isfield(Li, 'procOpt'), procOpt = Li.procOpt; end
                    return
                end
            end
        end

        function [S, f, r, arOrder] = modelSpectrum(obj, model)
        %MODELSPECTRUM  One-sided PSD, ACV, and AR order of a fitted noise model.
        %   S(w) = sigma^2 |B(e^jw)|^2 / |A(e^jw)|^2 on f in [0, pi]. For a diag
        %   (white) model the PSD is flat at sigma^2. Cheap; postModel only.
        %   INPUT  model struct with .type and (.ar/.ma/.sigma) or .sigma
        %   OUTPUT S Kx1 PSD; f Kx1 freqs; r Lx1 model ACV; arOrder 1x1
            S = [];  f = [];  r = [];  arOrder = 0;                      %#ok<NASGU>
            K = 256;
            f = linspace(0, pi, K).';
            sigma = 1;
            if isfield(model, 'sigma') && ~isempty(model.sigma), sigma = model.sigma(1); end
            a = 1;  b = 1;
            if isfield(model, 'ar') && ~isempty(model.ar), a = model.ar(:).'; end
            if isfield(model, 'ma') && ~isempty(model.ma), b = model.ma(:).'; end
            arOrder = numel(a) - 1;
            % frequency response of A and B at e^{jw}
            z = exp(-1i * f);                          % Kx1
            Aw = polyval(fliplr(a), z);                % A(e^{-jw})
            Bw = polyval(fliplr(b), z);
            Aw(Aw == 0) = eps;
            S = (sigma^2) * abs(Bw).^2 ./ abs(Aw).^2;  % one-sided PSD shape
            % model autocovariance (physical scale)
            try
                r = sigma^2 * lsqvce.armaACV(a, b, 64);
                r = r(:);
            catch
                r = [];
            end
        end

        function s = streamSigmaG(obj, key, vc)
        %STREAMSIGMAG  The VCE sigma of the group this stream belongs to (or NaN).
            s = NaN;
            if isempty(vc) || ~isfield(vc,'groups') || isempty(vc.groups), return, end
            gk = obj.vceGroupKeyList({key}, obj.vceGroupKeyFcn());
            vcKeys = cellstr(string(vc.groups));
            gi = find(strcmp(vcKeys, gk{1}), 1);
            if ~isempty(gi) && gi <= numel(vc.sigma), s = vc.sigma(gi); end
        end

        function m = epochCovHooks(obj)
        %EPOCHCOVHOOKS  Map stream 'satKey|obsType' -> dayFcn(sat,day), for
        %   MULTI-AXIS streams whose loader carries an epochCov hook. Only
        %   multi-axis streams are registered: LSQVCE.RESIDUALACV needs
        %   .Cep for those, whereas scalar streams use the path-(b) residual ACV
        %   (their descale is a single sigma). Absent hooks are simply omitted.
            m = containers.Map('KeyType','char','ValueType','any');
            L = obj.loaderList();
            for i = 1:numel(L)
                if ~isfield(L{i}, 'epochCov') || isempty(L{i}.epochCov), continue, end
                satKey = obj.satKeyOf(L{i}.sysID);
                for k = 1:numel(L{i}.obsTypes)
                    ot  = L{i}.obsTypes{k};
                    if obj.nAxesOf(ot) <= 1, continue, end   % scalar -> path (b)
                    key = sprintf('%s|%s', satKey, ot);
                    m(key) = makeDayFcn(L{i}, ot);
                end
            end
        end

        function k = satKeyOf(~, sysID)
        %SATKEYOF  A struct-field-safe key for a satellite id.
            k = matlab.lang.makeValidName(char(string(sysID)));
        end
        function lab = sysLabelOf(obj, satKey)
        %SYSLABELOF  A READABLE name for a satKey, for output only.
        %   satKey is a struct-field-safe MANGLING of the system id, so a
        %   NUMERIC id comes back prefixed: makeValidName('1') = 'x1'. That is
        %   correct internally (a field name cannot start with a digit) but it
        %   is meaningless on screen -- users asked what "x1" was, and the honest
        %   answer is "system 1, spelled so MATLAB can use it as a field name".
        %
        %   So for DISPLAY we undo it: prefer the loader's own 'name' if one was
        %   given to LSQVCE.MAKELOADER (e.g. 'demo-CD-pair'), else the ACTUAL
        %   sysID, and only fall back to the mangled key if neither is available.
        %
        %   INPUT   satKey 1xN char   the internal key
        %   OUTPUT  lab    1xN char   what to print
            lab = char(string(obj.sysIdOfKey(satKey)));
            L = obj.loaderList();
            for i = 1:numel(L)
                if ~strcmp(obj.satKeyOf(L{i}.sysID), satKey), continue, end
                if isstruct(L{i}) && isfield(L{i}, 'name') && ~isempty(L{i}.name)
                    lab = char(L{i}.name);
                end
                return
            end
        end

        function id = sysIdOfKey(obj, satKey)
        %SYSIDOFKEY  The ACTUAL sysID value (numeric or char) for a satKey.
        %   Reconstructing the id from the key string is unreliable for numeric
        %   ids (makeValidName('1')='x1'); instead we look it up against the
        %   loaders, so residual matching (isequal on the real sysID) is exact.
            L = obj.loaderList();
            id = satKey;                       % fallback: the key itself
            for i = 1:numel(L)
                if strcmp(obj.satKeyOf(L{i}.sysID), satKey)
                    id = L{i}.sysID;  return
                end
            end
        end
        function [satKey, ot] = splitStream(~, s)
            parts = split(string(s), '|');
            satKey = char(parts(1));  ot = char(parts(2));
        end
    end
end

% ==========================================================================
%  FILE-LOCAL HELPERS
% ==========================================================================
function s = mergeNV(s, nv)
%MERGENV  Merge Name/Value pairs into struct S, matching field names
%   CASE-INSENSITIVELY and canonicalises to the stored UpperCamelCase field
%   (WP-39 round a2, PINNED DEVIATION: the fluent configurators KEEP this
%   merge -- its case-insensitive-to-canonical matching is exactly MATLAB's
%   own name-value semantics, and the sentinel machinery an arguments
%   migration would need adds nothing a caller can observe. Canonical
%   casing is UpperCamelCase; typed validation stays with the downstream
%   consumers). Unknown names are added verbatim (validated downstream).
    fn = fieldnames(s);
    for i = 1:2:numel(nv)
        key = nv{i};
        hit = fn(strcmpi(fn, key));           % existing field, any case
        if ~isempty(hit), key = hit{1}; end   % canonicalise to stored case
        s.(key) = nv{i+1};
    end
end

function out = ternary(c, a, b)
%TERNARY  Inline conditional: return a if c is truthy, else b.
%   INPUT c 1x1 logical; a,b any   OUTPUT out  a if c else b
    if c, out = a; else, out = b; end
end

function [r, nSamp] = biasedAcv(x, maxLag)
%BIASEDACV  Biased (1/n) autocovariance r(0..maxLag) of a column vector.
%   Also returns the sample count, so callers can build an honest whiteness
%   band. Biased so it stays PSD and LEVINSONDURBIN never trips.
    x = x(:);
    x = x(isfinite(x));
    n = numel(x);
    nSamp = n;
    if n < 2, r = []; return, end
    maxLag = min(maxLag, n-1);
    x = x - mean(x);
    r = zeros(maxLag+1, 1);
    for h = 0:maxLag
        r(h+1) = sum(x(1:n-h) .* x(1+h:n)) / n;   % biased
    end
end

function [r, nSamp] = biasedAcvAxes(x, maxLag, nAx)
%BIASEDACVAXES  Biased ACV pooled over axes for epoch-major data.
%   x is [ax1 ax2 .. axd] for epoch 1, then epoch 2, ... The d axes share one
%   temporal law (the separable model), so estimate the ACV down epochs for
%   each axis and average. Also returns the pooled sample count. For nAx == 1
%   this reduces to BIASEDACV.
    if nargin < 3 || isempty(nAx) || nAx <= 1
        [r, nSamp] = biasedAcv(x, maxLag);  return
    end
    x = x(:);
    nE = floor(numel(x) / nAx);
    if nE < 2, r = []; nSamp = 0; return, end
    X = reshape(x(1:nE*nAx), nAx, nE).';          % epochs x axes
    maxLag = min(maxLag, nE-1);
    acc = zeros(maxLag+1, 1);  cnt = 0;
    for ax = 1:nAx
        ri = biasedAcv(X(:,ax), maxLag);
        if ~isempty(ri)
            if numel(ri) < maxLag+1, ri(end+1:maxLag+1) = 0; end
            acc = acc + ri(1:maxLag+1);  cnt = cnt + 1;
        end
    end
    if cnt == 0, r = []; nSamp = 0; else, r = acc / cnt; nSamp = nE; end
end

function s = setField(s, a, b, f, val)
%SETFIELD  s.(a).(b).(f) = val, creating the nesting as needed.
    if ~isfield(s, a),       s.(a)      = struct(); end
    if ~isfield(s.(a), b),   s.(a).(b)  = struct(); end
    s.(a).(b).(f) = val;
end

function s = setDiag(s, a, b, sigma)
%SETDIAG  Set s.(a).(b) to a clean 'diag' spec (drops any stale ar/ma fields).
    if ~isfield(s, a), s.(a) = struct(); end
    s.(a).(b) = struct('type', 'diag', 'sigma', sigma);
end

function tf = ljungBoxWhite(acv, nSamp, alpha)
%LJUNGBOXWHITE  Ljung-Box portmanteau whiteness test from an ACV vector.
%   tf = true means "cannot reject whiteness at level alpha" (i.e. treat the
%   stream as temporally white). Q = n(n+2) sum_{k=1..h} r_k^2/(n-k) ~ chi2(h)
%   under the white null; r_k = acv(k+1)/acv(1). Uses the same estimator as
%   LSQVCE.SYSTEMATIC. h = min(nLag, floor(n/5)).
    if nargin < 3 || isempty(alpha), alpha = 0.01; end
    acv = acv(:);
    n   = round(nSamp);
    if numel(acv) < 2 || acv(1) <= 0 || n < 10, tf = true; return, end
    h = min(numel(acv) - 1, floor(n/5));
    if h < 1, tf = true; return, end
    r = acv(2:h+1) / acv(1);
    Q = 0;
    for k = 1:h
        Q = Q + r(k)^2 / (n - k);
    end
    Q = n * (n + 2) * Q;
    p = 1 - lsqvce.stat.chi2Cdf(Q, h);
    tf = p > alpha;              % white unless we can reject at level alpha
end

function fcn = makeDayFcn(loader, obsType)
%MAKEDAYFCN  Wrap a loader's epochCov hook into the dayFcn residualAcv
%   wants: dayFcn(satId,day) -> struct with .arc(i).rows/.epochs and
%   .sigma (nAxes==1) or .Cep (nEpochs x nAxes x nAxes).
    fcn = @(satId, day) loader.epochCov(day, obsType);
end

function w = pluckResidual(res, ot, sysID)
%PLUCKRESIDUAL  Whitened residual vector for one (satellite,obsType) stream.
%   Thin wrapper over LSQVCE.COLLECTRESIDUALS (the same accessor applyShape
%   uses), so validate() reads residuals exactly as the loop does. The match is
%   on the real sysID (numeric or char), never the struct-field key.
    w = [];
    try
        w = lsqvce.collectResiduals(res, sysID, ot, 'w');   % whitened
    catch
        try
            w = lsqvce.collectResiduals(res, sysID, ot, 'v'); % fallback: raw
        catch
            w = [];
        end
    end
    w = w(:);
end

function [tf, p] = ljungBoxWhiteP(acv, nSamp, alpha)
%LJUNGBOXWHITEP  Ljung-Box whiteness test returning BOTH verdict and p-value.
%   Identical statistic to the file-local LJUNGBOXWHITE (which returns only the
%   verdict); validate() needs the p-value for its per-stream summary.
%   Q = n(n+2) sum_{k=1..h} r_k^2/(n-k) ~ chi2(h) under the white null.
    if nargin < 3 || isempty(alpha), alpha = 0.01; end
    acv = acv(:);
    n   = round(nSamp);
    if numel(acv) < 2 || acv(1) <= 0 || n < 10, tf = true; p = NaN; return, end
    h = min(numel(acv) - 1, floor(n/5));
    if h < 1, tf = true; p = NaN; return, end
    r = acv(2:h+1) / acv(1);
    Q = 0;
    for k = 1:h
        Q = Q + r(k)^2 / (n - k);
    end
    Q = n * (n + 2) * Q;
    p  = 1 - lsqvce.stat.chi2Cdf(Q, h);
    tf = p > alpha;
end

function share = contributionShare(r)
%CONTRIBUTIONSHARE  Level-1 relative information share from per-group redundancy.
%   share = tr(N^-1 N_g^cc) / (m-u) = r_g / sum_g r_g, which sums to 1 over the
%   groups. This is the SCALAR contribution of each system/type to the solution
%   (Ellmer 2018; contribution matrix C_g = N^-1 N_g^cc, share_g = tr(C_g)/u_red).
%
%   INPUT
%     r      1xnGroup double   per-group redundancy r_g = tr(N^-1 N_g^cc), as
%                              returned in vc.r by LSQVCE.VCE (aligned with
%                              vc.groups). May contain NaN for a group left
%                              below MinRedundancy.
%   OUTPUT
%     share  1xnGroup double   normalised share (nonnegative, sums to 1 over the
%                              finite entries). All-zero / all-NaN input yields an
%                              all-NaN row (no redundancy => share undefined),
%                              so it never poisons the trajectory with a silent 0.
%
%   Created 2026-07-17 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
    r = r(:).';
    share = nan(1, numel(r));
    fin = isfinite(r);
    s = sum(r(fin));
    if s > 0
        share(fin) = r(fin) / s;   % NaN groups stay NaN; finite groups normalised
    end
end

function str = contribLine(groups, share)
%CONTRIBLINE  Format 'key=xx.x%' pairs for the [CONTRIB] verbose banner.
%   INPUT
%     groups 1xnGroup   VCE group keys aligned with share; accepts a cellstr,
%                       string array, numeric vector, or a cell of scalars.
%     share  1xnGroup double            relative info share in [0,1] (or NaN)
%   OUTPUT
%     str    1xN char                   e.g. 'A=61.2%  B=38.8%' (NaN -> 'n/a')
%
%   Created 2026-07-17 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
    keys  = groupKeyStrings(groups);
    parts = cell(1, numel(share));
    for k = 1:numel(share)
        if isfinite(share(k))
            parts{k} = sprintf('%s=%.1f%%', keys{k}, 100*share(k));
        else
            parts{k} = sprintf('%s=n/a', keys{k});
        end
    end
    str = strjoin(parts, '  ');
end

function str = convergenceRow(satKey, obsType, sigma, isWhite, pval)
%CONVERGENCEROW  One stream's row of the level-3 convergence table.
%   '  -> GraceA "kinematicPosition"   sigma = 1.0342   white = yes  (p = 0.31)'
%
%   Both halves of the outer loop's stopping rule are shown on ONE row: the variance factor
%   (sigma, driven to 1 by the VCE) and whiteness of the residuals after
%   AR whitening). The gate stops only when BOTH hold for every stream, so a row
%   that fails one but not the other says immediately which half is stuck --
%   which the scalar max|sigma-1| cannot.
%
%   INPUT
%     satKey  1xN char     satellite/system key of the stream
%     obsType 1xN char     observation type, quoted in the output
%     sigma   1x1 double   the scale currently applied to this stream (1 = right)
%     isWhite 1x1 logical  did the residuals pass the Ljung-Box whiteness test
%     pval    1x1 double   that test's p-value ([] or NaN prints as n/a)
%   OUTPUT
%     str     1xN char     one row, newline-terminated
%
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
    if isWhite, wtxt = 'yes'; else, wtxt = 'NO '; end   % NO padded, columns align
    if isempty(pval) || ~isfinite(pval)
        ptxt = '(p = n/a)';
    else
        ptxt = sprintf('(p = %.3f)', pval);
    end
    % The SAME row shape as Stochastic model and Contribution. These rows used
    % to be indented by the 21-character stage-label column while the others
    % started at 2, so the three tables of a single pass stepped in and out.
    str = lsqvce.util.verbosity.streamLine(obsType, satKey, ...
              sprintf('sigma = %7.4f   white = %s  %s', sigma, wtxt, ptxt));
end

function str = sigmaLine(groups, sigma)
%SIGMALINE  Format 'key=sigma' pairs for the diagnostic [GATE] banner.
%   The gate reports max|sigma-1|, which names HOW FAR the worst group is from
%   unity but not WHICH group it is -- so a loop stalled on one stubborn stream
%   is indistinguishable from one drifting everywhere. This prints the per-group
%   sigmas behind that summary, which is the first thing wanted when a run will
%   not converge (see HISTORY_v1.md 1, "PER-STREAM SIGMAS in the convergence report").
%
%   INPUT
%     groups 1xnGroup   VCE group keys aligned with sigma; accepts a cellstr,
%                       string array, numeric vector, or a cell of scalars.
%     sigma  1xnGroup double   per-group variance factor (1 = correctly weighted)
%   OUTPUT
%     str    1xN char   e.g. 'posC=1.003  posD=0.982  rr=1.041' (NaN -> 'n/a')
%
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
    keys  = groupKeyStrings(groups);
    parts = cell(1, numel(sigma));
    for k = 1:numel(sigma)
        if k <= numel(keys), key = keys{k}; else, key = sprintf('g%d', k); end
        if isfinite(sigma(k))
            parts{k} = sprintf('%s=%.3f', key, sigma(k));
        else
            parts{k} = sprintf('%s=n/a', key);
        end
    end
    str = strjoin(parts, '  ');
end

function keys = groupKeyStrings(groups)
%GROUPKEYSTRINGS  Normalise any VCE group-key form to a cellstr.
%   VCE group keys arrive as a cellstr, a string array, a numeric vector, or a
%   cell of mixed scalars depending on how GroupBy was specified, and every
%   formatter needs brace-indexable char keys. Shared by CONTRIBLINE and
%   SIGMALINE so the two banners can never disagree about how a key is spelled.
%
%   INPUT   groups 1xnGroup  cellstr | string | numeric | cell of scalars
%   OUTPUT  keys   1xnGroup cell of char
%
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
    if isnumeric(groups)
        keys = cellfun(@num2str, num2cell(groups), 'UniformOutput', false);
    elseif isstring(groups)
        keys = cellstr(groups);
    elseif iscell(groups)
        keys = cellfun(@(x) char(string(x)), groups, 'UniformOutput', false);
    else
        keys = cellstr(string(groups));
    end
end

function y = padLags(v, n)
%PADLAGS  Pad a per-lag vector with zeros to length n.
%   Ellmer eqs. (6.5.22)-(6.5.23) define the contribution of a lag to be ZERO
%   when the lag exceeds the arc length, so short arcs simply contribute nothing
%   at the long lags rather than being excluded.
%   INPUT   v  kx1 double   the per-lag values actually computed
%           n  1x1 int      target length
%   OUTPUT  y  nx1 double   v padded (or truncated) to n, non-finite -> 0
    y = zeros(n, 1);
    k = min(numel(v), n);
    y(1:k) = v(1:k);
    y(~isfinite(y)) = 0;
end

function s = firstErrText(lastErr)
%FIRSTERRTEXT  Render the first captured Wfcn error for the skip warning.
%   INPUT   lastErr char   '' if nothing was captured
%   OUTPUT  s       char   ' First error: <text>.'  or  '' when empty
    if isempty(lastErr)
        s = '';
    else
        s = sprintf(' First error: %s.', lastErr);
    end
end

function fh = makeAcfFromAcv(c)
%MAKEACFFROMACV  Wrap a discrete autocovariance vector as an acf(lag) handle.
%   covModel's 'acf' type calls acf(t(idx) - t(idx).') with a LAG MATRIX, so the
%   discrete sequence c(0..L-1) is looked up by rounded |lag| and zero beyond the
%   last modelled lag (the model simply says "no correlation out there").
%   INPUT   c  Lx1 double   autocovariance, c(1) = variance at lag 0
%   OUTPUT  fh function_handle  @(lagMatrix, opt) -> covariance matrix
    c = c(:);
    L = numel(c);
    fh = @(lag, ~) acvLookup(c, L, lag);
end

function C = acvLookup(c, L, lag)
%ACVLOOKUP  c(|lag|) with zeros beyond the modelled range.
    k = round(abs(lag)) + 1;            % lag 0 -> index 1
    C = zeros(size(k));
    in = k >= 1 & k <= L;
    C(in) = c(k(in));
end
