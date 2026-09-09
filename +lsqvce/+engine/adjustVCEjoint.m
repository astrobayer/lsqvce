function out = adjustVCEjoint(accFcn, opt)
%LSQVCE.ENGINE.ADJUSTVCEJOINT  Joint (co-estimated) stochastic-model loop.
%
%   out = LSQVCE.ENGINE.ADJUSTVCEJOINT(accFcn)
%   out = LSQVCE.ENGINE.ADJUSTVCEJOINT(accFcn, 'Name', value, ...)
%
%   Runs Ellmer's (2018, §6.5.7) ten-step algorithm, in which ALL variance
%   components are updated TOGETHER from one residual pass and one shared probe:
%   the per-frequency PSD amplitudes Shat^j (the temporal law), the arc-wise
%   variance factors sigma2_m (the within-group variance factors), and -- when a prior is
%   present -- the geophysical factor sigma2_gpm.
%
%   JOINT vs SEPARATE -- WHY THIS IS ITS OWN DRIVER
%     LSQVCE.ENGINE.ADJUSTVCEJOINT alternates: fit the temporal whitening, then the
%     between-group variance components, each with its own criterion. That is the SEPARABLE model,
%     and it is the toolbox DEFAULT because the whiteness gate gives a clean
%     stopping rule and the alternation is robust.
%     Ellmer instead treats every PSD frequency and every arc factor as variance
%     components of ONE Helmert VCE and updates them SIMULTANEOUSLY (step 9).
%     This matters: alternating (block-coordinate) and joint updates share fixed
%     points only if the blocks are SEPARABLE, and here they are NOT -- the arc
%     factor and the frequency amplitudes both scale the same Sigma_ll (eqs.
%     6.5.7, 6.5.18), which is exactly why the gamma normalisation (eq. 6.5.28)
%     exists to prevent double counting. The two schemes can therefore converge
%     to DIFFERENT Sigma. Use 'joint' for Ellmer fidelity, 'separate' (the
%     default elsewhere) for the general separable model.
%
%   THE TEN STEPS (Ellmer §6.5.7), and where each lives
%     1-3  build Delta_l, A; form Sigma from the current model; decorrelate
%            -> the caller's accFcn does this (see INPUT below)
%     4-5  accumulate N, n; factor N = U'U and solve
%            -> acc.solve() and chol(acc.reducedNormalMatrix())
%     6    draw Zbar and form U^-1 Zbar ONCE for all arcs
%            -> LSQVCE.ENGINE.SHAREDPROBE
%     7    decorrelated residuals ebar = W^-T (l - A x - B y), the FULL
%          Table-6.2 vector (WP-33; the local fit is subtracted per system)
%            -> re-streamed per arc via the store's dataFcn
%     8    Rtilde, etilde, and from them Omega and s per lag
%            -> LSQVCE.KERNEL.REDUNDANCYLAGS
%     9    updated Shat^j, sigma2_m (and sigma2_gpm)
%            -> LSQVCE.ENGINE.PSDVCE and LSQVCE.ENGINE.ARCVARIANCEFACTORS
%     10   if not converged, rebuild Sigma and go to step 2
%
%   *** START MATTERS (verified; read this) ***
%     Ellmer's literal white start (eq. 6.5.62) is a DEGENERATE fixed point of
%     his own weighting -- see the warning in LSQVCE.ENGINE.PSDVCE. This driver
%     therefore defaults 'C0' to a mildly correlated covariance rather than a
%     delta, and passes any degeneracy straight through to psdVCE so the warning
%     is raised rather than hidden.
%
%   INPUT
%     accFcn  function_handle  @(model) -> a SOLVED-READY lsqvce.engine.BlockAccumulator
%             built under the given stochastic model. It is called once per outer
%             iteration and must apply the model when whitening (steps 1-3).
%             `model` is a struct with fields:
%               .C       nLag x 1  the current covariance function C^n_xx
%               .sigma2  M x 1     the current arc variance factors
%               .iter    1x1       the iteration number (1 on the first call)
%             On the FIRST call the caller may ignore the model and use its own
%             a-priori weighting; the driver then estimates from there.
%
%   OPTIONS (Name/Value)
%     'NumLags'   1x1 int    Nmax, number of lags/frequencies (default 60)
%     'C0'        nLag x 1   starting covariance function. Default: a mildly
%                            correlated [1; 0.3; 0.3; ...] rather than a delta,
%                            because the delta is degenerate (see above).
%     'Sigma2_0'  M x 1 | [] starting arc factors (default all ones, eq. 6.5.62)
%     'MaxIter'   1x1 int    outer iterations (default 10)
%     'Tol'       1x1        convergence on max|alpha2 - 1| (default 1e-2)
%     'Probes'    1x1 int    Z Hutchinson probes (default 100, Ellmer's value)
%     'Seed'      1x1 | []   RNG seed for reproducible probes (default [])
%     'Norm'      char       DCT normalisation, 'ortho' (default) | 'ellmer'
%     'MeanType'  char       arc-factor mean, 'robust' (default) | 'plain'
%     'Influence' string     "none" (default) | "huber" | "igg3" | "bisquare"
%              R6 WP-105: the ROBUST PASS INSIDE THE LOOP. Every pass ends
%              with one residual sweep; its LAW-whitened residuals (the
%              current law, level, arc factor and Cep divided out) go through
%              LSQVCE.KERNEL.ROBUSTWEIGHTS and the factors ride inside every
%              arc's whitener on the NEXT pass (model.rowWeights, which the
%              caller's factory forwards as UserOpt.rowWeights). Accumulation,
%              probe redundancies and the component update then all see ONE
%              re-weighted problem, so the Helmert equations stay consistent:
%              a row of weight w contributes w to Omega AND w to the
%              redundancy. NOT Ellmer (2018): he screens outliers in
%              preprocessing and down-weights whole ARCS through the variance
%              factors; this is Koch's iterated re-weighting inside a VCE,
%              and the levels it returns are those of the RETAINED
%              observations (no consistency correction is applied -- the
%              bias is measured by the oracle, test_vcePolicy). The residuals
%              are standardised by ONE MAD SCALE PER STREAM, never by the
%              estimated level (perStreamWeights says why). Prior streams
%              are never re-weighted. Prefer "igg3" / "bisquare" for gross
%              blunders: Huber keeps c/|u| of their weight. "none" is bitwise
%              the plain estimator.
%     'Tuning'    [] | 1xK   tuning constants of the influence (its defaults)
%     'Verbose'   logical    per-iteration banner (default false)
%     'StepSchedule' 'none' (default) | 'decay'.  WP-40 r2 (F-15): the
%              closed-loop alternation is a stochastic approximation; the
%              plain step is BISTABLE for short arcs (recovered-PSD corr can
%              collapse to negative values on unlucky seeds). 'decay' uses
%              the Robbins-Monro step Gamma_it = 2/(2+it) via PSDVCE's
%              'Gamma'. HONEST CLAIM (psd2acv-faithful mirror, WP-40 r6):
%              damping is DAMAGE LIMITATION -- it prevents the undamped
%              collapse (corr can hit ~0.2) but the closed loop does not
%              beat OPEN-LOOP pass-1 quality on the tested fixtures;
%              production guidance is open loop or <= 3 damped passes.
%              'none' is bitwise the old behaviour.
%     'AvgTail' 1x1 integer >= 0, default 0.  Polyak averaging of the
%              RETURNED spectrum over the last AvgTail passes (log domain);
%              (honest r6 claim: helps the damped loop settle; cannot
%              rescue an undamped collapse, and longer tails average over
%              late-pass drift -- keep it <= 4) and does nothing useful
%              without it (mirror-measured). The iteration itself always
%              runs on the raw C.
%
%   OUTPUT
%     out  1x1 struct
%       .C         nLag x 1   the converged covariance function
%       .Shat      nLag x 1   the converged PSD amplitudes
%       .sigma2    M x 1      the converged arc variance factors, ALIGNED TO
%                             .arcKeys (the pass-1 enumeration order; WP-34)
%       .robust    struct     .influence, .scale (1 x S, the MAD scale per
%                             stream of the last robust sweep), .streams
%                             (1 x S 'satKey|obsType', sorted), .factor
%                             (n x 1, the final per-observation factors in
%                             sweep order), .keys (A x 4 {sysID, obsType,
%                             day, arc}); empty under "none" (R6 WP-105)
%       .arcKeys   M x 1 cell arc keys s<sys>_<obsType>_d<day>_a<arc> pinned
%                             on pass 1 (LSQVCE.ENGINE.REDUNDANCY's 'arc' format) --
%                             the ONE identity .sigma2 and Sigma2_0 align to
%       .skips     S x 3 cell {pass, key, reason} for every arc skipped in
%                             any pass ('gap' | 'nonPD'); a skipped arc KEEPS
%                             its previous factor
%       .acc       1x1        the final accumulator
%       .sol       1x1        the final solution struct
%       .converged 1x1 logical whether max|alpha2-1| fell below Tol
%       .nIter     1x1        outer iterations actually run
%       .history   struct     .maxDev (nIter x 1), .Shat (nIter x nLag)
%
%   EXAMPLE
%     accFcn = @(model) buildAccumulatorWithModel(cfg, model);
%     r = lsqvce.engine.adjustVCEjoint(accFcn, 'NumLags', 60, 'MaxIter', 8);
%     semilogy(r.Shat);          % the estimated noise PSD
%     plot(r.sigma2);            % the arc-wise weights
%
%   ERRORS
%     lsqvce:adjustVCEjoint:needCfcn  the joint mode needs a cross-covariance
%                                     function and none was supplied
%     lsqvce:adjustVCEjoint:badC0     the starting cross-covariance is malformed
%     lsqvce:adjustVCEjoint:notPD     the joint covariance is not positive
%                                     definite, so it cannot be factored
%     lsqvce:adjustVCEjoint:skipThreshold  a pass skipped more than
%                                   SkipThreshold of the enumerated arcs (R5 M1)
%     lsqvce:adjustVCEjoint:badAxisRows  a multi-axis arc whose row count is not a
%                                   multiple of its axis count (R6 WP-90: an arc
%                                   carrying a per-epoch Cep, scalar or multi-axis,
%                                   is descaled per epoch before its lag statistics;
%                                   the former noCep refusal is gone -- the loader
%                                   composes D = I and says so once)
%     lsqvce:adjustVCEjoint:noArcs    a system declares no arcs
%     lsqvce:adjustVCEjoint:collapsed the influence function left fewer live
%                                     observations than parameters + 1 (beyond
%                                     the breakdown point; R6 WP-105)
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:adjustVCEjoint:localNotPD           System X
%     lsqvce:adjustVCEjoint:notConverged         [warning] The joint stochastic model did not settle in X outer iterations (last...
%     lsqvce:adjustVCEjoint:arcFactorsDiverging  [warning] arc-factor spread grew (F-01)
%
%   WHAT eq. 6.5.28 MULTIPLIES, AND WHY IT MATTERS (R6 WP-98, finding F-01)
%     gamma2_m = alpha2_m * sigma2_m is correct only when sigma2_m is the factor
%     the WEIGHTING used. It is not the per-arc factor: LSQVCE.ENGINE.JOINTMODELOPT
%     can express only a per-STREAM level in the covariance override, so that
%     level is what the data saw and that level is what multiplies here.
%     Multiplying by the per-arc factor instead made the update a geometric
%     recursion with NO FIXED POINT, because alpha2_m cancels its own arc's
%     factor exactly -- sweepArcs scales omegaT and sT by the same s2_m -- and
%     so never responds to it. MEASURED on a fixture with planted factors
%     1 / 9 / 1: the ratio ran 9.4, 79, 671, 5626, 4.7e4 over five passes, a
%     constant 8.4x per pass; it now settles at 9.4, 8.6, 8.6, 8.6, 8.6.
%     For a stream with arcFactors = false the two numbers are identical (every
%     arc carries the one stream level), so those runs are unchanged.
%     STILL OPEN: only the RATIO between a stream's arcs is pinned by evidence.
%     The absolute split between a stream's LEVEL and its arc factors is fixed
%     by gammaBar (eq. 6.5.30), and on that fixture the level settles at 0.678
%     where the undegraded arcs were planted at 0.5. HANDOVER 12.22.
%
%   See also LSQVCE.ENGINE.PSDVCE, LSQVCE.ENGINE.ARCVARIANCEFACTORS,
%            LSQVCE.KERNEL.REDUNDANCYLAGS, LSQVCE.ENGINE.SHAREDPROBE
%
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-105: the robust pass inside
%   the loop ('Influence', 'Tuning', model.rowWeights, out.robust, the
%   collapsed refusal) -- the one-loop composition of IRLS and the joint
%   estimator. With assistance from Claude (Anthropic), model Claude Fable 5.1.
%   Updated 2026-09-08 (CEST, UTC+02:00), R6 WP-98 (finding F-01 FIXED): eq.
%   6.5.28 now multiplies alpha2_m by the factor the WEIGHTING used -- the
%   per-stream level -- instead of the per-arc factor, which alpha2 cancels
%   and therefore never responded to. The arc-factor ratio settles instead of
%   growing 8.4x a pass; runs with arcFactors = false are unchanged, because
%   there the two numbers are identical. With assistance from Claude
%   (Anthropic), model Claude Opus 5.
%   Updated 2026-09-08 (CEST, UTC+02:00), R6 WP-97 (finding F-01): the loop
%   records the SPREAD of every arc-factor stream pass by pass and warns
%   (arcFactorsDiverging) when one keeps growing -- the per-arc factors are
%   estimated but never fed back into the weighting, so they compound. The
%   guard reports; it does not correct, because the mechanism is not yet
%   established. With assistance from Claude (Anthropic), model Claude Opus 5.
%   Updated 2026-09-06 (CEST, UTC+02:00), R6 WP-90: the per-axis descale runs
%   for every arc that carries a per-epoch Cep (scalar streams included); the
%   noCep refusal is retired. With assistance from Claude (Anthropic), model
%   Claude Fable 5.1.
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (a) (F-07): typed
%   `arguments` signature replaces inputParser. Option names were already
%   UpperCamelCase and are UNCHANGED. NOTE (premise corrected round a1-2,
%   measured on CI): MATLAB matches comma-syntax Name-Value NAMES
%   case-insensitively even under arguments blocks, so legacy casings
%   keep working; what tightened are the typed constraints, unknown-name
%   rejection, and enum VALUE casing (StepSchedule accepts exactly
%   "none"/"decay" -- strcmpi used to take 'Decay'). With assistance from
%   Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-40 r2 (F-15): optional
%   stabilised alternation -- 'StepSchedule','decay' (Robbins-Monro damped
%   psdVCE update) and 'AvgTail' (Polyak tail average of the returned
%   spectrum); defaults are bitwise the old path. With assistance from
%   Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-34 (fixes F-02, CRITICAL): arc
%   factors are matched by KEY (s<sys>_<obsType>_d<day>_a<arc>, pinned on
%   the pass-1 enumeration), never by position; skipped arcs keep their
%   previous factor and are counted in out.skips; a changed key set raises
%   arcMismatch. New outputs .arcKeys/.skips; Sigma2_0 aligns to .arcKeys.
%   With assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-18 (CEST, UTC+02:00), WP-33 (fixes F-03, CRITICAL): step
%   7 subtracts the per-system local fit B*y (recovered once per pass as
%   redundancy does) and the per-lag projector absorbs the local block
%   exactly -- shared implementation jointArcParts + redundancyLags
%   'LocalBbar'/'LocalLl'; evidence v10review.py T3. With assistance from
%   Claude (Anthropic), model Claude Fable 5.
%     lsqvce:adjustVCEjoint:arcMismatch  the arc KEY SET changed between
%                                        passes (the loader must be
%                                        deterministic; a covariance failure
%                                        is a SKIP, a vanished enumeration
%                                        entry is an identity error) (WP-34)
%     lsqvce:adjustVCEjoint:badSigma2_0  Sigma2_0 does not match the pass-1
%                                        arc count (align to out.arcKeys)
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    % WP-39 stage (a): typed `arguments` signature (rule 5) replaces
    % inputParser. Names were already UpperCamelCase here; NV name
    % matching stays case-insensitive (MATLAB semantics, measured round
    % a1-2) -- the typed constraints and enum-value casing are what
    % tightened. Release note: HANDOVER §9 WP-39.
arguments
    accFcn (1,1) function_handle
    opt.NumLags  (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(opt.NumLags, 2)} = 60
    opt.C0       double = []
    opt.Sigma2_0 double = []
    opt.MaxIter  (1,1) double {mustBeInteger, mustBePositive} = 10
    opt.Tol      (1,1) double {mustBePositive} = 1e-2
    opt.Probes   (1,1) double {mustBeInteger, mustBePositive} = 100
    opt.Seed     double {mustBeScalarOrEmpty} = []
    opt.Norm     (1,1) string = "ortho"
    opt.MeanType (1,1) string = "robust"
    opt.StepSchedule (1,1) string {mustBeMember(opt.StepSchedule, ["none","decay"])} = "none"
    opt.AvgTail  (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    opt.SkipThreshold (1,1) double {mustBeInRange(opt.SkipThreshold, 0, 1)} = 0.2
    % R5 M3 3d: per-stream policy (WP-83, decision R5-D4). Streams are keyed
    % 'satKey|obsType' (matlab.lang.makeValidName of the system id).
    opt.Pool (1,1) string {mustBeMember(opt.Pool, ["all","obsType","stream"])} = "all"
    opt.Estimate = []            % containers.Map key -> "none"|"scale"|"law"; [] = "law" for all
    opt.ArcFactors = []          % containers.Map key -> logical;              [] = true for all
    % R6 WP-105: the ROBUST pass inside the loop. "none" = the plain estimator
    % (bitwise the old behaviour). Otherwise every pass ends with one residual
    % sweep whose LAW-whitened residuals give the next pass's per-observation
    % factors (lsqvce.kernel.robustWeights); the factors ride inside every
    % arc's whitener (model.rowWeights -> the factory's UserOpt.rowWeights),
    % so the accumulation, the probe redundancies and the component update
    % of the next pass all see ONE re-weighted problem. Prior streams are
    % never re-weighted (a down-weighted datum is a rank defect, not an
    % outlier).
    opt.Influence (1,1) string {mustBeMember(opt.Influence, ["none","huber","igg3","bisquare"])} = "none"
    opt.Tuning double = []       % tuning constants of the influence ([] = its defaults)
    opt.OnPass = []              % R5 M5 5c': function handle fired after EVERY pass
                                 % with the ctx below (see lsqvce.util.onPass); return
                                 % false to stop the loop (stopReason "userStop")
    opt.OnPassContext struct = struct()   % extra ctx fields (fit: data, processingOptions)
    opt.Verbose  (1,1) logical = false
end
    o = opt;
    v = logical(o.Verbose);

    nLag = o.NumLags;

    % ---- starting model (step 1 / eq. 6.5.62, degeneracy-aware) -----------
    if isempty(o.C0)
        C = [1; 0.3*ones(nLag-1,1)];     % NOT a delta: see the header warning
    else
        C = o.C0(:);
        if numel(C) ~= nLag
            error('lsqvce:adjustVCEjoint:badC0', ...
                  'C0 has %d entries but NumLags is %d.', numel(C), nLag);
        end
    end
    sigma2 = o.Sigma2_0(:);              % may be empty; sized on first sweep

    hist = struct('maxDev', [], 'Shat', [], 'skips', {cell(0,3)}, ...
                  'streams', {{}}, 'sigma', [], 'white', [], 'whitePval', []);
    converged = false;  userStop = false;
    accF = []; solF = []; Shat = [];
    Cmap = containers.Map('KeyType', 'char', 'ValueType', 'any');   % R5 M3 3d: C per pool group
    pinnedKeys = {};                     % WP-34: pass-1 arc identity
    pinnedStream = {};                   % R6 WP-97: their stream keys, once
    s2Map = [];                          % keyed factors; built on pass 1
    spreadHist = containers.Map('KeyType', 'char', 'ValueType', 'any');   % R6 WP-97 (F-01)

    shatTail = {};                               % WP-40 r2: AvgTail buffer
    lastSol = [];                                % R5 M7 7b: the previous pass's solution
    robust = ~strcmp(o.Influence, "none");       % R6 WP-105
    W = containers.Map('KeyType', 'char', 'ValueType', 'any');   % per-arc row factors
    fac = [];  rkeys = cell(0, 4);  rscale = [];  rstreams = {};
    hist.robustScale = [];  hist.nDown = [];  hist.nZero = [];
    for it = 1:o.MaxIter
        % ---- steps 1-5: build under the current model, accumulate, solve --
        % R5 M3 3b: the caller's override needs a LEVEL per stream (the pooled
        % C is a shape): the mean arc factor of every 'satKey|obsType' stream.
        % R5 M7 7b: model.lastSol (x, Cxx, s0sq of the previous pass; [] on
        % pass 1) lets the factory update anything that has its own variance
        % component outside the streams -- the prior groups (Koch-Kusche).
        model = struct('C', C, 'sigma2', sigma2, 'iter', it, ...
                       'streamSigma2', streamMeans(pinnedKeys, sigma2), ...
                       'Cgroups', Cmap, 'Pool', char(o.Pool), 'Estimate', o.Estimate, ...
                       'lastSol', lastSol, 'rowWeights', W);   % R6 WP-105: the robust factors
        acc   = accFcn(model);
        sol   = acc.solve();
        lastSol = struct('x', sol.x(:), 'Cxx', sol.Cxx, 's0sq', sol.s0sq);
        accF  = acc;  solF = sol;

        N = acc.reducedNormalMatrix();
        [U, pd] = chol((N + N.')/2);      % upper: N = U'*U  (Ellmer's U)
        if pd > 0
            error('lsqvce:adjustVCEjoint:notPD', ...
                 ['The reduced normal matrix is not positive definite at outer ' ...
                  'iteration %d. Add constraints or a prior before estimating ' ...
                  'the stochastic model.'], it);
        end

        % ---- step 6: the shared probe, ONCE for every arc -----------------
        sp = lsqvce.engine.sharedProbe(U, 'numProbes', o.Probes, 'seed', o.Seed);

        % ---- steps 7-8: re-stream arcs, per-lag Omega and s ---------------
        % WP-34 (fixes F-02): arcs are KEYED, the factor map is pinned on the
        % PASS-1 ENUMERATION, and a skipped arc keeps its previous factor --
        % one skip can no longer shift every later arc's factor by position.
        [omegaT, sT, keys, keysAll, skips, wstats] = ...
            sweepArcs(acc, sol, sp.UiZ, nLag, C, s2Map, v);
        nArc = numel(keys);
        if nArc == 0
            error('lsqvce:adjustVCEjoint:noArcs', ...
                  'No arcs were enumerated -- nothing to estimate from.');
        end
        % R5 M1: a skip RATE above the threshold is a refusal, not a smaller
        % estimate. The skipped arcs keep their factors and the run continues
        % on the rest only while they are few; beyond SkipThreshold the law
        % estimated from the survivors would describe a different data set.
        nAll = numel(keysAll);
        if nAll > 0 && size(skips, 1) / nAll > o.SkipThreshold
            error('lsqvce:adjustVCEjoint:skipThreshold', ...
                 ['Pass %d skipped %d of %d arcs (%.0f %% > SkipThreshold %.0f %%): ' ...
                  '%s. A noise law estimated from the remaining arcs would not ' ...
                  'describe the data; fix the covariance model or raise ' ...
                  'Options.SkipThreshold deliberately.'], it, size(skips, 1), nAll, ...
                  100*size(skips, 1)/nAll, 100*o.SkipThreshold, skipSummary(skips));
        end
        if isempty(pinnedKeys)
            % Pass 1: pin the arc identity on the ENUMERATION (skipped arcs
            % included -- a pass-1 skip keeps the default factor 1).
            pinnedKeys = keysAll;
            pinnedStream = cellfun(@streamKeyOf, pinnedKeys(:), 'UniformOutput', false);
            s2Map = containers.Map('KeyType','char', 'ValueType','double');
            if isempty(o.Sigma2_0)
                for kk = 1:numel(pinnedKeys), s2Map(pinnedKeys{kk}) = 1; end
            else
                if numel(o.Sigma2_0) ~= numel(pinnedKeys)
                    error('lsqvce:adjustVCEjoint:badSigma2_0', ...
                        ['Sigma2_0 has %d entries but pass 1 enumerated %d ' ...
                         'arcs (see out.arcKeys for the pinned order).'], ...
                         numel(o.Sigma2_0), numel(pinnedKeys));
                end
                for kk = 1:numel(pinnedKeys)
                    s2Map(pinnedKeys{kk}) = o.Sigma2_0(kk);
                end
            end
        elseif ~isequal(sort(keysAll), sort(pinnedKeys))
            missing = setdiff(pinnedKeys, keysAll);
            newk    = setdiff(keysAll, pinnedKeys);
            error('lsqvce:adjustVCEjoint:arcMismatch', ...
                ['The arc key set changed between passes -- the loader must ' ...
                 'be deterministic across passes (an arc that merely fails ' ...
                 'its covariance is SKIPPED and kept, but one that vanishes ' ...
                 'from the enumeration is an identity error). Missing: %s. ' ...
                 'New: %s.'], ...
                strjoin([missing, {'(none)'}], ', '), ...
                strjoin([newk,    {'(none)'}], ', '));
        end
        if ~isempty(skips)
            for kk = 1:size(skips, 1)
                hist.skips(end+1, :) = {it, skips{kk,1}, skips{kk,2}}; %#ok<AGROW>
            end
            if v
                lsqvce.util.log.msg(['  [joint %2d] SKIPPED %d arc(s) this ' ...
                    'pass (factors kept): %s\n'], it, size(skips,1), ...
                    strjoin(skips(:,1).', ', '));
            end
        end

        % ---- step 9: co-update the PSD AND the arc factors ----------------
        s2Cols = cellfun(@(k) s2Map(k), keys(:));    % survived, keyed lookup
        % WP-40 r2 (F-15): the alternation is a stochastic approximation --
        % the plain step (Gamma = 1 every pass) re-amplifies per-bin
        % estimation noise and is BISTABLE for short arcs (corr collapse to
        % negative values on unlucky seeds; wp33_loop_mirror.py grid in the
        % round note). 'StepSchedule','decay' uses the Robbins-Monro step
        % Gamma_it = 2/(2+it): full first steps, averaging later. Measured
        % across the mirror grid: min corr rises from -0.66 to ~0.9 with
        % barely any ceiling loss. Default 'none' is bitwise the old path.
        if strcmpi(o.StepSchedule, 'decay')
            gammaIt = 2 / (2 + it);
        else
            gammaIt = 1;
        end
        % R5 M3 3d: the law is estimated PER GROUP of streams (Pool), from the
        % arcs of streams whose policy is "law"; "scale" streams get arc
        % factors only (their mean is the Helmert factor); "none" streams are
        % untouched. out.C stays the 'all'-group C (the first group's otherwise).
        streamOf = cellfun(@streamKeyOf, keys(:), 'UniformOutput', false);
        polOf    = cellfun(@(sk) policyOf(o.Estimate, sk), streamOf, 'UniformOutput', false);
        grpOf    = cellfun(@(sk) groupKeyOf(o.Pool, sk), streamOf, 'UniformOutput', false);
        lawCols  = strcmp(polOf, 'law');
        groups   = unique(grpOf(lawCols));
        pvMax = 0;  Shat = [];
        for g = 1:numel(groups)
            cols = lawCols & strcmp(grpOf, groups{g});
            Cg = C;
            if isKey(Cmap, groups{g}), Cg = Cmap(groups{g}); end
            pvg = lsqvce.engine.psdVCE(omegaT(:, cols), sT(:, cols), Cg, 'Norm', o.Norm, 'Gamma', gammaIt);
            Cmap(groups{g}) = pvg.Cnew;
            pvMax = max(pvMax, pvg.maxDev);
            if isempty(Shat), Shat = pvg.Shat; end         % the first group's spectrum for the record
        end
        if ~isempty(groups)
            if isKey(Cmap, 'all'), C = Cmap('all'); else, C = Cmap(groups{1}); end
        end
        pv = struct('maxDev', pvMax);
        if isempty(Shat), Shat = lsqvce.kernel.dctI(C, 'Norm', o.Norm); end
        % R6 WP-98 (finding F-01): eq. 6.5.28 multiplies alpha2_m by the factor
        % the WEIGHTING used. What the weighting used is the per-STREAM level --
        % LSQVCE.ENGINE.JOINTMODELOPT can express nothing finer -- so that is
        % what must multiply here. Multiplying by the per-ARC factor instead
        % makes the update a geometric recursion with no fixed point, because
        % alpha2_m cancels its own arc's factor exactly (both omegaT and sT are
        % scaled by it in sweepArcs) and so never responds to it: MEASURED
        % growth of 8.4x per pass over four decades. For a stream with
        % arcFactors = false every arc already carries the one stream level, so
        % these two are the same number and such runs cannot move.
        s2Applied = s2Cols;
        for sk = unique(streamOf(:)).'
            cols = strcmp(streamOf, sk{1});
            s2Applied(cols) = mean(s2Cols(cols));
        end
        av = lsqvce.engine.arcVarianceFactors(omegaT, sT, s2Applied, ...
                                       'Norm', o.Norm, 'MeanType', o.MeanType);
        % WP-40 r2 (F-15): Polyak tail averaging of the RETURNED spectrum
        % (log domain, last AvgTail passes). The loop keeps iterating on the
        % raw C; only the final Shat/C are averaged -- measured to add ~0.02
        % corr on top of the decay schedule (and to do NOTHING without it:
        % averaging cannot rescue an undamped collapse; mirror grid).
        if o.AvgTail > 0
            shatTail{end+1} = log(max(Shat, realmin));       %#ok<AGROW>
        end
        % R5 M4 4b: a stream WITHOUT arc factors still needs its LEVEL estimated
        % (the PSD shape is normalised; the level rode on the arc factors).
        % The level is the Helmert factor of the stream: the mean over its arcs
        % of gamma2_m = alpha2_m * sigma2_m (eq. 6.5.28, the arc's NEW absolute
        % factor) WITHOUT the unit-mean normalisation of the ellmer mode (that
        % normalisation exists because the PSD carries the level there; here
        % the shape is normalised, so the level must be absolute). One common
        % factor per stream. Found by the gravity demo: planted 1 / 2 / 0.5
        % came back 1 / 1 / 1 without it.
        fStream = containers.Map('KeyType', 'char', 'ValueType', 'double');
        for sk = unique(streamOf(:)).'
            cols = strcmp(streamOf, sk{1}) & ~strcmp(polOf, 'none') & ~arcFactorsOf(o.ArcFactors, sk{1});
            if ~any(cols), continue, end
            g2 = av.gamma2(cols);
            g2 = g2(isfinite(g2) & g2 > 0);
            if ~isempty(g2), fStream(sk{1}) = mean(g2); end
        end
        for kk = 1:nArc
            if strcmp(polOf{kk}, 'none'), continue, end
            if arcFactorsOf(o.ArcFactors, streamOf{kk})
                s2Map(keys{kk}) = av.sigma2(kk);     % keyed write-back, per arc
            elseif isKey(fStream, streamOf{kk})
                s2Map(keys{kk}) = fStream(streamOf{kk});    % one level per stream
            end
        end
        sigma2 = cellfun(@(k) s2Map(k), pinnedKeys(:));   % pinned order

        % R6 WP-97 (finding F-01): record the SPREAD of every arc-factor
        % stream, pass by pass. An arc factor is a RELATIVE weight inside its
        % stream, so its spread must settle; one that keeps growing is not an
        % estimate. arcSpreadGuard reads this after the loop.
        for sk = unique(streamOf(:)).'
            if ~arcFactorsOf(o.ArcFactors, sk{1}), continue, end
            f = sigma2(strcmp(pinnedStream, sk{1}));
            f = f(isfinite(f) & f > 0);
            if numel(f) < 2, continue, end
            if isKey(spreadHist, sk{1}), r0 = spreadHist(sk{1}); else, r0 = []; end
            spreadHist(sk{1}) = [r0; max(f)/min(f)];
        end

        % ---- step 9c (R6 WP-105): the robust pass -- next pass's factors ----
        % One residual sweep under THIS pass's accumulator; the factors of
        % this pass are divided out (lawWhitenedResiduals), so the influence
        % function sees the standardised residual under the current law,
        % level, arc factor and Cep. The factors are written for the NEXT
        % pass, one pass behind the model -- the same lag adjustRobust has.
        if robust
            res = lsqvce.engine.residuals(acc, sol, 'Whitened', true, 'Collect', true, 'Verbose', false);
            [vw, rkeys, sizes] = lawWhitenedResiduals(res, W);
            axes = arcAxes(res);
            [fac, rscale, rstreams] = perStreamWeights(vw, rkeys, sizes, axes, o.Influence, o.Tuning);
            nAlive = sum(fac > 0);
            if nAlive < size(sol.x, 1) + 1
                error('lsqvce:adjustVCEjoint:collapsed', ...
                     ['The influence function left %d observations with a non-zero ' ...
                      'weight, fewer than the %d parameters + 1: more than half the ' ...
                      'data is discrepant, which is beyond the breakdown point of ' ...
                      'any M-estimator. This is not an outlier problem -- check the model.'], ...
                      nAlive, size(sol.x, 1));
            end
            off = 0;
            for kk = 1:size(rkeys, 1)
                n = sizes(kk);
                W(lsqvce.util.rowWeightKey(rkeys{kk,1}, rkeys{kk,3}, rkeys{kk,2}, rkeys{kk,4})) = fac(off + (1:n));
                off = off + n;
            end
            hist.robustScale(it,1:numel(rscale)) = rscale;       %#ok<AGROW>  one MAD scale per stream
            hist.nDown(it,1) = sum(fac < 0.5 & fac > 0);         %#ok<AGROW>
            hist.nZero(it,1) = sum(fac == 0);                    %#ok<AGROW>
            if v
                lsqvce.util.log.msg('  [joint %2d] robust (%s): scales %s, %d down-weighted, %d at zero\n', ...
                    it, o.Influence, mat2str(rscale, 3), hist.nDown(it), hist.nZero(it));
            end
        end

        hist.maxDev(it,1)  = pv.maxDev;              %#ok<AGROW>
        hist.Shat(it,:)    = Shat(:).';              %#ok<AGROW>
        % R5 M4 4a: the per-stream pass record Results.convergence reads --
        % the stream's level factor and its Ljung-Box whiteness on the
        % law-whitened residuals of THIS pass (K = 20 lags, alpha 0.05).
        [hist.streams, sigRow, whiteRow, pRow] = whitenessRow(wstats, streamMeans(pinnedKeys, sigma2));
        hist.sigma(it, :)     = sigRow;              %#ok<AGROW>
        hist.white(it, :)     = whiteRow;            %#ok<AGROW>
        hist.whitePval(it, :) = pRow;                %#ok<AGROW>

        if v
            lsqvce.util.log.msg(['  [joint %2d] max|alpha2-1| = %.4e   ' ...
                'arcs = %d   gammaBar = %.4f\n'], it, pv.maxDev, nArc, av.gammaBar);
        end

        % ---- step 9b: the user's per-pass inspection hook (R5 M5 5c') -----
        ctx = struct('driver', 'adjustVCEjoint', 'pass', it, 'x', sol.x(:), ...
                     'Cxx', sol.Cxx, 'sigma0', sol.sigma0, 's0sq', sol.s0sq, 'C', C, 'Cgroups', Cmap, ...
                     'sigma2', sigma2, 'arcKeys', {pinnedKeys}, 'streams', {hist.streams}, ...
                     'sigma', hist.sigma(it, :), 'white', hist.white(it, :), ...
                     'whitePval', hist.whitePval(it, :), 'maxDev', pv.maxDev, ...
                     'skips', {skips}, 'acc', acc, ...
                     'factor', fac, 'keys', {rkeys}, 'scale', rscale);   % R6 WP-105 ([] when Influence "none")
        userStop = ~lsqvce.util.onPass(opt.OnPass, opt.OnPassContext, ctx);

        % ---- step 10: converged? -----------------------------------------
        if pv.maxDev <= o.Tol
            converged = true;
            break
        end
        if userStop, break, end
    end

    if ~converged && ~userStop
        lsqvce.util.log.warn('lsqvce:adjustVCEjoint:notConverged', ...
           ['The joint stochastic model did not settle in %d outer iterations ' ...
            '(last max|alpha2-1| = %.3e > tol %.1e). Inspect out.history.Shat: ' ...
            'a PSD that keeps moving usually means too few arcs per frequency, ' ...
            'or a starting covariance far from the data.'], ...
            o.MaxIter, hist.maxDev(end), o.Tol);
    end

    arcSpreadGuard(spreadHist);          % R6 WP-97 (F-01): report, never silence

    if o.AvgTail > 0 && numel(shatTail) >= o.AvgTail
        Shat = exp(mean(cat(2, shatTail{end-o.AvgTail+1:end}), 2));
        C    = lsqvce.kernel.dctI(Shat, 'Norm', o.Norm, 'Inverse', true);  % psdVCE's exact inverse
    end
    out = struct('C', C, 'Shat', Shat, 'sigma2', sigma2, 'acc', accF, ...
                 'userStop', userStop, ...
                 'Cgroups', Cmap, ...
                 'sol', solF, 'converged', converged, 'nIter', numel(hist.maxDev), ...
                 'history', hist, ...
                 'arcKeys', {pinnedKeys(:)}, ...
                 'skips', {hist.skips}, ...
                 'robust', struct('influence', char(o.Influence), 'scale', rscale, ...
                                  'streams', {rstreams}, 'factor', fac, 'keys', {rkeys}));   % R6 WP-105
end

%==========================================================================
function [fac, scale, streams] = perStreamWeights(vw, keys, sizes, axes, influence, tuning)
%PERSTREAMWEIGHTS  Robust factors PER STREAM (R6 WP-105, round 2).
%   The influence function standardises residuals by a MAD scale, and that
%   scale must be ONE PER STREAM and must NOT be the estimated level: the
%   first design used a global MAD and it drove the level of a clean stream
%   into a spiral (0.377 measured against 0.80 expected) -- a global scale
%   mis-standardises every stream whose level differs from the mix, the
%   next level is estimated from the mis-weighted residuals, and the next
%   scale sees that. With a MAD per stream the standardisation is
%   v / MAD(v), invariant to the level the whitener applied, so the weights
%   do not feed on the level they will shape: the level converges to the
%   retained-observation ratio, not to zero. Prior streams get factor 1.
%   A MULTI-AXIS arc (axes(kk) = d > 1, rows epoch-major) is weighted PER
%   EPOCH (round 4): the influence acts on the epoch's RMS over its d
%   law-whitened rows -- the multivariate M-estimate of a vector
%   observation -- and the d rows of an epoch share the factor. That is what
%   the sweep needs: its temporal covariance is per axis and per epoch, so a
%   factor must be per epoch to be expressible there.
%   OUTPUT fac n x 1; scale 1 x S (the MAD per stream); streams 1 x S cell
%          'satKey|obsType', sorted.
    n = numel(vw);
    fac = ones(n, 1);
    skey = cell(size(keys, 1), 1);
    for kk = 1:size(keys, 1)
        skey{kk} = sprintf('%s|%s', matlab.lang.makeValidName(char(string(keys{kk,1}))), ...
                           char(string(keys{kk,2})));
    end
    streams = unique(skey).';
    streams = streams(~startsWith(streams, 'prior_') & ~contains(streams, '|prior_'));
    scale = nan(1, numel(streams));
    ends = cumsum(sizes(:));  starts = ends - sizes(:) + 1;
    for si = 1:numel(streams)
        arcs = find(strcmp(skey, streams{si}));
        rows = [];  d = 1;
        for kk = arcs(:).'
            rows = [rows, starts(kk):ends(kk)]; %#ok<AGROW>
            d = max(d, axes(kk));
        end
        if d == 1
            [fac(rows), scale(si)] = lsqvce.kernel.robustWeights(vw(rows), influence, tuning);
        else
            % per-epoch RMS over the d rows; a factor per epoch, replicated
            u = vw(rows);
            uEp = sqrt(mean(reshape(u, d, []).^2, 1)).';
            [fEp, scale(si)] = lsqvce.kernel.robustWeights(uEp, influence, tuning);
            fac(rows) = repelem(fEp(:), d);
        end
    end
end

function axes = arcAxes(res)
%ARCAXES  The axis count of every arc, in the sweep order of lawWhitenedResiduals.
    axes = [];
    for k = 1:numel(res.sys)
        for j = 1:numel(res.sys(k).obs)
            axes(end+1) = res.sys(k).obs(j).nAxis; %#ok<AGROW>
        end
    end
end

%==========================================================================
function [omegaT, sT, keys, keysAll, skips, wstats] = sweepArcs(acc, sol, UiZ, nLag, C, s2Map, v)
%SWEEPARCS  Steps 7-8: re-stream every arc, returning the WEIGHTED per-lag
%   temporal vectors Omega^m_n and s^m_n (eqs. 6.5.22-6.5.23), one COLUMN per
%   arc. Mirrors the arc enumeration used by LSQVCE.ENGINE.REDUNDANCY and
%   LSQVCE.ENGINE.STREAMCONTRIBUTION, so grouping and whitening stay consistent.
%   WP-34 (fixes F-02): every arc is KEYED s<sys>_<obsType>_d<day>_a<arc>
%   (LSQVCE.ENGINE.REDUNDANCY's 'arc' key format, verbatim). OUTPUTS:
%     keys     1 x M   cell  keys of the arcs that produced columns, ALIGNED
%                            to the columns of omegaT/sT
%     keysAll  1 x Ma  cell  keys of every ENUMERATED arc (the identity the
%                            caller pins on pass 1)
%     skips    S x 2   cell  {key, reason} for every arc skipped THIS pass
%   The arc factor is looked up BY KEY in s2Map (containers.Map; missing or
%   empty map -> 1), never by position: one skip can no longer shift every
%   later arc's factor.
%
%   WHY ONE `C` WEIGHTS EVERY ARC, when Cgroups holds a law per pool group
%   (R6 WP-100 -- read this before "fixing" it, because it has been proposed
%   twice and measured twice). The laws ARE fitted per group: psdVCE runs once
%   per pool group and writes Cmap. What is shared is the WEIGHT in eqs.
%   6.5.22-6.5.23, and sharing it is what keeps the iteration from feeding on
%   itself. Threading each group's own freshly fitted law into its own arcs
%   was implemented and MEASURED on three two-stream fixtures: the second
%   group diverged every time -- planted lag-1 0.70 came back 0.852, planted
%   0.50 came back 1.331, planted 0.40 came back 1.365, and sigma0 fell to
%   0.709. With the shared C the same fixtures give 0.198 / 0.477 against
%   planted 0.20 / 0.50, sigma0 0.997. The behaviour is pinned by
%   test_vcePolicy>testPoolStreamGivesEachStreamItsOwnLaw.
    ids    = acc.systemIDs();
    x      = sol.x(:);
    omegaT = [];
    sT     = [];
    keys   = {};
    keysAll = {};
    skips  = cell(0, 2);
    nArc   = 0;
    wstats = containers.Map('KeyType', 'char', 'ValueType', 'any');   % stream -> [sums; m]

    for k = 1:numel(ids)
        st = acc.getStore(ids(k));
        uo = st.uo;
        % WP-33: recover this system's LOCAL solution and Schur pieces ONCE
        % per pass, exactly as LSQVCE.ENGINE.REDUNDANCY does -- the residual is the
        % full-vector one (Ellmer step 7, Table 6.2) and the projector gains
        % the exact local block via jointArcParts + redundancyLags.
        nlS = size(st.Nlc, 1);
        if nlS > 0
            ySys = st.Nll_r \ (st.nl_r - st.Nlc * x);
            Gsys = st.Nll_r \ st.Nlc;
            [LlS, plS] = chol((st.Nll_r + st.Nll_r.')/2, 'lower');
            if plS > 0
                error('lsqvce:adjustVCEjoint:localNotPD', ...
                     ['System %s: the local normal matrix is not positive ' ...
                      'definite -- the local block cannot be absorbed into ' ...
                      'the projector.'], string(st.sysID));
            end
        else
            ySys = [];  Gsys = [];  LlS = [];
        end
        for d = st.dayList(:).'
            meta = st.dataFcn('meta', d, [], uo);
            for ot = 1:numel(meta.obsTypes)
                obsType = meta.obsTypes{ot};
                for a = meta.arcs{ot}(:).'
                    key = sprintf('s%s_%s_d%s_a%s', string(st.sysID), ...
                                  string(obsType), string(d), string(a));
                    keysAll{end+1} = key;                    %#ok<AGROW>
                    [om, sn, why, wsum, wm] = oneArc(st, obsType, a, d, uo, UiZ, ...
                                           nLag, x, ySys, Gsys, LlS);
                    if isempty(om)
                        skips(end+1, :) = {key, why};        %#ok<AGROW>
                        continue
                    end
                    sk = sprintf('%s|%s', matlab.lang.makeValidName(char(string(st.sysID))), char(obsType));
                    if isKey(wstats, sk), w0 = wstats(sk); else, w0 = zeros(WHITE_LAGS() + 2, 1); end
                    wstats(sk) = w0 + [wsum; wm];
                    nArc = nArc + 1;
                    % weight by THIS ARC'S factor (keyed, F-02) and lag cov
                    if ~isempty(s2Map) && isKey(s2Map, key)
                        s2 = s2Map(key);
                    else
                        s2 = 1;                              % pass 1 default
                    end
                    keys{end+1}     = key;                   %#ok<AGROW>
                    omegaT(:, nArc) = s2 * C .* om;          %#ok<AGROW>
                    sT(:, nArc)     = s2 * C .* sn;          %#ok<AGROW>
                end
            end
        end
        if v
            lsqvce.util.log.msg('%sswept system %s\n', lsqvce.util.verbosity.banner('Joint estimation'), string(st.sysID));
        end
    end
end

%==========================================================================
function [om, sn, why, wsum, wm] = oneArc(st, obsType, arc, day, uo, UiZ, nLag, x, ySys, Gsys, LlS)
%ONEARC  Decorrelate one arc, form its residual, and return the per-lag Omega
%   and s from LSQVCE.KERNEL.REDUNDANCYLAGS. WP-34: a skipped arc REPORTS why --
%   `why` is '' on success, 'gap' (no data) or 'nonPD' (covariance not
%   positive definite this pass); silence was finding F-02's accomplice.
%   WP-33: ySys/Gsys/LlS are the system's local solution, Schur factor and
%   local Cholesky (computed once per pass by sweepArcs); the residual is
%   ebar = W^-T(l - A x - B y) and the projector absorbs the local block
%   exactly (jointArcParts + redundancyLags 'LocalBbar'/'LocalLl').
    uo.obsType = obsType;
    uo.day     = day;
    uo.arc     = arc;
    why = '';
    wsum = zeros(WHITE_LAGS + 1, 1);  wm = 0;     % R5 M4 4a: whiteness collector
    a = st.dataFcn(day, obsType, arc, uo);
    if isempty(a) || ~isfield(a, 'l') || isempty(a.l)
        om = [];  sn = [];  why = 'gap';  return
    end

    A = a.A;  l = a.l(:);
    m = size(A, 1);
    B = zeros(m, 0);
    if isfield(a, 'B') && ~isempty(a.B), B = a.B; end
    nAx = lsqvce.util.getfielddef(a, 'nAxis', 1);
    hasCep = isfield(a, 'Cep') && ~isempty(a.Cep);
    if nAx > 1 || hasCep
        % R5 M3 (finding H-15): a d-axis stream is d single-axis series that
        % share the TEMPORAL law once the FIXED per-epoch covariance has been
        % removed (the separable model). So: descale every epoch by Cep^-1/2
        % (rows stay epoch-major), then treat axis k = rows k:d:m as one arc of
        % m/d epochs against the epoch-indexed temporal covariance a.Cfcn, and
        % ADD the per-lag statistics over the axes -- the lag is a lag in
        % epochs, the dense factor is (m/d)^2 instead of m^2, and one arc
        % factor covers the arc. R6 WP-90: the branch is entered whenever the
        % arc CARRIES per-epoch error information -- a scalar stream with a
        % per-epoch sigma (d = 1) is descaled exactly like a 3-axis one; a
        % multi-axis arc without Cep (a raw loader) is D = I.
        mE = m / nAx;
        if mE ~= round(mE)
            error('lsqvce:adjustVCEjoint:badAxisRows', ...
                  'Arc (sys %s, %s, day %s, arc %s): %d rows are not a multiple of %d axes.', ...
                  string(st.sysID), string(obsType), string(day), string(arc), m, nAx);
        end
        Cep = ones(mE, nAx);
        if hasCep, Cep = a.Cep; end
        % pure descale: AR(0) (a = 1, s = 1, r = 1) through the separable kernel
        [Ad, Bd, ld] = lsqvce.kernel.separableWhiten(1, 1, 1, Cep, nAx, A, B, l, true);
        Sig = a.Cfcn(1:mE, uo);                       % TEMPORAL, over epochs
        if isvector(Sig), Sig = diag(Sig(:)); end
        [W, pd] = chol((Sig + Sig.')/2);
        if pd > 0
            om = [];  sn = [];  why = 'nonPD';  return
        end
        nLagA = min(nLag, mE);
        om = zeros(nLag, 1);  sn = zeros(nLag, 1);
        for k = 1:nAx
            rows = k:nAx:m;
            [AtilBar, Bbar, ebar] = jointArcParts(Ad(rows, :), Bd(rows, :), ld(rows), W, ...
                                                  x(1:size(A,2)), ySys, Gsys);
            r  = lsqvce.kernel.redundancyLags(AtilBar, W, UiZ, ebar, nLagA, ...
                                       'LocalBbar', Bbar, 'LocalLl', LlS);
            om = om + padTo(r.omega, nLag);
            sn = sn + padTo(r.s,     nLag);
            wsum = wsum + whiteLagSums(ebar);  wm = wm + numel(ebar);
        end
        return
    end
    if isfield(a, 'Cfcn') && ~isempty(a.Cfcn)
        Sig = a.Cfcn(1:m, uo);
        if isvector(Sig), Sig = diag(Sig(:)); end
        [W, pd] = chol((Sig + Sig.')/2);
        if pd > 0
            om = [];  sn = [];  why = 'nonPD';  return   % skipped, REPORTED
        end
    else
        % A Wfcn-only loader hands back a whitening OPERATOR, not the explicit
        % triangular W that redundancyLags needs for Rtilde (eq. 6.5.48) and
        % etilde (eq. 6.5.54). Rather than silently drop the arc -- which would
        % bias the estimate by an invisible amount -- say so.
        error('lsqvce:adjustVCEjoint:needCfcn', ...
             ['Arc (sys %s, %s, day %s, arc %s) supplies only a Wfcn whitening ' ...
              'operator. The joint spectral estimator needs the EXPLICIT arc ' ...
              'covariance to form W with Sigma = W''W, so the loader must also ' ...
              'provide Cfcn. Most LSQVCE.ENGINE.COVMODEL types do; a hand-written ' ...
              'loader may need one added.'], string(st.sysID), string(obsType), ...
              string(day), string(arc));
    end

    % Step 7 with the FULL Table-6.2 vector: subtract B*y, reduce the common
    % design (Atil = A - B G) so the shared probe matches the factored N, and
    % hand the exact local block to the projector (WP-33, fixes F-03).
    [AtilBar, Bbar, ebar] = jointArcParts(A, B, l, W, x(1:size(A,2)), ...
                                          ySys, Gsys);
    r  = lsqvce.kernel.redundancyLags(AtilBar, W, UiZ, ebar, min(nLag, m), ...
                               'LocalBbar', Bbar, 'LocalLl', LlS);
    om = padTo(r.omega, nLag);
    sn = padTo(r.s    , nLag);
    wsum = whiteLagSums(ebar);  wm = numel(ebar);
end

function K = WHITE_LAGS()
%WHITE_LAGS  Lags of the Ljung-Box whiteness statistic per stream (20).
    K = 20;
end

function s = whiteLagSums(e)
%WHITELAGSUMS  sum_i e(i+n) e(i), n = 0..K, of one arc's law-whitened residual
%   (zero-padded FFT, exact to rounding); aggregated per stream by sweepArcs.
    e = e(:);  K = WHITE_LAGS();
    n = numel(e);
    s = zeros(K + 1, 1);
    if n < 2, return, end
    nfft = 2^nextpow2(2*n - 1);
    F = fft(e, nfft);
    r = ifft(real(F .* conj(F)), 'symmetric');
    k = min(K + 1, n);
    s(1:k) = r(1:k);
end

%==========================================================================
function y = padTo(v, n)
%PADTO  Pad a per-lag vector with zeros to length n (eqs. 6.5.22-6.5.23 set the
%   contribution to ZERO for lags beyond the arc length).
    y = zeros(n, 1);
    k = min(numel(v), n);
    y(1:k) = v(1:k);
    y(~isfinite(y)) = 0;
end

function s = skipSummary(skips)
%SKIPSUMMARY  'gap: 3, nonPD: 12' from the {key, reason} skip list.
    if isempty(skips), s = 'none'; return, end
    [why, ~, j] = unique(skips(:, 2));
    n = accumarray(j(:), 1);
    parts = cell(1, numel(why));
    for k = 1:numel(why), parts{k} = sprintf('%s: %d', why{k}, n(k)); end
    s = strjoin(parts, ', ');
end

function M = streamMeans(arcKeys, sigma2)
%STREAMMEANS  containers.Map 'satKey|obsType' -> mean arc factor of the stream.
%   Keys are the pinned arc keys 's<sysID>_<obsType>_d<day>_a<arc>' (WP-34);
%   empty before pass 1 pins them (then every stream's level is its declared
%   sigma, factor 1).
    M = containers.Map('KeyType', 'char', 'ValueType', 'double');
    if isempty(arcKeys) || isempty(sigma2) || numel(arcKeys) ~= numel(sigma2), return, end
    acc = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for k = 1:numel(arcKeys)
        tok = regexp(arcKeys{k}, '^s(.+?)_(.+)_d[^_]+_a', 'tokens', 'once');   % R5-D8 6: obsType may carry underscores
        if isempty(tok), continue, end
        sk = sprintf('%s|%s', matlab.lang.makeValidName(tok{1}), tok{2});
        if ~isKey(acc, sk), acc(sk) = []; end
        v = sigma2(k);
        if isfinite(v) && v > 0, acc(sk) = [acc(sk); v]; end
    end
    ks = keys(acc);
    for k = 1:numel(ks)
        v = acc(ks{k});
        if ~isempty(v), M(ks{k}) = mean(v); end
    end
end

function sk = streamKeyOf(arcKey)
%STREAMKEYOF  's<sysID>_<obsType>_d..' -> 'satKey|obsType' (the override key).
    tok = regexp(arcKey, '^s(.+?)_(.+)_d[^_]+_a', 'tokens', 'once');   % R5-D8 6: obsType may carry underscores
    if isempty(tok), sk = arcKey; return, end
    sk = sprintf('%s|%s', matlab.lang.makeValidName(tok{1}), tok{2});
end

function arcSpreadGuard(spreadHist)
%ARCSPREADGUARD  Refuse to return a DIVERGING set of arc factors in silence.
%   R6 WP-97, kept as a RATCHET by WP-98 (HANDOVER 12.21, 12.22). An arc
%   variance factor is a RELATIVE weight inside its stream, so over the passes
%   its spread must settle. Finding F-01 was that it did not: eq. 6.5.28
%   multiplied by a factor the weighting never used, which made the update a
%   geometric recursion with no fixed point (measured 8.4x per pass). WP-98
%   fixed that, so on this tree nothing provokes this guard any more.
%
%   IT STAYS ANYWAY, and deliberately: it is what would catch this class
%   returning on a user's own data, where no test of ours runs. Its FIRING
%   direction was demonstrated against the pre-fix code and is on the record
%   with numbers -- 3/3 at spreads 4670 / 5600 / 6109 with a degraded arc,
%   silent 3/3 without one. Do not delete it because it is quiet: quiet is now
%   the correct reading, and dev/tests pins that reading.
%
%   THE RULE, and where its numbers come from: fire only when the spread BOTH
%   grew at every one of the last three passes AND ended at least 4x its first
%   pass. Measured silent on the shipped demos even before the fix:
%   demo_grace 4.13 -> 3.27 over six passes.
    ks = keys(spreadHist);
    for i = 1:numel(ks)
        r = spreadHist(ks{i});
        if numel(r) < 3, continue, end
        if all(diff(r(end-2:end)) > 0) && r(end) >= 4 * r(1)
            lsqvce.util.log.warnOnce('lsqvce:adjustVCEjoint:arcFactorsDiverging', ...
               ['Stream %s: the spread of its per-arc variance factors grew from %.3g to ' ...
                '%.3g over %d passes and was still growing at the last one. Per-arc factors ' ...
                'are estimated but never fed back into the weighting (finding F-01), so they ' ...
                'compound instead of settling. Take the ONE-PASS estimate (MaxPasses = 1) ' ...
                'for per-arc factors, or set arcFactors = false on this stream and read its ' ...
                'level instead. A stream with arcFactors = false is unaffected; the level of ' ...
                'THIS stream is the mean of these factors, so it drifts with them.'], ...
                ks{i}, r(1), r(end), numel(r));
        end
    end
end

function p = policyOf(E, sk)
%POLICYOF  'none' | 'scale' | 'law' for a stream; an empty map = 'law' for everybody.
    p = 'law';
    if isempty(E), return, end
    if isKey(E, sk), p = lower(char(E(sk))); end
end

function tf = arcFactorsOf(A, sk)
%ARCFACTORSOF  whether the stream's arcs get factors; an empty map = true.
    tf = true;
    if isempty(A), return, end
    if isKey(A, sk), tf = logical(A(sk)); end
end

function g = groupKeyOf(pool, sk)
%GROUPKEYOF  the pool group a stream's law belongs to.
    switch lower(char(pool))
        case 'all'
            g = 'all';
        case 'obstype'
            g = extractAfter(sk, '|');
        otherwise
            g = sk;                                  % 'stream'
    end
end

function [streams, sig, white, pval] = whitenessRow(wstats, sMeans)
%WHITENESSROW  One convergence-table row: per stream, the level factor
%   sqrt(mean arc factor) and the Ljung-Box test on the pass's law-whitened
%   residuals: Q = m (m+2) sum_k rho_k^2 / (m-k), k = 1..K, ~ chi2_K under
%   whiteness; white = p > 0.05 (the convention of the separate mode's
%   history, so Results.convergence reads both alike).
    streams = sort(keys(wstats));
    S = numel(streams);
    sig = ones(1, S);  white = false(1, S);  pval = nan(1, S);
    K = WHITE_LAGS();
    for j = 1:S
        w = wstats(streams{j});
        m = w(end);  r0 = w(1);
        if isKey(sMeans, streams{j}), sig(j) = sqrt(max(sMeans(streams{j}), eps)); end
        if m <= K + 1 || r0 <= 0, continue, end
        rho = w(2:K+1) / r0;
        Q = m * (m + 2) * sum(rho.^2 ./ (m - (1:K).'));
        pval(j)  = gammainc(Q/2, K/2, 'upper');      % 1 - chi2cdf(Q, K), base MATLAB
        white(j) = pval(j) > 0.05;
    end
end
