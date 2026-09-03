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
%     LSQVCE.ENGINE.ADJUSTVCE alternates: fit the temporal whitening, then the
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
%     lsqvce:adjustVCEjoint:noCep      (R5 M3, H-15) a multi-axis arc without the
%                                   per-epoch covariance the separable model needs
%     lsqvce:adjustVCEjoint:badAxisRows  a multi-axis arc whose row count is not a
%                                   multiple of its axis count
%     lsqvce:adjustVCEjoint:noArcs    a system declares no arcs
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:adjustVCEjoint:localNotPD           System X
%     lsqvce:adjustVCEjoint:notConverged         [warning] The joint stochastic model did not settle in X outer iterations (last...
%
%   See also LSQVCE.ENGINE.ADJUSTVCE (the separable default), LSQVCE.ENGINE.PSDVCE,
%            LSQVCE.ENGINE.ARCVARIANCEFACTORS, LSQVCE.KERNEL.REDUNDANCYLAGS,
%            LSQVCE.ENGINE.SHAREDPROBE
%
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
    s2Map = [];                          % keyed factors; built on pass 1

    shatTail = {};                               % WP-40 r2: AvgTail buffer
    lastSol = [];                                % R5 M7 7b: the previous pass's solution
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
                       'lastSol', lastSol);
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
        av = lsqvce.engine.arcVarianceFactors(omegaT, sT, s2Cols, ...
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
                     'skips', {skips}, 'acc', acc);
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
                 'skips', {hist.skips});
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
    if nAx > 1
        % R5 M3 (finding H-15): a d-axis stream is d single-axis series that
        % share the TEMPORAL law once the FIXED per-epoch covariance has been
        % removed (the separable model). So: descale every epoch by Cep^-1/2
        % (rows stay epoch-major), then treat axis k = rows k:d:m as one arc of
        % m/d epochs against the epoch-indexed temporal covariance a.Cfcn, and
        % ADD the per-lag statistics over the axes -- the lag is a lag in
        % epochs, the dense factor is (m/d)^2 instead of m^2, and one arc
        % factor covers the arc.
        if ~isfield(a, 'Cep') || isempty(a.Cep)
            error('lsqvce:adjustVCEjoint:noCep', ...
                 ['Arc (sys %s, %s, day %s, arc %s) is a %d-axis stream without a ' ...
                  'per-epoch covariance (Cep): declare the stream ''separable'' with ' ...
                  'd = %d and let buildDay return Cep, or supply the epochCovariance hook.'], ...
                  string(st.sysID), string(obsType), string(day), string(arc), nAx, nAx);
        end
        mE = m / nAx;
        if mE ~= round(mE)
            error('lsqvce:adjustVCEjoint:badAxisRows', ...
                  'Arc (sys %s, %s, day %s, arc %s): %d rows are not a multiple of %d axes.', ...
                  string(st.sysID), string(obsType), string(day), string(arc), m, nAx);
        end
        % pure descale: AR(0) (a = 1, s = 1, r = 1) through the separable kernel
        [Ad, Bd, ld] = lsqvce.kernel.separableWhiten(1, 1, 1, a.Cep, nAx, A, B, l, true);
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
        tok = regexp(arcKeys{k}, '^s(.+?)_(.+?)_d', 'tokens', 'once');
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
    tok = regexp(arcKey, '^s(.+?)_(.+?)_d', 'tokens', 'once');
    if isempty(tok), sk = arcKey; return, end
    sk = sprintf('%s|%s', matlab.lang.makeValidName(tok{1}), tok{2});
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
