function out = adjustVCEjoint(accFcn, varargin)
%LSQVCE.ADJUSTVCEJOINT  Joint (co-estimated) stochastic-model loop.
%
%   out = LSQVCE.ADJUSTVCEJOINT(accFcn)
%   out = LSQVCE.ADJUSTVCEJOINT(accFcn, 'Name', value, ...)
%
%   Runs Ellmer's (2018, §6.5.7) ten-step algorithm, in which ALL variance
%   components are updated TOGETHER from one residual pass and one shared probe:
%   the per-frequency PSD amplitudes Shat^j (the temporal law), the arc-wise
%   variance factors sigma2_m (the within-group variance factors), and -- when a prior is
%   present -- the geophysical factor sigma2_gpm.
%
%   JOINT vs SEPARATE -- WHY THIS IS ITS OWN DRIVER
%     LSQVCE.ADJUSTVCE alternates: fit the temporal whitening, then the
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
%            -> LSQVCE.SHAREDPROBE
%     7    decorrelated residuals ebar = W^-T (l - A x - B y), the FULL
%          Table-6.2 vector (WP-33; the local fit is subtracted per system)
%            -> re-streamed per arc via the store's dataFcn
%     8    Rtilde, etilde, and from them Omega and s per lag
%            -> LSQVCE.REDUNDANCYLAGS
%     9    updated Shat^j, sigma2_m (and sigma2_gpm)
%            -> LSQVCE.PSDVCE and LSQVCE.ARCVARIANCEFACTORS
%     10   if not converged, rebuild Sigma and go to step 2
%
%   *** START MATTERS (verified; read this) ***
%     Ellmer's literal white start (eq. 6.5.62) is a DEGENERATE fixed point of
%     his own weighting -- see the warning in LSQVCE.PSDVCE. This driver
%     therefore defaults 'C0' to a mildly correlated covariance rather than a
%     delta, and passes any degeneracy straight through to psdVCE so the warning
%     is raised rather than hidden.
%
%   INPUT
%     accFcn  function_handle  @(model) -> a SOLVED-READY lsqvce.BlockAccumulator
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
%                             on pass 1 (LSQVCE.REDUNDANCY's 'arc' format) --
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
%     r = lsqvce.adjustVCEjoint(accFcn, 'NumLags', 60, 'MaxIter', 8);
%     semilogy(r.Shat);          % the estimated noise PSD
%     plot(r.sigma2);            % the arc-wise weights
%
%   ERRORS
%     lsqvce:adjustVCEjoint:needCfcn  the joint mode needs a cross-covariance
%                                     function and none was supplied
%     lsqvce:adjustVCEjoint:badC0     the starting cross-covariance is malformed
%     lsqvce:adjustVCEjoint:notPD     the joint covariance is not positive
%                                     definite, so it cannot be factored
%     lsqvce:adjustVCEjoint:noArcs    a system declares no arcs
%
%   See also LSQVCE.ADJUSTVCE (the separable default), LSQVCE.PSDVCE,
%            LSQVCE.ARCVARIANCEFACTORS, LSQVCE.REDUNDANCYLAGS,
%            LSQVCE.SHAREDPROBE
%
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
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

    p = inputParser;  p.FunctionName = 'lsqvce.adjustVCEjoint';
    p.addParameter('NumLags' , 60, @(v) isscalar(v) && v >= 2 && v == floor(v));
    p.addParameter('C0'      , []);
    p.addParameter('Sigma2_0', []);
    p.addParameter('MaxIter' , 10, @(v) isscalar(v) && v >= 1);
    p.addParameter('Tol'     , 1e-2, @(v) isscalar(v) && v > 0);
    p.addParameter('Probes'  , 100, @(v) isscalar(v) && v >= 1);
    p.addParameter('Seed'    , []);
    p.addParameter('Norm'    , 'ortho');
    p.addParameter('MeanType', 'robust');
    p.addParameter('StepSchedule', 'none', ...
        @(s) any(strcmpi(s, {'none','decay'})));
    p.addParameter('AvgTail', 0, @(v) isscalar(v) && v >= 0 && v == floor(v));
    p.addParameter('Verbose' , false, @(v) isscalar(v));
    p.parse(varargin{:});
    o = p.Results;
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

    hist = struct('maxDev', [], 'Shat', [], 'skips', {cell(0,3)});
    converged = false;
    accF = []; solF = []; Shat = [];
    pinnedKeys = {};                     % WP-34: pass-1 arc identity
    s2Map = [];                          % keyed factors; built on pass 1

    shatTail = {};                               % WP-40 r2: AvgTail buffer
    for it = 1:o.MaxIter
        % ---- steps 1-5: build under the current model, accumulate, solve --
        model = struct('C', C, 'sigma2', sigma2, 'iter', it);
        acc   = accFcn(model);
        sol   = acc.solve();
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
        sp = lsqvce.sharedProbe(U, 'numProbes', o.Probes, 'seed', o.Seed);

        % ---- steps 7-8: re-stream arcs, per-lag Omega and s ---------------
        % WP-34 (fixes F-02): arcs are KEYED, the factor map is pinned on the
        % PASS-1 ENUMERATION, and a skipped arc keeps its previous factor --
        % one skip can no longer shift every later arc's factor by position.
        [omegaT, sT, keys, keysAll, skips] = ...
            sweepArcs(acc, sol, sp.UiZ, nLag, C, s2Map, v);
        nArc = numel(keys);
        if nArc == 0
            error('lsqvce:adjustVCEjoint:noArcs', ...
                  'No arcs were enumerated -- nothing to estimate from.');
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
        pv = lsqvce.psdVCE(omegaT, sT, C, 'Norm', o.Norm, 'Gamma', gammaIt);
        av = lsqvce.arcVarianceFactors(omegaT, sT, s2Cols, ...
                                       'Norm', o.Norm, 'MeanType', o.MeanType);
        Shat = pv.Shat;
        C    = pv.Cnew;
        % WP-40 r2 (F-15): Polyak tail averaging of the RETURNED spectrum
        % (log domain, last AvgTail passes). The loop keeps iterating on the
        % raw C; only the final Shat/C are averaged -- measured to add ~0.02
        % corr on top of the decay schedule (and to do NOTHING without it:
        % averaging cannot rescue an undamped collapse; mirror grid).
        if o.AvgTail > 0
            shatTail{end+1} = log(max(Shat, realmin));       %#ok<AGROW>
        end
        for kk = 1:nArc
            s2Map(keys{kk}) = av.sigma2(kk);         % keyed write-back
        end
        sigma2 = cellfun(@(k) s2Map(k), pinnedKeys(:));   % pinned order

        hist.maxDev(it,1)  = pv.maxDev;              %#ok<AGROW>
        hist.Shat(it,:)    = Shat(:).';              %#ok<AGROW>

        if v
            lsqvce.util.log.msg(['  [joint %2d] max|alpha2-1| = %.4e   ' ...
                'arcs = %d   gammaBar = %.4f\n'], it, pv.maxDev, nArc, av.gammaBar);
        end

        % ---- step 10: converged? -----------------------------------------
        if pv.maxDev <= o.Tol
            converged = true;
            break
        end
    end

    if ~converged
        lsqvce.util.log.warn('lsqvce:adjustVCEjoint:notConverged', ...
           ['The joint stochastic model did not settle in %d outer iterations ' ...
            '(last max|alpha2-1| = %.3e > tol %.1e). Inspect out.history.Shat: ' ...
            'a PSD that keeps moving usually means too few arcs per frequency, ' ...
            'or a starting covariance far from the data.'], ...
            o.MaxIter, hist.maxDev(end), o.Tol);
    end

    if o.AvgTail > 0 && numel(shatTail) >= o.AvgTail
        Shat = exp(mean(cat(2, shatTail{end-o.AvgTail+1:end}), 2));
        C    = lsqvce.dctI(Shat, 'Norm', o.Norm, 'Inverse', true);  % psdVCE's exact inverse
    end
    out = struct('C', C, 'Shat', Shat, 'sigma2', sigma2, 'acc', accF, ...
                 'sol', solF, 'converged', converged, 'nIter', numel(hist.maxDev), ...
                 'history', hist, ...
                 'arcKeys', {pinnedKeys(:)}, ...
                 'skips', {hist.skips});
end

%==========================================================================
function [omegaT, sT, keys, keysAll, skips] = sweepArcs(acc, sol, UiZ, nLag, C, s2Map, v)
%SWEEPARCS  Steps 7-8: re-stream every arc, returning the WEIGHTED per-lag
%   temporal vectors Omega^m_n and s^m_n (eqs. 6.5.22-6.5.23), one COLUMN per
%   arc. Mirrors the arc enumeration used by LSQVCE.REDUNDANCY and
%   LSQVCE.STREAMCONTRIBUTION, so grouping and whitening stay consistent.
%   WP-34 (fixes F-02): every arc is KEYED s<sys>_<obsType>_d<day>_a<arc>
%   (LSQVCE.REDUNDANCY's 'arc' key format, verbatim). OUTPUTS:
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

    for k = 1:numel(ids)
        st = acc.getStore(ids(k));
        uo = st.uo;
        % WP-33: recover this system's LOCAL solution and Schur pieces ONCE
        % per pass, exactly as LSQVCE.REDUNDANCY does -- the residual is the
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
                    [om, sn, why] = oneArc(st, obsType, a, d, uo, UiZ, ...
                                           nLag, x, ySys, Gsys, LlS);
                    if isempty(om)
                        skips(end+1, :) = {key, why};        %#ok<AGROW>
                        continue
                    end
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
function [om, sn, why] = oneArc(st, obsType, arc, day, uo, UiZ, nLag, x, ySys, Gsys, LlS)
%ONEARC  Decorrelate one arc, form its residual, and return the per-lag Omega
%   and s from LSQVCE.REDUNDANCYLAGS. WP-34: a skipped arc REPORTS why --
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
    a = st.dataFcn(day, obsType, arc, uo);
    if isempty(a) || ~isfield(a, 'l') || isempty(a.l)
        om = [];  sn = [];  why = 'gap';  return
    end

    A = a.A;  l = a.l(:);
    m = size(A, 1);

    % The arc's whitener. A loader may hand back either an explicit covariance
    % (Cfcn) or a whitening operator (Wfcn); redundancyLags needs the explicit
    % triangular W with Sigma = W'W, so build it from Cfcn when available.
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
              'provide Cfcn. Most LSQVCE.COVMODEL types do; a hand-written ' ...
              'loader may need one added.'], string(st.sysID), string(obsType), ...
              string(day), string(arc));
    end

    % Step 7 with the FULL Table-6.2 vector: subtract B*y, reduce the common
    % design (Atil = A - B G) so the shared probe matches the factored N, and
    % hand the exact local block to the projector (WP-33, fixes F-03).
    B = zeros(m, 0);
    if isfield(a, 'B') && ~isempty(a.B), B = a.B; end
    [AtilBar, Bbar, ebar] = jointArcParts(A, B, l, W, x(1:size(A,2)), ...
                                          ySys, Gsys);
    r  = lsqvce.redundancyLags(AtilBar, W, UiZ, ebar, min(nLag, m), ...
                               'LocalBbar', Bbar, 'LocalLl', LlS);
    om = padTo(r.omega, nLag);
    sn = padTo(r.s    , nLag);
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
