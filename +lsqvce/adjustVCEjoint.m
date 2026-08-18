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
%
%   OUTPUT
%     out  1x1 struct
%       .C         nLag x 1   the converged covariance function
%       .Shat      nLag x 1   the converged PSD amplitudes
%       .sigma2    M x 1      the converged arc variance factors
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
%   Updated 2026-08-18 (CEST, UTC+02:00), WP-33 (fixes F-03, CRITICAL): step
%   7 subtracts the per-system local fit B*y (recovered once per pass as
%   redundancy does) and the per-lag projector absorbs the local block
%   exactly -- shared implementation jointArcParts + redundancyLags
%   'LocalBbar'/'LocalLl'; evidence v10review.py T3. With assistance from
%   Claude (Anthropic), model Claude Fable 5.

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

    hist = struct('maxDev', [], 'Shat', []);
    converged = false;
    accF = []; solF = []; Shat = [];

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
        [omegaT, sT, nArc] = sweepArcs(acc, sol, sp.UiZ, nLag, C, sigma2, v);
        if nArc == 0
            error('lsqvce:adjustVCEjoint:noArcs', ...
                  'No arcs were enumerated -- nothing to estimate from.');
        end
        if isempty(sigma2) || numel(sigma2) ~= nArc
            sigma2 = ones(nArc, 1);       % eq. 6.5.62 default, now that M known
        end

        % ---- step 9: co-update the PSD AND the arc factors ----------------
        pv = lsqvce.psdVCE(omegaT, sT, C, 'Norm', o.Norm);
        av = lsqvce.arcVarianceFactors(omegaT, sT, sigma2, ...
                                       'Norm', o.Norm, 'MeanType', o.MeanType);
        Shat   = pv.Shat;
        C      = pv.Cnew;
        sigma2 = av.sigma2;

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

    out = struct('C', C, 'Shat', Shat, 'sigma2', sigma2, 'acc', accF, ...
                 'sol', solF, 'converged', converged, 'nIter', numel(hist.maxDev), ...
                 'history', hist);
end

%==========================================================================
function [omegaT, sT, nArc] = sweepArcs(acc, sol, UiZ, nLag, C, sigma2, v)
%SWEEPARCS  Steps 7-8: re-stream every arc, returning the WEIGHTED per-lag
%   temporal vectors Omega^m_n and s^m_n (eqs. 6.5.22-6.5.23), one COLUMN per
%   arc. Mirrors the arc enumeration used by LSQVCE.REDUNDANCY and
%   LSQVCE.STREAMCONTRIBUTION, so grouping and whitening stay consistent.
    ids    = acc.systemIDs();
    x      = sol.x(:);
    omegaT = [];
    sT     = [];
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
                    [om, sn] = oneArc(st, obsType, a, d, uo, UiZ, nLag, x, ...
                                      ySys, Gsys, LlS);
                    if isempty(om), continue, end
                    nArc = nArc + 1;
                    % weight by the CURRENT arc factor and lag covariance
                    if nArc <= numel(sigma2), s2 = sigma2(nArc); else, s2 = 1; end
                    omegaT(:, nArc) = s2 * C .* om;   %#ok<AGROW>
                    sT(:, nArc)     = s2 * C .* sn;   %#ok<AGROW>
                end
            end
        end
        if v
            lsqvce.util.log.msg('%sswept system %s\n', lsqvce.util.verbosity.banner('Joint estimation'), string(st.sysID));
        end
    end
end

%==========================================================================
function [om, sn] = oneArc(st, obsType, arc, day, uo, UiZ, nLag, x, ySys, Gsys, LlS)
%ONEARC  Decorrelate one arc, form its residual, and return the per-lag Omega
%   and s from LSQVCE.REDUNDANCYLAGS. Returns [] for a data gap.
%   WP-33: ySys/Gsys/LlS are the system's local solution, Schur factor and
%   local Cholesky (computed once per pass by sweepArcs); the residual is
%   ebar = W^-T(l - A x - B y) and the projector absorbs the local block
%   exactly (jointArcParts + redundancyLags 'LocalBbar'/'LocalLl').
    uo.obsType = obsType;
    uo.day     = day;
    uo.arc     = arc;
    a = st.dataFcn(day, obsType, arc, uo);
    if isempty(a) || ~isfield(a, 'l') || isempty(a.l)
        om = [];  sn = [];  return
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
            om = [];  sn = [];  return       % non-PD arc: skip rather than fail
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
