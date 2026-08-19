function VC = varComp(u, batchFcn, nBatch, groupFcn, nGroup, opt)
%LSQVCE.VARCOMP  Helmert variance-component estimation, streaming.
%
%   VC = LSQVCE.VARCOMP(u, batchFcn, nBatch, groupFcn, nGroup) estimates a
%   separate variance factor for each observation group, iterating until they
%   agree. Use this when the global test rejects but the residuals are clean
%   and structureless -- i.e. when the problem is the WEIGHTS, not the data.
%
%   This is the right tool when combining heterogeneous observation types whose
%   relative weighting you do not actually know a priori (the usual case when,
%   say, a range-rate instrument is combined with an interferometric one: both
%   vendors quote a noise figure, and at least one of them is wrong).
%
%   INPUTS
%     groupFcn  g = groupFcn(k) returns an nk x 1 vector of group indices in
%               1..nGroup for the observations in batch k.
%     nGroup    number of variance components.
%
%   METHOD
%     Rigorous Helmert requires the nGroup x nGroup matrix
%         H_ij = tr(Qxx N_i Qxx N_j),  N_i = A_i' P_i A_i
%     which needs only u x u accumulations per group -- perfectly streamable.
%     The right-hand side is q_i = v_i' P_i v_i. Solving H s = q gives the
%     variance components. We iterate: rescale P_i by 1/s_i, re-solve, repeat.
%
%     The simplified (Foerstner) estimator
%         s_i = v_i'P_i v_i / r_i,   r_i = n_i - tr(Qxx N_i)
%     is also returned, since it is what most textbooks actually use and it
%     converges to the same fixed point. Rigorous is the default.
%
%   NAME-VALUE
%     'Method'   'rigorous' (default) | 'simplified'
%     'MaxIter'  30
%     'Tol'      1e-6      convergence on max|s_i - 1|
%     'Verbose'  true
%     'MinObs'   10        below this a group's component is not estimable
%     'TraceMethod' 'exact' (default) | 'hutchinson'. Mirrors LSQVCE.VCE and
%                 LSQVCE.REDUNDANCY for entry-point parity. EXACT IS DEFAULT and
%                 authoritative: the redundancy traces tr(Qxx N_g) and the
%                 Helmert matrix tr(Qxx N_i Qxx N_j) are always formed EXACTLY
%                 (O(u^2) per group, streamable) and drive the estimate. With
%                 'hutchinson', those exact numbers are UNCHANGED; the estimator
%                 (Ellmer 2018 6.5.6; see LSQVCE.TRACEESTIMATE) additionally runs
%                 a MACHINERY CHECK on the final normal N -- tr(N^-1 N) = u
%                 exactly -- returned in VC.traceEstimate. This validates the
%                 estimator on THIS adjustment's conditioning WITHOUT yet
%                 replacing the exact per-group r_g. The full per-group
%                 Hutchinson substitution (a per-group scalar r_g in place of
%                 tr(Qxx N_g)) is a separate code path and is HISTORY_v1.md roadmap
%                 item 0; see LSQVCE.REDUNDANCY's help for the structural reason.
%     'Probes'    1x1  number of Hutchinson probes Z (default 100, Ellmer's).
%     'TraceSeed' 1x1  RNG seed for reproducible probes (default []: global stream).
%     'Backend'   "exact" (default) | "hutchinson" (VALUE case-sensitive,
%                 WP-39). Which estimator computes the
%                 VARIANCE COMPONENTS (distinct from 'TraceMethod', which only
%                 adds a machinery check to the exact backend). With 'hutchinson'
%                 the WHOLE estimation is delegated to LSQVCE.VCEHUTCHINSON -- the
%                 matrix-free Ellmer path that never forms the dense inverse Qxx.
%                 Use it only for LARGE u (thousands); below that the exact
%                 backend is faster AND exact AND gives the diagnostics. Strict
%                 Ellmer: the returned VC is then VCEHUTCHINSON's SUBSET struct
%                 (no per-obs fields; adds .redundSE, .estimator='hutchinson').
%                 If 'CheckD' is supplied it FALLS BACK to exact and warns (the
%                 contamination guard needs the per-observation path). This is
%                 HISTORY_v1.md roadmap item 0, tier T2 (an opt-in VCE backend).
%     'Seed'      1x1  RNG seed forwarded to the Hutchinson backend (falls back to
%                 'TraceSeed' if omitted). Ignored by the exact backend.
%     'Reprobe'   1x1  logical, forwarded to the Hutchinson backend (default true).
%
%   OUTPUT VC
%     .s          nGroup x 1  final variance components (converged: all ~1)
%     .scale      nGroup x 1  cumulative weight scaling applied: P_i -> P_i/scale_i
%     .sigma      nGroup x 1  implied sigma per group = sqrt(scale)
%     .nObs       nGroup x 1
%     .redund     nGroup x 1  partial redundancy per group
%     .iter, .converged, .history
%     .traceEstimate  present ONLY with 'TraceMethod','hutchinson': a struct
%                 (.value/.expected/.relErr/.se/.relSE/.probes/.note) from the
%                 machinery check tr(N^-1 N)=u on the final normal. NOT a
%                 per-group redundancy replacement (see 'TraceMethod' above).
%
%   EXAMPLE
%       % Two sensors observe the same arcs. Their quoted noise levels are
%       % almost never both right -- let the DATA say what the ratio is.
%       g    = [ones(n/2,1); 2*ones(n/2,1)];        % sensor 1 | sensor 2
%       gsrc = @(k) g(idx(:,k));
%
%       VC = lsqvce.varComp(u, src, 4, gsrc, 2);
%
%       % VC.sigma        the recovered sigma scale per group
%       % VC.resolution   its 1-sigma sampling noise (~0.5*sqrt(2/r_g))
%       % VC.significant  is the correction BIGGER than that noise?
%
%       % Do not "correct" your weights by 5%% when the resolution is 6%%.
%       if any(VC.significant)
%           lsqvce.util.log.msg('sigma ratio = %%.3f\n', VC.sigma(2)/VC.sigma(1));
%       end
%
%   ERRORS
%     lsqvce:varComp:badGroup      the group key is not a known name or handle
%     lsqvce:varComp:correlated    a batch carries a full covariance, so its
%                                  observations cannot be split between groups
%     lsqvce:varComp:singular      a group has no redundancy of its own, so its
%                                  variance component is not identifiable --
%                                  coarsen the grouping
%     lsqvce:varComp:contaminated  the contamination gauge is past the point
%                                  where the estimate is meaningful
%
%   See also LSQVCE.DIAGNOSE, LSQVCE.COMBINE, LSQVCE.VCE, LSQVCE.REDUNDANCY,
%            LSQVCE.TRACEESTIMATE, LSQVCE.VCEHUTCHINSON.
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) (F-07): typed
%   `arguments` signature, canonical UpperCamelCase option names (legacy
%   casings keep matching; enum VALUES case-sensitive lowercase). With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-40 (fixes F-11): a negative
%   component FREEZES its group for the pass (s = 1) instead of clamping
%   to 1e-6 -- the old clamp scaled the cumulative variance by 1e-6 and
%   fed a 1e6 reweighting oscillation. With assistance from Claude
%   (Anthropic), model Claude Fable 5.
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Updated 2026-07-17 (CEST, UTC+02:00): added 'TraceMethod'/'Probes'/
%   'TraceSeed' (entry-point parity with lsqvce.vce/redundancy; exact stays
%   default, adds a Hutchinson machinery check). Claude (Anthropic), Opus 4.8.
%   Updated 2026-07-17 (CEST, UTC+02:00): added 'backend','exact'|'hutchinson'
%   (item 0 tier T2) delegating the Hutchinson backend to LSQVCE.VCEHUTCHINSON;
%   falls back to exact + warn if 'CheckD' is given. Claude (Anthropic), Opus 4.8.
arguments
    u (1,1) double {mustBeInteger, mustBePositive}
    batchFcn (1,1) function_handle
    nBatch (1,1) double {mustBeInteger, mustBePositive}
    groupFcn (1,1) function_handle
    nGroup (1,1) double {mustBeInteger, mustBePositive}
    opt.Method (1,1) string {mustBeMember(opt.Method, ["rigorous","simplified"])} = "rigorous"
    opt.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 30
    opt.Tol (1,1) double {mustBePositive} = 1e-6
    opt.Verbose (1,1) logical = true
    opt.MinObs (1,1) double {mustBeGreaterThanOrEqual(opt.MinObs, 1)} = 10
    opt.CheckD = struct([])   % struct from lsqvce.diagnose, or [] (both accepted)
    opt.TraceMethod (1,1) string {mustBeMember(opt.TraceMethod, ["exact","hutchinson"])} = "exact"
    opt.Probes (1,1) double {mustBeInteger, mustBePositive} = 100
    opt.TraceSeed double {mustBeScalarOrEmpty} = []
    opt.Backend (1,1) string {mustBeMember(opt.Backend, ["exact","hutchinson"])} = "exact"
    opt.Seed double {mustBeScalarOrEmpty} = []
    opt.Reprobe (1,1) logical = true
end

% WP-39 stage (b): typed `arguments` signature (canonical UpperCamel
% names; NAME matching case-insensitive per MATLAB semantics; enum VALUES
% case-sensitive lowercase). Lowercase mapping keeps the body untouched.
o = struct('method', char(opt.Method), 'maxIter', opt.MaxIter, ...
           'tol', opt.Tol, 'verbose', opt.Verbose, 'minObs', opt.MinObs, ...
           'checkD', opt.CheckD, 'traceMethod', char(opt.TraceMethod), ...
           'probes', opt.Probes, 'traceSeed', opt.TraceSeed, ...
           'backend', char(opt.Backend), 'seed', opt.Seed, ...
           'reprobe', opt.Reprobe);
rigorous = strcmpi(o.method, 'rigorous');
useHutch = strcmpi(o.traceMethod, 'hutchinson');

% ---- BACKEND DISPATCH (item 0 tier T2) ------------------------------------
% If the Hutchinson backend is requested, hand the whole problem to
% LSQVCE.VCEHUTCHINSON (identical signature). Guard: the contamination check
% needs the exact per-obs path, so fall back to exact + warn if 'CheckD' given.
if strcmpi(o.backend, 'hutchinson')
    if ~isempty(o.checkD)
        lsqvce.util.log.warn('lsqvce:varComp:hutchNoCheckD', ...
            ['''backend'',''hutchinson'' cannot honour the ''checkD'' ' ...
             'contamination guard (the estimator has no per-observation ' ...
             'redundancy). Falling back to the EXACT backend so the guard ' ...
             'still protects you. Clean the data first, then re-run with the ' ...
             'Hutchinson backend if you need the speed.']);
    else
        % seed: prefer the vceHutchinson-style 'seed'; fall back to 'TraceSeed'.
        hseed = o.seed;
        if isempty(hseed), hseed = o.traceSeed; end
        VC = lsqvce.vceHutchinson(u, batchFcn, nBatch, groupFcn, nGroup, ...
                'method',  o.method, 'maxIter', o.maxIter, 'tol', o.tol, ...
                'verbose', o.verbose, 'minObs', o.minObs, ...
                'probes',  o.probes, 'seed', hseed, 'reprobe', o.reprobe);
        return
    end
end


% ---------------------------------------------------------------------------
% CONTAMINATION GUARD
% ---------------------------------------------------------------------------
% Variance-component estimation and outlier rejection interact viciously, and
% the toolbox must not let a user run them independently without warning.
%
% If a group contains blunders, they inflate that group's residual sum of
% squares, so Helmert assigns the group a LARGER variance, which DOWN-WEIGHTS
% the whole group, which shrinks its residuals, which hides the blunders
% further. It is the same self-masking pathology as Pope's tau -- one level up,
% and it corrupts the weights of every good observation in the group as well.
%
% The correct order is: clean the data first (or use a robust estimator), THEN
% estimate variance components. So if the caller hands us a diagnostic struct
% showing contamination, refuse.
if ~isempty(o.checkD) && isstruct(o.checkD) && isfield(o.checkD, 'contamLevel')
    switch o.checkD.contamLevel
    case {'refuse', 'breakdown', 'severe'}
        error('lsqvce:varComp:contaminated', ...
            ['varComp: REFUSING TO RUN on contaminated data (estimated %.0f%% ' ...
             'bad observations, contamination level "%s").\n' ...
             'Blunders in a group inflate that group''s variance component, ' ...
             'which DOWN-WEIGHTS the entire group, which shrinks its residuals, ' ...
             'which hides the blunders further -- and corrupts the weights of ' ...
             'every good observation in the group along the way.\n' ...
             'Clean the data first (lsqvce.snoop with ''test'',''mad'', or ' ...
             'lsqvce.robust), THEN estimate variance components. If you have a ' ...
             'reason to proceed anyway, omit ''checkD''.'], ...
             100 * o.checkD.fhat, o.checkD.contamLevel);
    case 'moderate'
        lsqvce.util.log.warn('lsqvce:varComp:someContamination', ...
            ['varComp: an estimated %.0f%% of observations are outlying. The ' ...
             'variance components below will be inflated by those blunders, ' ...
             'and the resulting weights will partly encode them as "noise". ' ...
             'Clean the data first if you can.'], 100 * o.checkD.fhat);
    end
end

scale   = ones(nGroup, 1);       % cumulative multiplier on the COVARIANCE
history = [];
converged = false;

if lsqvce.util.verbosity.progress(o.verbose)
    lsqvce.util.log.msg('\nHelmert variance components (%s), %d groups\n', o.method, nGroup);
    lsqvce.util.log.msg('%s\n', repmat('-', 1, 62));
end

for it = 1:o.maxIter

    % ---- pass 1: accumulate N, b, and per-group N_g with CURRENT weights
    N  = zeros(u);
    b  = zeros(u, 1);
    Ng = repmat({zeros(u)}, nGroup, 1);
    nObs = zeros(nGroup, 1);
    nTot = 0;

    for k = 1:nBatch
        [Ak, lk, pk] = batchFcn(k);
        gk = groupFcn(k);
        lk = lk(:);  gk = gk(:);
        nk = size(Ak, 1);
        if nk == 0, continue; end
        if isempty(pk), pk = ones(nk, 1); end
        pk = pk(:);
        if ~isvector(pk)
            error('lsqvce:varComp:correlated', ...
                  'varComp requires diagonal weights (batch %d gave a matrix).', k);
        end
        if any(gk < 1 | gk > nGroup | gk ~= fix(gk))
            error('lsqvce:varComp:badGroup', ...
                  'groupFcn(%d) returned an index outside 1..%d.', k, nGroup);
        end

        peff = pk ./ scale(gk);          % scaled weights
        sA   = sqrt(peff) .* Ak;
        N    = N + sA' * sA;
        b    = b + sA' * (sqrt(peff) .* lk);

        for g = 1:nGroup
            m = (gk == g);
            if any(m)
                sAg   = sqrt(peff(m)) .* Ak(m, :);
                Ng{g} = Ng{g} + sAg' * sAg;
                nObs(g) = nObs(g) + sum(m);
            end
        end
        nTot = nTot + nk;
    end

    [R, ok] = chol(N);
    if ok ~= 0
        error('lsqvce:varComp:singular', ...
              'Normal matrix became singular at iteration %d.', it);
    end
    x   = R \ (R' \ b);
    Rin = R \ eye(u);
    Qxx = Rin * Rin.';

    % ---- pass 2: q_g = v_g' P_g v_g   with the SAME scaled weights
    q = zeros(nGroup, 1);
    for k = 1:nBatch
        [Ak, lk, pk] = batchFcn(k);
        gk = groupFcn(k);
        lk = lk(:);  gk = gk(:);
        nk = size(Ak, 1);
        if nk == 0, continue; end
        if isempty(pk), pk = ones(nk, 1); end
        pk = pk(:);
        peff = pk ./ scale(gk);
        vk = Ak * x - lk;
        for g = 1:nGroup
            m = (gk == g);
            if any(m)
                q(g) = q(g) + sum(peff(m) .* vk(m).^2);
            end
        end
    end

    % ---- partial redundancies  r_g = n_g - tr(Qxx N_g)
    trQN = zeros(nGroup, 1);
    for g = 1:nGroup
        trQN(g) = sum(sum(Qxx .* Ng{g}));      % tr(Qxx*Ng) without the product
    end
    rg = nObs - trQN;

    estimable = nObs >= o.minObs & rg > 1e-8;

    % ---- solve for the components
    if rigorous
        % H_ij = tr(Qxx Ni Qxx Nj);  diagonal gets the extra n_i - 2 tr(Qxx Ni)
        QN = cell(nGroup, 1);
        for g = 1:nGroup, QN{g} = Qxx * Ng{g}; end
        H = zeros(nGroup);
        for i = 1:nGroup
            for j = 1:nGroup
                H(i, j) = sum(sum(QN{i}.' .* QN{j}));   % tr(QN_i * QN_j)
            end
        end
        H = H + diag(nObs - 2 * trQN);

        s = ones(nGroup, 1);
        ii = find(estimable);
        if numel(ii) >= 1
            Hs = H(ii, ii);
            if rcond(Hs) < 1e-12
                lsqvce.util.log.warn('lsqvce:varComp:illH', ...
                    ['The Helmert matrix is nearly singular: the groups are not ' ...
                     'separately estimable (they see the same parameters in the ' ...
                     'same way). Falling back to the simplified estimator.']);
                s(ii) = q(ii) ./ max(rg(ii), eps);
            else
                s(ii) = Hs \ q(ii);
            end
        end
    else
        s = ones(nGroup, 1);
        ii = find(estimable);
        s(ii) = q(ii) ./ max(rg(ii), eps);
    end

    % Guard: a negative variance component is a real possibility with the
    % rigorous estimator and is meaningless. Clamp and warn.
    bad = s <= 0 & estimable;
    if any(bad)
        lsqvce.util.log.warn('lsqvce:varComp:negative', ...
            ['Group(s) %s produced a NEGATIVE variance component. This is a ' ...
             'known pathology of rigorous Helmert when a group has little ' ...
             'partial redundancy. The component is FROZEN this pass (s = 1); ' ...
             'consider ''method'',''simplified'' or merging the group.'], ...
             mat2str(find(bad).'));
        % WP-40 (fixes F-11): FREEZE, never a small clamp. s multiplies the
        % CUMULATIVE scale, so the old s(bad) = 1e-6 divided that group's
        % variance by 1e6 -- and the next pass reweighted it 1e6 STRONGER:
        % oscillation fuel, not a guard. Freezing keeps the current scale,
        % the warning above flags the pathology, and the group re-enters
        % the estimation on the next pass with unchanged weight.
        s(bad) = 1;
    end
    s(~estimable) = 1;

    scale = scale .* s;
    dev   = max(abs(s(estimable) - 1));
    if isempty(dev), dev = 0; end
    history(end+1, :) = [it, s.']; %#ok<AGROW>

    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg('it %2d   ', it);
        lsqvce.util.log.msg('s=[%s]  max|s-1|=%.3e\n', ...
                strjoin(arrayfun(@(z) sprintf('%.4f', z), s.', 'uni', 0), ' '), dev);
    end

    if dev < o.tol
        converged = true;
        break;
    end
end

if lsqvce.util.verbosity.progress(o.verbose)
    lsqvce.util.log.msg('%s\n', repmat('-', 1, 62));
    if converged
        lsqvce.util.log.msg('Converged in %d iterations.\n', it);
    else
        lsqvce.util.log.msg(['NOT converged in %d iterations (max|s-1| = %.2e).\n' ...
                 'The groups may not be separately estimable.\n'], it, dev);
    end
    lsqvce.util.log.msg('\n  %-6s %8s %10s %12s %12s\n', ...
            'group', 'n', 'redund', 'var scale', 'sigma scale');
    lsqvce.util.log.msg('  %s\n', repmat('-', 1, 52));
    for g = 1:nGroup
        lsqvce.util.log.msg('  %-6d %8d %10.2f %12.5f %12.5f\n', ...
                g, nObs(g), rg(g), scale(g), sqrt(scale(g)));
    end
    lsqvce.util.log.msg('\n');

    % Resolution limit. A variance component estimated from r_g redundant
    % observations has a relative standard error of sqrt(2/r_g); the implied
    % SIGMA scale therefore has about half that. Quoting a weight correction
    % that is smaller than its own sampling noise is meaningless, and users
    % will otherwise happily "correct" their weights by 8% on the strength of
    % pure noise. Verified by Monte Carlo: with r_g = 500 the sigma ratio
    % scatters with sd ~ 0.045, so an apparent 8% discrepancy is a routine
    % 2-sigma fluctuation and means nothing.
    resol = zeros(nGroup, 1);
    for g = 1:nGroup
        if rg(g) > 0
            resol(g) = 0.5 * sqrt(2 / rg(g));   % rel. std err of sigma scale
        else
            resol(g) = Inf;
        end
    end
    lsqvce.util.log.msg('  Resolution limit (1-sigma) on each sigma scale:\n');
    for g = 1:nGroup
        lsqvce.util.log.msg('    group %-3d  +/- %.1f%%   (partial redundancy %.1f)\n', ...
                g, 100 * resol(g), rg(g));
    end
    lsqvce.util.log.msg('\n');

    rat  = sqrt(max(scale) / min(scale));
    noise = sqrt(sum(resol(estimable).^2));      % combined resolution
    if isempty(noise) || ~isfinite(noise), noise = Inf; end

    if (rat - 1) < noise
        lsqvce.util.log.msg(['  The spread between groups (%.1f%%) is SMALLER than the\n' ...
                 '  sampling noise of the estimator itself (%.1f%%). This is\n' ...
                 '  consistent with correct a priori weights. Do NOT apply\n' ...
                 '  these scale factors -- you would be fitting noise.\n\n'], ...
                 100 * (rat - 1), 100 * noise);
    elseif rat > 1.5
        lsqvce.util.log.msg(['  The groups disagree by a factor of %.2f in sigma, well\n' ...
                 '  beyond the %.1f%% resolution of the estimator. The a priori\n' ...
                 '  weighting was genuinely wrong. Apply the scale column to the\n' ...
                 '  covariances and re-run the adjustment.\n\n'], rat, 100 * noise);
    else
        lsqvce.util.log.msg(['  The groups differ by %.1f%%, which exceeds the %.1f%%\n' ...
                 '  resolution limit but is not dramatic. The correction is\n' ...
                 '  real but modest; applying it is defensible.\n\n'], ...
                 100 * (rat - 1), 100 * noise);
    end
end

% Resolution limit is returned regardless of 'Verbose', because a caller
% driving this from a script needs it just as much as a human reading the
% printout: it is the difference between a real weight correction and noise.
resolOut = zeros(nGroup, 1);
for g = 1:nGroup
    if rg(g) > 0
        resolOut(g) = 0.5 * sqrt(2 / rg(g));
    else
        resolOut(g) = Inf;
    end
end

VC = struct('s', s, 'scale', scale, 'sigma', sqrt(scale), ...
            'nObs', nObs, 'redund', rg, 'iter', it, ...
            'converged', converged, 'history', history, ...
            'method', o.method, 'x', x, 'estimable', estimable, ...
            'resolution', resolOut, ...
            'significant', abs(sqrt(scale) - 1) > resolOut);

% ---- optional Hutchinson trace estimator (Ellmer 2018 6.5.6) --------------
% Opt-in, EXPERIMENTAL, entry-point parity with lsqvce.vce / lsqvce.redundancy.
% The exact per-group traces above are UNCHANGED and remain authoritative; this
% only demonstrates the estimator on the final normal against a KNOWN answer:
% tr(N^-1 N) = u exactly, so the estimate must return ~u. It validates the
% estimator machinery end to end on THIS adjustment's normal (probe count,
% conditioning) WITHOUT yet replacing the exact per-group r_g. The full
% per-group substitution needs a group-level normal assembly and is HISTORY_v1.md
% roadmap item 0. R is the Cholesky of the FINAL normal N from the loop above.
if useHutch
    try
        % tr(N^-1 N) = tr(R^-1 (R'R) R^-T)?  We probe N^-1 N = I matrix-free:
        %   N^-1 N V = R \ (R' \ (N V)),  with N = R'*R (upper-triangular R).
        Xfun = @(V) R \ (R' \ (N * V));           % applies N^-1 N = I
        [tEst, tInfo] = lsqvce.traceEstimate(Xfun, 'Size', u, ...
                            'NumProbes', o.probes, 'Seed', o.traceSeed);
        VC.traceEstimate = struct('value', tEst, 'expected', u, ...
            'relErr', abs(tEst - u)/max(u,1), 'se', tInfo.se, ...
            'relSE', tInfo.relse, 'probes', o.probes, ...
            'note', ['machinery check on the final normal: tr(N^-1 N)=u; not ' ...
                     'a per-group redundancy (see help and roadmap item 0)']);
        if lsqvce.util.verbosity.progress(o.verbose)
            lsqvce.util.log.msg(['  Hutchinson check: tr(N^-1 N) est %.2f vs u = %d ' ...
                '(relErr %.2e, relSE %.2e, Z = %d)\n'], ...
                tEst, u, VC.traceEstimate.relErr, tInfo.relse, o.probes);
        end
    catch ME
        lsqvce.util.log.warn('lsqvce:varComp:hutchFailed', ...
            'Hutchinson trace check failed (%s); exact components are unaffected.', ...
            ME.identifier);
        VC.traceEstimate = struct('value', NaN, 'expected', u, ...
            'relErr', NaN, 'se', NaN, 'relSE', NaN, 'probes', o.probes, ...
            'note', 'estimator failed; see warning');
    end
end
end
