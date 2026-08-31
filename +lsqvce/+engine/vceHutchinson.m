function VC = vceHutchinson(u, batchFcn, nBatch, groupFcn, nGroup, opt)
%LSQVCE.ENGINE.VCEHUTCHINSON  Helmert variance components via the Ellmer (2018) 6.5.6
%   stochastic (Hutchinson) trace estimator -- a SPEED-ORIENTED sibling of
%   LSQVCE.ENGINE.VARCOMP for LARGE common blocks.
%
%   VC = LSQVCE.ENGINE.VCEHUTCHINSON(u, batchFcn, nBatch, groupFcn, nGroup) estimates a
%   separate variance factor per observation group, iterating until they agree,
%   exactly like LSQVCE.ENGINE.VARCOMP -- but it NEVER forms the dense inverse Qxx = N^-1
%   and NEVER forms the O(u^3) products Qxx*N_g. Instead it estimates every trace
%   it needs with Z random probes against the Cholesky factor of N (Ellmer eq.
%   6.5.40-6.5.48). It uses the SAME stochastic trace estimator as
%   LSQVCE.ENGINE.TRACEESTIMATE, but INLINES it so the shared factor Ubar = U^-1 raw
%   (eq. 6.5.48) and the per-group blocks W_g are formed once and reused across
%   every group and group pair -- the reuse that makes it fast. This is the "full
%   per-group Hutchinson" path (HISTORY_v1.md roadmap item 0), built as a SEPARATE
%   ENTRY POINT so the exact LSQVCE.ENGINE.VARCOMP / LSQVCE.ENGINE.ADJUSTVCE paths are untouched.
%
%   WHEN TO USE IT (and when NOT to)
%     Use it ONLY when u (the common-parameter count) is large -- thousands of
%     unknowns, e.g. a high-degree gravity field. The estimator's cost is
%     O(u^3 [one Cholesky] + G*Z*u^2), versus VARCOMP's O(G*u^3); the crossover is
%     u ~ Z, so the win is roughly u/Z and only materialises for u in the
%     thousands. For u up to a few thousand, LSQVCE.ENGINE.VARCOMP is BOTH faster AND
%     exact AND gives the per-observation diagnostics -- prefer it. This function
%     also trades accuracy for speed: every trace carries a Monte-Carlo sampling
%     error ~ 1/sqrt(Z), worse on low-redundancy or near-rank-deficient groups.
%
%   STRICT ELLMER -- WHAT THIS DELIBERATELY DOES *NOT* PRODUCE
%     Ellmer's estimator returns only PER-GROUP SCALAR traces. It therefore
%     CANNOT return the per-observation redundancies r_i = 1 - h_i that the exact
%     path derives, so this function does NOT provide, and must not be used for:
%       * per-observation redundancy r_i,
%       * data snooping / outlier flags (LSQVCE.ENGINE.SNOOP / LSQVCE.ENGINE.OUTLIERS),
%       * the (w,r) reliability / MDB numbers,
%       * the integrity guard sum(r_i) = m - u.
%     If you need any of those, use LSQVCE.ENGINE.VARCOMP (exact). This function returns
%     ONLY the per-group variance components, their scales, and the per-group
%     partial redundancies r_g (estimated, with a sampling error).
%
%   INPUTS  (identical to LSQVCE.ENGINE.VARCOMP)
%     u         1x1  number of common parameters (the normal is u x u).
%     batchFcn  [A,l,p] = batchFcn(k)  returns batch k: design A (nk x u),
%               observations l (nk x 1), diagonal weights p (nk x 1 or []).
%     nBatch    1x1  number of batches to stream over.
%     groupFcn  g = groupFcn(k)  returns an nk x 1 vector of group indices in
%               1..nGroup for the observations in batch k.
%     nGroup    1x1  number of variance components.
%
%   NAME-VALUE
%     'method'    'rigorous' (default) | 'simplified'. As LSQVCE.ENGINE.VARCOMP: rigorous
%                 solves H s = q with the estimated Helmert matrix; simplified
%                 uses s_g = q_g / r_g. Rigorous needs the estimated H_ij.
%     'maxIter'   30
%     'tol'       1e-6      convergence on max|s_g - 1|
%     'verbose'   true
%     'minObs'    10        below this a group's component is not estimable
%     'probes'    1x1  Z, number of Hutchinson probes (default 100, Ellmer's).
%                 Larger Z reduces the sampling error ~ 1/sqrt(Z) at linear cost.
%     'seed'      1x1  RNG seed for reproducible probes (default []: global stream).
%     'reprobe'   1x1  logical, draw a FRESH probe block each iteration (default
%                 true). Fresh probes avoid a fixed-probe bias correlating across
%                 iterations; set false to reuse one block (cheaper, deterministic
%                 given 'seed', but the components then share one probe draw).
%
%   OUTPUT VC  (a SUBSET of LSQVCE.ENGINE.VARCOMP's struct -- no per-obs fields)
%     .s          nGroup x 1  final variance components (converged: all ~1)
%     .scale      nGroup x 1  cumulative weight scaling: P_g -> P_g / scale_g
%     .sigma      nGroup x 1  implied sigma per group = sqrt(scale)
%     .nObs       nGroup x 1
%     .redund     nGroup x 1  ESTIMATED partial redundancy r_g (sampling error)
%     .redundSE   nGroup x 1  1-sigma Monte-Carlo standard error on each r_g
%     .iter, .converged, .history
%     .method     the estimator used
%     .estimable  nGroup x 1  logical
%     .resolution nGroup x 1  0.5*sqrt(2/r_g), the sigma-scale resolution limit
%     .probes     1x1  Z used
%     .estimator  'hutchinson'  (marks this as the estimated, not exact, path)
%
%   METHOD (matrix-free, Ellmer eq. 6.5.48 shared factor)
%     With N = U'U (upper-triangular U) and Ubar = U^-1 * raw for a raw u x Z
%     Rademacher block:
%       r_g   = n_g - tr(N^-1 N_g),   tr(N^-1 N_g) ~ tr(Ubar' N_g Ubar)
%       H_ij  = tr(N^-1 N_i N^-1 N_j) ~ sum(sum(W_i .* W_j)),  W_g = U^-T (N_g Ubar)
%     Ubar and each W_g are formed ONCE per iteration and reused for every group
%     pair, so no u x u inverse and no O(u^3) product ever appears.
%
%   ERRORS
%     lsqvce:vceHutchinson:badBatchFcn  the batch handle is not callable
%     lsqvce:vceHutchinson:badGroup     the group labelling is malformed
%     lsqvce:vceHutchinson:badGroupFcn  the grouping handle is not callable
%     lsqvce:vceHutchinson:correlated   a batch carries a full covariance and
%                                       cannot be split between groups
%     lsqvce:vceHutchinson:singular     a group has no redundancy of its own
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   See also LSQVCE.ENGINE.VARCOMP, LSQVCE.ENGINE.VCE, LSQVCE.ENGINE.TRACEESTIMATE, LSQVCE.ENGINE.ADJUSTVCE,
%            LSQVCE.ENGINE.REDUNDANCY.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-17 (CEST, UTC+02:00). Implements HISTORY_v1.md roadmap item 0, tier T1
%   (a separate strict-Ellmer entry point; the streaming block-arrow variant
%   remains future work).
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

arguments
    u
    batchFcn
    nBatch
    groupFcn
    nGroup
    opt.Method (1,1) string {mustBeMember(opt.Method, ["rigorous","simplified"])} = "rigorous"
    opt.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 30
    opt.Tol (1,1) double {mustBePositive} = 1e-6
    opt.Verbose (1,1) logical = true
    opt.MinObs (1,1) double {mustBeGreaterThanOrEqual(opt.MinObs,1)} = 10
    opt.Probes (1,1) double {mustBeInteger, mustBePositive} = 100
    opt.Seed double {mustBeScalarOrEmpty} = []
    opt.Reprobe (1,1) logical = true
end

    narginchk(5, inf);
    validateattributes(u,      {'numeric'}, {'scalar','integer','positive'}, mfilename, 'u', 1);
    validateattributes(nBatch, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'nBatch', 3);
    validateattributes(nGroup, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'nGroup', 5);
    if ~isa(batchFcn, 'function_handle')
        error('lsqvce:vceHutchinson:badBatchFcn', 'batchFcn must be a function handle.');
    end
    if ~isa(groupFcn, 'function_handle')
        error('lsqvce:vceHutchinson:badGroupFcn', 'groupFcn must be a function handle.');
    end

    % WP-39 stage (b) b2: typed `arguments` signature (canonical
    % UpperCamelCase; legacy casings keep matching -- MATLAB NV
    % semantics; enum VALUES case-sensitive lowercase).
    o = struct('method', char(opt.Method), 'maxIter', opt.MaxIter, 'tol', opt.Tol, 'verbose', opt.Verbose, 'minObs', opt.MinObs, 'probes', opt.Probes, 'seed', opt.Seed, 'reprobe', opt.Reprobe);
    rigorous = strcmpi(o.method, 'rigorous');
    Z = o.probes;

    scale   = ones(nGroup, 1);
    history = [];
    converged = false;

    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg('\nHelmert variance components -- HUTCHINSON estimator (%s), %d groups, Z=%d\n', ...
            o.method, nGroup, Z);
        lsqvce.util.log.msg('%s\n', repmat('-', 1, 66));
        lsqvce.util.log.msg(['NOTE: estimated (not exact); no per-obs r_i, snooping, or ' ...
            'sum(r_i)=m-u guard.\n']);
    end

    % A per-iteration RNG so 'seed' gives reproducible probe draws while still
    % allowing a fresh block each iteration (seed + iteration index).
    rg       = zeros(nGroup, 1);
    rgSE     = zeros(nGroup, 1);
    nObs     = zeros(nGroup, 1);

    for it = 1:o.maxIter

        % ---- pass 1: accumulate N, b, and per-group N_g with CURRENT weights
        % (identical to LSQVCE.ENGINE.VARCOMP: the per-group normals are assembled by
        %  re-streaming, so no accumulator core change is needed.)
        N  = zeros(u);
        b  = zeros(u, 1);
        Ng = repmat({zeros(u)}, nGroup, 1);
        nObs = zeros(nGroup, 1);

        for k = 1:nBatch
            [Ak, lk, pk] = batchFcn(k);
            gk = groupFcn(k);
            lk = lk(:);  gk = gk(:);
            nk = size(Ak, 1);
            if nk == 0, continue; end
            if isempty(pk), pk = ones(nk, 1); end
            pk = pk(:);
            if ~isvector(pk)
                error('lsqvce:vceHutchinson:correlated', ...
                      'vceHutchinson requires diagonal weights (batch %d gave a matrix).', k);
            end
            if any(gk < 1 | gk > nGroup | gk ~= fix(gk))
                error('lsqvce:vceHutchinson:badGroup', ...
                      'groupFcn(%d) returned an index outside 1..%d.', k, nGroup);
            end
            peff = pk ./ scale(gk);
            sA   = sqrt(peff) .* Ak;
            N    = N + sA' * sA;
            b    = b + sA' * (sqrt(peff) .* lk);
            for g = 1:nGroup
                mask = (gk == g);
                if any(mask)
                    sAg   = sqrt(peff(mask)) .* Ak(mask, :);
                    Ng{g} = Ng{g} + sAg' * sAg;
                    nObs(g) = nObs(g) + sum(mask);
                end
            end
        end

        [U, okc] = chol(N);            % N = U' * U, upper-triangular U
        if okc ~= 0
            error('lsqvce:vceHutchinson:singular', ...
                  'Normal matrix became singular at iteration %d.', it);
        end
        x = U \ (U' \ b);

        % ---- shared probe factor (Ellmer eq. 6.5.48): built ONCE, reused for
        %      every group and every group pair this iteration.
        raw  = drawProbes(u, Z, o.seed, o.reprobe, it);   % u x Z Rademacher
        Zbar = raw / sqrt(Z);          % eq. 6.5.43 normalisation
        Ubar = U \ Zbar;               % = U^-1 Zbar  (one triangular solve block)

        % ---- pass 2: q_g = v_g' P_g v_g  (same scaled weights, same as varComp)
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
                mask = (gk == g);
                if any(mask)
                    q(g) = q(g) + sum(peff(mask) .* vk(mask).^2);
                end
            end
        end

        % ---- per-group W_g = U^-T (N_g Ubar)  [u x Z], formed once per group.
        %      tr(N^-1 N_g)      ~ sum(sum(Zbar .* (N_g Ubar)))   [= Ubar'N_g Ubar]
        %      tr(N^-1 N_i N^-1 N_j) ~ sum(sum(W_i .* W_j))
        W    = cell(nGroup, 1);
        trQN = zeros(nGroup, 1);
        trSE = zeros(nGroup, 1);
        for g = 1:nGroup
            NgU   = Ng{g} * Ubar;              % u x Z
            % redundancy trace tr(N^-1 N_g) = tr(Ubar' N_g Ubar): the per-probe
            % quadratic forms are the column sums of Ubar .* (N_g Ubar); their
            % MEAN is the trace estimate (Zbar already carried the 1/Z, and Ubar
            % = U^-1 Zbar carries it too, so sum over columns IS the estimate).
            perCol   = sum(Ubar .* NgU, 1);    % 1 x Z, each column's z'(.)z / Z
            trQN(g)  = sum(perCol);
            if Z >= 2
                % un-normalised per-probe forms (multiply the 1/Z back out) give
                % the sample scatter -> standard error of the mean.
                trSE(g) = std(perCol * Z) / sqrt(Z);
            else
                trSE(g) = NaN;
            end
            if rigorous
                W{g} = U' \ NgU;              % U^-T (N_g Ubar) = W_g
            end
        end
        rg   = nObs - trQN;
        rgSE = trSE;                          % SE of r_g equals SE of the trace

        estimable = nObs >= o.minObs & rg > 1e-8;

        % ---- solve for the components
        if rigorous
            % H_ij = tr(N^-1 N_i N^-1 N_j) ~ sum(sum(W_i .* W_j)); diagonal adds
            % the exact n_i - 2 tr(N^-1 N_i) term (same closed form as varComp).
            H = zeros(nGroup);
            for i = 1:nGroup
                for j = i:nGroup
                    hij = sum(sum(W{i} .* W{j}));
                    H(i, j) = hij;
                    H(j, i) = hij;
                end
            end
            H = H + diag(nObs - 2 * trQN);
            s = ones(nGroup, 1);
            ii = find(estimable);
            if ~isempty(ii)
                Hs = H(ii, ii);
                if rcond(Hs) < 1e-12
                    lsqvce.util.log.warn('lsqvce:vceHutchinson:illH', ...
                        ['The estimated Helmert matrix is nearly singular: the groups ' ...
                         'are not separately estimable (or Z is too small to resolve ' ...
                         'them). Falling back to the simplified estimator this step.']);
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

        % Negative-component guard (same pathology and remedy as varComp).
        bad = s <= 0 & estimable;
        if any(bad)
            lsqvce.util.log.warn('lsqvce:vceHutchinson:negative', ...
                ['Group(s) %s produced a NEGATIVE variance component (rigorous ' ...
                 'Helmert with little partial redundancy, possibly aggravated by ' ...
                 'probe noise). Clamped; consider ''method'',''simplified'', more ' ...
                 '''probes'', or merging the group.'], mat2str(find(bad).'));
            s(bad) = 1e-6;
        end
        s(~estimable) = 1;

        scale = scale .* s;
        dev   = max(abs(s(estimable) - 1));
        if isempty(dev), dev = 0; end
        history(end+1, :) = [it, s.']; %#ok<AGROW>

        if lsqvce.util.verbosity.progress(o.verbose)
            lsqvce.util.log.msg('it %2d   s=[%s]  max|s-1|=%.3e\n', it, ...
                strjoin(arrayfun(@(z) sprintf('%.4f', z), s.', 'uni', 0), ' '), dev);
        end

        if dev < o.tol
            converged = true;
            break;
        end
    end

    % ---- resolution limit (same definition as varComp) + report -----------
    resolOut = zeros(nGroup, 1);
    for g = 1:nGroup
        if rg(g) > 0
            resolOut(g) = 0.5 * sqrt(2 / rg(g));
        else
            resolOut(g) = Inf;
        end
    end

    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg('%s\n', repmat('-', 1, 66));
        if converged
            lsqvce.util.log.msg('Converged in %d iterations.\n', it);
        else
            lsqvce.util.log.msg(['NOT converged in %d iterations (max|s-1| = %.2e). The ' ...
                'groups may not be separately estimable, or Z may be too small.\n'], it, dev);
        end
        lsqvce.util.log.msg('\n  %-6s %8s %12s %8s %12s %12s\n', ...
            'group', 'n', 'redund(est)', '+/-', 'var scale', 'sigma scale');
        lsqvce.util.log.msg('  %s\n', repmat('-', 1, 62));
        for g = 1:nGroup
            lsqvce.util.log.msg('  %-6d %8d %12.2f %8.2f %12.5f %12.5f\n', ...
                g, nObs(g), rg(g), rgSE(g), scale(g), sqrt(scale(g)));
        end
        lsqvce.util.log.msg('\n');
    end

    VC = struct('s', s, 'scale', scale, 'sigma', sqrt(scale), ...
                'nObs', nObs, 'redund', rg, 'redundSE', rgSE, ...
                'iter', it, 'converged', converged, 'history', history, ...
                'method', o.method, 'estimable', estimable, ...
                'resolution', resolOut, 'probes', Z, 'estimator', 'hutchinson');
end

%==========================================================================
function raw = drawProbes(u, Z, seed, reprobe, it)
%DRAWPROBES  A u x Z Rademacher (+/-1) probe block. If REPROBE, a fresh block is
%   drawn each iteration; SEED (if given) makes the whole sequence reproducible
%   by offsetting the stream per iteration. If ~REPROBE, the same block (given
%   SEED) is returned every iteration.
    if isempty(seed)
        rs = RandStream.getGlobalStream();
        raw = 2 * randi(rs, [0 1], u, Z) - 1;
        return
    end
    if reprobe
        rs = RandStream('twister', 'Seed', seed + it);   % per-iteration offset
    else
        rs = RandStream('twister', 'Seed', seed);
    end
    raw = 2 * randi(rs, [0 1], u, Z) - 1;
end
