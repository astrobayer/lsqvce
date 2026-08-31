function vc = vce(accFcn, opt)
%LSQVCE.ENGINE.VCE  Helmert variance component estimation, per arc.
%
%   vc = LSQVCE.ENGINE.VCE(accFcn, 'Name', value, ...)
%
%   Estimates one variance component per group by iterating
%
%       sigma_g^2 = vtv_g / r_g
%
%   re-whitening the data with the updated sigmas at each step. Converges in
%   2-4 iterations from a badly wrong start (verified: recovers the truth from
%   a 3x wrong initial weight in 2 iterations).
%
%   WHY YOU WANT THIS
%     The relative weighting of observation types -- 'pos' versus 'kbr' -- is
%     NOT known a priori, and getting it wrong biases the GRAVITY FIELD, not
%     merely the formal errors. Guessing sigma and hoping is the single most
%     common way to produce a confident wrong answer.
%
%   INPUT
%     accFcn   @(sigma) -> acc   a function that BUILDS AND RUNS a fresh
%                                accumulator with the given per-group sigmas.
%                                It must return the accumulator (already fully
%                                accumulated, not yet solved). On the FIRST
%                                pass it is called as accFcn([]) -- sigma is
%                                EMPTY, meaning unit weights; the function must
%                                accept that (WP-41: previously undocumented).
%
%   OPTIONS
%   (WP-39 stage (b): canonical UpperCamelCase names below; legacy
%   casings still match -- MATLAB NV semantics -- and enum VALUES are
%   case-sensitive lowercase.)
%     'Groups'      cellstr / numeric vector of group keys, in the order that
%                   sigma is indexed. If omitted, taken from the first pass.
%     'GroupBy'     as LSQVCE.ENGINE.REDUNDANCY (default 'arc')
%     'Sigma0'      initial sigmas (default all ones)
%     'MaxIter'     default 8
%     'Tol'         convergence on max|sigma_new/sigma_old - 1| (default 1e-3)
%     'MinRedundancy'  default 50
%     'FreezeR'  1x1 > 0, default 1.  Estimability threshold of the sigma
%                UPDATE (WP-40, fixes F-10; distinct from 'MinRedundancy',
%                which is LSQVCE.ENGINE.REDUNDANCY's per-observation testability
%                knob): a group with r_g < freezeR keeps its current sigma
%                this pass (upd = 1) and is warned about once -- the old code
%                silently floored r_g to 1, biasing sigma^2 = vtv/1 LOW.
%     'Verbose'     default true
%     'TraceMethod' 'exact' (default) | 'hutchinson'. Both estimate the SAME
%                   quantity, the group redundancy r_g = sum_{i in g}(1-h_i) =
%                   n_g - (local dof absorbed by g) - tr(N^-1 N_g^cc). 'exact'
%                   sums the per-arc leverages. 'hutchinson' streams it
%                   matrix-free (LSQVCE.ENGINE.STREAMCONTRIBUTION .r): the count and
%                   the local term are EXACT from the same sweep, only the
%                   common-block trace is a Monte-Carlo estimate -- an OPT-IN
%                   accelerator for very large nc where the exact O(nc^2)/arc
%                   leverage is dear. vtv (the residual sum of squares) stays
%                   EXACT. NOTE r_g is NOT tr(N^-1 N_g^cc): that trace is only
%                   the absorbed cc-attribution (it sums to nc, while
%                   sum_g r_g = m-u); using it as the Helmert denominator
%                   inflates sigma without bound (finding F-01, fixed WP-32;
%                   evidence validation/v10review.py T2). The final pass also
%                   fills vc.red.traceEstimate. Being stochastic it adds noise
%                   to sigma -- prefer 'exact' whenever affordable; see the
%                   accuracy note in LSQVCE.ENGINE.LEVERAGE.
%     'Probes'      1x1  number of Hutchinson probes Z (default 100, Ellmer's).
%     'TraceSeed'   1x1  RNG seed for reproducible probes (default []).
%
%   OUTPUT
%     vc.sigma      final variance components (one per group)
%     vc.groups     the group keys, aligned with sigma
%     vc.traceSeedUsed  probe seed actually used: the given TraceSeed, or
%                   the seed pinned up front from the private stream when
%                   TraceSeed was empty ('hutchinson', WP-36; replay with
%                   'TraceSeed', vc.traceSeedUsed reproduces vc bit-
%                   identically); [] on the 'exact' path
%     vc.r          final redundancy per group; in 'hutchinson' the streamed
%                   estimate n_g - local_g - tr(N^-1 N_g^cc), consistent with
%                   the loop's denominator (sum_g vc.r ~ m-u on both paths)
%     vc.history    nIter x nGroup, the sigma trajectory
%     vc.converged  logical
%     vc.nIter
%     vc.sol        the final solution (from the converged weighting)
%     vc.red        the final LSQVCE.ENGINE.REDUNDANCY output
%
%   COST
%     Each iteration = one accumulation pass + one redundancy pass, so roughly
%     DOUBLE the cost of a plain solve, times the number of iterations. See
%     LSQVCE.ENGINE.REDUNDANCY for measured timings; at global SH d/o 90 this is hours,
%     which is expected and unavoidable for exact per-arc components.
%
%   PARALLELISM (v4)
%     The redundancy passes inside each Helmert iteration honour the
%     accumulator's 'Parallel' option -- construct with
%     lsqvce.engine.BlockAccumulator(nc, 'Parallel','parfor') and VCE inherits it with no
%     change here. The arc tasks are independent, so this is where the hours
%     go down (target: ~4 h -> 30-40 min per iteration with 16 workers at
%     d/o 90). See LSQVCE.ENGINE.MAPREDUCE for the determinism rules.
%
%   INTERPRETING THE RESULT
%     sigma_g far from 1 means the a-priori covariance for that group was wrong
%     by that factor. A group whose sigma is inflated relative to its neighbours
%     usually contains outliers -- run LSQVCE.ENGINE.OUTLIERS before concluding that its
%     noise is genuinely higher.
%
%   EXAMPLE
%     accFcn = @(sig) buildAccumulator(cfg, sig);   % you write this
%     vc = lsqvce.engine.vce(accFcn, 'GroupBy','arc', 'MaxIter',6);
%
%     bar(vc.sigma);  xlabel('arc');  ylabel('\sigma_g');
%     yline(1, 'r--');
%     lsqvce.util.log.msg('converged in %d iterations\n', vc.nIter);
%
%   ERRORS
%     lsqvce:vce:badAccFcn           the accumulator handle is not callable
%     lsqvce:vce:badSigma0           the starting variance factor is not positive
%     lsqvce:vce:groupMismatch       the group labelling disagrees with the
%                                    accumulated batches
%     lsqvce:vce:hutchGroupMismatch  likewise, in the Hutchinson backend
%
%   See also LSQVCE.ENGINE.REDUNDANCY, LSQVCE.ENGINE.LEVERAGE, LSQVCE.ENGINE.OUTLIERS
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) (F-07): typed
%   `arguments` signature, canonical UpperCamelCase option names (legacy
%   casings keep matching; enum VALUES case-sensitive lowercase). With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-40 (fixes F-10): groups with
%   r_g < 'FreezeR' are FROZEN and warned about instead of silently
%   floored. With assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-36 (fixes F-06): an empty
%   TraceSeed is pinned ONCE up front from the package-private probe
%   stream (global RNG untouched), so the final vc.r probes with the
%   loop's identical Z; new output vc.traceSeedUsed. With assistance from
%   Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-18 (CEST, UTC+02:00), WP-32 (fixes F-01, CRITICAL):
%   'hutchinson' now divides by the streamed REDUNDANCY
%   r_g = n_g - local_g - tr(N^-1 N_g^cc) (LSQVCE.ENGINE.STREAMCONTRIBUTION .r)
%   instead of substituting the bare cc-trace, which sums to nc rather than
%   m-u and inflated sigma without bound; vc.r follows. Evidence:
%   validation/v10review.py T2. With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
arguments
    accFcn (1,1) function_handle
    opt.Groups = []
    opt.GroupBy = 'arc'
    opt.Sigma0 double = []
    opt.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 8
    opt.Tol (1,1) double {mustBePositive} = 1e-3
    opt.MinRedundancy (1,1) double {mustBePositive} = 50
    opt.FreezeR (1,1) double {mustBePositive} = 1
    opt.Verbose (1,1) logical = true
    opt.TraceMethod (1,1) string {mustBeMember(opt.TraceMethod, ["exact","hutchinson"])} = "exact"
    opt.Probes (1,1) double {mustBeInteger, mustBePositive} = 100
    opt.TraceSeed double {mustBeScalarOrEmpty} = []
end

    narginchk(1,inf);
    if ~isa(accFcn, 'function_handle')
        error('lsqvce:vce:badAccFcn', ...
              'accFcn must be a function handle @(sigma) -> accumulator.');
    end

    % WP-39 stage (b): typed `arguments` signature (canonical UpperCamel
    % names; NV NAME matching stays case-insensitive per MATLAB semantics;
    % enum VALUES case-sensitive lowercase). The lowercase o-struct keeps
    % the body untouched.
    o = struct('groups', {opt.Groups}, 'groupBy', opt.GroupBy, ...
               'sigma0', opt.Sigma0, 'maxIter', opt.MaxIter, ...
               'tol', opt.Tol, 'minRedundancy', opt.MinRedundancy, ...
               'freezeR', opt.FreezeR, 'verbose', opt.Verbose, ...
               'traceMethod', char(opt.TraceMethod), 'probes', opt.Probes, ...
               'traceSeed', opt.TraceSeed);
    % WP-36 (fixes F-06): with 'hutchinson' and an EMPTY TraceSeed, pin one
    % concrete seed up front (drawn from the package-private probe stream,
    % global RNG untouched). Every loop pass AND the final vc.r pass then
    % probe with IDENTICAL Z, so vc.r is exactly the converged loop
    % redundancy -- previously the final pass drew fresh probes and vc.r
    % was a different realisation than the denominators the loop used.
    if strcmpi(o.traceMethod, 'hutchinson') && isempty(o.traceSeed)
        [~, o.traceSeed] = probeMatrix(1, 1, 'rademacher', []);
    end

    % ---- pass 0: discover the groups --------------------------------------
    groups = o.groups;
    if isempty(groups)
        acc0 = accFcn([]);                          % [] -> unit weights
        sol0 = acc0.solve();
        red0 = lsqvce.engine.redundancy(acc0, sol0, 'GroupBy', o.groupBy, ...
                              'Collect', false, 'Verbose', false);
        groups = {red0.group.key};
    end
    nG = numel(groups);

    sig = o.sigma0;
    if isempty(sig), sig = ones(1, nG); end
    sig = sig(:).';
    if numel(sig) ~= nG
        error('lsqvce:vce:badSigma0', ...
              'Sigma0 has %d entries but there are %d groups.', numel(sig), nG);
    end

    hist = zeros(o.maxIter, nG);
    conv = false;

    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg('%s%d groups, grouping by ''%s''\n', lsqvce.util.verbosity.banner('Variance components'), nG, groupByLabel(o.groupBy));
    end

    warnedFreeze = false;                        % WP-40: one warning per run
    for it = 1:o.maxIter
        acc = accFcn(sig);
        sol = acc.solve();
        red = lsqvce.engine.redundancy(acc, sol, 'GroupBy', o.groupBy, ...
                             'MinRedundancy', o.minRedundancy, ...
                             'Collect', false, 'Verbose', false);

        % align the groups with sigma
        [tf, loc] = ismember(cellstr(string({red.group.key})), ...
                             cellstr(string(groups)));
        if ~all(tf)
            error('lsqvce:vce:groupMismatch', ...
                ['The groups changed between iterations. The loader must be ' ...
                 'deterministic and must not add or drop arcs.']);
        end

        rg  = zeros(1,nG);  vtv = zeros(1,nG);
        rg(loc)  = [red.group.r];
        vtv(loc) = [red.group.vtv];

        % ---- Level 1.1: streaming Hutchinson REDUNDANCY (WP-32) ----------
        % 'exact' (default) keeps the per-arc redundancies above. 'hutchinson'
        % estimates the SAME r_g matrix-free (LSQVCE.ENGINE.STREAMCONTRIBUTION .r):
        % r_g = n_g - local_g - tr(N^-1 N_g^cc), count and local term exact,
        % only the cc-trace Monte-Carlo -- the payoff regime is very large nc,
        % where the exact O(nc^2)/arc leverage is dear. vtv stays EXACT (the
        % residual pass is unchanged). NEVER substitute the bare .trace here:
        % it sums to nc, not m-u, and the Helmert update then divides by a
        % denominator that is orders of magnitude too small (F-01: measured
        % 8.4x inflation / divergence to 1e23; v10review.py T2b/c).
        if strcmpi(o.traceMethod, 'hutchinson')
            sc = lsqvce.engine.streamContribution(acc, sol, 'GroupBy', o.groupBy, ...
                     'Groups', groups, 'Probes', o.probes, 'Seed', o.traceSeed, ...
                     'WantDiag', false, 'Verbose', false);
            [tf2, loc2] = ismember(cellstr(string({sc.group.key})), ...
                                   cellstr(string(groups)));
            if ~all(tf2)
                error('lsqvce:vce:hutchGroupMismatch', ...
                    ['Streaming Hutchinson returned groups that do not match ' ...
                     'the redundancy grouping -- the loader must be ' ...
                     'deterministic across passes.']);
            end
            rg(loc2) = [sc.group.r];
        end

        % WP-40 (fixes F-10): varComp-style estimability instead of the old
        % SILENT rg = max(rg,1) -- which biased sigma^2 = vtv/1 LOW for any
        % group with true r_g < 1, without a word. A group below freezeR
        % redundant observations has no estimable factor this pass: its
        % sigma is FROZEN (upd = 1) and the group is warned about once.
        frozen = rg(:) < o.freezeR;
        if any(frozen) && ~warnedFreeze
            warnedFreeze = true;
            fk = find(frozen);
            lsqvce.util.log.warn('lsqvce:vce:lowRedundancy', ...
               ['%d group(s) have redundancy below freezeR = %g and their ' ...
                'sigma is FROZEN this pass (not estimable): %s (r_g = %s). ' ...
                'Merge tiny groups or lower freezeR knowingly.'], ...
                numel(fk), o.freezeR, ...
                strjoin(cellstr(string(groups(fk))), ', '), ...
                mat2str(round(rg(fk(:)).', 3)));
        end
        upd = ones(size(rg));
        upd(~frozen) = sqrt(vtv(~frozen) ./ rg(~frozen));
        sig = sig .* upd;
        hist(it,:) = sig;

        dmax = max(abs(upd - 1));
        if lsqvce.util.verbosity.progress(o.verbose)
            lsqvce.util.log.msg(['[vce %2d  ] max|dsigma| = %.3e   sigma in [%.4f, %.4f]  ' ...
                     'min r_g = %.0f\n'], it, dmax, min(sig), max(sig), min(rg));
        end

        if dmax < o.tol
            conv = true;
            break
        end
    end

    hist = hist(1:it, :);

    if ~conv
        lsqvce.util.log.warn('lsqvce:vce:notConverged', ...
            ['VCE did not converge in %d iterations (last change %.3e > tol ' ...
             '%.3e). If sigma is oscillating, some groups probably have too ' ...
             'little redundancy -- pool them.'], o.maxIter, dmax, o.tol);
    end

    % ---- final state, with the converged weighting ------------------------
    acc = accFcn(sig);
    sol = acc.solve();
    red = lsqvce.engine.redundancy(acc, sol, 'GroupBy', o.groupBy, ...
                         'MinRedundancy', o.minRedundancy, 'Verbose', o.verbose, ...
                         'TraceMethod', o.traceMethod, 'Probes', o.probes, ...
                         'TraceSeed', o.traceSeed);

    % Per-group redundancy aligned to the pinned group order. In 'hutchinson'
    % mode vc.r is the STREAMED REDUNDANCY (WP-32: n_g - local_g - trace, the
    % same estimator the loop divides by), so vc.r is consistent across both
    % paths and sum(vc.r) ~ m-u either way.
    rFinal = zeros(1, nG);
    [~, locF] = ismember(cellstr(string({red.group.key})), ...
                         cellstr(string(groups)));
    rFinal(locF) = [red.group.r];
    if strcmpi(o.traceMethod, 'hutchinson')
        scF = lsqvce.engine.streamContribution(acc, sol, 'GroupBy', o.groupBy, ...
                  'Groups', groups, 'Probes', o.probes, 'Seed', o.traceSeed, ...
                  'WantDiag', false, 'Verbose', false);
        [~, locS] = ismember(cellstr(string({scF.group.key})), ...
                             cellstr(string(groups)));
        rFinal(locS) = [scF.group.r];
    end

    vc = struct('sigma', sig, 'groups', {groups}, 'r', rFinal, ...
                'history', hist, 'converged', conv, 'nIter', it, ...
                'sol', sol, 'red', red, 'acc', acc, ...
                'traceSeedUsed', o.traceSeed);   % [] on 'exact'; the pinned
                                                 % probe seed otherwise
                                                 % (WP-36: replay with
                                                 % 'TraceSeed', this value
                                                 % reproduces vc bit-
                                                 % identically)

    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg('%s%s after %d iterations\n', lsqvce.util.verbosity.banner('Variance components'), ...
                ternary(conv, 'converged', 'STOPPED (not converged)'), it);
    end
end

function s = ternary(c, a, b)
%TERNARY  Inline conditional: return a if c is truthy, else b.
%   INPUT  c 1x1 logical/scalar; a,b any (the two branch values)
%   OUTPUT s   a if c else b
    if c, s = a; else, s = b; end
end

function s = groupByLabel(g)
%GROUPBYLABEL  Printable label for a GroupBy that may be a char/string or a
%   function handle (LSQVCE.ENGINE.REDUNDANCY accepts both; string() rejects handles).
    if isa(g, 'function_handle')
        s = func2str(g);
    else
        s = char(string(g));
    end
end
