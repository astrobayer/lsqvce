function R = adjustRobust(systems, days, opt)
%LSQVCE.ENGINE.ADJUSTROBUST  IRLS (Huber / IGG3 / bisquare) on the block path.
%   R = LSQVCE.ENGINE.ADJUSTROBUST(systems, days, 'Name', value, ...)
%
%   Robust M-estimation as a POLICY over the one accumulator (R5 M2, WP-81
%   round 2b; finding H-02): every pass accumulates the same systems with
%   per-observation weight factors riding INSIDE each arc's whitener
%   (lsqvce.engine.makeLoader, opt.rowWeights), solves, sweeps the residuals
%   once, and turns the law-whitened residuals into the next factors through
%   lsqvce.kernel.robustWeights. Because the factors live in the whitener,
%   the accumulation, the residual sweep and the final diagnostics all see
%   the same reweighted problem -- one definition, no second code path.
%   Works for every fit kind: the streaming kind directly, the in-memory and
%   batched kinds through fit's block adapter.
%
%   The result is NOT a least-squares result: its covariance is that of the
%   reweighted problem and is scaled by the reweighted sigma0 (report() says
%   so). The pass record is R.history (scale, max|dx|, counts); the final
%   factors are R.history(end).factor, one per observation in sweep order,
%   with R.history(end).keys naming the (sys, obsType, day, arc) rows. They
%   are the weights the RETURNED solve used (WP-132): R.x is the weighted
%   GLS solution on them, at every stop reason.
%
%   INPUT
%     systems  1xM cell | struct  as for LSQVCE.ENGINE.ADJUSTBLOCK
%     days     1xD double          the calendar
%   NAME-VALUE
%     'Backend'    "qr" | "normal"                          default "qr"
%     'PsdMode'    "toeplitz" | "ar" | "arma" | ""           how a declared 'psd'
%                  law whitens (default "": the covariance model's AR(30));
%                  LSQVCE.FIT passes Options.Noise.Whitener
%     'Influence'  "huber" | "igg3" | "bisquare"            default "huber"
%     'Tuning'     [] | constants of the influence          default [] (its default)
%     'MaxIter'    passes (>= 1)                            default 30
%     'Tol'        stop when max|dx| < Tol * max(max|x|, max sigma_x)
%                  -- relative to the problem's own scale     default 1e-8
%     'Parallel', 'DistributeDays', 'ParallelThreshold'     as adjustBlock
%     'Verbose'    logical                                  default false
%     'Provenance' struct   resolved user options for Results.meta (fit)
%     'Prior'      [] | lsqvce.engine.CommonPrior on the common block   default []
%     'Diagnostics' cell   Name-Value pairs for lsqvce.engine.verdict.diagnose default {}
%     'OnPass'     [] | @(ctx)  the per-pass inspection hook (lsqvce.util.onPass:
%                return false to stop after this pass)           default []
%     'OnPassContext' struct  fields merged into every ctx (lsqvce.fit
%                passes .data and .processingOptions)            default struct()
%   OUTPUT
%     R  1x1 lsqvce.Results   stopReason "converged" | "maxPasses"; history(k):
%        .pass .scale .dx .nDown (0 < factor < 0.5) .nZero (factor == 0)
%        (scale, nDown, nZero: what pass k derived from its residuals for
%        the next pass) .factor (last pass: n x 1, the weights the returned
%        solve used) .keys (last pass: n x 4 cell)
%
%   ERRORS
%     lsqvce:adjustRobust:collapsed   fewer live observations than nc + 1
%                                     after reweighting -- ALMOST EVERY row
%                                     rejected. MEASURED (R6 WP-110), this is
%                                     reached two ways, neither an ordinary
%                                     outlier problem: residuals offset
%                                     BODILY from zero by several MAD (a bias
%                                     the model did not absorb, a datum
%                                     defect -- the influence standardises
%                                     vw/s, not (vw - med)/s), or
%                                     contamination so near total that fewer
%                                     than nc + 1 clean rows are left. It is
%                                     NOT a 50 % contamination detector:
%                                     nAlive tracks the CLEAN count, so 60 %
%                                     contamination of 1000 rows still leaves
%                                     400 alive. Check the model, not the data.
%     + everything adjustBlock / makeLoader / robustWeights raise
%   WARNINGS
%     lsqvce:adjustRobust:notConverged  MaxIter reached; the last pass is
%                                       returned with stopReason "maxPasses"
%     lsqvce:robust:atBreakdown         the influence zeroed > 45 % of the
%                                       rows: at the breakdown point (once
%                                       per run; diagnostic only)
%
%   ACCURACY  on a diagonal law the fixed point is classic IRLS (same
%   influence, same MAD scale as the retired flat driver had); on a
%   correlated law the factors act on the LAW-WHITENED residuals (Yang 1994),
%   which the flat driver could not express.
%
%   EXAMPLE
%     R = lsqvce.engine.adjustRobust(systems, 1:3, 'Influence', "igg3");
%     R.report;  fac = R.history(end).factor;
%
%   See also LSQVCE.ENGINE.ADJUSTBLOCK, LSQVCE.KERNEL.ROBUSTWEIGHTS,
%            LSQVCE.ENGINE.RESIDUALS, LSQVCE.FIT
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-133: each pass writes its
%   factors into a new map, so the returned accumulator keeps the weights its
%   solve used; rejected rows (factor 0) leave the degrees of freedom
%   (makeLoader marks them excluded). With assistance from Claude (Anthropic),
%   model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-132: history(end).factor holds
%   the weights the RETURNED solve used (read back from the arc map), not the
%   next pass's -- R.x is the weighted solution on them at every stop. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-125: 'PsdMode' -- how a declared
%   psd law whitens (fit passes Options.Noise.Whitener). With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 4: the verdict comes
%   from lsqvce.engine.verdict.fromAccumulator (diagnoseAccum + diagnose
%   absorbed). With assistance from Claude (Anthropic), model Claude Opus
%   5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34, R7-35):
%   the ACCURACY note no longer cites the deleted flat driver; Prior /
%   Diagnostics / OnPass / OnPassContext documented. With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M2, WP-81 round 2b.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-115 (finding R7-17): the stop
%   rule's absolute floor of 1 on max|x| made it unit-dependent -- the same
%   fit in SH units (1e-9) stopped at pass 2 about one formal sigma from the
%   fixed point while reporting "converged". It is now relative to the
%   larger of max|x| and the largest formal sigma. With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-21 (CEST, UTC+02:00), R6 WP-110: the collapsed error said
%   it detected "more than half the data discrepant, beyond the breakdown
%   point". It does not, and the measurement is in HANDOVER 12.34: nAlive
%   tracks the CLEAN count, so the guard is reached by a bodily residual
%   offset or by near-total contamination, not at 50 %. The message and this
%   help now say what it detects. No behaviour changed -- the threshold is
%   the user's to set, and 12.34 puts that decision. With assistance from
%   Claude (Anthropic), model Claude Opus 5.
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-105: collectWhitened became
%   the engine-private lawWhitenedResiduals, the ONE definition this driver
%   and the joint estimator's robust pass share. With assistance from
%   Claude (Anthropic), model Claude Fable 5.1.

arguments
    systems
    days (1,:) double
    opt.Backend  (1,1) string {mustBeMember(opt.Backend, ["normal","qr"])} = "qr"
    opt.PsdMode  (1,1) string {mustBeMember(opt.PsdMode, ["","toeplitz","ar","arma"])} = ""
    opt.Influence (1,1) string {mustBeMember(opt.Influence, ["huber","igg3","bisquare"])} = "huber"
    opt.Tuning double = []
    opt.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 30
    opt.Tol (1,1) double {mustBePositive} = 1e-8
    opt.Parallel (1,1) string {mustBeMember(opt.Parallel, ["serial","parfor","spmd"])} = "serial"
    opt.DistributeDays (1,1) logical = false
    opt.ParallelThreshold (1,1) double {mustBeInteger, mustBeNonnegative} = 8
    opt.Verbose (1,1) logical = false
    opt.Prior = []                       % R5 M5 5b: CommonPrior on the common block, or []
    opt.Diagnostics cell = {}          % R5 M5 5b': Name-Value pairs for lsqvce.engine.verdict.diagnose
    opt.OnPass = []                    % R5 M5 5c': per-pass inspection hook (lsqvce.util.onPass)
    opt.OnPassContext struct = struct()
    opt.Provenance struct = struct()
end
    [specs, ~] = normalizeSystems(systems, days);       % engine-private, shared with adjustBlock
    nc = lsqvce.engine.resolveNc(systems);

    W = containers.Map('KeyType', 'char', 'ValueType', 'any');   % factors per arc
    x = [];  converged = false;  userStop = false;
    hist = struct('pass', {}, 'scale', {}, 'dx', {}, 'nDown', {}, 'nZero', {}, ...
                  'factor', {}, 'keys', {});
    if opt.Verbose
        lsqvce.util.log.msg('IRLS (%s) on the block path: %d system(s), %d day(s)\n', ...
                            opt.Influence, numel(specs), numel(days));
        lsqvce.util.log.msg('  %4s %14s %12s %10s %10s\n', 'pass', 'scale (MAD)', 'max|dx|', 'n down', 'n zero');
    end

    for it = 1:opt.MaxIter
        acc = newAcc(nc, specs, opt, W);
        acc.setPrior(opt.Prior);
        sol = acc.solve();
        res = lsqvce.engine.residuals(acc, sol, 'Whitened', true, 'Collect', true, 'Verbose', false);

        % law-whitened residuals of every observation, and the keys to put
        % the next factors back where they belong
        [vw, keys, sizes] = lawWhitenedResiduals(res, W);   % engine-private, THIS driver's alone (WP-108: adjustVCEjoint has its own producer and its own prior exclusion)
        [fac, s] = lsqvce.kernel.robustWeights(vw, opt.Influence, opt.Tuning);
        facUsed = appliedFactors(W, keys, sizes);   % WP-132: what THIS solve was weighted with

        nAlive = sum(fac > 0);
        if nAlive < nc + 1
            error('lsqvce:adjustRobust:collapsed', ...
                 ['The influence function left %d observations with a non-zero ' ...
                  'weight, fewer than the %d common parameters + 1: almost every ' ...
                  'row was rejected. Two things do that, and neither is an ' ...
                  'ordinary outlier problem. (1) The residuals are offset BODILY ' ...
                  'from zero by several MAD -- a bias the model did not absorb, ' ...
                  'a datum defect -- because the influence standardises vw/s, ' ...
                  'not (vw - med)/s. (2) Contamination so near total that fewer ' ...
                  'than %d clean rows are left; contamination merely past 50 %% ' ...
                  'does NOT reach here, since nAlive tracks the clean count. ' ...
                  'Check the model, not the data.'], ...
                  nAlive, nc, nc + 1);
        end
        warnAtBreakdown(fac, 'adjustRobust');      % engine-private: > 45 % zeroed, said once (R7 WP-115)

        dx = Inf;
        if ~isempty(x), dx = max(abs(sol.x(:) - x)); end
        x = sol.x(:);
        nDown = sum(fac < 0.5 & fac > 0);  nZero = sum(fac == 0);
        hist(end+1) = struct('pass', it, 'scale', s, 'dx', dx, 'nDown', nDown, ...
                             'nZero', nZero, 'factor', [], 'keys', {{}}); %#ok<AGROW>
        if opt.Verbose
            lsqvce.util.log.msg('  %4d %14.6e %12.3e %10d %10d\n', it, s, dx, nDown, nZero);
        end

        ctx = struct('driver', 'adjustRobust', 'pass', it, 'x', x, 'Cxx', sol.Cxx, ...
                     'sigma0', sol.sigma0, 's0sq', sol.s0sq, 'factor', fac, 'keys', {keys}, 'scale', s, ...
                     'dx', dx, 'acc', acc);
        userStop = ~lsqvce.util.onPass(opt.OnPass, opt.OnPassContext, ctx);
        % R7-17: relative to the problem's OWN scale -- the larger of the
        % parameters and their formal errors -- so the same fit in other
        % units stops on the same pass (the floor of 1 stopped SH-scale
        % coefficients, ~1e-9, at pass 2 ~1 sigma from the fixed point)
        if it > 1 && dx < opt.Tol * max(max(abs(x)), sqrt(max(diag(sol.Cxx))))
            converged = true;
            break
        end
        if userStop, break, end
        % write the factors for the next pass -- into a NEW map (WP-133, review
        % 3 finding 2): W is a handle that this pass's accumulator holds as
        % UserOpt.rowWeights and its loader reads on every call, so writing into
        % it re-weighted the RETURNED acc after the fact; at a maxPasses exit
        % R.diag, R.residuals, the redundancy and snoop replayed the next
        % pass's weights (sigma0 4.80 against 2.87 after one pass)
        W = copyMap(W);
        off = 0;
        for k = 1:size(keys, 1)
            n = sizes(k);
            W(weightKey(keys(k, :))) = fac(off + (1:n));
            off = off + n;
        end
    end
    % WP-132 (review low 24): the factors the RETURNED solve was weighted
    % with, read from W before the pass overwrote it. fac is the next pass's
    % (derived from this solve's residuals), one step ahead of sol / D at a
    % maxPasses or userStop exit, and a Lipschitz step ahead even at convergence.
    hist(end).factor = facUsed;
    hist(end).keys   = keys;

    if ~converged && ~userStop
        lsqvce.util.log.warn('lsqvce:adjustRobust:notConverged', ...
            ['IRLS reached MaxIter = %d (max|dx| = %.3e, scale %.3e); the last ' ...
             'pass is returned with stopReason "maxPasses".'], opt.MaxIter, dx, s);
    end

    % final diagnostics on the reweighted problem (the factors are in Wfcn)
    [D, dsol, batchFcn, nB] = lsqvce.engine.verdict.fromAccumulator(acc, sol, opt.Diagnostics{:});
    prov = opt.Provenance;
    if isempty(fieldnames(prov)), prov = rmfield(rmfield(rmfield(rmfield(rmfield(opt, 'Provenance'), 'Prior'), 'Diagnostics'), 'OnPass'), 'OnPassContext'); end
    R = lsqvce.Results(sol, D, acc, dsol, batchFcn, nB, specs, ...
                       lsqvce.engine.runMeta(struct(), prov));
    reason = "";  if userStop, reason = "userStop"; end
    R.setIterative(hist, struct('influence', char(opt.Influence), 'scale', s), converged, reason);
end

%==========================================================================
function f = appliedFactors(W, keys, sizes)
%APPLIEDFACTORS  The factors the last accumulation applied, in sweep order
%   (WP-132): W's entry per arc, 1 for an arc W does not hold (pass 1 is
%   unweighted). INPUT W containers.Map, keys n x 4 cell, sizes n x 1.
%   OUTPUT f sum(sizes) x 1 double.
    f = ones(sum(sizes), 1);  off = 0;
    for k = 1:size(keys, 1)
        n = sizes(k);  key = weightKey(keys(k, :));
        if W.isKey(key), f(off + (1:n)) = W(key); end
        off = off + n;
    end
end

%==========================================================================
function M = copyMap(W)
%COPYMAP  A new containers.Map with W's entries (WP-133): the accumulator a
%   pass built keeps its own, never written again.
    M = containers.Map('KeyType', 'char', 'ValueType', 'any');
    ks = keys(W);
    for i = 1:numel(ks), M(ks{i}) = W(ks{i}); end
end

%==========================================================================
function key = weightKey(row)
%WEIGHTKEY  W's key of one (sys, obsType, day, arc) row of keys.
    key = lsqvce.util.rowWeightKey(row{1}, row{3}, row{2}, row{4});
end

%==========================================================================
function acc = newAcc(nc, specs, opt, W)
%NEWACC  A fresh block accumulator whose loaders see the current factors.
    acc = lsqvce.engine.BlockAccumulator(nc, 'Backend', char(opt.Backend), ...
              'KeepLocal', true, 'Verbose', false, ...
              'Parallel', char(opt.Parallel), 'DistributeDays', opt.DistributeDays, ...
              'ParallelThreshold', opt.ParallelThreshold, ...
              'UserOpt', withPsd(struct('rowWeights', W), opt.PsdMode));
    acc.addSystems(specs);
end

%==========================================================================
function uo = withPsd(uo, md)
%WITHPSD  Add the psd whitening route to the loader UserOpt (WP-125); an
%   empty md leaves the covariance model's own default.
    if strlength(md) > 0, uo.psdMode = char(md); end
end
