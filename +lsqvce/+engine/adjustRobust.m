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
%   with R.history(end).keys naming the (sys, obsType, day, arc) rows.
%
%   INPUT
%     systems  1xM cell | struct  as for LSQVCE.ENGINE.ADJUSTBLOCK
%     days     1xD double          the calendar
%   NAME-VALUE
%     'Backend'    "qr" | "normal"                          default "qr"
%     'Influence'  "huber" | "igg3" | "bisquare"            default "huber"
%     'Tuning'     [] | constants of the influence          default [] (its default)
%     'MaxIter'    passes (>= 1)                            default 30
%     'Tol'        stop when max|dx| < Tol * max(max|x|, 1) default 1e-8
%     'Parallel', 'DistributeDays', 'ParallelThreshold'     as adjustBlock
%     'Verbose'    logical                                  default false
%     'Provenance' struct   resolved user options for Results.meta (fit)
%   OUTPUT
%     R  1x1 lsqvce.Results   stopReason "converged" | "maxPasses"; history(k):
%        .pass .scale .dx .nDown (0 < factor < 0.5) .nZero (factor == 0)
%        .factor (last pass: n x 1) .keys (last pass: n x 4 cell)
%
%   ERRORS
%     lsqvce:adjustRobust:collapsed   fewer live observations than parameters
%                                     + 1 after reweighting (beyond the
%                                     breakdown point: check the model)
%     + everything adjustBlock / makeLoader / robustWeights raise
%   WARNINGS
%     lsqvce:adjustRobust:notConverged  MaxIter reached; the last pass is
%                                       returned with stopReason "maxPasses"
%
%   ACCURACY  the fixed point equals lsqvce.engine.robust's on a diagonal
%   law (same influence, same MAD scale; pinned in test_fit to 1e-10); on a
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
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M2, WP-81 round 2b.
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-105: collectWhitened became
%   the engine-private lawWhitenedResiduals, the ONE definition this driver
%   and the joint estimator's robust pass share. With assistance from
%   Claude (Anthropic), model Claude Fable 5.1.

arguments
    systems
    days (1,:) double
    opt.Backend  (1,1) string {mustBeMember(opt.Backend, ["normal","qr"])} = "qr"
    opt.Influence (1,1) string {mustBeMember(opt.Influence, ["huber","igg3","bisquare"])} = "huber"
    opt.Tuning double = []
    opt.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 30
    opt.Tol (1,1) double {mustBePositive} = 1e-8
    opt.Parallel (1,1) string {mustBeMember(opt.Parallel, ["serial","parfor","spmd"])} = "serial"
    opt.DistributeDays (1,1) logical = false
    opt.ParallelThreshold (1,1) double {mustBeInteger, mustBeNonnegative} = 8
    opt.Verbose (1,1) logical = false
    opt.Prior = []                       % R5 M5 5b: CommonPrior on the common block, or []
    opt.Diagnostics cell = {}          % R5 M5 5b': Name-Value pairs for lsqvce.engine.diagnose
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

        nAlive = sum(fac > 0);
        if nAlive < nc + 1
            error('lsqvce:adjustRobust:collapsed', ...
                 ['The influence function left %d observations with a non-zero ' ...
                  'weight, fewer than the %d parameters + 1: more than half the ' ...
                  'data is discrepant, which is beyond the breakdown point of ' ...
                  'any M-estimator. This is not an outlier problem -- check the model.'], ...
                  nAlive, nc);
        end

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
        if it > 1 && dx < opt.Tol * max(max(abs(x)), 1)
            converged = true;
            break
        end
        if userStop, break, end
        % write the factors for the next pass
        off = 0;
        for k = 1:size(keys, 1)
            n = sizes(k);
            W(lsqvce.util.rowWeightKey(keys{k,1}, keys{k,3}, keys{k,2}, keys{k,4})) = fac(off + (1:n));
            off = off + n;
        end
    end
    hist(end).factor = fac;
    hist(end).keys   = keys;

    if ~converged && ~userStop
        lsqvce.util.log.warn('lsqvce:adjustRobust:notConverged', ...
            ['IRLS reached MaxIter = %d (max|dx| = %.3e, scale %.3e); the last ' ...
             'pass is returned with stopReason "maxPasses".'], opt.MaxIter, dx, s);
    end

    % final diagnostics on the reweighted problem (the factors are in Wfcn)
    [dsol, batchFcn, nB] = lsqvce.engine.diagnoseAccum(acc, sol);
    D = lsqvce.engine.diagnose(dsol, batchFcn, 'nBatch', nB, opt.Diagnostics{:});
    prov = opt.Provenance;
    if isempty(fieldnames(prov)), prov = rmfield(rmfield(rmfield(rmfield(rmfield(opt, 'Provenance'), 'Prior'), 'Diagnostics'), 'OnPass'), 'OnPassContext'); end
    R = lsqvce.Results(sol, D, acc, dsol, batchFcn, nB, specs, ...
                       lsqvce.engine.runMeta(struct(), prov));
    reason = "";  if userStop, reason = "userStop"; end
    R.setIterative(hist, struct('influence', char(opt.Influence), 'scale', s), converged, reason);
end

%==========================================================================
function acc = newAcc(nc, specs, opt, W)
%NEWACC  A fresh block accumulator whose loaders see the current factors.
    acc = lsqvce.engine.BlockAccumulator(nc, 'Backend', char(opt.Backend), ...
              'KeepLocal', true, 'Verbose', false, ...
              'Parallel', char(opt.Parallel), 'DistributeDays', opt.DistributeDays, ...
              'ParallelThreshold', opt.ParallelThreshold, ...
              'UserOpt', struct('rowWeights', W));
    acc.addSystems(specs);
end
