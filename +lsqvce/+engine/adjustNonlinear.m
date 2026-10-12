function R = adjustNonlinear(u, batchFcn, nBatch, x0, opt)
%LSQVCE.ENGINE.ADJUSTNONLINEAR  Gauss-Newton with Levenberg damping on the block path.
%   R = LSQVCE.ENGINE.ADJUSTNONLINEAR(u, batchFcn, nBatch, x0, 'Name', value, ...)
%
%   Re-linearisation as a POLICY over the one accumulator (R5 M2, WP-81
%   round 2c'; finding H-02): the batched source [A, w, p] = batchFcn(k, x)
%   (Jacobian and reduced observation l - f(x) at the current x) is read as
%   a local-free streaming system (LSQVCE.ENGINE.BATCHSYSTEM with a State
%   map holding x); every pass accumulates at the current x, solves for the
%   step, and the NEXT pass's accumulation at x + dx is also the step's
%   evaluation (its F decides acceptance), so nothing is read twice.
%
%   THE STEP (R7 WP-113). One least-squares solve on the accumulation at x
%   with two kinds of pseudo-observations on the COMMON block:
%     the prior, shifted to x:   K dx ~ x0 - K x            (R7-04)
%     Marquardt damping:         sqrt(lam d_j) dx_j ~ 0,  d = diag(J'PJ)  (R7-14)
%   Until 6.9 the prior entered as K dx ~ x0 -- a prior on the INCREMENT,
%   whose fixed point is not the declared optimum (b reported +0.60 against
%   the regularised +0.30 on the decay fixture) -- and the damping was the
%   accumulator's Reg "lm", which acts on a LOCAL block; a batch source has
%   none, so the damping was a no-op and a rejected step was re-solved
%   unchanged until MaxIter (review R7, dom_lsq F2 / F3).
%
%   STEP CONTROL. F(x) = l'Pl(x) + |(K x - x0)/sigma|^2, the DECLARED objective.
%   A step that raises F is rejected -- from the first pass on -- and lam
%   grows tenfold (0 -> 1e-3; above 1e12 -> diverged), the shorter step
%   re-solved on the SAME accumulation; an accepted step shrinks lam tenfold
%   (below 1e-12 -> 0, pure Gauss-Newton). Convergence when max|dx_j/sigma_j|
%   < TolX, or the relative change of F < TolF. A converged run returns x plus
%   the final UNDAMPED step; an unconverged one returns the best ACCEPTED point
%   x (never a rejected one). The covariance is that of the undamped
%   linearised problem at x. Mirror (mirrors/wp113/gauss_newton.py): the
%   exponential fit from [10 -5] / [0.1 1] and Rosenbrock from five starts
%   reach the optimum (|x - x*| <= 3e-9, 5e-16); the prior fixture lands on
%   the declared optimum to 1e-7.
%
%   INPUT
%     u         1x1 double          unknowns
%     batchFcn  1x1 function_handle [A, w(, p)] = batchFcn(k, x)
%     nBatch    1x1 double          batches
%     x0        u x 1 double        starting point
%   NAME-VALUE
%     'MaxIter' 20; 'TolX' 1e-8; 'TolF' 1e-10; 'Lambda0' 1e-6 (0 = start
%     undamped; a rejected step still turns the damping on);
%     'Backend' "qr"|"normal"; 'Parallel', 'DistributeDays', 'ParallelThreshold';
%     'Verbose' false; 'Provenance' struct (fit's resolved Options);
%     'Prior' [] | lsqvce.engine.CommonPrior on the common block; 'Diagnostics'
%     {} (Name-Value pairs for lsqvce.engine.verdict.diagnose); 'OnPass' [] | @(ctx)
%     the per-pass hook (lsqvce.util.onPass) with 'OnPassContext' struct()
%   OUTPUT
%     R  1x1 lsqvce.Results  stopReason "converged" | "maxPasses"; history(k):
%        .pass .vPv (F(x), the prior included) .dxRel .lambda .accepted;
%        R.noise = struct('x0', x0, 'reason')
%
%   ERRORS
%     lsqvce:adjustNonlinear:badX0       x0 is not u x 1
%     lsqvce:adjustNonlinear:nonFinite   the model returned non-finite residuals
%   WARNINGS
%     lsqvce:adjustNonlinear:diverged    damping reached 1e12 without a downhill step
%     lsqvce:adjustNonlinear:maxIter     MaxIter reached (stopReason "maxPasses")
%
%   ACCURACY  the Gauss-Newton fixed point of the dense reference in test_fit (x to 1e-8 on
%   the exponential-decay fixture); with a prior, the optimum of the declared
%   objective (test_fit, R7 WP-113).
%
%   EXAMPLE
%     bf = @(k, x) deal(J(x, k), l(k) - f(x, k));       % your model
%     R  = lsqvce.engine.adjustNonlinear(2, bf, 4, [0.5; 0]);
%
%   See also LSQVCE.ENGINE.BATCHSYSTEM,
%            LSQVCE.ENGINE.ADJUSTBLOCK, LSQVCE.FIT
%
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-127: a converged run's final
%   undamped step is kept only when F does not rise (on a plateau TolX met a
%   tiny damped step and the undamped one was never evaluated). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-113 (R7-14, R7-04): the step
%   damps the common block (Marquardt rows) and shifts the prior to x; the
%   descent check runs from pass 1 on the declared objective; a rejected step
%   is re-solved on the same accumulation; an unconverged run returns the best
%   accepted point. With assistance from Claude (Anthropic), model Claude Opus
%   5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 4: the verdict comes
%   from lsqvce.engine.verdict.fromAccumulator (diagnoseAccum + diagnose
%   absorbed). With assistance from Claude (Anthropic), model Claude Opus
%   5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34, R7-35):
%   the deleted SOLVENONLINEAR reference removed; Prior / Diagnostics /
%   OnPass / OnPassContext documented. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M2, WP-81 round 2c'.

arguments
    u (1,1) double {mustBeInteger, mustBePositive}
    batchFcn
    nBatch (1,1) double {mustBeInteger, mustBePositive}
    x0 double
    opt.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 20
    opt.TolX (1,1) double {mustBePositive} = 1e-8
    opt.TolF (1,1) double {mustBePositive} = 1e-10
    opt.Lambda0 (1,1) double {mustBeNonnegative} = 1e-6
    opt.Backend (1,1) string {mustBeMember(opt.Backend, ["normal","qr"])} = "qr"
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
    x0 = x0(:);
    if numel(x0) ~= u
        error('lsqvce:adjustNonlinear:badX0', 'x0 must have %d entries; it has %d.', u, numel(x0));
    end
    state = containers.Map('KeyType', 'char', 'ValueType', 'any');
    state('x') = x0;
    L = lsqvce.engine.batchSystem(u, batchFcn, nBatch, 'State', state, 'Name', 'nonlinear');
    days = 1:nBatch;

    x = x0;  lam = opt.Lambda0;  converged = false;  reason = 'maxIter';  userStop = false;
    hist = struct('pass', {}, 'vPv', {}, 'dxRel', {}, 'lambda', {}, 'accepted', {});
    if opt.Verbose
        lsqvce.util.log.msg('Gauss-Newton on the block path: u = %d, %d batch(es)\n', u, nBatch);
        lsqvce.util.log.msg('  %4s %16s %14s %12s\n', 'pass', 'F(x)', 'max|dx/sig|', 'lambda');
    end

    % R7 WP-113: one accumulation per point x (the evaluation of F and the
    % linearisation at once); the STEP is a solve on that accumulation with the
    % prior shifted to the current x and the damping on the common block, so
    % a rejected step is re-solved at the same x without re-reading the data.
    [acc, F, d] = accumulate(L, days, x, state, opt);
    sol = stepSolve(acc, x, lam, d, opt);
    for it = 1:opt.MaxIter
        dx   = sol.x(:);
        xNew = x + dx;
        [accNew, FNew, dNew] = accumulate(L, days, xNew, state, opt);  % evaluation AND next linearisation
        if ~isfinite(FNew) && ~isnan(FNew)
            error('lsqvce:adjustNonlinear:nonFinite', ...
                  'The model returned non-finite residuals at pass %d.', it);
        end
        if FNew > F                                           % R7-14: checked from pass 1
            if lam == 0, lam = 1e-3; else, lam = lam * 10; end   % rejected: do not move
            hist(end+1) = struct('pass', it, 'vPv', FNew, 'dxRel', NaN, 'lambda', lam, 'accepted', false); %#ok<AGROW>
            if opt.Verbose, lsqvce.util.log.msg('  %4d %16.8e %14s %12.1e   step rejected\n', it, FNew, '-', lam); end
            if lam > 1e12
                reason = 'diverged';
                lsqvce.util.log.warn('lsqvce:adjustNonlinear:diverged', ...
                    ['Levenberg damping reached %.1e without finding a downhill step. ' ...
                     'Either the starting point is far from the solution, the Jacobian ' ...
                     'is wrong (check it against a finite difference), or the model is ' ...
                     'not locally identifiable at this x.'], lam);
                break
            end
            sol = stepSolve(acc, x, lam, d, opt);             % same x, more damping: a SHORTER step
            continue
        end
        sigX  = sqrt(max(real(diag(sol.Cxx)), 0)) / max(sol.sigma0, realmin);
        if all(isfinite(sigX)) && all(sigX > 0)
            dxRel = max(abs(dx) ./ sigX);                 % the step in formal-error units
        else
            dxRel = max(abs(dx)) / max(max(abs(x)), 1);   % sigma0 NaN-honest: a relative step
        end
        Fold = F;
        x = xNew;  acc = accNew;  F = FNew;  d = dNew;
        lam = lam / 10;  if lam < 1e-12, lam = 0; end
        sol = stepSolve(acc, x, lam, d, opt);                 % the next step
        hist(end+1) = struct('pass', it, 'vPv', F, 'dxRel', dxRel, 'lambda', lam, 'accepted', true); %#ok<AGROW>
        if opt.Verbose, lsqvce.util.log.msg('  %4d %16.8e %14.3e %12.1e\n', it, F, dxRel, lam); end
        relF = abs(F - Fold) / max(abs(F), realmin);
        ctx = struct('driver', 'adjustNonlinear', 'pass', it, 'x', x, 'Cxx', sol.Cxx, ...
                     'sigma0', sol.sigma0, 's0sq', sol.s0sq, 'vPv', F, 'dxRel', dxRel, ...
                     'lambda', lam, 'acc', acc);
        userStop = ~lsqvce.util.onPass(opt.OnPass, opt.OnPassContext, ctx);
        if dxRel < opt.TolX
            converged = true;  reason = 'tolX';  break
        end
        if it > 1 && relF < opt.TolF
            converged = true;  reason = 'tolF';  break
        end
        if userStop, reason = 'userStop';  break, end
    end
    if ~converged && ~strcmp(reason, 'diverged') && ~userStop
        lastRel = NaN;  ok = [hist.accepted];
        if any(ok), lastRel = hist(find(ok, 1, 'last')).dxRel; end
        lsqvce.util.log.warn('lsqvce:adjustNonlinear:maxIter', ...
            ['Reached MaxIter = %d without meeting the tolerance (last max|dx/sigma| = %.3e, ' ...
             'target %.3e). The diagnostics describe an UNCONVERGED solution.'], ...
             opt.MaxIter, lastRel, opt.TolX);
    end

    prov = opt.Provenance;
    if isempty(fieldnames(prov)), prov = rmfield(rmfield(rmfield(rmfield(rmfield(opt, 'Provenance'), 'Prior'), 'Diagnostics'), 'OnPass'), 'OnPassContext'); end
    acc = finalAcc(L, days, x, state, opt);
    dsol = acc.solve();
    % R7-14: a converged run takes the final undamped step (below TolX, or
    % with a settled F) when it does not raise F (WP-127); an UNconverged one
    % returns the best ACCEPTED point x -- never x plus a step the loop had
    % just rejected. The covariance and the
    % diagnostics are those of the linearised problem at x either way.
    sol  = dsol;  sol.x = x;
    if converged
        % WP-127: the final undamped step is CHECKED like every other. On a
        % plateau sigX is huge, so a tiny damped step meets TolX, and the
        % undamped step from there took F from 28.8 to 49.6 in the mirror of
        % a logistic fit -- "converged" on a point never evaluated
        xF = x + dsol.x(:);
        [~, FF] = accumulate(L, days, xF, state, opt);
        state('x') = x;                       % the diagnostics replay the batches at x
        % tol: a rise below 1e-10 relative is the rounding of F's sum over the
        % rows (~n eps up to n ~ 1e6), not ascent -- at the optimum the step is ~0
        if isfinite(FF) && FF <= F * (1 + 1e-10)
            sol.x = xF;
        else
            reason = [reason '; final step rejected'];
        end
    end
    [D, dsol2, batchFcn, nB] = lsqvce.engine.verdict.fromAccumulator(acc, dsol, opt.Diagnostics{:});
    R = lsqvce.Results(sol, D, acc, dsol2, batchFcn, nB, [], lsqvce.engine.runMeta(struct(), prov));
    sr = "";  if userStop, sr = "userStop"; end
    R.setIterative(hist, struct('x0', x0, 'reason', reason, 'kind', 'nonlinear'), converged, sr);
end

function acc = finalAcc(L, days, x, state, opt)
%FINALACC  The undamped accumulation at x (KeepLocal, for the diagnostics),
%   the prior shifted to x (R7-04).
    state('x') = x;
    acc = lsqvce.engine.BlockAccumulator(L.nc, 'Backend', char(opt.Backend), ...
              'KeepLocal', true, 'Verbose', false, 'Parallel', char(opt.Parallel), ...
              'DistributeDays', opt.DistributeDays, 'ParallelThreshold', opt.ParallelThreshold);
    acc.addSystem(L.sysID, L.nl, L.closure, days);
    acc.setPrior(stepPrior(opt.Prior, x, 0, []));
end

function [acc, F, d] = accumulate(L, days, x, state, opt)
%ACCUMULATE  One accumulation at x: the linearisation for the next step and
%   the evaluation of the DECLARED objective F(x) = l'Pl(x) + |(K x - x0)/s|^2
%   (the data's whitened residual mass plus the prior's -- R7-04: the prior
%   is part of what is minimised, so it is part of the descent check), and
%   d = diag(J'PJ), Marquardt's scaling of the damping.
    state('x') = x;
    acc = lsqvce.engine.BlockAccumulator(L.nc, 'Backend', char(opt.Backend), 'KeepLocal', false, ...
              'Verbose', false, 'Parallel', char(opt.Parallel), 'DistributeDays', opt.DistributeDays, ...
              'ParallelThreshold', opt.ParallelThreshold);
    acc.addSystem(L.sysID, L.nl, L.closure, days);
    [N, ~, F] = acc.sys.commonBlocks();          % l'Pl at x (data)
    if ~isempty(opt.Prior), F = F + opt.Prior.vtv(x); end
    d = diag(N);
end

function sol = stepSolve(acc, x, lam, d, opt)
%STEPSOLVE  The Gauss-Newton step at x on the accumulation acc: the prior as
%   K dx ~ x0 - K x (R7-04 -- on the increment the prior's TARGET must move
%   with x) and, for lam > 0, Marquardt damping of the COMMON block as the
%   pseudo-observations sqrt(lam d_j) dx_j ~ 0 (R7-14 -- the accumulator's
%   'Reg' "lm" acts on a LOCAL block, and a batch source has none).
    acc.setPrior(stepPrior(opt.Prior, x, lam, d));
    sol = acc.solve();
end

function P = stepPrior(pri, x, lam, d)
%STEPPRIOR  The CommonPrior of one step: the declared prior shifted to x,
%   plus the damping rows at the same sigma (so their weight is lam d). [] when
%   there is neither.
    nc = numel(x);
    K = zeros(0, nc);  t = zeros(0, 1);  s = 1;  g = zeros(0, 1);
    if ~isempty(pri)
        K = pri.K;  t = pri.x0 - pri.K * x;  s = pri.sigma;  g = pri.group;
    end
    if lam > 0
        dd = max(d(:), eps * max([d(:); realmin]));   % a parameter J does not see is damped too
        K = [K; s * sqrt(lam) * diag(sqrt(dd))];
        t = [t; zeros(nc, 1)];
        g = [g; (max([g; 0]) + 1) * ones(nc, 1)];     % the damping is a group of its own (WP-119)
    end
    P = [];
    if ~isempty(K), P = lsqvce.engine.CommonPrior(K, t, s, g); end
end
