function R = adjustNonlinear(u, batchFcn, nBatch, x0, opt)
%LSQVCE.ENGINE.ADJUSTNONLINEAR  Gauss-Newton with Levenberg damping on the block path.
%   R = LSQVCE.ENGINE.ADJUSTNONLINEAR(u, batchFcn, nBatch, x0, 'Name', value, ...)
%
%   Re-linearisation as a POLICY over the one accumulator (R5 M2, WP-81
%   round 2c'; finding H-02): the batched source [A, w, p] = batchFcn(k, x)
%   (Jacobian and reduced observation l - f(x) at the current x) is read as
%   a local-free streaming system (LSQVCE.ENGINE.BATCHSYSTEM with a State
%   map holding x); every pass accumulates at the current x with the
%   accumulator's own Levenberg damping (Reg "lm"), solves for the step,
%   and the NEXT pass's accumulation at x + dx is also the step's
%   evaluation (its vPv decides acceptance), so nothing is read twice.
%
%   Step control mirrors LSQVCE.ENGINE.SOLVENONLINEAR exactly: a step that
%   raises vPv is rejected and the damping grows tenfold (up to 1e12 ->
%   diverged); an accepted step shrinks it tenfold; convergence when
%   max|dx_j / sigma_j| < TolX, or the relative vPv change < TolF. The final
%   Results is an UNDAMPED solve at the converged x, so its covariance is
%   that of the linearised problem, not of the damped one.
%
%   INPUT
%     u         1x1 double          unknowns
%     batchFcn  1x1 function_handle [A, w(, p)] = batchFcn(k, x)
%     nBatch    1x1 double          batches
%     x0        u x 1 double        starting point
%   NAME-VALUE
%     'MaxIter' 20; 'TolX' 1e-8; 'TolF' 1e-10; 'Lambda0' 1e-6 (0 = pure GN);
%     'Backend' "qr"|"normal"; 'Parallel', 'DistributeDays', 'ParallelThreshold';
%     'Verbose' false; 'Provenance' struct (fit's resolved Options)
%   OUTPUT
%     R  1x1 lsqvce.Results  stopReason "converged" | "maxPasses"; history(k):
%        .pass .vPv .dxRel .lambda .accepted; R.noise = struct('x0', x0, 'reason')
%
%   ERRORS
%     lsqvce:adjustNonlinear:badX0       x0 is not u x 1
%     lsqvce:adjustNonlinear:nonFinite   the model returned non-finite residuals
%   WARNINGS
%     lsqvce:adjustNonlinear:diverged    damping reached 1e12 without a downhill step
%     lsqvce:adjustNonlinear:maxIter     MaxIter reached (stopReason "maxPasses")
%
%   ACCURACY  the Gauss-Newton fixed point of the dense reference in test_fit (x to 1e-8 on
%   the exponential-decay fixture); the block path's vPv comes from the
%   accumulator's normal-form fold (NaN-honest on cancellation).
%
%   EXAMPLE
%     bf = @(k, x) deal(J(x, k), l(k) - f(x, k));       % your model
%     R  = lsqvce.engine.adjustNonlinear(2, bf, 4, [0.5; 0]);
%
%   See also LSQVCE.ENGINE.SOLVENONLINEAR, LSQVCE.ENGINE.BATCHSYSTEM,
%            LSQVCE.ENGINE.ADJUSTBLOCK, LSQVCE.FIT
%
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
    opt.Diagnostics cell = {}          % R5 M5 5b': Name-Value pairs for lsqvce.engine.diagnose
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

    x = x0;  lam = opt.Lambda0;  vPvPrev = Inf;  converged = false;  reason = 'maxIter';  userStop = false;
    hist = struct('pass', {}, 'vPv', {}, 'dxRel', {}, 'lambda', {}, 'accepted', {});
    if opt.Verbose
        lsqvce.util.log.msg('Gauss-Newton on the block path: u = %d, %d batch(es)\n', u, nBatch);
        lsqvce.util.log.msg('  %4s %16s %14s %12s\n', 'pass', 'vPv', 'max|dx/sig|', 'lambda');
    end

    [sol, acc] = accumulate(L, days, x, lam, state, opt);   % pass 1 at x0
    vPv = vpvOf(sol);
    for it = 1:opt.MaxIter
        dx   = sol.x(:);                                      % the (damped) step
        xNew = x + dx;
        [solNew, accNew] = accumulate(L, days, xNew, lam, state, opt);  % evaluation AND next linearisation
        vPvNew = vpvOf(solNew);
        if ~isfinite(vPvNew) && ~isnan(vPvNew)
            error('lsqvce:adjustNonlinear:nonFinite', ...
                  'The model returned non-finite residuals at pass %d.', it);
        end
        if opt.Lambda0 > 0 && it > 1 && vPvNew > vPv
            lam = lam * 10;                                   % rejected: do not move
            hist(end+1) = struct('pass', it, 'vPv', vPvNew, 'dxRel', NaN, 'lambda', lam, 'accepted', false); %#ok<AGROW>
            if opt.Verbose, lsqvce.util.log.msg('  %4d %16.8e %14s %12.1e   step rejected\n', it, vPvNew, '-', lam); end
            if lam > 1e12
                reason = 'diverged';
                lsqvce.util.log.warn('lsqvce:adjustNonlinear:diverged', ...
                    ['Levenberg damping reached %.1e without finding a downhill step. ' ...
                     'Either the starting point is far from the solution, the Jacobian ' ...
                     'is wrong (check it against a finite difference), or the model is ' ...
                     'not locally identifiable at this x.'], lam);
                break
            end
            [sol, acc] = accumulate(L, days, x, lam, state, opt);  % re-linearise at x, more damping
            continue
        end
        % accepted
        sigX  = sqrt(max(real(diag(sol.Cxx)), 0)) / max(sol.sigma0, realmin);
        if all(isfinite(sigX)) && all(sigX > 0)
            dxRel = max(abs(dx) ./ sigX);                 % the step in formal-error units
        else
            dxRel = max(abs(dx)) / max(max(abs(x)), 1);   % sigma0 NaN-honest: a relative step
        end
        x = xNew;  sol = solNew;  acc = accNew;
        lam = max(lam / 10, 1e-12);
        hist(end+1) = struct('pass', it, 'vPv', vPvNew, 'dxRel', dxRel, 'lambda', lam, 'accepted', true); %#ok<AGROW>
        if opt.Verbose, lsqvce.util.log.msg('  %4d %16.8e %14.3e %12.1e\n', it, vPvNew, dxRel, lam); end
        relF = abs(vPvNew - vPv) / max(abs(vPvNew), realmin);
        ctx = struct('driver', 'adjustNonlinear', 'pass', it, 'x', x, 'Cxx', sol.Cxx, ...
                     'sigma0', sol.sigma0, 's0sq', sol.s0sq, 'vPv', vPvNew, 'dxRel', dxRel, ...
                     'lambda', lam, 'acc', acc);
        userStop = ~lsqvce.util.onPass(opt.OnPass, opt.OnPassContext, ctx);
        if dxRel < opt.TolX
            converged = true;  reason = 'tolX';  break
        end
        if it > 1 && relF < opt.TolF
            converged = true;  reason = 'tolF';  break
        end
        if userStop, reason = 'userStop';  break, end
        vPv = vPvNew;  vPvPrev = vPv; %#ok<NASGU>
    end
    if ~converged && ~strcmp(reason, 'diverged') && ~userStop
        lsqvce.util.log.warn('lsqvce:adjustNonlinear:maxIter', ...
            ['Reached MaxIter = %d without meeting the tolerance (last max|dx/sigma| = %.3e, ' ...
             'target %.3e). The diagnostics describe an UNCONVERGED solution.'], ...
             opt.MaxIter, hist(end).dxRel, opt.TolX);
    end

    % final: an UNDAMPED linearised solve at x. Its x is the (converged) STEP
    % and its Cxx the covariance of the linearised problem; the estimate the
    % Results reports is x + step, while the diagnostics (dsol) stay in the
    % linearised frame where the reduced observations live.
    prov = opt.Provenance;
    if isempty(fieldnames(prov)), prov = rmfield(rmfield(rmfield(rmfield(rmfield(opt, 'Provenance'), 'Prior'), 'Diagnostics'), 'OnPass'), 'OnPassContext'); end
    acc = finalAcc(L, days, x, state, opt);
    dsol = acc.solve();
    sol  = dsol;  sol.x = x + dsol.x(:);
    [dsol2, batchFcn, nB] = lsqvce.engine.diagnoseAccum(acc, dsol);
    D = lsqvce.engine.diagnose(dsol2, batchFcn, 'nBatch', nB, opt.Diagnostics{:});
    R = lsqvce.Results(sol, D, acc, dsol2, batchFcn, nB, [], lsqvce.engine.runMeta(struct(), prov));
    sr = "";  if userStop, sr = "userStop"; end
    R.setIterative(hist, struct('x0', x0, 'reason', reason, 'kind', 'nonlinear'), converged, sr);
end

function acc = finalAcc(L, days, x, state, opt)
    state('x') = x;
    acc = lsqvce.engine.BlockAccumulator(L.nc, 'Backend', char(opt.Backend), ...
              'KeepLocal', true, 'Verbose', false, 'Parallel', char(opt.Parallel), ...
              'DistributeDays', opt.DistributeDays, 'ParallelThreshold', opt.ParallelThreshold);
    acc.setPrior(opt.Prior);
    acc.addSystem(L.sysID, L.nl, L.closure, days);
end

%==========================================================================
function [sol, acc] = accumulate(L, days, x, lam, state, opt)
    state('x') = x;
    nv = {'Backend', char(opt.Backend), 'KeepLocal', false, 'Verbose', false, ...
          'Parallel', char(opt.Parallel), 'DistributeDays', opt.DistributeDays, ...
          'ParallelThreshold', opt.ParallelThreshold};
    if lam > 0
        nv = [nv, {'Reg', 'lm', 'RegLambda', lam}];
    end
    acc = lsqvce.engine.BlockAccumulator(L.nc, nv{:});
    acc.setPrior(opt.Prior);
    acc.addSystem(L.sysID, L.nl, L.closure, days);
    sol = acc.solve();
    [~, ~, ltl] = acc.sys.commonBlocks();          % F(x): the whitened l'Pl at x
    sol.lPl = ltl;
end

function v = vpvOf(sol)
%VPVOF  The nonlinear OBJECTIVE F(x) = sum p_i w_i^2 at the linearisation point.
%   The reduced observations w = l - f(x) ARE the residuals at x, so F(x) is
%   the accumulated l'Pl of the whitened rows -- a plain sum of squares, never
%   a cancellation-prone difference. (The accumulator's vtv = l'l - n'x is
%   NaN-honest on purpose and went NaN at the converged trilateration point,
%   after which the step control compared NaNs for ever -- R5 M5 5b'.)
    v = NaN;
    if isfield(sol, 'lPl') && ~isempty(sol.lPl), v = sol.lPl; end
    if isnan(v) && isfield(sol, 'vtv') && ~isempty(sol.vtv), v = sol.vtv; end
end
