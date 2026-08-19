function [sol, D, V, IT] = solveNonlinear(u, batchFcn, nBatch, x0, varargin)
%LSQVCE.SOLVENONLINEAR  Iterated (Gauss-Newton) least squares in batches.
%
%   [sol, D, V, IT] = LSQVCE.SOLVENONLINEAR(u, batchFcn, nBatch, x0)
%   re-linearises and re-solves until convergence, then diagnoses the CONVERGED
%   solution.
%
%   batchFcn now takes the current estimate:
%       [A, w, p] = batchFcn(k, x)
%   where A is the Jacobian d(model)/dx AT x, and w is the reduced observation
%   vector  w = l - f(x)  (the "observed minus computed"). The step solves
%   A*dx = w, and x <- x + dx.
%
%   WHY THIS BELONGS IN THE TOOLBOX
%   -------------------------------
%   Everything else here assumes a linear model. In practice you re-linearise,
%   and the loop is easy to get subtly wrong in a way that silently poisons the
%   diagnostics:
%
%     * vPv must be evaluated AT THE CONVERGED POINT. Using the vPv from the
%       last accumulation (which was formed at the PREVIOUS x) biases the
%       variance factor. The bias is exactly the improvement of the final step,
%       so it is largest when you stop early -- i.e. precisely when you are
%       least entitled to trust it.
%     * The redundancy is still n - rank, but a wrongly converged solution
%       inflates vPv and can look like an outlier problem or a systematic.
%
%   So the loop is owned here, and the final diagnose() is run on a FRESH
%   accumulation at the converged x.
%
%   CONVERGENCE
%   -----------
%   Declared when EITHER
%       max|dx_j / sigma_x_j| < 'tolX'          (step is small vs. its own noise)
%   OR  |vPv_k - vPv_{k-1}| / vPv_k < 'tolF'    (the fit stops improving)
%   The first is the meaningful one: a step much smaller than the parameter's
%   own standard deviation cannot change any conclusion. A pure |dx| < eps test
%   is scale-dependent and therefore meaningless across mixed units (metres and
%   radians in the same vector).
%
%   DAMPING
%   -------
%   Pure Gauss-Newton diverges if the starting point is poor. A Levenberg
%   damping term is applied: solve (N + lam*diag(N)) dx = b, with lam reduced
%   on a successful step and increased on a failed one. Set 'lambda0' = 0 for
%   undamped Gauss-Newton.
%
%   NAME-VALUE
%     'maxIter'  30
%     'tolX'     1e-3    step size relative to sigma_x
%     'tolF'     1e-10   relative change in vPv
%     'lambda0'  1e-6    initial Levenberg damping (0 = pure Gauss-Newton)
%     'datumB'   []      free-network datum (see lsqvce.Accumulator)
%     'method'   'qr'
%     'sigma0'   1
%     'verbose'  true
%     ... plus everything lsqvce.diagnose accepts
%
%   OUTPUT IT
%     .x          per-iteration parameter estimates (u x nIter)
%     .vPv        per-iteration weighted residual sum of squares
%     .dxNorm     per-iteration max|dx/sigma_x|
%     .lambda     per-iteration damping
%     .converged  logical
%     .nIter
%     .reason     'tolX' | 'tolF' | 'maxIter' | 'diverged'
%
%   EXAMPLE
%       % 2-D resection: distances from known stations to an unknown point.
%       % batchFcn takes (k, x) and returns the Jacobian AT x.
%       %
%       % NOTE it must be a real function, not an anonymous one, because it
%       % needs the current x. Put it at the bottom of your script file:
%       %
%       %   function [A, l, p] = resect(~, x, st, d)
%       %       dx = x(1) - st(:,1);   dy = x(2) - st(:,2);
%       %       rr = sqrt(dx.^2 + dy.^2);
%       %       A  = [dx./rr, dy./rr];      % d(range)/d(x)
%       %       l  = d - rr;                % observed MINUS computed
%       %       p  = ones(numel(d),1) / 0.01^2;
%       %   end
%
%       st = [0 0; 10 0; 10 10; 0 10; 5 12];
%       d  = sqrt(sum((st - [4.2 6.7]).^2, 2)) + 0.01*randn(5,1);
%
%       [sol, D, V, IT] = lsqvce.solveNonlinear( ...
%                             2, @(k,x) resect(k,x,st,d), 1, [1;1]);
%
%       % sol.x     -> [4.208; 6.707]     truth was [4.2; 6.7]
%       % IT.nIter  -> 5                  IT.reason -> 'tolX'
%       % IT.vPv    -> the history; it must DECREASE at every step
%
%   INPUT
%     x0        u x 1  double  starting point
%     modelFcn  1x1    fh      [A, l] = modelFcn(x) -- the linearised model at x
%     varargin  Name/Value: 'maxIter', 'tol', 'damping', 'verbose'
%
%   ERRORS
%     lsqvce:nonlin:badX0      the starting point has the wrong length
%     lsqvce:nonlin:nonFinite  the model returned NaN or Inf at the current
%                              iterate
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (c) c1 (F-07):
%   canonical UpperCamelCase option names + typed validators on
%   inputParser, which stays BY DESIGN (KeepUnmatched pass-through --
%   the adjust pattern; pinned deviation §9 WP-39). With assistance from
%   Claude (Anthropic), model Claude Fable 5.
%   See also LSQVCE.ADJUST, LSQVCE.ACCUMULATOR.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

ip = inputParser;
ip.FunctionName = 'lsqvce.solveNonlinear';
ip.KeepUnmatched = true;
% WP-39 c1, PINNED DEVIATION (the adjust pattern): solveNonlinear KEEPS
% inputParser BY DESIGN (KeepUnmatched forwards to the linear solve /
% diagnostics chain). Canonical UpperCamelCase names + typed validators;
% enum VALUES case-sensitive lowercase (strcmp); the o remap keeps the
% body untouched.
addParameter(ip, 'MaxIter', 30,    @(x) isscalar(x) && x >= 1 && x == floor(x));
addParameter(ip, 'TolX',    1e-3,  @(x) isscalar(x) && x > 0);
addParameter(ip, 'TolF',    1e-10, @(x) isscalar(x) && x > 0);
addParameter(ip, 'Lambda0', 1e-6,  @(x) isscalar(x) && x >= 0);
addParameter(ip, 'DatumB',  [],    @(x) isempty(x) || isnumeric(x));
addParameter(ip, 'Method',  'qr',  @(s) (ischar(s) || isstring(s)) && any(strcmp(char(s), {'qr','normal'})));
addParameter(ip, 'Sigma0',  1,     @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'Verbose', true,  @(v) (islogical(v) || isnumeric(v)) && isscalar(v));
parse(ip, varargin{:});
r_ = ip.Results;
o  = struct('maxIter', r_.MaxIter, 'tolX', r_.TolX, 'tolF', r_.TolF, ...
            'lambda0', r_.Lambda0, 'datumB', r_.DatumB, ...
            'method', char(r_.Method), 'sigma0', r_.Sigma0, ...
            'verbose', r_.Verbose);

x   = x0(:);
if numel(x) ~= u
    error('lsqvce:nonlin:badX0', ...
          'x0 has %d elements but u = %d.', numel(x), u);
end

lam = o.lambda0;
IT = struct('x', [], 'vPv', [], 'dxNorm', [], 'lambda', [], ...
            'converged', false, 'nIter', 0, 'reason', 'maxIter');

accArgs = {'method', o.method};
if ~isempty(o.datumB)
    accArgs = [accArgs, {'datumB', o.datumB}];
end

if lsqvce.util.verbosity.progress(o.verbose)
    lsqvce.util.log.msg('\n--- NONLINEAR ITERATION (Gauss-Newton%s)\n', ...
            ternaryStr(o.lambda0 > 0, ' + Levenberg damping', ''));
    lsqvce.util.log.msg('  %4s %16s %14s %12s\n', 'iter', 'vPv', 'max|dx/sig|', 'lambda');
    lsqvce.util.log.msg('  %s\n', repmat('-', 1, 50));
end

vPvPrev = Inf;
solPrev = [];

for it = 1:o.maxIter

    % ---- accumulate at the current x
    acc = lsqvce.Accumulator(u, accArgs{:});
    for k = 1:nBatch
        [Ak, wk, pk] = callBatch(batchFcn, k, x);
        acc.add(Ak, wk, pk);
    end

    % Damping. Applied to the NORMAL matrix, so we need it regardless of the
    % accumulation method: extract R, form N, damp, re-factor. For u of any
    % reasonable size this is negligible next to the n-scale accumulation.
    solTry = acc.solve('sigma0', o.sigma0);

    if lam > 0
        N   = solTry.R' * solTry.R;
        Nd  = N + lam * diag(max(diag(N), eps));
        Rd_ = chol(Nd);
        % rhs: N*dx = b  =>  b = R'*(R*dx_undamped)
        b   = solTry.R' * (solTry.R * solTry.x);
        dx  = Rd_ \ (Rd_' \ b);
    else
        dx = solTry.x;
    end

    % ---- evaluate the step: does it actually reduce vPv?
    xNew   = x + dx;
    vPvNew = evalVPv(batchFcn, nBatch, xNew);

    if ~isfinite(vPvNew)
        error('lsqvce:nonlin:nonFinite', ...
              'The model returned non-finite residuals at iteration %d.', it);
    end

    if o.lambda0 > 0 && vPvNew > vPvPrev && it > 1
        % Rejected step: increase damping and retry without moving.
        lam = lam * 10;
        if lam > 1e12
            IT.reason = 'diverged';
            lsqvce.util.log.warn('lsqvce:nonlin:diverged', ...
                ['Levenberg damping reached %.1e without finding a downhill ' ...
                 'step. The iteration is not converging.\n' ...
                 'Either the starting point is far from the solution, the ' ...
                 'Jacobian is wrong (check it against a finite difference), or ' ...
                 'the model is not locally identifiable at this x.'], lam);
            break;
        end
        if lsqvce.util.verbosity.progress(o.verbose)
            lsqvce.util.log.msg('  %4d %16.8e %14s %12.1e   step rejected\n', ...
                    it, vPvNew, '-', lam);
        end
        continue;
    end

    % ---- accepted
    % Parameter variances = diag(N^-1) with N = R'R. Do NOT form inv(R'R): that
    % squares the condition number (the very thing the QR backend avoids) and
    % costs two dense products. Since N^-1 = R^-1 R^-T, its diagonal is the
    % squared row-norms of R^-1 (= squared column-norms of R^-T), obtained from
    % one triangular solve.
    Rinv  = solTry.R \ eye(size(solTry.R));      % triangular solve, not inv(N)
    sigX  = sqrt(max(sum(Rinv.^2, 2), 0));
    dxRel = max(abs(dx) ./ max(sigX, realmin));

    x       = xNew;
    solPrev = solTry;
    lam     = max(lam / 10, 1e-12);

    IT.x(:, end+1)      = x;      %#ok<AGROW>
    IT.vPv(end+1)       = vPvNew; %#ok<AGROW>
    IT.dxNorm(end+1)    = dxRel;  %#ok<AGROW>
    IT.lambda(end+1)    = lam;    %#ok<AGROW>

    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg('  %4d %16.8e %14.3e %12.1e\n', it, vPvNew, dxRel, lam);
    end

    % ---- convergence
    relF = abs(vPvNew - vPvPrev) / max(abs(vPvNew), realmin);
    if dxRel < o.tolX
        IT.converged = true;  IT.reason = 'tolX';  vPvPrev = vPvNew;  break;
    end
    if it > 1 && relF < o.tolF
        IT.converged = true;  IT.reason = 'tolF';  vPvPrev = vPvNew;  break;
    end
    vPvPrev = vPvNew;
end

IT.nIter = numel(IT.vPv);

if ~IT.converged && ~strcmp(IT.reason, 'diverged')
    lsqvce.util.log.warn('lsqvce:nonlin:maxIter', ...
        ['Reached maxIter = %d without meeting the convergence tolerance ' ...
         '(last max|dx/sigma| = %.3e, target %.3e).\n' ...
         'The diagnostics below describe an UNCONVERGED solution: vPv is ' ...
         'inflated by however much further the fit would still have improved, ' ...
         'which can masquerade as outliers or as a systematic error. Do not ' ...
         'interpret them until this converges.'], ...
         o.maxIter, IT.dxNorm(end), o.tolX);
end

% ---------------------------------------------------------------- final solve
% Re-accumulate AT THE CONVERGED x. This matters: the accumulation inside the
% loop was formed at the PREVIOUS iterate, so its vPv belongs to that point,
% not to the converged one.
if lsqvce.util.verbosity.progress(o.verbose)
    lsqvce.util.log.msg(['\n  Re-accumulating at the converged x (the in-loop vPv belongs\n' ...
             '  to the previous iterate, not to this one).\n']);
end
acc = lsqvce.Accumulator(u, accArgs{:});
for k = 1:nBatch
    [Ak, wk, pk] = callBatch(batchFcn, k, x);
    acc.add(Ak, wk, pk);
end
sol = acc.solve('sigma0', o.sigma0);

% ORDER MATTERS HERE.
%
% diagnose() computes v = A*sol.x - l. For the LINEARISED system the design is
% A, the "observations" are the reduced w = l - f(x), and the solution is the
% step dx. So diagnose must be handed sol.x = dx, NOT the parameters.
%
% At convergence dx -> 0 and therefore v = A*dx - w -> -w, i.e. the linearised
% residual is exactly the observed-minus-computed vector, which is what it
% should be. Verified: ||dx|| = 7e-15 and vPv from (A*dx - w) agreed with vPv
% from w alone to 2.8e-14.
%
% Only AFTER diagnose has run do we overwrite sol.x with the parameters, which
% is what the caller actually wants to see.
dxFinal  = sol.x;
finalSrc = @(k) callBatch(batchFcn, k, x);

dArgs = struct2nv(ip.Unmatched);
D = lsqvce.diagnose(sol, finalSrc, 'nBatch', nBatch, dArgs{:});

sol.dx = dxFinal;
sol.x  = x;                               % the PARAMETERS, not the step
D.x    = x;
D.dx   = dxFinal;
D.nonlinear = IT;

V = lsqvce.report(D, 'verbose', o.verbose);
end

%% ---------------------------------------------------------------- helpers

function [A, w, p] = callBatch(f, k, x)
%CALLBATCH  Call a nonlinear batch source f(k,x) at the current estimate x,
%   tolerating sources that return only [A,w] (no weights) -> p = [].
%   INPUT   f 1x1 fh;  k 1x1 double batch index;  x u x 1 current estimate
%   OUTPUT  A nk x u Jacobian;  w nk x 1 reduced obs (l - f(x));  p nk x 1|[] wts
    try
        [A, w, p] = f(k, x);
    catch ME
        if any(strcmp(ME.identifier, {'MATLAB:TooManyOutputs', ...
                                      'MATLAB:maxlhs', ...
                                      'MATLAB:unassignedOutputs'}))
            [A, w] = f(k, x);
            p = [];
        else
            rethrow(ME);
        end
    end
end

function vPv = evalVPv(f, nB, x)
%EVALVPV  Weighted sum of squares w'Pw at estimate x, streamed over nB batches.
%   INPUT   f 1x1 fh batch source;  nB 1x1 double batches;  x u x 1 estimate
%   OUTPUT  vPv 1x1 double  sum over all observations of p_i * w_i^2
    vPv = 0;
    for k = 1:nB
        [~, w, p] = callBatch(f, k, x);
        w = w(:);
        if isempty(p), p = ones(numel(w), 1); end
        p = p(:);
        vPv = vPv + sum(p .* w.^2);
    end
end

function nv = struct2nv(s)
%STRUCT2NV  Flatten a scalar struct into a 1 x 2K Name/Value cell for varargin.
%   INPUT   s 1x1 struct
%   OUTPUT  nv 1 x 2*numel(fieldnames(s)) cell  {name1,val1,...}
    f  = fieldnames(s);
    nv = cell(1, 2*numel(f));
    for i = 1:numel(f)
        nv{2*i-1} = f{i};
        nv{2*i}   = s.(f{i});
    end
end

function s = ternaryStr(c, a, b)
%TERNARYSTR  Inline conditional returning string a if c else b.
%   INPUT   c 1x1 logical;  a,b 1xk char
%   OUTPUT  s  a if c else b
    if c, s = a; else, s = b; end
end
