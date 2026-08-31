function out = demo_fit(opt)
%DEMO_FIT  The 2.0 user path end to end: a data struct, an Options, fit, report.
%
%   out = demo_fit()
%   out = demo_fit(Verbose=false)
%
%   The four names of the USER tier -- lsqvce.fit, lsqvce.Options,
%   lsqvce.Results, lsqvce.report -- on each of the three kinds of data
%   DATA_CONTRACT.md defines. Nothing here touches the ENGINE tier; every
%   engine call is made by fit.
%
%     (B) in-memory matrices   data.A / data.l / data.p       -> adjust, then
%         the same problem with a planted gross error and Robust="irls"
%     (C) batched callback     data.batchFcn / numBatches (numUnknowns on the
%                              linear form; with x0, u = numel(x0)), and
%         with data.x0 the nonlinear form (Gauss-Newton on exp(a t) + b)
%     (A) streaming            the GRACE-like physical oracle of
%         DEMO_GRAVITY_DATALOADER as a data struct: days, ragged arcs,
%         arc-wise unknowns, coloured noise, three systems -> Model "gls"
%         (one solve) and Model "vce" (stochastic model <-> components)
%
%   INPUT
%     opt.Verbose  1x1 logical  print the verdicts (default true)
%
%   OUTPUT
%     out  1x1 struct  .matrix  .robust  .batch  .nonlinear  .gls  .vce
%                      each an lsqvce.Results; .truth the planted values of
%                      the in-memory problems (x = [1; -2; 0.5], [a; b] =
%                      [0.7; -0.2]) and of the gravity oracle (cfg.xTrue)
%
%   ACCURACY
%     The in-memory problems are seeded (mt19937ar 17 / 5) and closed-form:
%     30 observations with noise ~ 1/sqrt(p), p in [0.5, 1.5], recover x
%     within 3 formal sigma (formal errors ~0.3-0.4); the gross error of 25 sigma is caught by IRLS (its
%     weight factor drops below 0.5) and the robust estimate returns to the
%     clean one within its formal error; the nonlinear form recovers [a; b]
%     to 5e-3 (dev/tests/test_fit.m pins both). The gravity oracle's recovery
%     quality is pinned in DEMO_GRAVITY's ACCURACY block (sigma0 0.9998,
%     stream ratios 0.996 / 1.000 / 1.019); this demo runs its "vce" loop
%     with MaxOuter = 2, enough to show the trajectory, not to converge --
%     the report says so (NOT CONVERGED is a finding, not an error).
%     dev/tests/test_fit.m testDemoFitRuns runs this function and checks every
%     Results it returns.
%
%   See also LSQVCE.FIT, LSQVCE.OPTIONS, LSQVCE.RESULTS, LSQVCE.REPORT,
%            DEMO_GRAVITY_DATALOADER, DEMO_GRAVITY
%
%   Updated 2026-08-31 (CEST, UTC+02:00), v4 SCOPED tables (Level L11,
%   WP-77): the declaration tables carry the scope column (commonParams
%   'common'|'day', localParams 'arc'|'day'|'system'); localsPerArc left the
%   USER surface. With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-74): the nonlinear batched declaration derives u = numel(x0)
%   (numUnknowns dropped there); the streaming declaration comes from
%   demo_gravity_dataLoader's [cfg, data] v4 form. With assistance from
%   Claude (Anthropic), model Claude Fable 5.
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5 -- R4 Level L7, WP-67.
arguments
    opt.Verbose (1,1) logical = true
end
    out = struct();
    say = @(varargin) sayIf(opt.Verbose, varargin{:});

    % ---- (B) in-memory matrices ------------------------------------------
    rs = RandStream('mt19937ar', 'Seed', 17);
    A  = [ones(30, 1), randn(rs, 30, 2)];
    xt = [1; -2; 0.5];
    p  = 0.5 + rand(rs, 30, 1);
    l  = A*xt + randn(rs, 30, 1) ./ sqrt(p);

    say('\n== (B) in-memory: fit(data) with data.A / data.l / data.p ==\n');
    data = struct('A', A, 'l', l, 'p', p);
    out.matrix = lsqvce.fit(data, lsqvce.Options(Verbose=opt.Verbose));
    show(out.matrix, opt.Verbose);

    say('\n== (B) the same with one gross error and Robust="irls" ==\n');
    lBad = l;  lBad(7) = lBad(7) + 25;
    data = struct('A', A, 'l', lBad, 'p', p);
    out.robust = lsqvce.fit(data, lsqvce.Options(Robust="irls", Influence="huber", Verbose=opt.Verbose));
    show(out.robust, opt.Verbose);
    say('   IRLS weight factor of the planted error: %.3f (clean rows ~1)\n', out.robust.source.factor(7));

    % ---- (C) batched callback, linear and nonlinear -------------------------
    say('\n== (C) batched: data.batchFcn over two batches ==\n');
    idx = {1:15, 16:30};
    bf  = @(k) deal(A(idx{k}, :), l(idx{k}), p(idx{k}));
    data = struct('batchFcn', bf, 'numUnknowns', 3, 'numBatches', 2);
    out.batch = lsqvce.fit(data, lsqvce.Options(Verbose=opt.Verbose));
    show(out.batch, opt.Verbose);

    say('\n== (C) nonlinear: l = exp(a t) + b, Jacobian at x, data.x0 ==\n');
    rs2 = RandStream('mt19937ar', 'Seed', 5);
    t   = linspace(0, 1, 40).';
    abt = [0.7; -0.2];
    lob = exp(abt(1)*t) + abt(2) + 1e-3*randn(rs2, 40, 1);
    idx = {1:20, 21:40};
    bfn = @(k, x) deal([t(idx{k}).*exp(x(1)*t(idx{k})), ones(20, 1)], ...
                       lob(idx{k}) - (exp(x(1)*t(idx{k})) + x(2)), []);
    data = struct('batchFcn', bfn, 'numBatches', 2, 'x0', [0.5; 0]);   % u = numel(x0), v3
    out.nonlinear = lsqvce.fit(data, lsqvce.Options(NonlinearMaxIter=15, TolX=1e-9, Verbose=opt.Verbose));
    show(out.nonlinear, opt.Verbose);
    say('   [a; b] = [%.4f; %.4f]  (truth [0.7; -0.2])\n', out.nonlinear.x);

    % ---- (A) streaming: the gravity oracle as a data struct -----------------
    say('\n== (A) streaming: the GRACE-like oracle, Model="gls" ==\n');
    [cfg, data] = demo_gravity_dataLoader();
    out.gls = lsqvce.fit(data, lsqvce.Options(Model="gls", Verbose=opt.Verbose));
    show(out.gls, opt.Verbose);

    say('\n== (A) streaming: Model="vce", two outer passes ==\n');
    w1 = warning('off', 'lsqvce:adjustVCE:notConverged');  % two passes cannot converge:
    w2 = warning('off', 'lsqvce:vce:notConverged');         % the warnings are expected here
    cleaner = onCleanup(@() warning([w1 w2]));              %#ok<NASGU>  restore both on exit
    out.vce = lsqvce.fit(data, lsqvce.Options(Model="vce", GroupBy="obsType", ...
                                              MaxOuter=2, Verbose=opt.Verbose));
    show(out.vce, opt.Verbose);
    say('   sigma per group after each pass:\n');
    if opt.Verbose, disp(out.vce.history.sigma); end

    out.truth = struct('x', xt, 'ab', abt, 'xGravity', cfg.xTrue);
end

%==========================================================================
function show(R, verbose)
    if ~verbose, return, end
    R.report('Verbose', true);
    fprintf('   sigma0 = %.4f\n', R.sigma0);
end

function sayIf(verbose, varargin)
    if verbose, fprintf(varargin{:}); end
end
