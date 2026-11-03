function out = demo_fit(opt)
%DEMO_FIT  The user path end to end: a TaskDefinition, an Options, fit, report.
%
%   out = demo_fit()
%   out = demo_fit(Verbose=false)
%
%   Five names of the USER tier -- lsqvce.TaskDefinition, lsqvce.fit,
%   lsqvce.Options, lsqvce.Results, lsqvce.report -- on each of the three
%   kinds of task the API reference defines. lsqvce.layout, which announces
%   the column layout of A and B so that priors can be declared in it, is not
%   needed here, where the prior below is written straight in the columns of A
%   (help lsqvce.layout, and the Manual, chapter Priors and regularisation).
%   Nothing here touches the ENGINE tier; every engine call is made by fit.
%
%     (B) in-memory matrices   lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p)
%         -> one solve, then the same problem with a planted gross error and
%         Robust="huber", and a free network whose DATUM is a declared prior
%         (task.addPrior with L spanning the null space): regularisation is
%         declared, never an option
%     (C) batched callback     lsqvce.TaskDefinition.fromBatches(@batchFcn,
%         NumBatches=K, NumUnknowns=u) on the linear form, and with X0=x0
%         (u = numel(x0)) the nonlinear form (Gauss-Newton on exp(a t) + b)
%     (A) streaming            the GRACE-like physical oracle, declared by
%         DEMO_GRAVITY_DATALOADER: days, ragged arcs, arc-wise unknowns,
%         coloured noise, three systems -> Estimate "none" (one solve) and
%         Estimate "noise" (the joint variance-component estimator, Pool
%         "stream", one level per stream), watched pass by pass through
%         Options.OnPass
%
%   INPUT
%     opt.Verbose  1x1 logical  print the verdicts (default true)
%
%   OUTPUT
%     out  1x1 struct  .matrix  .robust  .datum  .batch  .nonlinear  .gls  .vce
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
%     to 5e-3 (dev/tests/test_fit.m pins both). The free network's datum is
%     exact: with L spanning the null space, L*x = x0 holds to rounding for ANY
%     prior sigma and the solution is the minimum-constraint one
%     (dev/tests/test_fit.m testDatumIsADeclaredPrior pins both against the
%     pseudo-inverse). The gravity oracle's recovery quality is pinned in
%     DEMO_GRAVITY's ACCURACY block; no figure of it is repeated here, because
%     the levels that come back depend on which estimator runs. The counted
%     reference for the planted stream weights is dev/tests/test_vcePolicy.m
%     testRecoversPlantedRelativeWeights (planted 1 / 2 / 0.5, recovered to
%     25 % relative in six passes of the joint estimator). This demo runs its
%     "vce" loop with MaxPasses = 2, enough to show the trajectory, not to
%     converge -- the report says so (NOT CONVERGED is a finding, not an error).
%     dev/tests/test_fit.m testDemoFitRuns runs this function and checks the
%     six Results of the three kinds of task.
%
%   ERRORS
%     none of its own. A declaration error is raised by the add-call that
%     makes it (lsqvce:TaskDefinition:*). The free network WITHOUT its datum
%     is refused by the engines (lsqvce:qr:singularR); the demo catches that
%     refusal and prints its identifier.
%
%   See also LSQVCE.TASKDEFINITION, LSQVCE.FIT, LSQVCE.OPTIONS, LSQVCE.RESULTS,
%            LSQVCE.REPORT, DEMO_GRAVITY_DATALOADER, DEMO_GRAVITY
%
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: kinds B and C declared through lsqvce.TaskDefinition. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-95: the USER-tier census names
%   the fifth name (lsqvce.layout, not exercised here), the stale gravity
%   figures give way to the counted reference, and a free network's datum is
%   demonstrated as a declared prior.
%   With assistance from Claude (Anthropic), model Claude Opus 5.
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
%
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-124: references to the data
%   contract point at the API reference. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
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

    say('\n== (B) in-memory: fit(task) with TaskDefinition.fromMatrices(A, l, Weights=p) ==\n');
    task = lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p);
    out.matrix = lsqvce.fit(task, lsqvce.Options(Verbose=opt.Verbose));
    show(out.matrix, opt.Verbose);

    say('\n== (B) the same with one gross error and Robust="huber" ==\n');
    lBad = l;  lBad(7) = lBad(7) + 25;
    task = lsqvce.TaskDefinition.fromMatrices(A, lBad, Weights=p);
    out.robust = lsqvce.fit(task, lsqvce.Options(Robust="huber", Verbose=opt.Verbose));
    show(out.robust, opt.Verbose);
    say('   IRLS weight factor of the planted error: %.3f (clean rows ~1)\n', out.robust.history(end).factor(7));

    say('\n== (B) a free network: the datum DECLARED as a prior on the null space ==\n');
    % Regularisation is declared, never an option. task.addPrior(name, L=L,
    % Value=x0, Sigma=s) is one pseudo-observation group L*x ~ x0 of weight
    % 1/s^2; for a datum L spans the null space, and then any sigma gives the
    % same minimum-constraint x. An in-memory task takes its priors by L only.
    pairs = [1 2; 2 3; 3 4; 4 1; 1 3];                  % height differences only
    Af = zeros(size(pairs, 1), 4);
    for k = 1:size(pairs, 1), Af(k, pairs(k, 1)) = -1;  Af(k, pairs(k, 2)) = 1; end
    lf = Af*[0; 1.5; -0.7; 2.2] + 0.01*randn(rs, size(pairs, 1), 1);
    Ld = ones(1, 4);                                    % the null space: a common level
    task = lsqvce.TaskDefinition.fromMatrices(Af, lf);  % unit weights
    task.addPrior("datum", L=Ld, Value=0, Sigma=1);
    out.datum = lsqvce.fit(task, lsqvce.Options(Verbose=opt.Verbose));
    say('   the declared datum holds: L*x = %.1e\n', Ld*out.datum.x);
    try
        lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(Af, lf), lsqvce.Options(Verbose=false));
    catch ME
        say('   without it the same problem is refused by name: %s\n', ME.identifier);
    end

    % ---- (C) batched callback, linear and nonlinear -------------------------
    say('\n== (C) batched: TaskDefinition.fromBatches over two batches ==\n');
    idx = {1:15, 16:30};
    bf  = @(k) deal(A(idx{k}, :), l(idx{k}), p(idx{k}));
    task = lsqvce.TaskDefinition.fromBatches(bf, NumBatches=2, NumUnknowns=3);
    out.batch = lsqvce.fit(task, lsqvce.Options(Verbose=opt.Verbose));
    show(out.batch, opt.Verbose);

    say('\n== (C) nonlinear: l = exp(a t) + b, Jacobian at x, X0=x0 ==\n');
    rs2 = RandStream('mt19937ar', 'Seed', 5);
    t   = linspace(0, 1, 40).';
    abt = [0.7; -0.2];
    lob = exp(abt(1)*t) + abt(2) + 1e-3*randn(rs2, 40, 1);
    idx = {1:20, 21:40};
    bfn = @(k, x) deal([t(idx{k}).*exp(x(1)*t(idx{k})), ones(20, 1)], ...
                       lob(idx{k}) - (exp(x(1)*t(idx{k})) + x(2)), []);
    task = lsqvce.TaskDefinition.fromBatches(bfn, NumBatches=2, X0=[0.5; 0]);   % u = numel(x0)
    out.nonlinear = lsqvce.fit(task, lsqvce.Options(MaxPasses=15, Tolerance=1e-9, Verbose=opt.Verbose));
    show(out.nonlinear, opt.Verbose);
    say('   [a; b] = [%.4f; %.4f]  (truth [0.7; -0.2])\n', out.nonlinear.x);

    % ---- (A) streaming: the gravity oracle, as its loader declares it --------
    say('\n== (A) streaming: the GRACE-like oracle ==\n');
    % the declaration is DEMO_GRAVITY_DATALOADER's (its help says how it is
    % written); it goes to fit as it comes
    [cfg, gravity] = demo_gravity_dataLoader();   % its lsqvce.TaskDefinition
    out.gls = lsqvce.fit(gravity, lsqvce.Options(Verbose=opt.Verbose));
    show(out.gls, opt.Verbose);

    say('\n== (A) streaming: Estimate="noise", Pool="stream", ArcFactors=false, two passes ==\n');
    w1 = warning('off', 'lsqvce:adjustVCEjoint:notConverged');  % two passes cannot converge:
    cleaner = onCleanup(@() warning(w1));                       %#ok<NASGU>  the warning is expected
    % Watch the loop. OnPass fires after EVERY pass with the pass state, the
    % live accumulator, your declaration and every system's settings -- plot,
    % save, or return false to stop. Here: one line per pass.
    watch = @(ctx) say('   pass %d: max|alpha2-1| = %.3g, levels = %s\n', ...
                       ctx.pass, ctx.maxDev, mat2str(ctx.sigma, 3));
    out.vce = lsqvce.fit(gravity, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, ...
                                              MaxPasses=2, Verbose=opt.Verbose, OnPass=watch));
    show(out.vce, opt.Verbose);
    say('   sigma per stream after each pass:\n');
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
