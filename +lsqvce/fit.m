function R = fit(task, opts, nv)
%LSQVCE.FIT  Fit a least-squares adjustment (optionally with VCE) to a declared task.
%
%   R = lsqvce.fit(task)                        every default (one GLS solve)
%   R = lsqvce.fit(task, opts)                  a lsqvce.Options object
%   R = lsqvce.fit(task, Name, Value, ...)      the Options properties, straight on fit
%   R = lsqvce.fit(task, opts, Name, Value)     both: the pairs override opts
%
%   THE ONE VERB. You declare the problem as a lsqvce.TaskDefinition, say
%   what to do with lsqvce.Options, and get a lsqvce.Results back:
%
%       task = lsqvce.TaskDefinition(Days=59396:59425);
%       task.addCommon(coeffNames);
%       sys = task.addSystem("grace", Prepare=@prepareOrbit, Arcs=4);
%       sys.addLocal(["x0" "y0" "z0" "vx0" "vy0" "vz0"], Scope="arc");
%       sys.addStream("kbrr", Build=@buildKbrr, Noise=lsqvce.Noise.white(1e-7));
%       R = lsqvce.fit(task, Estimate="noise", Pool="obsType", OnPass=@myInspect);
%       R.report;  R.plot("psd");  T = R.toTable(task);
%
%   THREE KINDS OF TASK, made by the constructor that fits the data, all on
%   ONE accumulator:
%     streaming   lsqvce.TaskDefinition(Days=...) and its add-calls: arc-wise
%                 data read day by day through each stream's Build, local
%                 parameters eliminated per system, the noise declared per
%                 stream (LSQVCE.TASKDEFINITION, LSQVCE.SYSTEMDEFINITION)
%     in memory   lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p)
%     batched     lsqvce.TaskDefinition.fromBatches(@batchFcn, NumBatches=K,
%                 NumUnknowns=u), or X0=x0 for the nonlinear form
%   WHAT RUNS (opts.Estimate / opts.Robust / X0):
%     Estimate "none"    one generalised solve under the declared laws
%     Estimate "noise"   the joint variance-component estimator: one spectral
%                        shape per Pool group, one level per stream, per-arc
%                        factors when ArcFactors; a prior group declared with
%                        Estimate="scale" gets its sigma by Helmert /
%                        Koch-Kusche, one step per pass
%     Robust "huber" | "igg3" | "bisquare"   IRLS on the same loop (inside the
%                        joint loop with Estimate "noise")
%     X0                 Gauss-Newton with damping (fromBatches only)
%   REGULARISATION is declared on the task, never an option: task.addPrior on
%   the common parameters, sys.addPrior on a system's locals.
%   WATCHING A RUN: Options.OnPass fires after every pass with the pass state,
%   the live accumulator, ctx.data (the compiled declaration) and every
%   system's Settings.
%
%   INPUT
%     task   1x1 lsqvce.TaskDefinition   the declaration; fit compiles it with
%                          task.declaration() and fits that with
%                          LSQVCE.ENGINE.FITDECLARATION
%     opts   1x1 lsqvce.Options   default lsqvce.Options()
%     Name, Value          any property of lsqvce.Options
%   OUTPUT
%     R      1x1 lsqvce.Results   x, sigma0, sol (Cxx, dof, vtv), diag, acc, history,
%                                 noise (S, Sgroups, C, Cgroups, stream, model, priors,
%                                 checks), stopReason, meta (options as resolved,
%                                 warnings, seeds, version); methods report, summary,
%                                 convergence, residuals, contribution, contributions,
%                                 compare, plot, toTable, coefficients, reliability,
%                                 systematic, snoop, ellipse, redundancy (LSQVCE.RESULTS)
%
%   ACCURACY
%     fit adds no arithmetic of its own: it compiles the task (bookkeeping;
%     test_taskDefinition pins fit(task) BITWISE equal to the fit of the
%     equivalent hand-written struct) and returns LSQVCE.ENGINE.FITDECLARATION's
%     result -- see its ACCURACY for the estimate, the verdict and the plug-in
%     formal errors of an estimated noise model.
%
%   ERRORS
%     lsqvce:fit:structDeclaration   a struct was passed: since 7.0.0 the
%                                    declaration is a lsqvce.TaskDefinition
%                                    (RELEASE_NOTES, "Migrating a declaration
%                                    from v6 to v7"); the struct is the engine
%                                    tier's, LSQVCE.ENGINE.FITDECLARATION
%     MATLAB:validators:mustBeA      neither a struct nor a lsqvce.TaskDefinition
%     lsqvce:TaskDefinition:incomplete / badArcs / duplicateName   the task does
%                                    not compile (task.declaration())
%     lsqvce:fit:*, lsqvce:makeLoader:*, lsqvce:DayContext:*   raised by the
%                                    compiled declaration's checks and while the
%                                    days are read; LSQVCE.ENGINE.FITDECLARATION
%                                    lists them
%
%   See also LSQVCE.TASKDEFINITION, LSQVCE.OPTIONS, LSQVCE.RESULTS, LSQVCE.LAYOUT,
%            LSQVCE.ENGINE.FITDECLARATION
%
%   Created 2026-11-03 (CET, UTC+01:00) with assistance from Claude (Anthropic),
%   model Claude Opus 5.5 -- level v7, WP-162: the user tier of fit takes a
%   lsqvce.TaskDefinition only; the fit of a struct is
%   LSQVCE.ENGINE.FITDECLARATION, where the body of fit (2026-08-31, R4 L2
%   WP-56, to 2026-11-02, WP-160) now lives.
arguments
    task (1,1) {mustBeA(task, ["struct", "lsqvce.TaskDefinition"])}
    opts (1,1) lsqvce.Options = lsqvce.Options()
    nv.?lsqvce.Options
end
    if isstruct(task)
        error('lsqvce:fit:structDeclaration', ...
              ['lsqvce.fit takes a lsqvce.TaskDefinition since 7.0.0, not a data struct. ' ...
               'Declare the task with lsqvce.TaskDefinition (RELEASE_NOTES.md, "Migrating a ' ...
               'declaration from v6 to v7"), or fit a compiled struct on the engine tier with ' ...
               'lsqvce.engine.fitDeclaration(data, ...).']);
    end
    args = namedargs2cell(nv);
    R = lsqvce.engine.fitDeclaration(task.declaration(), opts, args{:});
end
