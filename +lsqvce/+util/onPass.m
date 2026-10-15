function go = onPass(hook, ctxExtra, ctx)
%ONPASS  Fire the user's per-pass inspection hook (R5 M5 5c').
%   go = lsqvce.util.onPass(hook, ctxExtra, ctx)
%
%   THE one place every iterative driver (the joint estimator, IRLS,
%   Gauss-Newton) calls Options.OnPass at the end of a pass. The context the
%   hook receives is the driver's ctx (pass, solution, noise / weight / step
%   state, and .acc -- the LIVE accumulator, so lsqvce.engine.residuals(ctx.acc,
%   ...) and every other engine call are yours inside the hook) merged with
%   ctxExtra, which lsqvce.fit fills with .data (the declaration as given) and
%   .processingOptions (per system, as declared). The hook is called through
%   lsqvce.util.callUser: a throwing hook is NOT swallowed (it arrives as
%   lsqvce:callUser:threw with the original as cause), and its output count is
%   asked once. Return nothing (or true) to continue, false to stop the loop
%   after this pass -- the driver then returns the pass with stopReason
%   "userStop".
%
%   INPUT   hook     [] | function_handle   Options.OnPass
%           ctxExtra struct                 fields merged into ctx (fit's data, processingOptions)
%           ctx      struct                 the driver's per-pass context
%   OUTPUT  go       1x1 logical            true = continue, false = stop
%   ACCURACY
%           Not applicable: it merges a context struct and calls the hook.
%           WP-139.
%   ERRORS  lsqvce:callUser:threw  the hook threw (its error is the cause)
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, "Not applicable" with its reason. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M5, WP-87 round 5c'
%   (the per-pass inspection seam of the project's opening brief, restored).
    go = true;
    if isempty(hook), return, end
    f = fieldnames(ctxExtra);
    for k = 1:numel(f), ctx.(f{k}) = ctxExtra.(f{k}); end
    [n, out] = lsqvce.util.callUser(hook, 0, 1, sprintf('OnPass (pass %d)', ctx.pass), ctx);
    if n == 1 && ~isempty(out) && (islogical(out) || isnumeric(out)) && isscalar(out) && ~out
        go = false;
    end
end
