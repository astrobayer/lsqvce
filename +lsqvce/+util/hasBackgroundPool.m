function tf = hasBackgroundPool()
%HASBACKGROUNDPOOL  True when MATLAB's backgroundPool can actually start.
%   tf = lsqvce.util.hasBackgroundPool()
%
%   The ONE availability probe for the background pool (rule 3): used by
%   lsqvce.engine.PrefetchLoader to decide prefetch-vs-passthrough and by the test
%   suite to guard pool-dependent tests. It PROBES BY USE - it calls
%   backgroundPool inside try/catch - because metadata probes are wrong
%   twice over: exist('backgroundPool','file') == 2 fails since the function
%   ships packaged/built-in (exist returns 5 or 0, never 2), and the
%   class+file exist pair is FALSE-NEGATIVE under `matlab -batch` (measured
%   2026-08-18, R2026a Update 4, Windows 11: license=1, ver=1, class=0,
%   file=0 while backgroundPool itself starts fine with 16 workers - WP-30
%   postflight finding). The silent effect of the false negative was
%   prefetch degrading to passthrough in EVERY batch run. Calling the
%   function is the only probe that answers the actual question; the lazy
%   singleton it starts is the very pool the caller is about to use, and on
%   machines without it the cost is one caught error.
%
%   INPUT
%     (none)
%
%   OUTPUT
%     tf   1x1 logical   true iff the Parallel Computing Toolbox is licensed
%                        and installed AND backgroundPool starts.
%
%   ERRORS
%     (none - a pool that cannot start IS the false branch)
%
%   EXAMPLE
%     if lsqvce.util.hasBackgroundPool()
%         fut = parfeval(backgroundPool, @myLoader, 1, d);
%     end
%
%   See also LSQVCE.ENGINE.PREFETCHLOADER, BACKGROUNDPOOL, PARFEVAL
%
%   Created 2026-08-18 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    tf = license('test', 'Distrib_Computing_Toolbox') && ...
         ~isempty(ver('parallel'));
    if ~tf
        return
    end
    try
        backgroundPool;
    catch ME                                     %#ok<NASGU>
        tf = false;      % capability probe: "no" IS the documented answer (R5 M1)
    end
end
