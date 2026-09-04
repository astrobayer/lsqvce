# Watching a run

`Options.OnPass` (chapter 4) is called after every pass with the pass state, the live
accumulator, your declaration and your `processingOptions` — plot, `drawnow`, save,
or return `false` to stop. Afterwards `R.plot("residuals" | "contribution" | "psd" |
"convergence", Stream=…, Visible=…)` returns one figure per question; nothing is exported
by the toolbox. `R.meta.warnings` lists every warning identifier of the run;
`R.stopReason` says why the loop stopped (`"single"`, `"converged"`, `"maxPasses"`,
`"userStop"`).

```matlab
function go = myInspect(ctx)
    res = lsqvce.engine.residuals(ctx.acc, struct('x', ctx.x, 'sigma0', ctx.sigma0, ...
                                                  's0sq', ctx.s0sq, 'Cxx', ctx.Cxx), 'Verbose', false);
    figure(1); plot(res.sys(1).obs(1).t, res.sys(1).obs(1).w, '.'); drawnow;
    go = ctx.maxDev > 1e-3;                 % or return nothing to continue
end
```
