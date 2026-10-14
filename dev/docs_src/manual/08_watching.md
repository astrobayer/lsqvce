# Watching a run

## During the run: `OnPass`

`Options.OnPass` is called after every pass of every iterative run (the noise loop,
IRLS, Gauss–Newton) with the pass state, the live accumulator, your declaration and your
`processingOptions`. Plot, `drawnow`, save — or return `false` to stop after the pass:

```matlab
function go = myInspect(ctx)
    res = lsqvce.engine.residuals(ctx.acc, struct('x', ctx.x, 'sigma0', ctx.sigma0, ...
              's0sq', ctx.s0sq, 'Cxx', ctx.Cxx), 'Verbose', false);
    figure(1); plot(res.sys(1).obs(1).t, res.sys(1).obs(1).w, '.'); drawnow;
    fprintf('pass %d: largest change %.3g\n', ctx.pass, ctx.maxDev);
end

R = lsqvce.fit(data, Estimate="noise", OnPass=@myInspect);
```

The fields of `ctx` are listed in the API reference (`lsqvce.Options`, *OnPass*). An
error inside the hook stops the run and arrives as `lsqvce:callUser:threw` with your error
as its cause. A stop by the hook is reported as `R.stopReason = "userStop"`, with the last
pass returned in full.

`Verbose = true` prints one progress line per stage; the default is silent.
`examples/demo_messaging.m` shows the output side by side.

## After the run: one figure per question

```matlab
R.plot("residuals")      % whitened residuals per stream
R.plot("psd")            % the estimated spectra per pool group
R.plot("convergence")    % the components pass by pass
R.plot("contribution")   % which stream determines which parameter
```

Each call returns one figure handle; nothing is exported by the toolbox. `R.stopReason`
says why the loop ended (`"single"`, `"converged"`, `"maxPasses"`, `"userStop"`,
`"diverged"`), and
`R.meta.warnings` lists every warning identifier the run raised.

## Reproducing a run

`R.meta` records the toolbox version, the MATLAB release, the creation time, every random
seed the run resolved and the options as resolved. Saving `R` keeps all of them; a run with
the same data, options and `Seed` reproduces the same numbers.
