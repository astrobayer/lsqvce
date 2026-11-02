# The workflow, step by step

This chapter is the checklist for your own data: the decisions in the order you make
them, each with the chapter that treats it in depth. If you have not yet run the tutorial
(*Your first adjustment*), do that first; it is the same path on a small case.
`examples/demo_workflow.m` walks the path as a script to run section by section.

| step | you decide | you write |
|---|---|---|
| 1 | which kind of problem | — |
| 2 | the parameters and their scopes | `task.addCommon`, `sys.addLocal` |
| 3 | the systems, their streams and their starting noise | `task.addSystem`, `sys.addStream` |
| 4 | how one day of one stream becomes rows | its `Build` function |
| 5 | whether the streams of a day share expensive work | the system's `Prepare` (optional) |
| 6 | — run once under the declared noise | |
| 7 | — estimate the noise | `lsqvce.Options` |
| 8 | regularisation, if needed | `addPrior` |
| 9 | — read the report and act on it | |
| 10 | — take the results out | |

## Step 1 — Decide which kind of problem you have

| your problem | kind | you write |
|---|---|---|
| the design matrix fits in memory, one block of unknowns | **matrices** | `lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p)` |
| too large for memory, or nonlinear, one block of unknowns | **batched** | `lsqvce.TaskDefinition.fromBatches(@batchFcn, NumBatches=K, NumUnknowns=u)` (or `X0=x0`) |
| data arrive per day and arc, several instruments, nuisance parameters per arc or day, correlated noise | **streaming** | `lsqvce.TaskDefinition(Days=...)` and its `add` calls, with your `Build` functions |

The first two are complete after one line:

```matlab
task = lsqvce.TaskDefinition.fromMatrices(A, l, Weights=1 ./ sigma.^2);   % weights, not sigmas
R = lsqvce.fit(task);
R.report
```

The rest of this chapter is the streaming kind, the one the toolbox is built for.

## Step 2 — List the parameters

Declare every parameter by name and scope. `addCommon` declares the parameters every
observation can see (for gravity recovery: the spherical-harmonic coefficients); a
system's `addLocal` declares the nuisance parameters that live on one arc, one day or the
whole system (initial states, accelerometer calibration, empirical terms):

```matlab
task = lsqvce.TaskDefinition(Days=1:30);
task.addCommon(myCoefficientNames);                    % Scope="common": one value per run
sys = task.addSystem("mission", Settings=struct('dataDir', '/data/mission'), Arcs=4);
sys.addLocal(["bias" "drift"], Scope="arc");           % one each per arc
sys.addLocal("accScale", Scope="day");                 % one per day
```

The toolbox derives every count and every column position from the names and scopes. See
*Declaring your problem* for the scopes, validity intervals and `Info` columns.

## Step 3 — Declare the streams and their noise

```matlab
sys.addStream("range", Build=@buildRange, Noise=lsqvce.Noise.white(2e-7));   % white, 2e-7
```

Each `addStream` is one observation type: its name, its `Build` function and its noise
law. Start with `lsqvce.Noise.white` and a realistic level; the estimator can find the
colour (*Choosing the stochastic model*). The system's `Settings` (step 2) is a struct of
your own that reaches your functions as `ctx.Settings`.

## Step 4 — Write `Build`

One function per stream returns the arcs of one day:

```matlab
function out = buildRange(day, ctx)
    arcs = readMyArcs(ctx.Settings.dataDir, day, "range");   % your reader, deterministic
    out = ctx.newDay();
    for a = 1:numel(arcs)
        m = numel(arcs(a).omc);
        out.addArc(myPartials(arcs(a)), arcs(a).omc, arcs(a).tSeconds, ...   % A, l, t [s]
                   Local={["bias" "drift"], [ones(m, 1), arcs(a).tRel]});    % by name
    end
end
```

Return raw, unwhitened rows; never scale them by your weights. `A` has one column per
common parameter, in `addCommon` order; the local partials go **by name** in `Local`, and
the toolbox places them in this arc's (or this day's) columns.

What the context gives you, in plain words:

| `ctx.` | what you get |
|---|---|
| `Day` | the day being built, one entry of the calendar |
| `System`, `Stream` | the names of the system and the stream |
| `Settings` | the system's `Settings`, unchanged |
| `Prepared` | what the system's `Prepare` returned for this day (step 5) |
| `ArcIds` | the scheduled arcs of the day |
| `NumCommon`, `NumLocal` | the number of columns of `A` and of the local design |
| `Common`, `Local` | tables of every parameter: name, scope, whether it exists on this day, `Info` |
| `cols(names)` | the columns of named parameters, for a local design you build as one matrix `B` |

`Build` must be **deterministic**: the toolbox calls it at least twice per day and pass,
and must get the same rows back each time.

## Step 5 — Share the expensive work of a day: `Prepare`

`Build` is called once per day **and stream**. If your system has several streams that
rest on one expensive computation — an orbit integration with its partials, a large file
both streams read, an Earth-rotation matrix for every epoch — that computation would run
once per stream. Move it into the system's `Prepare`:

```matlab
sys = task.addSystem("grace", Prepare=@myPrepare, Settings=struct('dataDir', '/data'));
sys.addStream("pos", Build=@buildPos, Noise=lsqvce.Noise.white(0.02), Axes=3);
sys.addStream("kbr", Build=@buildKbr, Noise=lsqvce.Noise.white(1e-7));

function P = myPrepare(day, ctx)
    P.orbit = integrateOrbit(ctx.Settings.dataDir, day);   % the expensive part, once
end

function out = buildKbr(day, ctx)
    out = rangeRateArcs(ctx.Prepared.orbit, ctx);          % what myPrepare returned today
end
```

What to know:

- It runs **once per build of a day**, before the first `Build` of that day, and its
  result `P` reaches every `Build` of that day as `ctx.Prepared`.
- Its context carries the day, the system's `Settings` and its scheduled arcs
  (`ctx.ArcIds`); `ctx.Stream` is empty.
- `Settings` must not have a field named `prepared` of its own: that name carries the
  result (`lsqvce:makeLoader:preparedField`, reported by `fit` under its own prefix).
- Only the current day's result is kept. A day that is built again (the default
  `CacheDays = 1` builds each day twice per pass) runs `Prepare` again, so it must be
  deterministic too. `CacheDays = Inf` builds every day only once per run, if memory
  allows.
- If it throws, the error is `lsqvce:makeLoader:prepareDay`, naming the day, with your
  error attached as the cause.

A system with one stream gains nothing from `Prepare`; leave it out. The tutorial (*Your
first adjustment*) and `template/kbr_variational.m` use it.

## Step 6 — First run: one solve under the declared model

```matlab
R = lsqvce.fit(task);          % Estimate = "none"
R.report
```

Check two things before anything else: that the run completes (a declaration mistake is
refused by the `add` call that makes it, or by `fit` naming the system and stream), and
that `R.sigma0` is plausible. A σ₀ far from 1 means the declared accuracies are wrong in
scale — fix that before estimating anything on top of it.

## Step 7 — Estimate the stochastic model

```matlab
opts = lsqvce.Options(Estimate="noise", Pool="stream", Seed=1);
R = lsqvce.fit(task, opts);
R.report
R.plot("psd"); R.plot("convergence");
```

The loop re-solves until the noise laws settle. `Pool` decides which streams share one
spectral shape; `ArcFactors` (default on) gives every arc its own variance factor; `Seed`
makes the probe-based redundancies reproducible. What a stream lets the estimator change
is its `Estimate` setting on `addStream` (`"law"`, `"scale"` or `"none"`).

## Step 8 — Add regularisation if the problem needs it

A free network needs a datum; a high-degree gravity field needs a Kaula rule. Both are
declared as priors on named parameters:

```matlab
task.addPrior("kaula", On="^[CS]", Sigma=1e-10, Estimate="scale");
```

See *Priors and regularisation*.

## Step 9 — Read the report and act on it

`R.report` prints a verdict and the evidence behind it: the global test, flagged
observations, the reliability, the whiteness of every stream. Each adverse verdict has a
standard next step — snooping, a richer noise law, a functional model correction. See
*Reading the results and the report*.

## Step 10 — Take the results out

```matlab
T = R.coefficients();                  % estimates, formal errors, confidence bounds
P = R.toTable(task);                   % every parameter by name, day and arc, with its estimate
save('run.mat', 'R');                  % R.meta records version, seeds and options
```
