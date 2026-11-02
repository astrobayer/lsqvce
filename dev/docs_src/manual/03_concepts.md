# Key concepts

This chapter introduces the vocabulary the rest of the Manual uses. It has no code to run.
Read it once before the tutorial. The words in **bold** are the names you will meet as
method names, options and chapter titles.

## The adjustment in one paragraph

You have observations `l` that depend on unknown parameters `x` through a linear (or
linearised) model, `l + v = A x`, where `v` are the errors. Least squares picks the `x`
that makes the weighted errors smallest. The weights are the inverse of the noise
covariance of the observations. Usually you do not know that covariance well. *Variance
component estimation* (VCE) estimates it from the residuals, re-weights, and solves
again until the weights settle. lsqvce does all of that in one call, `lsqvce.fit`.

## Three kinds of problem

Every problem is declared as one `lsqvce.TaskDefinition`, made in one of three ways:

| kind | made with | use it when |
|---|---|---|
| **matrices** | `TaskDefinition.fromMatrices(A, l, Weights=p)` | the whole design matrix fits in memory |
| **batched** | `TaskDefinition.fromBatches(@batchFcn, NumBatches=K, NumUnknowns=u)` (or `X0=x0`) | it does not fit, or the model is nonlinear; your function returns one block of rows at a time |
| **streaming** | `TaskDefinition(Days=...)` and its `add` calls | the data come per day and per arc, from several instruments, with nuisance parameters and correlated noise |

The first two are ordinary least squares with diagnostics. The streaming kind is what the
toolbox is built for, and the rest of this chapter describes it.

## Parameters: common and local

Ask of every unknown: *over what span is it constant, and who sees it?*

- **Common parameters** (`task.addCommon`) are seen by all data, for example the
  gravity-field coefficients or station coordinates. They are what you want to know.
- **Local parameters** (`sys.addLocal` of a system) belong to one **arc**, one **day** or
  one **system**, for example an initial orbit state, an accelerometer bias, or a
  per-arc offset. They are nuisance parameters. The toolbox eliminates them arc by arc
  while it reads the data, so their number does not cost memory.

You declare both **by name and scope**: `Scope="common"`, `"day"`, `"arc"` or `"system"`
says over which span the parameter is constant. A declared name is a template: the
toolbox instantiates it once per scope (an `"arc"` local once per arc) and counts the
columns itself. Your code refers to parameters by name; you never compute a column index.

## Data: systems, days, arcs and streams

```
task = lsqvce.TaskDefinition(Days=1:30)       the calendar
├── addCommon(names)                            what everybody estimates
├── addPrior(name, ...)                         priors on the common parameters
└── addSystem(name, Prepare=@p)                 one per satellite, receiver, instrument group
     ├── Prepare                                (optional) once per day: what its streams share
     ├── addLocal(names, Scope=...)             its nuisance parameters
     ├── addPrior(name, ...)                    priors on its local parameters
     └── addStream(name, Build=@b, Noise=...)   one per observation type  ──►  the streams
          └── Build                             once per day: the rows of this stream
```

- A **system** is one source of data with its own nuisance parameters, for example one
  satellite, or a satellite pair whose inter-satellite ranging depends on both orbits.
- A **day** is the unit of reading. The toolbox asks for the data of one day at a time.
- An **observation type** is a kind of measurement within a system, for example `"pos"`
  and `"kbr"`. Each `addStream` declares one; each (system, observation type) pair is one
  **stream**, with its own `Build` function and its own noise law.
- An **arc** is one continuous piece of a stream on a day. The noise is correlated
  within an arc and independent between arcs. `Build` adds the day's arcs one by one.

## The two functions you write

You do not hand the toolbox your data. You hand it functions that produce the data on
demand. This is what keeps memory small, but it has one consequence: **the toolbox calls
your functions more than once for the same day**, so they must return the same result
every time.

**`Build`** (one per stream, required) returns the arcs of one day of its stream:

```matlab
out = build(day, ctx)            % given to addStream(..., Build=@build)
%   out = ctx.newDay();
%   out.addArc(A, l, t, Local={names, P})   once per arc:
%      A  design of the common parameters   (m x nc, columns in addCommon order)
%      l  observed minus computed           (m x 1)
%      t  epochs in seconds                 (m x 1)
%      Local  the partials of local parameters, by name
```

**`Prepare`** (one per system, optional) does the work that every stream of a day
shares, *once* per day, before the first `Build` of that day. Its result reaches every
`Build` of the day as `ctx.Prepared`:

```matlab
P = prepare(day, ctx)            % given to addSystem(..., Prepare=@prepare)
```

Without `Prepare`, a system with two streams computes anything they share twice per day:
an orbit integration, a large file read, a coordinate rotation. With it, the shared work
happens once per build of the day. The tutorial uses it; chapter *Declaring your problem*
shows when it pays.

Both receive the same two arguments: the day and a `lsqvce.DayContext`, `ctx`, which
carries the names of the system and stream, your `Settings`, the arc schedule of the day
(`ctx.ArcIds`) and the parameter names, and finds columns by name (`ctx.cols`).

## What `fit` does with them

For a run with `Estimate = "noise"`, `fit` repeats one **pass** until the noise settles:

```
for every pass
    for every system and every day
        P = Prepare(day, ctx)                     once per day (if declared)
        for every stream
            out = Build(day, ctx)                 the day's arcs, with ctx.Prepared = P
            whiten every arc with the stream's current noise law
            eliminate the local parameters of the arc
            add the arc to the normal equations of the common parameters
    solve for the common parameters
    read the days again: residuals, redundancies
    update the noise law of every stream, its level, and every arc factor
until no variance component changes by more than Tolerance
diagnose: outlier tests, reliability, whiteness  ->  R.report
```

With `Estimate = "none"` (the default) the loop runs once: one solve under the noise you
declared. Each day is read at least twice per pass, once to accumulate and once for the
residuals. `Options.CacheDays` decides how many built days are kept in memory between
the two reads.

## The noise model

Each stream has a **noise law** (`addStream(..., Noise=lsqvce.Noise.<law>(...))`): how
its errors correlate from epoch to epoch (the law) and how large they are (its `sigma`).
On top of that, every arc can have its own **arc factor**, which catches one bad pass
without spoiling the whole stream. With `Estimate = "noise"` the toolbox estimates all of
them, as far as the stream's `Estimate` setting allows. You start from an honest guess,
usually white noise (`lsqvce.Noise.white(sigma)`) with a realistic `sigma`.

The result of all weights together is summarised by **σ₀** (`R.sigma0`), the a-posteriori
standard deviation of unit weight. A σ₀ near 1 means the declared and the observed noise
agree.

## The result and the verdict

`fit` returns one `lsqvce.Results` object, `R`. It holds the estimate `R.x`, its
covariance, the estimated noise model, the per-observation diagnostics, the history of
the passes, and a record of how the run was made (`R.meta`). `R.report` reads all of this
and prints a verdict in five levels:

| level | question |
|---|---|
| 0 integrity | does the arithmetic hold? |
| 1 global | are model and noise consistent with the data? |
| 2 localise | which observations are suspect? |
| 3 reliability | which errors could go undetected? |
| 4 parameters | how well is each parameter known? |

The verdict is one of `PASS`, `PASS WITH WARNINGS`, `SUSPECT` and `INVALID` (see
*Reading the results and the report*).

## Glossary

| term | meaning |
|---|---|
| arc | a continuous piece of one stream on one day; one `addArc` of its `Build` |
| arc factor | a variance factor of one arc, estimated on top of the stream's level |
| Build | the function of a stream that returns its arcs of one day (`addStream`) |
| common parameter | an unknown seen by all data (`task.addCommon`) |
| context | the `lsqvce.DayContext` your Prepare and Build receive: names, settings, schedule, columns |
| day | the unit in which the toolbox reads data |
| local parameter | a nuisance unknown of one arc, day or system, eliminated on the fly |
| observation type | a kind of measurement within a system; the name of an `addStream` |
| pass | one solve plus one update of the noise model |
| pool | streams that share one estimated spectral shape (`Options.Pool`) |
| Prepare | the function of a system that computes, once per day, what its streams share |
| prior | pseudo-observations on the parameters (`addPrior`): datum, regularisation |
| redundancy | how much of an observation is checked by the others (0 to 1) |
| `Sigma` | per-epoch accuracy delivered with an arc (`addArc(..., Sigma=)`), relative to the law's `sigma` |
| σ₀ | a-posteriori standard deviation of unit weight; about 1 when the weights are right |
| stream | one (system, observation type) pair, with one noise law |
| system | one source of data with its own local parameters (`task.addSystem`) |
| task | the whole declaration: one `lsqvce.TaskDefinition` |
| whitening | transforming correlated errors into uncorrelated ones before the solve |
