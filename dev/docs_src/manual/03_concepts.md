# Key concepts

This chapter introduces the vocabulary the rest of the Manual uses. It has no code to run.
Read it once before the tutorial. The words in **bold** are the names you will meet as
field names, options and chapter titles.

## The adjustment in one paragraph

You have observations `l` that depend on unknown parameters `x` through a linear (or
linearised) model, `l + v = A x`, where `v` are the errors. Least squares picks the `x`
that makes the weighted errors smallest. The weights are the inverse of the noise
covariance of the observations. Usually you do not know that covariance well. *Variance
component estimation* (VCE) estimates it from the residuals, re-weights, and solves
again until the weights settle. lsqvce does all of that in one call, `lsqvce.fit`.

## Three kinds of problem

`fit` recognises the kind of your problem from the fields of `data`:

| kind | fields you give | use it when |
|---|---|---|
| **matrices** | `A`, `l`, `p` | the whole design matrix fits in memory |
| **batched** | `batchFcn`, `numBatches`, `numUnknowns` (or `x0`) | it does not fit, or the model is nonlinear; your function returns one block of rows at a time |
| **streaming** | `commonParams`, `days`, `systems` | the data come per day and per arc, from several instruments, with nuisance parameters and correlated noise |

The first two are ordinary least squares with diagnostics. The streaming kind is what the
toolbox is built for, and the rest of this chapter describes it.

## Parameters: common and local

Ask of every unknown: *over what span is it constant, and who sees it?*

- **Common parameters** (`data.commonParams`) are seen by all data, for example the
  gravity-field coefficients or station coordinates. They are what you want to know.
- **Local parameters** (`localParams` of a system) belong to one **arc**, one **day** or
  one **system**, for example an initial orbit state, an accelerometer bias, or a
  per-arc offset. They are nuisance parameters. The toolbox eliminates them arc by arc
  while it reads the data, so their number does not cost memory.

You declare both as tables with one row per *kind* of parameter. Column 1 is the scope
(`'common'`, `'day'`, `'arc'` or `'system'`); the other columns are yours (a name, a
degree and order, whatever helps you). The toolbox instantiates each row once per scope
and counts the columns itself: you never compute a column index by hand.

## Data: systems, days, arcs and streams

```
data
├── commonParams          what everybody estimates
├── days                  the calendar, e.g. 1:30
└── systems(k)            one per satellite, receiver, or instrument group
     ├── localParams      its nuisance parameters
     ├── covariance.<type>  one noise law per observation type  ──►  the streams
     ├── prepareDay       (optional) the work all its streams of a day share
     └── buildDay         the rows of one day and one observation type
```

- A **system** is one source of data with its own nuisance parameters, for example one
  satellite, or a satellite pair whose inter-satellite ranging depends on both orbits.
- A **day** is the unit of reading. The toolbox asks for the data of one day at a time.
- An **observation type** is a kind of measurement within a system, for example `'pos'`
  and `'kbr'`. The observation types are the field names of the system's `covariance`
  struct. Each (system, observation type) pair is one **stream**.
- An **arc** is one continuous piece of a stream on a day. The noise is correlated
  within an arc and independent between arcs. `buildDay` returns one cell per arc.

## The two functions you write

You do not hand the toolbox your data. You hand it functions that produce the data on
demand. This is what keeps memory small, but it has one consequence: **the toolbox calls
your functions more than once for the same day**, so they must return the same result
every time.

**`buildDay`** (required) returns the design rows and observations of one day for one
observation type:

```matlab
D = buildDay(day, obsType, common, local, procOpt, info)
%   D.A{a}  design of the common parameters, arc a      (m_a x nc)
%   D.B{a}  design of the local parameters,  arc a      (m_a x nl)
%   D.l{a}  observed minus computed                     (m_a x 1)
%   D.t{a}  epochs in seconds                           (m_a x 1)
```

**`prepareDay`** (optional) does the work that every stream of a day shares, *once* per
day, before the first `buildDay` of that day. Its result reaches every `buildDay` of the
day as `procOpt.prepared`:

```matlab
P = prepareDay(day, common, local, procOpt, info)
```

Without `prepareDay`, a system with two observation types computes anything they share
twice per day: an orbit integration, a large file read, a coordinate rotation. With it,
the shared work happens once per build of the day. The tutorial uses it; chapter
*Declaring your problem* shows when it pays.

Both functions take *as many arguments as they declare*. `@(day, obsType)` is a valid
`buildDay` when you need nothing else. `info`, the 6th argument, tells you where each
block of local parameters sits in `B`.

## What `fit` does with them

For a run with `Estimate = "noise"`, `fit` repeats one **pass** until the noise settles:

```
for every pass
    for every system and every day
        P = prepareDay(day, ...)                  once per day (if declared)
        for every observation type
            D = buildDay(day, obsType, ...)       the day's rows, with procOpt.prepared = P
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

Each stream has a **noise law** (`covariance.<type>`): how its errors correlate from
epoch to epoch (`type`), and how large they are (`sigma`). On top of that, every arc can
have its own **arc factor**, which catches one bad pass without spoiling the whole stream.
With `Estimate = "noise"` the toolbox estimates all of them. You start from an honest
guess, usually white noise (`type = 'diag'`) with a realistic `sigma`.

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
| arc | a continuous piece of one stream on one day; one cell of `buildDay`'s output |
| arc factor | a variance factor of one arc, estimated on top of the stream's level |
| common parameter | an unknown seen by all data (`data.commonParams`) |
| `Cep` | per-epoch accuracy delivered with the data (`D.Cep`), relative to the law's `sigma` |
| day | the unit in which the toolbox reads data |
| local parameter | a nuisance unknown of one arc, day or system, eliminated on the fly |
| observation type | a kind of measurement within a system; a field name of `covariance` |
| pass | one solve plus one update of the noise model |
| pool | streams that share one estimated spectral shape (`Options.Pool`) |
| prior | pseudo-observations on the parameters (`data.priors`): datum, regularisation |
| redundancy | how much of an observation is checked by the others (0 to 1) |
| σ₀ | a-posteriori standard deviation of unit weight; about 1 when the weights are right |
| stream | one (system, observation type) pair, with one noise law |
| system | one source of data with its own local parameters (`data.systems(k)`) |
| whitening | transforming correlated errors into uncorrelated ones before the solve |
