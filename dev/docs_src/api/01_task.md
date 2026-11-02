# The task definition

Everything the toolbox knows about a problem enters through one object, a
`lsqvce.TaskDefinition`, handed to `lsqvce.fit(task, opts)`. It is filled with flat calls
that take name–value arguments; every call checks its arguments when it is made and names
the one at fault. Nothing is read before `fit` runs. `help lsqvce.TaskDefinition` and
`doc lsqvce.TaskDefinition` show the same contract.

| how it is made | for |
|---|---|
| `lsqvce.TaskDefinition(Days=days)` and its `add` calls | **streaming**: arc-wise data read day by day through your `Build` functions, local unknowns eliminated on the fly, estimated stochastic models |
| `lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p)` | **matrices**: a design matrix that fits in memory, no local unknowns, diagonal weights |
| `lsqvce.TaskDefinition.fromBatches(@batchFcn, NumBatches=K, NumUnknowns=u)` or `X0=x0` | **batched**: a design too large for one matrix, or a nonlinear model |

## `lsqvce.TaskDefinition(Name=Value)`

| name | type / size | default | meaning |
|---|---|---|---|
| `Days` | double 1×D, integers, each once | `[]` (then every system sets its own `Days`) | the default calendar (e.g. MJD); the loader cache is keyed on the days |
| `Name` | string | `""` | a label for the run, shown by `disp(task)` |

| method | returns | meaning |
|---|---|---|
| `addCommon(names, Name=Value)` | — | common parameters |
| `sys = addSystem(name, Name=Value)` | `lsqvce.SystemDefinition` | an observing system |
| `addPrior(name, Name=Value)` | — | a prior group on the common parameters (chapter *Priors*) |
| `lay = layout()` | struct | the columns of A and B, with names (below) |
| `data = declaration()` | struct | the compiled engine-tier struct (chapter *The compiled declaration*) |
| `disp(task)` | — | parameter templates per scope, systems, streams and priors |

## `addCommon(names, Name=Value)` and `sys.addLocal(names, Name=Value)`

One call adds one **template** per name, all of one scope. A template becomes one
parameter per instance of its scope.

| name | type / size | default | meaning |
|---|---|---|---|
| `names` | string array | required | unique within the task (commons) or within the system (locals); a local may not share a name with a common |
| `Scope` | `"common"`, `"day"` (commons); `"arc"`, `"day"`, `"system"` (locals) | `"common"` / `"arc"` | over which span the parameter is constant |
| `Valid` | double 1×2 `[start end]`, fractional MJD, or `[]` | `[]` (always) | a `"day"` or `"arc"` template exists only on days its interval overlaps (start < d+1 and end > d); `"common"` and `"system"` templates are always present |
| `Info` | table, one row per name | none | your own columns; `Prepare` and `Build` read each template's row as a struct, `ctx.Common.Info{k}` / `ctx.Local.Info{k}` |

| scope | instantiated as | typical use |
|---|---|---|
| `"common"` | one column of A, always | spherical-harmonic coefficients |
| `"day"` (common) | one column of A per calendar day that carries at least one arc in at least one system | a daily common parameter |
| `"arc"` | one column of the system's B per (day, arc) | initial state, per-arc empirical parameters; needs an arc schedule |
| `"day"` (local) | one column per day of the system that carries ≥ 1 arc — every calendar day without a schedule | accelerometer scale / bias / drift |
| `"system"` | one column for the whole adjustment | instrument constants |

**Column layout of A.** The `"common"` templates first, in `addCommon` order; then one day
block per calendar day (the sorted union of every system's days), holding that day's
`"day"` commons. **Column layout of B (per system).** Per day, in calendar order: the arc
blocks in sorted arc-id order (each the `"arc"` templates in `addLocal` order), then the
day block; after all days one trailing `"system"` block. Your code never computes these
positions: `Build` gives partials by name, or takes columns from `ctx.cols`.

## `sys = addSystem(name, Name=Value)`

| name | type / size | default | meaning |
|---|---|---|---|
| `name` | string, unique in the task | required | the system's name in reports |
| `Prepare` | function handle `@(day, ctx)` | none | the work every stream of a day shares, run once per build of the day (chapter *Prepare and Build*) |
| `Settings` | 1×1 struct | `struct()` | anything your functions need (paths, force models); `ctx.Settings`. The field name `prepared` is reserved |
| `Days` | double 1×D, integers, each once | the task's `Days` | this system's calendar |
| `Arcs` | count, id vector, or `@(day, stream)` | none | the arc schedule (chapter *Prepare and Build*); required by `"arc"` locals |
| `Id` | finite double scalar or single char | its position | an identifier, unique across systems; give every system an id or none |

A system is complete with one `addStream`. Its methods:

| method | meaning |
|---|---|
| `addLocal(names, Name=Value)` | local parameters (above) |
| `addStream(name, Name=Value)` | an observation stream (below) |
| `addPrior(name, Name=Value)` | a prior group on this system's locals (chapter *Priors*) |

## `sys.addStream(name, Name=Value)`

| name | type / size | default | meaning |
|---|---|---|---|
| `name` | string, a MATLAB identifier, unique in the system | required | the observation type; streams are processed in the order added |
| `Build` | function handle `@(day, ctx)` | required | returns the stream's arcs of one day (chapter *Prepare and Build*) |
| `Noise` | 1×1 `lsqvce.Noise` | required | the temporal law and its level (chapter *Noise laws*) |
| `Estimate` | `"law"`, `"scale"`, `"none"` | `"law"` (`"scale"` for a per-row `sigma`) | what `Estimate = "noise"` may change: shape and level, the level only, nothing |
| `ArcFactors` | logical | `Options.ArcFactors` | per-arc variance factors for this stream |
| `Axes` | integer ≥ 1 | 1 | axes per epoch; rows are epoch-major |
| `Arcs` | as `addSystem` | the system's | this stream's own schedule; a system gives every stream a schedule, or none |

## `fromMatrices` and `fromBatches`

| call | arguments | meaning |
|---|---|---|
| `fromMatrices(A, l, Weights=p)` | A m×u, l m×1, p m×1 | weights are the diagonal of P (`1./sigma.^2`); absent = unit weights |
| `fromBatches(@batchFcn, NumBatches=K, NumUnknowns=u)` | `[A, l, p] = batchFcn(k)`, k = 1…K | one batch of design rows, observations and weights at a time |
| `fromBatches(@batchFcn, NumBatches=K, X0=x0)` | `[J, w, p] = batchFcn(k, x)` | the nonlinear form: the Jacobian at x and `w = l − f(x)`; damped Gauss–Newton; u = `numel(x0)` |

Both take priors by `L` only (they have no parameter names) and run with
`Estimate = "none"` (variance components need streams, `lsqvce:fit:modelKind`); `Robust`
selects M-estimation, except on the nonlinear form (`lsqvce:fit:notYet`). Their
per-observation diagnostics are in the whitened frame, `R.diag.v = sqrt(p).*v`.

## `lay = task.layout()` — the columns announced

```matlab
lay = task.layout();                           % the same as lsqvce.layout(task)
lay.nc, lay.common(j)                          % .name, .row, .day of column j of A
lay.systems(k).name, lay.systems(k).nl
lay.systems(k).local(c)                        % .name, .row, .day, .arc of column c of B
c = lay.systems(k).col(row, day, arc);         % lookup by template row, 0 = not instantiated
c = lay.colCommon(row, day);
```

This is the column map a prior operator `L` is written in, and the one `R.toTable(task)`
uses. `row` is the position of a name in `addCommon` / `addLocal` order.

## Errors while declaring

| identifier | raised by | meaning |
|---|---|---|
| `lsqvce:TaskDefinition:duplicateName` | `addCommon`, `addLocal`, `addSystem`, `addStream`, `addPrior`, `declaration` | a name used twice where it must be unique (also a local named like a common) |
| `lsqvce:TaskDefinition:badName` | `addStream` | a stream name that is not a MATLAB identifier |
| `lsqvce:TaskDefinition:badValid` | `addCommon`, `addLocal` | `Valid` is not `[]` or `[start end]` with start ≤ end |
| `lsqvce:TaskDefinition:infoRows` | `addCommon`, `addLocal` | `Info` has not one row per name |
| `lsqvce:TaskDefinition:badHandle` | `addSystem`, `addStream` | `Prepare` or `Build` is not a function handle taking `(day, ctx)` |
| `lsqvce:TaskDefinition:missingArgument` | `addStream`, `addPrior`, `fromBatches` | `Build`, `Noise`, `Sigma` or `NumBatches` missing |
| `lsqvce:TaskDefinition:priorTarget` | `addPrior` | both `On` and `L`, or neither |
| `lsqvce:TaskDefinition:wrongKind` | `addCommon`, `addSystem`, `addPrior`, `layout` | a streaming call on an in-memory task |
| `lsqvce:TaskDefinition:badArcs` | `declaration` | some streams of a system have a schedule and others none |
| `lsqvce:TaskDefinition:incomplete` | `declaration` | no common parameter, no system, a system without a stream, or no calendar |
| `MATLAB:validators:*` | every call | a value of the wrong kind (a `Scope` that does not exist, a negative `Sigma`, …) |

The checks the engines make on the compiled struct (chapter *The compiled declaration*)
still run when `fit` starts; they carry `lsqvce:fit:` identifiers.
