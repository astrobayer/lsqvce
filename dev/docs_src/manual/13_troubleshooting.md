# Troubleshooting

Every error the toolbox raises has an identifier (`ME.identifier`) and a message that
names what is wrong: the field, the system, the day, the stream. Read the message first;
this chapter adds the usual cause and the usual fix. The full list of identifiers is in
the API reference (chapter *Error and warning identifiers*).

```matlab
try
    R = lsqvce.fit(task, opts);
catch ME
    disp(ME.identifier)
    disp(ME.message)
    if ~isempty(ME.cause), disp(ME.cause{1}.message), end   % YOUR error, when your function threw
end
```

## While you declare: the `add` calls

A `lsqvce.TaskDefinition` checks every call when you make it, so most mistakes are
refused at the line that makes them, before `fit` runs. What only the whole task can know
is checked by `task.declaration()`, which `fit` calls first. A few are raised while `fit`
builds the loaders; `fit` reports those under its own prefix: the table lists them by the
`makeLoader` identifier, and the `ME.identifier` you see is `lsqvce:fit:` followed by its
last part (`localPlan`, say), with the original attached as the cause.

| identifier | usual cause | fix |
|---|---|---|
| `MATLAB:validators:*` | a setting of the wrong kind: `Scope="arcs"`, `Estimate="all"`, a negative `Sigma` | the message names the setting and its allowed values |
| `lsqvce:TaskDefinition:duplicateName` | a name used twice: two commons, two locals of a system, a local named like a common, two systems, streams or priors | every name is unique where it lives |
| `lsqvce:TaskDefinition:badName` | a stream name that is not a MATLAB identifier (`"pos C"`) | letters, digits and underscores, starting with a letter |
| `lsqvce:TaskDefinition:badValid` | `Valid` not `[start end]` with start ≤ end | fractional MJD, or `[]` for always |
| `lsqvce:TaskDefinition:infoRows` | an `Info` table without one row per name | one row per name |
| `lsqvce:TaskDefinition:badHandle` | `Prepare` or `Build` is not a function handle, or takes fewer than two arguments | `@myBuild`, with `function out = myBuild(day, ctx)` |
| `lsqvce:TaskDefinition:missingArgument` | `addStream` without `Build` or `Noise`; a prior without `Sigma`; `fromBatches` without `NumBatches` | add it |
| `lsqvce:TaskDefinition:priorTarget` | a prior with both `On` and `L`, or neither | one of the two |
| `lsqvce:TaskDefinition:wrongKind` | `addCommon`, `addSystem` or a prior with `On` on a task made by `fromMatrices` / `fromBatches` | in-memory tasks take priors by `L` only |
| `lsqvce:TaskDefinition:badArcs` | some streams of a system have an arc schedule and others none | give the system one `Arcs`, or every stream its own |
| `lsqvce:TaskDefinition:incomplete` | no `addCommon`, no `addSystem`, a system without a stream, or no calendar | complete the task |
| `lsqvce:makeLoader:localPlan` | `"arc"` locals without an arc schedule | `addSystem(..., Arcs=...)`: a count, an id vector or `@(day, stream)` |
| `lsqvce:makeLoader:preparedField` | the system's `Settings` has a field `prepared` | rename your field; the name carries `Prepare`'s result |
| `lsqvce:fit:badPrior` | a pattern `On` that matches no name, or `L` of the wrong width | check the names; write `L` in the columns `task.layout()` announces |
| `lsqvce:fit:modelKind` | `Estimate = "noise"` on a task made by `fromMatrices` / `fromBatches` | the noise estimation needs streams: use the streaming kind |

## While reading: your `Prepare` and `Build`

These come from inside the run, on the first day where something is wrong. The message
names the day and the stream.

| identifier | usual cause | fix |
|---|---|---|
| `lsqvce:makeLoader:buildDay` | your `Build` threw an error (the message says `buildDay`: the stream's generated reader) | read `ME.cause{1}`: it is your own error, unchanged |
| `lsqvce:makeLoader:prepareDay` | your `Prepare` threw an error | as above |
| `lsqvce:DayContext:unknownName` | a name in `Local` or `ctx.cols` that is not a parameter of the system, or not on this day (`Valid`) | check `ctx.Local` (its `OnDay` column) |
| `lsqvce:DayContext:arcRequired` | an `"arc"` local in `ctx.cols` without `Arc=id`, or with an id the day does not schedule | `ctx.cols(names, Arc=ctx.ArcIds(a))` |
| `lsqvce:DayContext:mixedTables` | one `ctx.cols` call names commons and locals | one call for each: commons are columns of `A`, locals of `B` |
| `lsqvce:DayContext:localAndB` | one arc gives both `Local` and `B` | one of the two |
| `lsqvce:DayContext:badLocal` | `Local` not `{names, partials}` rows, partials not m × numel(names), a name twice | as the message says |
| `lsqvce:DayContext:badArc` | `Arc=` without a schedule, an id not scheduled, an arc added twice | follow `ctx.ArcIds` |
| `lsqvce:DayContext:badReturn` | `Build` returned something other than `ctx.newDay()` | `out = ctx.newDay(); out.addArc(...);` and return `out` |
| `lsqvce:makeLoader:dayShape` | an arc whose `A`, `l`, `t` disagree in their numbers of rows | one row count per arc |
| `lsqvce:accum:arcBadNc`, `arcBadNl` | `A` does not have `ctx.NumCommon` columns, or `B` not `ctx.NumLocal` | build `A` with one column per `addCommon` name; give locals by name |
| `lsqvce:covModel:denseArcTooLong` | an arc longer than 20 000 epochs under a dense noise law or under `Estimate = "noise"` | cut the arcs shorter, or see *Large problems* |

A frequent mistake that raises nothing: a `Build` that is **not deterministic**, for
example one drawing fresh random noise on every call. The toolbox reads each day at least
twice per pass and then sees two different data sets; the report's integrity check
(Level 0) is the place where this shows. Seed any random numbers on the day and the
stream, as the tutorial does.

## During the solve

| message | meaning | what to do |
|---|---|---|
| `lsqvce:qr:singularR` / `lsqvce:qr:rankDeficient` | the common parameters are not all determined by the data: a free network without a datum, an offset that a per-arc bias absorbs | add a datum or a regularisation as a prior (*Priors and regularisation*), or remove the parameter the data cannot see. The toolbox never guesses a datum |
| `lsqvce:adjustVCEjoint:notConverged` (warning) | the noise loop used up `MaxPasses` before it settled; `R.stopReason` is `"maxPasses"` | the last pass is returned in full. Look at `R.plot("convergence")`; raise `MaxPasses`, or see *When the loop did not settle* in *Reading the results and the report* |

## The run finishes, but the numbers look wrong

| symptom | likely cause | where to read |
|---|---|---|
| σ₀ far from 1 in a run with `Estimate = "none"` | the declared `sigma` (or `p`, or `Cep`) is wrong in scale; often an information matrix used where a covariance was expected | *Choosing the stochastic model*, "Run one pass with `Estimate = "none"` first" |
| σ₀ fine, but the formal errors look too small | the noise is correlated and was declared white | run `Estimate = "noise"`; the tutorial shows the effect |
| `R.sigma0` is `NaN` | the quadratic form vᵀPv cancelled in double precision, and the toolbox reports that instead of clamping it | a `NaN` is a result, not a bug; the theory (section 1.3) explains when it happens |
| the verdict is `SUSPECT` | at least one finding at `ERROR` level | work through the `[ERROR]` lines of the report, top level first |
| the verdict is `INVALID` | the integrity check failed; nothing else was evaluated | most often a non-deterministic `Build` |
| `R.summary()` or `R.plot("psd")` raises `lsqvce:Results:noHistory` | the run was a single solve | these need `Estimate = "noise"` (or a robust run) |
| `R.residuals()` raises `lsqvce:Results:noAccumulator` | the matrices or batched kind | read `R.diag.v` instead |

## Asking for help

When you report a problem, include `R.meta` (toolbox version, MATLAB release, the options
as resolved, the seeds) or, if `fit` did not return, the full output of
`getReport(ME, 'extended')`. Never attach measurement data you are not allowed to share;
a case with synthetic data of the same shape is better.
