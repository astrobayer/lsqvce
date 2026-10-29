# Troubleshooting

Every error the toolbox raises has an identifier (`ME.identifier`) and a message that
names what is wrong: the field, the system, the day, the stream. Read the message first;
this chapter adds the usual cause and the usual fix. The full list of identifiers is in
the API reference (chapter *Error and warning identifiers*).

```matlab
try
    R = lsqvce.fit(data, opts);
catch ME
    disp(ME.identifier)
    disp(ME.message)
    if ~isempty(ME.cause), disp(ME.cause{1}.message), end   % YOUR error, when your function threw
end
```

## Before any data is read: the declaration

`fit` checks the whole declaration before it calls your functions. These errors cost
nothing to fix and to re-run. A few of them are raised while `fit` builds the loaders;
`fit` reports those under its own prefix: the table lists them by the `makeLoader`
identifier, and the `ME.identifier` you see is the prefix `lsqvce:fit:` followed by its
last part (`prepareUnused`, say), with the original attached as the cause.

| identifier | usual cause | fix |
|---|---|---|
| `lsqvce:fit:unknownField` | a typo in a field name, e.g. `covaraince`; the message lists the allowed names | correct the name. Unknown fields are errors on purpose: a typo is never ignored |
| `lsqvce:fit:missingField` | a system without `buildDay` or without `covariance` | add it |
| `lsqvce:fit:noCommon` | `data.commonParams` missing, or not a cell table | one row per parameter, scope token in column 1 |
| `lsqvce:fit:badScope` | a scope token that does not exist (`'arcs'` instead of `'arc'`) | `commonParams`: `'common'` or `'day'`; `localParams`: `'arc'`, `'day'` or `'system'` |
| `lsqvce:fit:noDays`, `badDays`, `duplicateDay` | no day list, a non-integer day, or a day listed twice | `data.days = 1:30` (or a system's own `days`) |
| `lsqvce:fit:builderKind` | `buildDay` or `prepareDay` is a name (`'myBuildDay'`) instead of a handle | `@myBuildDay` |
| `lsqvce:makeLoader:localPlan` | `'arc'` parameters without an arc schedule | declare `s.arcs` (a count, an id vector or `@(day, obsType)`) |
| `lsqvce:makeLoader:prepareUnused` | a `prepareDay` with a `buildDay` that declares fewer than five arguments, so it could never see the result | declare `buildDay(day, obsType, common, local, procOpt[, info])` |
| `lsqvce:makeLoader:preparedField` | your `processingOptions` already has a field `prepared` | rename your field; the name is taken by `prepareDay`'s result |
| `lsqvce:fit:duplicateSystem` | two systems with the same `id` | give each system its own `id` |
| `lsqvce:fit:badPrior` | `L` not `p x nc`, or `x0` of the wrong length | write the prior in the columns `lsqvce.layout(data)` announces |
| `lsqvce:fit:modelKind` | `Estimate = "noise"` on the matrices or batched kind | the noise estimation needs streams: use the streaming kind |

## While reading: your `buildDay` and `prepareDay`

These come from inside the run, on the first day where something is wrong. The message
names the day and the stream.

| identifier | usual cause | fix |
|---|---|---|
| `lsqvce:makeLoader:buildDay` | your `buildDay` threw an error | read `ME.cause{1}`: it is your own error, unchanged. Call your `buildDay` for that day by hand |
| `lsqvce:makeLoader:prepareDay` | your `prepareDay` threw an error | as above, with `prepareDay` |
| `lsqvce:makeLoader:dayShape` | `A`, `B`, `l`, `t` are not cells of one per arc, or the arcs disagree in their numbers of rows | return `struct('A', {cellOfArcs}, ...)`: note the braces around the cell |
| `lsqvce:accum:arcBadNc`, `arcBadNl` | `A` does not have `info.nc` columns, or `B` not `info.nl` | build `B = zeros(m, info.nl)` and place each block by `info.offsetOf`, `info.dayOffset` |
| `lsqvce:covModel:denseArcTooLong` | an arc longer than 20 000 epochs under a dense noise law or under `Estimate = "noise"` | cut the arcs shorter, or see *Large problems* |

A frequent mistake that raises nothing: a `buildDay` that is **not deterministic**, for
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
| the verdict is `INVALID` | the integrity check failed; nothing else was evaluated | most often a non-deterministic `buildDay` |
| `R.summary()` or `R.plot("psd")` raises `lsqvce:Results:noHistory` | the run was a single solve | these need `Estimate = "noise"` (or a robust run) |
| `R.residuals()` raises `lsqvce:Results:noAccumulator` | the matrices or batched kind | read `R.diag.v` instead |

## Asking for help

When you report a problem, include `R.meta` (toolbox version, MATLAB release, the options
as resolved, the seeds) or, if `fit` did not return, the full output of
`getReport(ME, 'extended')`. Never attach measurement data you are not allowed to share;
a case with synthetic data of the same shape is better.
