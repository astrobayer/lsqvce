# Prepare and Build — the functions you write

Two kinds of function produce the data on demand. Each has one home:

```matlab
sys = task.addSystem(name, Prepare=@prepare, ...);   % P   = prepare(day, ctx)
sys.addStream(name, Build=@build, Noise=noise, ...);  % out = build(day, ctx)
```

| function | runs | receives | returns |
|---|---|---|---|
| `Prepare` (optional, per system) | once per build of a day, before the first stream of that day | the day and a `lsqvce.DayContext` | anything; it reaches every `Build` of the day as `ctx.Prepared` |
| `Build` (required, per stream) | once per build of a day, for its stream | the day and a `lsqvce.DayContext` | the `lsqvce.DayData` of `ctx.newDay()`, filled with `addArc` |

Both take exactly two arguments. Every pass over the data reads each day twice
(accumulation, then the residual sweep), and an iterative run multiplies that by its
passes. Hence two rules:

1. **Determinism.** The same day must give byte-identical output: no random numbers,
   clocks or mutable globals (seed a simulator on the day and the stream).
2. **Share the day's common work through `Prepare`.** The toolbox keeps
   `Options.CacheDays` built days resident (default one). Work that several streams of a
   day need — reading a file, integrating the orbits and their partials — goes into the
   system's `Prepare`, which runs once for all of them each time the day is built.

| `Prepare` rule | |
|---|---|
| how often | once each time the day is built: once per run with `Options.CacheDays = Inf`; with the default `CacheDays = 1` a day is built twice per pass, so it runs as often; with `Options.Parallel` every worker holds its own copy of the loader, so it runs on every worker that builds the day |
| memory | the result of the current day only, whatever `CacheDays` |
| refused | `Settings` with a field `prepared` (`lsqvce:makeLoader:preparedField`, reported by `fit` under its own prefix, the original attached as the cause) |
| when it fails | `lsqvce:makeLoader:prepareDay` names the day; your error is the cause (an `lsqvce:*` error passes unchanged) |

## `lsqvce.DayContext` — what the functions know

| property | type | meaning |
|---|---|---|
| `Day` | double | the day being built |
| `System`, `Stream` | string | the names of the system and the stream (`Stream` is `""` in `Prepare`) |
| `Settings` | struct | the system's `Settings` |
| `Prepared` | any | what `Prepare` returned for this day (`[]` in `Prepare`, or without one) |
| `NumCommon` | double | nc, the columns of A |
| `NumLocal` | double | nl, the columns of the system's B (0 without locals) |
| `ArcIds` | double 1×A | the stream's scheduled arcs of the day; in `Prepare`, the system's; `[]` without a schedule |
| `Common`, `Local` | table | one row per template: `Name`, `Scope`, `OnDay` (the template exists on this day), `Info` (a cell; the template's `Info` row as a struct, `struct()` without one) |

| method | returns | meaning |
|---|---|---|
| `newDay()` | `lsqvce.DayData` | an empty day of this stream |
| `cols(names)` | double 1×n | the columns of named parameters on this day: of `"common"` and `"day"` commons in A, of `"day"` and `"system"` locals in B |
| `cols(names, Arc=id)` | double 1×n | the same for `"arc"` locals of scheduled arc `id` |

One `cols` call names commons or locals, not both.

## `out.addArc(A, l, t, Name=Value)` — the arcs of the day

| argument | type / size | default | meaning |
|---|---|---|---|
| `A` | m × nc | required | design of the common parameters, columns in A's layout (chapter *The task definition*) |
| `l` | m × 1 | required | observed minus computed, raw: never whitened or weighted |
| `t` | m × 1 | required | epoch of each row in **seconds** (the whitener needs the lags) |
| `Local` | n × 2 cell `{names, P; ...}` | none | local partials by name, P m × numel(names); each name's scope decides its column (this arc's, this day's or the system's block) |
| `B` | m × nl, may be sparse | none | the local design as one matrix, every column from `ctx.cols`. Give `Local` or `B`, never both; neither = an arc that sees no local parameter |
| `Sigma` | e×1, e×d or e×d×d | none | the per-epoch accuracy of the arc (below) |
| `Excluded` | m × 1 logical | none | rows of **weight 0**: zero in A, the local partials and l, kept for their index (residuals, diagnostics) but not observations — out of the degrees of freedom, `rejected` in the verdict. Excluded rows are deleted from the law's covariance exactly, so they do not condition their neighbours |
| `Arc` | a scheduled arc id | the next | places this arc in the schedule |

Arc lengths may differ freely. An arc with 0 rows is a **data gap** and is skipped.

### `Sigma` — the per-epoch accuracy

| shape | meaning | the numbers are |
|---|---|---|
| e × 1 | scalar stream with per-epoch accuracy | standard deviations |
| e × d | multi-axis stream, axes uncorrelated within an epoch | standard deviations per axis |
| e × d × d | multi-axis stream, full per-epoch covariance (d ≥ 2) | covariance matrices |

A scalar stream's `Sigma` is always **standard deviations**. MATLAB stores an e × 1 × 1
array as e × 1, so a per-epoch *variance* array would be read as sigmas: pass the square
root. A full `Sigma` is composed through its Cholesky factor in the given axis order; when
it varies along the arc, the composed covariance therefore depends on the frame and the
axis order (XYZ against NEU: about 2 % of C). A constant or diagonal `Sigma` does not.

The axis count d is `size(Sigma, 2)`, the same for the whole day, and must agree with the
stream's `Axes` when one above 1 is declared; e·d must equal the arc's row count. Rows of a multi-axis stream are
**epoch-major** (x₁ y₁ z₁ x₂ …) and every arc holds whole epochs. `Sigma` is relative to
the law's `sigma`: declare `lsqvce.Noise.white(1)` when `Sigma` already carries the
absolute accuracy. A multi-axis day without `Sigma` is whitened with D = I and reported
once per run (`lsqvce:makeLoader:noCep`); a `var` law takes no `Sigma`.

## The arc schedule — `Arcs`

| form | meaning |
|---|---|
| omitted | the arcs of a day are those `Build` adds (ids 1…n); days may differ in their count. `"day"` locals live on every calendar day; a calendar day with no arc in any stream is refused (`lsqvce:makeLoader:emptyDay`) |
| `4` or `[1 2 4]` | static: arcs 1…4, or these ids, on every day |
| `@(day, stream) ids` | a resolver returning the ids of a day *without reading data* (ragged schedules, gaps, manoeuvres) |

With a schedule, `ctx.ArcIds` lists the day's arcs; `addArc` fills them in order, `Arc=id`
places an arc out of order, and an arc the day schedules but `Build` never adds is a data
gap. `"arc"` locals require a schedule. A stream may carry its own `Arcs` on `addStream`;
then every stream of the system needs one (from the system or its own).

## Checks while reading

| identifier | meaning |
|---|---|
| `lsqvce:makeLoader:buildDay` | a `Build` threw (the message names the day and the stream); your error is attached as the cause (an `lsqvce:*` error passes unchanged) |
| `lsqvce:DayContext:badReturn` | `Build` returned something other than its `ctx.newDay()` |
| `lsqvce:DayContext:unknownName` | a name in `Local` or `cols` that is not a parameter of the system (or task), or not instantiated on this day |
| `lsqvce:DayContext:arcRequired` | an `"arc"` local without `Arc=id`, or with an id the day does not schedule |
| `lsqvce:DayContext:mixedTables` | one `cols` call names commons and locals |
| `lsqvce:DayContext:localAndB` | one arc gives both `Local` and `B` |
| `lsqvce:DayContext:badLocal` | `Local` is not `{names, P}` rows with P m × numel(names), or names a parameter twice in one arc |
| `lsqvce:DayContext:badArc` | `Arc` without a schedule, an unscheduled id, an arc added twice, more arcs than scheduled |
| `lsqvce:makeLoader:dayShape` | an arc's `A`, `B`, `l`, `t` disagree in rows |
| `lsqvce:accum:arcBadNc` / `arcBadNl` | `A` or `B` has the wrong number of columns |
| `lsqvce:makeLoader:arcCount` | a declared schedule names n arcs and the day carries a different number |
| `lsqvce:makeLoader:emptyDay` | without a schedule, a calendar day returned no arc in any stream |
| `lsqvce:makeLoader:badArc` | an arc of a multi-axis day that does not hold whole epochs |
| `lsqvce:makeLoader:badDayCep` | `Sigma` is malformed or does not match the arc's rows |
| `lsqvce:makeLoader:axisConflict` | `Axes` disagrees with `Sigma` or does not divide the row count |
| `lsqvce:covModel:cepOnVar` | a `var` law with a `Sigma` |
| `lsqvce:makeLoader:noCep` | *(warning, once per run)* a multi-axis day without `Sigma` |
