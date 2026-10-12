# `buildDay` — the function you write

## Signature

`fit` passes as many of these arguments as the handle declares (`nargin`); a
`varargin` builder receives all six:

```matlab
D = buildDay(day, obsType)
D = buildDay(day, obsType, common, local, procOpt)
D = buildDay(day, obsType, common, local, procOpt, info)
```

| argument | what it is |
|---|---|
| `day` | one entry of the system's calendar (numeric) |
| `obsType` | one entry of `fieldnames(covariance)` (char) |
| `common` | `data.commonParams` |
| `local` | this system's `localParams` (or `{}`) |
| `procOpt` | this system's `processingOptions` |
| `info` | this day's expansion: every column offset of A and B (below); `[]` when nothing is scoped |

`fit` calls `buildDay` once per (day, observation type). Every pass over the data reads
each day twice (accumulation, then the residual sweep), and an iterative run multiplies
that by its passes. Hence two rules:

1. **Determinism.** The same arguments must give byte-identical output: no random
   numbers, clocks or mutable globals.
2. **Cache expensive reads yourself.** The toolbox keeps `Options.CacheDays` built days
   resident (default one); a raw file read shared by several observation types belongs
   in your own cache, in the closure that binds `buildDay`.

## The `info` argument

```matlab
info.arcId            % 1xA  ids of this day's arcs
info.offsetOf         % @(arcId) -> column offset of that arc's block in B
info.count            % 1xA  columns of each arc block (= number of 'arc' rows)
info.dayOffset        % 1x1  this day's 'day' block offset in B ([] if none)
info.dayCount         % 1x1  its width (0 if none)
info.sysOffset        % 1x1  the 'system' block offset in B
info.sysCount         % 1x1  its width
info.nl               % 1x1  total local count of the system
info.nc               % 1x1  total common count
info.commonCount      % 1x1  width of A's leading 'common' block
info.commonDayOffset  % 1x1  this day's common day-block offset in A ([] if none)
info.commonDayCount   % 1x1  its width (0 if none)
info.arcRows, info.dayRows, info.sysRows, info.commonDayRows
                      %      template rows that exist on this day (with paramColumns)
```

Fill `B(rows, info.offsetOf(id) + (1:k))` and `A(:, info.commonDayOffset + (1:w))`.

## The return value: one cell per arc

Every field is a **1 × nArc cell**; cell a is arc a. A day with one arc returns `{A}`;
a day without data returns empty cells.

| field | each cell | meaning |
|---|---|---|
| `A` | m_a × nc | design matrix of the common unknowns |
| `B` | m_a × nl | design matrix of the local unknowns; may be `sparse` |
| `l` | m_a × 1 | observed minus computed (the reduced observation) |
| `t` | m_a × 1 | epoch of each row in **seconds** (the whitener needs the lags) |
| `Cep` | n_a×1, n_a×d or n_a×d×d | *optional* — the per-epoch accuracy of the arc (below) |
| `nAxis` | scalar | *optional* — the axis count of a multi-axis stream without `Cep` |

Arc lengths may differ freely. An arc whose cells have 0 rows is a **data gap** and is
skipped. Return **raw, unwhitened** matrices: the toolbox whitens every arc from the
covariance spec and reweights it during the variance-component loop. No other field is
allowed (`lsqvce:makeLoader:badDay`).

### `Cep` — the per-epoch accuracy

| shape | meaning | the numbers are |
|---|---|---|
| n_a × 1 | scalar stream with per-epoch accuracy | standard deviations |
| n_a × d | multi-axis stream, axes uncorrelated within an epoch | standard deviations per axis |
| n_a × d × d | multi-axis stream, full per-epoch covariance | covariance matrices |

The axis count d is `size(Cep, 2)`, the same for the whole day; n_a·d must equal the
arc's row count. Rows of a multi-axis stream are **epoch-major** (x₁ y₁ z₁ x₂ …) and
every arc holds whole epochs. `Cep` is relative to the level `sigma` of the covariance
spec: declare `sigma = 1` when `Cep` already carries the absolute accuracy. A
multi-axis day without `Cep` is whitened with D = I and reported once per run
(`lsqvce:makeLoader:noCep`); a `var` law takes no `Cep`.

## The arc schedule — `systems(k).arcs`

| form | meaning |
|---|---|
| omitted | the arcs of a day are the cells `buildDay` returns (ids 1…nArc); days may differ in their count. `'day'` rows live on every calendar day; a calendar day returning no arc in any stream is refused (`lsqvce:makeLoader:emptyDay`). |
| `4` or `[1 2 4]` | static: these arc ids on every day |
| `struct('kbr', 2, 'pos', 4)` | static, per observation type |
| `@(day, obsType) ids` | a resolver returning the ids of a day *without reading data* (ragged schedules, gaps, manoeuvres) |

With a declared schedule cell a is the a-th scheduled id, and the number of cells must
match (`lsqvce:makeLoader:arcCount`). An arc a scheduled day does not carry is returned
as an empty cell. `'arc'` template rows require a declared schedule.

## Checks while reading

| identifier | meaning |
|---|---|
| `lsqvce:makeLoader:buildDay` | `buildDay` threw; its error is attached as the cause |
| `lsqvce:makeLoader:badDay` | a required field is missing or not a 1×nArc cell, or an unknown field is present |
| `lsqvce:makeLoader:dayShape` | the cells disagree in number, or an arc's `A`, `B`, `l`, `t` disagree in rows |
| `lsqvce:accum:arcBadNc` / `arcBadNl` | `A` or `B` has the wrong number of columns |
| `lsqvce:makeLoader:arcCount` | a declared schedule names n arcs and the day returned a different number of cells |
| `lsqvce:makeLoader:emptyDay` | without a schedule, a calendar day returned no arc in any stream |
| `lsqvce:makeLoader:badArc` | an arc id outside the declared schedule, or an arc of a multi-axis day that does not hold whole epochs |
| `lsqvce:makeLoader:badDayCep` | `Cep` is malformed or does not match the arc's rows |
| `lsqvce:makeLoader:axisConflict` | `nAxis` disagrees with `Cep` or does not divide the row count |
| `lsqvce:makeLoader:cepTwice` | the covariance spec carries a `Cep` and the day returns one |
| `lsqvce:covModel:cepOnVar` | a `var` law with a `Cep` |
| `lsqvce:makeLoader:noCep` | *(warning, once per run)* a multi-axis day without `Cep` |

Errors in the declaration itself are raised by `fit` before any data is read and carry
`lsqvce:fit:` identifiers (`unknownField`, `noSystems`, `noCommon`, `missingField`,
`noDays`, `badId`, `duplicateSystem`, `builderKind`, `badArcs`, `badScope`, …).
