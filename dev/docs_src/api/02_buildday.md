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
| `procOpt` | this system's `processingOptions`; with a `prepareDay`, also the field `prepared` (its result — the name is then reserved) |
| `info` | this day's expansion: every column offset of A and B (below); `[]` when nothing is scoped |

`fit` calls `buildDay` once per (day, observation type). Every pass over the data reads
each day twice (accumulation, then the residual sweep), and an iterative run multiplies
that by its passes. Hence two rules:

1. **Determinism.** The same arguments must give byte-identical output: no random
   numbers, clocks or mutable globals.
2. **Share the day's common work through `prepareDay`.** The toolbox keeps
   `Options.CacheDays` built days resident (default one). Work that several observation
   types of a day need — reading a file, integrating the orbits and their partials — goes
   into the system's `prepareDay`, which runs once for all of them each time the day is
   built (next section).

## Where the handle goes — one per system

`buildDay` is not an option of `fit`; it is part of the data. Every entry of
`data.systems` carries its own handle in the field `buildDay`, next to its own
`covariance`, `localParams`, `days`, `arcs` and `processingOptions`, and `fit` calls each
system's handle on that system's calendar only. `lsqvce.Options` says how one run
proceeds (`OnPass` among it); the `data` struct says what is read and by whom.

| your situation | what to do |
|---|---|
| several instruments or satellites in **one** adjustment | one entry of `data.systems` each, each with its own handle; they share only the common parameters |
| **one reader**, several setups (data directory, background models, editing thresholds) | pass the setup in as data: plain values in `processingOptions` (the 5th argument), or a closure that binds a configuration object |
| data sets that do **not** share common parameters | one `data` struct each, one `fit` call each |

```matlab
% two instruments, two readers, one adjustment
data.systems(1).buildDay = @buildDayKbr;
data.systems(1).covariance.kbr = struct('type', 'diag', 'sigma', 1e-7);
data.systems(2).buildDay = @buildDayPos;
data.systems(2).covariance.pos = struct('type', 'diag', 'sigma', 0.02);

% one reader, the setup as data: processingOptions is the 5th argument ...
s.buildDay          = @myBuildDay;     % D = myBuildDay(day, obsType, common, local, procOpt)
s.processingOptions = struct('dataDir', '/data/setupA', 'gravityModel', 'GOCO06s');

% ... or bound by a closure (the examples do this)
cfg = mySetup('setupB');
s.buildDay = @(day, obsType, common, local, procOpt, info) ...
                 myBuildDay(cfg, day, obsType, common, local, procOpt, info);
```

A closure declares its own argument count, and that count is what `fit` passes
(`nargin`, above): declare `info` when your rows need the column offsets of local or
`'day'` parameters. What the handle captures must not change during the run — it is
called again on every pass (rule 1). It is also the place for a cache of your own that
outlives a pass — the raw files of a day, say, when they are smaller than the built rows;
what the streams of one day share goes into `prepareDay` (rule 2).
`buildDay` that is not a function handle is refused before any data is read
(`lsqvce:fit:builderKind`).

## What all streams of a day share — `prepareDay`

`buildDay` is called once per (day, observation type). When the streams of a day rest on
one expensive computation — the reference orbits and the variational partials of both
satellites, used by the kinematic orbits **and** the range-rates — give the system a
`prepareDay` handle. `fit` calls it **once per build of a day**, before the first stream of
that day is built, and hands its result to every `buildDay` of that day as
`procOpt.prepared`:

```matlab
s.prepareDay = @prepareGrace;       % P = prepareGrace(day, common, local, procOpt, info)
s.buildDay   = @buildDayGrace;      % D = buildDayGrace(day, obsType, common, local, procOpt, info)

function P = prepareGrace(day, ~, ~, procOpt, info)
    P.C = integrateWithPartials('C', day, procOpt, info);   % r, v and d[r;v]/dp, whole day
    P.D = integrateWithPartials('D', day, procOpt, info);
end

function D = buildDayGrace(day, obsType, ~, ~, procOpt, info)
    P = procOpt.prepared;               % what prepareGrace returned for this day
    switch obsType
        case 'posC', D = kinematicRows(P.C, day, info);
        case 'posD', D = kinematicRows(P.D, day, info);
        case 'kbr',  D = rangeRateRows(P.C, P.D, day, info);
    end
end
```

| rule | |
|---|---|
| arguments | as many of `(day, common, local, procOpt, info)` as the handle declares, none to all five; `info` is the same day slice `buildDay` gets, `procOpt` the plain `processingOptions` |
| how often | once each time the day is built, for all its streams: once per run with `Options.CacheDays = Inf`; with the default `CacheDays = 1` a day is built twice per pass (accumulation, residual sweep), so the hook runs as often; with `Options.Parallel` every worker holds its own copy of the loader, so it runs on every worker that builds the day |
| memory | the result of the current day only, whatever `CacheDays`: the toolbox keeps the built streams, not the result |
| determinism | like `buildDay`: it runs again whenever its day is built again |
| refused | a handle of another class (`lsqvce:fit:builderKind`); a `buildDay` that declares fewer than five arguments and so never sees `procOpt` (`lsqvce:makeLoader:prepareUnused`); a `processingOptions` that already has a field `prepared` (`lsqvce:makeLoader:preparedField`) — `fit` reports both under its own prefix, the original attached as the cause |
| when it fails | `lsqvce:makeLoader:prepareDay` names the day; your error is the cause (an `lsqvce:*` error passes unchanged); a handle that does not resolve is reported the same way |

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
| `excluded` | m_a × 1 logical | *optional* — rows of **weight 0**: zero in `A`, `B` and `l`, kept for their index (residuals, diagnostics) but not observations — out of the degrees of freedom, `rejected` in the verdict |

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
| `lsqvce:makeLoader:buildDay` | `buildDay` threw; its error is attached as the cause (an `lsqvce:*` error passes unchanged) |
| `lsqvce:makeLoader:badDay` | a required field is missing or not a 1×nArc cell, or an unknown field is present, or `excluded` marks a row that is not zero |
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
`noDays`, `badDays`, `duplicateDay`, `badId`, `mixedIdTypes`, `duplicateSystem`, `builderKind`,
`badArcs`, `badScope`, …).
