# The `data` struct

Everything the toolbox knows about a problem enters through one plain MATLAB struct,
`data`, handed to `lsqvce.fit(data, opts)`. `fit` tells three **kinds** apart by the
fields present:

| kind | fields | for |
|---|---|---|
| **streaming** | `commonParams`, `days`, `systems` (+ `priors`, `paramColumns`) | arc-wise data read day by day through a `buildDay` function, local unknowns eliminated on the fly, estimated stochastic models |
| **matrices** | `A`, `l`, `p` (+ `priors`) | a design matrix that fits in memory, no local unknowns, diagonal weights |
| **batched** | `batchFcn`, `numBatches`, `numUnknowns` or `x0` (+ `priors`) | a design too large for one matrix, or a nonlinear model |

`data` is a **scalar struct**; field names are case-sensitive (lowerCamelCase). An
unknown field is an error (`lsqvce:fit:unknownField`), so a misspelt field cannot be
ignored silently.

## The streaming kind

| field | type / size | required | meaning |
|---|---|---|---|
| `commonParams` | cell R×C, C ≥ 1 | yes | the common parameter **templates** (below) |
| `days` | double 1×D | unless every system has its own `days` | the default calendar: integer day identifiers, each once (the loader cache is keyed on them; a repeated day is refused, `lsqvce:fit:duplicateDay`) |
| `systems` | struct 1×M, M ≥ 1 | yes | one entry per observing system (below) |
| `priors` | struct array | no | prior groups on the common parameters (chapter Priors) |
| `paramColumns` | struct | no | which table columns hold a row's validity and its name (below) |

### `commonParams` — the common parameter table

One row per **template**. **Column 1 is the scope token**, columns 2…C are your own
convention (degree, order, name, epoch …) and are passed through to `buildDay`:

| token | instantiated as | typical use |
|---|---|---|
| `'common'` | one column of A, always | spherical-harmonic coefficients |
| `'day'` | one column per calendar day that carries at least one arc in at least one system | a daily common parameter |

`'arc'` is not allowed in this table (`lsqvce:fit:badScope`): arcs belong to a system.
nc, the number of columns of A, is **derived** from the expansion and never declared.

**Column layout of A.** The `'common'` rows first, in table order; then one day block
per calendar day, in calendar order (the sorted union of every system's day list).
`buildDay` receives the position of the current day's block as
`info.commonDayOffset` / `info.commonDayCount`. A minimal table:
`[repmat({'common'}, nc, 1), num2cell((1:nc).')]`.

### `systems(k)` — one entry per observing system

All entries share the same fields; fill an unused optional field with `[]`.

| field | type / size | default | meaning |
|---|---|---|---|
| `buildDay` | function handle | *(required)* | your day builder (chapter `buildDay`) |
| `covariance` | struct | *(required)* | one field per observation type holding a covariance spec (chapter Covariance specifications). **The fieldnames are the observation types**, in processing order; `buildDay` is called once per (day, type). |
| `localParams` | cell R×C | `{}` | the local parameter templates (below); not declared → nl = 0 |
| `arcs` | double, struct or function handle | *(none)* | the arc schedule (chapter `buildDay`); required only by `'arc'` template rows |
| `days` | double 1×D | `data.days` | this system's calendar: integers, each once (`lsqvce:fit:duplicateDay`) |
| `id` | finite double scalar or single char | the index k | identifier, unique across systems; all ids of one kind -- numbers or characters (the index default is a number: give every system an id, or none; `lsqvce:fit:badId`) |
| `name` | char / string | `''` | label used in reports |
| `processingOptions` | struct | `struct()` | anything `buildDay` needs (paths, force models); its 5th argument |
| `priors` | struct array | none | prior groups on this system's local parameters (chapter Priors) |

`buildDay` + `covariance` is a complete system.

### `localParams` — the local parameter table

Column 1 is the scope token, columns 2…C your convention:

| token | one block per … | typical use |
|---|---|---|
| `'arc'` | (day, arc) of the system | initial state, per-arc empirical parameters |
| `'day'` | day of the system that carries ≥ 1 arc — every calendar day when no schedule is declared | accelerometer scale / bias / drift |
| `'system'` | whole adjustment | instrument constants |

`'arc'` rows require a declared arc schedule; `'day'` rows without a schedule live on
every calendar day of the system; `'system'`-only tables need neither. nl is derived.

**Column layout of B (per system).** Per day, in calendar order: the arc blocks in
sorted arc-id order (each block = the `'arc'` rows in table order), then the `'day'`
block; after all days one trailing `'system'` block. A day whose declared schedule
names no arc gets no block at all. `buildDay` never computes an offset: every offset
arrives in its `info` argument.

A template has the same rows for every arc. Arcs with genuinely different parameter
counts are an engine-tier feature (`lsqvce.engine.makeLoader` with `LocalsPerArc`).

### `paramColumns` — validity and names

| field | meaning |
|---|---|
| `start`, `end` | the columns of both tables that hold a row's validity (`[MJD seconds]` pairs or fractional MJD; `[]` = always valid; any other entry, or a start after its end, is refused: `lsqvce:fit:badValidity`). A `'day'` / `'arc'` row is then instantiated only on calendar days its validity overlaps (start < d+1 and end > d); `'common'` / `'system'` rows are always present. |
| `name` | the column holding a row's name, matched by the `group` field of a prior (default 2) |

Without `paramColumns` every row is valid always. `buildDay`'s `info` then lists which
template rows exist on the day (`info.arcRows`, `info.dayRows`, `info.sysRows`,
`info.commonDayRows`).

### `lsqvce.layout(data)` — the columns announced

```matlab
lay = lsqvce.layout(data);
lay.nc, lay.common(j)                          % template row (and day) of column j of A
lay.systems(k).nl, lay.systems(k).local(c)     % (row, day, arc) of column c of B
c = lay.systems(k).col(row, day, arc);          % lookup, 0 = not instantiated
c = lay.colCommon(row, day);
```

This is the column map a prior operator `L` is written in, and the one `R.toTable`
uses to write estimates back into your tables.

## The matrices kind

| field | size | meaning |
|---|---|---|
| `A` | m×u | design matrix |
| `l` | m×1 | observed minus computed |
| `p` | m×1 | weights (the diagonal of P); absent = unit weights |
| `priors` | struct array | prior groups (chapter Priors) |

## The batched kind

| field | size | meaning |
|---|---|---|
| `batchFcn` | function handle | `[A, l, p] = batchFcn(k)`, k = 1…`numBatches`: one batch of design rows, observations and weights. With `x0`: `[A, w, p] = batchFcn(k, x)` returns the Jacobian at `x` and `w = l − f(x)`. |
| `numBatches` | scalar | the number of batches |
| `numUnknowns` | scalar | u, linear form only |
| `x0` | u×1 | the starting point; makes the problem **nonlinear** (damped Gauss–Newton). u = `numel(x0)`; `numUnknowns` must then be absent (`lsqvce:fit:uCount`). |
| `priors` | struct array | prior groups (chapter Priors) |

Both in-memory kinds run with `Estimate = "none"` (variance components need the
streams of a streaming system, `lsqvce:fit:modelKind`); `Robust` selects
M-estimation. Their per-observation diagnostics are in the whitened frame,
`R.diag.v = sqrt(p).*v`.
