# Declaring the data

This chapter is the binding data contract of the toolbox (formerly the file `DATA_CONTRACT.md`, retired 2026-09-04: the Manual is the single authority). Everything the toolbox knows about your problem enters through the `data` struct described here.


**Status.** Binding specification, **version 6** (adopted 2026-09-09, Review R6 WP-102:
ONE shape changed — `buildDay` returns one CELL per arc, §0c — and the arc schedule became
optional; everything the v5 declaration says about tables, priors and `Cep` is unchanged.
Option names are those of `lsqvce.Options` v5, the v5 → v6 map is in `RELEASE_NOTES.md`).
`lsqvce.fit` is live for all three data kinds, `Estimate = "none" | "noise"`, `Robust`, the
nonlinear form and `OnPass`. §9 shows the equivalent direct engine calls. This file is the
source of the `fit` help text and of the HTML guide: one authority, no second copy.

Created 2026-08-31 (CEST, UTC+02:00) with assistance from Claude (Anthropic), model Claude
Fable 5. Rewritten for v3 and again for v4 (scoped tables) the same day (Claude Fable 5);
v6 (the cells are the arcs) 2026-09-09, Claude Fable 5.1.

---

### 0a. The v5 delta (R5 M5): three additions, no shape changes

1. **`data.priors(k) = struct('group' | 'L', 'x0', 'sigma', 'from', 'name', 'estimate')`** — declared
   pseudo-observation groups `L·x ~ x0` with weight `1/sigma²` on the COMMON parameters,
   folded into the solve AND into the diagnostic sweep (they carry a redundancy). A free
   network's datum is the prior `L = B'` on the null space: the exact minimum-constraint
   solution, for any weight, because the data's range and the null space are orthogonal
   (`lsqvce:fit:badPrior` on malformed fields). A Kaula rule or a background model is the
   same declaration with a scientific `L`. With `estimate = "scale"` (R5 M7 7b) the group has
   its OWN variance component: under `Estimate="noise"` its sigma is updated every pass by the
   Helmert / Koch–Kusche rule (Ω_k = ‖(L x − x0)/σ‖², r_k = p − tr(L N⁻¹ L')/σ², σ² ← σ²·Ω_k/r_k;
   the fixed point Ω_k = r_k is "the prior's weighted residual sum equals its redundancy"), so
   λ = s0²/σ_k² is the regularisation the data support; `R.noise.priors` and
   `R.noise.priorHistory` record it. `from` is the pass at which the group becomes active
   (default 1; a single-pass run reports `lsqvce:fit:priorNotActive` and skips it). A group is
   declared EITHER by `group` — a regular expression on the NAME column of `commonParams`
   (`data.paramColumns.name`, default 2); fit builds L, one row per matched row per instance —
   OR by `L` (p × nc in the columns of A); never both (R5-D8 round 6). One `sigma` per group:
   the rows of L are the members, and several rows are what make `estimate="scale"` estimable
   (a single row has redundancy ≤ 1).
1b. **`data.systems(k).priors(l) = struct('group' | 'L', 'x0', 'sigma', 'from', 'name', 'estimate')`**
   (R5-D8 round 5/6) — the priors on system k's LOCAL parameters, the literal counterpart of
   `data.priors(j)`: one struct per GROUP, declared EITHER by **`group`** — a regular expression
   matched against the NAME column of the table (column 2; `data.paramColumns.name` to say
   otherwise); the toolbox builds L with one row per matched template row per instance — OR by
   **`L`** (q × nl in the columns of that system's B as **`lsqvce.layout(data)`** announces them,
   used as is); never both. Then `x0`, ONE `sigma` per group (relative weights inside a group go
   into L), `from`, `name` (default the pattern), `estimate = "none" | "scale"`. Every row must lie inside one day's local block
   (couple days through common parameters); the width of L is checked by name
   (`lsqvce:fit:badPrior`). Each group rides as its own observation stream `prior_<name>` of the
   system (pseudo-observations A = 0, B = L/σ, l = x0/σ, unit law), so accumulation, elimination,
   the sweeps, `residuals()` and the redundancy count treat it as data — and with
   `estimate = "scale"` its σ is estimated every pass of `Estimate="noise"` as that stream's level,
   the Helmert / Koch–Kusche factor the joint estimator computes for every stream: IDENTICAL
   behaviour to the common groups (`R.noise.model{k}.sigma` for stream `satKey|prior_<name>`).
   Common priors and every system's priors carry their own settings. THE RULE (R5-D8): the
   parameters and their columns are yours (tables, `buildDay`, plots, `toTable`); regularisation is
   the toolbox's (mechanism, redundancy, activation, VCE per group); the interface between them is
   an operator in your columns, announced by `lsqvce.layout`, checked by width.
   `Options.Reg / RegLambda / RegMatrix` are gone.
1c. **`data.paramColumns = struct('start', c1, 'end', c2)`** (R5-D8 round 3) — names the columns
   of BOTH tables that hold a row's validity (`[MJD seconds]` pairs or fractional MJD; `[]` =
   always). A `'day'` / `'arc'` row is then instantiated only on the calendar days its validity
   overlaps (half-open day: start < d+1 and end > d); `'common'` / `'system'` rows are always
   present. nc and nl are derived accordingly, and `buildDay`'s 6th argument tells you which
   template rows exist on the day: `info.arcRows`, `info.dayRows`, `info.sysRows`,
   `info.commonDayRows` (indices into the tables, in block order). Without `paramColumns`
   every row is valid always — the v4 behaviour. This is how a one-day manoeuvre parameter and
   a "literal" row become the same mechanism as everything else.
2. **Per-stream policy on the covariance spec** — `covariance.<type>.estimate = "none" |
   "scale" | "law"` (default `"law"` under `Estimate="noise"`) and `.arcFactors` (logical,
   default the Options policy). A law whose `sigma` is a VECTOR (per-row weights) is a fixed
   structure: it can carry a level (`"scale"`) but not a spectral shape (`"law"` is refused,
   `lsqvce:fit:lawOnFixedStructure`). `Options.Pool` decides which streams SHARE one
   estimated spectral shape; every stream keeps its own level. To read an arc's standard deviation back, use all
   three terms: **declared sigma × sqrt(arc factor) × `R.sigma0`**. Since 6.1.6 the factor
   is the arc's absolute Helmert factor and **weights the solve** (review R7-13), so
   `sigma0` settles near 1 and carries only what no component absorbed, as in any GLS
   adjustment (§7.7, §7.8 item 2). (R6 finding F-01, fixed in WP-98: eq. 6.5.28 multiplies
   by the factor the weighting used, so the factors settle across passes instead of growing
   8.4x a pass. Runs with `arcFactors = false` are unchanged.)
3. **The per-pass inspection hook** — `Options.OnPass = @(ctx)` receives, after every pass
   of every iterative run, the pass state, the live accumulator, `ctx.data` (this
   declaration as given) and `ctx.processingOptions{k}` (system k's `processingOptions` as
   declared); return `false` to stop. Nothing else about `buildDay` or the tables changed.

Semantics made explicit in v5: a full SPD weight matrix returned by a batch source is
whitened exactly (`chol(P)`), and per-observation diagnostics of a batch problem are
reported in the WHITENED frame (`R.diag.v = sqrt(p).*v`; only their sums are
frame-invariant).

### 0b. The v5.3 delta (R6 WP-90): the error information travels with the day

What the preprocessing knows about the accuracy of every observation — a per-epoch
standard deviation, per-axis standard deviations, or a full 3 × 3 covariance per epoch — is
DATA, and it comes back from `buildDay` WITH the day, as **`D.Cep`** (§3.4, §7), for ANY
stream: `nEpoch × 1` for a scalar stream, `nEpoch × d` or `nEpoch × d × d` for a multi-axis
one. The axis count is derived from it. The covariance spec (§5) declares the TEMPORAL law
only — `diag` (white), `arma`, `psd`, `acf` — and the toolbox composes
C = D (R_t ⊗ I_d) D per arc. Retired with this: the `separable` type (it was `arma` + Cep),
the `epochCovariance` hook (an unknown field now), and a static `Cep` on the spec (refused
when the day carries one: one fact, one authority). Nothing else about the declaration
changed.

### 0c. The v6 delta (R6 WP-102): the cells are the arcs, and the schedule is optional

`buildDay` returns every field as a **1 × nArc cell, one cell per arc of the day**: `A`, `B`,
`l`, `t` and, when given, `Cep` (§3.3, §3.4). A day with one arc returns `{A}`; a day
without data returns empty cells. The v5 row-index fields `arcIdx` / `arcId` are gone — the
cells say what they said — and an **empty cell is a data gap** (§4.3). The driver was the
GRAVFIRE preprocessing, which delivers every day as a cell array, one cell per
noise-independence piece: the day's cells now ARE the arcs, with no bookkeeping in between.

`arcs` is **optional**. Without it, the arcs of a day are its cells (ids 1 … nArc) and the
toolbox learns their number when the engine's per-day `meta` probe reads the day — which the
engine reads right after, so the day cache makes the discovery free (`CacheDays` counts DAYS:
one resident day is one slot per stream). `'day'` template rows then live on **every**
calendar day, and a calendar day that returns no arc in any stream is refused at its probe
(`lsqvce:makeLoader:emptyDay`) — take it off the calendar, or declare a schedule. A schedule
is **required** only by `'arc'` template rows (nl counts one block per (day, arc) and must be
known before any day is read); declared, cell a is the a-th scheduled id and the count must
match (`lsqvce:makeLoader:arcCount`). The v5 string form `'fromDay'` is retired: omitting
`arcs` is what it meant.

### 0. The v4 principle: state nothing the data already states — and scope, don't number

The contract declares **structure**, never **counts**, and since v4 it declares **scope,
never numbering**: column 1 of both parameter tables is a **scope token**, the rows are
**templates**, and the toolbox expands them over the calendar. Every count is derived:

| quantity | derived from |
|---|---|
| nc (common unknowns) | the expansion of `data.commonParams` (`'common'` rows once + `'day'` rows per calendar day with arcs) |
| nl (local unknowns) | the expansion of `localParams` (`'arc'` rows per (day, arc) + `'day'` rows per day with arcs + `'system'` rows once); **not declared → nl = 0** |
| observation types | `fieldnames(covariance)` |
| u (batched nonlinear) | `numel(x0)` |

The v2 count fields are gone (v3), and so is v3's hand-numbering split (`localsPerArc` XOR
free `localParams`): **one mechanism**, the scoped table. Every column offset reaches
`buildDay` through its 6th argument (`info`, §6.3) — no offset is ever computed in user
code. The driver was the GRACE parameterization: arc-wise ICs, day-wise ACC
scale/bias/drift, per-arc KBR empiricals — three scopes, one table.

### 1. The mental model in one paragraph

You describe **one adjustment problem** with a plain MATLAB struct called `data`. It answers
three questions: *what is estimated for everybody* (`data.commonParams`, one row per common
unknown), *which days are processed* (`data.days`), and *which observing systems deliver
data* (`data.systems`, one entry per satellite or instrument). Each system points to **one
function you write**, `buildDay`, that returns the design matrices and observations of a
single day. Everything else — reading days in the right order, cutting them into arcs,
whitening with your covariance model, iterating the variance components — is the toolbox's
job:

```matlab
opts = lsqvce.Options(Estimate="noise");  % what to do
R    = lsqvce.fit(data, opts);           % the one verb
R.report                                 % the verdict
```

You never call `makeLoader` yourself; `fit` builds the loaders from `data`. `R.systems`
(the existing Results field) records what ran: per system its id, nl, the loader closure,
the day list and the name.

---

### 2. The `data` struct — field by field

`data` is a **scalar struct**. Field names are **case-sensitive** and use lowerCamelCase.
Unknown field names are an **error** (`lsqvce:fit:unknownField`), which is deliberate: a
typo in a struct field would otherwise be silently ignored. The retired v2 names (`common`,
`numCommon`, `numLocal`, `observationTypes`, `systemId`, `cache`, `loader`, the plot hooks)
are unknown fields now — the error names the field, so a v2 declaration fails loudly, not
subtly.

#### 2.1 `data.commonParams` — the shared unknowns (required)

| field | type / size | meaning |
|---|---|---|
| `commonParams` | cell, R×C (C ≥ 1) | one row per common **template**. **Column 1 is the scope token**: `'common'` (one column of A, always) or `'day'` (one column per calendar day that carries at least one scheduled arc in at least one system). Columns 2…C are *your* convention (degree, order, epoch, …). nc is **derived** from the expansion; the table is piped to `buildDay` as its 3rd argument. `'arc'` is forbidden here — arc sets are system-specific (`lsqvce:fit:badScope`). |

**Binding column layout of A:** the `'common'` rows first (table order), then one day
block per calendar day (calendar order — the sorted union of every declared day list).
`buildDay` receives this day's block position as `info.commonDayOffset` /
`info.commonDayCount` (§6.3). Minimal table when you have nothing to say:
`[repmat({'common'}, nc, 1), num2cell((1:nc).')]`.

#### 2.2 `data.days` — the default schedule

| field | type / size | meaning |
|---|---|---|
| `days` | double, 1×D | day identifiers, e.g. `1:30` or `datenum` values cast to integers. Every system uses this list unless it carries its own `days`. Required unless every system has one. |

Days must be **numeric integers**: the loader cache is keyed on the number.

#### 2.3 `data.systems(k)` — one entry per observing system

`data.systems` is a **1×M struct array** (M ≥ 1). All entries have the same fields; fill an
unused optional field with `[]`.

**Required — two fields**

| field | type / size | meaning |
|---|---|---|
| `buildDay` | function handle | your day builder, §3. |
| `covariance` | struct | one field per observation type holding a covariance spec (§5), e.g. `cov.kbr = struct('type','diag','sigma',0.2)`. **The fieldnames ARE the observation types** (their order is the processing order); `buildDay` is called once per (day, type). |

**Optional**

| field | type / size | default | meaning |
|---|---|---|---|
| `arcs` | double or fh | *(none)* | the arc schedule, §4: omitted, the cells `buildDay` returns are the arcs of the day; required only by `'arc'` template rows. |
| `localParams` | cell, R×C (C ≥ 1) | `{}` | the local **templates**, §6: column 1 is the scope token `'arc'` \| `'day'` \| `'system'`, columns 2…C your convention; nl is **derived** from the expansion. Passed to `buildDay` as its 4th argument. Not declared → **nl = 0** (§6.4). |
| `days` | double, 1×D | `data.days` | per-system schedule override. |
| `id` | double 1×1 **or** single char | the index k | system identifier, unique across systems (explicit ids and index defaults share one namespace). Longer labels go in `name`. |
| `name` | char / string | `''` | human-readable label used in reports. |
| `processingOptions` | struct | `struct()` | anything your `buildDay` needs (paths, force models); passed as its 5th argument. |

**Rule of thumb for beginners:** `buildDay` + `covariance` is a complete system — return one
cell per arc and the arcs are declared. Add `'arc'` template rows in `localParams` when arcs
carry their own unknowns (and then `arcs`, the schedule those rows need), the rest when you
need it. (The v5.2 `epochCovariance` hook is gone —
the per-epoch error information comes back from `buildDay` as `D.Cep`, §7.)

#### 2.4 A pitfall you will hit once: `struct()` with cell values

```matlab
% WRONG: struct() expands cell values into a struct ARRAY
s = struct('localParams', {'bias','drift'});          % s is 1x2
% RIGHT: double braces, or assign field by field
s = struct('localParams', {{'bias';'drift'}});        % s is 1x1
s.localParams = {'bias'; 'drift'};                    % clearest
```

---

### 3. The `buildDay` contract

#### 3.1 Who calls it, and how often

`fit` calls `buildDay` **once per (day, observationType)**; the cells it returns are the
day's arcs. Every pass over the data reads each day **twice** (once to accumulate the
normal equations, once for the residuals), and the variance-component loop multiplies that
by the number of passes. Two consequences:

1. **Determinism is mandatory.** Same arguments → byte-identical output. Random numbers,
   clocks, mutable globals are forbidden inside `buildDay`. (The residuals are formed as
   `v = l − A·x − B·y` from the second read; any difference lands in σ₀.)
2. **Cache expensive reads yourself.** The toolbox caches *one built day* (every stream of
   it; `Options.CacheDays` for more); if the raw file read is the expensive part and serves
   several observation types, keep your own cache in the closure that binds `buildDay`.

#### 3.2 Signature — declare only what you need

`fit` passes as many of these six arguments as your handle **declares** (`nargin`):

```matlab
D = buildDay(day, obsType)                                        % classic
D = buildDay(day, obsType, common, local, procOpt)                % parameter-driven
D = buildDay(day, obsType, common, local, procOpt, info)          % with the expansion info
```

| argument | what it is |
|---|---|
| `day` | one entry of the day list (numeric) |
| `obsType` | one entry of `fieldnames(covariance)` (char) |
| `common` | `data.commonParams` |
| `local` | `localParams` of this system (or `{}`) |
| `procOpt` | `processingOptions` of this system |
| `info` | this day's expansion slice — every offset of A and B (§6.3); `[]` only when nothing is scoped |

A `varargin` builder receives all six.

#### 3.3 Return value — required fields: THE CELLS ARE THE ARCS

Every field is a **1 × nArc cell**, one cell per arc of the day; cell a is arc a. A day with
one arc returns `{A}`, `{B}`, `{l}`, `{t}`; a day without data returns empty cells (nArc = 0).

| field | size | meaning |
|---|---|---|
| `A` | cell, 1×nArc; each m_a×nc | design matrix of the **common** unknowns of arc a |
| `B` | cell, 1×nArc; each m_a×nl | design matrix of the **local** unknowns; may be `sparse` |
| `l` | cell, 1×nArc; each m_a×1 | **observed minus computed** (the reduced observation). Not the raw observation. |
| `t` | cell, 1×nArc; each m_a×1 | epoch of each row in **seconds** (the whitener needs the lags) |

Arc lengths may differ freely; the whitener is built per arc. An arc whose four entries have
**0 rows is a data gap** (§4.3) and is skipped. Return **raw, unwhitened** matrices. Whitening
happens per arc inside the toolbox from your covariance spec, and the VCE loop reweights there
too — a builder that scales its own rows would be reweighted twice.

#### 3.4 Return value — optional fields

| field | size | meaning |
|---|---|---|
| `Cep` | cell, 1×nArc; each nEpoch_a×1, nEpoch_a×d or nEpoch_a×d×d | **the error information of each arc** (§7): the per-epoch standard deviation of a scalar stream, the per-axis standard deviations, or the full per-epoch covariance of a multi-axis stream. nEpoch_a·d must equal the arc's row count; the axis count d is derived from it and is one for the day. |
| `nAxis` | 1×1 | the axis count of a multi-axis stream WITHOUT `Cep` (rows epoch-major, m_a = nEpoch_a·nAxis, §7). Derived from `Cep` when that is present; given as well, it must agree. |

No other field is allowed — `arcIdx` / `arcId` (the v5 way of saying what the cells say) or a
misspelt `cep` are refused by name (`lsqvce:makeLoader:badDay`), so a typo cannot be ignored
silently.

#### 3.5 What is checked, and the error you get

| error id | meaning |
|---|---|
| `lsqvce:makeLoader:buildDay` | your `buildDay` threw; the original error is attached as the cause |
| `lsqvce:makeLoader:badDay` | a required field (`A`, `B`, `l`, `t`) is missing or is not a 1×nArc cell, or the day carries a field that is not `A` `B` `l` `t` `Cep` `nAxis` |
| `lsqvce:makeLoader:dayShape` | the cells do not agree in number, or an arc's `A`, `B`, `l`, `t` do not share a row count |
| `lsqvce:accum:arcBadNc` / `arcBadNl` | `A` or `B` has the wrong number of columns — the derived nc/nl are enforced against **every** arc |
| `lsqvce:makeLoader:arcCount` | a declared schedule names n arcs for the day and `buildDay` returned a different number of cells |
| `lsqvce:makeLoader:emptyDay` | without a schedule, a calendar day that owns a `'day'` block returned no arc in any stream |
| `lsqvce:makeLoader:badArc` | an arc id outside the declared schedule was asked for (a static schedule declared wrongly), or an arc of a multi-axis day does not hold whole epochs |
| `lsqvce:makeLoader:badDayCep` | `Cep` is not a 1×nArc cell, or an arc's entry is not nEpoch×1 / nEpoch×d / nEpoch×d×d, or nEpoch·d is not the arc's row count |
| `lsqvce:makeLoader:axisConflict` | `nAxis` disagrees with `Cep` (or with a `d` on the spec), or does not divide the row count |
| `lsqvce:makeLoader:cepTwice` | the covariance spec carries a `Cep` AND the day returns one |
| `lsqvce:covModel:cepOnVar` | a `var` law (its own cross-axis structure) with a `Cep` |
| `lsqvce:makeLoader:noCep` | *(warning, once per run)* a multi-axis day without `Cep`: the axes are weighted by the temporal law alone (D = I) |

These are raised **while reading**, by the engine, and keep their engine ids. Errors in the
*declaration* itself are raised by `fit` before any data is read and carry `lsqvce:fit:` ids
(`unknownField`, `noSystems`, `noCommon`, `missingField`, `localCount`, `noDays`, `badId`,
`duplicateSystem`, `builderKind`, `badArcs`, and every `makeLoader` construction id
relabelled, e.g. `lsqvce:fit:badScope`).

---

### 4. Arc schedules

#### 4.0 None — the cells are the schedule (the default since v6)

Omit `arcs`. The arcs of a day are the cells `buildDay` returns, ids 1 … nArc, and days may
differ freely in their count. The toolbox learns the count when the engine's per-day `meta`
probe reads the day (the engine reads it right after; the cache keeps a whole day resident,
so nothing is read twice). Consequences:

- `'day'` template rows (§6) and `'day'` rows of `commonParams` are instantiated on **every**
  calendar day of the system — the calendar is the declaration.
- A calendar day that returns **no** arc in any stream is refused at its probe
  (`lsqvce:makeLoader:emptyDay`) before any solve could go rank-deficient on an unobserved
  block. Take the day off the system's calendar (`days`), or declare a schedule that names no
  arc for it.
- `'arc'` template rows need a schedule (§4.1 / §4.2): nl counts one block per (day, arc)
  and must be known before any day is read (the `localPlan` construction error otherwise).

#### 4.1 Static — every day looks the same

```matlab
s.arcs = 4;            % arcs 1..4 on every day, every observation type
s.arcs = [1 2 4];      % these ids on every day
s.arcs.kbr = 2;  s.arcs.pos = 4;    % per observation type
```

Cell a of the returned day is arc ids(a); the day must return exactly as many cells as the
schedule names (`lsqvce:makeLoader:arcCount`), which catches a mis-declared schedule.

#### 4.2 Ragged — days differ (gaps, manoeuvres, short months)

Declare a **resolver** that returns the arc ids of a day *without reading data*:

```matlab
s.arcs = @(day, obsType) arcTable{day};      % 1xA vector of ids for that day
```

and return one cell per id, in that order. This is the form `'arc'` template rows need: a day
whose list is empty gets no local blocks and no common day block (§6.4).

#### 4.3 Data gaps — the empty cell

Under a declared schedule, an arc the day does not carry is returned as an **empty cell**
(0 rows in `A`, `B`, `l`, `t`): every engine skips it, no padding, no error. The same run
declared with a resolver that simply omits the id gives bitwise the same result. (An empty
cell of an arc that carries `'arc'` template rows leaves those unknowns unobserved — say the
gap in the resolver instead.)

The v5 string form `arcs = 'fromDay'` is retired (`lsqvce:fit:badArcs`): omitting `arcs` is
what it meant.

---

### 5. Covariance specs — the `covariance.(obsType)` structs

One struct per observation type; `type` selects the **temporal law** — how the errors of
one stream are correlated down the epochs. The most common:

| `type` | fields | use for |
|---|---|---|
| `'diag'` | `sigma` (scalar std) | uncorrelated noise, e.g. `struct('type','diag','sigma',0.2)` |
| `'arma'` | `ar` `[1 a1 … ap]`, `ma` (default 1), `sigma` (innovation std) | coloured noise with a known filter |
| `'acf'` | `acf` `@(dt,opt)`, `sigma` | correlation given as a function of the lag |
| `'psd'` | `f` (nf×1, Hz), `S` (nf×1) | noise given as a one-sided PSD |
| `'var'` | `ar` (d×d×p), `sigma` (d×d) | vector-AR for multi-axis streams with cross-axis temporal correlation (no `Cep`) |

`sigma` is always a **standard deviation**, never a variance, and it is the LEVEL of the law.
The accuracy of the individual observations — what varies from epoch to epoch — is not
declared here: it comes back from `buildDay` with the day as `D.Cep` (§7), and the toolbox
composes the two. The VCE loop scales and reshapes these specs at run time; you only give
the starting point. Full reference: `doc lsqvce.engine.covModel`.

---

### 6. Local unknowns — the scoped `localParams` table

One mechanism: the table's rows are templates, column 1 says **how often** each row
instantiates.

#### 6.1 The three scopes

```matlab
s.localParams = [ ...
    [repmat({'arc'}, 6, 1),  num2cell((1:6).')];   % 6 ICs   -> once per (day, arc)
    {'arc','kbrBias'; 'arc','kbrCos1cpr'};         % 2 more  -> once per (day, arc)
    {'day', 'accScale'; 'day', 'accBias'};         % 2       -> once per day with arcs
    {'system', 'antennaOffset'}];                  % 1       -> once, full adjustment
```

| token | one block per … | typical use |
|---|---|---|
| `'arc'` | (day, arc) of this system | initial conditions, per-arc empiricals |
| `'day'` | day whose schedule carries ≥ 1 arc — every calendar day when no schedule is declared | ACC scale/bias/drift |
| `'system'` | whole adjustment | instrument constants; also the v3 free table (see migration) |

Requires `days` whenever `'arc'`/`'day'` rows are present, and a declared schedule (static or
resolver, §4) for `'arc'` rows; `'day'` rows without a schedule live on every calendar day
(§4.0); `'system'`-only tables need neither. nl is derived; declaring `numLocal` too is the
conflict error.

**Binding column layout of B (per system):** per day (declared order) — the arc blocks in
**sorted arc-id order** (each block = the `'arc'` rows, table order), then the `'day'`
block — and after all days one trailing `'system'` block. A day whose declared schedule has
no arcs gets **no blocks at all** (§6.4).

#### 6.2 Heterogeneous per-arc counts stay ENGINE-tier

A template has the same rows for every arc. If arcs genuinely carry *different* counts,
that is the ENGINE-tier `lsqvce.engine.makeLoader(..., 'LocalsPerArc', @(day, arcId) n)`
(the WP-52 plan) — deliberately not on the `fit` surface.

#### 6.3 Using `info` (buildDay's 6th argument)

```matlab
info.arcId            % 1xA  ids of this day's arcs
info.offsetOf         % @(arcId) -> column offset of that arc's block in B
info.count            % 1xA  columns of each arc block (= number of 'arc' rows)
info.dayOffset        % 1x1  this day's 'day' block offset ([] if none)
info.dayCount         % 1x1  its width (0 if none)
info.sysOffset        % 1x1  the 'system' block offset
info.sysCount         % 1x1  its width
info.nl               % 1x1  total local count of the adjustment
info.nc               % 1x1  total common count
info.commonCount      % 1x1  width of A's leading 'common' block
info.commonDayOffset  % 1x1  this day's common day-block offset in A ([] if none)
info.commonDayCount   % 1x1  its width (0 if none)
```

Fill `B(rows, info.offsetOf(id) + (1:k))`, `A(:, info.commonDayOffset + (1:w))` — and
**never compute an offset yourself**.

#### 6.4 Days without arcs, systems without locals

A calendar day on which **no system** schedules an arc gets no common day block (nc
shrinks); a system day whose declared schedule has no arcs gets no arc and no `'day'` blocks
(nl shrinks) — offsets stay gap-safe on ragged schedules and manoeuvre days. A system
**without** a schedule claims every day of its calendar (§4.0). A system without local
unknowns is still the **default**: declare no `localParams` and return `B = zeros(m, 0)`
(the engine skips the pre-elimination; empty diagnostics per HANDOVER G-07).

---

### 7. The error information — `D.Cep` on the returned day

The covariance of a stream has two parts. The **temporal law** R_t (white, AR, a PSD) is
declared in §5 and estimated by `Estimate="noise"`. The **per-epoch accuracy** D — what the
preprocessing knows about each observation: a PPP covariance, a per-epoch sigma from a
fit, a weight from an earlier screening — is DATA, and you return it with the day:

```matlab
D.Cep{a} = sigmaPerEpoch;                  % nEpoch_a x 1      scalar stream
D.Cep{a} = [sx sy sz];                     % nEpoch_a x 3      per-axis sigmas (diagonal D)
D.Cep{a} = permute(reshape(Q.', 3, 3, []), [3 1 2]);   % nEpoch_a x 3 x 3  full per-epoch covariance
```

(one entry per arc, a 1×nArc cell like `A`; v6)

The toolbox composes C = D (R_t ⊗ I_d) D per arc: every epoch is de-scaled by D^-1/2 (for
d = 1 that is a division by the sigma), then the temporal whitener runs down the epochs of
each axis. The axis count d is `size(Cep, 2)`; `nAxis` is derived from it. Rows of a
multi-axis stream are **epoch-major** (x₁ y₁ z₁ x₂ y₂ z₂ …), so nEpoch·d must equal the row
count, and an arc must hold whole epochs. Under `Estimate="noise"` the estimator de-scales
the residuals of every arc by the same D before it forms the lag statistics, so the
estimated law is the law of the *normalised* errors — a scalar stream with a per-epoch
sigma is handled exactly like a 3-axis one.

Rules:

- One authority: `Cep` lives on the day. A `Cep` on the covariance spec together with one
  on the day is refused (`cepTwice`); a `d` on the spec that disagrees with the day is
  refused (`axisConflict`).
- A multi-axis day **without** `Cep` (set `nAxis`) is composed with D = I and reported once
  per run (`lsqvce:makeLoader:noCep`) — the omission is distinguishable from a choice.
- `sigma` on the spec stays the LEVEL of the law; `Cep` is relative to it. Declare `sigma`
  1 when `Cep` already carries the absolute accuracy.
- A `var` law carries its own cross-axis structure and takes no `Cep`.

---

### 8. Worked example — the GRACE parameterization: three scopes, one table

Arc-wise initial conditions, **day-wise** accelerometer scale/bias/drift, per-arc KBR
empiricals — the case v3 could only hand-number, and the driver of v4:

```matlab
% ---- the one function you write --------------------------------------
function D = buildDayGrace(day, obsType, common, ~, procOpt, info)
    raw = readGraceDay(procOpt.dataDir, day, obsType);   % your reader, deterministic;
                                                         % one entry per arc of the schedule
    nArc = numel(info.arcId);                            % the resolver's ids, in order
    D = struct('A', {cell(1, nArc)}, 'B', {cell(1, nArc)}, 'l', {cell(1, nArc)}, 't', {cell(1, nArc)});
    for a = 1:nArc                                       % THE CELLS ARE THE ARCS (v6)
        arc = raw.arc(a);                                % empty on a data gap -> empty cell
        D.A{a} = shPartials(arc, common);                % m_a x info.nc; the 'common' block
        D.l{a} = arc.observed - arc.computed;            % O - C
        D.t{a} = arc.tSeconds;
        B = sparse(numel(D.l{a}), info.nl);
        off = info.offsetOf(info.arcId(a));
        B(:, off + (1:info.count(a))) = icAndKbrPartials(arc);
        if info.dayCount > 0                             % this day's ACC block
            B(:, info.dayOffset + (1:info.dayCount)) = accPartials(arc);
        end
        D.B{a} = B;
    end
end

% ---- the declaration ---------------------------------------------------
data.commonParams = [repmat({'common'}, 3717, 1), shParameterTable(60)];
data.days = 1:31;
cov.kbr = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',1e-7);
cov.pos = struct('type','diag','sigma',0.02);             % fieldnames = the two types

s = struct();
s.buildDay = @buildDayGrace;
s.covariance = cov;
s.id = 'C';
s.name = 'GRACE-FO C';
s.arcs = @(day, obsType) arcSchedule('C', day);            % required by the 'arc' rows below
s.localParams = [ ...
    [repmat({'arc'}, 6, 1), {'x0';'y0';'z0';'vx0';'vy0';'vz0'}];  % ICs per arc
    {'arc', 'kbrBias'; 'arc', 'kbrCos1cpr'};                      % KBR per arc
    [repmat({'day'}, 3, 1), {'accScale';'accBias';'accDrift'}]];  % ACC per day
s.processingOptions = struct('dataDir', '/data/gracefo');
data.systems(1) = s;

s.id = 'D';  s.name = 'GRACE-FO D';  s.arcs = @(day, obsType) arcSchedule('D', day);
data.systems(2) = s;

% nl per satellite = SUM over its days-with-arcs of (arcs*8 + 3) -- derived,
% gap-safe on manoeuvre days, never stated.

% ---- the run -----------------------------------------------------------
opts = lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, Parallel=struct('Mode', "parfor"));
R = lsqvce.fit(data, opts);
R.report
```

---

### 9. What `fit` does internally — the equivalent direct engine calls

```matlab
L = cell(1, numel(data.systems));
for k = 1:numel(data.systems)
    s = data.systems(k);
    L{k} = lsqvce.engine.makeLoader(s.buildDay, ...
        ObservationTypes=fieldnames(s.covariance).', Covariance=s.covariance, ...
        SystemId=s.id, Name=s.name, Arcs=s.arcs, ...                 % Arcs only when declared
        LocalParams=s.localParams, Scoped=true, CommonPlan=plan, ...  % fit builds `plan`
        Days=data.days, CommonParams=data.commonParams, ...           % over the global calendar
        ProcessingOptions=s.processingOptions);
end
R = lsqvce.engine.adjustBlock(L, data.days);                        % Estimate "none"
R = lsqvce.engine.adjustRobust(L, data.days, 'Influence', "huber");  % Robust "huber"
% Estimate "noise" is lsqvce.engine.adjustVCEjoint over an accumulator factory that
% lsqvce.fit builds (jointModelOpt + the block accumulator); use fit for it.
```

These are ENGINE-tier calls (experts); a user needs none of them. The ENGINE tier keeps the
v2/v3 conveniences (`numCommon`/`numLocal` counts, free un-scoped tables with
`Scoped=false` — the default — `LocalsPerArc` incl. the heterogeneous `@(day,arcId)` form,
`cache`, raw 4-arg loader closures) — they left the USER surface, not the toolbox.

Field names map 1:1 to `makeLoader`'s Name-Value names (lowerCamelCase field ↔
UpperCamelCase name); `fit` performs exactly this mapping internally, plus the three
derivations of §0.

---

### 10. Checklist before the first run

- [ ] `buildDay` returns `A`, `B`, `l`, `t` as one cell per arc, each arc with one row count;
      `l` is O − C; nothing is whitened; no field beyond `Cep` / `nAxis`.
- [ ] Calling `buildDay` twice with the same arguments gives identical output.
- [ ] `t` is in seconds; the sampling matches the covariance spec you declared.
- [ ] `data.commonParams` rows carry the scope token in column 1; `size(A, 2)` must equal
      `info.nc` (use `zeros(m, info.nc)` and place blocks by the `info` offsets).
- [ ] Every observation type is a field of `covariance` — the fieldnames are the type list.
- [ ] Arc schedule: none (the cells are the arcs) unless `'arc'` template rows need one;
      declared, one cell per scheduled id, an empty cell for a gap.
- [ ] Local unknowns: scoped `localParams` templates, `B` filled from `info` offsets only;
      not declared → `B = zeros(m, 0)`.
- [ ] The error information: `Cep` on the day, one cell per arc (nEpoch×1 for a scalar
      stream; nEpoch×d or nEpoch×d×d, rows epoch-major, for a multi-axis one); the spec
      declares the temporal law and its level only.
- [ ] No unknown field names in `data` (a typo is an error, on purpose — and so are the v2
      count fields and, since v4, `localsPerArc`).

---

### 11. The two in-memory kinds (no local unknowns, diagonal weights)

Not every problem streams. When A fits in memory and there are no arc-wise unknowns, `data`
has no `systems` at all — `fit` tells the kinds apart by the fields present:

**Kind B — matrices**

| field | size | meaning |
|---|---|---|
| `A` | m×u | design matrix |
| `l` | m×1 | observed minus computed |
| `p` | m×1 | weights (diagonal of P), optional; absent = unit weights |

**Kind C — batched callback** (A too large for one matrix, or naturally batched)

| field | size | meaning |
|---|---|---|
| `batchFcn` | fh | `[A, l, p] = batchFcn(k)`, k = 1…nBatch — the `lsqvce.engine.adjust` contract |
| `numBatches` | 1×1 | nBatch |
| `numUnknowns` | 1×1 | u — **linear form only** |
| `x0` | u×1 | optional. With it the problem is **nonlinear**: `[A, w, p] = batchFcn(k, x)` returns the Jacobian at `x` and `w = l − f(x)`; solved by Gauss-Newton (`lsqvce.engine.solveNonlinear`). **u = numel(x0)** — `numUnknowns` must then be absent (`lsqvce:fit:uCount`) |

Rules: `opts.Estimate` must be `"none"` (variance components need the streams of a
streaming system — `lsqvce:fit:modelKind`); `opts.Robust = "huber" | "igg3" | "bisquare"`
switches to M-estimation (`RobustOpts` in Options); shapes are checked before anything runs
(`lsqvce:fit:badMatrix`, `lsqvce:fit:missingField`, `lsqvce:fit:unknownField`). Both kinds
run on the ONE accumulator as a local-free streaming system (batch k = day k), so
`R.residuals()`, `R.contribution()`, `R.plot(..)` and `OnPass` work on them too; the
per-observation diagnostics are in the whitened frame (`R.diag.v = sqrt(p).*v`).

What comes back is the same `Results`: `R.x`, `R.sigma0`, `R.coefficients()`, `R.report()`,
`R.diag`, and for IRLS / Gauss–Newton the pass record on `R.history` (the final IRLS factors:
`R.history(end).factor`) with `R.stopReason`.

```matlab
data = struct('A', A, 'l', l, 'p', p);
R = lsqvce.fit(data);                                   % adjustBlock on the batch adapter
R = lsqvce.fit(data, Robust="huber");                   % adjustRobust (IRLS)
data = struct('batchFcn', @(k, x) myJacobian(k, x), 'numBatches', 20, 'x0', x0);
R = lsqvce.fit(data, MaxPasses=15);                     % adjustNonlinear (Gauss-Newton)
```

Stream names carry no meaning: the axis count of a stream comes from the day it returns
(`Cep` or `nAxis`, §7), never from its name.
