# lsqvce — The data contract: what you give `fit`, and how to write `buildDay`

**Status.** Binding specification, **version 5** (adopted 2026-09-03, Review R5 level M5:
the v4 SHAPES are unchanged — every v4 declaration still fits — and three things were
ADDED, listed in §0a below; option names are those of `lsqvce.Options` v5, the v4 → v5
map is in `RELEASE_NOTES.md`).
`lsqvce.fit` is live for all three data kinds, `Estimate = "none" | "noise"`, `Robust`, the
nonlinear form and `OnPass`. §9 shows the equivalent direct engine calls. This file is the
source of the `fit` help text and of the HTML guide: one authority, no second copy.

Created 2026-08-31 (CEST, UTC+02:00) with assistance from Claude (Anthropic), model Claude
Fable 5. Rewritten for v3 and again for v4 (scoped tables) the same day (Claude Fable 5).

---

## 0a. The v5 delta (R5 M5): three additions, no shape changes

1. **`data.priors(k) = struct('L', pxnc, 'x0', px1, 'sigma', s, 'name', .., 'estimate', "none"|"scale")`** — declared
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
   `R.noise.priorHistory` record it. R5-D8 round 1 adds `from` (the pass at which the group
   becomes active, default 1; a single-pass run reports `lsqvce:fit:priorNotActive` and skips
   it) and `group` (a label: priors sharing it are merged into ONE group -- rows divided by
   their declared sigmas, the group sigma a common factor -- so `estimate="scale"` estimates
   one component for the members).
1b. **`data.systems(k).priors(l) = struct('L', 'x0', 'sigma', 'from', 'name', 'estimate')`** (R5-D8
   round 5) — the priors on system k's LOCAL parameters, the literal counterpart of `data.priors(j)`:
   one struct per GROUP, `L` (q × nl) in the columns of that system's B as **`lsqvce.layout(data)`**
   announces them before fit (`lay.systems(k).col(row, day, arc)`), `x0`, one `sigma` per group,
   `from`, `name`, `estimate = "none" | "scale"`. Every row must lie inside one day's local block
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
   estimated spectral shape; every stream keeps its own level.
3. **The per-pass inspection hook** — `Options.OnPass = @(ctx)` receives, after every pass
   of every iterative run, the pass state, the live accumulator, `ctx.data` (this
   declaration as given) and `ctx.processingOptions{k}` (system k's `processingOptions` as
   declared); return `false` to stop. Nothing else about `buildDay` or the tables changed.

Semantics made explicit in v5: a full SPD weight matrix returned by a batch source is
whitened exactly (`chol(P)`), and per-observation diagnostics of a batch problem are
reported in the WHITENED frame (`R.diag.v = sqrt(p).*v`; only their sums are
frame-invariant). Multi-axis streams must be declared `separable` (a `diag` law ignores the
per-epoch covariance).

## 0. The v4 principle: state nothing the data already states — and scope, don't number

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

## 1. The mental model in one paragraph

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

## 2. The `data` struct — field by field

`data` is a **scalar struct**. Field names are **case-sensitive** and use lowerCamelCase.
Unknown field names are an **error** (`lsqvce:fit:unknownField`), which is deliberate: a
typo in a struct field would otherwise be silently ignored. The retired v2 names (`common`,
`numCommon`, `numLocal`, `observationTypes`, `systemId`, `cache`, `loader`, the plot hooks)
are unknown fields now — the error names the field, so a v2 declaration fails loudly, not
subtly.

### 2.1 `data.commonParams` — the shared unknowns (required)

| field | type / size | meaning |
|---|---|---|
| `commonParams` | cell, R×C (C ≥ 1) | one row per common **template**. **Column 1 is the scope token**: `'common'` (one column of A, always) or `'day'` (one column per calendar day that carries at least one scheduled arc in at least one system). Columns 2…C are *your* convention (degree, order, epoch, …). nc is **derived** from the expansion; the table is piped to `buildDay` as its 3rd argument. `'arc'` is forbidden here — arc sets are system-specific (`lsqvce:fit:badScope`). |

**Binding column layout of A:** the `'common'` rows first (table order), then one day
block per calendar day (calendar order — the sorted union of every declared day list).
`buildDay` receives this day's block position as `info.commonDayOffset` /
`info.commonDayCount` (§6.3). Minimal table when you have nothing to say:
`[repmat({'common'}, nc, 1), num2cell((1:nc).')]`.

### 2.2 `data.days` — the default schedule

| field | type / size | meaning |
|---|---|---|
| `days` | double, 1×D | day identifiers, e.g. `1:30` or `datenum` values cast to integers. Every system uses this list unless it carries its own `days`. Required unless every system has one. |

Days must be **numeric integers**: the loader cache is keyed on the number.

### 2.3 `data.systems(k)` — one entry per observing system

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
| `arcs` | double or fh | `1` | which arcs each day carries, §4. (`'fromDay'` is ENGINE-tier only.) |
| `localParams` | cell, R×C (C ≥ 1) | `{}` | the local **templates**, §6: column 1 is the scope token `'arc'` \| `'day'` \| `'system'`, columns 2…C your convention; nl is **derived** from the expansion. Passed to `buildDay` as its 4th argument. Not declared → **nl = 0** (§6.4). |
| `days` | double, 1×D | `data.days` | per-system schedule override. |
| `id` | double 1×1 **or** single char | the index k | system identifier, unique across systems (explicit ids and index defaults share one namespace). Longer labels go in `name`. |
| `name` | char / string | `''` | human-readable label used in reports. |
| `processingOptions` | struct | `struct()` | anything your `buildDay` needs (paths, force models); passed as its 5th argument. |
| `epochCovariance` | function handle | `[]` | `@(day,obsType)` per-epoch D hook for multi-axis streams — the exact separable-descale path of the VCE Shape step, §7. |

**Rule of thumb for beginners:** `buildDay` + `covariance` is a complete system. Add `arcs`
when a day has more than one arc, `'arc'` template rows in `localParams` when arcs carry
their own unknowns, the rest when you need it.

### 2.4 A pitfall you will hit once: `struct()` with cell values

```matlab
% WRONG: struct() expands cell values into a struct ARRAY
s = struct('localParams', {'bias','drift'});          % s is 1x2
% RIGHT: double braces, or assign field by field
s = struct('localParams', {{'bias';'drift'}});        % s is 1x1
s.localParams = {'bias'; 'drift'};                    % clearest
```

---

## 3. The `buildDay` contract

### 3.1 Who calls it, and how often

`fit` calls `buildDay` **once per (day, observationType)** and cuts the returned day into
arcs itself. Every pass over the data reads each day **twice** (once to accumulate the
normal equations, once for the residuals), and the variance-component loop multiplies that
by the number of passes. Two consequences:

1. **Determinism is mandatory.** Same arguments → byte-identical output. Random numbers,
   clocks, mutable globals are forbidden inside `buildDay`. (The residuals are formed as
   `v = l − A·x − B·y` from the second read; any difference lands in σ₀.)
2. **Cache expensive reads yourself.** The toolbox caches *one built day*; if the raw file
   read is the expensive part and serves several observation types, keep your own cache
   in the closure that binds `buildDay`.

### 3.2 Signature — declare only what you need

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

### 3.3 Return value — required fields

| field | size | meaning |
|---|---|---|
| `A` | m×nc | design matrix of the **common** unknowns, all observations of the day |
| `B` | m×nl | design matrix of the **local** unknowns; may be `sparse` (recommended when arc-wise) |
| `l` | m×1 | **observed minus computed** (the reduced observation). Not the raw observation. |
| `t` | m×1 | epoch of each row in **seconds** (the whitener needs the lags) |

Return **raw, unwhitened** matrices. Whitening happens per arc inside the toolbox from your
covariance spec, and the VCE loop reweights there too — a builder that scales its own rows
would be reweighted twice.

### 3.4 Return value — optional fields

| field | size | meaning |
|---|---|---|
| `arcIdx` | cell, 1×nArc | row-index vector of each arc (lengths may differ). Absent → the whole day is one arc. |
| `arcId` | double, 1×nArc | the **identity** of each arc, when the schedule is ragged (§4.2). Must name every `arcIdx` entry exactly once with unique positive integers. |
| `nAxis` | 1×1 | for multi-axis streams (positions: 3), rows are epoch-major, m = nEpoch·nAxis (§7). |
| `Cep` | nEpoch×d or nEpoch×d×d | per-epoch scaling of a multi-axis stream for this day (§7). **This is the standard D source.** |

### 3.5 What is checked, and the error you get

| error id | meaning |
|---|---|
| `lsqvce:makeLoader:buildDay` | your `buildDay` threw; the original error is attached as the cause |
| `lsqvce:makeLoader:badDay` | a required field (`A`, `B`, `l`, `t`) is missing |
| `lsqvce:makeLoader:dayShape` | `A`, `B`, `l`, `t` do not share a row count |
| `lsqvce:accum:arcBadNc` / `arcBadNl` | `A` or `B` has the wrong number of columns — the derived nc/nl are enforced against **every** day |
| `lsqvce:makeLoader:badArc` | a positional arc id is out of range (static schedule declared wrongly) |
| `lsqvce:makeLoader:badArcId` | `arcId` does not name every `arcIdx` entry exactly once |

These are raised **while reading**, by the engine, and keep their engine ids. Errors in the
*declaration* itself are raised by `fit` before any data is read and carry `lsqvce:fit:` ids
(`unknownField`, `noSystems`, `noCommon`, `missingField`, `localCount`, `noDays`, `badId`,
`duplicateSystem`, `builderKind`, `badArcs`, and every `makeLoader` construction id
relabelled, e.g. `lsqvce:fit:paramType`).

---

## 4. Arc schedules

### 4.1 Static — every day looks the same

```matlab
s.arcs = 4;            % arcs 1..4 on every day, every observation type
s.arcs = [1 2 4];      % these ids on every day
s.arcs.kbr = 2;  s.arcs.pos = 4;    % per observation type
```

Arcs are then addressed **by position** in `arcIdx`; an id beyond `numel(arcIdx)` is an
error, which catches a mis-declared schedule.

### 4.2 Ragged — days differ (gaps, manoeuvres, short months)

Declare a **resolver** that returns the arc ids of a day *without reading data*:

```matlab
s.arcs = @(day, obsType) myArcTable{day};      % 1xA vector of ids for that day
```

and let `buildDay` return `arcId` naming the arcs it actually carries. Arcs are then
addressed **by id**; an id that a day does not carry is treated as a **data gap** and
skipped — no padding, no empty arcs.

(The v2 string form `arcs = 'fromDay'` is not part of this contract: it reads data before
the run starts and cannot be combined with `localsPerArc`. It survives as an ENGINE-tier
convenience on `lsqvce.engine.makeLoader`; `fit` rejects it with `lsqvce:fit:badArcs`.)

---

## 5. Covariance specs — the `covariance.(obsType)` structs

One struct per observation type; `type` selects the model. The most common:

| `type` | fields | use for |
|---|---|---|
| `'diag'` | `sigma` (scalar or m×1 std) | uncorrelated noise, e.g. `struct('type','diag','sigma',0.2)` |
| `'arma'` | `ar` `[1 a1 … ap]`, `ma` (default 1), `sigma` (innovation std) | coloured noise with a known filter |
| `'acf'` | `acf` `@(dt,opt)`, `sigma` | correlation given as a function of the lag |
| `'psd'` | `f` (nf×1, Hz), `S` (nf×1) | noise given as a one-sided PSD |
| `'var'` | `ar` (d×d×p), `sigma` (d×d) | vector-AR for multi-axis streams |

`sigma` is always a **standard deviation**, never a variance. The VCE loop scales these
specs at run time; you only give the starting point. Full reference: `doc lsqvce.engine.covModel`.

---

## 6. Local unknowns — the scoped `localParams` table

One mechanism: the table's rows are templates, column 1 says **how often** each row
instantiates.

### 6.1 The three scopes

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
| `'day'` | day whose schedule carries ≥ 1 arc | ACC scale/bias/drift |
| `'system'` | whole adjustment | instrument constants; also the v3 free table (see migration) |

Requires `days` and a schedule resolvable **without reading data** (static or resolver)
whenever `'arc'`/`'day'` rows are present; `'system'`-only tables need neither. nl is
derived; declaring `numLocal` too is the conflict error.

**Binding column layout of B (per system):** per day (declared order) — the arc blocks in
**sorted arc-id order** (each block = the `'arc'` rows, table order), then the `'day'`
block — and after all days one trailing `'system'` block. A day whose schedule has no
arcs gets **no blocks at all** (§6.4).

### 6.2 Heterogeneous per-arc counts stay ENGINE-tier

A template has the same rows for every arc. If arcs genuinely carry *different* counts,
that is the ENGINE-tier `lsqvce.engine.makeLoader(..., 'LocalsPerArc', @(day, arcId) n)`
(the WP-52 plan) — deliberately not on the `fit` surface.

### 6.3 Using `info` (buildDay's 6th argument)

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

### 6.4 Days without arcs, systems without locals

A calendar day on which **no system** schedules an arc gets no common day block (nc
shrinks); a system day without arcs gets no arc and no `'day'` blocks (nl shrinks) —
offsets stay gap-safe on ragged schedules and manoeuvre days. A system without local
unknowns is still the **default**: declare no `localParams` and return `B = zeros(m, 0)`
(the engine skips the pre-elimination; empty diagnostics per HANDOVER G-07).

---

## 7. Multi-axis streams (positions)

A kinematic position has three components per epoch that share one temporal law. Return
the rows **epoch-major** (x₁ y₁ z₁ x₂ y₂ z₂ …), set `nAxis = 3`, and provide the per-epoch
scaling `D`:

- **Standard: `Cep` on the returned day** (nEpoch×d or nEpoch×d×d) — you know D when you
  build the day. Alternative for a *static* scaling: `Cep` on the covariance spec. A
  multi-axis stream with no `D` is an error (`noD`).
- **`epochCovariance` hook (optional, VCE only).** Whitening never needs it — but the
  **Shape step** of the VCE iteration (re-estimating the temporal AR law from residuals)
  does not read the day's `Cep`: with the hook it takes the *exact separable de-scale* path
  (`lsqvce.kernel.residualAcv` with the full per-epoch D, fetched per day **without**
  rebuilding the day); without it, it falls back to descaling by the stream's current
  scalar σ — an approximation that biases the estimated AR law when D varies strongly
  across epochs (changing geometry/PDOP). Declare the hook when you run `Estimate="noise"` on
  multi-axis streams and have D available independently of the day build.

---

## 8. Worked example — the GRACE parameterization: three scopes, one table

Arc-wise initial conditions, **day-wise** accelerometer scale/bias/drift, per-arc KBR
empiricals — the case v3 could only hand-number, and the driver of v4:

```matlab
% ---- the one function you write --------------------------------------
function D = buildDayGrace(day, obsType, common, ~, procOpt, info)
    raw = readGraceDay(procOpt.dataDir, day, obsType);   % your reader, deterministic
    D.A = shPartials(raw, common);                       % m x info.nc; the 'common' block
    D.l = raw.observed - raw.computed;                   % O - C
    D.t = raw.tSeconds;
    D.arcIdx = raw.arcRows;                              % cell, one row vector per arc
    D.arcId  = raw.arcIds;                               % this day's arc identities
    D.B = sparse(numel(D.l), info.nl);
    for a = 1:numel(D.arcId)
        off = info.offsetOf(D.arcId(a));
        D.B(D.arcIdx{a}, off + (1:info.count(a))) = icAndKbrPartials(raw, a);
    end
    if info.dayCount > 0                                 % this day's ACC block
        D.B(:, info.dayOffset + (1:info.dayCount)) = accPartials(raw);
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
s.arcs = @(day, obsType) arcSchedule('C', day);
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

## 9. What `fit` does internally — the equivalent direct engine calls

```matlab
L = cell(1, numel(data.systems));
for k = 1:numel(data.systems)
    s = data.systems(k);
    L{k} = lsqvce.engine.makeLoader(s.buildDay, ...
        ObservationTypes=fieldnames(s.covariance).', Covariance=s.covariance, ...
        SystemId=s.id, Name=s.name, Arcs=s.arcs, ...
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
`'fromDay'`, `cache`, the plot hooks, raw 4-arg loader closures) — they left the USER
surface, not the toolbox.

Field names map 1:1 to `makeLoader`'s Name-Value names (lowerCamelCase field ↔
UpperCamelCase name); `fit` performs exactly this mapping internally, plus the three
derivations of §0.

---

## 10. Checklist before the first run

- [ ] `buildDay` returns `A`, `B`, `l`, `t` with one row count; `l` is O − C; nothing is whitened.
- [ ] Calling `buildDay` twice with the same arguments gives identical output.
- [ ] `t` is in seconds; the sampling matches the covariance spec you declared.
- [ ] `data.commonParams` rows carry the scope token in column 1; `size(A, 2)` must equal
      `info.nc` (use `zeros(m, info.nc)` and place blocks by the `info` offsets).
- [ ] Every observation type is a field of `covariance` — the fieldnames are the type list.
- [ ] Ragged schedule → resolver + `arcId`; static schedule → count/vector, no `arcId`.
- [ ] Local unknowns: scoped `localParams` templates, `B` filled from `info` offsets only;
      not declared → `B = zeros(m, 0)`.
- [ ] Multi-axis streams: epoch-major rows, `nAxis`, `Cep` on the day; the `epochCovariance`
      hook when Estimate="noise" should descale per epoch (the separable law).
- [ ] No unknown field names in `data` (a typo is an error, on purpose — and so are the v2
      count fields and, since v4, `localsPerArc`).

---

## 11. The two in-memory kinds (no local unknowns, diagonal weights)

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

**One sentence on stream names (from the L2 CI catch):** a stream whose observation type is
`pos` is treated as **3 components per epoch** by default (`Options.ComponentsPerEpoch`), which
selects the separable model and needs a per-epoch scaling (§7). A scalar stream that happens
to be called `pos` must declare `ComponentsPerEpoch = struct('pos', 1)`.
