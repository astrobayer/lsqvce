# lsqvce — The data contract: what you give `fit`, and how to write `buildDay`

**Status.** Binding specification, **version 3** (adopted 2026-08-31, Level L10 WP-72; the
v2 contract of 2026-08-31/R4-D4 is superseded — migration table in `RELEASE_NOTES.md`).
`lsqvce.fit` is live for all three data kinds and Models `gls`, `vce`, `pooled`, `ellmer`.
§9 shows the equivalent direct engine calls. This file is the source of the `fit` help text
and of the HTML guide: one authority, no second copy.

Created 2026-08-31 (CEST, UTC+02:00) with assistance from Claude (Anthropic), model Claude
Fable 5. Rewritten for v3 the same day (Claude Fable 5).

---

## 0. The v3 principle: state nothing the data already states

The contract declares **structure**, never **counts**. Every count is derived from the one
place it already lives:

| quantity | derived from |
|---|---|
| nc (common unknowns) | `size(data.commonParams, 1)` |
| nl (local unknowns) | `size(localParams, 1)`, or the `localsPerArc` plan; **neither declared → nl = 0** |
| observation types | `fieldnames(covariance)` |
| u (batched nonlinear) | `numel(x0)` |

The v2 count fields `numCommon`, `numLocal` and `numUnknowns`-with-`x0` are gone, and so is
the duplicate type list `observationTypes`. One source each — nothing left to disagree.

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
opts = lsqvce.Options(Model="vce");      % what to do
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
| `commonParams` | cell, R×C | one row per common unknown. **nc = R.** The column count and meaning are *your* convention (degree, order, epoch, …) — the toolbox reads only the row count and pipes the table to `buildDay` as its 3rd argument. |

The table **is** the declaration: describing the unknowns is your job, counting them is not.
Minimal table when you have nothing to say about them: `num2cell((1:nc).')`.

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
| `localsPerArc` | double or fh | `[]` | the toolbox owns the per-arc numbering, §6.1. **XOR** `localParams`. |
| `localParams` | cell, R×C | `{}` | the local unknowns, one row each (**nl = R**), you own the numbering, §6.2. Passed to `buildDay` as its 4th argument. Neither this nor `localsPerArc` → **nl = 0** (§6.4). |
| `days` | double, 1×D | `data.days` | per-system schedule override. |
| `id` | double 1×1 **or** single char | the index k | system identifier, unique across systems (explicit ids and index defaults share one namespace). Longer labels go in `name`. |
| `name` | char / string | `''` | human-readable label used in reports. |
| `processingOptions` | struct | `struct()` | anything your `buildDay` needs (paths, force models); passed as its 5th argument. |
| `epochCovariance` | function handle | `[]` | `@(day,obsType)` per-epoch D hook for multi-axis streams — the exact separable-descale path of the VCE Shape step, §7. |

**Rule of thumb for beginners:** `buildDay` + `covariance` is a complete system. Add `arcs`
when a day has more than one arc, `localsPerArc` when arcs carry their own unknowns, the
rest when you need it.

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
D = buildDay(day, obsType, common, local, procOpt, arcInfo)       % with local-block plan
```

| argument | what it is |
|---|---|
| `day` | one entry of the day list (numeric) |
| `obsType` | one entry of `fieldnames(covariance)` (char) |
| `common` | `data.commonParams` |
| `local` | `localParams` of this system (or `{}`) |
| `procOpt` | `processingOptions` of this system |
| `arcInfo` | the local-block plan for this day when `localsPerArc` is declared, else `[]` (§6.3) |

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
relabelled, e.g. `lsqvce:fit:noCov`).

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

## 6. Local (arc-wise) unknowns

nl has two declared sources and one default — never more than one applies:

### 6.1 `localsPerArc` — the toolbox owns the numbering

```matlab
s.localsPerArc = 6;                        % six per arc, every arc
s.localsPerArc = @(day, arcId) 6;          % or per arc
```

The toolbox numbers every (day, arc) over the whole adjustment, derives nl, and hands
each day its slice as `arcInfo` (§6.3). Requires `data.days` (or per-system `days`) and a
schedule that is resolvable without reading data (static or resolver).

### 6.2 `localParams` — you own the numbering

An R×C cell, one row per local unknown (**nl = R**, columns your convention), passed as
`buildDay`'s 4th argument. You then place every local column yourself, consistently across
days. Declaring **both** `localsPerArc` and `localParams` is the error `lsqvce:fit:localCount`.

### 6.3 Using `arcInfo`

```matlab
arcInfo.arcId      % 1xA  ids of this day's arcs
arcInfo.offset     % 1xA  column offset of each arc's block in B
arcInfo.count      % 1xA  columns of each block
arcInfo.nl         % 1x1  total local count of the adjustment
arcInfo.offsetOf   % @(arcId) -> offset
```

Fill `B(rows, off + (1:k))` and **never compute an offset yourself**.

### 6.4 A system without local unknowns (nl = 0)

The **default**: declare neither `localsPerArc` nor `localParams` and return
`B = zeros(m, 0)` from `buildDay`. The engine then skips the pre-elimination for that
system; local diagnostics come back as the consistent empty shapes (decided 2026-08-31,
HANDOVER G-07 — in v3 no `numLocal = 0` declaration is needed any more).

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
  across epochs (changing geometry/PDOP). Declare the hook when you run `Model="vce"` on
  multi-axis streams and have D available independently of the day build.

---

## 8. Worked example — two satellites, ragged arcs, per-arc locals

```matlab
% ---- the one function you write --------------------------------------
function D = buildDayGrace(day, obsType, common, ~, procOpt, arcInfo)
    raw = readGraceDay(procOpt.dataDir, day, obsType);   % your reader, deterministic
    D.A = shPartials(raw, common);                       % m x nc
    D.l = raw.observed - raw.computed;                   % O - C
    D.t = raw.tSeconds;
    D.arcIdx = raw.arcRows;                              % cell, one row vector per arc
    D.arcId  = raw.arcIds;                               % this day's arc identities
    D.B = sparse(numel(D.l), arcInfo.nl);
    for a = 1:numel(D.arcId)
        off = arcInfo.offsetOf(D.arcId(a));
        D.B(D.arcIdx{a}, off + (1:arcInfo.count(a))) = localPartials(raw, a);
    end
end

% ---- the declaration ---------------------------------------------------
data.commonParams = shParameterTable(60);                 % nc = 3717 rows, C your columns
data.days = 1:31;
cov.kbr = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',1e-7);
cov.pos = struct('type','diag','sigma',0.02);             % fieldnames = the two types

s = struct();
s.buildDay = @buildDayGrace;
s.covariance = cov;
s.id = 'C';
s.name = 'GRACE-FO C';
s.arcs = @(day, obsType) arcSchedule('C', day);
s.localsPerArc = 6;
s.processingOptions = struct('dataDir', '/data/gracefo');
data.systems(1) = s;

s.id = 'D';  s.name = 'GRACE-FO D';  s.arcs = @(day, obsType) arcSchedule('D', day);
data.systems(2) = s;

% ---- the run -----------------------------------------------------------
opts = lsqvce.Options(Model="vce", GroupBy="satObs", Parallel="parfor");
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
        LocalsPerArc=s.localsPerArc, Days=data.days, ...
        CommonParams=data.commonParams, ...
        ProcessingOptions=s.processingOptions);
end
R = lsqvce.engine.adjustBlock(L, data.days);                        % Model "gls"
R = lsqvce.engine.adjustVCE(L, data.days, 'GroupBy', 'satObs', ...   % Model "vce"
        'Shape', struct(...), 'Scale', struct(...)).run();          % (constructor-only since 2.0)
```

These are ENGINE-tier calls (experts); a user needs none of them. The ENGINE tier keeps the
v2 conveniences (`numCommon`/`numLocal` counts, `'fromDay'`, `cache`, the plot hooks, raw
4-arg loader closures) — they left the USER surface, not the toolbox.

Field names map 1:1 to `makeLoader`'s Name-Value names (lowerCamelCase field ↔
UpperCamelCase name); `fit` performs exactly this mapping internally, plus the three
derivations of §0.

---

## 10. Checklist before the first run

- [ ] `buildDay` returns `A`, `B`, `l`, `t` with one row count; `l` is O − C; nothing is whitened.
- [ ] Calling `buildDay` twice with the same arguments gives identical output.
- [ ] `t` is in seconds; the sampling matches the covariance spec you declared.
- [ ] `data.commonParams` has one row per common unknown — `size(A, 2)` must equal it.
- [ ] Every observation type is a field of `covariance` — the fieldnames are the type list.
- [ ] Ragged schedule → resolver + `arcId`; static schedule → count/vector, no `arcId`.
- [ ] Local unknowns: `localsPerArc` (fill `B` from `arcInfo`) or `localParams` (you number);
      neither → `B = zeros(m, 0)`.
- [ ] Multi-axis streams: epoch-major rows, `nAxis`, `Cep` on the day; the `epochCovariance`
      hook when Model="vce" should take the exact Shape path.
- [ ] No unknown field names in `data` (a typo is an error, on purpose — and so is every v2
      count field).

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

Rules: `opts.Model` must be `"gls"` (variance components need the groups of a streaming
system — `lsqvce:fit:modelKind`); `opts.Robust = "irls"` switches to M-estimation
(`Influence`, `Tuning`, `RobustMaxIter`, `RobustTol` in Options); shapes are checked before
anything runs (`lsqvce:fit:badMatrix`, `lsqvce:fit:missingField`, `lsqvce:fit:unknownField`).

What comes back is the same `Results`: `R.x`, `R.sigma0`, `R.coefficients()`, `R.report()`,
and the per-observation diagnostics in `R.diag` (residuals `R.diag.v`). Two methods need the
streaming engine and say so: `R.residuals()` and `R.contribution()` refuse with
`lsqvce:Results:noAccumulator`. The driver's own record (verdict, IRLS weights, or the
nonlinear iteration log) is on `R.source`.

```matlab
data = struct('A', A, 'l', l, 'p', p);
R = lsqvce.fit(data);                                   % = lsqvce.engine.adjust, one batch
R = lsqvce.fit(data, lsqvce.Options(Robust="irls"));    % = lsqvce.engine.robust
data = struct('batchFcn', @(k, x) myJacobian(k, x), 'numBatches', 20, 'x0', x0);
R = lsqvce.fit(data, lsqvce.Options(NonlinearMaxIter=15));   % = lsqvce.engine.solveNonlinear
```

**One sentence on stream names (from the L2 CI catch):** a stream whose observation type is
`pos` is treated as **3 components per epoch** by default (`Options.ComponentsPerEpoch`), which
selects the separable model and needs a per-epoch scaling (§7). A scalar stream that happens
to be called `pos` must declare `ComponentsPerEpoch = struct('pos', 1)`.
