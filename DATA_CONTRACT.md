# lsqvce — The data contract: what you give `fit`, and how to write `buildDay`

**Status.** Binding specification adopted 2026-08-31 (R4 decision R4-D4, HANDOVER §11.9).
`lsqvce.fit` is live for the streaming data path (Level L2, WP-56; Models `gls`, `vce`,
`pooled`); the in-memory path `fit(A, l, …)`, `Robust`, and Model `ellmer` follow in WP-56b.
§9 shows the equivalent direct engine calls. This file is the source of the `fit` help text
and of the HTML guide (WP-56, WP-60): one authority, no second copy.

Created 2026-08-31 (CEST, UTC+02:00) with assistance from Claude (Anthropic), model Claude
Fable 5.

---

## 1. The mental model in one paragraph

You describe **one adjustment problem** with a plain MATLAB struct called `data`. It answers
three questions: *what is estimated for everybody* (`data.common`), *which days are
processed* (`data.days`), and *which observing systems deliver data* (`data.systems`, one
entry per satellite or instrument). Each system points to **one function you write**,
`buildDay`, that returns the design matrices and observations of a single day. Everything
else — reading days in the right order, cutting them into arcs, whitening with your
covariance model, iterating the variance components — is the toolbox's job:

```matlab
opts = lsqvce.Options(Model="vce");      % what to do (Level L1)
R    = lsqvce.fit(data, opts);           % the one verb (Level L2)
R.report                                 % the verdict
```

You never call `makeLoader` yourself; `fit` builds the loaders from `data`. `R.systems`
(the existing Results field) records what ran: per system its id, nl, the loader closure,
the day list and the name.

---

## 2. The `data` struct — field by field

`data` is a **scalar struct**. Field names are **case-sensitive** and use lowerCamelCase.
Unknown field names are an **error** (`lsqvce:fit:unknownField`), which is deliberate: a
typo in a struct field would otherwise be silently ignored.

### 2.1 `data.common` — the shared unknowns (nc stated exactly once)

Exactly **one** of the two fields:

| field | type / size | meaning |
|---|---|---|
| `commonParams` | cell, R×14 | one row per common unknown; the 14 columns follow *your* convention (degree, order, epoch, …). nc = R. Passed to `buildDay` as its 3rd argument. |
| `numCommon` | double, 1×1 | just the count nc, when you do not want a parameter table. |

*Why here and not per system?* The common block is shared by every system; stating nc in
two places would let the two disagree, and that disagreement was a real bug class.

### 2.2 `data.days` — the default schedule

| field | type / size | meaning |
|---|---|---|
| `days` | double, 1×D | day identifiers, e.g. `1:30` or `datenum` values cast to integers. Every system uses this list unless it carries its own `days`. Required unless every system has one. |

Days must be **numeric integers**: the loader cache is keyed on the number.

### 2.3 `data.systems(k)` — one entry per observing system

`data.systems` is a **1×M struct array** (M ≥ 1). All entries have the same fields; fill an
unused optional field with `[]`.

**Required**

| field | type / size | meaning |
|---|---|---|
| `buildDay` | function handle | your day builder, §3. |
| `observationTypes` | cellstr / string, 1×T | the observation types this system delivers, e.g. `{'kbr'}` or `{'pos','kbr'}`. `buildDay` is called once per (day, type). |
| `covariance` | struct | one field per observation type holding a covariance spec (§5), e.g. `cov.kbr = struct('type','diag','sigma',0.2)`. |
| `systemId` | double 1×1 **or** a single char | unique across systems, e.g. `1` and `2`, or `'C'` and `'D'`. Longer labels go in `name`. |
| one of `localsPerArc` \| `numLocal` \| `localParams` | see §6 | how many local (arc-wise) unknowns nl this system has. **0 is valid** (§6.4). |

**Optional**

| field | type / size | default | meaning |
|---|---|---|---|
| `arcs` | double, fh, or `'fromDay'` | `1` | which arcs each day carries, §4. |
| `days` | double, 1×D | `data.days` | per-system schedule override. |
| `name` | char / string | `''` | human-readable label used in reports. |
| `cache` | logical | `true` | keep the last day in memory (needed for the arc reads of that day; switch off only for memory emergencies). |
| `epochCovariance` | function handle | `[]` | `@(day,obsType)` returning the per-epoch scaling of a multi-axis stream, §7. |
| `processingOptions` | struct | `struct()` | anything your `buildDay` needs (paths, force models); passed as its 5th argument. |
| `plotPostSolve`, `plotPostModel` | fh or containers.Map | `[]` | plot hooks fired after the solve / after the noise model (advanced). |

**Rule of thumb for beginners:** start with `buildDay`, `observationTypes`, `covariance`,
`systemId`, `localsPerArc`, and `arcs` as a resolver. Add the rest when you need it.

### 2.4 A pitfall you will hit once: `struct()` with cell values

```matlab
% WRONG: struct() expands cell values into a struct ARRAY
s = struct('observationTypes', {'pos','kbr'});      % s is 1x2
% RIGHT: double braces, or assign field by field
s = struct('observationTypes', {{'pos','kbr'}});    % s is 1x1
s.observationTypes = {'pos','kbr'};                 % clearest
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
| `obsType` | one entry of `observationTypes` (char) |
| `common` | `data.common.commonParams` (or `{}`) |
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
| `Cep` | nEpoch×d or nEpoch×d×d | per-epoch scaling of a multi-axis stream for this day (§7). |

### 3.5 What is checked, and the error you get

| error id | meaning |
|---|---|
| `lsqvce:makeLoader:buildDay` | your `buildDay` threw; the original error is attached as the cause |
| `lsqvce:makeLoader:badDay` | a required field (`A`, `B`, `l`, `t`) is missing |
| `lsqvce:makeLoader:dayShape` | `A`, `B`, `l`, `t` do not share a row count |
| `lsqvce:accum:arcBadNc` / `arcBadNl` | `A` or `B` has the wrong number of columns |
| `lsqvce:makeLoader:badArc` | a positional arc id is out of range (static schedule declared wrongly) |
| `lsqvce:makeLoader:badArcId` | `arcId` does not name every `arcIdx` entry exactly once |

These are raised **while reading**, by the engine, and keep their engine ids. Errors in the
*declaration* itself are raised by `fit` before any data is read and carry `lsqvce:fit:` ids
(`unknownField`, `noSystems`, `missingField`, `localCount`, `commonCount`, `noDays`,
`duplicateSystem`, `builderKind`, and every `makeLoader` construction id relabelled, e.g.
`lsqvce:fit:noCov`).

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

`s.arcs = 'fromDay'` is the lazy alternative: the day is read during the schedule probe and
`arcId` (or `1:numel(arcIdx)`) is taken from it. It works, but it reads data before the
run starts and cannot be combined with `localsPerArc`.

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
specs at run time; you only give the starting point. Full reference: `doc lsqvce.covModel`.

---

## 6. Local (arc-wise) unknowns

### 6.1 `localsPerArc` — recommended

```matlab
s.localsPerArc = 6;                        % six per arc, every arc
s.localsPerArc = @(day, arcId) 6;          % or per arc
```

The toolbox numbers every (day, arc) over the whole adjustment, derives nl, and hands
each day its slice as `arcInfo` (§6.3). Requires `data.days` (or per-system `days`) and a
schedule that is resolvable without reading data (static or resolver — not `'fromDay'`).

### 6.2 `numLocal` / `localParams` — you own the numbering

`numLocal = nl` (a count) or `localParams` (R×14 cell, nl = R, passed as `buildDay`'s 4th
argument). You then place every local column yourself, consistently across days.

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

Valid (decided 2026-08-31, HANDOVER G-07). Declare `numLocal = 0` (or `localsPerArc = 0`)
and return `B = zeros(m, 0)` from `buildDay`. The engine then skips the pre-elimination for
that system; local diagnostics come back as the consistent empty shapes.

---

## 7. Multi-axis streams (positions)

A kinematic position has three components per epoch that share one temporal law. Return
the rows **epoch-major** (x₁ y₁ z₁ x₂ y₂ z₂ …), set `nAxis = 3`, and provide the per-epoch
scaling `D` in one of three ways (first found wins): `Cep` on the covariance spec, `Cep` on
the returned day, or the system's `epochCovariance` hook. A multi-axis stream with no `D`
is an error (`noD`).

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
data.common.commonParams = shParameterTable(60);          % nc = 3717 rows
data.days = 1:31;
cov.kbr = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',1e-7);
cov.pos = struct('type','diag','sigma',0.02);

s = struct();
s.buildDay = @buildDayGrace;
s.observationTypes = {'kbr','pos'};
s.covariance = cov;
s.systemId = 'C';
s.name = 'GRACE-FO C';
s.arcs = @(day, obsType) arcSchedule('C', day);
s.localsPerArc = 6;
s.processingOptions = struct('dataDir', '/data/gracefo');
data.systems(1) = s;

s.systemId = 'D';  s.name = 'GRACE-FO D';  s.arcs = @(day, obsType) arcSchedule('D', day);
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
    L{k} = lsqvce.makeLoader(s.buildDay, ...
        ObservationTypes=s.observationTypes, Covariance=s.covariance, ...
        SystemId=s.systemId, Name=s.name, Arcs=s.arcs, ...
        LocalsPerArc=s.localsPerArc, Days=data.days, ...
        CommonParams=data.common.commonParams, ...
        ProcessingOptions=s.processingOptions);
end
R = lsqvce.adjustBlock(L, data.days);                       % one GLS solve
R = lsqvce.adjustVCE(L, data.days).groupComponentsBy('satObs').run();   % VCE loop
```

Field names map 1:1 to `makeLoader`'s Name-Value names (lowerCamelCase field ↔
UpperCamelCase name); `fit` performs exactly this mapping internally.

---

## 10. Checklist before the first run

- [ ] `buildDay` returns `A`, `B`, `l`, `t` with one row count; `l` is O − C; nothing is whitened.
- [ ] Calling `buildDay` twice with the same arguments gives identical output.
- [ ] `t` is in seconds; the sampling matches the covariance spec you declared.
- [ ] Every `systemId` is unique; nc is stated once in `data.common`.
- [ ] Ragged schedule → resolver + `arcId`; static schedule → count/vector, no `arcId`.
- [ ] Local unknowns via `localsPerArc` (and `B` filled from `arcInfo`); nl = 0 is fine with `B = zeros(m, 0)`.
- [ ] Multi-axis streams: epoch-major rows, `nAxis`, and a source for `D`.
- [ ] No unknown field names in `data` (a typo is an error, on purpose).
