# GRAVFIRE → `lsqvce.fit`: preparing the data, step by step

*R5 M0 close-out, 2026-09-02. Companion of `template/gravfire_data.m`, which
implements every step below and was verified on `testdata.mat` (2 days,
degree 10) on R2026a: the declaration passes every `lsqvce:fit` check and a
dry run with random design matrices accumulates both days in 2 s. Read this
guide with that file open; the STEP numbers match its comments. Created
2026-09-02 (CEST, UTC+02:00) with assistance from Claude (Anthropic), model
Claude Fable 5.1.*

---

## What you are building

One MATLAB struct, `data`, that `lsqvce.fit` reads. It has three top-level
fields and nothing else:

```
data.commonParams   cell R×15   the shared unknowns (your 117 Stokes rows)
data.days           double 1×D  the calendar (MJD integers)
data.systems        struct 1×1  ONE system: GraceA + GraceB + KBR together
```

`fit` **derives** every count from this (nc = 117, nl = 144, arcs per day,
observation types) and enforces them against every day it reads. You never
type a count.

`data.systems` has these fields:

| field | what it holds | in the template |
|---|---|---|
| `id` | `1` (numeric, or ONE character — not a word) | STEP 7 |
| `name` | `'GRACE'` (any label; used in reports) | STEP 7 |
| `buildDay` | your day builder, see Part B | STEP 7 |
| `covariance` | a struct whose **field names are the observation types** | STEP 3 |
| `arcs` | a function `@(day, obsType)` returning the piece ids of that day | STEP 5 |
| `localParams` | one day's local templates, 72 rows | STEP 4 |
| `processingOptions` | anything `buildDay` needs (your bundles, the design function) | STEP 6 |

---

## Part A — the declaration (what `gravfire_data.m` does today; what GRAVFIRE should emit tomorrow)

### STEP 1 — the common table: put the scope token in column 1

Your `commonData.parameter` is 117 × 14 (`'c'`/`'s'`, degree, order, `'even'`, …,
a-priori sigma 1e-9 in column 13). Prepend one column with the word `'common'`:

```matlab
data.commonParams = [repmat({'common'}, 117, 1), commonData.parameter];   % 117 x 15
```

(In the 15-column table -- token in column 1 -- your column 13 (the column scale) is column 14 and your column 11 (regularisation) is column 12; see Part C and STEP 6.)

`'common'` means "one column of A for the whole adjustment". Columns 2…15 are
yours; `fit` never reads them, it hands the table to `buildDay` as its third
argument. **In GRAVFIRE:** emit the 15-column table directly; keep column 13 (the
column scale) and column 11 (regularisation) in it -- STEP 6 and Part C say what happens to them.

### STEP 6 — column 11 (the constraining factor) becomes `data.priors` (R5-D8 round 1)

`{'Tikhonov', k, sigma}` on row j means: constrain parameter j to its a-priori value
(column 10; `[]` = 0) with standard deviation `sigma`, from pass `k` on. The template
translates every such COMMON row into `data.priors(i) = struct('L', e_j', 'x0', ..,
'sigma', sigma, 'from', k, 'name', ..)`. An optional 4th entry `{'Tikhonov', k, sigma,
'label'}` is a GROUP label: all rows with the same label share ONE variance component
estimated by VCE (`estimate="scale"`, Koch-Kusche; a single parameter cannot estimate its
own sigma -- its redundancy is at most 1 -- a labelled group can). `'none'` or `[]` = no
prior. In your bundle every common row is `[]` and the constraints sit on LOCAL rows
(`*accscale`): those become `data.systems.priors` in STEP 7 (R5-D8 round 5) -- one struct per group, L in the
columns of B as `lsqvce.layout` announces them; the same mechanism and the same Koch-Kusche
estimation as `data.priors` on the common parameters.

### STEP 2 — the calendar

```matlab
data.days = double(myOptions.days(:).');       % [59396 59397]
```

Integers (MJD is fine). Every day listed must be a day the bundle has.

### STEP 3 — observation types become FIELD NAMES

Today the three groups are addressed by number (`observations{day, g}`,
`variance(g)`). `fit` addresses them by **name**: the fieldnames of
`covariance`, in that order. Choose short valid identifiers:

| group | `name` / `observationQuantity` | field name |
|---|---|---|
| 1 | GraceA / kinematicPosition | `posA` |
| 2 | GraceB / kinematicPosition | `posB` |
| 3 | GraceX / range-rate | `kbr` |

and declare the **starting** temporal law of each as white, with the level as
a standard deviation (not the variance):

```matlab
cov.posA = struct('type','diag','sigma', sqrt(sd.variance(1)));   % 0.1 m
cov.posB = struct('type','diag','sigma', sqrt(sd.variance(2)));   % 0.1 m
cov.kbr  = struct('type','diag','sigma', sqrt(sd.variance(3)));   % 1e-6 m/s
```

The VCE loop estimates the temporal law (`fftmodel` in your options) **from
these starting points**; nothing about the FFT model is declared here, and
nothing about the per-epoch accuracy either — that is the day's `Cep`
(STEP 10), which the toolbox composes with the law (v5.3; the `separable`
type and its `d` are gone).
**In GRAVFIRE:** store the bundle as `bundles.posA`, `.posB`, `.kbr` (or keep
the `days × groups` cell and a name→column map, as the template does). Never
rely on group order alone.

### STEP 4 — the local table: ONE day's templates, all `'day'`-scoped

Your `systemsData{1}.parameter` is the **expanded** table: 144 rows = 2 days ×
2 satellites × 36 unknowns, with absolute times `[MJD sec]` and the day in the
flags. `fit` wants the opposite: the **36 rows per satellite for one day**,
and it expands them over the calendar itself. So:

1. keep only the rows of the first day (72 rows);
2. column 1 = `'day'` for every row ("one instance per calendar day");
3. absolute spans `[59396 21600]` → seconds of day `21600` (columns 3, 4 after
   the token; the spline windows are the same every day);
4. the knot table `[MJD.frac value]` → `[secondsOfDay value]`;
5. drop the day from the flag vector `[59396 1 1 1 0]` → `[1 1 1 0]`.

The template's `templateFromExpanded` does exactly this; **in GRAVFIRE**, emit
the 72-row template instead of the 144-row expansion and the helper disappears.
The result (`data.systems.localParams`, 72 × 15) looks like:

```
{'day'}  {'GraceA_x0'}          {0}      {86390}  {[]}          {8640}  {[1 1 1 1]} ...
{'day'}  {'GraceA_xaccspline'}  {21600}  {86399}  {2161x2}      {6480}  {[1 1 1 0]} ...
```

Why no `'arc'` rows: your initial conditions span the whole day (`numarc = 1`),
so nothing is per piece. If one day you integrate per piece, those rows
become `'arc'` — one line changes, nothing else.

**How the `'day'` mechanism works (from `makeLoader.buildScopePlan`).** Before any
data is read, the engine walks the calendar once: for every day it asks the arc
resolver (STEP 5) which pieces exist; a day with at least one piece gets ONE
`'day'` block = all `'day'` rows in table order; a day without pieces gets no
block. Column offsets follow in day order, nl is the total. For your table:

```
B = [ day 59396: columns 1..72 | day 59397: columns 73..144 ]       nl = 144
      within a block: GraceA_x0, GraceA_y0, ..., GraceB_zaccspline#7 (table order)
```

`buildDay` for day 59397 receives `info.dayOffset = 72`, `info.dayCount = 72`;
the engine eliminates each day's block when the day is closed. Two consequences:
the ORDER of instances is the order of your rows (that is where your control
lives), and every day gets the SAME set — a template cannot say "this window
exists on day 59397 only" or "this scale is fixed on day 59396".

**Full control AND automation — the v5 proposal (R5-D8, to be decided).** Keep
the template form as the default and add an instance SELECTOR column:

```
column 1  scope   'arc' | 'day' | 'system'
column 2  days    []           = every day (today's behaviour)
                  [59396 ...]  = this row instantiates on these days only
columns 3…  yours
```

Two rows with the same name and different selectors are two instances (how
"the same parameter, regularised differently on day 2" is written); the block
of a day = the instantiated rows in table order; `info.dayRows` tells `buildDay`
which template row sits in column `info.dayOffset + k`; per-day counts may
differ; nl is still derived without reading; a selector naming a day outside
the calendar is refused at declaration. Your 144-row list converts mechanically:
group by name → 72 rows with `[]`; rows that differ on one day → two rows with
explicit selectors. (A resolver `localParams = @(day) rows`, like `arcs`, gives
the same power but hides the structure in a function; it stays engine-tier.)

### STEP 5 — the arc schedule = the pieces, resolvable WITHOUT reading data

`indexOfPiece` (1…4 for positions, 1…8 for KBR) are the Toeplitz blocks. `fit`
must know the schedule before it reads anything (it sizes the covariance grid
and the local blocks from it), so declare a resolver from `numberOfPieces`:

```matlab
nPieces = sd.numberOfPieces;                                   % days x groups
s.arcs  = @(day, obsType) 1:nPieces(day == data.days, groupIndex(obsType));
```

and let `buildDay` return `arcId = unique(indexOfPiece)` (STEP 9). The ragged
last KBR piece (2111 rows on 59396) needs nothing special: pieces are
addressed by id, not by position.

### STEP 6 — bind what `buildDay` needs

```matlab
po.days      = data.days;
po.obsTypes  = {'posA','posB','kbr'};
po.bundles   = sd.observations;                % the (days x groups) cell
po.designFcn = @gravfireDesign;                % your A/B builder (Part B)
s.processingOptions = po;
```

`fit` passes `processingOptions` untouched as `buildDay`'s fifth argument.
Nothing here is copied per call.

### STEP 7 — assemble the system

```matlab
s.id = 1;  s.name = 'GRACE';  s.buildDay = @buildGravfireDay;
s.covariance = cov;  s.localParams = localParams;  s.arcs = ...;  s.processingOptions = po;
data.systems = s;
```

`id` must be numeric or a single character (`'GRACE'` is rejected by
`lsqvce:fit:badId` — that was the one rule the first dry run hit). Why ONE
system: the KBR rows depend on the locals of **both** satellites; a system in
`fit` has one local block B shared by all its observation types, which is
exactly that coupling. Three systems would sever it.

---

## Part B — `buildDay` (what you write; the template has everything except STEP 11)

`fit` calls `D = buildDay(day, obsType, common, local, po, info)` **once per
(day, observation type)** and reads each day twice per pass (accumulate,
residuals). It must return raw, un-whitened matrices — the toolbox whitens,
weights, eliminates and cuts. Required fields: `A`, `B`, `l`, `t`; optional:
`arcIdx`, `arcId`, `nAxis`, `Cep`.

### STEP 8 — observations, epoch-major

`lobs` is `nEpoch × nAxis`; `fit` wants one column, epoch-major
(x₁ y₁ z₁ x₂ y₂ z₂ …), and the epoch of every row in seconds:

```matlab
D.l = reshape(lobs.', [], 1);           % m = nEpoch * nAxis
D.t = repelem(time(:,2), nAxis);        % seconds of day, one per row
```

Drop the epochs flagged in `nanIndexOfPiece` from `lobs`, `time`, `indexOfPiece`
and `iSpp` **together** before this (the template does; the file has none).

### STEP 9 — pieces → `arcIdx` + `arcId`

```matlab
pieceRow = repelem(indexOfPiece, nAxis);
ids      = unique(indexOfPiece).';                                   % [1 2 3 4]
D.arcIdx = arrayfun(@(k) find(pieceRow == k).', ids, 'UniformOutput', false);
D.arcId  = ids;
```

`arcIdx{k}` are the ROW numbers of piece `k`; `arcId` says which piece each
entry is. Every id the resolver (STEP 5) announced must appear exactly once.

### STEP 10 — the error information: `Cep` on the day (v5.3, R6 WP-90)

What the preprocessing knows about the accuracy of each observation comes
back **with the day**, as `D.Cep`; the covariance spec (STEP 3) declares only
the temporal law and its level. For the positions, `iSpp` holds the 3×3 epoch
covariance of the kinematic orbit as 9 columns per epoch. Hand it over as
`nEpoch × 3 × 3` — the axis count is derived from it:

```matlab
D.Cep = permute(reshape(iSpp.', 3, 3, []), [3 1 2]);     % nEpoch x 3 x 3
```

The toolbox composes C = D (R_t ⊗ I_3) D per piece: every epoch is de-scaled
by its 3×3 factor, then the temporal whitener runs down the epochs of each
axis; under `Estimate="noise"` the estimator de-scales the residuals the same
way before it forms the lag statistics. The KBR stream returns no `Cep` when
its epochs are equally accurate — or `D.Cep = sigmaPerEpoch` (nEpoch × 1)
when they are not: a scalar stream's per-epoch sigma travels the same way
(this is also where a per-epoch weight from an earlier screening would go,
instead of the `weights` of Part C). Nothing is declared twice: a `Cep` or
`d` on the spec that disagrees with the day is refused by name. **In
GRAVFIRE:** store `iSpp` as `nEpoch × 3 × 3` once, at preparation, instead of
reshaping it per pass.

### STEP 11 — the design matrices (yours)

```matlab
[D.A, D.B] = gravfireDesign(o, day, obsType, common, local, po, info);
```

with `A` of size `m × info.nc` (117) and `B` of size `m × info.nl` (144 for
two days), rows in the same epoch-major order as `l`. Where each local's
column is: **ask `info`, never compute it** —

```matlab
info.dayOffset          % column offset of THIS day's 72-wide 'day' block in B
info.dayCount           % 72
info.nc, info.nl        % 117, 144
```

so `B(:, info.dayOffset + (1:72))` is where GraceA's and GraceB's day-locals
go (in template row order: GraceA_x0 … GraceB_zaccspline#7); the other day's
72 columns stay zero (use `sparse` for B if you like). KBR rows fill both
satellites' columns; posA rows fill only GraceA's. `A` is the same 117
columns for every stream. The dry run replaces this step with seeded random
matrices so the rest of the pipeline can be checked before it exists.

### What `buildDay` must NOT do

Whiten, apply `Sll`/`weights`, regularise, or split the day into pieces for
the covariance. Return raw A/B/l/t and the piece index; the toolbox owns the
rest. A builder that scales its own rows is scaled twice.

---

## Part C — the things that are not (yet) in `data`

| GRAVFIRE item | where it goes now (v4) | where it goes in v5 (plan, HANDOVER §12) |
|---|---|---|
| **column 13 = column SCALE** (1e-9 on the Stokes rows, 1 / 1e-7 on the locals): the factor that balances the columns of A and B so that N has ≈ unit diagonal — *not* a prior (corrected 2026-09-02) | v4 has **no** column scaling: apply it in your design function (`A(:,j)*s_j`, `B(:,j)*s_j`) and undo it on `R.x` (`x = S*x̃`) and `R.Cxx` (`S*C̃*S`) yourself | a declared per-row `scale`: the toolbox solves the equilibrated system and unscales the results and the reported condition number (R5-D8) |
| **column 11 = regularisation** (`'none'` or `{'Tikhonov', 1, 0.01}`) — API parked by the user (R5-D8) | `data.priors(k) = struct('L', ..., 'sigma', ...)` on the common block (v5: a declared pseudo-observation group; `{'Tikhonov', 1, 0.01}` on row j is `L = e_j'`, `sigma = 0.01`, stacked over rows); `lsqvce.Options(Reg=..., RegMatrix=...)` remains the block-path regularisation | CLOSED (M7 7b): one `data.priors(k)` struct per label with `estimate = "scale"` is the `priorGroup` — its rows share ONE variance component, updated every pass by the Helmert / Koch–Kusche rule (`R.noise.priors`, `R.noise.priorHistory`) |
| `weights` (Huber, from the previous outer iteration) | `Robust="huber"` runs IRLS on the streaming kind (R5 M2) — iteration state stays the toolbox's; a FIXED per-epoch accuracy from an earlier screening goes into `D.Cep` (nEpoch × 1, STEP 10, v5.3) | unchanged |
| `covarianceModel 'fftmodel'`, `correlationLength 2160` | `Options(Estimate="noise", Noise=struct('Grid', 2160))` | per-stream `estimate="law"`, `Noise.Grid` from the schedule |
| `maximumIterations 7`, `minimumIterations 3` | `MaxPasses=7`; **no** minimum knob | `MaxPasses=7`, `MinPasses=3` (added on this evidence) |
| `minimumNumberObservations 100` | `MinRedundancy=100` | `Noise.FreezeRedundancy` |
| `limitConditionNumber 1e13` | not a knob: rank/condition refusals + reported condition | unchanged |
| `lest`, `iteridx`, `Sacl`, `facl`, rotations, `niPhi`/`niZhi`, `nX12`/`ndX12`/`drho` | inputs of your `gravfireDesign` only | unchanged |

---

## Part D — running it

```matlab
S    = load('testdata.mat');
data = gravfire_data(S.commonData, S.myOptions, S.systemsData, 'DryRun', true);
R    = lsqvce.fit(data, lsqvce.Options(Verbose=1));        % 2 s: pipeline check

data = gravfire_data(S.commonData, S.myOptions, S.systemsData, 'DesignFcn', @gravfireDesign);
R    = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=7, Noise=struct('FreezeRedundancy', 100), OnPass=@(ctx) fprintf('pass %d: max|alpha2-1| = %.3g\n', ctx.pass, ctx.maxDev)));
R.report;  T = R.coefficients;  C = R.convergence;
```

> **Multi-axis streams on the joint estimator (R5 M3, 2026-09-02):** the 3-axis
> position pieces are processed per axis after the per-epoch de-scale by `Cep`; both
> days at 2160 lags take ~40 s (HANDOVER §12.5, round 3c). Since v5.3 the positions
> declare `diag` and return `Cep` (STEP 10) — no `separable` type, no axis count in the
> spec.

If `fit` rejects the declaration, the identifier names the rule
(`lsqvce:fit:unknownField`, `:badId`, `:badScope`, …) and no data has been
read. If it fails while reading, the identifier is `lsqvce:makeLoader:*` and
the day and observation type are in the message; your own error inside the
design function arrives as the **cause** of `lsqvce:makeLoader:buildDay`.

## Part E — the checklist for GRAVFIRE's next release of the bundle

1. Emit `data` (STEPS 1–7) directly; keep the per-day bundles under
   `processingOptions` or in the closure.
2. Local table as 72 templates with seconds-of-day (STEP 4), not 144 expanded rows.
3. Observation types by name (STEP 3).
4. `iSpp` as `nEpoch × 3 × 3` (STEP 10).
5. No `weights`, `lest`, `iteridx` in the bundle (iteration state is the toolbox's).
6. `gravfireDesign` deterministic: same `(day, obsType)` → byte-identical A/B.
7. Priors stay numbers in the tables (Part C), never applied in the builder.
