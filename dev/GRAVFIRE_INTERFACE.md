# GRAVFIRE → `lsqvce.fit`: preparing the data, step by step

*WP-161, 2026-11-02: rewritten for `lsqvce.TaskDefinition` (level v7). The declaration is a
task filled by add-calls; each stream's `Build` is `out = f(day, ctx)` and adds one arc per
piece; the `DesignFcn` is `[A, B] = f(o, day, ctx)`; `gravfire_data` returns `[task, data]`.
Updated 2026-11-02 (CET, UTC+01:00) with assistance from Claude (Anthropic), model Claude
Opus 5.5.*

*R5 M0 close-out, 2026-09-02. Companion of `template/gravfire_data.m`, which implements
every step below. Read this guide with that file open; the STEP numbers match its comments.
Measured on `testdata.mat` (2 days, degree 10) on R2026a, 2026-09-02, with the hand-written
declaration struct of that time: it passed every `lsqvce:fit` check and a dry run with random
design matrices accumulated both days in 2 s. Since v7 the task compiles to that engine struct
(`task.declaration()`; WP-160 pins the compilation bitwise on its own fixture, not on
`testdata.mat`). Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
(Anthropic), model Claude Fable 5.1. Updated 2026-09-09 for v6 (R6 WP-102): THE PIECES ARE
THE ARCS — no arc schedule is declared (STEP 5 is empty; STEPS 8–11 renumbered). Model
Claude Fable 5.1. Updated 2026-10-17 (WP-148): the flagged epochs stay on the grid with
weight 0 (STEP 8). Model Claude Opus 5.5.*

---

## What you are building

One `lsqvce.TaskDefinition`, `task`, that `lsqvce.fit` reads. It is filled with flat calls;
each call checks its arguments at once and names the one at fault, before any file is read:

| call | what it declares | STEP |
|---|---|---|
| `task = lsqvce.TaskDefinition(Days=days)` | the calendar | 2 |
| `task.addCommon(names, Info=T)` | the shared unknowns: your 117 Stokes rows | 1 |
| `sys = task.addSystem("GRACE", Settings=po, Id=1)` | ONE system: GraceA + GraceB + KBR | 6 (`po`), 7 |
| `sys.addLocal(names, Scope="day", Info=T)` | one day's local templates, 72 rows | 4, 7 |
| `sys.addStream(type, Build=@b, Noise=n)` | the three observation types `posA`, `posB`, `kbr` | 3, 7 |
| `task.addPrior(...)`, `sys.addPrior(...)` | column 11 (Tikhonov) of the two tables | 6, 7 |
| *(no `Arcs`)* | the pieces `Build` adds are the arcs of the day | 5 |

`fit` **derives** every count from this (nc = 117, nl = 144, arcs per day, observation types)
and enforces them against every day it reads. You never type a count. `disp(task)` lists the
declaration; `task.layout()` names every column of A and B.

```matlab
[task, data] = gravfire_data(commonData, myOptions, systemsData, 'DesignFcn', @gravfireDesign);
```

`data = task.declaration()` is the compiled engine-tier struct (`commonParams`, `days`,
`priors`, `systems`), for engine-tier work and for inspecting the compiled tables;
`fit(task)` and `fit(data)` are the same fit. That struct is documented in
`dev/docs_src/api/07_engine_declaration.md`; nothing in this guide needs it.

---

## Part A — the declaration (what `gravfire_data.m` does today; what GRAVFIRE should emit tomorrow)

The code makes the task first (STEP 2), because `addCommon` is one of its methods; the
STEPs are numbered as before.

### STEP 1 — the common table: one unique name per row, your columns as `Info`

Your `commonData.parameter` is 117 × 14 (`'c'`/`'s'`, degree, order, `'even'`, …, a-priori
sigma 1e-9 in column 13). Every parameter needs a **unique name**; the template's one rule
(`commonName`) is `type_degree_order`, and the 14 columns ride along unchanged as `Info`:

```matlab
P = commonData.parameter;                                   % 117 x 14 cell
task.addCommon(commonNames(P), Info=gravfireInfo(P));       % 'c_2_0', 's_2_1', ...
```

`gravfireInfo` turns the cell table into a MATLAB table with one row per name and the
variables `col1` … `col14` — GRAVFIRE's own column numbers, so column 11 (regularisation) is
`col11` and column 13 (the column scale) is `col13`. `fit` never reads `Info`; your `Build` and
`DesignFcn` read row j as `ctx.Common.Info{j}` (Part B). All rows are `Scope="common"` (the
default): one column of A each, in `addCommon` order. Two rows with one type, degree and order
are refused by the call (`lsqvce:TaskDefinition:duplicateName`). **In GRAVFIRE:** emit the
names with the table; keep column 13 and column 11 in it — STEP 6 and Part C say what happens
to them.

### STEP 2 — the calendar

```matlab
days = double(myOptions.days(:).');             % [59396 59397]
task = lsqvce.TaskDefinition(Days=days);
```

Integers (MJD is fine), each once. Every day listed must be a day the bundle has.

### STEP 3 — observation types become STREAM names

Today the three groups are addressed by number (`observations{day, g}`, `variance(g)`). `fit`
addresses them by **name**: one `addStream` per type, processed in the order added. Choose
short valid MATLAB identifiers (anything else is `lsqvce:TaskDefinition:badName`); the
template's `obsTypeOf` does:

| group | `name` / `observationQuantity` | stream name |
|---|---|---|
| 1 | GraceA / kinematicPosition | `posA` |
| 2 | GraceB / kinematicPosition | `posB` |
| 3 | GraceX / range-rate | `kbr` |

(another quantity is `gravfire:unknownQuantity`), and declares the **starting** temporal law
of each as white, with the level as a standard deviation (not the variance):

```matlab
Noise = lsqvce.Noise.white(sqrt(sd.variance(g)))    % posA, posB: 0.1 m; kbr: 1e-6 m/s
```

The stream's `Estimate` is left at its default `"law"`: under `Options(Estimate="noise")` the
estimator estimates the temporal law (`fftmodel` in your options) **from this starting level**,
the shape from a white start; nothing about the FFT model is declared here, and nothing about
the per-epoch accuracy either — that travels with each arc as `Sigma` (STEP 9), which the
toolbox composes with the law (C = D (R_t ⊗ I_3) D). **In GRAVFIRE:** store the bundle as
`bundles.posA`, `.posB`, `.kbr` (or keep the `days × groups` cell and look the column up by
name, as the template's `Build` does). Never rely on group order alone.

### STEP 4 — the local table: ONE day's templates, all `"day"`-scoped

Your `systemsData{1}.parameter` is the **expanded** table: 144 rows = 2 days × 2 satellites ×
36 unknowns, with absolute times `[MJD sec]` and the day in the flags. `fit` wants the
opposite: the **36 rows per satellite for one day**, and it instantiates them on every calendar
day itself. So `templateFromExpanded`:

1. keeps only the rows of the first day (72 rows);
2. converts the absolute spans `[59396 21600]` → seconds of day `21600` (columns 2, 3; the
   spline windows are the same every day);
3. converts the knot table `[MJD.frac value]` → `[secondsOfDay value]` (column 4);
4. drops the day from the flag vector `[59396 1 1 1 0]` → `[1 1 1 0]` (column 6).

The scope is not a column any more: it is the `Scope="day"` of the `addLocal` call (STEP 7),
column 1 is the name. **In GRAVFIRE**, emit the 72-row template instead of the 144-row
expansion and the helper disappears. The result (72 × 14, its columns again as `Info`) looks
like:

```
{'GraceA_x0'}          {0}      {86390}  {[]}          {8640}  {[1 1 1 1]} ...
{'GraceA_xaccspline'}  {21600}  {86399}  {2161x2}      {6480}  {[1 1 1 0]} ...
```

The spans are **not** declared as validity (`addLocal`'s `Valid`): they travel as
`Info.col2` / `Info.col3` for your design function, and the templates of the first day serve
every calendar day. (A template limited to some days is `Valid=[start end]` in fractional MJD,
API reference, chapter *The task definition*; the template does not use it.)

Why no `"arc"` locals: your initial conditions span the whole day (`numarc = 1`), so nothing
is per piece. If one day you integrate per piece, those rows become `Scope="arc"` — and then
need an arc schedule (STEP 5).

**How the `"day"` scope works.** Without an arc schedule (STEP 5) EVERY calendar day gets
ONE `"day"` block = all 72 templates in `addLocal` order — the calendar is the declaration, and
a day that then returns no arc in any stream is refused when it is read
(`lsqvce:makeLoader:emptyDay`: take it off the calendar). The blocks follow in calendar order,
nl is the total. For your table:

```
B = [ day 59396: columns 1..72 | day 59397: columns 73..144 ]       nl = 144
      within a block: GraceA_x0, GraceA_y0, ..., GraceB_zaccspline#7 (addLocal order)
```

Your code never computes these positions: the `DesignFcn` takes the columns of today's locals
from `ctx.cols` by name (STEP 10), and `task.layout()` announces them
(`lay.systems(1).local(c).name`, `.day`). The engine eliminates each day's block when the day
is closed. Two consequences: the ORDER of the columns within a block is the order of your rows,
and every day gets the SAME set of 72.

### STEP 5 — the arc schedule: none (v6, R6 WP-102)

`indexOfPiece` (1…4 for positions, 1…8 for KBR) are the Toeplitz blocks, and THE PIECES ARE
THE ARCS: `Build` adds one arc per piece (STEP 11) and the toolbox takes them as that day's
arcs, ids 1…nPiece. Nothing is declared — `numberOfPieces` is not restated, no resolver, no
`myArcTable`. The ragged last KBR piece (2111 rows on 59396) needs nothing special: arc lengths
may differ freely.

A schedule (`addSystem(..., Arcs=...)`: a count, an id vector or `@(day, stream) ids`) is
REQUIRED only by `Scope="arc"` locals, because nl then counts one block per (day, arc) and must
be known before any day is read. Your locals are all `"day"`-scoped, so none is given. (With a
schedule, `ctx.ArcIds` lists the day's arcs, `addArc` fills them in order, and a scheduled id
that `Build` never adds is that id's data gap.)

### STEP 6 — bind what `Build` needs: the system's `Settings`

```matlab
po.days      = days;
po.obsTypes  = {'posA','posB','kbr'};
po.bundles   = sd.observations;                % the (days x groups) cell
po.designFcn = opt.DesignFcn;                  % 'DesignFcn': your A/B builder (Part B)
po.dryRun    = opt.DryRun;                     % 'DryRun', true: seeded random A/B instead
po.myOptions = myOptions;
```

STEP 7 makes it the system's `Settings` (`addSystem(..., Settings=po)`). `fit` hands
`Settings` untouched to every `Build` as `ctx.Settings` (the field name `prepared` is reserved).
The template declares no `Prepare`: each stream's `Build` reads its own (day, group) bundle.
Work that all streams of a day share — GRAVFIRE's orbit integration and its partials — would go
into the system's `Prepare` (`addSystem(..., Prepare=@f)`, run once per build of the day, its
result `ctx.Prepared`).

**Column 11 on the COMMON rows → `task.addPrior`** (STEP 6 in the template's comments,
R5-D8 round 1). `{'Tikhonov', k, sigma}` on row j means: constrain parameter j to its a-priori
value (column 10; `[]` = 0) with standard deviation `sigma`, from pass `k` on. An optional 4th
entry `{'Tikhonov', k, sigma, 'label'}` is a GROUP label. `priorsFromColumn11` makes ONE
prior group per label, declared by `L` in the columns of A (`addCommon` order):

```matlab
task.addPrior(name, L=L, Value=x0, Sigma=sigma, Estimate=est, From=k);
```

- rows sharing a label: one group named after the label, `L` stacks their unit rows e_j', `Value`
  their column-10 values (q × 1), one shared `sigma` estimated by VCE (`Estimate="scale"`,
  Koch-Kusche; several rows are what make a sigma estimable);
- a row without a label: its own group named after the row (`'c_2_0'`), `L = e_j'`, `sigma`
  fixed (`Estimate="none"`; a single row has redundancy at most 1);
- `'none'` or `[]`: no prior. Members of one label that disagree on `sigma` or on `fromPass`
  are `gravfire:badConstraint`.

In your bundle every common row is `[]` and the constraints sit on LOCAL rows (`*accscale`):
those become `sys.addPrior` in STEP 7 — the same mechanism and the same Koch-Kusche estimation.

### STEP 7 — assemble the system: its local templates and its streams

```matlab
sys = task.addSystem("GRACE", Settings=po, Id=1);
sys.addLocal(string(T(:, 1)), Scope="day", Info=gravfireInfo(T));       % 72 templates
for g = 1:G
    sys.addStream(obsTypes{g}, Build=@buildGravfireDay, ...
                  Noise=lsqvce.Noise.white(sqrt(sd.variance(g))));
end
```

`"GRACE"` is the system's name; its `Id` must be numeric or a single character (a word as
`Id` is refused when `fit` starts, `lsqvce:fit:badId` — that was the one rule the first dry run
hit). Why ONE system: the KBR rows depend on the locals of **both** satellites; a system has
one local block B shared by all its streams, which is exactly that coupling. Three systems
would sever it. All three streams share one `Build`; it tells them apart by `ctx.Stream`. No
`Axes` is declared: the default (1) leaves the axis count to the `Sigma` the positions return
with each arc (3, STEP 9).

**Column 11 on the LOCAL rows → `sys.addPrior`** (STEP 7 in the template's comments,
R5-D8 round 6). Regularisation of local parameters is DECLARED, never computed.
`gravfireLocalPriors` makes one group per label in the PATTERN form: a regexp on the local
names, e.g. `'^(GraceA_xaccscale|GraceB_xaccscale)$'`, from which the toolbox builds L with
one row per matched template per day — the group spans the calendar:

```matlab
sys.addPrior(name, On=pattern, Value=x0, Sigma=sigma, Estimate=est, From=k);
```

Labelled groups get `Estimate="scale"`, unlabelled templates one fixed-sigma group each. The
members of a labelled local group must share `sigma`, `fromPass` **and** the a-priori value
(column 10, WP-144): the pattern form builds L per member and per day in the toolbox's order,
so one value per member cannot be written without restating that order (`gravfire:badConstraint`
names the remedy: separate labels, or the group by `L=`). The alternative to `On=` is `L=`,
written in the columns of B as `task.layout()` announces them — never both
(`lsqvce:TaskDefinition:priorTarget`). A local group is reported as the stream `prior_<name>`
of the system (API reference, chapter *Priors*).

---

## Part B — `Build` (what you write; the template has everything except STEP 10)

`fit` calls `out = buildGravfireDay(day, ctx)` **once per (day, stream)** each time the day is
built; with the default `Options.CacheDays = 1` a day is built twice per pass (accumulate,
residuals). It reads what it needs from the `lsqvce.DayContext`:

```matlab
po = ctx.Settings;                             % STEP 6
obsType = char(ctx.Stream);                    % 'posA' | 'posB' | 'kbr'
o  = po.bundles{find(po.days == day, 1), groupIndex(po.obsTypes, obsType)};
```

and returns the `lsqvce.DayData` of `ctx.newDay()`, filled with one `out.addArc` per piece
(STEP 11); anything else is `lsqvce:DayContext:badReturn`. Every row is raw and un-whitened —
the toolbox whitens, weights and eliminates. Per arc: `A`, `l`, `t` (required), the local
design `B=`, the per-epoch accuracy `Sigma=`, the rows of weight 0 `Excluded=`. Same day, same
stream → byte-identical output.

### STEP 8 — observations, epoch-major

`lobs` is `nEpoch × nAxis`; `fit` wants one column, epoch-major (x₁ y₁ z₁ x₂ y₂ z₂ …), and
the epoch of every row in seconds:

```matlab
l = reshape(lobs.', [], 1);             % m = nEpoch * nAxis, the whole day
t = repelem(time(:,2), nAxis);          % seconds of day, one per row
```

(the cut into one arc per piece is STEP 11).

The epochs flagged in `nanIndexOfPiece` **stay on the grid with weight 0** (WP-148). Zero
their rows in `l`, `A` and `B` (STEP 10), give them any valid 3×3 in `Sigma` (the template
takes the mean of the others), and pass the rows as `Excluded` (STEP 11):

```matlab
flag = logical(nanIndexOfPiece(:));                  % nEpoch x 1
lobs(flag, :) = 0;                                   % NaN in the file: zero here
ex   = reshape(repmat(flag.', nAxis, 1), [], 1);     % m x 1, epoch-major
```

The toolbox deletes such an epoch from the law's covariance: it is out of the degrees of
freedom and its values never reach the solution; an excluded row that is not zero in A, B and
l is refused (`lsqvce:makeLoader:badDay`). Do **not** remove the flagged epochs: a removed
epoch leaves a gap inside its piece, and a `psd` law (`Estimate="noise"`), which counts lags in
samples, takes the gap for ONE step (`lsqvce:covModel:unevenGrid` warns). Your design function
therefore returns the rows of **every** epoch of the day, the flagged ones included.

### STEP 9 — the error information: `Sigma` (v5.3, R6 WP-90)

What the preprocessing knows about the accuracy of each observation comes back **with each
arc**, as `Sigma=`; the stream's `Noise` (STEP 3) declares only the temporal law and its level.
For the positions, `iSpp` holds the 3×3 epoch covariance of the kinematic orbit as 9 columns
per epoch. Hand it over as `nEpoch × 3 × 3` — the axis count is derived from it:

```matlab
Cep = permute(reshape(iSpp.', 3, 3, []), [3 1 2]);       % nEpoch x 3 x 3, whole day
```

(cut per piece in STEP 11: `Sigma=Cep(epochsOfPiece_k, :, :)`).

The toolbox composes C = D (R_t ⊗ I_3) D per arc: every epoch is de-scaled by its 3×3
factor, then the temporal whitener runs down the epochs of each axis; under
`Estimate="noise"` the estimator de-scales the residuals the same way before it forms the lag
statistics. `Sigma` is relative to the law's `sigma`: the starting covariance of a position arc
carries the declared level `sqrt(variance(g))` times the `iSpp` factor, and under
`Estimate="noise"` the level is estimated. The KBR stream passes no `Sigma` when its epochs are
equally accurate — or `Sigma=sigmaPerEpoch` (nEpoch × 1, standard deviations, not variances)
when they are not: a scalar stream's per-epoch accuracy travels the same way (this is also
where a per-epoch weight from an earlier screening would go, instead of the `weights` of
Part C). **In GRAVFIRE:** store `iSpp` as `nEpoch × 3 × 3` once, at preparation, instead of
reshaping it per pass.

### STEP 10 — the design matrices (yours)

```matlab
[A, B] = po.designFcn(o, day, ctx);
```

with `A` of size `m × ctx.NumCommon` (117, columns in `addCommon` order) and `B` of size
`m × ctx.NumLocal` (144 for two days), rows in the same epoch-major order as `l` — one row per
axis of **every** epoch, the flagged ones too (the template zeroes those afterwards). Another
size is `gravfire:designShape`. Where each local's column is: **ask `ctx.cols` by name, never
compute an offset** —

```matlab
B = sparse(m, ctx.NumLocal);
B(:, ctx.cols(["GraceA_x0" "GraceA_y0" "GraceA_z0"])) = PxyzA;   % this day's columns
```

The locals are `"day"`-scoped, so a name gives its column on THIS day; the other day's 72
columns stay zero (`sparse` is accepted). `ctx.cols("c_2_0")` gives a common's column of A;
one call names commons or locals, not both (`lsqvce:DayContext:mixedTables`), and a name that
is not a parameter is `lsqvce:DayContext:unknownName`. The GRAVFIRE columns of every template
are at hand as `ctx.Common.Info{j}` and `ctx.Local.Info{k}` (fields `col1` … `col14`: the
spans in `col2` / `col3`, the knots in `col4`, the column scale in `col13`). KBR rows fill both
satellites' columns; posA rows fill only GraceA's. `A` is the same 117 columns for every
stream. The template takes B as one matrix (`B=`); the by-name form `Local={names, P}` is the
alternative (API reference, chapter *Prepare and Build*). The dry run replaces this step with
random matrices, seeded per day and stream, so the rest of the pipeline can be checked before it
exists; without `DesignFcn` and without `DryRun` the read fails with `gravfire:designMissing`.

### STEP 11 — the cut: one arc per piece (v6)

```matlab
ids = unique(indexOfPiece).';                                  % [1 2 3 4]
out = ctx.newDay();
for k = 1:numel(ids)
    ep   = find(indexOfPiece == ids(k));                       % this piece's epochs
    rows = reshape(((ep(:) - 1) * nAxis + (1:nAxis)).', [], 1); % their epoch-major rows
    sigma = [];
    if nAxis > 1, sigma = Cep(ep, :, :); end                   % positions only
    out.addArc(A(rows, :), l(rows), t(rows), B=B(rows, :), Sigma=sigma, Excluded=ex(rows));
end
```

Arc k is piece k of the day, in piece order; the toolbox reads the count from the arcs added
(and the cache keeps a whole day resident, so nothing is read twice within a build). Every arc
holds whole epochs (`lsqvce:makeLoader:badArc` otherwise). Whole-day matrices with row lists
(`arcIdx`) were the v5 form and are gone.

### What `Build` must NOT do

Whiten, apply `Sll`/`weights`, regularise, or split the day into pieces for the covariance.
Return raw A/B/l/t per piece; the toolbox owns the rest. A builder that scales its own rows is
scaled twice.

---

## Part C — the things that are not (yet) in the declaration

- **column 13 = column SCALE** (1e-9 on the Stokes rows, 1 / 1e-7 on the locals): the factor
  that balances the columns of A and B so that N has ≈ unit diagonal — *not* a prior
  (corrected 2026-09-02). The task has **no** per-parameter scale: the factor reaches your
  design function as `Info.col13`; apply it there (`A(:,j)*s_j`, `B(:,j)*s_j`) and undo it on
  `R.x` (`x = S*x̃`) and `R.sol.Cxx` (`S*C̃*S`) yourself. (The declared per-row scale of the
  R5-D8 plan, HANDOVER §12, was not built.)
- **column 11 = regularisation** (`'none'` or `{'Tikhonov', 1, 0.01}`): `task.addPrior` on the
  commons (STEP 6), `sys.addPrior` on the locals (STEP 7). A labelled group with
  `Estimate="scale"` shares ONE variance component, updated every pass of `Estimate="noise"` by
  the Helmert / Koch–Kusche rule (commons: `R.noise.priors`, `R.noise.priorHistory`; locals:
  the stream `prior_<name>`). `lsqvce.Options` carries no regularisation: everything about the
  problem is declared on the task.
- **`weights`** (Huber, from the previous outer iteration): `Options(Robust="huber")` — the
  iteration state stays the toolbox's; a FIXED per-epoch accuracy from an earlier screening
  goes into `Sigma=` (nEpoch × 1 per arc, STEP 9).
- **`covarianceModel 'fftmodel'`, `correlationLength 2160`**:
  `Options(Estimate="noise", Noise=struct('Grid', 2160))`; what the estimator may change per
  stream is the stream's `Estimate` (`"law"`, the default; STEP 3).
- **`maximumIterations 7`, `minimumIterations 3`**: `MaxPasses=7`; **no** minimum knob.
  `MinPasses` was added on this evidence and RETIRED in R6 WP-101, where it was found to have
  had no consumer since the day it was added — help line, property, `disp` table, and nothing
  else. Say the word and it comes back as a floor on the pass loop.
- **`minimumNumberObservations 100`**: the nearest knob is
  `Options(Noise=struct('FreezeRedundancy', 100))` — a spectral bin whose redundancy is below it
  keeps its value instead of being updated.
- **`limitConditionNumber 1e13`**: not a knob — rank/condition refusals and the reported
  condition.
- **`lest`, `iteridx`, `Sacl`, `facl`, rotations, `niPhi`/`niZhi`, `nX12`/`ndX12`/`drho`**:
  inputs of your `gravfireDesign` only (through the bundle `o` or `ctx.Settings`).

---

## Part D — running it

```matlab
S    = load('testdata.mat');
task = gravfire_data(S.commonData, S.myOptions, S.systemsData, 'DryRun', true);
R    = lsqvce.fit(task, lsqvce.Options(Verbose=1));        % pipeline check

[task, data] = gravfire_data(S.commonData, S.myOptions, S.systemsData, 'DesignFcn', @gravfireDesign);
R    = lsqvce.fit(task, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, ...
           MaxPasses=7, Noise=struct('FreezeRedundancy', 100), ...
           OnPass=@(ctx) fprintf('pass %d: max|alpha2-1| = %.3g\n', ctx.pass, ctx.maxDev)));
R.report;  T = R.coefficients;  C = R.convergence;
```

> **Multi-axis streams on the joint estimator (R5 M3, 2026-09-02):** the 3-axis position
> pieces are processed per axis after the per-epoch de-scale by the epoch covariance; both days
> at 2160 lags took ~40 s (HANDOVER §12.5, round 3c). Since v5.3 the positions declare a white
> law and return their per-epoch covariance with the day (STEP 9; `Sigma=` since v7) — no
> separable type, no axis count in the declaration.

If a declaration call is wrong, it fails at once with `lsqvce:TaskDefinition:*` and names the
argument (`:duplicateName`, `:badName`, `:priorTarget`, …); the checks on the compiled
declaration run when `fit` starts (`lsqvce:fit:badId`, `:badPrior`, …). No data has been read
in either case. If it fails while reading, the identifier is `lsqvce:makeLoader:*` or
`lsqvce:DayContext:*` and the day and stream are in the message; your own error inside `Build`
or the design function arrives as the **cause** of `lsqvce:makeLoader:buildDay`.

## Part E — the checklist for GRAVFIRE's next release of the bundle

1. Emit what the add-calls take (STEPS 1–7) directly — unique names, the 14 columns as `Info`;
   keep the per-day bundles in the system's `Settings` (or in the closure of `Build`).
2. Local table as 72 templates with seconds-of-day (STEP 4), not 144 expanded rows.
3. Observation types by name (STEP 3).
4. `iSpp` as `nEpoch × 3 × 3` (STEP 9), cut per piece into `Sigma=` (STEP 11).
5. No `weights`, `lest`, `iteridx` in the bundle (iteration state is the toolbox's).
6. `gravfireDesign` deterministic: same `(day, stream)` → byte-identical A/B.
7. Priors stay numbers in the tables (column 11 → `addPrior`, Part C), never applied in the
   builder.
