# lsqvce release notes

## 6.1 — 2026-09-09 (R6 WP-105: the robust pass inside the joint loop)

Additive on 6.0; no declaration changes.

- **`Robust` together with `Estimate="noise"` runs.** Every pass of the joint estimator ends
  with one residual sweep; its law-whitened residuals (the current law, level, arc factor
  and `Cep` divided out) go through the influence function (`huber` | `igg3` | `bisquare`,
  `RobustOpts.Tuning`) and the factors ride inside every arc's whitener on the next pass, so
  accumulation, probe redundancies and the component update see ONE re-weighted problem.
  `MaxPasses` is the budget; `RobustOpts.MaxIter` / `Tol` belong to the stand-alone IRLS
  only. Prior streams are never re-weighted.
- What comes back: `R.noise.robust` (`influence`, `scale` and `streams` — one MAD scale per
  stream, `factor` per observation in sweep order, `keys`), `R.history.robustScale` (one
  column per stream) / `nDown` / `nZero` per pass, and the same fields on the `OnPass`
  context. The residuals are standardised by the stream's own MAD, never by the estimated
  level (a global scale, or the level itself, sends a clean stream's level into a spiral —
  measured). Prefer `igg3` / `bisquare` for gross blunders: Huber keeps c/|u| of their
  weight.
- The price, stated: the levels are those of the RETAINED observations — no consistency
  correction is applied. Under Huber (c = 1.345) on Gaussian residuals the effective level
  ratio is E[w u²]/E[w] ≈ 0.80 (σ ≈ 0.89); the oracle in `test_vcePolicy` measures it.
  This is Koch's iterated re-weighting inside a VCE, not Ellmer (2018), whose robustness is
  preprocessing screening plus the arc-wise variance factors (Manual chapter 7, §7.8).
- New identifier `lsqvce:adjustVCEjoint:collapsed` (beyond the breakdown point);
  `lsqvce:fit:notYet` now fires only for `Robust` on a nonlinear (`x0`) batch problem.

## 6.0 — 2026-09-09 (R6 WP-102: the cells are the arcs)

BREAKING for every `buildDay` (the sixth contract change: 2.0 namespace, 3.0
counts derived, 4.0 scoped tables, 5.0 one loop, 6.0 this); everything else
about the declaration — tables, priors, `Cep`, `Options` — is unchanged.

- **`buildDay` returns ONE CELL PER ARC:** `D.A`, `D.B`, `D.l`, `D.t` are
  1 × nArc cells, cell a = arc a (`{A}` for a single arc, empty cells for a
  day without data); `D.Cep` is a matching cell (nEpoch_a × 1 / × d / × d × d
  per arc); `D.nAxis` stays scalar. Arc lengths differ freely. The v5 fields
  `arcIdx` / `arcId` are GONE — a returned `arcIdx`, `arcId` or any field
  beyond `A B l t Cep nAxis` is refused by name (`lsqvce:makeLoader:badDay`),
  as is a matrix where a cell belongs. Driver: the GRAVFIRE preprocessing
  delivers every day as a cell array, one cell per noise-independence piece
  — the day's cells now ARE the arcs, no bookkeeping in between.
- **`arcs` is OPTIONAL.** Omitted, the arcs of a day are its cells (ids
  1 … nArc), learned when the engine's per-day `meta` probe reads the day
  (which the engine reads right after; the cache holds a whole day). Then
  `'day'` template rows and common `'day'` rows live on EVERY calendar day
  of the system, and a calendar day returning no arc in any stream is
  refused at its probe (`lsqvce:makeLoader:emptyDay`) — take it off the
  calendar. `arcs` is REQUIRED only by `'arc'` template rows and the
  engine-tier `LocalsPerArc` (`localPlan`: nl must be known before a read).
- **Declared schedule:** cell a is the a-th scheduled id and the count must
  match (`lsqvce:makeLoader:arcCount`); an **EMPTY cell is the data gap** of
  a scheduled arc (replaces the v5 `arcId` subset). An id outside the
  schedule stays `badArc`. `badArcId` is retired.
- **BREAKING (no aliases):** the string form `arcs = 'fromDay'` is retired
  everywhere (`lsqvce:fit:badArcs` / `lsqvce:makeLoader:badArcs`): omitting
  `arcs` is what it meant.
- `CacheDays` now means what it says: one resident DAY is one slot per data
  stream (the key is day|stream; with one slot a 3-stream day evicted what
  the engine was about to read).
- Test instrument `dev/tests/private/arcCells.m`: the ONE cut from the v5
  row-index day shape to the v6 cells (schedule order, gaps, per-arc `Cep`),
  so every fixture keeps generating whole-day matrices. 558 points (+5).
- Docs: the Manual chapter 3 (§0c, §2.3, §3, §4.0–4.3, §6, §7, §8, §10),
  `makeLoader` / `fit` help, both templates (`gravfire_data` STEP 5 gone,
  STEPS 8–11 = the cut), `dev/GRAVFIRE_INTERFACE.md`, every example
  builder; record: `dev/HANDOVER.md` §12.26.

**Migration v5 → v6**

| v5 | v6 |
|---|---|
| `D.A = A; D.B = B; D.l = l; D.t = t; D.arcIdx = {rows1, rows2}` | `D.A = {A(rows1,:), A(rows2,:)}` … one cell per arc for `A`, `B`, `l`, `t` |
| a single-arc day: matrices, no `arcIdx` | `D.A = {A}`, `D.B = {B}`, `D.l = {l}`, `D.t = {t}` |
| `D.Cep` per day-epoch | `D.Cep = {Cep1, Cep2}` per arc |
| `D.arcId = [1 3]` with `s.arcs = @(day,~) [1 2 3]` | cells `{arc1, {}, arc3}` — the empty cell is arc 2's gap — or a resolver returning `[1 3]` |
| `s.arcs = 4` on a `'system'`/`'day'`-scoped system | omit `arcs`: the four cells are the arcs |
| `s.arcs = 'fromDay'` | omit `arcs` |
| `'arc'` template rows | unchanged: keep the schedule (count, vector or resolver) |

## 5.3 — 2026-09-06 (R6 WP-90: the error information travels with the day)

Additive on 5.2 for every declaration that used `diag` / `arma` / `psd` laws;
BREAKING for two seams, named below.

- **`buildDay` returns the error information with the day, for ANY stream:**
  `D.Cep` is nEpoch × 1 (the per-epoch standard deviation of a scalar stream),
  nEpoch × d (per-axis standard deviations) or nEpoch × d × d (the full
  per-epoch covariance, e.g. from a PPP solution). The axis count is derived
  from it (`nAxis` stays only for a multi-axis stream without `Cep`; given as
  well, it must agree). The covariance spec declares the TEMPORAL law only
  (`diag`, `arma`, `psd`, `acf`) and its level; the toolbox composes
  C = D (R_t ⊗ I_d) D per arc in ONE place (`covModel`), for every law,
  and the joint estimator de-scales the residuals of every arc that carries
  `Cep` — a scalar stream with a per-epoch sigma is estimated correctly
  (before 5.3 that information could not travel with the day at all, and a
  `Cep` on a `diag` stream was silently ignored).
- Validated at the day, with day and stream in the message:
  `lsqvce:makeLoader:badDayCep` (shape), `:axisConflict` (`nAxis` vs `Cep`,
  or a `d` on the spec), `:cepTwice` (a `Cep` on the spec and on the day),
  `lsqvce:covModel:cepOnVar`; a multi-axis day without `Cep` is composed with
  D = I and reported once per run (`lsqvce:makeLoader:noCep`).
- Fixed: the day's `Cep` was sliced with the arc ID instead of the arc's
  position in `arcId` — a ragged day whose ids are not 1:n whitened an arc
  with another arc's epochs (R6 finding E-01).
- **BREAKING (no aliases):** the covariance type `separable` is retired —
  declare the temporal law (`arma` / `diag` / `psd`) and return `Cep` with
  the day (`covModel:unknownType` names the replacement); the
  `epochCovariance` field of `data.systems(k)` is retired (an unknown field
  now — its only consumer, `kernel.residualAcv`, had no engine caller since
  5.0 and is deleted); `makeLoader`'s `EpochCovariance` option and the
  loader's `epochCov` field are gone; `jointModelOpt`'s `Separable` option
  is gone (the loader composes).
- Docs: the Manual chapter 3 (§0b, §3.4, §3.5, §5, §7 "the error information",
  §11), chapter 7, the GRAVFIRE template and guide (positions: `diag` +
  `Cep` from iSpp), `lsqvce_template` (one stub), `demo_grace_dataLoader`.
  Review R6 and the test audit of the chain: `dev/HANDOVER.md` §12.14.

## 5.2 — 2026-09-03 (R5-D8: the parameter-table interface)

Additive on 5.1 except one retirement (`Options.Reg / RegLambda / RegMatrix`).

- **Regularisation is declared, in your columns.** `data.priors(j)` (common
  parameters, L in the columns of A) and **`data.systems(k).priors(l)`**
  (system k's local parameters, L in the columns of B) are one mechanism:
  one struct per group, declared by **`group`** (a regular expression on the
  table's name column — the toolbox builds L, one row per matched row per
  instance) or by **`L`** (used as is), plus `x0`, one `sigma`, `from`,
  `name`, `estimate`; pseudo-observations counted in the redundancy; `estimate = "scale"`
  estimates the group's sigma every pass by Helmert / Koch–Kusche — for local
  groups as the level of their own stream `prior_<name>`, identical to the
  common groups. `from` activates a group at a pass; a label (`group` on
  common priors, the struct on local ones) merges members.
- **`lsqvce.layout(data)`** announces the column layout before `fit`, so L
  can be written in it (`lay.systems(k).col(row, day, arc)`).
- **`data.paramColumns = struct('start', c1, 'end', c2)`** names the validity
  columns of both tables: a `'day'` / `'arc'` row is instantiated only on the
  days it overlaps; `buildDay`'s info carries the rows present.
- **`R.toTable(data)`** writes estimates and standard deviations back into the
  expanded n × 15 tables (columns 9/10 by default) with the instance map.
- The 15-column convention: the scope token is column 1; your 14 columns
  follow. Column 11 `{'Tikhonov', fromPass, sigma[, label]}` is translated by
  the GRAVFIRE template into the declarations above.
- Retired: `Options.Reg / RegLambda / RegMatrix` (19 Options fields);
  `localPriors(row)` and the `priorDay` handle of the intermediate rounds.
- Fixed on the way: an observation type with an underscore was cut by the arc
  key parsers; the prior arc now takes its law through the override path.

## 5.1 — 2026-09-03 (Review R5, levels M6–M7: additive on 5.0)

Additive; no v5.0 declaration or option changes meaning.

- **`data.priors(k).estimate = "scale"`** — a prior group gets its own variance
  component (Koch–Kusche): under `Estimate="noise"` its sigma is updated every
  pass by the Helmert rule σ² ← σ²·Ω_k/r_k (fixed point Ω_k = r_k, "the prior's
  weighted residual sum equals its redundancy"); `R.noise.priors`,
  `R.noise.priorHistory`. This closes the GRAVFIRE guide's `priorGroup` question.
- **`R.compare(Rred)`** — the nested-model F test between two fits
  (`lsqvce.engine.ftest` re-homed onto `Results`).
- **`R.plot("residuals" | "contribution" | "psd" | "convergence")`** in its
  final form; `Options.OnPass` unchanged since 5.0.
- **`examples/demo_gravfire.m`** — a GRAVFIRE bundle through
  `template/gravfire_data.m` and `fit` (pipeline check or `Estimate="noise"`
  with `OnPass` and `CacheDays=Inf`); every demo explains v5 in its prose.
- **Test tree**: every public unit has its seven categories tested or argued
  (`dev/tests/category_exemptions.csv` is the review surface); new suites
  `test_metamorphic` (invariances) and `test_speedBudget` (ratio budgets
  pinned to their own measurements). 523 points, 26 suites.
- **Docs**: `lsqvce_theory.pdf` at the root is the R5 theory (54 pp.); the
  HTML manuals carry a v5 banner naming the binding authorities;
  `DATA_CONTRACT.md` v5 documents priors, per-stream policy and `OnPass`.

Development record: `dev/HANDOVER.md` §12.11 (the R5 closing summary).

## 5.0 — 2026-09-03 (Review R5, levels M0–M5: the one-loop redesign)

**Breaking release.** The toolbox has ONE accumulator and ONE stochastic-model
estimator; the option space is 14 knobs + 5 typed groups; every retired name
has a successor in the migration table below (the authority is
`dev/HANDOVER.md` §12.10). Nothing on the data contract changed shape — v4
declarations still fit — but three things were added to it.

### What changed on the user surface

- **`Options` v5** — `Estimate="none"|"noise"` replaces `Model`; the v4 model
  words are explicit policy (`"vce"` = `Estimate="noise", Pool="stream",
  ArcFactors=false`; `"pooled"` = `Pool="all", ArcFactors=false`; `"ellmer"` =
  `Estimate="noise"`). One budget: `MaxPasses` / `Tolerance` (`MinPasses` was
  retired in R6 WP-101: it never had a consumer).
  `Robust` names the influence (`"huber"|"igg3"|"bisquare"`); groups
  `Noise`, `RobustOpts`, `Diagnostics`, `Parallel`, `GaussNewton` merge partial
  structs onto defaults. `CacheDays` (resident built days) and `OnPass` (see
  below) are new. `fit(data, Name=Value)` is the same mechanism as the object.
- **`Options.OnPass = @(ctx)`** — your inspection hook after EVERY pass of every
  iterative run (variance components, IRLS, Gauss–Newton): the pass state, the
  live accumulator, your declaration and every system's `processingOptions`;
  return `false` to stop (`stopReason = "userStop"`). Successor of the v4
  `plotPostModel` / `plotPostSolve` hooks and of `ShowSeconds` / `ExportDir`;
  post-mortem plotting is `R.plot("residuals"|"contribution"|"psd"|"convergence")`.
- **`data.priors(k) = struct('L', 'x0', 'sigma', 'name')`** — declared
  pseudo-observation groups on the common parameters. A free network's datum
  is the prior `L = B'` on the null space (the exact minimum-constraint
  solution, weight-independent); `DatumB` / `Sigma0` / `Sparse` are gone.
- **Per-stream policy on the covariance spec** — `estimate = "none" | "scale"
  | "law"` and `arcFactors` (logical); `Pool = "all" | "obsType" | "stream"`
  decides which streams share one estimated spectral shape; every stream keeps
  its own level. A law with a per-row `sigma` is a fixed structure (`"scale"`
  only).
- **Every fit kind rides the block path.** In-memory (`A, l, p`) and batched
  (`batchFcn`) problems become one local-free streaming system; `Robust`
  works on the streaming kind; Gauss–Newton (`x0`) is a policy on the same
  loop. A full SPD weight matrix per batch is whitened exactly. Per-observation
  diagnostics of a batch problem are reported in the WHITENED frame
  (`R.diag.v = sqrt(p).*v`); only their sums are frame-invariant.
- **`Results`**: `stopReason` (`"single" | "converged" | "maxPasses" |
  "userStop"`), `meta.warnings` (every warning identifier of the run),
  `noise.Cgroups` / `noise.stream` / `noise.model{k}.sigma`, `history.streams /
  sigma / white / whitePval` on every iterative run (so `convergence()` and
  `summary()` work for `Estimate="noise"`), `plot(kind)`. `R.source.verdict` →
  `R.report().verdict`; `history.nOuter` → `history.nIter`.

### What changed underneath (why the numbers may differ from 4.0)

- The joint (Ellmer) estimator is the one estimator (finding H-03): the
  separate mode's AR law fitted to the biased raw-residual ACF is retired
  (R5-D4). Acceptance is the physical oracle (planted fields and stream sigmas
  recovered within their formal errors), not equality to 4.0.
- Multi-axis streams are handled per axis after the per-epoch descale (H-15);
  the joint law lives in the PSD domain, positive by construction and resampled
  to every arc's length (H-17); the pooled shape is normalised and each stream
  carries its own level; per-lag sums are FFT autocorrelations (a pass on a
  two-day GRACE-FO bundle at 2160 lags: 12 min → 39 s).
- The QR backend refuses rank deficiency by `rank()`'s tolerance (not at the
  roundoff boundary — the refusal was platform-dependent) and folds a
  rank-deficient contribution in its range space (no `\` on a singular factor,
  no warning ahead of the refusal) (H-16).
- Error doctrine: every catch refuses, reports or counts; `util.callUser`
  wraps user functions with the original error as cause; the identifier
  catalogue is closed (`dev/docs_src/error_ids.csv`); skip RATES refuse
  (`Noise.SkipThreshold`).

### Migration table (v4 → v5)

| v4 | v5 |
|---|---|
| `Model="gls"` | `Estimate="none"` (default) |
| `Model="ellmer"` | `Estimate="noise"` (+ `Pool`, `ArcFactors=true` default) |
| `Model="vce"` | `Estimate="noise", Pool="stream", ArcFactors=false` |
| `Model="pooled"` | `Estimate="noise", Pool="all", ArcFactors=false` |
| `Robust="irls"` + `Influence` | `Robust="huber"|"igg3"|"bisquare"` |
| `Tuning`, `RobustMaxIter`, `RobustTol` | `RobustOpts.Tuning / .MaxIter / .Tol` |
| `MaxOuter`, `JointMaxIter`, `VceMaxIter`, `NonlinearMaxIter` | `MaxPasses` |
| `JointTol`, `OuterTol`, `VceTol`, `TolX` | `Tolerance` |
| `TolF`, `Lambda0` | `GaussNewton.TolF / .Lambda0` |
| `NumLags` | `Noise.Grid` |
| `PsdMode` | `Noise.Whitener` (`"toeplitz"` default, as before) |
| `MeanType`, `Smooth`, `MinRedundancy`, `SkipThreshold` | `Noise.MeanType / .Smooth / .FreezeRedundancy / .SkipThreshold` |
| `Parallel`, `DistributeDays`, `ParallelThreshold` | `Parallel.Mode / .DistributeDays / .Threshold` |
| `Alpha`, `Bonferroni`, `Beta`, `RMin`, `RWeak`, `Test`, `RobustCheck`, `Store`, `NWorst` | `Diagnostics.<same>` |
| `LogFid`, `LogFlag`, `ShowSeconds`, `ExportDir`, `plotPostSolve`, `plotPostModel` | `OnPass`; `R.plot(kind)` (`LogFile` was offered here and retired in R6 WP-101: it never opened a file) |
| `DatumB` | `data.priors = struct('L', B.', 'sigma', 1, 'name', 'datum')` |
| `Sigma0`, `Sparse`, `Verify`, `Progress` | retired |
| `GroupBy`, `ComponentsPerEpoch`, `MaxOrder`, `AcvTol`, `MaxLag`, `WhiteAlpha` | retired with the separate mode (groups are streams; the axis count comes from the builder; whiteness is fixed at K = 20, α = 0.05) |
| `Norm`, `C0`, `StepSchedule`, `AvgTail` | engine-tier only (`lsqvce.engine.adjustVCEjoint`) |
| `lsqvce.engine.adjust / robust / solveNonlinear / Accumulator / adjustVCE / Combine / BatchIndex / varComp / vceHutchinson` | gone; `fit` on the block path; `lsqvce.engine.adjustBlock / adjustRobust / adjustNonlinear / adjustVCEjoint / batchSystem` at engine tier |
| `R.source.verdict`, `history.nOuter` | `R.report().verdict`, `history.nIter` |

Full ledger, decisions and round notes: `dev/HANDOVER.md` §12; theory:
`lsqvce_theory.pdf` (R6 WP-101: it shipped under a second name in dev/ as well; that duplicate is gone -- one document, one name).

## 4.0 — 2026-08-31 (Level L11, WP-75: SCOPED parameter tables)

**Breaking: both parameter tables gain a binding scope column.** Column 1 of
`commonParams` (`'common' | 'day'`) and `localParams` (`'arc' | 'day' |
'system'`) is now the scope token; rows are **templates** the toolbox expands
over the calendar — nc and nl are derived from the expansion, and every column
offset reaches `buildDay` through its 6th argument. Driver: the GRACE
parameterization (arc-wise ICs, day-wise ACC scale/bias/drift, per-arc KBR
empiricals), which under v3 forced day-wise locals into hand-numbering.

### Migration table (v3 declaration → v4)

| v3 | v4 |
|---|---|
| `commonParams = T` (R×C, free) | prefix every row: `[repmat({'common'}, R, 1), T]` |
| `localParams = T` (R×C, free) | prefix every row: `[repmat({'system'}, R, 1), T]` |
| `localsPerArc = n` | delete the field; declare n `'arc'` rows in `localParams` |
| `localsPerArc = @(day, arcId) n` | not expressible as a template (constant rows per arc); heterogeneous per-arc counts stay ENGINE-tier (`lsqvce.engine.makeLoader(..., 'LocalsPerArc', fh)`) |
| day-wise parameters hand-packed into `commonParams` or hand-numbered in `localParams` | one `'day'` template row each; the toolbox numbers them |

### Layout (binding)
A: `['common' rows, table order | day blocks, calendar order]` — a calendar
day gets a block iff at least one system's **schedule** carries an arc that
day. B per system: per day `['arc' blocks, sorted ids | 'day' block]`, then
one trailing `'system'` block; a system day without arcs gets no blocks.
`info` (buildDay's 6th argument): `.nl .arcId .offsetOf(arcId) .count
.dayOffset .dayCount .sysOffset .sysCount .nc .commonCount .commonDayOffset
.commonDayCount`.

### Error-id changes
New `lsqvce:fit:badScope` (bad scope token; the localParams case is
makeLoader's `badScope`, relabelled). `lsqvce:fit:localCount` retired (one
local mechanism). `localsPerArc` on the USER surface is `unknownField`.

### Unchanged
`buildDay`'s signature and return contract, covariance specs, arc schedules,
`epochCovariance`, `lsqvce.Options`, `lsqvce.Results`, all engine numerics;
the ENGINE tier keeps v3 free-form tables and `LocalsPerArc` (default
`Scoped=false` on `lsqvce.engine.makeLoader`).

---

## 3.0 — 2026-08-31 (Level L10, WP-72: the v3 data contract)

**Breaking: the `data` contract states nothing the data already states.** Every
count is derived from the one place it already lives — nc from the
`commonParams` rows, nl from the `localParams` rows or the `localsPerArc` plan
(neither → nl = 0), the observation types from `fieldnames(covariance)`, u from
`numel(x0)` on the nonlinear batched kind. The retired v2 fields are **unknown
fields** now, so a v2 declaration fails loudly with `lsqvce:fit:unknownField`
naming the field. `lsqvce.engine.*` keeps the v2 conveniences (counts,
`'fromDay'`, `cache`, plot hooks, raw loader closures) — they left the USER
surface, not the toolbox.

### Migration table (v2 declaration → v3)

| v2 | v3 |
|---|---|
| `data.common.commonParams = T` | `data.commonParams = T` |
| `data.common.numCommon = nc` | `data.commonParams = num2cell((1:nc).')` — the table is the declaration |
| `s.observationTypes = {'kbr','pos'}` | delete — the types are `fieldnames(s.covariance)`, in declaration order |
| `s.systemId = 'C'` | `s.id = 'C'` (optional; defaults to the system's index) |
| `s.numLocal = nl` | `s.localParams = num2cell((1:nl).')` — or `s.localsPerArc` |
| `s.numLocal = 0` | delete — nl = 0 is the default (neither `localsPerArc` nor `localParams`) |
| `s.localParams` (R×14) | unchanged name; the fixed 14 columns are gone — R×C, C your convention (same for `commonParams`) |
| `s.arcs = 'fromDay'` | a resolver `@(day,obsType)` (`'fromDay'` is ENGINE-tier only: `lsqvce:fit:badArcs`) |
| `s.cache`, `s.plotPostSolve`, `s.plotPostModel` | delete — ENGINE-tier (`lsqvce.engine.makeLoader`), not data declaration |
| `s.loader = @(d,ot,a,uo) ...` | ENGINE tier: call `lsqvce.engine.adjustBlock` / `adjustVCE` with the struct-array form directly |
| `data.numUnknowns` with `data.x0` (kind C) | delete — u = `numel(x0)` (`lsqvce:fit:uCount` if both); the linear form keeps `numUnknowns` |

### Error-id changes
`lsqvce:fit:commonCount` → `lsqvce:fit:noCommon` (commonParams missing or not a
2-D cell). `lsqvce:fit:localCount` now means *both* local ways declared. New:
`lsqvce:fit:badId`, `lsqvce:fit:badArcs`, `lsqvce:fit:uCount`.
`lsqvce:makeLoader:paramCols` retired (no fixed column count to violate).

### Unchanged
`buildDay`'s six-argument signature and return contract, the covariance specs,
arc schedules (static and resolver), `arcInfo`, the multi-axis rules incl. the
`epochCovariance` hook (§7 documents *why* it exists: the exact
separable-descale path of the VCE Shape step), `lsqvce.Options`,
`lsqvce.Results`, all engine numerics.

---

## 2.0 — 2026-09-01 (R4 track, Levels L1–L6)

**The one breaking release of the R4 redesign.** Version 1.x users: read the
migration table at the end; everyone else starts at `doc lsqvce` → Getting
started.

### The user surface is four names
`data` struct (DATA_CONTRACT.md) + `lsqvce.Options` + `lsqvce.fit` + `lsqvce.Results` (`.report`).
`fit` covers three kinds of data (streaming with days/arcs/local unknowns,
in-memory matrices, batched callbacks incl. nonlinear) and four models
(`gls`, `vce`, `pooled`, `ellmer`), with `Robust="irls"` on the in-memory kinds.

### BREAKING
1. **Namespace layering (R4-D1, WP-64a).** Engine and kernel functions moved
   into the sub-packages `lsqvce.engine` and `lsqvce.kernel`; the root package
   holds the USER tier only. **No aliases** (toolbox rule 10): old names stop
   working. The toolbox's own error and warning identifiers
   (`lsqvce:<unit>:<rule>`) are unchanged; only the MATLAB-generated ids of
   `validateattributes` (`MATLAB:lsqvce:engine:<fn>:<check>`) now carry the
   qualified name, because they are built from it.
2. **`Results` from the in-memory kinds**: `R.sigma0` is the a-posteriori std
   of unit weight and `R.sol.Cxx` the covariance scaled by it (the streaming
   contract); the driver's a-priori `Sigma0` is kept as `sol.sigma0Prior`.
3. `adjust`'s help describes its `sol` as it is (no `.Cxx`; `sigma0` is the
   a-priori value) — the documentation changed, the numbers did not.

### Fixed
- **G-06** the raw-closure template ignored `adjustVCE`'s per-satellite reweighting.
- **G-07** a system with no local unknowns (nl = 0) is valid.
- **G-08** `jointModelOpt`'s covariance lookup assumed unit sampling; for
  dt ≠ 1 s every arc was skipped as non-PD.
- The batch→`Results` bridge establishes the streaming `Results` contract.

### Added
`lsqvce.fit`, `lsqvce.Options` (59 typed knobs, engine defaults marked in
`disp`), `Results.coefficients()`, `DATA_CONTRACT.md`, the Getting-started page
and examples gallery, `lsqvce_theory.pdf` (the ONE theory authority, R4-D3 —
the manual's theory chapter now points there), `DistributeDays` /
`ParallelThreshold` through `fit`, audit check 18 (tier census + Windows case
guard). `buildfile.m` (`buildtool`) existed L5-L8 and was removed in L9: an
mltbx is not foreseen; the gate is `cd('dev/tests'); ci_gate`, the audit
`python3 dev/docs_src/audit_tree.py`.

### Fluent configurators removed (R4-D2, WP-64b, Level L6)
`lsqvce.engine.adjustVCE` is **constructor-only**. The six chained
configurators are gone (no aliases); their option groups are Name-Value
**structs** on the constructor, merged into the engine defaults with the same
case-insensitive key matching the methods used:

| 1.x fluent call | 2.0 constructor Name-Value |
|---|---|
| `.groupComponentsBy(key)` | `'GroupBy', key` |
| `.stochasticModel('MaxOrder',p, ...)` | `'Shape', struct('MaxOrder',p, ...)` |
| `.varianceComponents('MaxIter',n, ...)` | `'Scale', struct('MaxIter',n, ...)` |
| `.iteration('MaxOuter',k,'Tol',t)` | `'Outer', struct('MaxOuter',k,'Tol',t)` |
| `.joint('Probes',m, ...)` | `'VceMode', 'joint', 'Joint', struct('Probes',m, ...)` — the method used to switch the mode; the group alone does not |
| `.componentsPerEpoch(map)` | `'Axes', struct(obsType, count)` |

`P.plots(...)`, `run()`, `empiricalStochastic()`, `validate()` are unchanged.
Through the USER front door nothing changes: `lsqvce.fit` passes exactly these
groups from `lsqvce.Options`. The rewrite of the toolbox's own tests and
examples was scripted (`dev/docs_src/refactor_fluent.py`, kept as protocol).

### Migration table (1.x → 2.0)
| 1.x | 2.0 |
|---|---|
| the direct driver calls (`lsqvce.adjust(u, batchFcn, n)`, …) | `R = lsqvce.fit(data, opts)` (recommended), or the qualified engine name below |
| `lsqvce.Accumulator` | `lsqvce.engine.Accumulator` |
| `lsqvce.BatchIndex` | `lsqvce.engine.BatchIndex` |
| `lsqvce.BlockAccumulator` | `lsqvce.engine.BlockAccumulator` |
| `lsqvce.Combine` | `lsqvce.engine.Combine` |
| `lsqvce.CommonPrior` | `lsqvce.engine.CommonPrior` |
| `lsqvce.NormalEqSystem` | `lsqvce.engine.NormalEqSystem` |
| `lsqvce.PrefetchLoader` | `lsqvce.engine.PrefetchLoader` |
| `lsqvce.QRSystem` | `lsqvce.engine.QRSystem` |
| `lsqvce.ReducedContribution` | `lsqvce.engine.ReducedContribution` |
| `lsqvce.adjust` | `lsqvce.engine.adjust` |
| `lsqvce.adjustBlock` | `lsqvce.engine.adjustBlock` |
| `lsqvce.adjustVCE` | `lsqvce.engine.adjustVCE` |
| `lsqvce.adjustVCEjoint` | `lsqvce.engine.adjustVCEjoint` |
| `lsqvce.applyRegularisation` | `lsqvce.engine.applyRegularisation` |
| `lsqvce.ar2arma` | `lsqvce.kernel.ar2arma` |
| `lsqvce.arcVarianceFactors` | `lsqvce.engine.arcVarianceFactors` |
| `lsqvce.armaACV` | `lsqvce.kernel.armaACV` |
| `lsqvce.armaWhiten` | `lsqvce.kernel.armaWhiten` |
| `lsqvce.checkPSD` | `lsqvce.engine.checkPSD` |
| `lsqvce.collectResiduals` | `lsqvce.engine.collectResiduals` |
| `lsqvce.contribution` | `lsqvce.engine.contribution` |
| `lsqvce.covComponent` | `lsqvce.engine.covComponent` |
| `lsqvce.covModel` | `lsqvce.engine.covModel` |
| `lsqvce.crossCorr` | `lsqvce.kernel.crossCorr` |
| `lsqvce.datum` | `lsqvce.engine.datum` |
| `lsqvce.dctI` | `lsqvce.kernel.dctI` |
| `lsqvce.decorrFactor` | `lsqvce.kernel.decorrFactor` |
| `lsqvce.diagnose` | `lsqvce.engine.diagnose` |
| `lsqvce.diagnoseAccum` | `lsqvce.engine.diagnoseAccum` |
| `lsqvce.ellipse` | `lsqvce.engine.ellipse` |
| `lsqvce.ftest` | `lsqvce.engine.ftest` |
| `lsqvce.gcvLambda` | `lsqvce.engine.gcvLambda` |
| `lsqvce.iterateNoiseModel` | `lsqvce.engine.iterateNoiseModel` |
| `lsqvce.jointModelOpt` | `lsqvce.engine.jointModelOpt` |
| `lsqvce.leverage` | `lsqvce.engine.leverage` |
| `lsqvce.leverageBlocks` | `lsqvce.engine.leverageBlocks` |
| `lsqvce.levinsonDurbin` | `lsqvce.kernel.levinsonDurbin` |
| `lsqvce.localDiagnostics` | `lsqvce.engine.localDiagnostics` |
| `lsqvce.makeLoader` | `lsqvce.engine.makeLoader` |
| `lsqvce.mapreduce` | `lsqvce.engine.mapreduce` |
| `lsqvce.normalizeSatList` | `lsqvce.engine.normalizeSatList` |
| `lsqvce.outliers` | `lsqvce.engine.outliers` |
| `lsqvce.poolGroups` | `lsqvce.engine.poolGroups` |
| `lsqvce.psd2acv` | `lsqvce.kernel.psd2acv` |
| `lsqvce.psdVCE` | `lsqvce.engine.psdVCE` |
| `lsqvce.redundancy` | `lsqvce.engine.redundancy` |
| `lsqvce.redundancyLags` | `lsqvce.kernel.redundancyLags` |
| `lsqvce.residualAcv` | `lsqvce.kernel.residualAcv` |
| `lsqvce.residuals` | `lsqvce.engine.residuals` |
| `lsqvce.resolveNc` | `lsqvce.engine.resolveNc` |
| `lsqvce.resolveSatSpecs` | `lsqvce.engine.resolveSatSpecs` |
| `lsqvce.robust` | `lsqvce.engine.robust` |
| `lsqvce.satKeyOf` | `lsqvce.engine.satKeyOf` |
| `lsqvce.selectArAcv` | `lsqvce.engine.selectArAcv` |
| `lsqvce.selectArOrder` | `lsqvce.engine.selectArOrder` |
| `lsqvce.separableWhiten` | `lsqvce.kernel.separableWhiten` |
| `lsqvce.sharedProbe` | `lsqvce.engine.sharedProbe` |
| `lsqvce.snoop` | `lsqvce.engine.snoop` |
| `lsqvce.solveNonlinear` | `lsqvce.engine.solveNonlinear` |
| `lsqvce.stat` | `lsqvce.engine.stat` |
| `lsqvce.streamContribution` | `lsqvce.engine.streamContribution` |
| `lsqvce.synthLoader` | `lsqvce.engine.synthLoader` |
| `lsqvce.systematic` | `lsqvce.engine.systematic` |
| `lsqvce.toeplitzCofactor` | `lsqvce.kernel.toeplitzCofactor` |
| `lsqvce.traceEstimate` | `lsqvce.engine.traceEstimate` |
| `lsqvce.varACV` | `lsqvce.kernel.varACV` |
| `lsqvce.varComp` | `lsqvce.engine.varComp` |
| `lsqvce.varWhiten` | `lsqvce.kernel.varWhiten` |
| `lsqvce.vce` | `lsqvce.engine.vce` |
| `lsqvce.vceHutchinson` | `lsqvce.engine.vceHutchinson` |
| `lsqvce.whitenessCorr` | `lsqvce.kernel.whitenessCorr` |
| `lsqvce.whittle` | `lsqvce.engine.whittle` |

Created 2026-09-01 (CEST, UTC+02:00) with assistance from Claude (Anthropic), model Claude Fable 5 — R4 Level L5, WP-65.
