# lsqvce release notes

## 5.2 — 2026-09-03 (R5-D8: the parameter-table interface)

Additive on 5.1 except one retirement (`Options.Reg / RegLambda / RegMatrix`).

- **Regularisation is declared, in your columns.** `data.priors(j)` (common
  parameters, L in the columns of A) and **`data.systems(k).priors(l)`**
  (system k's local parameters, L in the columns of B) are one mechanism:
  one struct per group (`L`, `x0`, `sigma`, `from`, `name`, `estimate`),
  pseudo-observations counted in the redundancy; `estimate = "scale"`
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
  `Estimate="noise"`). One budget: `MaxPasses` / `MinPasses` / `Tolerance`.
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
| `LogFid`, `LogFlag`, `ShowSeconds`, `ExportDir`, `plotPostSolve`, `plotPostModel` | `LogFile`; `OnPass`; `R.plot(kind)` |
| `DatumB` | `data.priors = struct('L', B.', 'sigma', 1, 'name', 'datum')` |
| `Sigma0`, `Sparse`, `Verify`, `Progress` | retired |
| `GroupBy`, `ComponentsPerEpoch`, `MaxOrder`, `AcvTol`, `MaxLag`, `WhiteAlpha` | retired with the separate mode (groups are streams; the axis count comes from the builder; whiteness is fixed at K = 20, α = 0.05) |
| `Norm`, `C0`, `StepSchedule`, `AvgTail` | engine-tier only (`lsqvce.engine.adjustVCEjoint`) |
| `lsqvce.engine.adjust / robust / solveNonlinear / Accumulator / adjustVCE / Combine / BatchIndex / varComp / vceHutchinson` | gone; `fit` on the block path; `lsqvce.engine.adjustBlock / adjustRobust / adjustNonlinear / adjustVCEjoint / batchSystem` at engine tier |
| `R.source.verdict`, `history.nOuter` | `R.report().verdict`, `history.nIter` |

Full ledger, decisions and round notes: `dev/HANDOVER.md` §12; theory:
`dev/lsqvce_theory_r5.pdf`.

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
