# lsqvce release notes

## 6.26.2 — 2026-10-17 (WP-149: a refused v'Pv no longer reaches the fit's Results)

A correction. Estimates are unchanged; only numbers that were NaN change.

- `Backend = "normal"` on a very tight fit: when the backend refused a cancelled v'Pv
  (`lsqvce:accum:cancellation`, `lsqvce:normalEq:cancellation`), the Results that `lsqvce.fit`
  returned carried `R.sigma0 = NaN` and NaN formal errors in `R.sol.Cxx`, although the warning
  said the residual pass supersedes the refused number. `R.sigma0`, `R.sol.vtv` and `R.sol.Cxx`
  now come from the verdict's own residual pass (no data is read again), for every estimate and
  every producer of a Results, and agree with the QR backend's to 1e-9 or better (measured
  2.7e-11 at most).
- A robust noise run (`Robust = "igg3"` / `"bisquare"`) on the same refusal now keeps the
  truncation correction of its hard solution (measured 1.3 % off QR's sigma0 without it).

## 6.26.1 — 2026-10-17 (WP-144: the package and template findings of the repository review)

Corrections. Estimates change only where they were wrong (the prior in the second point).

- `lsqvce.engine.synthLoader`, the loader template: its one-day cache was keyed on six named
  fields, so a run that differed only in `mArc`, `x_c`, `cov`, `dt` or another configuration
  field received the previous run's day, silently. It is now keyed on the whole configuration,
  compared by value, and its help states the rule for a copy.
- A common prior group with `estimate = "scale"` on `Backend = "normal"`: when the backend
  refused a cancelled v'Pv (a very tight fit), the group's update divided by NaN, the group kept
  its declared sigma, and the loop reported "converged". It now takes N⁻¹ from the cofactor, which
  is finite, and estimates the group as the QR backend does.
- `R.plot("psd")` draws every pool group's spectrum on a logarithmic axis (with several groups
  the axis was linear).
- `template/gravfire_data.m`: the members of a labelled local Tikhonov group must share
  `fromPass` and the a-priori value (a later member's were dropped); `gravfire:badConstraint`
  names the disagreement. The template now runs in the test suite on a generated bundle, and
  `demo_gravfire` runs its pipeline there instead of its no-bundle message.

## 6.26.0 — 2026-10-16 (WP-147: a stop rule that low-information spectral bins can meet)

The noise loop's stop rule. The estimator is unchanged; a run may now stop on an earlier pass,
and its estimate then differs from the later stop's by the step that remained (measured: within
0.1 of a standard error).

- A spectral bin has settled when it moved by at most `Tolerance·max(1, sqrt(200/r))`, where
  `r` is the bin's redundancy. A bin of redundancy 200 or more is held to `Tolerance` as
  before; a bin the data barely determine (the zero-frequency bin that per-arc parameters
  absorb, or many bins on short arcs) is held, at the default `Tolerance` 1e-2, to a tenth of
  its own standard error `sqrt(2/r)`. Before, every bin was held to `Tolerance`; on four test
  fixtures that took 14–168 passes, while the solution was within 0.1–0.2 of its standard
  error of the end point after the default 10 and the not-converged warning fired on every
  one. Now 8–25 passes, each stop within 0.1 of a standard error of the end point. A run never
  stops later than before, and a smaller `Tolerance` still tightens every bin.
- The not-converged warning names the bins that still move, with their redundancy against the
  median of their pool group, and how far the solution moved in the last pass in units of its
  standard error.
- `R.history` and `R.convergence()` carry `binRatio` (the bins' stop quantity, ≤ 1 = settled)
  and `xStep` (the solution's last step in standard errors); `maxDev` stays as a reading. The
  `OnPass` context carries both.
- `lsqvce.engine.psdVCE` returns each bin's deviation `.dev` and redundancy `.redundancy`.

## 6.25.2 — 2026-10-16 (WP-145: the examples converge, and the gate hears them)

Corrections to four examples and to the test tooling. Nothing in `+lsqvce` changes, so no
estimate of the toolbox changes.

- `demo_gravity`, `demo_grace`, `demo_noisePolicy` and `demo_sensors` stopped on their pass
  budget without converging (`lsqvce:adjustVCEjoint:notConverged`). Measured cause: spectral
  bins that the data barely determine — the zero-frequency bin that per-arc locals absorb, or
  60 bins on 60-sample arcs of white noise — settle slowly or decay toward zero, so the worst
  bin never meets the tolerance. They now use `Noise.Smooth = 1` (gravity, grace, noisePolicy)
  or `Noise.Grid = 8` (sensors) with the toolbox budget of 10 passes, and converge in 3–7 passes
  (`demo_sensors`' per-arc run in 14 of 20) on the same estimates. Their help states the measured values.
- `demo_sensors` estimated its per-sensor noise levels on data that still held the bad day it
  plants, and reported them as recovered while its table showed LRI and GPS 32 % off. It now
  finds the bad day first (one variance factor per sensor-day), leaves that day out, and shows
  every stream against its truth (all within 4 %).
- The test gate: a warning raised inside an example is now part of the gate's warning
  inventory (it used to be discarded with the example's printed text); the gate's log must be
  written by the run it judges; `dev/tests/run_lsqvce` finds the suites in `dev/tests` again,
  and `LogFile=''` switches its log off.

## 6.25.1 — 2026-10-16 (WP-142: the code review of WP-141)

Corrections to `prepareDay` (6.25.0). No estimate changes.

- A `prepareDay` that declares no argument (`@() load(...)`) is called with none; 6.25.0 passed
  it one and failed with "Too many input arguments".
- Its result is kept for the **current day only**, whatever `Options.CacheDays`: the engines
  read every stream of a day together, and `CacheDays` keeps the built streams. 6.25.0 kept
  `CacheDays` results, never read again once the day's streams were built — with
  `CacheDays = Inf` every day's result (whole-day partials, say) stayed in memory to the end.
- How often it runs, said exactly: once each time a day is built, for all its streams — once
  per run with `CacheDays = Inf`, twice per pass with the default `CacheDays = 1`, and on every
  parallel worker that builds the day. 6.25.0 said "once per day".
- A `buildDay` or `prepareDay` handle that does not resolve is `lsqvce:makeLoader:buildDay` /
  `lsqvce:makeLoader:prepareDay` with the day, like any error they throw (it escaped as a
  raw MATLAB error); with a `prepareDay`, an unresolvable `buildDay` is refused while the
  loader is built. Through `fit`, the refusals name the fields you wrote (`prepareDay`,
  `processingOptions`), not makeLoader's option names.
- Since 6.25.0, `template/kbr_variational`'s `kit.buildDay` called directly needs
  `procOpt.prepared` set to what `kit.prepareDay` returns for the day.
- Documentation: the API's buildDay chapter, the Manual (the template's design, the memory of
  large runs) and `Options.CacheDays` say the above; the identifier appendix names the
  `lsqvce:fit:<reason>` family relabelled from makeLoader; the README's theory reference is
  section 1.3.

## 6.25.0 — 2026-10-15 (WP-141: prepareDay, what all streams of a day share)

New, optional: `data.systems(k).prepareDay`. No estimate changes; nothing changes for a
declaration without it.

- `buildDay` runs once per (day, observation type), so work that several streams of a day
  need was repeated for each of them, and again in every pass — or cached by hand, as
  `template/kbr_variational` did. `prepareDay(day, common, local, procOpt, info)` runs
  **once per day**, before the first stream of that day is built; every `buildDay` of the day
  receives its result as `procOpt.prepared`. It stays resident as long as the day
  (`Options.CacheDays`, at least the current day): with `CacheDays = Inf`, once per day per run.
- Refused while building the loader: a `prepareDay` that is not a function handle
  (`lsqvce:fit:builderKind`), a `buildDay` with fewer than five arguments, which never sees
  `procOpt` (`lsqvce:makeLoader:prepareUnused`), and a `processingOptions` that already has a
  field `prepared` (`lsqvce:makeLoader:preparedField`) -- `fit` reports both under its
  own prefix, the original as the cause. A `prepareDay` that throws is
  `lsqvce:makeLoader:prepareDay`, naming the day, with your error as the cause.
- `template/kbr_variational` integrates the C/D pair in `prepareDay` and picks each stream in
  `buildDay`; its own day cache is gone.

## 6.24.2 — 2026-10-15 (WP-140: the observations of WP-139)

Documentation and test bookkeeping only: no estimate, interface or identifier changes.

- Help text that said what the code does not: `applyRegularisation`'s `'tsvd'` is the accumulator's
  truncation of the unfloored local block (since 6.20.0), not a floor; `mapreduce` says why
  `'parfor'` matches `'serial'` to 1e-14 rather than bitwise (single-threaded workers); the claims
  "4/4 injected outliers detected" (`outliers`) and "2 iterations from a 3x wrong weight" (`vce`)
  had nothing behind them and give way to what a test pins.
- Help structure: stamps that sat inside a help block (six files), one OPTIONS block in
  `arcVarianceFactors`, one INPUT block in `separableWhiten` and `varWhiten`, every public unit's
  H1 line in full (`%LSQVCE.ENGINE.RUNMETA`, ...).
- Tests: a check that compared CommonPrior's absorbed DoF with itself now compares it with a dense
  trace; copied or wrong tolerance comments say their own source; the checkPSD verdict test counts
  as a precision test. The R2 Python mirror's spectral check uses an absolute tolerance only.

## 6.24.1 — 2026-10-15 (WP-139: every unit states its ACCURACY)

Documentation only: no estimate, interface or identifier changes.

- The 48 public units that 6.24.0 froze without an ACCURACY block (audit check 34) state one
  now: the numerical claim, the test that pins it with its tolerance, and the reference it is
  checked against (a dense solve, a closed form, a Python mirror's recorded measurement), or
  "Not applicable" with the reason where a unit computes nothing (settings, keys, bookkeeping).
  No number was written without a source in a test, the unit's own help or the project record.
- The audit's list of frozen units is empty: every public unit must state its ACCURACY.

## 6.24.0 — 2026-10-15 (WP-138: the code review's judgement calls)

No estimate changes. Identifiers change where marked **ids**.

- **ids** A day list that is not numeric integers is `lsqvce:fit:badDays` (it was `noDays`, the
  identifier of a missing list). System ids of two kinds are `lsqvce:fit:mixedIdTypes` (it was
  `badId`, which keeps an id that is no number or single letter), the accumulator's reason name.
  The accumulator judges the kinds before its duplicate test: `addSystem(67)` after `'C'` was
  called a duplicate (`char(67)` is `'C'`); it is `lsqvce:accum:mixedIdTypes`.
- **ids** `gcvLambda`'s two errors sit under the segment of its warnings: `lsqvce:gcv:badData`,
  `lsqvce:gcv:nullUnobserved` (6.22.0–6.23.0: `lsqvce:gcvLambda:*`).
- **ids** `lsqvce.report` takes `Iteration.loop` — `"noise"` (the default), `"irls"`,
  `"gaussNewton"` — and owns the wording; the display label `Iteration.what` it branched on is
  retired. Another value, or `.what`, is refused (`lsqvce:report:badIteration`). `R.report()`
  passes the code.
- One definition each where the code had several: `lsqvce.util.priorStream` makes and
  recognises the prior pseudo-streams' observation type (nine sites spelled the prefix, and `fit`
  parsed stream keys for it); the parameter templates' validity rule (`declaration` and
  `makeLoader` each had one); the rules "ids of one kind, each day once" (`declaration` and the
  accumulator); the verdict's redundancy formula (four copies) and its list of prior batches.
  `R.snoop` re-diagnoses with the thresholds of the `R.diag` it was handed, read rather than
  restated. The contamination gauge has its own file.
- A local prior batch with a full weight matrix — reachable only by a direct caller of
  `lsqvce.engine.verdict.diagnose` — is judged by the generalised redundancy; it was skipped
  without a word and its datum row counted as uncontrolled.
- Every method of `lsqvce.engine.verdict` states its ACCURACY, `ftest`'s p-value and critical
  value now pinned against the closed form of F(2, r); so do the other units this release
  touched. The audit freezes the units still without one (check 34): a new unit must state its
  accuracy.

## 6.23.0 — 2026-10-15 (WP-137: the code review of WP-135)

Only `R.snoop` behaves differently (**verdict**); no estimate changes.

- **`R.snoop` warns at the contamination level `indeterminate` (verdict).** WP-135 added the
  level — the residuals are far larger than sigma0 promised, the robust scale sees no excess of
  outliers, and below n ≈ 450 a sigma0 declared too small cannot be told from a robust scale
  captured by more than half the data — but `R.snoop` had no case for it and ran on in silence.
  Measured on the mirror: with 55 % blunders, n = 30…200 lands there in 100 % of draws, so a
  real breakdown went unchallenged. It now warns (`lsqvce:snoop:indeterminate`), names both
  causes, and proceeds under its `MaxRejectFrac` cap (user decision: warn, not refuse — a clean
  run with a mis-declared sigma0 lands in the same level). The report's `SNOOP_INDETERMINATE`
  says so.
- **Tests for what WP-135 changed without one:** `R.snoop` refusing at `breakdown`
  (`lsqvce:snoop:breakdown`, never fired by a test before) and proceeding under `Force` with a
  warning; snoop's re-solve keeping a regularised local block's `Reg`; snoop testing only the
  rows `R.diag` tests (`RMin`); `R.compare` reading the solve's rank and redundancy;
  `Noise.Whitener` in any case on the noise path; the advice of an unsettled Gauss-Newton run.
- `gcvLambda`'s ACCURACY states the measured and the worst-case error and cites its authority;
  its test no longer mutes warnings that could not fire, and its tolerance sits at its stated
  source. The contamination mirror is the code's ladder in full.

## 6.22.1 — 2026-10-15 (WP-136: where `buildDay` goes)

No change to any computation.

- **The API reference, chapter `buildDay`, says where the handle goes** (*Where the handle goes —
  one per system*): every entry of `data.systems` carries its own `buildDay`; several readers in
  one adjustment are several systems; one reader for several setups takes the setup as data,
  through `processingOptions` (its 5th argument) or a closure; data sets that share no common
  parameters are separate `fit` calls. The chapter gave the signature only, so a reader looking
  for the handle next to `Options.OnPass` did not find it (user question, 2026-10-15).
- `fit`'s help and the Manual (*Declaring your problem*, *Several systems*) say the same in one
  sentence each.

## 6.22.0 — 2026-10-14 (WP-135: the front end, the loader and the diagnostics of the third whole-code review)

Numbers change for the cases marked **numbers**; the diagnostics' verdict changes where
marked **verdict**.

- **A repeated day is refused.** `data.days = [1:15, 15:30]` (or an overlap in a system's
  `days`) accumulated day 15 twice, silently: x pulled toward it, σ too small, dof inflated
  (`lsqvce:fit:duplicateDay`; the accumulator refuses it too, `lsqvce:accum:duplicateDay`).
- **System ids are of one kind.** A character id next to a number (an index default is a
  number) made the accumulator's id list a character row: robust factors went under keys the
  loader never read and were reported as applied (`lsqvce:fit:badId`,
  `lsqvce:accum:mixedIdTypes`). Ids `'i'` / `'j'` stay characters (they parsed as 1i).
- **`lsqvce.layout` works with a system without `localParams`** (it stopped; the Manual's
  own example).
- **Local prior groups follow the common prior's rules (verdict).** A local group's row the
  data cannot see is a datum row — not tested or counted (it made the verdict SUSPECT) — and
  `R.snoop()` never rejects a prior row, common or local. `R.systematic()` leaves the priors'
  rows out of its sequences; a covariate may be given over the data rows only.
- **The contamination ladder is calibrated for every n (verdict).** Its null band was
  measured at n ≈ 300; at n = 30 clean data raised `SNOOP_CONTAMINATED` / `SNOOP_REFUSE`
  (ERRORs) in 3 % of fits, and a σ0 declared 3.5× too small read "breakdown" in half the
  fits, with `R.snoop` refusing. The band now follows the measured null law of `fhat` at
  every n (false alarms ≤ 0.25 %), "breakdown" needs `fhat` implausibly low for clean data,
  and where n is too small to tell a captured robust scale from a σ0 declared too small the
  report says so (`SNOOP_INDETERMINATE`, a WARN). The price: below n ≈ 200, a 10 %
  contamination is no longer flagged by the gauge alone (the outlier test still sees the
  blunders themselves).
- **GCV is the GCV of the data (numbers).** `Reg = "tikhonov"` with `RegLambda = 0` chose
  λ from a functional of the normal equations, not of the data (‖y − y_true‖ up to 2770
  against 1.0); a `RegMatrix` with a null space (a first difference) put every candidate λ
  above 10³‖N‖. Both fixed; `gcvLambda` takes m and l'l.
- **Diagnostics in the far tail.** Below a Bonferroni level of ~1e-16 the critical values
  came back Inf / NaN and the test switched itself off; they are formed from the upper tail
  now (`lsqvce.kernel.stat.tCrit`). `SNOOP_VACUOUS` fires when a residual must carry ≥ 90 %
  of v'Pv to be flagged (the old criterion could never fire). Under `Test = "mad"` the MDB
  is labelled the Baarda (known-σ) value.
- **A joint robust run's report quotes `R.sigma0` (numbers).** Its s0, χ² test and scaled
  parameter column read the hard solution's vPv without the truncation correction.
- **Smaller:** a throwing `buildDay` is reported as `lsqvce:makeLoader:buildDay` with the
  day and stream (documented, never raised; a toolbox `lsqvce:*` error raised inside it
  passes unchanged); a local prior group's `from` and name are
  validated; a "scale" group activated on the last pass is scaled in the closing solve; a
  system without own days claims only `data.days`; `Noise.Whitener` is validated in any
  case (`lsqvce:Options:badWhitener`); validity entries of another type and static arc-id
  vectors with repeats are refused; `synthLoader` refuses robust factors it cannot apply;
  `R.convergence().converged` reads `R.converged`; an unsettled IRLS run is advised on
  `RobustOpts`; `R.compare` uses the solve's rank and redundancy; snoop's re-diagnosis keeps
  `Store` / `NWorst` / `RobustCheck` and its re-solve the solve's `Reg`;
  `outliers.mask` is a logical mask.

## 6.21.0 — 2026-10-13 (WP-134: the noise estimator and the kernel of the third whole-code review)

Numbers change for the cases marked **numbers**.

- **A common prior with `estimate = "scale"` is in the stop rule (numbers).** The joint loop
  stopped on the law and the levels alone, and a run that settled on pass 1 — before the
  group's first update — returned the declared σ unestimated (σ 100 against 0.92, the
  solution 154 against 5.4 off in the review's case) with `R.converged` true. The loop now
  waits for the group's update and converges only when its change is below `Tolerance`
  (`R.history.priorDev`).
- **The noise estimator no longer squares the normals.** Its probe used chol of the reduced
  normal matrix and refused, as `lsqvce:adjustVCEjoint:notPD`, problems the QR solve takes
  (cond(A) from ~3e8); where chol succeeded the probe term was 4–25 % off. It uses the
  solve's own factor now (numbers on ill-conditioned problems).
- **A long arc is refused under `Estimate = "noise"`, not run out of memory.** The estimator
  forms each arc's covariance densely whatever `Noise.Whitener` says; the `maxDenseArc`
  limit (20 000 epochs by default) now guards that too (`lsqvce:covModel:denseArcTooLong`).
  The Manual's cost claim for `"ar"` is corrected: it speeds up the solve, not the estimator.
- **One blunder no longer switches off robustness (numbers).** A single residual beyond
  ~4.5e12 MAD switched the robust scale to the standard deviation it inflates (Huber's
  estimate moved by 2e8, 8σ rows kept full weight). Tuning constants are validated.
- **`SkipThreshold` is per stream**, as documented: a small stream can no longer be skipped
  whole among large ones and keep its declared level silently.
- **An arc without redundancy keeps its factor (numbers).** A probe-estimated redundancy ≤ 0
  (short arcs with locals) gave the arc 1e8 times the weight of the others; it is held now
  (`lsqvce:arcVarianceFactors:held`).
- **A `var` law of order 1 can be declared** (`ar` d×d was refused); a VAR(0) is white noise.
- **Smaller:** `armaACV` follows the impulse response to its decay for any lag count (a
  near-unit root came out up to 27 % low, depending on n; beyond 2^20 samples the
  truncation is reported); an `acf` law honours a per-row σ (it was averaged); an uneven
  epoch grid under a `psd` / `arma` law is reported (`lsqvce:covModel:unevenGrid`);
  `crossCorr` refuses a partial epoch; `ar2arma` normalises its input; the robust joint
  sweep leaves local prior rows out; `Sigma2_0` values are validated; `AvgTail` averages
  every pool group's law; `jointModelOpt` carries the pass and the robust factors for engine
  callers; `streamContribution` returns its probe seed and no standard error for one probe;
  `redundancy`'s Hutchinson check compares the per-group traces.

## 6.20.0 — 2026-10-13 (WP-133: the solve core of the third whole-code review)

Numbers change for the cases marked **numbers**.

- **A row of weight 0 is not an observation (numbers).** A zero data weight (`data.p`, a
  batch source's `p`), a robust rejection (`Robust` "igg3" / "bisquare", factor 0) and a
  rejected epoch of the joint estimator deleted the row's information but kept it in the
  degrees of freedom — σ₀ 0.753 against 0.949 with 10 of 30 rows at weight 0 — and the
  verdict called it "verified". Such a row now keeps its index but leaves the dof, Σr, the
  tests and Bonferroni's n; the verdict counts it in its own quadrant, **"rejected"**
  (`R.diag.quad.rejected`, `R.diag.excludedIdx`, `R.reliability()`). Weights must be ≥ 0 and
  finite (`lsqvce:batchSystem:weightMatrix`); a day may mark rows with `.excluded`.
- **A robust fit's diagnostics replay the returned solve.** At a `maxPasses` exit the
  verdict, `residuals` and `snoop` replayed the next pass's weights (the factor map was
  shared and overwritten).
- **Regularised local blocks count their effective redundancy (numbers).** `Reg` "tsvd" now
  truncates (the dof counts the kept rank, as its help promised); "tikhonov" / "lm" count
  m − u + Σ tr(Q_ll Λ) and the data's residual mass without the penalty y'Λy. The verdict's
  Σr identity read INVALID for every regularised solve.
- **Gauss-Newton:** `R.residuals()` measures the residuals at the estimate (it took the
  absolute x as a step: σ₀ off by orders); a trial point where the model is non-finite is a
  rejected step, not an abort (`lsqvce:adjustNonlinear:nonFinite` only at x0).
- **The joint robust hard solution:** its σ₀ correction (truncation) reaches `residuals`,
  which no longer warns `sigma0Mismatch` on every igg3 / bisquare noise run.
- **Smaller:** a sparse local design no longer stops the rank decision or `snoop`'s
  downdate; a local block the data never see is a zero block (it gave NaN locals); a prior
  that leaves no redundancy is refused (`lsqvce:accum:noRedundancy`, it gave Cxx = 0); an
  all-zero prior row is refused (`lsqvce:prior:zeroRow`); a `LocalBlockSize` whose blocks are
  coupled is warned and solved whole (`lsqvce:accum:notBlockDiagonal`).

## 6.19.0 — 2026-10-13 (WP-132: the five recorded lows of the second whole-code review)

Numbers change for the cases marked **numbers**.

- **`chi2Inv` converges in the far lower tail (numbers).** Below p ≈ 1e-40 the iteration ran
  out of steps and returned wherever it stopped (r = 1, p = 1e-60: 9e31 times too large). It
  now starts from the small-x expansion where that is closer, iterates on log F below the
  median, and brackets its steps; a quantile below `realmin` comes back as 0 or a subnormal.
  The quantiles the toolbox uses (p ≥ 0.005) move by at most 1e-14.
- **Robust IRLS returns the factors its solution used.** `R.history(end).factor` was the
  next pass's factor, one step ahead of the returned x, Cxx and diagnostics. It is now the
  weights the returned solve used, so R.x is the weighted solution on them at every stop.
- **The local covariances survive a refused v'Pv.** After the normal backend refused a
  cancellation-destroyed v'Pv, `residuals` repaired σ₀ but returned every `Cyy` as NaN; they
  now carry the repaired σ₀² and the solve's cofactor (`BlockAccumulator.commonCofactor`).
- **`Sigma2_0` weights the first pass per arc** (engine tier, `adjustVCEjoint`). Pass 1 was
  weighted with the mean of the start factors while the update multiplied each arc's own; the
  arcs are now enumerated and pinned first (one extra accumulation when `Sigma2_0` is given).
- **The arc factor takes the whole law as its component, for any law (numbers).** It assumed
  the arc's covariance Toeplitz; for a non-stationary law (a per-row sigma, say) that is a
  different, less efficient component (sd 0.62 against 0.32 on a 1..10 sigma ramp). The pair
  is now in closed form (`redundancyLags` `.sLaw` / `.omegaLaw`); Toeplitz laws are unchanged
  to rounding.

## 6.18.0 — 2026-10-13 (WP-131: noise model, kernel and contributions of the second whole-code review)

Numbers change for the cases marked **numbers**.

- **Level-2 contributions on the Level-1 split (numbers).** `contribution(..., 'Level', 2)`
  split the reduced normal matrix differently from Level 1 whenever a system's locals spanned
  several groups (`GroupBy` "obsType" or "day"): a group's diagonal did not sum to its share
  and could be negative. Both levels now use the leverage split N_g = Z_g'Z_g (user decision).
- **A declared PSD law with a physical frequency axis works under `Estimate = "noise"`.** The
  estimator's law (on normalised frequency) kept the declared axis and failed on pass 1, or at
  equal lengths integrated on the wrong axis.
- **`Noise.Smooth > 0` can converge.** The stop test read the unsmoothed factors, which are
  not 1 at a smoothed fixed point; it now reads the change of the spectrum's shape.
- **The lag sums of a non-stationary covariance are exact (numbers).** A per-row sigma, a
  down-weighted epoch or a gap after mid-arc was taken for a Toeplitz covariance.
- **`armaACV` of a sparse AR polynomial (numbers).** A polynomial such as `[1 0 0.998]` or a
  seasonal one stopped its impulse response early (r(0) 36 % low).
- **`crossCorr`'s verdict is scale-free.** It tested covariances against a correlation band.
- **A non-robust `Estimate = "noise"` run is no longer reported as "a robust (none) fit".**
- **A prior beyond the last pass is reported** (`lsqvce:fit:priorNotActive`) instead of being
  dropped silently when the joint loop ends before its `from`.
- Smaller fixes: the MDB of a correlated batch (direct `diagnose` calls, **numbers**);
  `residuals(..., 'Whitened', false)` keeps the whitened σ₀; `outliers` under `'spmd'`;
  `fit`'s Name=Value groups override the given Options field by field; `toTable` without
  `localParams`; a one-argument arc resolver with day-scoped common parameters.
- **A data stream named `prior_*` is refused** (`lsqvce:makeLoader:priorType`): the name is
  reserved for the prior pseudo-streams, and such a stream was dropped silently (user decision).
- **`R.snoop` tests at `R.diag`'s Alpha, Bonferroni and Test** unless told otherwise (user
  decision): an observation `R.diag` flagged at Alpha 0.01 could pass `R.snoop` as clean.
- Tests: 673 (+11).

## 6.17.0 — 2026-10-13 (WP-130: the solve core and the diagnostics of the second whole-code review)

Numbers change for the cases marked **numbers**.

- **A rank-deficient local block is solved on its factor (numbers).** One accidentally
  duplicated local column made the QR pivots fail, and the fallback squared the block: the
  pseudo-inverse cut a weak but determined direction together with the null one (x 31 formal
  sd off for a quadratic in raw seconds), the degrees of freedom counted every local, and the
  local solution and covariance were solved without the pseudo-inverse. The rank is now
  decided on the factor; the solution, σ₀, the local covariance, the diagnostics and snoop all
  use the same pseudo-inverse, and the degrees of freedom count the rank.
- **Gauss-Newton reports the point it returns (numbers).** σ₀, the covariance and the
  diagnostics belonged to the linear prediction at x + dx, even when that step was rejected
  (σ₀ 0.19 against 0.77: formal errors up to 4× too small). They are now evaluated at the
  returned x; `R.diag.x` is the estimate (the report's x column and the ellipse centre read
  the final increment, ~0); a diverged run says `stopReason = "diverged"`.
- **The datum test does not depend on parameter units (numbers).** A parameter in small units
  with a prior group of its own was taken for a datum and reported sd 0.
- **`Test = "mad"` leaves the datum rows out of its scale (numbers).** Rounding gave them
  near-zero scores about half the time, which pulled the MAD down.
- **`R.snoop` after rejection.** With `Reaccumulate = false` it returned the input `D` and the
  pre-rejection σ₀; both now describe the survivors. Its Bonferroni count leaves the datum
  rows out, as `diagnose` does; `REJ.history` carries the critical value.
- Tests: 662 (+4).

## 6.16.1 — 2026-10-13 (WP-129: snoop's downdate of a local block on its QR factor)

- **`R.snoop` keeps a badly conditioned local block exact after a rejection (numbers).**
  6.16.0 built the diagnostics cofactor of a local block from its QR factor, but snoop's
  rank-one downdate rebuilt that block from the squared `Nll` again: from cond(B) ≈ 7e7 one
  rejection lost the weak direction (a quadratic local polynomial in raw seconds over 10⁴ s:
  the leverage sum over the survivors 6 instead of 7). The downdate now runs on the factor.
- Tests: 658 (+1).

## 6.16.0 — 2026-10-13 (WP-128: the remaining findings of the whole-code review)

Numbers change for `R.snoop` after a rejection, for the MAD-scaled test, and for the
diagnostics of a local block with a badly conditioned design (marked **numbers**).

- **`R.snoop` keeps a prior a prior (numbers).** The clean re-solve after a rejection put the
  prior's rows in as data, so a datum's covariance grew with its σ (sigX 250 instead of
  0.164 at σ = 1e3 on a levelling net). The prior's rows are never rejection candidates.
- **`R.snoop` returns one index space.** The re-diagnosis numbered the surviving rows, while
  `D.rejectedIdx` numbered the original ones. All index fields of the returned `D` are now
  original observation numbers; `D.obsIdx` names the original number of each per-observation
  entry.
- **The MAD-scaled test uses its own critical value (numbers).** Pope's τ holds only when
  the scale contains the tested residual. The MAD test now uses Student t with r/2 degrees
  of freedom (`lsqvce.kernel.stat.madCrit`), calibrated to the nominal family-wise level
  (user decision). The τ value gave 63 times the nominal false-alarm rate at n = 30, 5.3
  times at n = 200. On short networks under heavy contamination one pass finds fewer
  blunders than before (n = 50, 10 % at 25σ: 60 % of runs find ≥ 4 of 5, against 93 %).
- **The diagnostics of a badly conditioned local block are exact (numbers).** A local block
  the QR backend eliminates exactly (cond(B) up to ~4.5e12) was re-inverted from the squared
  `Nll` for the cofactor, which cut eigenvalues from cond(B) ≈ 7e7 on: a quadratic local
  polynomial in raw seconds over 10⁴ s lost one redundancy unit. The cofactor, the local
  solution and the local covariance now come from the QR factor.
- **A late local prior is reported in a single solve.** `lsqvce:fit:priorNotActive` counted
  `data.priors` only; a local group with `from > 1` vanished silently.
- Tests: 657 (+3).

## 6.15.0 — 2026-10-12 (WP-127: findings of the whole-code review)

Numbers change for ARMA laws with an MA part, for split datum groups and for the report
of any fit with a datum prior (marked **numbers**).

- **The ARMA whitener with an MA part is exact (numbers).** The recursion 1/B(z) carries
  state; only the first p rows were replaced by the exact start, so the rows after them were
  not chol(C)⁻¹ (F'F off inv(C) by 1.6 % at a = [1 −0.5], b = [1 0.4]; whitened variance 10.7
  instead of 1 at row 2 for b = [1 0.8]). The exact start now spans the filter memory and the
  recursion continues from its exact innovations: chol(C)⁻¹ to 3e-13.
- **A datum may span several prior groups (numbers).** Translation and rotation declared as
  two groups were both refused as datum (they overlap off-centre), and the covariance grew
  with the datum σ. Groups the data cannot see now need orthogonality only to the groups it
  can see.
- **Datum rows are not observations (numbers, report).** Their redundancy is 0 by
  construction; they were counted as uncontrolled observations with MDB = ∞, so every free
  network's report said SUSPECT. The diagnostics keep them in Σr and vPv but do not test or
  count them; new fields `D.nTests`, `D.datumRows`; quadrant `"datum"` in `R.reliability()`.
- **The joint estimator runs with a local prior that starts late** (`from ≥ 2`). Its
  pseudo-arc was a "gap" skip before `from` (counted against `SkipThreshold`) and one more
  stream after it, which broke the per-pass record. Inactive arcs are no longer skips; the
  record has one column per stream of pass 1.
- **Gauss–Newton checks its final step.** A converged run applied one undamped step without
  evaluating it; near a plateau that took F from 28.8 to 49.6 (a logistic fit). The step is
  kept only when F does not rise. Deeper on the plateau the undamped linearisation is singular
  and the run is refused by name (`lsqvce:qr:singularR`), as before.
- **A NaN or Inf in A or B is refused** (`lsqvce:accum:arcNotFinite`), as one in l already was.
- `.gitignore` covers `dev/references/local/` (the folder moved to `dev/`) and every
  `dev/validation/*_measure.txt`.
- Tests: 654 (+3).

## 6.14.0 — 2026-10-12 (WP-126: the review findings on WP-125)

Numbers change for `Noise.Whitener = "arma"` with an MA order (marked **numbers**).

- **A dense law refuses an arc above `maxDenseArc` epochs.** Since 6.13.0 an `acf`, `custom` or
  Toeplitz `psd` law is factorised over the whole arc (O(m³), 8·m² bytes) and the default
  `Noise.Whitener = "toeplitz"` reaches declared `psd` laws, so a long arc could exhaust memory.
  Above 20 000 epochs — the largest dense factor the toolbox formed before 6.13.0 — it is now
  refused (`lsqvce:covModel:denseArcTooLong`) with the routes out: `Whitener = "ar"`, shorter
  arcs, or `covariance.<type>.maxDenseArc` raised deliberately.
- **`lsqvce.kernel.ar2arma` fits the AR part forwards (numbers).** The Hankel regression ran
  the recursion backwards and returned the reciprocal roots (0.94 → 1.064), so the stationarity
  guard fired on almost every input and `Whitener = "arma"` silently used the long AR. The true
  ARMA is now recovered (1e-14 on the round-trip test).
- **A non-invertible MA part falls back to the AR** with the warning
  `lsqvce:ar2arma:nonInvertibleMA`, as a non-stationary AR part already did, instead of
  failing the run in `covModel` (`maNotInvertible` stays for a declared `arma` law).
- **`arcVarianceFactors` returns `.pooled`,** the pooled Helmert factor that is the joint
  estimator's stream level since 6.13.0 — now computed in one place and unit-tested.
- Tests: 651 (+4).

## 6.13.0 — 2026-10-12 (WP-125: the code findings of the documentation review)

Numbers change for three of the items below (marked **numbers**).

- **Dense laws are whitened over the whole arc (numbers).** An `acf`, `custom` or Toeplitz
  `psd` law was factorised per batch of 20 000 rows, so an arc longer than one batch lost the
  correlation between batches. The arc is now factorised whole (exact; O(m³), m² memory).
  A filter law whose overlap does not fit in a batch is also whitened as one batch.
- **`Noise.Whitener` governs every `psd` law (numbers).** A declared `psd` law was whitened by
  a fixed AR(30) filter in single solves and the stand-alone IRLS; it now follows
  `Noise.Whitener` like the estimated law (default `"toeplitz"`: exact).
- **The level of a stream without arc factors is the pooled Helmert factor (numbers),**
  sigma² · Σq / Σt over its arcs, instead of the mean of the arcs' ratios. Measured in
  dev/validation/wp125_level_mirror.py: on mixed arc lengths the scatter drops from 0.047 to
  0.010, and the bias from probe noise from +1 % to +0.1 %; equal arcs are unchanged.
- **ARMA laws with an MA part overlap their batches by the decay of 1/B(z).** The overlap used
  to count the AR order only, so the restart transient of the recursive part leaked into every
  batch after the first. A non-invertible MA part is refused
  (`lsqvce:covModel:maNotInvertible`).
- **`R.contribution()` shares the common solution by its trace** t_g = tr(N⁻¹N_g), computed
  exactly in the redundancy sweep (`red.group.tcc`). It used to divide the groups'
  REDUNDANCIES, which rank the groups nearly upside down on the information they carry.
- **The report:** the contamination cross-check runs whatever the `Bonferroni` setting; the
  PASS text no longer promises "half a sigma"; the remedies name `Backend = "qr"`; an IRLS or
  Gauss–Newton run reports its own pass count and names its loop.
- **`R.convergence()`** adds the two stop-rule quantities `maxDev` and `levelDev`;
  `maxSigmaOffset` is described as what it is (the stream levels, not the stop rule).

## 6.12.2 — 2026-10-12 (WP-124: three documents, each with one job)

No change to any computation.

- **The Manual** (`doc lsqvce`, `lsqvce_manual.pdf`) is rewritten as a step-by-step user guide:
  purpose, installation, the workflow, declaring your problem, choosing the stochastic model,
  priors, reading the results and the report (every level, every finding, what to do), watching a
  run, the examples (purpose, logic, outline and implementation of each, and the templates), and
  large problems. It no longer carries version history or development notes.
- **The API reference** (`html/lsqvce_api.html`, `lsqvce_api.pdf`) is new: the `data` struct,
  `buildDay`, covariance specifications, `Options`, priors, `Results`, the help of the user-tier
  units, the engine / kernel / util units and the identifier catalogue. The interface tables that
  were in the Manual moved here.
- **The theory** (`lsqvce_theory.pdf`) now states only the mathematical framework, matched to the
  code: square-root information as the default backend, elimination per system, the datum per
  prior group, the hat-spectrum components, absolute arc factors, the Koch-convention robust
  nugget with the truncation correction, the exact two-sided MDB. The review and redesign
  chapters of September 2026 are archived under `dev/reviews/R5_2026-09-01/`.
- The help blocks of `fit`, `Options`, `Results`, `report` and `layout` describe the current
  interface only, and every pointer to "the Manual, chapter 3" now names the API reference.

## 6.12.1 — 2026-10-12 (WP-123: the kernel whiteners take a sparse design matrix)

- **A sparse `B` no longer breaks a whitened stream.** The Manual (chapter 3) allows a sparse
  local design matrix, and `synthLoader` recommends it, but the kernel whiteners reshape to N-D
  arrays or call `filter`, which MATLAB refuses for sparse input: a multi-axis stream with
  `Cep` stopped in `lsqvce.kernel.separableWhiten` with `MATLAB:reshape:sparseNDmatrix`
  (found by the WP-122 template). `separableWhiten`, `varWhiten` and `armaWhiten` now run
  their dense filter on the non-zero columns of a sparse block only and return it sparse, at
  its full width (memory: batch rows x active columns). Dense input takes the unchanged path.
- `template/kbr_variational.m` returns its position streams' `B` sparse again (6.12.0 had to
  return it dense).

## 6.12.0 — 2026-10-12 (WP-122: a variational-equation range-rate builder for real data)

- **`template/kbr_variational.m`** — a runnable `buildDay` for the GRACE-FO C/D pair:
  - Reads a day of Level-1B files (GNV1B, KBR1B, ACT1B, SCA1B; the v04 ASCII layouts, all
    overridable) and a reference field from an ICGEM `.gfc`.
  - Integrates every arc's orbit together with its variational equations
    Psi = d[r v]/d[r0 v0 | C_nm S_nm | ACC bias | ACC scale] (Adams–Bashforth–Moulton 8,
    exact coefficients).
  - Returns the range-rate rows d rho_dot = [h' e'] · (Psi_D - Psi_C) and the position rows as
    `lsqvce` streams, with arc-wise initial states and KBR empiricals and day-wise ACC calibration.
  - Acceleration, gravity gradient and da/dC_nm are analytic: normalized solid harmonics by
    the Cartesian (pole-safe) recursion, and exact sparse derivative operators degree n -> n+1.
  - Every column carries a declared unit (1e-10 for C_nm, S_nm, ...), so the QR rank decision
    sees O(1) columns; physical value = unit .* estimate (`kit.unitCommon .* R.x`).
  - Validated against oracles not built here (a Python numerical twin in dev/validation and
    four new gate tests): Legendre functions, closed-form J2, Kepler, re-integration, finite
    differences, and a closed loop through `lsqvce.fit` on synthetic L1B files.
  - What it does not model is stated in its help: the Earth-rotation stand-in (use the IERS
    2010 transformation on real data), forces beyond the static field and the accelerometer
    (`cfg.background`), and the product conventions to check against the handbook.

## 6.11.1 — 2026-10-12 (WP-121: diagnose reads the data once)

- **`diagnose` sweeps the observations once** (the residue §12.48 of review R7 recorded).
  The robust check (the MAD scale and the contamination gauge `fhat`, on by default) had
  its own sweep: every batch was re-read through the loader and its O(nc²) `halfRows`
  computed twice. Its raw scores (P v)ᵢ / sqrt((P Qvv P)ᵢᵢ) are now formed in the main pass
  from the same quantities -- one read and one `halfRows` per row, results bitwise
  unchanged. Only `'Test', "mad"`, which needs the robust scale before it can classify,
  still sweeps first.
- Corrected record: snoop's clean re-solve (`Reaccumulate=true`) has rebuilt the
  block-arrow system since 6.4.0 (WP-86 part 3); §12.48 had listed it as still dense.

## 6.11.0 — 2026-10-12 (WP-119: a datum beside a regulariser)

- **The datum test runs per declared prior group** (the residue §12.48 of review R7
  recorded). In 6.10.0 the stacked prior was a datum only as a whole, so minimum
  constraints declared beside a regularising group kept a σ-dependent covariance:
  se 0.09 / 41.6 / 4.2e4 on a levelled height at datum σ 1 / 1e3 / 1e6, with x unchanged.
  - A group (`data.priors(k)`) is a datum when the data and every other group are blind
    to its rows. `R.sol.Cxx` and the diagnostics' parameter cofactor are then the
    minimum-constraint ones of the data plus the other groups, s0²·pinv(N + Ki'Ki), at
    every σ, on both backends.
  - A group only partly in the null space, or one another group overlaps (a ridge on every
    parameter beside a datum), is information: its x moves with its σ, and it keeps the
    solve's covariance. Declare datum rows as a group of their own.
  - `lsqvce.engine.CommonPrior` carries the group of each row (`group`, 4th argument,
    default one group; `lsqvce:prior:badGroup`). Gauss–Newton's damping rows are a group
    of their own.

## 6.10.0 — 2026-10-11 (WP-114: priors without a second path)

- **A datum's covariance no longer depends on the datum's σ** (review R7, R7-05). x was
  already the minimum-constraint solution for any weight. The covariance, however, was
  s0²·inv(N + BB'/σ²), so it carried σ² on the datum directions: se 0.005 at σ = 1 and
  2.48 at σ = 1e3 on the levelling network.
  - `BlockAccumulator` now recognises a datum: a prior whose directions the data do not
    see, |N K'| ≤ 1e-10 |N||K|.
  - For a datum, `sol.Cxx` and the parameter cofactor of the diagnostics (`D.sigX`,
    `D.Q.block`, correlations) are the minimum-constraint ones, s0²·pinv(N).
  - They are computed stably: the datum rows are re-weighted to the data's scale, then
    inv(Rm'Rm) − Dm·Dm'. Subtracting the σ² term from the solve's cofactor would lose
    11 digits at σ = 1e6.
  - Leverage keeps the solve's cofactor, so the datum rows stay absorbed. Data rows are
    unaffected.
  - A prior only partly in the null space is genuine information and keeps the solve's
    covariance.
- **`residuals()` and `redundancy()` count a common prior** (R7-15).
  - `residuals()` adds the prior's residual to vtv and its rows to dof. New field:
    `res.prior`.
  - `redundancy()` closes the identity as sum(r_i) + r_prior = m + nk − u. New field:
    `red.priorRedundancy`.
  - Before, every fit with `data.priors` warned `lsqvce:residuals:sigma0Mismatch`, which
    blamed the loader, and `lsqvce:redundancy:identityViolated`. The wrong dof was marked
    authoritative.

## 6.9.0 — 2026-10-11 (WP-113: Gauss–Newton damps the parameters and keeps the prior on x)

- **The damping now acts on the common block** (review R7, R7-14). `adjustNonlinear` used to
  damp through the accumulator's `Reg "lm"`. That option acts on a local block, and a batch
  source has none, so the damping did nothing:
  - a rejected step was re-solved unchanged until `MaxIter`;
  - `R.x` was then the last accepted point plus that rejected step (F 36.8 against the best
    accepted 4.9 on the review's exponential fit).

  What changed:
  - The step is now one solve with Marquardt rows `sqrt(λ·d_j)·dx_j ≈ 0`, d = diag(J'PJ).
  - A rejection re-solves on the same accumulation with ten times the damping.
  - The descent check runs from the first pass, on the declared objective with the prior
    included.
  - An unconverged run returns the best accepted point.
- **A prior on a nonlinear problem acts on x, not on the step** (R7-04). Each step now
  solves `K·dx ≈ x0 − K·x`. The old form `K·dx ≈ x0` converged away from the declared
  optimum: b = +0.60 against +0.30 with σ = 1e-6, and the unconstrained answer when x0 = 0.
- `Lambda0 = 0` now means "start undamped": a rejected step still turns the damping on,
  starting at 1e-3.
- Mirror (the R7 review's mirrors, `wp113`):
  - the exponential fit reaches the optimum from [10, −5] and from [0.1, 1], with 30 passes;
  - Rosenbrock converges from five starts, the classic one included;
  - the prior fixtures land on the declared optimum to 1e-7.
- `R.history(k).vPv` is now F(x) with the prior term included.

## 6.8.0 — 2026-10-11 (WP-112: the QR backend's accumulator path is square-root throughout)

- **The default backend no longer squares the condition number on the accumulator path**
  (review R7, R7-03). `BlockAccumulator` used to read each system's QR factor back as
  N = R1'R1, Schur-eliminate the locals in normal form, and re-root the contribution with
  `chol` in `merge`. It did so also for systems without locals, and the prior-folded solve
  did the same.
  - Now the locals are eliminated on the factor: `QRSystem.reducedRows` re-triangularises
    with the local columns first.
  - The reduced contribution is merged as square-root rows (`ReducedContribution.Rows`).
  - A common prior is stacked as rows (K/σ, x0/σ).
  - Measured in a mirror at cond(A) = 1e7–1e9, with and without locals: x_c is 1e-8–1e-4
    formal sd from the dense `[A B] \ l`, against a median of 2–3 sd before. No full-rank
    problem is refused; before, 4–11 of 20 were refused at 1e8.5–1e9.
- **The rank decision is the code's own** (R7-16): `QRSystem.RankTol` = 1e3·eps on the
  pivots of the final factor.
  - Every free levelling network is refused (`lsqvce:qr:singularR`), and every one with a
    datum is solved. Before, `chol` decided by rounding: 61 % refused, 39 % solved with an
    arbitrary datum and sd ~1e7.
  - Problems with cond(A) above about 4e12 are now refused.
- At the WP-49/50 deep regime the QR path now returns an honest sigma0 at the fold (no
  `lsqvce:qr:cancellation`): the residual row rides along with the merged rows. A
  regularised local block still folds in normal form, with the WP-50 guard.
- Numbers on the QR backend change in the last digits for well-conditioned problems, and
  substantially for ill-conditioned ones. The normal backend is unchanged.

## 6.7.0 — 2026-10-11 (WP-111 part 3: the joint law as positive components)

- **The joint estimator's law is a non-negative spectrum on `Noise.Grid` bins** (review
  R7, R7-02 and R7-20; user decision 2026-10-11). The spectrum is linear between the bins
  and its covariance at every lag of every arc is the exact integral r = K·S, with K from
  the new kernel `lsqvce.kernel.hatACV`; its projection form `hatACV(L, m, X)` gives
  K'·X by one FFT, without forming K. Each bin is a variance component whose matrix is
  the covariance of a "hat" spectrum, which is positive semi-definite. Its Helmert factor
  (the projections of the arc's lag sums over **all** its lags, summed per pool group) is
  therefore never negative, and the law needs no floor and no zero-frequency rule.
  - Ellmer's components were the columns of the inverse type-I DCT of a truncated ACF.
    They are indefinite, so a bin's factor could go negative, and a floored bin left the
    declared covariance near-singular.
  - Measured on the runner (HANDOVER §12.44), on AR(1) streams at φ 0.8 / 0.9:
    - the shipped estimator's level ran away to 5e7 × the truth, with a negative lag-1
      correlation and coefficient z-scores up to 9.4;
    - at φ 0.5 it never converged;
    - on the same data the positive-component update converged 18/18, with the level
      0.91–1.07 × the truth and |z| ≤ 2.33.
- **Ellmer's white start is regular now:** S0 = ones (white noise) is the default.
  - Option `'C0'` is replaced by `'S0'` (the starting spectrum), with error
    `lsqvce:adjustVCEjoint:badS0`.
  - `'Norm'` is retired from `adjustVCEjoint` and `arcVarianceFactors`.
  - The warning `lsqvce:psdVCE:degenerateStart` is retired: its fixed point belonged to the
    DCT weighting.
- **`psdVCE` is rewritten** as `psdVCE(Q, T, Sprior)` on per-bin statistics.
  - `'Gamma'`, `'Smooth'` and `'FreezeRedundancy'` stay. A bin now freezes on its
    redundancy, `Sprior_j · sum T`.
  - `'Norm'`, `'MinAmplitude'`, `'ZeroFreqRule'`, `'AllowDegenerate'`, `.Cnew` and
    `.degenerate` are gone.
- **The arc factor is the Helmert factor of the arc's whole covariance** (r'Ω / r's over
  all lags). It used to be the DCT-summed, C-weighted lag vector truncated at Grid.
- **Each arc is weighted by its own law.** The bin statistics of an arc are scaled by the
  arc's own variance r(0), read from the covariance it was weighted with. The one pooled C
  that weighted every arc (R6 WP-100) is gone, because K does not depend on the law.
- **The law is reported as the bin spectrum.**
  - `out.S` / `out.Sgroups` and `R.noise.S` / `R.noise.Sgroups`, next to its covariance
    `C` / `Cgroups` at lags 0..Grid-1, which is derived and kept for the record.
  - `out.Shat`, `R.noise.Shat` and `history.Shat` are renamed `S`.
  - `R.plot("psd")` draws the bins.
- **`jointModelOpt` takes `model.S` / `model.Sgroups`** and emits the 'psd' spec with
  ordinates 2·S on fNorm, scaled to r(0) = 1, so covModel reproduces K·S exactly.
  - Error `badC` is now `badS`.
  - The cosine transform of a truncated ACF, floored at eps·max, is gone.
- `Noise.Grid` is the spectral resolution. It no longer has to cover the longest arc:
  the law is exact at every lag whatever the Grid.
- `redundancyLags`: the general (non-Toeplitz) Σ⁻¹ lag sums run in BLAS-3 column blocks
  (memory N × 256), because the joint sweep now asks for every lag.
- **Every estimated law changes.** Runs with `Estimate = "noise"` give different numbers
  from 6.6.2.

## 6.6.2 — 2026-10-11 (WP-111 part 2: the PSD law's covariance without a lag period)

- **A declared PSD law no longer becomes a circulant covariance** (review R7, R7-01).
  `lsqvce.kernel.psd2acv` integrates the piecewise-linear spectrum through the given
  ordinates EXACTLY (a closed form per segment) instead of by trapezoids, and `covModel`
  'psd' 'toeplitz' no longer resamples the law onto floor(m/2)+1 points. Trapezoidal
  quadrature on nf points is periodic in the lag with period 2(nf-1): with the joint
  estimator's 60-bin law, epoch 1 of every arc longer than 118 epochs was correlated with
  epoch m like a neighbour (r(m-1) = r(1) on even arcs, = r(0) on odd ones, Toeplitz cond
  3.3e8 at φ 0.5 -- measured on the runner, HANDOVER §12.44). Now r(m-1)/r(0) ~ 1e-8 and
  cond 9 on every arc length; r(0) = trapz(f, S) as before. Numbers change for every
  'psd' law whose grid is coarse against the arc (the joint estimator's laws above all);
  fine grids (thousands of points) move by the quadrature error only.

## 6.6.1 — 2026-10-11 (WP-115 part 4: one stream key, one arc key)

- **The stream key and the arc key each have one definition** (review R7, R7-43; coding
  rule R1). `lsqvce.util.streamKey(sysID, obsType)` returns `'satKey|obsType'` (and the
  satKey); `lsqvce.util.arcKey(sysID, obsType, day, arc)` returns
  `'s<sys>_<obsType>_d<day>_a<arc>'`. They were composed inline in eight places and the
  stream was recovered from the arc key by a regex in three. The joint estimator now
  carries each arc's stream next to it: new output `out.arcStreams`, aligned to
  `out.arcKeys`. Every key keeps its format; no number changes.
- `lsqvce.engine.satKeyOf` is retired: nothing in the toolbox called it, and its `'s_'`
  prefix was a third convention for the same name.

## 6.6.0 — 2026-10-11 (WP-115 part 3b: the robust pass as an additive nugget; an unbiased robust level)

- **The joint estimator's robust factors are an additive NUGGET** (review R7, R7-10; user
  decision (b)). A factor w of an epoch raises its temporal variance to Σ_ee / w and leaves
  its covariances with the neighbours to the law (`covModel`'s new `cm.nugget`,
  q = 1/w − 1, composed before `Cep`). Until 6.5.0 the factor was folded into `Cep`
  (sigma / sqrt(w)), which scaled those covariances too: under a correlated law a
  zero-weight epoch stayed a phantom observation, and at φ = 0.9 the formal standard
  deviation was 0.31 of the true one (measured on the runner, HANDOVER §12.43). With the
  nugget, formal/true is 1.000 and a zero weight is a deletion. A nugget arc is whitened in
  ONE batch (the law's filter plus an exact Woodbury correction over the down-weighted
  epochs, 0.1–4 ms per 4096-epoch arc); `kernel.redundancyLags` gains `'Nugget'` /
  `'LawW'` to keep its structured Σ⁻¹ sums on such an arc.
- **The level of a robust joint run is no longer biased low** (R7-12). Measured over 20
  seeds with the true law declared: plain 0.997, IGG3 0.81 / 0.86. The influence selects
  the epochs with large residuals; counting the redundancy with the weights does not undo
  that (a mirror moved 0.74 → 0.76). Now (user decision) the component update sees only
  the REJECTED epochs (deleted), partial factors weight the solution only, and each
  stream's level is divided by the Gaussian truncation factor of the rejection point
  (0.9733 for IGG3's k1 = 3; 1 for Huber), returned as `R.noise.robust.truncation`.
  The returned solution is hard as well (rejected epochs deleted, partial factors 1), its
  `sigma0` / `Cxx` corrected for the rejected rows and τ: with the partial factors in the
  final solve, σ·σ0 stayed 0.83 on the runner even with the corrected levels.
  Every robust `Estimate = "noise"` fit changes its numbers.
- `kernel.robustWeights` returns the influence's rejection point as a fifth output.
- New identifiers: `lsqvce:covModel:badNugget`, `lsqvce:covModel:nuggetNeedsWholeArc`,
  `lsqvce:redundancyLags:badNugget`.

## 6.5.0 — 2026-10-11 (WP-86 part 4: one verdict unit; the verdict's methods on R)

- **`lsqvce.engine.verdict` is the ONE engine unit that judges a fit** (the R5 M4
  "verdict from the sweep" row; user decision: the full absorption). `diagnoseAccum`,
  `diagnose`, `snoop`, `systematic`, `ftest` and `ellipse` are its static methods —
  `verdict.source`, `verdict.diagnose`, `verdict.snoop`, `verdict.systematic`,
  `verdict.ftest`, `verdict.ellipse` — with their bodies unchanged; the new
  `verdict.fromAccumulator(acc, sol, ...)` is the source + diagnose pair every driver
  ran (fit, adjustBlock, adjustRobust, adjustNonlinear). The distributions moved to the
  KERNEL tier: **`lsqvce.kernel.stat`** (was `lsqvce.engine.stat`). Breaking for code
  that called the six engine units or `engine.stat` directly: prefix `verdict.`
  (`lsqvce.engine.diagnose(...)` → `lsqvce.engine.verdict.diagnose(...)`,
  `diagnoseAccum` → `verdict.source`) and `engine.stat` → `kernel.stat`. Error
  identifiers are unchanged.
- **Everything the verdict tells you to do is a method on `R`.** New `Results` methods:
  `R.systematic(...)` (Ljung–Box / runs / binned means on the residuals),
  `[D, REJ, sol] = R.snoop(...)` (iterative snooping; R unchanged), `T = R.reliability()`
  (per observation: v, r, w, MDB and the (w, r) quadrant; `lsqvce:Results:noPerObservation`
  without the per-observation diagnostics), `E = R.ellipse(i, j)`. The report's remedies
  and snoop's runtime messages name `R.systematic()` / `R.snoop()`, never an engine
  function (finding H-05).
- `R.diag.store` records what the diagnostics kept per observation.

## 6.4.0 — 2026-10-11 (WP-86 part 3: one block formula, one commons factor)

- **The VCE redundancy pass runs on `lsqvce.engine.BlockCofactor`.** `engine.redundancy`
  (behind `engine.vce`, `outliers`, `contribution`) built, per system, its own copy of the
  commons factor and the nc × nl products (`leverageBlocks`) -- nSys copies of an
  nc × nc matrix, about 0.7 GB each at a GRAVFIRE month -- and evaluated the same
  block-arrow leverage formula a second time (`leverage`). It now builds ONE cofactor
  from the accumulator (`BlockCofactor.fromAccumulator`, the same one the diagnostics use)
  and broadcasts it to every arc task. **`lsqvce.engine.leverage` and
  `lsqvce.engine.leverageBlocks` are retired** with their five identifiers
  (`lsqvce:leverage:*`, `lsqvce:leverageBlocks:*`); `BlockCofactor.leverage(s, Ac, B)`
  is the per-system form (`h` of system s's rows from their commons and local parts).
- **`snoop`'s clean re-solve keeps the block structure.** With `Reaccumulate=true` the
  surviving rows used to be re-solved with every parameter as a common one (a dense
  u_full factor). Each original system is now rebuilt from its surviving batches with its
  own locals (`engine.batchSystem` gained `'LocalColumns'` / `'SystemId'`), the batches
  that touch no locals as one local-free system: the original structure and column
  layout, and the returned `sol.Q` is block-arrow again.
- `BlockCofactor` gained `systemOf` and the static `fromAccumulator`; a local-free system
  keeps its index in `Q.sys`.

## 6.3.0 — 2026-10-10 (WP-86 part 2: the diagnostics' cofactor is kept as blocks)

- **No dense factor of the full parameter set (R7-24).** The diagnostics used a dense
  u × u factor over commons plus all locals (13 729 at a GRAVFIRE month, 1.5 GB per copy;
  `diagnose` held four: R, its inverse, Qxx and the correlation matrix). The new engine unit
  `lsqvce.engine.BlockCofactor` answers every question they ask of Qxx — leverage, Qxx·a,
  diag(Qxx), a 2 × 2 block, the strongest correlation, a conditioning bound, snoop's
  one-row removal — from the reduced commons factor and each system's local blocks, the
  ones the solve used: about 1 GB at that scale. Verified against the dense inverse to
  ~1e-16 (the Python mirror wp86/blockcofactor.py in the R7 review folder) and pinned by
  seven new tests.
- **`diagnoseAccum` reads no data itself.** One batch is one arc; the index comes from
  the loaders' memoised meta answers and the row count from the accumulator. A VCE run with
  `CacheDays = 1` builds each day 2n+3 times (2n+4 in 6.2.0).
- **`snoop` downdates the blocks.** The rank-one removal is a rank-one update of the
  commons factor plus the observation's own local block, and the estimate follows in
  closed form (x + p·v·Qxx~·a): no row-search pass, no right-hand-side re-accumulation pass
  per rejection. The clean re-solve reads one batch at a time (it concatenated every
  surviving row at full width), but still treats all parameters as commons.
- **Changed result fields (breaking for code that read them):** `D.R`, `D.Qxx` and
  `D.corr` are gone, and `R.source.R` with them; the cofactor is `D.Q` / `R.source.Q`
  (`D.Q.block(idx)` is `Qxx(idx, idx)`, `D.Q.diagonal()` its diagonal). `D.condN` /
  `D.condR` are now a LOWER bound (`D.condIsBound`), formed without an SVD; the report says
  so. Retired identifiers: `lsqvce:diagnose:largeU` (the dense matrices it warned about are
  gone), `lsqvce:snoop:badBatch`, `lsqvce:snoop:badIndex` (their helpers are gone). New:
  `lsqvce:BlockCofactor:badBlocks`, `lsqvce:BlockCofactor:multiSystemRow`.
- **What it does not change:** the per-observation cost stays O(nc²) — the commons
  dominate (about 1.7× fewer flops than the dense factor at GRAVFIRE scale).

## 6.2.0 — 2026-10-10 (WP-86 part 1: the diagnostics no longer store the data)

- **Every fit's diagnostic source streams (R7-24).** `lsqvce.engine.diagnoseAccum` stored
  every whitened row at full width — n × (commons + all locals) doubles, kept alive in
  `R.batchFcn` — about 228 GB for one GRAVFIRE-scale month. It now folds each batch into
  the full factor during one enumeration pass, keeps only the batch's address (system,
  day, obsType, arc, sub-batch), and rebuilds a batch on demand through the system's
  loader and the same whitening iterator, caching one arc. The arithmetic is unchanged
  (same batches, same order); the price is one re-read of the data per diagnostic pass,
  and the rebuild relies on the loader contract that an arc is returned identically on
  every call.
- **What remains (WP-86 part 2):** the diagnostics still build a dense factor over the
  FULL parameter set (≈ 1.5 GB at 13 729 parameters) and cost ~u²/2 flops per row and
  pass (~2·10¹⁴ for a GRAVFIRE month). The README says so.

## 6.1.9 — 2026-10-10 (R7 WP-118 part 3: the memory claim corrected; long code on the record)

No behaviour change.

- **The README no longer says "nothing of size n is ever formed" (R7-24).** The solve
  streams; the diagnostics stage does not yet: every fit rebuilds the whitened rows at full
  width, n × (commons + all locals), and keeps them in the `Results` — about 228 GB for one
  GRAVFIRE-scale month, which runs out of memory after the solve has succeeded. Streaming
  that stage is WP-86 (open). `diagnose`'s help and `dev/DEVELOPMENT.md` say the same.
- **Long functions and long lines are on the record (R7-44).** The five functions over 400
  lines (`diagnose`, `report`, `adjustVCEjoint`, `snoop`, `makeLoader`) carry a written
  reason and a code-line ceiling in `dev/docs_src/long_functions.csv`, and the 74 lines over 120
  characters a per-file ceiling in `long_lines.csv`; audit check 33 refuses growth, a new
  unrecorded case, and a stale row. `makeLoader`'s 569-character struct literal is wrapped.

## 6.1.8 — 2026-10-10 (R7 WP-118 part 2: four reported numbers corrected; MeanType retired)

- **`R.compare` counts every parameter (R7-19).** The nested-model F test took the
  redundancy from the common parameters alone; on a streaming fit the eliminated locals are
  parameters too, so the full model's s0² was too small and F too large: 7.4 % false
  positives at α = 0.05. It now uses `R.diag.u` (commons + locals). Kinds B/C (no locals)
  are unchanged.
- **The whiteness test is sized on short arcs (R7-38).** The per-stream Ljung–Box
  statistic of the joint estimator (`R.history.whitePval`, `white`) pooled the arcs' lag
  products but divided by (N − k) as for one long series. It now divides by the number of
  products each lag actually holds. Measured size under a white null: 64 arcs × 30 epochs
  0.001 → 0.048; 200 × 12 0.000 → 0.051; long arcs and single arcs unchanged.
- **Kind B/C warnings reach `R.meta.warnings` (R7-23).** The run's warning inventory was
  reset after the batch path had already warned (`lsqvce:fit:priorNotActive`), so the id
  was lost — and after a run that fired it, the next one warned nothing at all. `fit` now
  resets once, first.
- **Worker warnings reach the run inventory (R7-41).** With `Parallel.Mode = "parfor"` or
  `"spmd"`, the log sink never reached the workers: their warnings were missing from
  `R.meta.warnings` and a once-per-run warning printed once per worker. The queue is now
  handed to every task, the tasks return their warning ids, and the client merges them in
  task order (the order a serial run fires them in); a once-per-run warning shows once.
  New public helpers `lsqvce.util.log.attachWorker` / `mergeInventory` for callers
  running their own parfor.
- **`Noise.MeanType` is retired (user decision).** Since 6.1.6 the arc factors are
  absolute and gamma-bar no longer scales them; setting the field is refused
  (`lsqvce:Options:badGroupField`), and `adjustVCEjoint` no longer takes `'MeanType'`. The
  kernel `lsqvce.engine.arcVarianceFactors` keeps both of the thesis' means.
- **Engine contract (user decision, documented in 6.1.6):** an accFcn driving
  `adjustVCEjoint` must apply `model.arcFactor`.

## 6.1.7 — 2026-10-10 (R7 WP-118 part 1: docs truth, one default per fact, honest numbers)

- **The "95 %" error ellipse covers 95 % (R7-39).** `lsqvce.engine.ellipse` with the
  estimated s0² (the default, `'aposteriori'`) now scales by 2·F(2, r) instead of χ²(2),
  which covered 86 % at redundancy 5 and 90.5 % at 10. With σ0 known (`'apriori'`) χ²(2)
  stays; the two meet at large r. The used constant is returned as `E.k2`.
- **Unseeded joint runs are replayable (R7-37).** With `Seed = []` the joint estimator
  pins one probe seed up front (from the private probe stream, the global RNG untouched),
  probes with it on every pass, returns it as `out.seed`, and `fit` records it in
  `R.meta.seeds.jointSeed`. Before, every pass drew a fresh seed and none was recorded.
  Behaviour change: an unseeded run now uses the same probes on every pass, as a seeded
  run always did.
- **Diagnostics no longer advance your RNG (R7-36).** Past 5000 observations the MDB
  reservoir of `lsqvce.engine.diagnose` drew from the global stream; it has its own fixed
  stream now, so the reported MDB quantiles are reproducible.
- **One default per fact (R7-35).** `diagnose`'s `rMin` defaults to 0.01, the value
  `Options.Diagnostics.RMin` always had (it was 1e-6), and `snoop`'s re-diagnosis inherits
  the thresholds of the D it was given. The documented-but-dead `diagnose` option
  `'verify'` is retired (passing it is an error now). `adjustBlock`'s help no longer
  documents `ShowSeconds` / `ExportDir`, which it never accepted.
- **Every accepted option is documented, every named unit exists (R7-34, R7-35).** Two
  new audit checks: 31 (each Name-Value option an `arguments` block accepts is named in
  its help; 17 were not) and 32 (each `lsqvce.*` unit named in help text, runtime
  messages, the Manual, README, templates and examples exists; 33 named deleted units,
  among them `snoop`'s runtime remedy).
- **The verbose locals line reports the right size (R7-40).** `localCovariance` takes the
  solve's s0² (`acc.localCovariance(sysID, sol.Cxx, sol.s0sq)`, now one formula shared
  with `residuals`); the log line used a pre-solve hint and read ~300× too large. API
  change on this engine method: the third argument is required.

## 6.1.6 — 2026-10-10 (R7 WP-116 part 2: the arc factors weight the solve)

Behaviour change for every `Estimate = "noise"` run with `ArcFactors = true`; no declaration
changes.

- **Per-arc factors weight the solve (R7-13, user decision: Ellmer's convention).** Until
  now an arc's estimated factor was reported but only its stream's mean reached the
  weighting, so a bad arc kept its stream's weight. Each arc now enters the next pass with
  its own variance (stream level × arc factor), and the update (eq. 6.5.28) multiplies by
  that same factor. Measured on one arc 10× in sigma among nine, beside two other streams
  (16 seeds): the true scatter of the estimated parameters is 13–47 % smaller.
- **Arc factors are absolute — a 6.1.5 regression fixed.** The 6.1.5 per-stream form,
  (gamma2_m / gamma-bar_s²) · mean(gamma2_s), counted a stream's level twice as soon as its
  arcs differ: on the fixture above the factors came back 1176 / 11.65 against the planted
  100 / 1, and the other stream's level 30 % low. A factor is now the arc's Helmert factor
  gamma2_m itself (100.0 / 1.00). 6.1.5's own fixtures had homogeneous arcs and could not
  see it. `Noise.MeanType` (gamma-bar) now sets only the scale of the clamp floor.
- A caller driving `lsqvce.engine.adjustVCEjoint` with its own accFcn must apply
  `model.arcFactor` (`lsqvce.engine.jointModelOpt` does, as `covOverride.byArc`); one that
  drops it gets compounding factors and the `arcFactorsDiverging` warning.

## 6.1.5 — 2026-10-10 (R7 WP-116 part 1: the joint-loop policy)

Behaviour changes in `Estimate = "noise"` runs; no declaration changes.

- **One stream's level no longer depends on another's (R7-09).** The arc factors were
  normalised by gamma-bar over every stream's arcs, so an arc-factor stream's reconstructed
  sigma moved with an unrelated stream's truth (1.01 / 0.79 / 1.15). gamma-bar is now formed
  per stream, and an arc-factor stream keeps its own absolute level (the Helmert mean the
  other streams already get). Note: the review's proposed remedy, per-stream gamma-bar alone,
  was measured on its own mirror and loses that level (a stream truly 3× its declared sigma
  came back 1.93); the shipped form recovers 3.02.
- **The joint loop's stop rule covers every free component (R7-11).** It used to watch the
  law alone, so a run with no `"law"` stream stopped "converged" on pass 1 under the
  DECLARED weights and prior groups were never updated. It now also needs every level and
  arc factor to settle (`R.history.levelDev`), and the returned solution is accumulated
  under the returned model — before, x, Cxx and sigma0 were one pass behind the noise model
  handed back next to them.
- **`Noise.Smooth` and `Noise.FreezeRedundancy` do what they say (R7-22).** Both were
  documented and reached nothing. `Smooth` is `psdVCE`'s moving average; `FreezeRedundancy`
  is a real freeze: a frequency whose summed redundancy is below it holds its value. The
  GRAVFIRE recipe's `FreezeRedundancy = 100` now takes effect.
- **A `var` law under `Estimate = "noise"` is refused by name (R7-07)**
  (`lsqvce:fit:varUnderNoise`): the joint sweep has no vector-AR path and used a wrong
  covariance for it silently.
- Per-arc factors still do not weight the solve; that is the next change (R7-13, user decision).

## 6.1.4 — 2026-10-10 (R7 WP-115: robust conventions — wrong-day weights, calibrated multi-axis epochs, honest scales)

Two silently-wrong-numbers fixes, three behaviour corrections, one new diagnostic and
documentation truth. Of the design questions review R7 raised, the user decided all three on
2026-10-10; two ship here, and the third (an additive nugget instead of the multiplicative
Koch fold under correlated laws, R7-10) is a larger change to the loader, the sweep and the
lag kernel and follows in its own release.

- **Clean multi-axis epochs are no longer rejected (R7-08).** The joint loop's robust pass
  reduced a d-axis epoch to its RMS — a χ-shaped value with median near 1 — and the
  influence standardises uncentred, so `igg3` zeroed about 24 % of CLEAN 3-axis epochs (62 %
  at d = 6). The epoch statistic is now χ-calibrated (user decision): the sum of squares is
  mapped through χ²(d) to the equivalent normal deviate with the sign of the rows' sum, so a
  clean epoch is standardised exactly like a clean scalar residual (0.3 % zeroed at any d,
  measured on a mirror); blunders are still rejected. Scalar streams are unchanged.
- **A breakdown warning (user decision).** `lsqvce:robust:atBreakdown`, once per run, when a
  robust pass (stand-alone or joint) zeroes more than 45 % of its observations. The collapse
  refusal is unchanged — it remains the model-error detector 6.1.3 described. Diagnostic only;
  Huber never zeroes a factor, so it cannot fire under Huber.

- **Robust weights could be applied to the wrong day (R7-06).** `lsqvce.util.rowWeightKey`
  printed day ids with `%g`, which keeps six significant digits, so seven-digit day ids —
  YYYYMMDD dates, Julian Day Numbers — collapsed onto one key: with `Robust` on, one day's
  blunder factors were applied to another day of equal length, or a misleading
  `makeLoader:rowWeightSize` blamed the builder. Integer ids now print with `%d` (exact,
  the loader cache's own format), others with `%.17g`. MJD and datenum days (< 1e6) were
  never affected.
- **The robust scale's fallback is no longer silent (R7-21).** When half or more of the
  whitened residuals tie (quantised data, zero-redundancy rows) the MAD is degenerate and
  `lsqvce.kernel.robustWeights` falls back to the std; it now warns once per run
  (`lsqvce:robustWeights:scaleFallback`), returns the path taken as a fourth output, and
  treats a near-tie MAD (a 1e-15 spread) as a tie instead of a scale — which had given every
  ordinary row a factor of ~1e-14 that the collapse guard could not see.
- **The IRLS stop rule is unit-invariant (R7-17).** `adjustRobust` stopped when
  `max|dx| < Tol·max(max|x|, 1)`; the floor of 1 made SH-scale coefficients (~1e-9) stop at
  pass 2, "converged", about one formal sigma from the fixed point. The scale is now the
  larger of `max|x|` and the largest formal sigma: unchanged for parameters of order one and
  larger, and the same pass count in any units.
- **`R.report()` says when a fit is robust (R7-18).** `adjustRobust`'s help promised it and
  the report said nothing; a robust fit's sigma0 and covariance are those of the reweighted
  problem, not an M-estimator's scale or sandwich covariance. One level-1 NOTE
  (`ROBUST_FIT`) now says so.
- **Claims withdrawn.** The joint loop's `collapsed` refusal no longer says "beyond the
  breakdown point" (R7-25: the claim 6.1.3 retired from `adjustRobust` survived there, wrapped
  across two help lines and inside the runtime message), and counts epochs, which is what it
  compares. "The level a robust pass returns is UNBIASED" is withdrawn to an open question
  (R7-12) in the help, the Manual and the test that pinned it; it is settled by measurement
  in a later part. Audit check 28 now reads wrapped help lines and runtime messages, so
  neither route can carry a retired claim again.

## 6.1.3 — 2026-09-21 (R6 WP-110: what the collapse guard really detects, and a 6.1.2 claim corrected)

Documentation, a corrected error message and tests; no estimator behaviour changed.

- **6.1.2 said the collapse guard "can fire only for `n < 2(nc+1)` … it cannot fire at any
  contamination". That is wrong, and this release withdraws it.** The half-sample bound
  6.1.2 proved is real, but it is conditional on `|med| ≤ 3.4·MAD`, and that premise is the
  50 % breakdown point wearing a different hat: measured on deterministic samples, it holds
  at 49 % contamination (ratio 0.50) and fails at 51 % (ratio 518.89). Above breakdown
  `nAlive` falls with the clean fraction and the guard does fire. The bound therefore says
  no more than "we are below breakdown"; 6.1.2 drew a consequence its own premise forbids.
- **What the guard actually detects, measured.** `nAlive` tracks the CLEAN row count to the
  row, and the threshold is an absolute `nc + 1`, so contamination reaches it only when
  almost nothing clean is left (99.8 % on a 3600-row fit against `nc = 6`; 60 % contamination
  of 1000 rows still leaves 400 alive). What does reach it on a full-sized problem is
  residuals offset **bodily** from zero by several MAD — an unabsorbed bias, a datum defect —
  because the influence standardises `vw/s`, not `(vw − med)/s`. That is a model error, not
  an outlier problem, which is what the message always meant to say and now does.
- **`lsqvce:adjustRobust:collapsed`'s message and help are corrected** to name those two
  routes instead of claiming a 50 % breakdown detector. The phrase it used to carry is in
  `superseded_claims.csv`, so it cannot come back. **The threshold itself is unchanged** —
  turning the guard into a genuine breakdown test means choosing a policy, which is the
  user's call, not a release note's.
- **Issue #110's test 6 is closed**: the arc-scoped local fixture with priors that WP-109
  said it needed now exists, and a degenerate arc — every row rejected — is carried by its
  local prior alone, with plain GLS on the same data as the control.

## 6.1.2 — 2026-09-20 (R6 WP-109: the PDF Manual is gated, and issue #110 is closed out)

Documentation and tests; no code in `+lsqvce` changed, no declaration changes.

- **`lsqvce_manual.pdf` carried 6.1.0 into the 6.1.1 release, because nothing read it.**
  It is rebuilt, and audit check 29 now holds its document title to `+lsqvce/Contents.m`
  — the same single version source check 17 (a) holds the README title to. The PDF cannot
  be gated the way the HTML Manual is (check 23 compares a fresh build byte for byte):
  rebuilding the identical source changes the font subset tags, the document `/ID` and
  the xref offsets, which is measured here rather than assumed — 21 of 141 inflated
  streams differ, every one of them accounted for, and not one content byte. The title is
  what can be gated and the title is what drifted. The check fails CLOSED: a missing file
  or an unreadable Info dictionary is a finding, not silence.
- **Issue #110's remaining tests (3–6) are in the suite**, each with its margin measured
  in the same run rather than assumed: planted blunders against the fixture's own
  clean-data error under `igg3` and `bisquare`; the IRLS fixed point on kind C with an
  explicit `p`; and `arma` against the PSD evaluated from it, where the AR-order gap is
  measured by the GLS pair so the robust pair can be asked the only question IRLS can
  answer — that re-weighting does not amplify it.
- **A finding, recorded not fixed: the collapse guard is close to unreachable.**
  `adjustRobust` refuses a pass when `nAlive < nc + 1`, but `robustWeights` standardises
  by the normalised MAD, so at least half of any sample survives a redescending influence
  at its default tuning: `nAlive >= ceil(n/2)` always, and the guard can only fire for
  `n < 2*(nc + 1)`. On a real-sized fit it cannot fire at any contamination. 6.1.1
  reported that prior pseudo-rows padded `nAlive` — a guard failing open; the guard was
  nearly unreachable to begin with. The tuning is the user's, so nothing is changed here.

## 6.1.1 — 2026-09-20 (R6 WP-108: the stand-alone IRLS stops re-weighting priors)

Fix and tests only; no declaration changes, no new capability.

- **`Robust` alone no longer lets a local prior group into the robust sweep.** 6.1 said
  "Prior streams are never re-weighted", and `adjustVCEjoint` did exactly that. The
  STAND-ALONE driver did not: a local group rides as its own stream `prior_<name>`, and
  `adjustRobust`'s residual collection took every stream. Three consequences, all
  measured on CI before the fix (issue #110): the pseudo-residuals entered the MAD scale
  that sets every DATA row's factor; they were counted in `nAlive`, `nDown` and `nZero`,
  so the collapse guard could pass on data already past the breakdown point — a guard
  failing OPEN; and because `priorArc` never applies `rowWeights`, the next pass divided
  a prior residual by a `sqrt(factor)` that had never been applied, so that residual grew
  and its bookkeeping factor fell, pass after pass. `private/lawWhitenedResiduals` now
  skips `prior_` streams, which is the whole fix — `adjustRobust`'s write-back iterates
  the keys that function returns, so nothing is written under a prior key any more.
  The two halves of one rule are now written in the same terms.
- **Not changed, deliberately:** `priorArc` still does not apply `rowWeights`. Making the
  bookkeeping consistent the other way would down-weight a prior for disagreeing with the
  data, which is what a prior exists to resist.
- **`adjustRobust`'s `% shared with adjustVCEjoint` comment was wrong** and is corrected:
  `lawWhitenedResiduals` is called only from `adjustRobust`. The joint loop has its own
  producer and its own exclusion.

## 6.1.0 — 2026-09-09, tagged 2026-09-19 (R6 WP-105: the robust pass inside the joint loop)

Tagged as `v6.1.0` on 2026-09-19 (WP-107): the version string gains its PATCH field
(ci-actions `VERSIONING.md`: `MAJOR.MINOR.PATCH`, one annotated tag per release). The
content is the 6.1 of 2026-09-09; no code changed.

Additive on 6.0; no declaration changes.

- **`Robust` together with `Estimate="noise"` runs — in the Koch convention.** Every pass of
  the joint estimator ends with one residual sweep; the RAW residuals, descaled by each
  epoch's declared accuracy and reduced to one value per epoch, go through the influence
  function (`huber` | `igg3` | `bisquare`, `RobustOpts.Tuning`) per stream, and the factors
  are folded into the per-epoch error information of the next pass (sigma / sqrt(w),
  BEFORE the temporal law) — the one composition C = D (R_t ⊗ I_d) D carries them
  everywhere. Measured: weighting the law-whitened innovation instead (Yang, the
  stand-alone IRLS's convention) neutralises the filter tail around every blunder and a
  spectral estimator cannot recover the colour (0.054 against 0.20 in 16 passes); with the
  same blunders declared through `Cep` the plain estimator recovers it.
  `MaxPasses` is the budget; `RobustOpts.MaxIter` / `Tol` belong to the stand-alone IRLS
  only. Prior streams are never re-weighted.
- What comes back: `R.noise.robust` (`influence`, `scale` and `streams` — one MAD scale per
  stream, `factor` per observation in sweep order, `keys`), `R.history.robustScale` (one
  column per stream) / `nDown` / `nZero` per pass, and the same fields on the `OnPass`
  context. The residuals are standardised by the stream's own MAD, never by the estimated
  level (a global scale, or the level itself, sends a clean stream's level into a spiral —
  measured). Prefer `igg3` / `bisquare` for gross blunders: Huber keeps c/|u| of their
  weight.
- The level is UNBIASED under the robust pass (measured 1.007 on a clean stream under IGG3;
  the retained-observation ratio E[w u²]/E[w] that was expected was refuted by the oracle):
  Ellmer's per-lag redundancy is the model expectation under the same weighted whitening,
  so the weights cancel in alpha². Multi-axis streams are weighted per EPOCH (the RMS over
  the d law-whitened rows, one factor per epoch). This is Koch's iterated re-weighting
  inside a VCE, not Ellmer (2018), whose robustness is preprocessing screening plus the
  arc-wise variance factors (Manual chapter 7, §7.8).
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
