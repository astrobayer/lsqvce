# Reading the results and the report

`lsqvce.fit` returns a `lsqvce.Results` object. This chapter covers how to read it:
first the object itself, then the five-level report that judges it, then the per-stream
evidence of an `Estimate = "noise"` run, which the report does not print. For each
adverse finding it says what to do. The full field lists are in the API reference, and
the derivations of every statistic are in the theory.

```matlab
R = lsqvce.fit(data, opts);
R.report                      % the verdict, levels 0..4
V = R.report(Verbose=false);  % the same verdict, silent, as a struct
```

## The `Results` object

### Properties

All properties are read-only.

| property | what it holds |
|---|---|
| `R.x` | the common-parameter estimate, `nc x 1` (the eliminated local parameters are not in it; `R.toTable` and `R.residuals` give them back) |
| `R.sigma0` | the a-posteriori standard deviation of unit weight |
| `R.sol` | the solver's raw output: `.x`, `.sigma0`, `.Cxx` (the covariance, **already scaled by `sigma0^2`**), `.vtv`, ... |
| `R.diag` | the per-observation diagnostics the report reads (next sections) |
| `R.history` | the pass record of an iterative run; empty for a single solve |
| `R.noise` | the estimated stochastic model (`Estimate = "noise"`), or the robust record |
| `R.converged`, `R.stopReason` | whether the loop settled, and why it stopped |
| `R.meta` | run provenance: `.toolboxVersion`, `.matlabRelease`, `.createdUTC`, `.seeds`, `.options` (as resolved), `.warnings` |
| `R.acc`, `R.source`, `R.batchFcn`, `R.nBatch`, `R.systems` | the accumulator, the diagnosed full-parameter solution and its replayable batch source (for expert use) |

`R.stopReason` takes one of five values:

| value | meaning |
|---|---|
| `"single"` | one solve, nothing iterated (`Estimate = "none"` without `Robust` or `x0`) |
| `"converged"` | every stop criterion was met |
| `"maxPasses"` | the pass budget ran out. The last pass is returned with its covariance and one warning: this is a result, not an error |
| `"userStop"` | your `Options.OnPass` hook returned `false` |
| `"diverged"` | Gauss–Newton found no downhill step; more passes cannot help |

### Methods

| method | answers |
|---|---|
| `R.report(...)` | Is the adjustment acceptable? Prints the five levels and returns the verdict |
| `R.coefficients(Alpha=0.05)` | the estimates as a table: `index`, `estimate`, `se`, `lower`, `upper` (normal bounds, `se = sqrt(diag(Cxx))`) |
| `R.toTable(data, Estimate=9, Std=10)` | your parameter tables with the estimate and its standard deviation written into the named columns, one row per instantiated parameter (commons and locals) |
| `R.reliability()` | one row per observation: `obs`, `v`, `r`, `w`, `mdb`, `quadrant` |
| `R.snoop(...)` | iterative data snooping on this fit; `R` itself is not changed |
| `R.systematic(...)` | runs, Ljung–Box and binned-mean tests for unmodelled signal in the residuals |
| `R.compare(Rred)` | nested-model F test: does the full model `R` earn its extra parameters over `Rred`? Both fits must share one fixed weighting. After `Estimate = "noise"` each model's v'Pv is normalised to its own redundancy, so F would be 1 whatever the data, and such fits are refused. Declare the noise model the full fit estimated, with `estimate = "none"`, refit both, and compare those. |
| `R.ellipse(i, j, ...)` | confidence ellipse of the parameter pair `(i, j)` |
| `R.residuals(...)` | the collected residuals per system and stream (re-reads the arcs) |
| `R.contribution(...)` | the information share of every stream in the common parameters |
| `R.summary()` | one row per **stream**: what each ended up as (iterative runs) |
| `R.convergence(...)` | one row per **pass**: when the loop settled (iterative runs) |
| `R.plot(kind, ...)` | `"residuals"`, `"contribution"`, `"psd"`, `"convergence"`. Returns one figure handle and exports nothing |
| `R.redundancy()` | `sum(r_i) = n − u` |
| `disp(R)` | a one-screen summary |

Some methods cannot run on every result:

* `summary`, `convergence`, `plot("psd")` and `plot("convergence")` need a variance-component history.
  On a single solve they raise `lsqvce:Results:noHistory`.
* `residuals`, `contribution`, `toTable`, `plot("residuals")` and `plot("contribution")` need the
  streaming accumulator. On the two in-memory kinds they raise
  `lsqvce:Results:noAccumulator`; read `R.diag.v` there instead.
* `reliability` and `systematic` need every residual (`Diagnostics.Store = "all"`, the
  default). Under `"flagged"` or `"none"` they refuse (`lsqvce:Results:noPerObservation`,
  `lsqvce:systematic:noResiduals`). The report still runs, because it reads only aggregates
  and the flagged list.

**Standard deviations.** `R.sol.Cxx`, `R.coefficients` and `R.toTable` all carry the
a-posteriori scaling. An arc's standard deviation in an `Estimate = "noise"` run has three
factors: declared `sigma` × `sqrt(arc factor)` × `R.sigma0`. The chapter *Choosing the stochastic model* has the details.

## How the report is built

`R.report` is `lsqvce.report(R.diag)`. It also passes along the run provenance and the
loop state, so the printed report gets one provenance line at Level 0 and one
loop-settled line at Level 1 (described below). A bare `lsqvce.report(R.diag)` prints
neither line.

| name-value | default | effect |
|---|---|---|
| `Verbose` | `true` | print the report; `false` only returns `V` |
| `Width` | `78` | line width, ≥ 40 |

The report always evaluates and prints all five levels. There is no option to print only
one level; to filter, use the `level` field of the findings. The text goes through the
toolbox's message log, like every other status line.

The returned struct:

```matlab
V.verdict       % 'PASS' | 'PASS WITH WARNINGS' | 'SUSPECT' | 'INVALID'
V.message       % the sentence printed under the verdict
V.findings      % cell of struct(level, severity, code, msg)
V.nError, V.nWarn, V.nNote, V.integrityOK
```

Each finding has a severity, `ERROR`, `WARN` or `NOTE`, and a stable **code**. A pipeline
should branch on the codes, not on the wording:

```matlab
codes = cellfun(@(f) f.code, V.findings, 'UniformOutput', false);
if any(ismember(codes, {'SNOOP_CONTAMINATED', 'SNOOP_REFUSE', 'SNOOP_BREAKDOWN'}))
    R2 = lsqvce.fit(data, opts, Robust="igg3");
end
```

### The tests behind it: `Options.Diagnostics`

The diagnostics are computed by `fit` itself, in a second sweep over the data. Their
settings belong to the `Diagnostics` option group:

| field | default | meaning |
|---|---|---|
| `Alpha` | `1e-3` | significance of each per-observation test |
| `Bonferroni` | `true` | divide `Alpha` by `n`: the family-wise rate over all `n` tests |
| `Beta` | `0.2` | type-II error for the MDB (power 0.8) |
| `RMin` | `0.01` | redundancy number at or below which an observation is **uncontrolled**: its `w` is set to 0 and its MDB to `Inf` |
| `RWeak` | `0.1` | redundancy number below which an observation is **weakly controlled** |
| `Test` | `"auto"` | `"pope"` (studentised by the estimated `s0`), `"baarda"` (by the a-priori sigma0), `"mad"` (by a robust MAD scale); `"auto"` selects `"pope"` |
| `RobustCheck` | `true` | also form the MAD scale and the contamination gauge |
| `Store` | `"all"` | keep every per-observation vector, or only the `"flagged"` ones, or `"none"` |
| `NWorst` | `1000` | how many largest-`|w|` observations are retained (`R.diag.worstW`, `.worstIdx`) |

```matlab
opts = lsqvce.Options(Diagnostics=struct('Test', "mad"));   % fields merge into the defaults
```

**What the numbers refer to.**

* **Observations.** An observation index is a row in **sweep order**: the arcs in the
  order they were accumulated, then the rows of a declared prior. Prior rows are counted in
  `n` and in the redundancy.
* **Units.** On the streaming kind every row is whitened by its declared law before it is
  tested. `v`, `w`, `r_i` and the MDB therefore belong to whitened rows, and an MDB is in
  units of that row's own sigma. For correlated data only sums such as `sum(r_i)` and `vPv`
  are the same in both frames.
* **Parameters.** A parameter index in the report refers to the **full** parameter vector:
  the commons first (indices `1..numel(R.x)`), then each system's local parameters in
  block order. The `u` in the header counts both.
* **The a-priori sigma0.** For whitened rows the a-priori sigma0 is 1. The report's `s0`
  is the a-posteriori value recomputed from the residual sweep and should match `R.sigma0`;
  Level 0 checks this.

## The verdict

The verdict is printed last, between two rules:

```
  VERDICT: <verdict>        (<e> error, <w> warning, <n> note)
  <message>
```

| verdict | rule | what it means | what to do |
|---|---|---|---|
| `INVALID` | a Level-0 integrity check failed | The arithmetic or the setup is broken (for example the second pass did not replay the batches the solve saw). Levels 1–4 are **not evaluated**. | Fix the setup. Do not read or publish any number from this run. |
| `SUSPECT` | integrity holds, at least one `ERROR` | The adjustment computed, but at least one finding is serious enough that the result should not be published as it stands. | Work through the `[ERROR]` items, top level first. |
| `PASS WITH WARNINGS` | no `ERROR`, at least one `WARN` | Nothing disqualifies the result, but a reviewer will ask about the warnings. | Address each `[WARN ]` item or be ready to defend it. |
| `PASS` | no `ERROR`, no `WARN` | Every check passed. | Still read the `REL_EXTERNAL` line at Level 3. It is a `NOTE` whenever the reliability is normal for a network of this size, even when its absolute value is above one sigma. |

The levels are ordered so that a failure at one level makes the levels below it
unreliable. Read the report from the top, and fix the highest failing level first.

## Level 0 — INTEGRITY (does the arithmetic hold?)

These checks test the computation, not the data.

| code | severity | condition | action |
|---|---|---|---|
| `PROVENANCE` | NOTE | always (through `R.report`) | Prints the toolbox version, the MATLAB release and the creation time. |
| `SUM_R` | NOTE / **ERROR** | `sum(r_i)` against `n − u`; ERROR if the relative error > 1e-6 | The second pass does not replay the design or the weights the solve used. Make `buildDay` / `batchFcn` deterministic: the same rows in the same order on every call. |
| `VPV_SIGN` | **ERROR** | `vPv` negative or non-finite | A numerical failure. Check the declared covariances for non-SPD blocks. |
| `VPV_MISMATCH` | WARN / **ERROR** | the accumulator's `vPv` and the recomputed `sum p v²` differ by > 1e-6 (WARN) or > 1e-2 (ERROR, invalid) | The recomputed value is the one the report uses. A large gap means cancellation in the normal-equation path; re-run with `Backend = "qr"`. |
| `NO_REDUNDANCY` | **ERROR** | `r = n − u = 0` | Nothing is testable. Add observations or remove parameters. |
| `S0_FALLBACK` | WARN | the estimated variance factor was unusable | The `w` values were studentised with the a-priori sigma0. Read them as Baarda-type, not Pope-type. |

## Level 1 — GLOBAL (is the model consistent?)

Level 1 checks whether the functional and stochastic model are consistent overall.

**Loop state and robust fits.** If the run iterated (any `stopReason` other than
`"single"`), one of these lines comes first:

| code | severity | meaning | action |
|---|---|---|---|
| `NOISE_SETTLED` | NOTE | the loop met its stop rule; the line names the loop (the stochastic model, the robust re-weighting or the Gauss–Newton iteration) and its number of passes | — |
| `NOISE_UNSETTLED` | WARN | the loop stopped on `"maxPasses"` or `"userStop"` before it settled, or a Gauss-Newton fit stopped on `"diverged"` (no downhill step). Every number below rests on weights or a point that were still moving | Treat the result as provisional. The noise loop and Gauss–Newton: raise `MaxPasses` or loosen `Tolerance`; the robust re-weighting: raise `RobustOpts.MaxIter` or loosen `RobustOpts.Tol`. For the noise loop `R.convergence()` shows what was still changing (see "When the loop did not settle" below); a robust or Gauss–Newton run lists its passes in `R.history`. For `"diverged"` more passes cannot help: check the start point and the Jacobian. |
| `ROBUST_FIT` | NOTE | the run used a robust influence (`Robust = "huber"`, `"igg3"` or `"bisquare"`) | `sigma0` and the covariance belong to the **reweighted** problem. They are not an M-estimator's scale or sandwich covariance, so the global test and the parameter bounds are approximate in both directions. |

**The global test.** Two lines always appear:

* `GLOBAL_STAT` prints `s0^2/sigma0^2` and `T = vPv/sigma0^2`, together with the
  chi-square acceptance band for `r` degrees of freedom. The band is fixed at a two-sided
  **α = 1 %** (quantiles 0.005 and 0.995, divided by `r`); it does not use
  `Diagnostics.Alpha`.
* `GLOBAL_CI` gives `sigma0_hat` with its **95 % confidence interval** and the smallest
  relative deviation this redundancy can tell apart from sampling noise. **Read the
  interval, not the point estimate.**

What follows depends on whether the ratio falls inside the band:

| code | severity | condition | meaning and action |
|---|---|---|---|
| `GLOBAL_PASS` | NOTE | ratio inside the band | The model is consistent. |
| `GLOBAL_SHARP` | NOTE | ratio outside a band narrower than ±1 % | At this redundancy the test rejects any imperfection, however harmless. Judge the **ratio** and the confidence interval, not the rejection itself. |
| `GLOBAL_HIGH` | WARN | ratio above the band | The residuals are larger than the weights promised. There are three causes, and the global test cannot tell them apart: **(a) outliers**, see Level 2. **(b) Unmodelled systematics**: run `R.systematic(Covariates=...)` and plot `w` against time, elevation, orbit position and similar covariates; neither this test nor data snooping can see a smooth bias. **(c) Over-optimistic weights**: if the residuals are clean and structureless, estimate the weights with `Estimate = "noise"`, or correct the declared `sigma`. Do not rescale `Cxx` by `s0²` before you have ruled out (a) and (b). |
| `GLOBAL_LOW` | WARN | ratio below the band | The residuals are smaller than predicted. Either the weights are pessimistic (the declared `sigma` is too large; correct it or use `Estimate = "noise"`), or the model is **over-parameterised** and absorbs the noise. Test whether each parameter group is justified with `R.compare(Rred)`, and look at the correlations at Level 4. |

## Level 2 — LOCALISE (where is it broken?)

Level 2 locates individual bad observations by data snooping: one test per observation
on its normalised residual

  `w_i = (P v)_i / (scale · sqrt((P Qvv P)_ii))`.

The `scale` is the estimated `s0` (Pope), the a-priori sigma0 (Baarda) or the robust MAD
scale. With `Bonferroni = true` the level of each test is `Alpha / n`. The critical value
follows the scale: Pope's tau, the normal distribution, and for the MAD a Student t with
r/2 degrees of freedom (calibrated to the nominal level; Pope's tau there gave up to 60 times
the nominal false alarms on small networks). The price of the valid level is power on short
networks under heavy contamination: with 50 observations and 10 % blunders of 25 sigma, four
or more of five are found in about 60 % of cases (snooping with `Test="mad"` does not change
that); from about 200 observations all of them.

### The contamination cross-check

Pope's statistic divides by `s0`, and `s0` is estimated from the same data that contain
the blunders. A group of blunders therefore inflates the scale it is tested against, and
the blunders mask each other. To catch this, the report compares `s0` with a robust MAD
scale of the raw standardised residuals and estimates the contaminated fraction `fhat`:
the share of residuals beyond 2.5 robust sigmas. On clean data `fhat` sits near
`P(|Z| > 2.5) = 0.0124`. As a backstop it also prints `s0/sigma0`, which contamination
cannot suppress but which is only as good as your declared sigma0.

The null band is measured, for every `n` (the scored residuals): on clean data the count
`n·fhat` is negative binomial with mean `n · 0.0124 · (1 + 27/n)` (a small sample's MAD is
biased low) and variance 1.69 times that mean. `D.fhatTol` is its 99.75 % point (false
alarms ≤ 0.25 %), and `fhat` counts as *collapsed* only where it is implausibly low for
clean data (`P(K ≤ n·fhat) ≤ 1 %`). Below `n ≈ 450` even `fhat = 0` is common on clean
data (59 % at `n = 30`), so there a captured robust scale cannot be told from a sigma0
declared too small (`D.collapseDetectable = false`).

The cross-check appears whenever `RobustCheck = true`, whatever the `Bonferroni` setting. Exactly
one of these lines is printed:

| code | severity | condition | action |
|---|---|---|---|
| `SNOOP_SCALECHECK` | NOTE | `clean`: `fhat` within its null band `D.fhatTol` | `s0` and the robust scale agree, so the flag list can be believed. (This line is omitted under `Test = "mad"`.) |
| `SNOOP_INDETERMINATE` | WARN | `fhat` within its null band, `s0/sigma0 > max(3, 2·(1 + 3/sqrt(2r)))`, and a collapse undetectable at this `n` | Either sigma0 is declared too small (the common case: then the flag list can be believed) or more than half the data are bad and the robust scale is captured (then it cannot). Check sigma0 first; if it is trustworthy, fix the functional model before rejecting anything. `R.snoop` runs here with a warning (`lsqvce:snoop:indeterminate`), its rejections capped at `MaxRejectFrac`. |
| `SNOOP_MODERATE` | WARN | `fhat` above the null band, below 0.10 | Pope's test is starting to lose power. Re-diagnose with `Test = "mad"`, and reject one observation at a time (`R.snoop(Test="mad")`). |
| `SNOOP_CONTAMINATED` | **ERROR** | `0.10 ≤ fhat < 0.25` | Pope's test is **blind**: an empty flag list from a Pope run is not evidence of clean data, and iterating does not help. Use `Test = "mad"`, or refit with `Robust = "huber"` (`"igg3"` / `"bisquare"` for gross blunders). |
| `SNOOP_REFUSE` | **ERROR** | `fhat ≥ 0.25` | So large a fraction of "outliers" points to a model error, not to outliers. Rejecting them would produce a confidently wrong result, and the flag list is incomplete anyway. Run `R.systematic()` and fix the functional model before rejecting anything. |
| `SNOOP_BREAKDOWN` | **ERROR** | `fhat` collapsed below its null band while `s0/sigma0 > max(3, 2·(1 + 3/sqrt(2r)))` | The robust scale has adopted the blunders as its reference (more than about half the data are bad). **Reject nothing** and trust no flag list in this report. Either the declared sigma or the functional model is badly wrong; fix that first. |
| `SNOOP_NOCHECK` | NOTE | no robust scale could be formed (too few usable residuals, or `RobustCheck = false`) | Pope's test is used unchecked. Above about 5 % contamination it loses power without warning. |

In addition:

* `SNOOP_BACKSTOP` (NOTE) appears whenever `s0/sigma0` exceeds its clean-data tolerance
  `1 + 3/sqrt(2r)`. It is not printed in the breakdown and indeterminate cases, which
  already report the backstop. If your sigma0 is trustworthy, the residuals really are that much larger than
  promised.
* `SNOOP_MAD` (NOTE) states that the robust scale replaced `s0` (`Test = "mad"`).

### The snooping result

* `SNOOP_SETUP` (NOTE) prints the test, `Alpha`, the corrected level, `n`, the critical
  value `|w| > k` and the number of false alarms expected by chance (`Alpha` under
  Bonferroni).
* `SNOOP_NOBONF` (WARN) replaces it when `Bonferroni = false`. It prints the uncorrected
  level and roughly how many good observations will be flagged by chance. Set
  `Bonferroni = true`.
* `SNOOP_VACUOUS` (**ERROR**, Pope only) appears when the critical value comes so close to
  the upper bound `sqrt(r)` of Pope's statistic that a residual must carry at least 90 % of
  `vPv` to be flagged (`k^2 / r >= 0.9`, `D.vacuityShare`). Only a blunder that dominates
  the whole adjustment can then be found, and an empty list says nothing. Accept a weaker
  `Alpha` or a weaker family-wise correction, or accept that this network cannot detect its
  own blunders.
* `SNOOP_CLEAN` (NOTE) means nothing exceeds the critical value. It prints the largest
  `|w|` and its observation index.
* `SNOOP_FLAGS` (WARN) gives the number flagged, the number expected by chance, and so
  about how many are real. `SNOOP_TABLE` (NOTE) lists up to 15 of them, worst first, with
  columns `obs`, `w`, `r_i`, `v`, `MDB`. Read `w` and `r_i` **together**.
* `SNOOP_LOWR` (WARN): some flagged observations have `r_i < RWeak`. A large residual on
  a weakly controlled observation may be error absorbed from a neighbour. Rejecting it is
  close to a coin flip; prefer to re-measure.
* `SNOOP_ITER` (NOTE) prescribes the remedy. Reject **one** observation at a time (the
  largest `|w|`), re-solve and test again, because blunders mask each other and a single
  pass misses some of them:

```matlab
[D2, REJ, sol2] = R.snoop();               % Pope, Bonferroni, one at a time
[D2, REJ, sol2] = R.snoop(Test="mad");     % under contamination
REJ.idx      % rejected, in order;  REJ.refused  % tripped but r_i < RMinReject (0.10)
```

`R.snoop` takes `Test` (`"pope"` or `"mad"`), `Alpha` and `Bonferroni` name-values; left out,
they are `R.diag`'s, so `Options.Diagnostics` reaches them. It caps the number of rejections at
`max(10, 0.01·n)` and stops at a fraction of 0.25 (`MaxRejectFrac`). It refuses to run at
`SNOOP_BREAKDOWN` and `SNOOP_REFUSE` (`Force = true` proceeds with a warning), warns at
`SNOOP_INDETERMINATE` (and at `SNOOP_CONTAMINATED` under Pope's test), and tests only what
`R.diag` tests (`r_i > RMin`). It refuses
observations with `r_i < RMinReject`, re-solves cleanly
over the survivors at the end, and leaves `R` unchanged. A prior's rows, common or local,
are never rejected, and the re-solve keeps the prior as a prior (a datum stays a datum);
the re-diagnosis keeps `R.diag`'s thresholds and storage (`Store`, `NWorst`,
`RobustCheck`). The returned
`D2` counts in the original observation numbers (`D2.flagIdx`, `D2.rejectedIdx`, ...); its
per-observation vectors cover the survivors, and `D2.obsIdx` names their original numbers.
It cannot reject single rows of a
correlated batch (`lsqvce:snoop:correlatedBatch`): there, remove the epoch in `buildDay`
or use `Robust`.

## Level 3 — RELIABILITY (what can you not check?)

Level 3 shows which observations the tests could not check at all.

* `REL_HEAD` (NOTE) prints the smallest redundancy number and its observation, the mean
  redundancy `r/n`, and the MDB non-centrality `delta0`. `delta0` is computed at the
  **Bonferroni-corrected** level and `Beta`: 4.13 at a level of 1e-3 with `Beta = 0.2`,
  and larger as `n` grows. It is the normal (Baarda, known-sigma) value. Under
  `Test = "mad"` the threshold is Student t with `r/2` degrees of freedom, so the bias that
  can hide from that test is larger (about +50 % at `n = 30`, +25 % at `n = 50`, +6 % at
  `n = 200`); the line then says that the MDB, the envelope and the external reliability
  are the Baarda values.
* `REL_QUADRANTS` (NOTE) counts the observations in each cell of the (|w|, r) table.
  `|w|` is high when it exceeds the critical value; `r` is low when it is below `RWeak`.

| quadrant | \|w\| | r | meaning |
|---|---|---|---|
| clean outlier | high | high | reject with confidence |
| **AMBIGUOUS** | high | low | could be a neighbour's error; re-measure |
| verified | low | high | tested, and passed |
| **UNTESTED** | low | low | **not** "fine": the test had no power |

A residual shows only the fraction `r_i` of an error in its observation. A small `|w|` on a
low-`r` observation is therefore no evidence of health. `R.reliability()` gives the same
quadrant per observation; the rows of a datum prior carry the label `"datum"` and are not
counted in any quadrant.

| code | severity | condition | action |
|---|---|---|---|
| `REL_SILENT` | WARN | observations in the UNTESTED quadrant | These observations did not pass the test; the test could not see them. Strengthen the geometry (more or better-placed observations), or validate them externally. |
| `REL_AMBIGUOUS` | WARN | flagged observations with `r_i < RWeak` | Re-measure rather than reject. |
| `REL_ENVELOPE` | NOTE | always, when MDBs are finite | The hidden-blunder envelope: the median, 90th and 99th percentile and the worst MDB, in units of each observation's sigma. "Half of my observations could hide an undetected error larger than the median." The percentiles come from a fixed-seed reservoir of 5000 samples. |
| `REL_MDBINF` | **ERROR** | observations with infinite MDB | No error in these observations, however large, can be detected. They drop out of the percentiles above, so the envelope looks better than it is. |
| `REL_UNCONTROLLED` | **ERROR** | observations with `r_i ≤ RMin` | Their residual is forced to zero, so a blunder in them passes into the parameters unseen. The first indices are listed. Add observations that constrain the same parameters, or remove the parameter that only those observations determine. |
| `REL_WEAK` | WARN | `RMin < r_i < RWeak` | A blunder in these is only partly visible. Data snooping has little power here. |
| `REL_GOOD` | NOTE | none weak, none uncontrolled | Snooping has power everywhere. |

**External reliability (`REL_EXTERNAL`)** answers one question: how far could the worst
undetectable blunder (one MDB in the worst observation) move a parameter, in units of that
parameter's own sigma? The report names the observation and the parameter. Under a healthy
geometry this number shrinks like `1/sqrt(n)`, so it is judged against the value expected
for a network of this size, `4.0 · delta0 / sqrt(n)`:

| ratio to expected | severity | action |
|---|---|---|
| ≤ 1.8 | NOTE | Normal for this size. |
| 1.8 – 3 | WARN | The geometry is weaker than it should be. |
| > 3 | **ERROR** | One observation carries far more leverage than the rest. Find it through the minimum `r_i` and the named observation, and do not trust the formal sigmas. |

If the absolute value exceeds one sigma, the line adds that an undetectable blunder could
move the parameter by more than its stated precision. When that happens at a normal ratio,
the remedy is more observations, not a different estimator.

## Level 4 — PARAMETERS (what do you actually know?)

Level 4 checks conditioning and correlations, and lists the parameter precisions.

* `PAR_COND` (NOTE) prints `cond(R)` and `cond(N)`. On the `fit` path these are **lower
  bounds** formed from the blocks of the cofactor (no SVD), and they are marked
  `[LOWER BOUND]`. The true conditioning can be worse than printed, never better, so the
  thresholds below are optimistic.
* `PAR_ILLCOND`: WARN when `cond(N) > 1e10`, with roughly `log10(cond(N))` of 16 digits
  lost. **ERROR** when `cond(N) > 1e14`: the printed sigmas are meaningless. Rescale the
  parameters so that the columns of the design have comparable norms, and use
  `Backend = "qr"`, which avoids squaring the condition number. A rank defect (a free
  network) needs a datum, which you declare as a prior (chapter *Priors and regularisation*).
* `PAR_CORR` reports the strongest correlation and its parameter pair. It is a NOTE at
  most 0.95, a WARN above 0.95 (nearly inseparable: use `R.ellipse(i, j)` rather than the
  individual sigmas), and an **ERROR** above 0.999. In the ERROR case only a combination of
  the two parameters is determined. Fix one of them, reparameterise to estimate the
  combination, or add observations that separate them.
* `PAR_TABLE` (NOTE) lists the first 30 parameters of the full vector: `x`, the a-priori
  sigma (`sigma0 · sqrt(Qxx_jj)`) and the scaled sigma (`s0 · sqrt(Qxx_jj)`). The scaled
  column is what `R.coefficients` reports for the commons.
* `PAR_SCALE` (WARN) appears when `|s0²/sigma0² − 1| > 0.1`, so that the two columns
  differ. Which column is honest depends on **why** the variance factor deviates. If the
  weights were wrong, the scaled column is right. If there are outliers or systematics,
  **neither** is right. Resolve Levels 1 and 2 before quoting either.

## Beyond the report: the stochastic model

The report judges one solve under the final weights. Whether those weights describe the
noise, stream by stream, is a separate question, and the report does not answer it. For an
`Estimate = "noise"` run, read it from `R.history`, `R.summary()` and `R.convergence()`.

### Whiteness per stream

On every pass the loop tests the law-whitened residuals of each stream with a Ljung–Box
test. The test pools over the stream's arcs, uses **K = 20 lags** and calls a stream white
when **p > 0.05**.

```matlab
T = R.summary();             % stream | sigma | arOrder | white | pValue  (final pass)
T(~T.white, :)               % the streams that never whitened
R.history.streams            % stream keys 'satKey|obsType', the column order of
R.history.white(end, :)      %   these per-pass matrices (nPass x nStream)
R.history.whitePval(end, :)
```

In `summary()`, `sigma` is the stream's level factor `sqrt(mean arc factor)`. It is near
1 when the declared sigma was right, and it tells you by how much the declared sigma was
off. `arOrder` is `NaN`, because the joint estimator fits a spectrum, not an AR order.

**A stream that is not white** has residual correlation that its law does not describe.
Possible causes and remedies, roughly in order:

1. **The law was not estimated for it.** If the stream's covariance spec has
   `estimate = "none"` or `"scale"`, the declared law is used unchanged (its
   `R.noise.model{k}.C` is then empty). Set `estimate = "law"` for that stream.
2. **It shares a law with streams that behave differently.** Use a finer `Pool`
   (`"obsType"`, or `"stream"` for one law per stream).
3. **The spectral grid is too coarse** for its correlation: raise `Noise.Grid`.
4. **The functional model is incomplete.** An unmodelled signal (a drift, a
   once-per-revolution term) appears as coloured noise. Run
   `R.systematic(Covariates=struct('time', t, ...))`: the binned-mean test points to the
   covariate. A correlation that lives only at a long lag (an echo, a repeat pattern) is
   invisible at K = 20; test it with `R.systematic(Lags=L)`.

`R.plot("psd")` draws the estimated spectral shape of each pool group, and
`R.plot("residuals", Stream="...")` draws a stream's whitened residuals against time.

### When the loop did not settle

The loop stops when two conditions hold. Every arc factor changed by at most
`Tolerance`. Every spectral bin j changed by at most `Tolerance·max(1, sqrt(200/r_j))`,
where `r_j` is the bin's redundancy and `sqrt(2/r_j)` the relative standard error of its
estimate. A bin with a redundancy of 200 or more is held to `Tolerance`. A bin with less,
such as the zero-frequency bin that per-arc parameters absorb, may move further, in
proportion to its own standard error: at the default `Tolerance` 1e-2, by a tenth of it.
A finer step would ask for more precision than its estimate has. A smaller `Tolerance`
tightens every bin. (Before 6.26.0 every bin was held to `Tolerance`. On the fixtures of the
development record that took 14–168 passes, while the solution was already within
0.1–0.2 of its standard error of the end point after 10.)

Each pass records `R.history.binRatio` (the largest change of a bin divided by its bound,
so ≤ 1 means every bin settled), `R.history.levelDev`, `R.history.maxDev` (the largest
`|α² − 1|`, a reading) and `R.history.xStep` (the largest change of the solution since the
pass before, in units of its own standard error). `R.convergence()` gives one row per pass
with `pass`, `maxSigmaOffset` (the largest `|sigma − 1|` of the per-stream level factors:
a reading, not the stop rule), `nWhite`, `minWhitePval`, `converged`, `maxDev`,
`levelDev`, `binRatio` and `xStep`. `R.convergence(Plot=true)` and
`R.plot("convergence")` draw the same record.

With `stopReason = "maxPasses"`, the warning names the bins that still move, with their
redundancy against the median of their group, and the solution's last step:

* If `xStep` is a few hundredths or less, the solution has settled to well within its
  errors, and the bins that still move are ones the data barely determine.
* If `binRatio` is still falling, the loop needs more passes: raise `MaxPasses`.
* If it oscillates, there is too little redundancy per component: coarsen `Pool`
  (`"obsType"` or `"all"`), use a smaller `Noise.Grid`, or smooth the spectrum
  (`Noise.Smooth`). A spectrum that keeps moving also points to too few arcs per frequency,
  or a start far from the data.
* If one stream stays non-white on every pass, see "Whiteness per stream" above.

### Skipped arcs and the warnings inventory

An arc that cannot be used in a pass is **skipped and keeps its previous factor**. The
reason is `'gap'` or `'nonPD'` (its covariance could not be factored). `R.history.skips`
lists every skip as `{pass, arcKey, reason}`. A pass that skips more than
`Noise.SkipThreshold` (0.2) of the arcs is refused with
`lsqvce:adjustVCEjoint:skipThreshold`. A few `nonPD` skips usually mean a `Cep` or a law
that is not positive definite on those arcs, so inspect that day's `buildDay` output.

Warnings are raised once per run. Muted or not, every warning identifier the run raised is
listed in `R.meta.warnings` (first-fired order). Check it in batch jobs. The ones that
qualify a result:

| identifier | meaning |
|---|---|
| `lsqvce:adjustVCEjoint:notConverged` | the stochastic model did not settle (the same fact as `NOISE_UNSETTLED`) |
| `lsqvce:adjustVCEjoint:arcFactorsDiverging` | the spread of a stream's arc factors kept growing over the passes: they are not an estimate. Take the one-pass estimate (`MaxPasses = 1`) or set `arcFactors = false` on that stream and read its level |
| `lsqvce:robust:atBreakdown` | a robust pass zeroed more than 45 % of the epochs (rows); the data are not a few blunders on a good model |
| `lsqvce:adjustRobust:notConverged`, `lsqvce:adjustNonlinear:maxIter` | IRLS or Gauss–Newton exhausted its budget (`stopReason = "maxPasses"`) |
| `lsqvce:systematic:lbUnreliable` | `R.systematic`'s Ljung–Box is anti-conservative at this `u/n`; lead with the binned test |

### Robust runs

A robust run reports `ROBUST_FIT` at Level 1. Its weights come back in `R.noise.robust`
when the robust pass runs inside the noise loop: `influence`, `scale` per stream, `factor`
per epoch in sweep order, `keys`, `truncation` (the factor the hard σ0 is corrected by) and
`frozen` (the pass a cycling rejected set was frozen at, 0 if none). For the stand-alone IRLS they are in
`R.noise.influence` and `R.noise.scale`, with per-pass `scale`, `dx`, `nDown` and `nZero`
in `R.history`. Down-weighted epochs are the robust counterpart of a flag list. Compare
them with a Pope or MAD snooping run before you remove anything permanently.

## Adverse findings at a glance

| you see | first suspect | do |
|---|---|---|
| `INVALID` (`SUM_R`, `VPV_*`, `NO_REDUNDANCY`) | setup or arithmetic | Make the data source deterministic, use `Backend = "qr"`, add redundancy. Read nothing else. |
| `NOISE_UNSETTLED`, `stopReason = "maxPasses"` | budget or poorly separable components | `R.convergence()`, then raise `MaxPasses`, loosen `Tolerance` or coarsen `Pool` / `Noise.Grid` |
| `NOISE_UNSETTLED`, `stopReason = "diverged"` | Gauss-Newton found no downhill step | check the start point `x0` and the Jacobian against a finite difference |
| `GLOBAL_HIGH`, flags present | outliers | `R.snoop()`, one at a time |
| `GLOBAL_HIGH`, no flags | systematics or optimistic weights | `R.systematic(Covariates=...)`; if clean and structureless, `Estimate = "noise"` or correct `sigma` |
| `GLOBAL_LOW` | pessimistic weights or over-parameterisation | correct `sigma` / `Estimate = "noise"`; `R.compare(Rred)`; Level-4 correlations |
| `GLOBAL_HIGH`/`LOW` with `GLOBAL_SHARP` | very large `r` | judge the ratio and the confidence interval, not the rejection |
| `SNOOP_MODERATE` / `SNOOP_CONTAMINATED` | contamination masks Pope | `Diagnostics.Test = "mad"`, `R.snoop(Test="mad")`, or `Robust` |
| `SNOOP_REFUSE` / `SNOOP_BREAKDOWN` | model error, or wrong sigma0 | reject nothing; fix the functional model or the declared sigma |
| `SNOOP_VACUOUS` | redundancy too small for the corrected level | weaker `Alpha`, or accept that blunders cannot be detected |
| `SNOOP_LOWR`, `REL_AMBIGUOUS` | weak geometry | re-measure, do not reject |
| `REL_UNCONTROLLED`, `REL_MDBINF`, `REL_SILENT`, `REL_WEAK` | untestable observations | strengthen the geometry; validate those observations externally |
| `REL_EXTERNAL` WARN/ERROR | one observation with excessive leverage | find it (minimum `r_i`, named observation); distrust the formal sigmas |
| `PAR_ILLCOND` | scaling or rank defect | rescale columns, `Backend = "qr"`, declare a datum prior |
| `PAR_CORR` ERROR | inseparable parameters | fix or reparameterise, or add separating observations |
| `PAR_SCALE` | variance factor ≠ 1 | resolve Levels 1–2 before quoting either sigma column |
| a stream not white in `R.summary()` | stream's law or the functional model | `estimate = "law"`, finer `Pool` / `Noise.Grid`, `R.systematic` |

## A worked reading

Every report has the same layout. The numbers below are placeholders; `demo_lsqvce` runs
six small adjustments, each with one failure planted (clean, outliers, a systematic drift,
an uncontrolled observation, wrong weights, collinear parameters), and prints the verdict
of each, so you can see the real output for each case.

```
==============================================================================
  LEAST-SQUARES ADJUSTMENT REPORT
  n = <n> observations,  u = <u> parameters,  r = n-u = <r>
  <b> batches,  vPv source: accumulator (v'v)
==============================================================================

--- LEVEL 0: INTEGRITY  (does the arithmetic hold?)
       run: lsqvce <version>, MATLAB R<release>, created <UTC time>
       sum(r_i) = <value> vs n-u = <r>  (rel. err <e>)  [identity holds]

--- LEVEL 1: GLOBAL     (is the model consistent?)
       The stochastic model settled in <k> pass(es).
       s0^2/sigma0^2 = <ratio>     T = vPv/sigma0^2 = <T>
       chi2 acceptance band (alpha=1%, r=<r>): [<lo>, <hi>]  width <w>
       sigma0_hat = <s0>,  95% CI [<a>, <b>]   (a priori sigma0 = 1.00000)
       ...
       Global test PASSES: functional and stochastic model are consistent.

--- LEVEL 2: LOCALISE   (where is it broken?)
       Scale cross-check: s0 = <s0>, robust scale = <m> (ratio <q>), ...
       Data snooping: POPE test, alpha = 1.0e-03 per test, Bonferroni-
       corrected to <alpha/n> for n = <n> simultaneous tests.
       Critical value |w| > <k>.   Expected false alarms: 0.00
       No observation exceeds the critical value. Largest |w| = <w> ...

--- LEVEL 3: RELIABILITY(what can you not check?)
       Redundancy numbers: min = <rmin> (observation <i>), mean = <r/n>
       MDB non-centrality delta0 = <d0>  (alpha_eff = <..>, beta = 0.20)
       <quadrant table>
       Hidden-blunder envelope -- ...  median <m>  90th pct <..>  ...
       Every observation has r_i >= 0.10: ...
       External reliability: the worst undetectable bias (the MDB of
       observation <i>) would shift parameter <j> by <e> of its own sigma.
       That is normal for a network of this size ...

--- LEVEL 4: PARAMETERS (what do you actually know?)
       cond(R) >= <..>     cond(N) >= <..>   [LOWER BOUND]
       ...
       Maximum parameter correlation is <c> (parameters <i> and <j>).
       <parameter table: par | x | sigma(a priori) | sigma(scaled)>

==============================================================================
  VERDICT: PASS        (0 error, 0 warning, <n> note)
  All checks passed: ...
==============================================================================
```

How to read it, in order:

1. **Level 0.** Check that `[identity holds]` is printed and that no `[ERROR]` appears.
   If one does, stop: the verdict is `INVALID`.
2. **Level 1, loop line.** If the run iterated, check that the stochastic model settled.
   If not, every later line is provisional.
3. **Level 1, the ratio.** Read the confidence interval on sigma0, not the pass/fail.
   A ratio well above 1 sends you to Level 2 (outliers) and to `R.systematic`.
4. **Level 2.** Read the scale cross-check before the flag list: if the check reports
   contamination, the list is not to be believed. With flags, read each `w` next to its
   `r_i`, and reject with `R.snoop`, one observation at a time.
5. **Level 3.** Check that the quadrant table has no UNTESTED or AMBIGUOUS counts and that
   no uncontrolled or infinite-MDB observations are reported. Check that the external
   reliability is normal for this `n`.
6. **Level 4.** Check that the conditioning is far from 1e10 and that no correlation is
   above 0.95. If `PAR_SCALE` fires, decide which sigma column applies only after Levels
   1 and 2 are resolved.
7. **Outside the report.** For a `noise` run, check `R.summary()` for non-white streams,
   `R.history.skips`, and `R.meta.warnings`.

Then quote the result: `R.coefficients()` for the commons, `R.toTable(data)` for every
parameter in your own table layout, and `R.meta` for the provenance.
