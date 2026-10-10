# Options — the complete reference

`lsqvce.Options` is a typed value class: every knob has a type, a range and a default,
checked at assignment (`MATLAB:validators:*` identifiers) or at construction
(`lsqvce:Options:*`). It says *how to run*; everything about the data, including
regularisation, is declared on `data` (chapter 3). Three equivalent spellings:

```matlab
opts = lsqvce.Options(Estimate="noise", MaxPasses=6);   R = lsqvce.fit(data, opts);
opts.Noise.Grid = 2160;                                  % change one later (validated)
R = lsqvce.fit(data, Estimate="noise", MaxPasses=6);     % straight on fit
```

`disp(opts)` lists every knob, a `*` marking the defaults. The option space is
**14 top-level knobs and 5 groups**; the count is pinned by a test so it cannot grow
unnoticed.

## What to estimate

**`Estimate`** — `"none"` (default) | `"noise"`.
`"none"` is one generalised least-squares solve under the declared covariance laws.
`"noise"` runs the joint variance-component estimator (Ellmer 2018, theory ch. 6): for
every stream whose spec allows it (`covariance.<type>.estimate ~= "none"`), the spectral
*shape* of its noise is estimated per pool group from Hutchinson-probed redundancies, one
*level* per stream, and (with `ArcFactors`) one factor per arc; prior groups with
`estimate="scale"` get their sigma by Helmert / Koch–Kusche in the same passes. The v4
words map as: `"ellmer"` → `Estimate="noise"`; `"vce"` → `Estimate="noise", Pool="stream",
ArcFactors=false`; `"pooled"` → `Estimate="noise", Pool="all", ArcFactors=false`.

**`Pool`** — `"all"` (default) | `"obsType"` | `"stream"`.
Which streams *share* one estimated spectral shape: all streams of the run; the streams
of one observation type across systems (both satellites' KBR, say); or every
(system, type) stream on its own. Every stream keeps its own level regardless. Fewer
groups = more arcs per shape = a steadier estimate; more groups = more physics.

**`ArcFactors`** — logical, default `true`.
Per-arc variance factors (theory 6.5.4): every arc of a stream gets a keyed factor on top
of the stream's level, normalised to unit mean over the run. `false` = one level per
stream (the Helmert factor of its arcs). Turn off when arcs are few or short.

**`Robust`** — `"none"` (default) | `"huber"` | `"igg3"` | `"bisquare"`.
M-estimation on the same loop, every data kind, in two deliberately different conventions
(chapter 7, §7.8). With `Estimate="none"` it is the stand-alone IRLS (`adjustRobust`): the
influence function re-weights the law-whitened *innovations* (Yang) with a MAD scale each
pass, `RobustOpts.MaxIter` / `Tol` are its budget, the final factors are
`R.history(end).factor`. With `Estimate="noise"` (since 6.1) the robust pass runs *inside* the
joint estimator and weights *observations* (Koch): the raw residuals, descaled by each
epoch's declared accuracy, one factor per epoch, folded into `Cep` before the temporal law;
`MaxPasses` is the budget, the factors are `R.noise.robust.factor`. Prefer `"igg3"` or
`"bisquare"` for gross blunders — Huber keeps c/|u| of their weight.

## The loop

**`MaxPasses`** — integer ≥ 1, default 10. The pass budget of the variance-component
loop and of Gauss–Newton. **IRLS has its own**, `RobustOpts.MaxIter` (30), as it has its
own tolerance — `engineArgs` hands `adjustRobust` that one and never `MaxPasses`
(R6 WP-101). Exhausting the budget returns the last pass with
`R.stopReason = "maxPasses"` and one warning.

**`Tolerance`** — positive, default 1e-2. The stop rule of the variance-component loop:
max |α² − 1| over frequencies and arcs below it; for Gauss–Newton the relative step in
units of the formal error; IRLS uses `RobustOpts.Tol`.

**`Probes`** — integer ≥ 1, default 100 (Ellmer's value). Hutchinson probes for the
per-lag redundancies of the joint estimator. There is **no exact route** in that
estimator, so 0 is refused at assignment (R6 WP-101; the `'TraceMethod','exact'` route
exists only on the engine-tier `lsqvce.engine.redundancy`, which `fit` does not call).

**`Seed`** — `[]` (default) or a double. RNG seed for the probes; recorded in
`R.meta.seeds` so a run is reproducible.

## Run control

**`Backend`** — `"qr"` (default) | `"normal"`. Square-root (Householder QR, rank refusal
by `rank()`'s tolerance) or normal equations (Cholesky). Same numbers to rounding; QR is
the safer default on ill-conditioned common blocks.

**`CacheDays`** — positive, default 1. Built days kept resident per system. A two-stage
pass builds every day twice (accumulate + sweep) and a VCE run 2n+1 times; `Inf` keeps the
calendar and builds each day once per run — your memory decision (a built type-day is
m × nc doubles, ≈190 MB at degree 60 for a 6-h three-axis piece).

**`Verbose`** — logical, 0 or 1, default `false`. Screen messages: 0 silent, 1 progress.
This surface is a FLAG: `Options.engineArgs` hands every driver `logical(Verbose)`, so
level 2 has never been reachable through `fit` — it is an engine-tier setting of
`lsqvce.util.verbosity` (R6 WP-101). A 2 is still accepted here and behaves as 1.
(Level 3 was merged into 2 in 2026-07; `lsqvce.util.verbosity` refuses it, and since v5.3 so does
`Options` itself, by name.)

**`OnPass`** — `[]` (default) or `@(ctx)`. Your inspection hook, called at the end of
**every pass** of every iterative run with one struct: `driver`, `pass`, `x`, `Cxx`,
`sigma0`, `s0sq`; the loop state (noise: `C`, `Cgroups`, `sigma2`, `arcKeys`, `streams`,
`sigma`, `white`, `whitePval`, `maxDev`, `skips`; IRLS: `factor`, `keys`, `scale`, `dx`;
Gauss–Newton: `vPv`, `dxRel`, `lambda`); `acc` — the **live accumulator** (call
`lsqvce.engine.residuals(ctx.acc, …)` or any engine function inside the hook); `data` —
your declaration as passed to `fit`; `processingOptions{k}` — system k's, as declared.
Plot, `drawnow`, save wherever you like. Return nothing to continue, `false` to stop after
the pass (`R.stopReason = "userStop"`, the last pass returned with its covariance). A
throwing hook is not swallowed: it arrives as `lsqvce:callUser:threw` with your error as
the cause. This replaces the v4 `plotPostSolve` / `plotPostModel` hooks and `ShowSeconds`
/ `ExportDir`.

## The groups

A group is a struct; assign a whole (partial) struct or one field. Partial structs are
merged onto the defaults; an unknown field is refused (`lsqvce:Options:badGroupField`).

### `Noise` — the estimated stochastic model

| field | default | meaning |
|---|---|---|
| `Grid` | 60 | frequencies (lags) of the estimated shape; an integer ≥ 2, and ≥ the longest arc for an exact law. `"auto"` was documented until v5.3 and never implemented; it is refused now (R6 WP-101) |
| `Whitener` | `"toeplitz"` | how a PSD law whitens an arc: exact Toeplitz Cholesky, or an AR / ARMA filter fitted to it (`"ar"`, `"arma"`: O(m·p) instead of O(m²)) |
| `MeanType` | `"robust"` | gamma-bar, the mean of the arc factors: trimmed (`"robust"`) or plain. Since 6.1.6 the factors are absolute and gamma-bar sets only the scale of their clamp floor (§7.8 item 3) |
| `Smooth` | 0 | moving-average half-width on the estimated PSD (0 = none) |
| `FreezeRedundancy` | 1 | a frequency whose redundancy is below this holds its value instead of being updated |
| `SkipThreshold` | 0.2 | the skip *rate* (arcs the sweep could not process / arcs of the stream) above which the pass is refused rather than silently thinned |

### `RobustOpts` — IRLS

| field | default | meaning |
|---|---|---|
| `Tuning` | `[]` | the influence function's tuning constant (its textbook default when empty: Huber 1.345, IGG3 (1.5, 2.5), bisquare 4.685) |
| `MaxIter` | 30 | IRLS passes |
| `Tol` | 1e-8 | relative change of x that ends the iteration |

### `Diagnostics` — the per-observation tests (`R.diag`, `R.report`)

| field | default | meaning |
|---|---|---|
| `Alpha` | 1e-3 | significance level of the outlier test |
| `Bonferroni` | `true` | correct α for the number of tests |
| `Beta` | 0.2 | power (1 − β) for the minimal detectable bias |
| `RMin` | 0.01 | redundancy number below which an observation is "untestable" |
| `RWeak` | 0.1 | below which it is "weakly controlled" |
| `Test` | `"auto"` | `"pope"` (σ₀ estimated) / `"baarda"` (σ₀ known) / `"mad"` (robust scale) / auto by the a-priori backstop |
| `RobustCheck` | `true` | run the MAD-scale backstop that flags a contaminated fit |
| `Store` | `"all"` | keep per-observation arrays for all / `"flagged"` / `"none"` |
| `NWorst` | 1000 | how many worst observations `R.report` keeps |

### `Parallel`

| field | default | meaning |
|---|---|---|
| `Mode` | `"serial"` | `"parfor"` or `"spmd"` with the Parallel Computing Toolbox; serial fallback with one warning without it |
| `DistributeDays` | `false` | the day loop of the accumulation over workers (bitwise-identical fold order) |
| `Threshold` | 8 | days below which distribution is not attempted |

### `GaussNewton` — problems with `data.x0`

| field | default | meaning |
|---|---|---|
| `TolF` | 1e-10 | relative change of the objective that ends the iteration |
| `Lambda0` | 1e-6 | initial Levenberg–Marquardt damping; a rejected step multiplies it |

## Retired knobs

Regularisation is no longer an option (`Reg`, `RegLambda`, `RegMatrix`, `DatumB`,
`Sigma0`, `Sparse` are gone): declare `data.priors` / `data.systems(k).priors` (chapter 6).
The full v4 → v5 map is chapter 12.
