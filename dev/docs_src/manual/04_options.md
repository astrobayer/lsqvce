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
IRLS M-estimation on the same loop, every data kind: the influence function named here
reweights the standardised whitened residuals with a MAD scale each pass; the final
factors are `R.history(end).factor`. Not combinable with `Estimate="noise"` in one run
(refused by name: run the robust pass with `Estimate="none"`).

## The loop

**`MaxPasses`** — integer ≥ 1, default 10. The one pass budget of every iterative run
(variance components, IRLS, Gauss–Newton). Exhausting it returns the last pass with
`R.stopReason = "maxPasses"` and one warning.

**`MinPasses`** — integer ≥ 1, default 1. Passes run even if the tolerance is already met
(GRAVFIRE's `minimumIterations`).

**`Tolerance`** — positive, default 1e-2. The stop rule of the variance-component loop:
max |α² − 1| over frequencies and arcs below it; for Gauss–Newton the relative step in
units of the formal error; IRLS uses `RobustOpts.Tol`.

**`Probes`** — integer ≥ 0, default 100. Hutchinson probes for the per-lag redundancies
of the joint estimator; `0` = the exact route (measured to be the better default up to
u ≈ 4000; the probe route is faster only for the dense inverse at small nc).

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

**`Verbose`** — logical or 0..3, default `false`. Screen messages; 0 is silent.

**`LogFile`** — string, default `""`. Mirror the messages to this file.

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
| `Grid` | 60 | frequencies (lags) of the estimated shape; must be ≥ the longest arc for an exact law (`"auto"` = the longest arc) |
| `Whitener` | `"toeplitz"` | how a PSD law whitens an arc: exact Toeplitz Cholesky, or an AR / ARMA filter fitted to it (`"ar"`, `"arma"`: O(m·p) instead of O(m²)) |
| `MeanType` | `"robust"` | the normalising mean of the arc factors: trimmed (`"robust"`) or plain |
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
