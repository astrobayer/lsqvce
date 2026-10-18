# `lsqvce.Options`

`lsqvce.Options` is a typed value class: every knob has a type, a range and a default,
checked at assignment (`MATLAB:validators:*`) or at construction (`lsqvce:Options:*`).
It says *how to run*; everything about the problem, including regularisation, is
declared on `data`. Three equivalent spellings:

```matlab
opts = lsqvce.Options(Estimate="noise", MaxPasses=6);   R = lsqvce.fit(data, opts);
opts.Noise.Grid = 120;                                   % change one later (validated)
R = lsqvce.fit(data, Estimate="noise", MaxPasses=6);     % name-value pairs straight on fit
```

`disp(opts)` lists every knob, marking the defaults with `*`.

## What to estimate

| knob | values | default | meaning |
|---|---|---|---|
| `Estimate` | `"none"`, `"noise"` | `"none"` | `"none"`: one generalised least-squares solve under the declared laws. `"noise"`: the joint variance-component estimator — per pool group the spectral shape of the noise, per stream its level, per arc a factor (with `ArcFactors`), and the σ of every prior group with `estimate = "scale"`. |
| `Pool` | `"all"`, `"obsType"`, `"stream"` | `"all"` | which streams share one estimated spectral shape: all streams of the run, the streams of one observation type across systems, or every (system, type) stream on its own. Every stream keeps its own level. |
| `ArcFactors` | logical | `true` | one variance factor per arc on top of the stream's level; the factor weights the solve. `false`: one level per stream. Turn off when arcs are few or short. |
| `Robust` | `"none"`, `"huber"`, `"igg3"`, `"bisquare"` | `"none"` | M-estimation. With `Estimate = "none"`: stand-alone IRLS on the law-whitened innovations with a MAD scale per pass (budget `RobustOpts`, factors in `R.history(end).factor`). With `Estimate = "noise"`: a robust weight per epoch inside the joint estimator, on the residuals descaled by the declared accuracy (budget `MaxPasses`, factors in `R.noise.robust.factor`). Prefer `"igg3"` or `"bisquare"` for gross blunders. |

## The loop

| knob | values | default | meaning |
|---|---|---|---|
| `MaxPasses` | integer ≥ 1 | 10 | pass budget of the variance-component loop and of Gauss–Newton (IRLS uses `RobustOpts.MaxIter`). Exhausting it returns the last pass with `R.stopReason = "maxPasses"` and one warning. |
| `Tolerance` | > 0 | 1e-2 | stop rule of the variance-component loop: every arc factor's change ≤ `Tolerance`, every spectral bin's change ≤ `Tolerance`·max(1, sqrt(200/r)), r its redundancy (at the default a tenth of the bin's standard error sqrt(2/r) where that is larger); for Gauss–Newton the relative step in units of the formal error |
| `Probes` | integer ≥ 1 | 100 | Hutchinson probes for the per-lag redundancies of the joint estimator |
| `Seed` | `[]` or double | `[]` | seed for the probes, recorded in `R.meta.seeds` |

## Run control

| knob | values | default | meaning |
|---|---|---|---|
| `Backend` | `"qr"`, `"normal"` | `"qr"` | square-root information (Householder QR) or normal equations (Cholesky); the same numbers to rounding, QR the safer on ill-conditioned common blocks |
| `CacheDays` | > 0 | 1 | built days kept resident per system; `Inf` keeps the calendar and builds each day once per run (a built type-day is m × nc doubles); a system's `prepareDay` result is the current day's only, whatever `CacheDays` |
| `Verbose` | logical | `false` | screen messages: 0 silent, 1 progress |
| `OnPass` | `[]` or `@(ctx)` | `[]` | inspection hook called after every pass of every iterative run (below) |

### `OnPass`

The hook receives one struct: `driver`, `pass`, `x`, `Cxx`, `sigma0`, `s0sq`; the loop
state (noise loop: `C`, `Cgroups`, `sigma2`, `arcKeys`, `streams`, `sigma`, `white`,
`whitePval`, `maxDev`, `binRatio`, `xStep`, `skips`; IRLS: `factor`, `keys`, `scale`, `dx`; Gauss–Newton:
`vPv`, `dxRel`, `lambda`); `acc`, the live accumulator; `data`, the declaration as
passed to `fit`; `processingOptions{k}`, system k's as declared. Return nothing to
continue, `false` to stop after the pass (`R.stopReason = "userStop"`). An error in the
hook arrives as `lsqvce:callUser:threw` with your error as its cause.

## The groups

A group is a struct; assign a whole or partial struct, or one field. A partial struct
is merged onto the defaults; an unknown field is refused
(`lsqvce:Options:badGroupField`).

### `Noise` — the estimated stochastic model

| field | default | meaning |
|---|---|---|
| `Grid` | 60 | spectral bins of the estimated law between 0 and Nyquist (its resolution), an integer ≥ 2; the law is exact at every lag of every arc whatever the Grid |
| `Whitener` | `"toeplitz"` | how every `psd` law (the estimated one and the declared ones) whitens an arc: exact Toeplitz Cholesky over the whole arc, or a fitted AR / ARMA filter (`"ar"`, `"arma"`; any case; another value is refused, `lsqvce:Options:badWhitener`) |
| `Smooth` | 0 | moving-average half-width on the estimated spectrum (0 = none) |
| `FreezeRedundancy` | 1 | a bin whose redundancy is below this keeps its value instead of being updated |
| `SkipThreshold` | 0.2 | the share of skipped arcs of a stream above which a pass is refused |

### `RobustOpts` — IRLS

| field | default | meaning |
|---|---|---|
| `Tuning` | `[]` | tuning constant of the influence function; empty = Huber 1.345, IGG3 (1.5, 2.5), bisquare 4.685 |
| `MaxIter` | 30 | IRLS passes |
| `Tol` | 1e-8 | IRLS stops when every parameter moved by less than `Tol` of its own formal sigma between two passes |

### `Diagnostics` — the per-observation tests

| field | default | meaning |
|---|---|---|
| `Alpha` | 1e-3 | significance level of the outlier test |
| `Bonferroni` | `true` | divide α by the number of tests |
| `Beta` | 0.2 | type-II error; the power of the minimal detectable bias is 1 − β |
| `RMin` | 0.01 | redundancy number below which an observation is untestable |
| `RWeak` | 0.1 | redundancy number below which it is weakly controlled |
| `Test` | `"auto"` | `"pope"` (σ₀ estimated), `"baarda"` (σ₀ known), `"mad"` (robust scale), or automatic |
| `RobustCheck` | `true` | run the MAD-scale check that flags a contaminated fit |
| `Store` | `"all"` | keep the per-observation arrays for all observations, only `"flagged"` ones, or `"none"` |
| `NWorst` | 1000 | how many worst observations `R.report` keeps |

### `Parallel`

| field | default | meaning |
|---|---|---|
| `Mode` | `"serial"` | `"parfor"` or `"spmd"` with the Parallel Computing Toolbox; serial with one warning without it |
| `DistributeDays` | `false` | distribute the day loop of the accumulation over workers (bitwise-identical fold order) |
| `Threshold` | 8 | number of days below which distribution is not attempted |

### `GaussNewton` — problems with `data.x0`

| field | default | meaning |
|---|---|---|
| `TolF` | 1e-10 | relative change of the objective that ends the iteration |
| `Lambda0` | 1e-6 | initial Marquardt damping (rows `sqrt(λ·diag(JᵀPJ))·dx ≈ 0`); a rejected step multiplies it by 10 (from 0: 1e-3), an accepted one divides it (below 1e-12: 0) |
