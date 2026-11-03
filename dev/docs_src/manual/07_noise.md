# Choosing the stochastic model

The weights of a least-squares solution are the inverse of the noise covariance. This
chapter says what to declare, what to let the data estimate, and how to read the
estimate back. The laws and their arguments are listed in the API reference (*Noise
laws*) and in `help lsqvce.Noise`; the mathematics is in the theory.

## Three facts, three places

| fact | declared in | estimated? |
|---|---|---|
| the **temporal law** — how the errors of a stream correlate from epoch to epoch | `addStream(..., Noise=lsqvce.Noise.<law>(...))` | yes, with `Estimate = "noise"`, as the stream's `Estimate` allows |
| the **level** of that law | the law's `sigma` (a standard deviation) | yes |
| the **per-epoch accuracy** — what the preprocessing knows about each observation | `addArc(..., Sigma=)` in `Build` | never: it is data |

The toolbox combines them per arc. A filter law (`white`, `arma`, an AR-whitened `psd`, `var`)
is applied in batches of up to 20 000 rows that overlap by the memory of the filter (for an
`arma` law with an MA part, until its impulse response has decayed to rounding), so the
batching does not change the result. A dense law (`acf`, `custom`, a `psd` whitened as
Toeplitz) is factorised over the whole arc: exact, at a cost of O(m³) and m² memory for an
arc of m rows — for arcs of many thousand epochs prefer a filter law. Above 20 000 epochs a
dense law is refused (`lsqvce:covModel:denseArcTooLong`, 3.2 GB for the factor alone at that
length); `MaxDenseArc=` on the law raises the limit deliberately. Under
`Estimate = "noise"` the estimator's sweep forms every arc's covariance densely whatever the
whitener, so there the same limit applies to every law, a filter law included.

## What to declare

| law | constructor | use when | cost per arc of m rows |
|---|---|---|---|
| white | `lsqvce.Noise.white(sigma)` | the colour is unknown (the default start), or truly white | O(m) |
| ARMA filter | `lsqvce.Noise.arma(ar, sigma, MA=ma)` | a whitening filter is known | O(m·p) |
| spectrum | `lsqvce.Noise.psd(f, S)` | a measured one-sided spectrum, f from 0 Hz | O(m³) Toeplitz, O(m·p) with `Noise.Whitener="ar"` |
| correlation | `lsqvce.Noise.acf(@rho, sigma)` | the correlation is known as a function of the lag | O(m³) |
| vector AR | `lsqvce.Noise.var(ar, Sigma)` | axes correlated in time across axes | O(m·p·d²) |
| any | `lsqvce.Noise.custom(@covFcn)` | anything else | O(m³) |

| stream's `Estimate` | under `Estimate = "noise"` | use for |
|---|---|---|
| `"law"` (default; `"scale"` for a per-row `sigma`) | shape and level estimated; the declared law gives only its `sigma`, as the starting level | the colour is unknown |
| `"scale"` | the declared shape kept, its level estimated | a trusted spectrum or filter of uncertain level |
| `"none"` | used exactly as declared | a reference stream, a stream you trust completely |

- **You do not know the colour of the noise.** Declare `lsqvce.Noise.white(sigma)` with a
  realistic `sigma` and run `Estimate = "noise"`: the estimator starts from white noise and
  finds the spectral shape. Under the stream's `Estimate = "law"` (the default) this is
  what happens to every stream, whatever law it declares. (For an `arma` law `sigma` is an
  innovation standard deviation; as a level it is read as the marginal one.)
- **You know the noise spectrum.** Declare `lsqvce.Noise.psd(f, S)` or an `arma` filter
  and set the stream's `Estimate = "scale"` to keep the shape and estimate only the level,
  or `"none"` to use it as given.
- **Your axes are correlated in time across axes.** Declare `lsqvce.Noise.var` (no
  per-arc `Sigma`). It is used as declared: `Estimate = "noise"` has no path for it and
  refuses it (`lsqvce:fit:varUnderNoise`).
- **A stream you trust completely** — a prior stream, a reference — takes
  `Estimate = "none"` on its `addStream`.

## Running the estimator

```matlab
opts = lsqvce.Options(Estimate="noise", Pool="obsType", ArcFactors=true, ...
                      MaxPasses=7, Seed=1, Noise=struct('Grid', 120, 'Whitener', "ar"));
R = lsqvce.fit(task, opts);
```

Per pass the toolbox solves under the current laws, sweeps every arc once for its
whitened residuals and their redundancies, and updates the spectral shape of every pool
group, the level of every stream, the factor of every arc and the σ of every prior with
`Estimate = "scale"`. It stops when no component changes by more than `Tolerance` — the
σ of a prior group included, which is first updated on the pass after it enters the solve,
so the loop does not stop before — or after `MaxPasses`. A spectral bin the data barely
determine may change by more, in proportion to its own standard error (chapter
*Reading the results and the report*, "When the loop did not settle").

The choices that matter:

| choice | guidance |
|---|---|
| `Pool` | the default `"stream"` gives every stream its own shape. Streams of the same physics may share one (`"obsType"`: both satellites' range-rate), but measured on identical spectra that gains nothing for x. Never pool streams of different physics or units: positions with range-rate drives the normal matrix towards singularity, and a shared law biases both levels (`poolNotWhite` says so) |
| `ArcFactors` | on when arcs are numerous and long enough to carry their own variance; off for few or short arcs |
| `Noise.Grid` | the spectral resolution of the estimated law; finer grids resolve narrow spectral features but give each bin less redundancy per pass. Keep 60 with the default `Smooth = 1`. Fewer bins cut the estimator's variance too, but bias the level low when the spectrum has a sharp low-frequency peak (16 bins: −12 % at AR(0.98)) |
| `Noise.Whitener` | how every `psd` law — the estimated one and any you declare — whitens an arc: `"ar"` for long arcs (cost O(m·p) instead of O(m³) in the solve); `"toeplitz"` (the default, exact) when the spectrum has narrow features a low-order filter cannot follow. Under `Estimate = "noise"` `"ar"` saves the solve only — the sweep forms each arc's covariance densely (above) — and the estimated spectrum is then whitened by an AR(30) that approximates it: exact enough for short-memory noise (an AR(1) law to about 1 %), biased where AR(30) cannot follow the shape (a 1/f floor: tens of percent per bin, and the loop need not settle) |
| `Probes`, `Seed` | the redundancies are estimated with random probes: 100 probes are the standard; a fixed `Seed` makes the run reproducible |

**Run one pass with `Estimate = "none"` first.** If σ₀ is near 1 the declared weights are
consistent. A σ₀ far from 1 says the `Sigma` scale is wrong — most often an information
matrix used where a covariance was expected — and estimating a law on top will not fix it.

## Reading the estimate

```matlab
R.noise.stream           % the stream keys, in the order of model{}
R.noise.model{k}.sigma   % the estimated level of stream k
R.noise.Sgroups          % the estimated spectrum of each pool group
R.summary()              % one row per stream: level, whiteness
R.plot("psd")            % the spectra
```

The standard deviation of one arc is the product of three terms: the declared `sigma`, the
square root of the arc factor, and `R.sigma0`. With `sigma = 1` and a `Sigma` carrying the
absolute accuracy, that product reads directly as *how much too optimistic the
preprocessing's covariance was on this arc*.

After convergence the whitened residuals of every stream should be white. The loop tests
that in every pass, per stream, with a Ljung–Box test over 20 lags at α = 0.05; the result
is in `R.summary()` and `R.history.white` (it is not a stop criterion). It detects short-lag
correlation, periodic terms and excess low-frequency noise; a correlation that exists only
at one long lag (an echo, a repeat-track pattern) is invisible to it — test for that with
`R.systematic('Lags', L)`.

Three assumptions of the model can be contradicted by the data, and the loop checks them
after convergence (`R.noise.checks`, one warning each):

| warning | what the data say | what to do |
|---|---|---|
| `lsqvce:adjustVCEjoint:poolNotWhite` | a stream of a pool group is not white under the group's shared law. Either the group mixes spectra, which biases both levels and the formal errors, or the stream has something unmodelled in it (an event, a drift) | compare `Pool = "stream"`; regroup, or find the event |
| `lsqvce:adjustVCEjoint:levelVaries` | a stream's arcs scatter far beyond their sampling error around one level: the level changes from arc to arc or day to day | arc factors for that stream |
| `lsqvce:adjustVCEjoint:memoryVsArc` | the estimated correlation time exceeds half the median arc: the arcs' local parameters absorb the long memory, and the level comes out low (−12 % at a correlation time of 200 against 256-epoch arcs) | longer arcs, or a law that models the low frequencies |

The formal errors of an estimated noise model are *plug-in*: `R.sol.Cxx` treats the
estimated covariance as known. On short or few arcs they are 3–13 % optimistic in σ; on long
records the optimism is negligible.

Arcs can be skipped in a pass: `gap` (no data) or `nonPD` (the arc's covariance is not
positive definite — a mis-scaled `Sigma` produces exactly this). A skipped arc keeps its
previous factor; when more than `Noise.SkipThreshold` of a stream's arcs are skipped, the
pass is refused rather than estimated on a fraction of the data.

## Outliers

Gross outliers belong to preprocessing. Inside the adjustment there are three tools:

1. **Arc factors** down-weight a whole bad arc automatically (`ArcFactors = true`).
2. **Robust estimation** down-weights single epochs: `Robust = "igg3"` or `"bisquare"`
   together with `Estimate = "noise"` weights every epoch of every stream inside the
   joint estimator, standardised by the stream's own robust scale; a rejected epoch is
   removed from the arc exactly. The statistic is Baarda's w for correlated observations
   — what deleting the epoch removes from v'Pv — tested one epoch at a time, the largest
   first (data snooping), and the truncation of the epochs that are kept is accounted
   for lag by lag. So the level, the law and σ0 of a robust run equal those of a plain
   run on clean data, on a coloured stream too (within 0.25 % in variance, measured),
   and a small blunder on a strongly correlated stream is found: at a lag-1 correlation
   of 0.98 a blunder of 2 standard deviations stands out by about 14 in w. Prefer a
   redescending function (`igg3`, `bisquare`) for gross blunders: Huber keeps part of
   their weight. The factors are in `R.noise.robust.factor`; a pass that rejects more
   than 45 % of the epochs warns that it has reached the breakdown point. A rejected set
   that flips back and forth between two passes is frozen at the union of the two
   (`R.noise.robust.frozen` names the pass).
3. **Diagnose, then decide.** Read the outlier tests in `R.report`, repair or remove the
   observations in `Build` and run again. Weights taken from a robust run can be made
   data: zero the rejected epochs' rows and mark them `Excluded` in the arc `Build`
   adds — the toolbox deletes such rows from a correlated law exactly, whole epochs or
   single axes (removing them
   leaves a gap that a `psd` or `arma` law, which counts lags in samples, takes for one
   step — `lsqvce:covModel:unevenGrid` says so); store the weights, never recompute them
   inside `Build`.

With `Estimate = "none"`, `Robust` runs a stand-alone iteratively re-weighted solve on the
whitened innovations, with its own budget `RobustOpts`.
