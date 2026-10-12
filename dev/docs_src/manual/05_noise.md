# Choosing the stochastic model

The weights of a least-squares solution are the inverse of the noise covariance. This
chapter says what to declare, what to let the data estimate, and how to read the
estimate back. The laws and their fields are listed in the API reference (*Covariance
specifications*); the mathematics is in the theory.

## Three facts, three places

| fact | declared in | estimated? |
|---|---|---|
| the **temporal law** — how the errors of a stream correlate from epoch to epoch | `covariance.<type>.type` and its parameters | yes, with `Estimate = "noise"` |
| the **level** of that law | `covariance.<type>.sigma` (a standard deviation) | yes |
| the **per-epoch accuracy** — what the preprocessing knows about each observation | `D.Cep`, returned by `buildDay` | never: it is data |

The toolbox combines them per arc. An arc is whitened in batches of up to 20 000 rows. For
a filter law (`diag`, `arma`, an AR-whitened `psd`, `var`) consecutive batches overlap by the
memory of the filter, so batching does not change the result. A dense law (`acf`,
`custom`, a `psd` whitened as Toeplitz) is factorised per batch: on an arc longer than one
batch the correlation between batches is dropped — cut such arcs, or use a filter law.

## What to declare

- **You do not know the colour of the noise.** Declare `type = 'diag'` with a realistic
  `sigma` and run `Estimate = "noise"`: the estimator starts from white noise and finds
  the spectral shape. Under `estimate = "law"` (the default) this is what happens to every
  stream, whatever law it declares: the declared law contributes only its `sigma`, as the
  starting level, and the shape is estimated from white noise. (For an `arma` law `sigma`
  is an innovation standard deviation; as a level it is read as the marginal one.)
- **You know the noise spectrum.** Declare `type = 'psd'` with `f` and `S`, or an `arma`
  filter; set `estimate = "scale"` to keep the shape and estimate only the level, or
  `"none"` to use it as given.
- **Your axes are correlated in time across axes.** Declare a `var` law (no `Cep`).
- **A stream you trust completely** — a prior stream, a reference — takes
  `estimate = "none"`.

## Running the estimator

```matlab
opts = lsqvce.Options(Estimate="noise", Pool="obsType", ArcFactors=true, ...
                      MaxPasses=7, Seed=1, Noise=struct('Grid', 120, 'Whitener', "ar"));
R = lsqvce.fit(data, opts);
```

Per pass the toolbox solves under the current laws, sweeps every arc once for its
whitened residuals and their redundancies, and updates the spectral shape of every pool
group, the level of every stream, the factor of every arc and the σ of every prior with
`estimate = "scale"`. It stops when no component changes by more than `Tolerance`, or
after `MaxPasses`.

The choices that matter:

| choice | guidance |
|---|---|
| `Pool` | streams of the same physics may share a shape (`"obsType"`: both satellites' range-rate); never pool streams of different physics or units — positions with range-rate drives the normal matrix towards singularity |
| `ArcFactors` | on when arcs are numerous and long enough to carry their own variance; off for few or short arcs |
| `Noise.Grid` | the spectral resolution of the estimated law; finer grids resolve narrow spectral features but give each bin less redundancy per pass |
| `Noise.Whitener` | `"ar"` for long arcs (cost O(b·p) instead of O(b³)); `"toeplitz"` when the spectrum has narrow features a low-order filter cannot follow |
| `Probes`, `Seed` | the redundancies are estimated with random probes: 100 probes are the standard; a fixed `Seed` makes the run reproducible |

**Run one pass with `Estimate = "none"` first.** If σ₀ is near 1 the declared weights are
consistent. A σ₀ far from 1 says the `Cep` scale is wrong — most often an information
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
square root of the arc factor, and `R.sigma0`. With `sigma = 1` and a `Cep` carrying the
absolute accuracy, that product reads directly as *how much too optimistic the
preprocessing's covariance was on this arc*.

After convergence the whitened residuals of every stream should be white. The loop tests
that in every pass, per stream, with a Ljung–Box test over 20 lags at α = 0.05; the result
is in `R.summary()` and `R.history.white` (it is not a stop criterion). It detects short-lag
correlation, periodic terms and excess low-frequency noise; a correlation that exists only
at one long lag (an echo, a repeat-track pattern) is invisible to it — test for that with
`R.systematic('Lags', L)`.

Arcs can be skipped in a pass: `gap` (no data) or `nonPD` (the arc's covariance is not
positive definite — a mis-scaled `Cep` produces exactly this). A skipped arc keeps its
previous factor; when more than `Noise.SkipThreshold` of a stream's arcs are skipped, the
pass is refused rather than estimated on a fraction of the data.

## Outliers

Gross outliers belong to preprocessing. Inside the adjustment there are three tools:

1. **Arc factors** down-weight a whole bad arc automatically (`ArcFactors = true`).
2. **Robust estimation** down-weights single epochs: `Robust = "igg3"` or `"bisquare"`
   together with `Estimate = "noise"` weights every epoch of every stream inside the
   joint estimator, standardised by the stream's own robust scale; a rejected epoch is
   removed from the arc exactly. Prefer a redescending function (`igg3`, `bisquare`) for
   gross blunders: Huber keeps part of their weight. The factors are in
   `R.noise.robust.factor`; a pass that rejects more than 45 % of the epochs warns that
   it has reached the breakdown point.
3. **Diagnose, then decide.** Read the outlier tests in `R.report`, repair or remove the
   observations in `buildDay` and run again. Weights taken from a robust run can be made
   data: remove rejected epochs (a gap) in `buildDay`; store the weights, never recompute
   them inside `buildDay`.

With `Estimate = "none"`, `Robust` runs a stand-alone iteratively re-weighted solve on the
whitened innovations, with its own budget `RobustOpts`.
