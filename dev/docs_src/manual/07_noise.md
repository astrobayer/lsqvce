# Stochastic models and their estimation

This chapter is the reference for the covariance specification: what you declare, what
travels with the data, how the two are composed, what the estimator does with them, and how
to read the result back. Chapter 3 lists the fields as a contract; this chapter says what
they mean and how to choose them.

## 7.1 The division of labour — three facts, three places

The noise of one observation stream is three separate facts, and each is stated in exactly
one place:

| fact | where it is stated | who owns it |
|---|---|---|
| the **temporal law** R_t — how the errors of a stream correlate down the epochs | `covariance.<stream>.type` and its parameters | you declare a starting point; the estimator refines it |
| the **level** of that law | `covariance.<stream>.sigma` (always a standard deviation) | you declare it; the estimator scales it |
| the **per-epoch accuracy** D — what the preprocessing knows about each observation | `D.Cep`, returned by `buildDay` with the day (chapter 3, §7) | you, always; it is never estimated |

The toolbox composes them per arc as **C = D (R_t ⊗ I_d) D**: every epoch is de-scaled by
D^-1/2, then the temporal whitener runs down the epochs of each axis. Batches are
overlapped by the filter memory automatically, so cutting an arc into batches never changes
the result.

`Cep` is **relative to** `sigma`. A `diag` law of level s with a per-epoch sigma c whitens
as `l ./ (s*c)`. If `Cep` already carries the absolute accuracy — a PPP covariance in m²,
for instance — declare `sigma = 1`, or the accuracy is counted twice.

## 7.2 `Cep` — the three shapes, and the asymmetry that bites

`D.Cep` is a **cell like `A`**, one entry per arc (v6). Each entry may be:

| shape | meaning | what the numbers are |
|---|---|---|
| nEpoch × 1 | a scalar stream with per-epoch accuracy | **standard deviations** |
| nEpoch × d | a multi-axis stream, axes uncorrelated within an epoch | **standard deviations** per axis |
| nEpoch × d × d | a multi-axis stream, full per-epoch covariance | **covariance matrices** (the toolbox takes their Cholesky factor) |

The vector forms are sigmas and the matrix form is a covariance. That asymmetry is the
commonest way to be wrong by a square: check one epoch by hand before a campaign.

Rules:

- the axis count is `d = size(Cep, 2)`, derived, and the same for the whole day; a `nAxis`
  given as well must agree (`axisConflict`);
- a `Cep` on the covariance spec *and* on the day is refused (`cepTwice`) — one fact, one
  authority. The spec field is the engine-tier carrier the loader fills from the day;
- a `var` law takes **no** `Cep` (`cepOnVar`): it carries its own cross-axis structure;
- a multi-axis day **without** `Cep` is composed with D = I and reported once per run
  (`noCep`), so an omission is distinguishable from a choice;
- rows of a multi-axis stream are **epoch-major** (x₁ y₁ z₁ x₂ y₂ z₂ …) and every arc holds
  whole epochs. The temporal law is placed on the **epoch** grid, `t(1:d:end)`.

## 7.3 The temporal laws

| `type` | fields | when to use it | cost per batch |
|---|---|---|---|
| `diag` | `sigma` | white noise — and the honest **starting point** when the colour is unknown: `Estimate="noise"` discovers it | O(b) |
| `arma` | `ar = [1 a₁ … a_p]`, `ma` (default 1), `sigma` (innovation std) | a known whitening filter | O(b·p) |
| `psd` | `f` (Hz, one-sided), `S` (unit²/Hz) | noise given as a measured spectrum | O(b·p) as AR, O(b³) as Toeplitz |
| `acf` | `acf = @(dt, opt)`, `sigma` | correlation as a function of the lag | O(b³) |
| `var` | `ar` (d×d×p), `sigma` (d×d) | genuine **cross-axis temporal** correlation | O(b·p·d²) |
| `custom` | `fcn = @(t, opt)` | anything else, dense | O(b³) |

`sigma` is a standard deviation, never a variance.

**A PSD is best whitened as an AR filter, and that is exact.** Fitting an AR(p) to the
spectrum and whitening with the FIR filter is not an approximation of the dense Cholesky:
the filter *is* the Cholesky inverse factor, and both are lower-triangular square roots of
the inverse covariance, unique up to sign, agreeing to ~1e-14 row for row. Set
`Noise.Whitener = "ar"` and pay O(b·p) instead of O(b³). Validate the fit once per stream
with `lsqvce.engine.checkPSD` (aim for `acvRelErr < 1e-3`; raise `arOrder`, or fall back to
`"toeplitz"` if it does not reach that).

## 7.4 The two policy fields on a spec

Beside the law, a spec may carry:

| field | values | meaning |
|---|---|---|
| `estimate` | `"law"` (default), `"scale"`, `"none"` | `"law"`: the spectral shape **and** the level are estimated; `"scale"`: your declared shape is kept, only the level is estimated; `"none"`: the stream is left exactly as declared |
| `arcFactors` | logical | per-arc variance factors for this stream (Ellmer §6.5.4); default: the `Options.ArcFactors` policy |

One refusal to know about: a **per-row `sigma` vector** makes the spec a *fixed structure*.
It can carry a level but not a spectral shape, so `estimate = "law"` on it is
`lsqvce:fit:lawOnFixedStructure`. A `Cep` does **not** trigger this — per-epoch accuracy
together with an estimated temporal law is exactly the supported combination, and the
reason the error information travels with the day at all.

`Options.Pool` decides which streams **share one estimated spectral shape**; every stream
keeps its own level regardless. `"all"` (the default) pools every stream, `"obsType"`
groups by stream name, `"stream"` gives each (system, stream) its own shape. With a single
system, `"obsType"` and `"stream"` are the same grouping. Pooling streams of different
physics is a mistake worth naming: pooling positions with range-rate was measured to drive
the conditioning to RCOND 1e-50 on the gravity demo.

## 7.5 What `Estimate="noise"` does

The joint estimator (theory chapter 6; Ellmer 2018 §6.5.7) treats every PSD frequency and
every arc factor as variance components of **one** Helmert VCE and updates them together.
Per pass it

1. accumulates and solves under the current laws,
2. sweeps every arc once for its whitened residual and its Hutchinson-probed redundancies
   at every lag (all lags from one FFT per column),
3. updates the spectral shape per pool group in the PSD domain (positive by construction,
   resampled to each arc's length), the level per stream, the per-arc factors, and the sigma
   of any prior group with `estimate = "scale"`,

and stops at `Tolerance` (max |alpha² − 1|) or `MaxPasses`. `Noise.Grid` is the number of
lags/frequencies of the estimated shape: at least the longest arc for an exact law.
`Probes` is the Hutchinson probe count (100 is Ellmer's value); the joint estimator has no
exact route, so `Probes = 0` is refused.

Arcs can be **skipped** for a pass with reason `gap` (no data) or `nonPD` (the arc's
covariance is not positive definite — a mis-scaled `Cep` produces exactly this). A skipped
arc keeps its previous factor and is listed in the run's skip record; beyond
`Noise.SkipThreshold` (0.2) the pass is refused rather than estimated on a fraction of the
data.

For multi-axis streams the residuals are de-scaled by D **before** the lag statistics are
formed, per axis over epochs. The estimated law is therefore the law of the *normalised*
errors, which is what makes a 3-axis stream with a per-epoch 3×3 and a scalar stream with a
per-epoch sigma behave identically.

## 7.6 Worked example — kinematic orbits with a 3×3 per epoch, KBR with one variance

The GRACE-FO case: two position streams whose preprocessing delivers a full 3×3 covariance
per epoch, one range-rate stream with a single scalar variance, all in one system (the
range-rate rows depend on both satellites' local parameters), and Ellmer's estimator.

```matlab
% ---- the specs: the TEMPORAL law and its level, nothing else ----------
% positions: Cep carries the absolute accuracy -> level 1; white start, so
% the estimator finds the colour; per-arc factors on
cov.posA = struct('type','diag', 'sigma', 1, 'estimate', "law");
cov.posB = struct('type','diag', 'sigma', 1, 'estimate', "law");
% range-rate: the scalar variance IS the level; no Cep on the day
cov.kbr  = struct('type','diag', 'sigma', sqrt(varianceKbr), 'estimate', "law");
s.covariance = cov;

% ---- buildDay: the error information travels with the day ------------
%   D.Cep{k} = nEp_k x 3 x 3     the per-epoch covariance of arc k
%   (rows epoch-major; the kbr stream returns no Cep field at all)

% ---- the run ---------------------------------------------------------
opts = lsqvce.Options( ...
        Estimate   = "noise", ...          % the joint (Ellmer) estimator
        Pool       = "obsType", ...        % one shape per stream, NOT one for all three
        ArcFactors = true, ...             % Ellmer eq. 6.5.4
        MaxPasses  = 7, ...
        Probes     = 100, Seed = 1, ...    % reproducible probes
        Noise      = struct('Grid', 2160, 'Whitener', "ar"), ...
        CacheDays  = Inf, ...
        Verbose    = 1);
R = lsqvce.fit(data, opts);
```

Why each line: `sigma = 1` on the positions because the 3×3 already carries m²;
`Pool = "obsType"` because positions and range-rate differ in units and in spectrum;
`Grid = 2160` because that is the length of the longest arc in epochs (6 h at 10 s, 3 h at
5 s); `Whitener = "ar"` because a 2160-epoch arc on the dense Toeplitz path is O(b³).

**Before trusting the numbers, run one pass with `Estimate = "none"`.** If sigma0 is near 1
the declared weights are consistent. A sigma0 far from 1 says the `Cep` scale is wrong —
most often an information matrix used where a covariance was expected — and no amount of
estimating a law on top of it will fix that.

## 7.7 Reading the result back

```matlab
R.noise.Cgroups          % one covariance function per pool group
R.noise.model{k}.sigma   % the estimated LEVEL of stream k
R.noise.model{k}.C       % its law; EMPTY where no law was estimated
R.noise.stream           % the stream keys, in the order of model{}
R.history(end).white     % Ljung-Box whiteness per stream (K = 20, alpha = 0.05)
R.convergence(); R.plot("psd"); R.plot("convergence");
```

An arc's standard deviation is **all three terms**: declared `sigma` × sqrt(arc factor) ×
`R.sigma0`. The level and the factor split one scale by a convention (Ellmer eq. 6.5.30's
gamma-bar normalisation) and sigma0 carries the rest, as in any GLS adjustment. Measured on
a fixture planted at 0.5 / 1.5 / 0.5: 0.517 / 1.499 / 0.507. With `sigma = 1` and a `Cep`
carrying the absolute accuracy, that product reads directly as "by how much the
preprocessing's covariance was too optimistic on this arc".

The acceptance of the estimator is physical, and its two halves are pinned differently. The
planted FIELD is recovered within its formal errors (`demo_gravity` reports z-scores per
coefficient; `test_physOracle` asserts them). The planted stream SIGMAS carry no formal
error, and what the suite asserts is a RELATIVE bound: the ratios between streams to 25 %
(measured 1.18 / 2.10 / 0.51 against a planted 1 / 2 / 0.5).

## 7.8 Robustness: what Ellmer does, and what `Robust` adds inside the loop

Three things are worth separating, because Ellmer's own robustness is not per-observation
M-estimation and all three of his are available with `Estimate = "noise"` alone:

1. **Gross outlier screening happens in preprocessing**, before the adjustment (Ellmer
   §6.2, after Klinger 2018: the level-1B data is screened for outliers and periods of
   reduced quality, and segments shorter than one hour are discarded). Nothing in the
   adjustment replaces that.
2. **The arc-wise variance factors ARE the down-weighting** at arc granularity: a bad arc
   receives a large variance factor and loses weight automatically. That is
   `ArcFactors = true`.
3. **The normalisation of those factors is itself robust**: gamma-bar (eq. 6.5.29) is
   computed in GROOPS from the central half of the sorted factors, ignoring both tails.
   That is `Noise.MeanType = "robust"`, the default here.

**Per-observation down-weighting inside the same loop (since 6.1).** `Robust = "huber" |
"igg3" | "bisquare"` together with `Estimate = "noise"` runs the robust pass *inside* the
joint estimator: every pass ends with one residual sweep, its law-whitened residuals (the
current law, level, arc factor and `Cep` divided out — the standardised residual under the
current stochastic model) go through the influence function (`RobustOpts.Tuning` for the
constants), and the resulting factors ride inside every arc's whitener on the next pass.
Accumulation, the probe redundancies and the component update therefore all see **one**
re-weighted problem: a row of weight w contributes w to Omega *and* w to the redundancy, so
the Helmert equations stay consistent — the same thing that keeps the amplitudes and the arc
factors from double counting. `MaxPasses` is the budget (`RobustOpts.MaxIter` / `Tol` belong
to the stand-alone IRLS). Prior streams are never re-weighted: a down-weighted datum is a
rank defect, not an outlier. What comes back: `R.noise.robust` (`influence`, `scale`,
`factor` per observation in sweep order, `keys`), and `R.history.robustScale` / `nDown` /
`nZero` per pass.

**The price, stated.** The levels are those of the *retained* observations: no consistency
correction is applied. Under Huber (c = 1.345) on Gaussian residuals a row of weight
w = min(1, c/|u|) contributes w·u² and w, so the effective level ratio is
E[w u²]/E[w] ≈ 0.80 (a sigma about 11 % low), and a stream with no outliers at all is still
estimated at that ratio. The oracle in `test_vcePolicy` measures exactly this. If you want
the plain GLS level with outliers removed, use `R.noise.robust.factor` to identify them, fix
the data, and re-run without `Robust`. This composition is Koch's iterated re-weighting
inside a variance-component estimation, not something Ellmer (2018) does.

Two further routes remain useful:

- **Diagnose, then decide.** Run `Estimate = "noise"` and read `R.report()` — the w-tests,
  redundancies and MDBs identify the observations carrying the outliers; remove or repair
  them in your `buildDay` and re-run.
- **Fold fixed weights into `Cep`.** Take factors from any robust run and let `buildDay`
  return `Cep` divided by sqrt(factor) on those rows. The weights are then data; the noise
  estimation runs on the re-weighted stream without the retained-observation bias moving
  with every pass. The weights must be stored, never recomputed inside `buildDay` (the day
  builder must be deterministic).
