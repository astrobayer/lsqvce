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
3. updates the spectral shape per pool group in the PSD domain (the loader turns it into
   its exact ACV at every lag of each arc -- since 6.6.2 no longer a resampled, circulant
   one, review R7-01), the level per stream, the per-arc factors, and the sigma
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

The whiteness test pools a stream's arcs: the lag products of every arc are summed, and the
statistic divides each lag by the number of products it actually holds (N_k, the sum over arcs
of m_a − k), so it is sized on short arcs too (since 6.1.8; before, 64 white arcs of 30
epochs were flagged with probability 0.001 instead of 0.05). K stays 20 lags: on arcs of
thousands of epochs the test sees only the shortest lags.

An arc's standard deviation is **all three terms**: declared `sigma` × sqrt(arc factor) ×
`R.sigma0`. Since 6.1.6 an arc factor is the arc's **absolute** Helmert factor (Ellmer eq.
6.5.28) — the same number a stream without arc factors averages into its level — so sigma0
settles near 1 and a stream's level is the mean of its factors. Ellmer normalises the factors
by gamma-bar (eqs. 6.5.29–30) because his PSD carries the level; here the law is a shape, so
the level has to stay in the factors. Two earlier forms were wrong: gamma-bar over every
stream's arcs tied one stream's sigma to *another* stream's truth (review R7-09: 1.01 / 0.79 /
1.15 as only the other stream changed), and the 6.1.5 per-stream form counted the level twice
as soon as a stream's arcs differ (one arc 10× in sigma among nine: 1176 / 11.65 against the
planted 100 / 1). With `sigma = 1` and a `Cep`
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
2. **The arc-wise variance factors** are Ellmer's down-weighting at arc granularity: a bad
   arc receives a large variance factor, and since 6.1.6 that factor **weights the solve**
   (review R7-13): each arc enters the next pass with its own variance, level × factor,
   and the same factor multiplies the update (eq. 6.5.28). Measured on one arc 10× in sigma
   among nine, beside two other streams (16 seeds): the factors come back 100.0 / 1.00, and
   the true scatter of the estimated parameters is 13–47 % smaller than when only the
   stream's mean factor reached the weighting. That is `ArcFactors = true`; the factors are
   reported in `R.noise.model{k}.sigma2`.
3. **The normalisation of those factors**: gamma-bar (eq. 6.5.29) is computed in GROOPS from
   the central half of the sorted factors, ignoring both tails. In this joint estimator the
   factors are absolute (§7.7), so gamma-bar no longer scales them; it sets only the scale of
   the floor below which a factor is clamped (`lsqvce:arcVarianceFactors:clamped`), with the
   robust mean. The option that chose plain or robust, `Noise.MeanType`, is retired since
   6.1.8; the kernel `lsqvce.engine.arcVarianceFactors` keeps both forms of the thesis.

**Per-epoch down-weighting inside the same loop (since 6.1), in the Koch convention.**
`Robust = "huber" | "igg3" | "bisquare"` together with `Estimate = "noise"` runs the robust
pass *inside* the joint estimator. Every pass ends with one residual sweep; the **raw**
residuals, descaled by each epoch's own declared accuracy (`Cep`) and reduced to one value
per epoch, go through the influence function per stream (`RobustOpts.Tuning` for the
constants). For a multi-axis epoch that value is **χ-calibrated** (since 6.1.4): the sum of
squares of its d descaled rows is mapped through the χ²(d) distribution to the equivalent
normal deviate, with the sign of the rows' sum, so a clean 3-axis epoch is standardised
exactly like a clean scalar residual — the RMS used before was χ-shaped and, standardised
uncentred, zeroed about a quarter of clean 3-axis epochs under `igg3`. When a pass zeroes
more than 45 % of the epochs the toolbox warns once that it is at the breakdown point
(`lsqvce:robust:atBreakdown`; diagnostic only — the stand-alone IRLS says the same about
rows). The factors enter the next pass as an **additive nugget** on the temporal law
(since 6.6.0): epoch e's temporal variance Σ_ee becomes Σ_ee / w_e, its covariances with the
neighbours stay the law's, **before** the composition with `Cep`:
C = D ((R_t + diag(q ∘ Σ_ee)) ⊗ I_d) D with q = 1/w − 1. Until 6.5.0 the factors were
*folded* into `Cep` (sigma / sqrt(w)), which scales the neighbours' covariances too: under a
correlated law a zero-weight epoch then stayed a phantom observation the filter leaned on —
at φ = 0.9 the formal standard deviation was 0.31 of the true one, and the estimate 1.2–1.5
times worse than deleting the epoch (review R7, R7-10; measured at tier A, HANDOVER §12.43).
With the nugget a zero weight *is* a deletion (formal/true 1.000). A nugget arc is whitened
in one batch: the law's filter, then an exact rank-k correction over the k down-weighted
epochs (Woodbury), O(m k + k³) — 0.1–4 ms for a 4096-epoch arc. A weight
is a statement about the *observation*; the stand-alone IRLS (`Estimate = "none"`) weights
the law-whitened *innovation* instead (Yang 1994), and the two are deliberately different:
the control measurement showed that an innovation weight of zero neutralises the whitening
filter's whole tail around a blunder, so a spectral estimator loses exactly the lag
correlations it needs there (a planted colour of 0.20 crawled to 0.054 in 16 passes), while
the same blunders declared through `Cep` let the plain estimator recover it. `MaxPasses` is
the budget (`RobustOpts.MaxIter` / `Tol` belong to the stand-alone IRLS). Prior streams are
never re-weighted: a down-weighted datum is a rank defect, not an outlier. What comes back:
`R.noise.robust` (`influence`, `scale` per stream, `factor` per **epoch** in sweep order,
`keys`, `truncation` per stream), and `R.history.robustScale` / `nDown` / `nZero` per pass.

Two things the mechanism learned by measurement, both binding on how to use it. **The
standardisation is the stream's own MAD scale, never the estimated level**: a global scale
mis-standardises every stream whose level differs from the mix, and using the level the
whitener just applied makes the weights feed on the level they shape — the clean stream's
level spiralled to 0.38 that way. With `v / MAD(v)` per stream the weights are invariant to
the level, and the level converges. **Prefer a redescending influence for gross blunders**:
Huber is monotone and keeps c/|u| of a blunder's weight, so 2 % of a stream at 25 sigma
retain 0.02 · 1.345 · 25 ≈ 0.67 of white variance and dilute a lag-1 colour of 0.20 to
about 0.12; `igg3` (zero beyond 3 sigma) and `bisquare` reject them outright and recover
the colour.

**The level (since 6.6.0).** Up to 6.5.0 the factors also weighted the component update,
and the level of a *clean* stream came back **14–19 % low** — measured on the runner over 20
seeds with the true law declared (plain 0.997, `igg3` 0.81 / 0.86; HANDOVER §12.43; the one
earlier reading of 1.007 had been read against 1 while the truth was the process sd 1.155).
The cause is selection: the influence down-weights exactly the epochs with large residuals,
which no redundancy count undoes (counting the redundancy with the weights, exact for fixed
weights, moved a mirror only from 0.74 to 0.76). So the component update now sees only the
**rejected** epochs (factor 0, deleted by their nugget); partial factors weight the solution
only; and each stream's level is divided by the Gaussian truncation factor of the rejection
point, τ = F_{d+2}(q_k) / F_d(q_k) with P(χ²_d > q_k) = P(|N(0,1)| > k) — 0.9733 for
`igg3`'s k1 = 3 at d = 1, practically 1 for `bisquare`, exactly 1 for Huber, which rejects
nothing (`R.noise.robust.truncation`). A mirror of the loop reads 0.978 / 0.976 (φ = 0.5 /
white) before the factor. The returned solution is hard too: after the loop one accumulation
with the rejected epochs deleted and every partial factor set back to 1, whose `sigma0` and
`Cxx` are corrected for the rejected rows and for τ. With the partial factors in that last
solve, σ·σ0 had stayed at 0.83, because `sigma0`, recomputed under the corrected levels,
gave the selection bias back (measured on the runner). Under a correlated law τ is the white-noise value, an approximation
of the same order as that residual. This composition is Koch's iterated re-weighting inside
a variance-component estimation, not something Ellmer (2018) does.

Two further routes remain useful:

- **Diagnose, then decide.** Run `Estimate = "noise"` and read `R.report()` — the w-tests,
  redundancies and MDBs identify the observations carrying the outliers; remove or repair
  them in your `buildDay` and re-run.
- **Declare fixed weights as data.** Take factors from any robust run and let `buildDay`
  remove the rejected epochs (a data gap) or return `Cep` divided by sqrt(factor) for a partial
  one. The weights are then data and stay put from pass to pass. Under a correlated law a
  `Cep` division is the fold described above, so a zero factor belongs in the gap, not in
  `Cep`. The weights must be stored, never recomputed inside `buildDay` (the day
  builder must be deterministic).
