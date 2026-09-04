# Stochastic models and their estimation

Every observation stream carries a covariance spec in `data.systems(k).covariance.<type>`
(chapter 3 lists the fields): `diag` (weights), `arma` / `ar` (a whitening filter),
`acf` / `psd` (an exact Toeplitz law), `separable` (a temporal law times a fixed per-epoch
covariance — the declaration for multi-axis streams; a `diag` law would ignore the
per-epoch covariance). Batches are overlapped by the filter memory automatically, so
batching never changes the result.

With `Estimate="noise"` the joint estimator (theory ch. 6, Ellmer 2018) runs on every
stream whose spec allows it (`estimate = "law"` default; `"scale"` = level only, the
declared shape kept; `"none"` = untouched). Per pass it accumulates and solves under the
current laws, sweeps every arc once for its whitened residual and Hutchinson-probed
redundancies at every lag (all lags from one FFT per column), updates the spectral shape
per pool group in the PSD domain (positive by construction, resampled to every arc's
length), the level per stream, the per-arc factors, and the prior groups' σ; it stops at
`Tolerance` or `MaxPasses`. A law with a per-row sigma is a fixed structure: it can carry
a level but not a spectral shape.

What comes back: `R.noise.Cgroups` (one covariance function per pool group),
`R.noise.stream` / `R.noise.model{k}` (per stream: the level, its law, its arc factors),
`R.history.streams / sigma / white / whitePval` per pass (the Ljung–Box whiteness of every
stream, K = 20, α = 0.05), `R.convergence()`, `R.plot("psd")`, `R.plot("convergence")`.

The acceptance of the estimator is physical: a planted field and planted stream sigmas
are recovered within their formal errors (`demo_gravity`, counted in the test suite).
