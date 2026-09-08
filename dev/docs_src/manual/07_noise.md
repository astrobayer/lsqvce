# Stochastic models and their estimation

Every observation stream carries a covariance spec in `data.systems(k).covariance.<type>`
(chapter 3 lists the fields) — its TEMPORAL law: `diag` (white), `arma` / `ar` (a whitening
filter), `acf` / `psd` (an exact Toeplitz law). The per-epoch accuracy of the observations
(a sigma per epoch, or a full d × d covariance per epoch for a multi-axis stream) is not
declared: `buildDay` returns it with the day as `D.Cep` (chapter 3, §7) and the toolbox
composes C = D (R_t ⊗ I_d) D per arc — for every temporal law, scalar streams included.
Batches are overlapped by the filter memory automatically, so batching never changes the
result.

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

The acceptance of the estimator is physical, and the two halves are pinned differently
(R6 WP-101 — the sentence here used to claim one bound for both). The planted FIELD is
recovered within its formal errors: `demo_gravity` reports z-scores per coefficient and
`test_physOracle` asserts them. The planted stream SIGMAS carry no formal error, and what
the suite asserts is a RELATIVE bound: the ratios between streams, to 25 %
(`test_vcePolicy>testRecoversPlantedRelativeWeights`; measured 1.18 / 2.10 / 0.51 against
a planted 1 / 2 / 0.5). A per-arc standard deviation is read back as declared sigma ×
sqrt(arc factor) × `R.sigma0` (chapter 3 §2).
