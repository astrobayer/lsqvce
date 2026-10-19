# `lsqvce.Results`

Every run of `lsqvce.fit` returns one `lsqvce.Results`. Its full help block, with every
property and method, is in the chapter *Help of the user-tier units*; this chapter is
the map.

## Properties

| property | size | content |
|---|---|---|
| `x` | nc×1 | the common-parameter estimate |
| `sigma0` | scalar | a-posteriori standard deviation of unit weight |
| `stopReason` | string | `"single"`, `"converged"`, `"maxPasses"`, `"userStop"` or `"diverged"` (Gauss–Newton found no downhill step) |
| `converged` | logical | whether the iteration converged (true for a single solve) |
| `sol` | struct | the solve: `.x`, `.sigma0`, `.Cxx` (covariance scaled by σ₀²), … |
| `diag` | struct | per-observation diagnostics: residuals in the whitened frame, redundancy numbers, test statistics, minimal detectable biases, the robust check |
| `history` | 1×I struct | one record per pass of an iterative run (empty for a single solve) |
| `noise` | struct | the estimated stochastic model: per stream `.stream` keys and `.model{k}` (`.sigma` the level, `.C` the law's autocovariance); `.S` / `.Sgroups` the bin spectra per pool group; `.priors`, `.priorHistory`; `.robust` (influence, scale per stream, one factor per epoch) |
| `acc` | `BlockAccumulator` | the accumulator: re-solve, inspect blocks, residuals |
| `systems` | 1×M struct | what ran per system: id, nl, day list, name |
| `meta` | struct | provenance: `.toolboxVersion`, `.matlabRelease`, `.createdUTC`, `.seeds`, `.options` (the resolved options), `.warnings` |

## Methods

| method | returns |
|---|---|
| `R.report()` | the verdict and its evidence, levels 0 to 4 (`lsqvce.report`) |
| `R.summary()` | one row per stream: final σ, law, whiteness (iterative runs) |
| `R.convergence()` | one row per pass: the stream levels, whiteness, `maxDev`, and the stop-rule quantities `binRatio` and `levelDev`, with the solution's step `xStep` (iterative runs) |
| `R.coefficients(Alpha=a)` | estimates with formal errors and (1 − a) confidence bounds |
| `R.toTable(data)` | your parameter tables with estimate and standard deviation written in |
| `R.residuals(...)` | whitened (and raw) residuals per stream, re-read from the arcs |
| `R.contributions()` | every (system, day)'s reduced contribution under the run's final weights, for `lsqvce.solveContributions` (one more pass over the data) |
| `R.contribution(...)` | the share of every group in the common solution, t_g = tr(N⁻¹N_g) (Level 2: per parameter) |
| `R.redundancy()` | the total redundancy n − u |
| `R.reliability()` | per observation: r, w, the minimal detectable bias and its quadrant |
| `R.systematic(...)` | unmodelled signal in the residuals: Ljung–Box, runs, binned means |
| `R.snoop(...)` | iterative data snooping with a clean re-solve; R is not changed |
| `R.compare(Rred)` | the nested-model F test against a reduced fit, under one fixed weighting: fits with estimated weights (`Estimate = "noise"`, `Robust`) are refused (`lsqvce:Results:compareEstimatedWeights`) |
| `R.ellipse(i, j, ...)` | the confidence ellipse of a parameter pair |
| `R.plot(kind)` | `"residuals"`, `"contribution"`, `"psd"`, `"convergence"`: one figure per question |

The readers that need an iteration (`summary`, `convergence`) raise
`lsqvce:Results:noHistory` on a single solve; those that need the accumulator
(`residuals`, `contribution`) raise `lsqvce:Results:noAccumulator` on the in-memory
kinds.

## `lsqvce.solveContributions` — a bootstrap over days

`C = R.contributions()` streams every day once more and eliminates each day's own local
parameters, so that the day's reduced normal equations remain. `out =
lsqvce.solveContributions(C, c)` solves any weighted combination of the days:

    N* = sum_d c_d N_d,   n* = sum_d c_d n_d,   l'l* = sum_d c_d l'l_d,   dof* = sum_d c_d dof_d

`c` holds one multiplicity per day of `C.days`: how often a replicate drew the day (0, 1,
2, ...), or any non-negative real for a Bayesian bootstrap. On the QR backend the days'
square-root rows are stacked (times √c_d); the common prior is added once, and the solve is
fit's own. One replicate costs O(D·nc²), independent of how many observations a day held.

| output | meaning |
|---|---|
| `out.x`, `out.Cxx`, `out.sigma0`, `out.s0sq`, `out.dof`, `out.vtv` | the combination's solution |
| `out.weights` | `.estimate`, `.robust`, `.fixed = true`: the weights every day carries |
| `out.locals` | with `'Locals', true`: per drawn day `.sysID`, `.day`, `.multiplicity`, `.localIdx`, `.y`, `.Cyy` |

The weights are the run's final ones, fixed: a noise run's estimated variance components
and a robust run's hard weights are not re-estimated per combination (a bootstrap of the
estimator under fixed weights); a robust run's σ0 is divided by its truncation factor, as
fit's is. Contributions exist only where the days are separable:

| identifier | cause |
|---|---|
| `lsqvce:accum:localsCoupleDays` | a local parameter takes observations on two days (declare it per day or per arc) |
| `lsqvce:accum:regularisedLocals` | the accumulation regularised a local block |
| `lsqvce:accum:dayLocalsRankDeficient` | a day leaves its own locals undetermined |
| `lsqvce:Results:contributionsNonlinear` | a Gauss–Newton run (its rows are linearised about x) |
| `lsqvce:Results:contributionsInconsistent` | the days' dof do not add up to the run's |
| `lsqvce:solveContributions:badMultiplicity` | `c` does not fit `C.days`, or is negative, or draws nothing |

