# `lsqvce.Results`

Every run of `lsqvce.fit` returns one `lsqvce.Results`. Its full help block, with every
property and method, is in the chapter *Help of the user-tier units*; this chapter is
the map.

## Properties

| property | size | content |
|---|---|---|
| `x` | nc×1 | the common-parameter estimate |
| `sigma0` | scalar | a-posteriori standard deviation of unit weight |
| `stopReason` | string | `"single"`, `"converged"`, `"maxPasses"` or `"userStop"` |
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
| `R.convergence()` | one row per pass: when the loop settled (iterative runs) |
| `R.coefficients(Alpha=a)` | estimates with formal errors and (1 − a) confidence bounds |
| `R.toTable(data)` | your parameter tables with estimate and standard deviation written in |
| `R.residuals(...)` | whitened (and raw) residuals per stream, re-read from the arcs |
| `R.contribution(...)` | the information share of every stream in the common parameters |
| `R.redundancy()` | the total redundancy n − u |
| `R.reliability()` | per observation: r, w, the minimal detectable bias and its quadrant |
| `R.systematic(...)` | unmodelled signal in the residuals: Ljung–Box, runs, binned means |
| `R.snoop(...)` | iterative data snooping with a clean re-solve; R is not changed |
| `R.compare(Rred)` | the nested-model F test against a reduced fit |
| `R.ellipse(i, j, ...)` | the confidence ellipse of a parameter pair |
| `R.plot(kind)` | `"residuals"`, `"contribution"`, `"psd"`, `"convergence"`: one figure per question |

The readers that need an iteration (`summary`, `convergence`) raise
`lsqvce:Results:noHistory` on a single solve; those that need the accumulator
(`residuals`, `contribution`) raise `lsqvce:Results:noAccumulator` on the in-memory
kinds.
