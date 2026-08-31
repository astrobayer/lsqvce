# lsqvce release notes

## 2.0 — 2026-09-01 (R4 track, Levels L1–L5)

**The one breaking release of the R4 redesign.** Version 1.x users: read the
migration table at the end; everyone else starts at `doc lsqvce` → Getting
started.

### The user surface is four names
`data` struct (DATA_CONTRACT.md) + `lsqvce.Options` + `lsqvce.fit` + `lsqvce.Results` (`.report`).
`fit` covers three kinds of data (streaming with days/arcs/local unknowns,
in-memory matrices, batched callbacks incl. nonlinear) and four models
(`gls`, `vce`, `pooled`, `ellmer`), with `Robust="irls"` on the in-memory kinds.

### BREAKING
1. **Namespace layering (R4-D1, WP-64a).** Engine and kernel functions moved
   into the sub-packages `lsqvce.engine` and `lsqvce.kernel`; the root package
   holds the USER tier only. **No aliases** (toolbox rule 10): old names stop
   working. Error and warning identifiers (`lsqvce:<unit>:<rule>`) are unchanged.
2. **`Results` from the in-memory kinds**: `R.sigma0` is the a-posteriori std
   of unit weight and `R.sol.Cxx` the covariance scaled by it (the streaming
   contract); the driver's a-priori `Sigma0` is kept as `sol.sigma0Prior`.
3. `adjust`'s help describes its `sol` as it is (no `.Cxx`; `sigma0` is the
   a-priori value) — the documentation changed, the numbers did not.

### Fixed
- **G-06** the raw-closure template ignored `adjustVCE`'s per-satellite reweighting.
- **G-07** a system with no local unknowns (nl = 0) is valid.
- **G-08** `jointModelOpt`'s covariance lookup assumed unit sampling; for
  dt ≠ 1 s every arc was skipped as non-PD.
- The batch→`Results` bridge establishes the streaming `Results` contract.

### Added
`lsqvce.fit`, `lsqvce.Options` (59 typed knobs, engine defaults marked in
`disp`), `Results.coefficients()`, `DATA_CONTRACT.md`, the Getting-started page
and examples gallery, `lsqvce_theory.pdf` (the ONE theory authority, R4-D3 —
the manual's theory chapter now points there), `DistributeDays` /
`ParallelThreshold` through `fit`, audit check 18 (tier census + Windows case
guard), `buildfile.m` (`buildtool audit`, `test`, `package`).

### Still fluent (R4-D2, deferred to Level L6)
`lsqvce.engine.adjustVCE` keeps its configurator methods for now; `fit` is the
supported way to reach the loop. Their deletion (constructor Name-Values
`Shape` / `Scale` / `Outer` / `Joint` / `Axes` instead) is the next level.

### Migration table (1.x → 2.0)
| 1.x | 2.0 |
|---|---|
| the direct driver calls (`lsqvce.adjust(u, batchFcn, n)`, …) | `R = lsqvce.fit(data, opts)` (recommended), or the qualified engine name below |
| `lsqvce.Accumulator` | `lsqvce.engine.Accumulator` |
| `lsqvce.BatchIndex` | `lsqvce.engine.BatchIndex` |
| `lsqvce.BlockAccumulator` | `lsqvce.engine.BlockAccumulator` |
| `lsqvce.Combine` | `lsqvce.engine.Combine` |
| `lsqvce.CommonPrior` | `lsqvce.engine.CommonPrior` |
| `lsqvce.NormalEqSystem` | `lsqvce.engine.NormalEqSystem` |
| `lsqvce.PrefetchLoader` | `lsqvce.engine.PrefetchLoader` |
| `lsqvce.QRSystem` | `lsqvce.engine.QRSystem` |
| `lsqvce.ReducedContribution` | `lsqvce.engine.ReducedContribution` |
| `lsqvce.adjust` | `lsqvce.engine.adjust` |
| `lsqvce.adjustBlock` | `lsqvce.engine.adjustBlock` |
| `lsqvce.adjustVCE` | `lsqvce.engine.adjustVCE` |
| `lsqvce.adjustVCEjoint` | `lsqvce.engine.adjustVCEjoint` |
| `lsqvce.applyRegularisation` | `lsqvce.engine.applyRegularisation` |
| `lsqvce.ar2arma` | `lsqvce.kernel.ar2arma` |
| `lsqvce.arcVarianceFactors` | `lsqvce.engine.arcVarianceFactors` |
| `lsqvce.armaACV` | `lsqvce.kernel.armaACV` |
| `lsqvce.armaWhiten` | `lsqvce.kernel.armaWhiten` |
| `lsqvce.checkPSD` | `lsqvce.engine.checkPSD` |
| `lsqvce.collectResiduals` | `lsqvce.engine.collectResiduals` |
| `lsqvce.contribution` | `lsqvce.engine.contribution` |
| `lsqvce.covComponent` | `lsqvce.engine.covComponent` |
| `lsqvce.covModel` | `lsqvce.engine.covModel` |
| `lsqvce.crossCorr` | `lsqvce.kernel.crossCorr` |
| `lsqvce.datum` | `lsqvce.engine.datum` |
| `lsqvce.dctI` | `lsqvce.kernel.dctI` |
| `lsqvce.decorrFactor` | `lsqvce.kernel.decorrFactor` |
| `lsqvce.diagnose` | `lsqvce.engine.diagnose` |
| `lsqvce.diagnoseAccum` | `lsqvce.engine.diagnoseAccum` |
| `lsqvce.ellipse` | `lsqvce.engine.ellipse` |
| `lsqvce.ftest` | `lsqvce.engine.ftest` |
| `lsqvce.gcvLambda` | `lsqvce.engine.gcvLambda` |
| `lsqvce.iterateNoiseModel` | `lsqvce.engine.iterateNoiseModel` |
| `lsqvce.jointModelOpt` | `lsqvce.engine.jointModelOpt` |
| `lsqvce.leverage` | `lsqvce.engine.leverage` |
| `lsqvce.leverageBlocks` | `lsqvce.engine.leverageBlocks` |
| `lsqvce.levinsonDurbin` | `lsqvce.kernel.levinsonDurbin` |
| `lsqvce.localDiagnostics` | `lsqvce.engine.localDiagnostics` |
| `lsqvce.makeLoader` | `lsqvce.engine.makeLoader` |
| `lsqvce.mapreduce` | `lsqvce.engine.mapreduce` |
| `lsqvce.normalizeSatList` | `lsqvce.engine.normalizeSatList` |
| `lsqvce.outliers` | `lsqvce.engine.outliers` |
| `lsqvce.poolGroups` | `lsqvce.engine.poolGroups` |
| `lsqvce.psd2acv` | `lsqvce.kernel.psd2acv` |
| `lsqvce.psdVCE` | `lsqvce.engine.psdVCE` |
| `lsqvce.redundancy` | `lsqvce.engine.redundancy` |
| `lsqvce.redundancyLags` | `lsqvce.kernel.redundancyLags` |
| `lsqvce.residualAcv` | `lsqvce.kernel.residualAcv` |
| `lsqvce.residuals` | `lsqvce.engine.residuals` |
| `lsqvce.resolveNc` | `lsqvce.engine.resolveNc` |
| `lsqvce.resolveSatSpecs` | `lsqvce.engine.resolveSatSpecs` |
| `lsqvce.robust` | `lsqvce.engine.robust` |
| `lsqvce.satKeyOf` | `lsqvce.engine.satKeyOf` |
| `lsqvce.selectArAcv` | `lsqvce.engine.selectArAcv` |
| `lsqvce.selectArOrder` | `lsqvce.engine.selectArOrder` |
| `lsqvce.separableWhiten` | `lsqvce.kernel.separableWhiten` |
| `lsqvce.sharedProbe` | `lsqvce.engine.sharedProbe` |
| `lsqvce.snoop` | `lsqvce.engine.snoop` |
| `lsqvce.solveNonlinear` | `lsqvce.engine.solveNonlinear` |
| `lsqvce.stat` | `lsqvce.engine.stat` |
| `lsqvce.streamContribution` | `lsqvce.engine.streamContribution` |
| `lsqvce.synthLoader` | `lsqvce.engine.synthLoader` |
| `lsqvce.systematic` | `lsqvce.engine.systematic` |
| `lsqvce.toeplitzCofactor` | `lsqvce.kernel.toeplitzCofactor` |
| `lsqvce.traceEstimate` | `lsqvce.engine.traceEstimate` |
| `lsqvce.varACV` | `lsqvce.kernel.varACV` |
| `lsqvce.varComp` | `lsqvce.engine.varComp` |
| `lsqvce.varWhiten` | `lsqvce.kernel.varWhiten` |
| `lsqvce.vce` | `lsqvce.engine.vce` |
| `lsqvce.vceHutchinson` | `lsqvce.engine.vceHutchinson` |
| `lsqvce.whitenessCorr` | `lsqvce.kernel.whitenessCorr` |
| `lsqvce.whittle` | `lsqvce.engine.whittle` |

Created 2026-09-01 (CEST, UTC+02:00) with assistance from Claude (Anthropic), model Claude Fable 5 — R4 Level L5, WP-65.
