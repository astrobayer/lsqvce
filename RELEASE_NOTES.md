# lsqvce release notes

## 6.10.0 — 2026-10-11 (WP-114: priors without a second path)

- **A datum's covariance no longer depends on the datum's σ** (review R7, R7-05). x was
  already the minimum-constraint solution for any weight. The covariance, however, was
  s0²·inv(N + BB'/σ²), so it carried σ² on the datum directions: se 0.005 at σ = 1 and
  2.48 at σ = 1e3 on the levelling network.
  - `BlockAccumulator` now recognises a datum: a prior whose directions the data do not
    see, |N K'| ≤ 1e-10 |N||K|.
  - For a datum, `sol.Cxx` and the parameter cofactor of the diagnostics (`D.sigX`,
    `D.Q.block`, correlations) are the minimum-constraint ones, s0²·pinv(N).
  - They are computed stably: the datum rows are re-weighted to the data's scale, then
    inv(Rm'Rm) − Dm·Dm'. Subtracting the σ² term from the solve's cofactor would lose
    11 digits at σ = 1e6.
  - Leverage keeps the solve's cofactor, so the datum rows stay absorbed. Data rows are
    unaffected.
  - A prior only partly in the null space is genuine information and keeps the solve's
    covariance.
- **`residuals()` and `redundancy()` count a common prior** (R7-15).
  - `residuals()` adds the prior's residual to vtv and its rows to dof. New field:
    `res.prior`.
  - `redundancy()` closes the identity as sum(r_i) + r_prior = m + nk − u. New field:
    `red.priorRedundancy`.
  - Before, every fit with `data.priors` warned `lsqvce:residuals:sigma0Mismatch`, which
    blamed the loader, and `lsqvce:redundancy:identityViolated`. The wrong dof was marked
    authoritative.

## 6.9.0 — 2026-10-11 (WP-113: Gauss–Newton damps the parameters and keeps the prior on x)

- **The damping now acts on the common block** (review R7, R7-14). `adjustNonlinear` used to
  damp through the accumulator's `Reg "lm"`. That option acts on a local block, and a batch
  source has none, so the damping did nothing:
  - a rejected step was re-solved unchanged until `MaxIter`;
  - `R.x` was then the last accepted point plus that rejected step (F 36.8 against the best
    accepted 4.9 on the review's exponential fit).

  What changed:
  - The step is now one solve with Marquardt rows `sqrt(λ·d_j)·dx_j ≈ 0`, d = diag(J'PJ).
  - A rejection re-solves on the same accumulation with ten times the damping.
  - The descent check runs from the first pass, on the declared objective with the prior
    included.
  - An unconverged run returns the best accepted point.
- **A prior on a nonlinear problem acts on x, not on the step** (R7-04). Each step now
  solves `K·dx ≈ x0 − K·x`. The old form `K·dx ≈ x0` converged away from the declared
  optimum: b = +0.60 against +0.30 with σ = 1e-6, and the unconstrained answer when x0 = 0.
- `Lambda0 = 0` now means "start undamped": a rejected step still turns the damping on,
  starting at 1e-3.
- Mirror (the R7 review's mirrors, `wp113`):
  - the exponential fit reaches the optimum from [10, −5] and from [0.1, 1], with 30 passes;
  - Rosenbrock converges from five starts, the classic one included;
  - the prior fixtures land on the declared optimum to 1e-7.
- `R.history(k).vPv` is now F(x) with the prior term included.

## 6.8.0 — 2026-10-11 (WP-112: the QR backend's accumulator path is square-root throughout)

- **The default backend no longer squares the condition number on the accumulator path**
  (review R7, R7-03). `BlockAccumulator` used to read each system's QR factor back as
  N = R1'R1, Schur-eliminate the locals in normal form, and re-root the contribution with
  `chol` in `merge`. It did so also for systems without locals, and the prior-folded solve
  did the same.
  - Now the locals are eliminated on the factor: `QRSystem.reducedRows` re-triangularises
    with the local columns first.
  - The reduced contribution is merged as square-root rows (`ReducedContribution.Rows`).
  - A common prior is stacked as rows (K/σ, x0/σ).
  - Measured in a mirror at cond(A) = 1e7–1e9, with and without locals: x_c is 1e-8–1e-4
    formal sd from the dense `[A B] \ l`, against a median of 2–3 sd before. No full-rank
    problem is refused; before, 4–11 of 20 were refused at 1e8.5–1e9.
- **The rank decision is the code's own** (R7-16): `QRSystem.RankTol` = 1e3·eps on the
  pivots of the final factor.
  - Every free levelling network is refused (`lsqvce:qr:singularR`), and every one with a
    datum is solved. Before, `chol` decided by rounding: 61 % refused, 39 % solved with an
    arbitrary datum and sd ~1e7.
  - Problems with cond(A) above about 4e12 are now refused.
- At the WP-49/50 deep regime the QR path now returns an honest sigma0 at the fold (no
  `lsqvce:qr:cancellation`): the residual row rides along with the merged rows. A
  regularised local block still folds in normal form, with the WP-50 guard.
- Numbers on the QR backend change in the last digits for well-conditioned problems, and
  substantially for ill-conditioned ones. The normal backend is unchanged.

## 6.7.0 — 2026-10-11 (WP-111 part 3: the joint law as positive components)

- **The joint estimator's law is a non-negative spectrum on `Noise.Grid` bins** (review
  R7, R7-02 and R7-20; user decision 2026-10-11). The spectrum is linear between the bins
  and its covariance at every lag of every arc is the exact integral r = K·S, with K from
  the new kernel `lsqvce.kernel.hatACV`; its projection form `hatACV(L, m, X)` gives
  K'·X by one FFT, without forming K. Each bin is a variance component whose matrix is
  the covariance of a "hat" spectrum, which is positive semi-definite. Its Helmert factor
  (the projections of the arc's lag sums over **all** its lags, summed per pool group) is
  therefore never negative, and the law needs no floor and no zero-frequency rule.
  - Ellmer's components were the columns of the inverse type-I DCT of a truncated ACF.
    They are indefinite, so a bin's factor could go negative, and a floored bin left the
    declared covariance near-singular.
  - Measured on the runner (HANDOVER §12.44), on AR(1) streams at φ 0.8 / 0.9:
    - the shipped estimator's level ran away to 5e7 × the truth, with a negative lag-1
      correlation and coefficient z-scores up to 9.4;
    - at φ 0.5 it never converged;
    - on the same data the positive-component update converged 18/18, with the level
      0.91–1.07 × the truth and |z| ≤ 2.33.
- **Ellmer's white start is regular now:** S0 = ones (white noise) is the default.
  - Option `'C0'` is replaced by `'S0'` (the starting spectrum), with error
    `lsqvce:adjustVCEjoint:badS0`.
  - `'Norm'` is retired from `adjustVCEjoint` and `arcVarianceFactors`.
  - The warning `lsqvce:psdVCE:degenerateStart` is retired: its fixed point belonged to the
    DCT weighting.
- **`psdVCE` is rewritten** as `psdVCE(Q, T, Sprior)` on per-bin statistics.
  - `'Gamma'`, `'Smooth'` and `'FreezeRedundancy'` stay. A bin now freezes on its
    redundancy, `Sprior_j · sum T`.
  - `'Norm'`, `'MinAmplitude'`, `'ZeroFreqRule'`, `'AllowDegenerate'`, `.Cnew` and
    `.degenerate` are gone.
- **The arc factor is the Helmert factor of the arc's whole covariance** (r'Ω / r's over
  all lags). It used to be the DCT-summed, C-weighted lag vector truncated at Grid.
- **Each arc is weighted by its own law.** The bin statistics of an arc are scaled by the
  arc's own variance r(0), read from the covariance it was weighted with. The one pooled C
  that weighted every arc (R6 WP-100) is gone, because K does not depend on the law.
- **The law is reported as the bin spectrum.**
  - `out.S` / `out.Sgroups` and `R.noise.S` / `R.noise.Sgroups`, next to its covariance
    `C` / `Cgroups` at lags 0..Grid-1, which is derived and kept for the record.
  - `out.Shat`, `R.noise.Shat` and `history.Shat` are renamed `S`.
  - `R.plot("psd")` draws the bins.
- **`jointModelOpt` takes `model.S` / `model.Sgroups`** and emits the 'psd' spec with
  ordinates 2·S on fNorm, scaled to r(0) = 1, so covModel reproduces K·S exactly.
  - Error `badC` is now `badS`.
  - The cosine transform of a truncated ACF, floored at eps·max, is gone.
- `Noise.Grid` is the spectral resolution. It no longer has to cover the longest arc:
  the law is exact at every lag whatever the Grid.
- `redundancyLags`: the general (non-Toeplitz) Σ⁻¹ lag sums run in BLAS-3 column blocks
  (memory N × 256), because the joint sweep now asks for every lag.
- **Every estimated law changes.** Runs with `Estimate = "noise"` give different numbers
  from 6.6.2.

## 6.6.2 — 2026-10-11 (WP-111 part 2: the PSD law's covariance without a lag period)

- **A declared PSD law no longer becomes a circulant covariance** (review R7, R7-01).
  `lsqvce.kernel.psd2acv` integrates the piecewise-linear spectrum through the given
  ordinates EXACTLY (a closed form per segment) instead of by trapezoids, and `covModel`
  'psd' 'toeplitz' no longer resamples the law onto floor(m/2)+1 points. Trapezoidal
  quadrature on nf points is periodic in the lag with period 2(nf-1): with the joint
  estimator's 60-bin law, epoch 1 of every arc longer than 118 epochs was correlated with
  epoch m like a neighbour (r(m-1) = r(1) on even arcs, = r(0) on odd ones, Toeplitz cond
  3.3e8 at φ 0.5 -- measured on the runner, HANDOVER §12.44). Now r(m-1)/r(0) ~ 1e-8 and
  cond 9 on every arc length; r(0) = trapz(f, S) as before. Numbers change for every
  'psd' law whose grid is coarse against the arc (the joint estimator's laws above all);
  fine grids (thousands of points) move by the quadrature error only.

## 6.6.1 — 2026-10-11 (WP-115 part 4: one stream key, one arc key)

- **The stream key and the arc key each have one definition** (review R7, R7-43; coding
  rule R1). `lsqvce.util.streamKey(sysID, obsType)` returns `'satKey|obsType'` (and the
  satKey); `lsqvce.util.arcKey(sysID, obsType, day, arc)` returns
  `'s<sys>_<obsType>_d<day>_a<arc>'`. They were composed inline in eight places and the
  stream was recovered from the arc key by a regex in three. The joint estimator now
  carries each arc's stream next to it: new output `out.arcStreams`, aligned to
  `out.arcKeys`. Every key keeps its format; no number changes.
- `lsqvce.engine.satKeyOf` is retired: nothing in the toolbox called it, and its `'s_'`
  prefix was a third convention for the same name.

## 6.6.0 — 2026-10-11 (WP-115 part 3b: the robust pass as an additive nugget; an unbiased robust level)

- **The joint estimator's robust factors are an additive NUGGET** (review R7, R7-10; user
  decision (b)). A factor w of an epoch raises its temporal variance to Σ_ee / w and leaves
  its covariances with the neighbours to the law (`covModel`'s new `cm.nugget`,
  q = 1/w − 1, composed before `Cep`). Until 6.5.0 the factor was folded into `Cep`
  (sigma / sqrt(w)), which scaled those covariances too: under a correlated law a
  zero-weight epoch stayed a phantom observation, and at φ = 0.9 the formal standard
  deviation was 0.31 of the true one (measured on the runner, HANDOVER §12.43). With the
  nugget, formal/true is 1.000 and a zero weight is a deletion. A nugget arc is whitened in
  ONE batch (the law's filter plus an exact Woodbury correction over the down-weighted
  epochs, 0.1–4 ms per 4096-epoch arc); `kernel.redundancyLags` gains `'Nugget'` /
  `'LawW'` to keep its structured Σ⁻¹ sums on such an arc.
- **The level of a robust joint run is no longer biased low** (R7-12). Measured over 20
  seeds with the true law declared: plain 0.997, IGG3 0.81 / 0.86. The influence selects
  the epochs with large residuals; counting the redundancy with the weights does not undo
  that (a mirror moved 0.74 → 0.76). Now (user decision) the component update sees only
  the REJECTED epochs (deleted), partial factors weight the solution only, and each
  stream's level is divided by the Gaussian truncation factor of the rejection point
  (0.9733 for IGG3's k1 = 3; 1 for Huber), returned as `R.noise.robust.truncation`.
  The returned solution is hard as well (rejected epochs deleted, partial factors 1), its
  `sigma0` / `Cxx` corrected for the rejected rows and τ: with the partial factors in the
  final solve, σ·σ0 stayed 0.83 on the runner even with the corrected levels.
  Every robust `Estimate = "noise"` fit changes its numbers.
- `kernel.robustWeights` returns the influence's rejection point as a fifth output.
- New identifiers: `lsqvce:covModel:badNugget`, `lsqvce:covModel:nuggetNeedsWholeArc`,
  `lsqvce:redundancyLags:badNugget`.

## 6.5.0 — 2026-10-11 (WP-86 part 4: one verdict unit; the verdict's methods on R)

- **`lsqvce.engine.verdict` is the ONE engine unit that judges a fit** (the R5 M4
  "verdict from the sweep" row; user decision: the full absorption). `diagnoseAccum`,
  `diagnose`, `snoop`, `systematic`, `ftest` and `ellipse` are its static methods —
  `verdict.source`, `verdict.diagnose`, `verdict.snoop`, `verdict.systematic`,
  `verdict.ftest`, `verdict.ellipse` — with their bodies unchanged; the new
  `verdict.fromAccumulator(acc, sol, ...)` is the source + diagnose pair every driver
  ran (fit, adjustBlock, adjustRobust, adjustNonlinear). The distributions moved to the
  KERNEL tier: **`lsqvce.kernel.stat`** (was `lsqvce.engine.stat`). Breaking for code
  that called the six engine units or `engine.stat` directly: prefix `verdict.`
  (`lsqvce.engine.diagnose(...)` → `lsqvce.engine.verdict.diagnose(...)`,
  `diagnoseAccum` → `verdict.source`) and `engine.stat` → `kernel.stat`. Error
  identifiers are unchanged.
- **Everything the verdict tells you to do is a method on `R`.** New `Results` methods:
  `R.systematic(...)` (Ljung–Box / runs / binned means on the residuals),
  `[D, REJ, sol] = R.snoop(...)` (iterative snooping; R unchanged), `T = R.reliability()`
  (per observation: v, r, w, MDB and the (w, r) quadrant; `lsqvce:Results:noPerObservation`
  without the per-observation diagnostics), `E = R.ellipse(i, j)`. The report's remedies
  and snoop's runtime messages name `R.systematic()` / `R.snoop()`, never an engine
  function (finding H-05).
- `R.diag.store` records what the diagnostics kept per observation.

## 6.4.0 — 2026-10-11 (WP-86 part 3: one block formula, one commons factor)

- **The VCE redundancy pass runs on `lsqvce.engine.BlockCofactor`.** `engine.redundancy`
  (behind `engine.vce`, `outliers`, `contribution`) built, per system, its own copy of the
  commons factor and the nc × nl products (`leverageBlocks`) -- nSys copies of an
  nc × nc matrix, about 0.7 GB each at a GRAVFIRE month -- and evaluated the same
  block-arrow leverage formula a second time (`leverage`). It now builds ONE cofactor
  from the accumulator (`BlockCofactor.fromAccumulator`, the same one the diagnostics use)
  and broadcasts it to every arc task. **`lsqvce.engine.leverage` and
  `lsqvce.engine.leverageBlocks` are retired** with their five identifiers
  (`lsqvce:leverage:*`, `lsqvce:leverageBlocks:*`); `BlockCofactor.leverage(s, Ac, B)`
  is the per-system form (`h` of system s's rows from their commons and local parts).
- **`snoop`'s clean re-solve keeps the block structure.** With `Reaccumulate=true` the
  surviving rows used to be re-solved with every parameter as a common one (a dense
  u_full factor). Each original system is now rebuilt from its surviving batches with its
  own locals (`engine.batchSystem` gained `'LocalColumns'` / `'SystemId'`), the batches
  that touch no locals as one local-free system: the original structure and column
  layout, and the returned `sol.Q` is block-arrow again.
- `BlockCofactor` gained `systemOf` and the static `fromAccumulator`; a local-free system
  keeps its index in `Q.sys`.

## 6.3.0 — 2026-10-10 (WP-86 part 2: the diagnostics' cofactor is kept as blocks)

- **No dense factor of the full parameter set (R7-24).** The diagnostics used a dense
  u × u factor over commons plus all locals (13 729 at a GRAVFIRE month, 1.5 GB per copy;
  `diagnose` held four: R, its inverse, Qxx and the correlation matrix). The new engine unit
  `lsqvce.engine.BlockCofactor` answers every question they ask of Qxx — leverage, Qxx·a,
  diag(Qxx), a 2 × 2 block, the strongest correlation, a conditioning bound, snoop's
  one-row removal — from the reduced commons factor and each system's local blocks, the
  ones the solve used: about 1 GB at that scale. Verified against the dense inverse to
  ~1e-16 (the Python mirror wp86/blockcofactor.py in the R7 review folder) and pinned by
  seven new tests.
- **`diagnoseAccum` reads no data itself.** One batch is one arc; the index comes from
  the loaders' memoised meta answers and the row count from the accumulator. A VCE run with
  `CacheDays = 1` builds each day 2n+3 times (2n+4 in 6.2.0).
- **`snoop` downdates the blocks.** The rank-one removal is a rank-one update of the
  commons factor plus the observation's own local block, and the estimate follows in
  closed form (x + p·v·Qxx~·a): no row-search pass, no right-hand-side re-accumulation pass
  per rejection. The clean re-solve reads one batch at a time (it concatenated every
  surviving row at full width), but still treats all parameters as commons.
- **Changed result fields (breaking for code that read them):** `D.R`, `D.Qxx` and
  `D.corr` are gone, and `R.source.R` with them; the cofactor is `D.Q` / `R.source.Q`
  (`D.Q.block(idx)` is `Qxx(idx, idx)`, `D.Q.diagonal()` its diagonal). `D.condN` /
  `D.condR` are now a LOWER bound (`D.condIsBound`), formed without an SVD; the report says
  so. Retired identifiers: `lsqvce:diagnose:largeU` (the dense matrices it warned about are
  gone), `lsqvce:snoop:badBatch`, `lsqvce:snoop:badIndex` (their helpers are gone). New:
  `lsqvce:BlockCofactor:badBlocks`, `lsqvce:BlockCofactor:multiSystemRow`.
- **What it does not change:** the per-observation cost stays O(nc²) — the commons
  dominate (about 1.7× fewer flops than the dense factor at GRAVFIRE scale).

## 6.2.0 — 2026-10-10 (WP-86 part 1: the diagnostics no longer store the data)

- **Every fit's diagnostic source streams (R7-24).** `lsqvce.engine.diagnoseAccum` stored
  every whitened row at full width — n × (commons + all locals) doubles, kept alive in
  `R.batchFcn` — about 228 GB for one GRAVFIRE-scale month. It now folds each batch into
  the full factor during one enumeration pass, keeps only the batch's address (system,
  day, obsType, arc, sub-batch), and rebuilds a batch on demand through the system's
  loader and the same whitening iterator, caching one arc. The arithmetic is unchanged
  (same batches, same order); the price is one re-read of the data per diagnostic pass,
  and the rebuild relies on the loader contract that an arc is returned identically on
  every call.
- **What remains (WP-86 part 2):** the diagnostics still build a dense factor over the
  FULL parameter set (≈ 1.5 GB at 13 729 parameters) and cost ~u²/2 flops per row and
  pass (~2·10¹⁴ for a GRAVFIRE month). The README says so.

## 6.1.9 — 2026-10-10 (R7 WP-118 part 3: the memory claim corrected; long code on the record)

No behaviour change.

- **The README no longer says "nothing of size n is ever formed" (R7-24).** The solve
  streams; the diagnostics stage does not yet: every fit rebuilds the whitened rows at full
  width, n × (commons + all locals), and keeps them in the `Results` — about 228 GB for one
  GRAVFIRE-scale month, which runs out of memory after the solve has succeeded. Streaming
  that stage is WP-86 (open). `diagnose`'s help and `dev/DEVELOPMENT.md` say the same.
- **Long functions and long lines are on the record (R7-44).** The five functions over 400
  lines (`diagnose`, `report`, `adjustVCEjoint`, `snoop`, `makeLoader`) carry a written
  reason and a code-line ceiling in `dev/docs_src/long_functions.csv`, and the 74 lines over 120
  characters a per-file ceiling in `long_lines.csv`; audit check 33 refuses growth, a new
  unrecorded case, and a stale row. `makeLoader`'s 569-character struct literal is wrapped.

## 6.1.8 — 2026-10-10 (R7 WP-118 part 2: four reported numbers corrected; MeanType retired)

- **`R.compare` counts every parameter (R7-19).** The nested-model F test took the
  redundancy from the common parameters alone; on a streaming fit the eliminated locals are
  parameters too, so the full model's s0² was too small and F too large: 7.4 % false
  positives at α = 0.05. It now uses `R.diag.u` (commons + locals). Kinds B/C (no locals)
  are unchanged.
- **The whiteness test is sized on short arcs (R7-38).** The per-stream Ljung–Box
  statistic of the joint estimator (`R.history.whitePval`, `white`) pooled the arcs' lag
  products but divided by (N − k) as for one long series. It now divides by the number of
  products each lag actually holds. Measured size under a white null: 64 arcs × 30 epochs
  0.001 → 0.048; 200 × 12 0.000 → 0.051; long arcs and single arcs unchanged.
- **Kind B/C warnings reach `R.meta.warnings` (R7-23).** The run's warning inventory was
  reset after the batch path had already warned (`lsqvce:fit:priorNotActive`), so the id
  was lost — and after a run that fired it, the next one warned nothing at all. `fit` now
  resets once, first.
- **Worker warnings reach the run inventory (R7-41).** With `Parallel.Mode = "parfor"` or
  `"spmd"`, the log sink never reached the workers: their warnings were missing from
  `R.meta.warnings` and a once-per-run warning printed once per worker. The queue is now
  handed to every task, the tasks return their warning ids, and the client merges them in
  task order (the order a serial run fires them in); a once-per-run warning shows once.
  New public helpers `lsqvce.util.log.attachWorker` / `mergeInventory` for callers
  running their own parfor.
- **`Noise.MeanType` is retired (user decision).** Since 6.1.6 the arc factors are
  absolute and gamma-bar no longer scales them; setting the field is refused
  (`lsqvce:Options:badGroupField`), and `adjustVCEjoint` no longer takes `'MeanType'`. The
  kernel `lsqvce.engine.arcVarianceFactors` keeps both of the thesis' means.
- **Engine contract (user decision, documented in 6.1.6):** an accFcn driving
  `adjustVCEjoint` must apply `model.arcFactor`.

## 6.1.7 — 2026-10-10 (R7 WP-118 part 1: docs truth, one default per fact, honest numbers)

- **The "95 %" error ellipse covers 95 % (R7-39).** `lsqvce.engine.ellipse` with the
  estimated s0² (the default, `'aposteriori'`) now scales by 2·F(2, r) instead of χ²(2),
  which covered 86 % at redundancy 5 and 90.5 % at 10. With σ0 known (`'apriori'`) χ²(2)
  stays; the two meet at large r. The used constant is returned as `E.k2`.
- **Unseeded joint runs are replayable (R7-37).** With `Seed = []` the joint estimator
  pins one probe seed up front (from the private probe stream, the global RNG untouched),
  probes with it on every pass, returns it as `out.seed`, and `fit` records it in
  `R.meta.seeds.jointSeed`. Before, every pass drew a fresh seed and none was recorded.
  Behaviour change: an unseeded run now uses the same probes on every pass, as a seeded
  run always did.
- **Diagnostics no longer advance your RNG (R7-36).** Past 5000 observations the MDB
  reservoir of `lsqvce.engine.diagnose` drew from the global stream; it has its own fixed
  stream now, so the reported MDB quantiles are reproducible.
- **One default per fact (R7-35).** `diagnose`'s `rMin` defaults to 0.01, the value
  `Options.Diagnostics.RMin` always had (it was 1e-6), and `snoop`'s re-diagnosis inherits
  the thresholds of the D it was given. The documented-but-dead `diagnose` option
  `'verify'` is retired (passing it is an error now). `adjustBlock`'s help no longer
  documents `ShowSeconds` / `ExportDir`, which it never accepted.
- **Every accepted option is documented, every named unit exists (R7-34, R7-35).** Two
  new audit checks: 31 (each Name-Value option an `arguments` block accepts is named in
  its help; 17 were not) and 32 (each `lsqvce.*` unit named in help text, runtime
  messages, the Manual, README, templates and examples exists; 33 named deleted units,
  among them `snoop`'s runtime remedy).
- **The verbose locals line reports the right size (R7-40).** `localCovariance` takes the
  solve's s0² (`acc.localCovariance(sysID, sol.Cxx, sol.s0sq)`, now one formula shared
  with `residuals`); the log line used a pre-solve hint and read ~300× too large. API
  change on this engine method: the third argument is required.

## 6.1.6 — 2026-10-10 (R7 WP-116 part 2: the arc factors weight the solve)

Behaviour change for every `Estimate = "noise"` run with `ArcFactors = true`; no declaration
changes.

- **Per-arc factors weight the solve (R7-13, user decision: Ellmer's convention).** Until
  now an arc's estimated factor was reported but only its stream's mean reached the
  weighting, so a bad arc kept its stream's weight. Each arc now enters the next pass with
  its own variance (stream level × arc factor), and the update (eq. 6.5.28) multiplies by
  that same factor. Measured on one arc 10× in sigma among nine, beside two other streams
  (16 seeds): the true scatter of the estimated parameters is 13–47 % smaller.
- **Arc factors are absolute — a 6.1.5 regression fixed.** The 6.1.5 per-stream form,
  (gamma2_m / gamma-bar_s²) · mean(gamma2_s), counted a stream's level twice as soon as its
  arcs differ: on the fixture above the factors came back 1176 / 11.65 against the planted
  100 / 1, and the other stream's level 30 % low. A factor is now the arc's Helmert factor
  gamma2_m itself (100.0 / 1.00). 6.1.5's own fixtures had homogeneous arcs and could not
  see it. `Noise.MeanType` (gamma-bar) now sets only the scale of the clamp floor.
- A caller driving `lsqvce.engine.adjustVCEjoint` with its own accFcn must apply
  `model.arcFactor` (`lsqvce.engine.jointModelOpt` does, as `covOverride.byArc`); one that
  drops it gets compounding factors and the `arcFactorsDiverging` warning.

## 6.1.5 — 2026-10-10 (R7 WP-116 part 1: the joint-loop policy)

Behaviour changes in `Estimate = "noise"` runs; no declaration changes.

- **One stream's level no longer depends on another's (R7-09).** The arc factors were
  normalised by gamma-bar over every stream's arcs, so an arc-factor stream's reconstructed
  sigma moved with an unrelated stream's truth (1.01 / 0.79 / 1.15). gamma-bar is now formed
  per stream, and an arc-factor stream keeps its own absolute level (the Helmert mean the
  other streams already get). Note: the review's proposed remedy, per-stream gamma-bar alone,
  was measured on its own mirror and loses that level (a stream truly 3× its declared sigma
  came back 1.93); the shipped form recovers 3.02.
- **The joint loop's stop rule covers every free component (R7-11).** It used to watch the
  law alone, so a run with no `"law"` stream stopped "converged" on pass 1 under the
  DECLARED weights and prior groups were never updated. It now also needs every level and
  arc factor to settle (`R.history.levelDev`), and the returned solution is accumulated
  under the returned model — before, x, Cxx and sigma0 were one pass behind the noise model
  handed back next to them.
- **`Noise.Smooth` and `Noise.FreezeRedundancy` do what they say (R7-22).** Both were
  documented and reached nothing. `Smooth` is `psdVCE`'s moving average; `FreezeRedundancy`
  is a real freeze: a frequency whose summed redundancy is below it holds its value. The
  GRAVFIRE recipe's `FreezeRedundancy = 100` now takes effect.
- **A `var` law under `Estimate = "noise"` is refused by name (R7-07)**
  (`lsqvce:fit:varUnderNoise`): the joint sweep has no vector-AR path and used a wrong
  covariance for it silently.
- Per-arc factors still do not weight the solve; that is the next change (R7-13, user decision).

## 6.1.4 — 2026-10-10 (R7 WP-115: robust conventions — wrong-day weights, calibrated multi-axis epochs, honest scales)

Two silently-wrong-numbers fixes, three behaviour corrections, one new diagnostic and
documentation truth. Of the design questions review R7 raised, the user decided all three on
2026-10-10; two ship here, and the third (an additive nugget instead of the multiplicative
Koch fold under correlated laws, R7-10) is a larger change to the loader, the sweep and the
lag kernel and follows in its own release.

- **Clean multi-axis epochs are no longer rejected (R7-08).** The joint loop's robust pass
  reduced a d-axis epoch to its RMS — a χ-shaped value with median near 1 — and the
  influence standardises uncentred, so `igg3` zeroed about 24 % of CLEAN 3-axis epochs (62 %
  at d = 6). The epoch statistic is now χ-calibrated (user decision): the sum of squares is
  mapped through χ²(d) to the equivalent normal deviate with the sign of the rows' sum, so a
  clean epoch is standardised exactly like a clean scalar residual (0.3 % zeroed at any d,
  measured on a mirror); blunders are still rejected. Scalar streams are unchanged.
- **A breakdown warning (user decision).** `lsqvce:robust:atBreakdown`, once per run, when a
  robust pass (stand-alone or joint) zeroes more than 45 % of its observations. The collapse
  refusal is unchanged — it remains the model-error detector 6.1.3 described. Diagnostic only;
  Huber never zeroes a factor, so it cannot fire under Huber.

- **Robust weights could be applied to the wrong day (R7-06).** `lsqvce.util.rowWeightKey`
  printed day ids with `%g`, which keeps six significant digits, so seven-digit day ids —
  YYYYMMDD dates, Julian Day Numbers — collapsed onto one key: with `Robust` on, one day's
  blunder factors were applied to another day of equal length, or a misleading
  `makeLoader:rowWeightSize` blamed the builder. Integer ids now print with `%d` (exact,
  the loader cache's own format), others with `%.17g`. MJD and datenum days (< 1e6) were
  never affected.
- **The robust scale's fallback is no longer silent (R7-21).** When half or more of the
  whitened residuals tie (quantised data, zero-redundancy rows) the MAD is degenerate and
  `lsqvce.kernel.robustWeights` falls back to the std; it now warns once per run
  (`lsqvce:robustWeights:scaleFallback`), returns the path taken as a fourth output, and
  treats a near-tie MAD (a 1e-15 spread) as a tie instead of a scale — which had given every
  ordinary row a factor of ~1e-14 that the collapse guard could not see.
- **The IRLS stop rule is unit-invariant (R7-17).** `adjustRobust` stopped when
  `max|dx| < Tol·max(max|x|, 1)`; the floor of 1 made SH-scale coefficients (~1e-9) stop at
  pass 2, "converged", about one formal sigma from the fixed point. The scale is now the
  larger of `max|x|` and the largest formal sigma: unchanged for parameters of order one and
  larger, and the same pass count in any units.
- **`R.report()` says when a fit is robust (R7-18).** `adjustRobust`'s help promised it and
  the report said nothing; a robust fit's sigma0 and covariance are those of the reweighted
  problem, not an M-estimator's scale or sandwich covariance. One level-1 NOTE
  (`ROBUST_FIT`) now says so.
- **Claims withdrawn.** The joint loop's `collapsed` refusal no longer says "beyond the
  breakdown point" (R7-25: the claim 6.1.3 retired from `adjustRobust` survived there, wrapped
  across two help lines and inside the runtime message), and counts epochs, which is what it
  compares. "The level a robust pass returns is UNBIASED" is withdrawn to an open question
  (R7-12) in the help, the Manual and the test that pinned it; it is settled by measurement
  in a later part. Audit check 28 now reads wrapped help lines and runtime messages, so
  neither route can carry a retired claim again.

## 6.1.3 — 2026-09-21 (R6 WP-110: what the collapse guard really detects, and a 6.1.2 claim corrected)

Documentation, a corrected error message and tests; no estimator behaviour changed.

- **6.1.2 said the collapse guard "can fire only for `n < 2(nc+1)` … it cannot fire at any
  contamination". That is wrong, and this release withdraws it.** The half-sample bound
  6.1.2 proved is real, but it is conditional on `|med| ≤ 3.4·MAD`, and that premise is the
  50 % breakdown point wearing a different hat: measured on deterministic samples, it holds
  at 49 % contamination (ratio 0.50) and fails at 51 % (ratio 518.89). Above breakdown
  `nAlive` falls with the clean fraction and the guard does fire. The bound therefore says
  no more than "we are below breakdown"; 6.1.2 drew a consequence its own premise forbids.
- **What the guard actually detects, measured.** `nAlive` tracks the CLEAN row count to the
  row, and the threshold is an absolute `nc + 1`, so contamination reaches it only when
  almost nothing clean is left (99.8 % on a 3600-row fit against `nc = 6`; 60 % contamination
  of 1000 rows still leaves 400 alive). What does reach it on a full-sized problem is
  residuals offset **bodily** from zero by several MAD — an unabsorbed bias, a datum defect —
  because the influence standardises `vw/s`, not `(vw − med)/s`. That is a model error, not
  an outlier problem, which is what the message always meant to say and now does.
- **`lsqvce:adjustRobust:collapsed`'s message and help are corrected** to name those two
  routes instead of claiming a 50 % breakdown detector. The phrase it used to carry is in
  `superseded_claims.csv`, so it cannot come back. **The threshold itself is unchanged** —
  turning the guard into a genuine breakdown test means choosing a policy, which is the
  user's call, not a release note's.
- **Issue #110's test 6 is closed**: the arc-scoped local fixture with priors that WP-109
  said it needed now exists, and a degenerate arc — every row rejected — is carried by its
  local prior alone, with plain GLS on the same data as the control.

## 6.1.2 — 2026-09-20 (R6 WP-109: the PDF Manual is gated, and issue #110 is closed out)

Documentation and tests; no code in `+lsqvce` changed, no declaration changes.

- **`lsqvce_manual.pdf` carried 6.1.0 into the 6.1.1 release, because nothing read it.**
  It is rebuilt, and audit check 29 now holds its document title to `+lsqvce/Contents.m`
  — the same single version source check 17 (a) holds the README title to. The PDF cannot
  be gated the way the HTML Manual is (check 23 compares a fresh build byte for byte):
  rebuilding the identical source changes the font subset tags, the document `/ID` and
  the xref offsets, which is measured here rather than assumed — 21 of 141 inflated
  streams differ, every one of them accounted for, and not one content byte. The title is
  what can be gated and the title is what drifted. The check fails CLOSED: a missing file
  or an unreadable Info dictionary is a finding, not silence.
- **Issue #110's remaining tests (3–6) are in the suite**, each with its margin measured
  in the same run rather than assumed: planted blunders against the fixture's own
  clean-data error under `igg3` and `bisquare`; the IRLS fixed point on kind C with an
  explicit `p`; and `arma` against the PSD evaluated from it, where the AR-order gap is
  measured by the GLS pair so the robust pair can be asked the only question IRLS can
  answer — that re-weighting does not amplify it.
- **A finding, recorded not fixed: the collapse guard is close to unreachable.**
  `adjustRobust` refuses a pass when `nAlive < nc + 1`, but `robustWeights` standardises
  by the normalised MAD, so at least half of any sample survives a redescending influence
  at its default tuning: `nAlive >= ceil(n/2)` always, and the guard can only fire for
  `n < 2*(nc + 1)`. On a real-sized fit it cannot fire at any contamination. 6.1.1
  reported that prior pseudo-rows padded `nAlive` — a guard failing open; the guard was
  nearly unreachable to begin with. The tuning is the user's, so nothing is changed here.

## 6.1.1 — 2026-09-20 (R6 WP-108: the stand-alone IRLS stops re-weighting priors)

Fix and tests only; no declaration changes, no new capability.

- **`Robust` alone no longer lets a local prior group into the robust sweep.** 6.1 said
  "Prior streams are never re-weighted", and `adjustVCEjoint` did exactly that. The
  STAND-ALONE driver did not: a local group rides as its own stream `prior_<name>`, and
  `adjustRobust`'s residual collection took every stream. Three consequences, all
  measured on CI before the fix (issue #110): the pseudo-residuals entered the MAD scale
  that sets every DATA row's factor; they were counted in `nAlive`, `nDown` and `nZero`,
  so the collapse guard could pass on data already past the breakdown point — a guard
  failing OPEN; and because `priorArc` never applies `rowWeights`, the next pass divided
  a prior residual by a `sqrt(factor)` that had never been applied, so that residual grew
  and its bookkeeping factor fell, pass after pass. `private/lawWhitenedResiduals` now
  skips `prior_` streams, which is the whole fix — `adjustRobust`'s write-back iterates
  the keys that function returns, so nothing is written under a prior key any more.
  The two halves of one rule are now written in the same terms.
- **Not changed, deliberately:** `priorArc` still does not apply `rowWeights`. Making the
  bookkeeping consistent the other way would down-weight a prior for disagreeing with the
  data, which is what a prior exists to resist.
- **`adjustRobust`'s `% shared with adjustVCEjoint` comment was wrong** and is corrected:
  `lawWhitenedResiduals` is called only from `adjustRobust`. The joint loop has its own
  producer and its own exclusion.

## 6.1.0 — 2026-09-09, tagged 2026-09-19 (R6 WP-105: the robust pass inside the joint loop)

Tagged as `v6.1.0` on 2026-09-19 (WP-107): the version string gains its PATCH field
(ci-actions `VERSIONING.md`: `MAJOR.MINOR.PATCH`, one annotated tag per release). The
content is the 6.1 of 2026-09-09; no code changed.

Additive on 6.0; no declaration changes.

- **`Robust` together with `Estimate="noise"` runs — in the Koch convention.** Every pass of
  the joint estimator ends with one residual sweep; the RAW residuals, descaled by each
  epoch's declared accuracy and reduced to one value per epoch, go through the influence
  function (`huber` | `igg3` | `bisquare`, `RobustOpts.Tuning`) per stream, and the factors
  are folded into the per-epoch error information of the next pass (sigma / sqrt(w),
  BEFORE the temporal law) — the one composition C = D (R_t ⊗ I_d) D carries them
  everywhere. Measured: weighting the law-whitened innovation instead (Yang, the
  stand-alone IRLS's convention) neutralises the filter tail around every blunder and a
  spectral estimator cannot recover the colour (0.054 against 0.20 in 16 passes); with the
  same blunders declared through `Cep` the plain estimator recovers it.
  `MaxPasses` is the budget; `RobustOpts.MaxIter` / `Tol` belong to the stand-alone IRLS
  only. Prior streams are never re-weighted.
- What comes back: `R.noise.robust` (`influence`, `scale` and `streams` — one MAD scale per
  stream, `factor` per observation in sweep order, `keys`), `R.history.robustScale` (one
  column per stream) / `nDown` / `nZero` per pass, and the same fields on the `OnPass`
  context. The residuals are standardised by the stream's own MAD, never by the estimated
  level (a global scale, or the level itself, sends a clean stream's level into a spiral —
  measured). Prefer `igg3` / `bisquare` for gross blunders: Huber keeps c/|u| of their
  weight.
- The level is UNBIASED under the robust pass (measured 1.007 on a clean stream under IGG3;
  the retained-observation ratio E[w u²]/E[w] that was expected was refuted by the oracle):
  Ellmer's per-lag redundancy is the model expectation under the same weighted whitening,
  so the weights cancel in alpha². Multi-axis streams are weighted per EPOCH (the RMS over
  the d law-whitened rows, one factor per epoch). This is Koch's iterated re-weighting
  inside a VCE, not Ellmer (2018), whose robustness is preprocessing screening plus the
  arc-wise variance factors (Manual chapter 7, §7.8).
- New identifier `lsqvce:adjustVCEjoint:collapsed` (beyond the breakdown point);
  `lsqvce:fit:notYet` now fires only for `Robust` on a nonlinear (`x0`) batch problem.

## 6.0 — 2026-09-09 (R6 WP-102: the cells are the arcs)

BREAKING for every `buildDay` (the sixth contract change: 2.0 namespace, 3.0
counts derived, 4.0 scoped tables, 5.0 one loop, 6.0 this); everything else
about the declaration — tables, priors, `Cep`, `Options` — is unchanged.

- **`buildDay` returns ONE CELL PER ARC:** `D.A`, `D.B`, `D.l`, `D.t` are
  1 × nArc cells, cell a = arc a (`{A}` for a single arc, empty cells for a
  day without data); `D.Cep` is a matching cell (nEpoch_a × 1 / × d / × d × d
  per arc); `D.nAxis` stays scalar. Arc lengths differ freely. The v5 fields
  `arcIdx` / `arcId` are GONE — a returned `arcIdx`, `arcId` or any field
  beyond `A B l t Cep nAxis` is refused by name (`lsqvce:makeLoader:badDay`),
  as is a matrix where a cell belongs. Driver: the GRAVFIRE preprocessing
  delivers every day as a cell array, one cell per noise-independence piece
  — the day's cells now ARE the arcs, no bookkeeping in between.
- **`arcs` is OPTIONAL.** Omitted, the arcs of a day are its cells (ids
  1 … nArc), learned when the engine's per-day `meta` probe reads the day
  (which the engine reads right after; the cache holds a whole day). Then
  `'day'` template rows and common `'day'` rows live on EVERY calendar day
  of the system, and a calendar day returning no arc in any stream is
  refused at its probe (`lsqvce:makeLoader:emptyDay`) — take it off the
  calendar. `arcs` is REQUIRED only by `'arc'` template rows and the
  engine-tier `LocalsPerArc` (`localPlan`: nl must be known before a read).
- **Declared schedule:** cell a is the a-th scheduled id and the count must
  match (`lsqvce:makeLoader:arcCount`); an **EMPTY cell is the data gap** of
  a scheduled arc (replaces the v5 `arcId` subset). An id outside the
  schedule stays `badArc`. `badArcId` is retired.
- **BREAKING (no aliases):** the string form `arcs = 'fromDay'` is retired
  everywhere (`lsqvce:fit:badArcs` / `lsqvce:makeLoader:badArcs`): omitting
  `arcs` is what it meant.
- `CacheDays` now means what it says: one resident DAY is one slot per data
  stream (the key is day|stream; with one slot a 3-stream day evicted what
  the engine was about to read).
- Test instrument `dev/tests/private/arcCells.m`: the ONE cut from the v5
  row-index day shape to the v6 cells (schedule order, gaps, per-arc `Cep`),
  so every fixture keeps generating whole-day matrices. 558 points (+5).
- Docs: the Manual chapter 3 (§0c, §2.3, §3, §4.0–4.3, §6, §7, §8, §10),
  `makeLoader` / `fit` help, both templates (`gravfire_data` STEP 5 gone,
  STEPS 8–11 = the cut), `dev/GRAVFIRE_INTERFACE.md`, every example
  builder; record: `dev/HANDOVER.md` §12.26.

**Migration v5 → v6**

| v5 | v6 |
|---|---|
| `D.A = A; D.B = B; D.l = l; D.t = t; D.arcIdx = {rows1, rows2}` | `D.A = {A(rows1,:), A(rows2,:)}` … one cell per arc for `A`, `B`, `l`, `t` |
| a single-arc day: matrices, no `arcIdx` | `D.A = {A}`, `D.B = {B}`, `D.l = {l}`, `D.t = {t}` |
| `D.Cep` per day-epoch | `D.Cep = {Cep1, Cep2}` per arc |
| `D.arcId = [1 3]` with `s.arcs = @(day,~) [1 2 3]` | cells `{arc1, {}, arc3}` — the empty cell is arc 2's gap — or a resolver returning `[1 3]` |
| `s.arcs = 4` on a `'system'`/`'day'`-scoped system | omit `arcs`: the four cells are the arcs |
| `s.arcs = 'fromDay'` | omit `arcs` |
| `'arc'` template rows | unchanged: keep the schedule (count, vector or resolver) |

## 5.3 — 2026-09-06 (R6 WP-90: the error information travels with the day)

Additive on 5.2 for every declaration that used `diag` / `arma` / `psd` laws;
BREAKING for two seams, named below.

- **`buildDay` returns the error information with the day, for ANY stream:**
  `D.Cep` is nEpoch × 1 (the per-epoch standard deviation of a scalar stream),
  nEpoch × d (per-axis standard deviations) or nEpoch × d × d (the full
  per-epoch covariance, e.g. from a PPP solution). The axis count is derived
  from it (`nAxis` stays only for a multi-axis stream without `Cep`; given as
  well, it must agree). The covariance spec declares the TEMPORAL law only
  (`diag`, `arma`, `psd`, `acf`) and its level; the toolbox composes
  C = D (R_t ⊗ I_d) D per arc in ONE place (`covModel`), for every law,
  and the joint estimator de-scales the residuals of every arc that carries
  `Cep` — a scalar stream with a per-epoch sigma is estimated correctly
  (before 5.3 that information could not travel with the day at all, and a
  `Cep` on a `diag` stream was silently ignored).
- Validated at the day, with day and stream in the message:
  `lsqvce:makeLoader:badDayCep` (shape), `:axisConflict` (`nAxis` vs `Cep`,
  or a `d` on the spec), `:cepTwice` (a `Cep` on the spec and on the day),
  `lsqvce:covModel:cepOnVar`; a multi-axis day without `Cep` is composed with
  D = I and reported once per run (`lsqvce:makeLoader:noCep`).
- Fixed: the day's `Cep` was sliced with the arc ID instead of the arc's
  position in `arcId` — a ragged day whose ids are not 1:n whitened an arc
  with another arc's epochs (R6 finding E-01).
- **BREAKING (no aliases):** the covariance type `separable` is retired —
  declare the temporal law (`arma` / `diag` / `psd`) and return `Cep` with
  the day (`covModel:unknownType` names the replacement); the
  `epochCovariance` field of `data.systems(k)` is retired (an unknown field
  now — its only consumer, `kernel.residualAcv`, had no engine caller since
  5.0 and is deleted); `makeLoader`'s `EpochCovariance` option and the
  loader's `epochCov` field are gone; `jointModelOpt`'s `Separable` option
  is gone (the loader composes).
- Docs: the Manual chapter 3 (§0b, §3.4, §3.5, §5, §7 "the error information",
  §11), chapter 7, the GRAVFIRE template and guide (positions: `diag` +
  `Cep` from iSpp), `lsqvce_template` (one stub), `demo_grace_dataLoader`.
  Review R6 and the test audit of the chain: `dev/HANDOVER.md` §12.14.

## 5.2 — 2026-09-03 (R5-D8: the parameter-table interface)

Additive on 5.1 except one retirement (`Options.Reg / RegLambda / RegMatrix`).

- **Regularisation is declared, in your columns.** `data.priors(j)` (common
  parameters, L in the columns of A) and **`data.systems(k).priors(l)`**
  (system k's local parameters, L in the columns of B) are one mechanism:
  one struct per group, declared by **`group`** (a regular expression on the
  table's name column — the toolbox builds L, one row per matched row per
  instance) or by **`L`** (used as is), plus `x0`, one `sigma`, `from`,
  `name`, `estimate`; pseudo-observations counted in the redundancy; `estimate = "scale"`
  estimates the group's sigma every pass by Helmert / Koch–Kusche — for local
  groups as the level of their own stream `prior_<name>`, identical to the
  common groups. `from` activates a group at a pass; a label (`group` on
  common priors, the struct on local ones) merges members.
- **`lsqvce.layout(data)`** announces the column layout before `fit`, so L
  can be written in it (`lay.systems(k).col(row, day, arc)`).
- **`data.paramColumns = struct('start', c1, 'end', c2)`** names the validity
  columns of both tables: a `'day'` / `'arc'` row is instantiated only on the
  days it overlaps; `buildDay`'s info carries the rows present.
- **`R.toTable(data)`** writes estimates and standard deviations back into the
  expanded n × 15 tables (columns 9/10 by default) with the instance map.
- The 15-column convention: the scope token is column 1; your 14 columns
  follow. Column 11 `{'Tikhonov', fromPass, sigma[, label]}` is translated by
  the GRAVFIRE template into the declarations above.
- Retired: `Options.Reg / RegLambda / RegMatrix` (19 Options fields);
  `localPriors(row)` and the `priorDay` handle of the intermediate rounds.
- Fixed on the way: an observation type with an underscore was cut by the arc
  key parsers; the prior arc now takes its law through the override path.

## 5.1 — 2026-09-03 (Review R5, levels M6–M7: additive on 5.0)

Additive; no v5.0 declaration or option changes meaning.

- **`data.priors(k).estimate = "scale"`** — a prior group gets its own variance
  component (Koch–Kusche): under `Estimate="noise"` its sigma is updated every
  pass by the Helmert rule σ² ← σ²·Ω_k/r_k (fixed point Ω_k = r_k, "the prior's
  weighted residual sum equals its redundancy"); `R.noise.priors`,
  `R.noise.priorHistory`. This closes the GRAVFIRE guide's `priorGroup` question.
- **`R.compare(Rred)`** — the nested-model F test between two fits
  (`lsqvce.engine.ftest` re-homed onto `Results`).
- **`R.plot("residuals" | "contribution" | "psd" | "convergence")`** in its
  final form; `Options.OnPass` unchanged since 5.0.
- **`examples/demo_gravfire.m`** — a GRAVFIRE bundle through
  `template/gravfire_data.m` and `fit` (pipeline check or `Estimate="noise"`
  with `OnPass` and `CacheDays=Inf`); every demo explains v5 in its prose.
- **Test tree**: every public unit has its seven categories tested or argued
  (`dev/tests/category_exemptions.csv` is the review surface); new suites
  `test_metamorphic` (invariances) and `test_speedBudget` (ratio budgets
  pinned to their own measurements). 523 points, 26 suites.
- **Docs**: `lsqvce_theory.pdf` at the root is the R5 theory (54 pp.); the
  HTML manuals carry a v5 banner naming the binding authorities;
  `DATA_CONTRACT.md` v5 documents priors, per-stream policy and `OnPass`.

Development record: `dev/HANDOVER.md` §12.11 (the R5 closing summary).

## 5.0 — 2026-09-03 (Review R5, levels M0–M5: the one-loop redesign)

**Breaking release.** The toolbox has ONE accumulator and ONE stochastic-model
estimator; the option space is 14 knobs + 5 typed groups; every retired name
has a successor in the migration table below (the authority is
`dev/HANDOVER.md` §12.10). Nothing on the data contract changed shape — v4
declarations still fit — but three things were added to it.

### What changed on the user surface

- **`Options` v5** — `Estimate="none"|"noise"` replaces `Model`; the v4 model
  words are explicit policy (`"vce"` = `Estimate="noise", Pool="stream",
  ArcFactors=false`; `"pooled"` = `Pool="all", ArcFactors=false`; `"ellmer"` =
  `Estimate="noise"`). One budget: `MaxPasses` / `Tolerance` (`MinPasses` was
  retired in R6 WP-101: it never had a consumer).
  `Robust` names the influence (`"huber"|"igg3"|"bisquare"`); groups
  `Noise`, `RobustOpts`, `Diagnostics`, `Parallel`, `GaussNewton` merge partial
  structs onto defaults. `CacheDays` (resident built days) and `OnPass` (see
  below) are new. `fit(data, Name=Value)` is the same mechanism as the object.
- **`Options.OnPass = @(ctx)`** — your inspection hook after EVERY pass of every
  iterative run (variance components, IRLS, Gauss–Newton): the pass state, the
  live accumulator, your declaration and every system's `processingOptions`;
  return `false` to stop (`stopReason = "userStop"`). Successor of the v4
  `plotPostModel` / `plotPostSolve` hooks and of `ShowSeconds` / `ExportDir`;
  post-mortem plotting is `R.plot("residuals"|"contribution"|"psd"|"convergence")`.
- **`data.priors(k) = struct('L', 'x0', 'sigma', 'name')`** — declared
  pseudo-observation groups on the common parameters. A free network's datum
  is the prior `L = B'` on the null space (the exact minimum-constraint
  solution, weight-independent); `DatumB` / `Sigma0` / `Sparse` are gone.
- **Per-stream policy on the covariance spec** — `estimate = "none" | "scale"
  | "law"` and `arcFactors` (logical); `Pool = "all" | "obsType" | "stream"`
  decides which streams share one estimated spectral shape; every stream keeps
  its own level. A law with a per-row `sigma` is a fixed structure (`"scale"`
  only).
- **Every fit kind rides the block path.** In-memory (`A, l, p`) and batched
  (`batchFcn`) problems become one local-free streaming system; `Robust`
  works on the streaming kind; Gauss–Newton (`x0`) is a policy on the same
  loop. A full SPD weight matrix per batch is whitened exactly. Per-observation
  diagnostics of a batch problem are reported in the WHITENED frame
  (`R.diag.v = sqrt(p).*v`); only their sums are frame-invariant.
- **`Results`**: `stopReason` (`"single" | "converged" | "maxPasses" |
  "userStop"`), `meta.warnings` (every warning identifier of the run),
  `noise.Cgroups` / `noise.stream` / `noise.model{k}.sigma`, `history.streams /
  sigma / white / whitePval` on every iterative run (so `convergence()` and
  `summary()` work for `Estimate="noise"`), `plot(kind)`. `R.source.verdict` →
  `R.report().verdict`; `history.nOuter` → `history.nIter`.

### What changed underneath (why the numbers may differ from 4.0)

- The joint (Ellmer) estimator is the one estimator (finding H-03): the
  separate mode's AR law fitted to the biased raw-residual ACF is retired
  (R5-D4). Acceptance is the physical oracle (planted fields and stream sigmas
  recovered within their formal errors), not equality to 4.0.
- Multi-axis streams are handled per axis after the per-epoch descale (H-15);
  the joint law lives in the PSD domain, positive by construction and resampled
  to every arc's length (H-17); the pooled shape is normalised and each stream
  carries its own level; per-lag sums are FFT autocorrelations (a pass on a
  two-day GRACE-FO bundle at 2160 lags: 12 min → 39 s).
- The QR backend refuses rank deficiency by `rank()`'s tolerance (not at the
  roundoff boundary — the refusal was platform-dependent) and folds a
  rank-deficient contribution in its range space (no `\` on a singular factor,
  no warning ahead of the refusal) (H-16).
- Error doctrine: every catch refuses, reports or counts; `util.callUser`
  wraps user functions with the original error as cause; the identifier
  catalogue is closed (`dev/docs_src/error_ids.csv`); skip RATES refuse
  (`Noise.SkipThreshold`).

### Migration table (v4 → v5)

| v4 | v5 |
|---|---|
| `Model="gls"` | `Estimate="none"` (default) |
| `Model="ellmer"` | `Estimate="noise"` (+ `Pool`, `ArcFactors=true` default) |
| `Model="vce"` | `Estimate="noise", Pool="stream", ArcFactors=false` |
| `Model="pooled"` | `Estimate="noise", Pool="all", ArcFactors=false` |
| `Robust="irls"` + `Influence` | `Robust="huber"|"igg3"|"bisquare"` |
| `Tuning`, `RobustMaxIter`, `RobustTol` | `RobustOpts.Tuning / .MaxIter / .Tol` |
| `MaxOuter`, `JointMaxIter`, `VceMaxIter`, `NonlinearMaxIter` | `MaxPasses` |
| `JointTol`, `OuterTol`, `VceTol`, `TolX` | `Tolerance` |
| `TolF`, `Lambda0` | `GaussNewton.TolF / .Lambda0` |
| `NumLags` | `Noise.Grid` |
| `PsdMode` | `Noise.Whitener` (`"toeplitz"` default, as before) |
| `MeanType`, `Smooth`, `MinRedundancy`, `SkipThreshold` | `Noise.MeanType / .Smooth / .FreezeRedundancy / .SkipThreshold` |
| `Parallel`, `DistributeDays`, `ParallelThreshold` | `Parallel.Mode / .DistributeDays / .Threshold` |
| `Alpha`, `Bonferroni`, `Beta`, `RMin`, `RWeak`, `Test`, `RobustCheck`, `Store`, `NWorst` | `Diagnostics.<same>` |
| `LogFid`, `LogFlag`, `ShowSeconds`, `ExportDir`, `plotPostSolve`, `plotPostModel` | `OnPass`; `R.plot(kind)` (`LogFile` was offered here and retired in R6 WP-101: it never opened a file) |
| `DatumB` | `data.priors = struct('L', B.', 'sigma', 1, 'name', 'datum')` |
| `Sigma0`, `Sparse`, `Verify`, `Progress` | retired |
| `GroupBy`, `ComponentsPerEpoch`, `MaxOrder`, `AcvTol`, `MaxLag`, `WhiteAlpha` | retired with the separate mode (groups are streams; the axis count comes from the builder; whiteness is fixed at K = 20, α = 0.05) |
| `Norm`, `C0`, `StepSchedule`, `AvgTail` | engine-tier only (`lsqvce.engine.adjustVCEjoint`) |
| `lsqvce.engine.adjust / robust / solveNonlinear / Accumulator / adjustVCE / Combine / BatchIndex / varComp / vceHutchinson` | gone; `fit` on the block path; `lsqvce.engine.adjustBlock / adjustRobust / adjustNonlinear / adjustVCEjoint / batchSystem` at engine tier |
| `R.source.verdict`, `history.nOuter` | `R.report().verdict`, `history.nIter` |

Full ledger, decisions and round notes: `dev/HANDOVER.md` §12; theory:
`lsqvce_theory.pdf` (R6 WP-101: it shipped under a second name in dev/ as well; that duplicate is gone -- one document, one name).

## 4.0 — 2026-08-31 (Level L11, WP-75: SCOPED parameter tables)

**Breaking: both parameter tables gain a binding scope column.** Column 1 of
`commonParams` (`'common' | 'day'`) and `localParams` (`'arc' | 'day' |
'system'`) is now the scope token; rows are **templates** the toolbox expands
over the calendar — nc and nl are derived from the expansion, and every column
offset reaches `buildDay` through its 6th argument. Driver: the GRACE
parameterization (arc-wise ICs, day-wise ACC scale/bias/drift, per-arc KBR
empiricals), which under v3 forced day-wise locals into hand-numbering.

### Migration table (v3 declaration → v4)

| v3 | v4 |
|---|---|
| `commonParams = T` (R×C, free) | prefix every row: `[repmat({'common'}, R, 1), T]` |
| `localParams = T` (R×C, free) | prefix every row: `[repmat({'system'}, R, 1), T]` |
| `localsPerArc = n` | delete the field; declare n `'arc'` rows in `localParams` |
| `localsPerArc = @(day, arcId) n` | not expressible as a template (constant rows per arc); heterogeneous per-arc counts stay ENGINE-tier (`lsqvce.engine.makeLoader(..., 'LocalsPerArc', fh)`) |
| day-wise parameters hand-packed into `commonParams` or hand-numbered in `localParams` | one `'day'` template row each; the toolbox numbers them |

### Layout (binding)
A: `['common' rows, table order | day blocks, calendar order]` — a calendar
day gets a block iff at least one system's **schedule** carries an arc that
day. B per system: per day `['arc' blocks, sorted ids | 'day' block]`, then
one trailing `'system'` block; a system day without arcs gets no blocks.
`info` (buildDay's 6th argument): `.nl .arcId .offsetOf(arcId) .count
.dayOffset .dayCount .sysOffset .sysCount .nc .commonCount .commonDayOffset
.commonDayCount`.

### Error-id changes
New `lsqvce:fit:badScope` (bad scope token; the localParams case is
makeLoader's `badScope`, relabelled). `lsqvce:fit:localCount` retired (one
local mechanism). `localsPerArc` on the USER surface is `unknownField`.

### Unchanged
`buildDay`'s signature and return contract, covariance specs, arc schedules,
`epochCovariance`, `lsqvce.Options`, `lsqvce.Results`, all engine numerics;
the ENGINE tier keeps v3 free-form tables and `LocalsPerArc` (default
`Scoped=false` on `lsqvce.engine.makeLoader`).

---

## 3.0 — 2026-08-31 (Level L10, WP-72: the v3 data contract)

**Breaking: the `data` contract states nothing the data already states.** Every
count is derived from the one place it already lives — nc from the
`commonParams` rows, nl from the `localParams` rows or the `localsPerArc` plan
(neither → nl = 0), the observation types from `fieldnames(covariance)`, u from
`numel(x0)` on the nonlinear batched kind. The retired v2 fields are **unknown
fields** now, so a v2 declaration fails loudly with `lsqvce:fit:unknownField`
naming the field. `lsqvce.engine.*` keeps the v2 conveniences (counts,
`'fromDay'`, `cache`, plot hooks, raw loader closures) — they left the USER
surface, not the toolbox.

### Migration table (v2 declaration → v3)

| v2 | v3 |
|---|---|
| `data.common.commonParams = T` | `data.commonParams = T` |
| `data.common.numCommon = nc` | `data.commonParams = num2cell((1:nc).')` — the table is the declaration |
| `s.observationTypes = {'kbr','pos'}` | delete — the types are `fieldnames(s.covariance)`, in declaration order |
| `s.systemId = 'C'` | `s.id = 'C'` (optional; defaults to the system's index) |
| `s.numLocal = nl` | `s.localParams = num2cell((1:nl).')` — or `s.localsPerArc` |
| `s.numLocal = 0` | delete — nl = 0 is the default (neither `localsPerArc` nor `localParams`) |
| `s.localParams` (R×14) | unchanged name; the fixed 14 columns are gone — R×C, C your convention (same for `commonParams`) |
| `s.arcs = 'fromDay'` | a resolver `@(day,obsType)` (`'fromDay'` is ENGINE-tier only: `lsqvce:fit:badArcs`) |
| `s.cache`, `s.plotPostSolve`, `s.plotPostModel` | delete — ENGINE-tier (`lsqvce.engine.makeLoader`), not data declaration |
| `s.loader = @(d,ot,a,uo) ...` | ENGINE tier: call `lsqvce.engine.adjustBlock` / `adjustVCE` with the struct-array form directly |
| `data.numUnknowns` with `data.x0` (kind C) | delete — u = `numel(x0)` (`lsqvce:fit:uCount` if both); the linear form keeps `numUnknowns` |

### Error-id changes
`lsqvce:fit:commonCount` → `lsqvce:fit:noCommon` (commonParams missing or not a
2-D cell). `lsqvce:fit:localCount` now means *both* local ways declared. New:
`lsqvce:fit:badId`, `lsqvce:fit:badArcs`, `lsqvce:fit:uCount`.
`lsqvce:makeLoader:paramCols` retired (no fixed column count to violate).

### Unchanged
`buildDay`'s six-argument signature and return contract, the covariance specs,
arc schedules (static and resolver), `arcInfo`, the multi-axis rules incl. the
`epochCovariance` hook (§7 documents *why* it exists: the exact
separable-descale path of the VCE Shape step), `lsqvce.Options`,
`lsqvce.Results`, all engine numerics.

---

## 2.0 — 2026-09-01 (R4 track, Levels L1–L6)

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
   working. The toolbox's own error and warning identifiers
   (`lsqvce:<unit>:<rule>`) are unchanged; only the MATLAB-generated ids of
   `validateattributes` (`MATLAB:lsqvce:engine:<fn>:<check>`) now carry the
   qualified name, because they are built from it.
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
guard). `buildfile.m` (`buildtool`) existed L5-L8 and was removed in L9: an
mltbx is not foreseen; the gate is `cd('dev/tests'); ci_gate`, the audit
`python3 dev/docs_src/audit_tree.py`.

### Fluent configurators removed (R4-D2, WP-64b, Level L6)
`lsqvce.engine.adjustVCE` is **constructor-only**. The six chained
configurators are gone (no aliases); their option groups are Name-Value
**structs** on the constructor, merged into the engine defaults with the same
case-insensitive key matching the methods used:

| 1.x fluent call | 2.0 constructor Name-Value |
|---|---|
| `.groupComponentsBy(key)` | `'GroupBy', key` |
| `.stochasticModel('MaxOrder',p, ...)` | `'Shape', struct('MaxOrder',p, ...)` |
| `.varianceComponents('MaxIter',n, ...)` | `'Scale', struct('MaxIter',n, ...)` |
| `.iteration('MaxOuter',k,'Tol',t)` | `'Outer', struct('MaxOuter',k,'Tol',t)` |
| `.joint('Probes',m, ...)` | `'VceMode', 'joint', 'Joint', struct('Probes',m, ...)` — the method used to switch the mode; the group alone does not |
| `.componentsPerEpoch(map)` | `'Axes', struct(obsType, count)` |

`P.plots(...)`, `run()`, `empiricalStochastic()`, `validate()` are unchanged.
Through the USER front door nothing changes: `lsqvce.fit` passes exactly these
groups from `lsqvce.Options`. The rewrite of the toolbox's own tests and
examples was scripted (`dev/docs_src/refactor_fluent.py`, kept as protocol).

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
