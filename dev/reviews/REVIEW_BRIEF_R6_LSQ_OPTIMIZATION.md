# lsqvce — Reviewer's brief R6: least-squares adjustment and optimization

*Created 2026-09-06 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Fable 5.1 — review track R6, round 0 (the brief). Tree reviewed
for this brief: `main` at `af842e0`, toolbox version 5.2 (03-Sep-2026),
31 suites under `dev/tests/`, CI green, audit exit 0. No MATLAB was run while
writing it (tier B): every code statement below was READ, none was executed.*

This file is a BRIEF, not a status document. Status lives in `dev/HANDOVER.md`
only (§13 carries the R6 ledger); this brief tells the reviewer what to read,
what the code claims, what to check, and how to report.

---

## 1. Your role

You are the **domain-scientist reviewer** for the lane *least-squares adjustment
and optimization*: estimation theory, the variance-component estimator, the
iterative policies (VCE, IRLS, Gauss–Newton), convergence, identifiability,
numerical conditioning, and the plausibility of every number the toolbox
prints. A separate **software-engineering** reviewer covers architecture,
tests-as-tests, layering and documentation mechanics; do not cross into that
lane, and do not expect it to cover yours.

Rules that bind every finding (from the project's review discipline):

1. **Findings only.** You modify no code. A finding is a statement with a
   symptom; the implementation delegate turns it into a work package.
2. **Reproduced or inferred, on every line.** A finding you computed
   (ran a test, differenced a Jacobian, evaluated an estimator on a fixture)
   is REPRODUCED. A finding you derived from reading is INFERRED. Never mix.
3. **Falsifiable.** Each finding names what is wrong, why (the estimation
   theory), and the observable symptom that would prove or disprove it.
4. **Rank by severity.** A wrong estimator that silently changes every result
   outranks a suboptimal stopping rule. Lead with the worst.
5. **Say what checks out.** An area examined and found sound is reported as
   such — one line, what was verified and how. A report with only defects is
   indistinguishable from one that did not look.
6. **A recorded decision is not a finding.** Re-proposing something the record
   already decided, with its reasoning, is noise. CONTRADICTING a recorded
   decision is a finding, stated as a contradiction (§7 lists them).
7. **A green suite proves nothing about the theory.** A test whose expected
   value was computed with the same defect passes. Use the oracles (§6).

---

## 2. Read first, in this order

| # | Source | What it settles |
|---|---|---|
| 1 | `dev/HANDOVER.md` §0, §12.2, §12.4, §12.11 | standing rules; what R5 confirmed as RIGHT; the v5 design in one paragraph; the closing summary of the last track |
| 2 | `dev/HANDOVER.md` §12.3 (H-01…H-17), §12.6 (R5-D1…D8) | every finding of the last review with its disposition; the decisions in force |
| 3 | `lsqvce_theory.pdf` (root, 54 pp.; LaTeX in `dev/docs_src/theory_r5/`) | THE authority for every formula: Part I ch. 1 linear model (GM model, whitening, NEQ vs QR, block-arrow elimination, datum/priors, nonlinear models), ch. 2 VCE (Helmert, redundancy, per-lag components, arc factors, group weights, Koch–Kusche, convergence/identifiability), ch. 3 stochastic modelling (Toeplitz, ACV↔PSD, whitening, multi-axis), ch. 4 quality (Baarda, Pope, reliability, whiteness); Appendix C the option mapping as implemented |
| 4 | `html/lsqvce_manual.html` ch. 4 (Options), 6 (priors), 7 (stochastic models), 8 (diagnostics), 11 (theory pointer) | what the USER is promised |
| 5 | Ellmer (2018) §6.5.1–6.5.7 (eqs. 6.5.7, 6.5.16–6.5.28, 6.5.48, 6.5.54, 6.5.62); Kvas (2019) §2.2–2.3, §3.3 | the method the joint estimator claims to implement; equation numbers in help texts refer to Ellmer unless said otherwise (`dev/references/README.md`) |
| 6 | `dev/DEVELOPMENT.md` §"Diagnostics reference", §"Design decisions, and the evidence", §"WP-79 measurements" | Pope's tau under contamination; QR vs NEQ decided by cond(N); the measurements deletions cite |

Where the theory PDF and a help block disagree, the PDF is the claim to
test and the help block is a documentation finding (hand it to the SWE lane).

---

## 3. The tree you are reviewing

**User surface (4 names):** `lsqvce.fit` (the one verb), `lsqvce.Options`
(14 knobs + 5 groups), `lsqvce.Results` (x, sigma0, sol.Cxx/dof/vtv, diag,
history, noise, stopReason, meta; methods report / summary / convergence /
residuals / contribution / compare / plot / toTable), `lsqvce.report`.

**One loop, three policies** — `fit` dispatches on `Options.Estimate`,
`Options.Robust`, `data.x0`:

| policy | driver | what iterates | stop rule (read from code) |
|---|---|---|---|
| `Estimate="none"` | `engine.adjustBlock` | nothing: one generalised solve under the declared laws | — |
| `Estimate="noise"` | `engine.adjustVCEjoint` | PSD amplitudes per pool group, per-stream level, per-arc factors, prior-group sigmas | `max_j |alpha2_j − 1| <= Tol` over the PSD amplitudes (adjustVCEjoint.m:429–431), `Tol = Options.Tolerance = 1e-2`, budget `MaxPasses = 10` |
| `Robust="huber"|"igg3"|"bisquare"` | `engine.adjustRobust` | IRLS: factors from `kernel.robustWeights` (MAD scale) on the law-whitened residuals | `max|dx| < Tol·max(max|x|, 1)`, `Tol = 1e-8`, `MaxIter = 30` (adjustRobust.m:129) |
| `data.x0` present | `engine.adjustNonlinear` | Gauss–Newton with Levenberg damping through the accumulator's `Reg="lm"` | `max|dx_j/sigma_j| < TolX (1e-8)` or relative vPv change `< TolF (1e-10)`; reject on vPv increase, λ×10 up to 1e12 → diverged, λ/10 on accept (adjustNonlinear.m:93–141) |

`Robust ≠ "none"` together with `Estimate="noise"`, and `Robust` together with
`x0`, are REFUSED by name (`lsqvce:fit:notYet`) — the one-loop composition is
a later level, not a defect.

**One accumulator:** `engine.BlockAccumulator` (block-arrow N with per-scope
elimination of local parameters, `preeliminate` at :1195; day-order fold
`foldDayPart` :1026; `reducedNormalMatrix` :683; `solve` :450) over two
backends, `engine.NormalEqSystem` and `engine.QRSystem` (default `"qr"`;
range-space square root of a rank-deficient contribution in `QRSystem.merge`,
H-16). Regularisation is DECLARED, never an option: `data.priors(j)` and
`data.systems(k).priors(l)` are pseudo-observation groups `L·x ~ x0` with
weight `1/sigma²`, ridden as one more observation stream (`prior_<name>`); a
free network's datum is the prior on its null space.

**The estimator's kernels:** `kernel.redundancyLags` (per-lag Ω_n, s_n from
one shared probe block, all lags by one zero-padded FFT per column),
`engine.psdVCE` (eqs. 6.5.19–6.5.27, DCT-I domain), `engine.arcVarianceFactors`
(eq. 6.5.28, gamma normalisation), `engine.sharedProbe` (Z Hutchinson probes,
seeded), `engine.jointModelOpt` (the PSD-domain law on normalised frequency,
floored at eps·max), `engine.covModel` (places the law on the arc's physical
axis with its dt; AR/ARMA/Toeplitz/separable whiteners), `kernel.dctI`,
`kernel.psd2acv`, `kernel.levinsonDurbin`, `kernel.separableWhiten`.

**Diagnostics:** `engine.diagnose` / `diagnoseAccum` (hat matrix, redundancy
numbers, sigma0 tests), `engine.stat`, `engine.snoop` (Baarda / Pope),
`engine.systematic`, `engine.ellipse`, whiteness (Ljung–Box K = 20 collected
in the joint sweep), `engine.ftest` reached as `Results.compare`.

---

## 4. What the code claims — and the derivation to repeat

Every row is a claim you can re-derive from first principles and compare with
the implementation. Compare with the CODE, not with the docstring.

### 4.1 The linear model and the solve
- **Gauss–Markov with a general Σ.** `x̂ = (AᵀΣ⁻¹A)⁻¹AᵀΣ⁻¹l`, `Σ_x̂ = σ₀² N⁻¹`,
  `σ̂₀² = vᵀΣ⁻¹v / (m − u)`. Whitening `W⁻ᵀ` per arc from `chol` of the arc's
  Toeplitz (or separable) covariance; the toolbox accumulates whitened rows.
  Check: the whitened residual `R.diag.v` on kind B equals `sqrt(p).*v_flat`
  (pinned; verify the pin is what you would define).
- **NEQ vs QR.** Same x̂ to roundoff; QR is the default because the decision is
  about `cond(N) = cond(A)²` (DEVELOPMENT §"The decision is about cond(N)").
  Check: on an ill-conditioned fixture (`demo_gravity_design`, polar gap)
  difference the two backends' x̂ against the formal error, not against each
  other's tolerance.
- **Block-arrow elimination.** Local parameters y per (system, day, arc)
  eliminated into the common block: `N_red = N_cc − N_cl N_ll⁻¹ N_lc`. Check
  the local solution y and `Cyy` are RECOVERED, not re-solved, and that the
  datum of a local block with `nl > rank` is refused, not ridged.
- **Rank deficiency.** `QRSystem.merge` takes the square root of a
  rank-deficient contribution in its range space (eigenvalues above
  `rank()`'s tolerance); `solve()` refuses `lsqvce:qr:singularR`; the NEQ
  backend `lsqvce:acc:rankDeficient`. Check the refusal threshold is a CHOICE
  of the code (`max(size)*eps(max(dR))`), never LAPACK's behaviour at the
  boundary (H-16, the two-platform lesson) — and that a legitimately
  rank-deficient single system whose SUM is full rank is accepted.
- **Priors and datum.** A prior group is a pseudo-observation stream: it
  enters accumulation, elimination, the sweep, `residuals()`, and the
  redundancy count with NO second path (D8-2 principle). Minimum-constraint
  datum `L = Bᵀ` on the null space gives the exact minimum-norm solution for
  ANY weight when B spans the null space. Check weight-independence, and that
  a prior with `from > 1` in a single-pass run is reported as not applied.

### 4.2 The joint variance-component estimator (Ellmer §6.5.7, ten steps)
- **Per-lag components.** Ω_n = σ²_m C_n (êᵀΣ⁻ᵀV_nΣ⁻¹ê), s_n = σ²_m C_n
  tr(R V_n), the redundancy R through Z probes (Hutchinson). α²_j = Ω̃_j/s̃_j
  in the DCT-I domain, Ŝ_j = α²_j S_j. Re-derive the trace estimator's
  expectation and its variance in Z; check the probe block is shared across
  arcs (one draw per pass, seeded into `R.meta`) and that the FFT lag sums
  equal the loop form (pinned to 1e-10 in `test_ellmer`).
- **Arc factors and the level (eq. 6.5.28).** γ²_m = α²_m σ²_m is the arc's
  new absolute factor; the stream LEVEL is the Helmert factor = mean over its
  arcs of γ²_m written to every arc (`adjustVCEjoint.m:380–391`); the pooled
  shape is normalised to `C(1) = 1` so level and shape do not double count.
  Check: the product level×shape×arcFactor is the ONLY place the variance
  enters the whitener (3b/4b defects were exactly a missing factor here).
- **Fixed point.** Alternating (separate) and joint updates share a fixed
  point only when the blocks are separable, and here they are not (help
  block, adjustVCEjoint.m:13–26). Check on a fixture with a KNOWN law: does
  the joint loop return to the truth when STARTED at the truth (corr ≈ 1),
  and where does it land from the default start (`C0 = [1; 0.3; 0.3; …]`,
  a deliberately non-degenerate start — eq. 6.5.62's white start is a
  degenerate fixed point, `psdVCE` warns `degenerateStart`)?
- **Convergence rule.** The loop stops on the PSD amplitudes only
  (`max|α²−1| <= Tol`); arc factors, levels and prior sigmas are NOT in the
  criterion. Decide whether that is a correct stopping rule for the joint
  problem or a partial one, and what symptom a partial one produces
  (levels still moving at "converged").
- **Step schedule.** `StepSchedule="decay"` = Robbins–Monro Γ_it = 2/(2+it)
  through `psdVCE 'Gamma'`; the help block calls damping "damage limitation"
  and says the closed loop does not beat open-loop pass-1 quality on the
  tested fixtures. Judge that claim against the theory of stochastic
  approximation for this estimator — is it a property of the estimator or of
  the fixtures?
- **Multi-axis streams.** Axis k = rows k:d:m, descaled by Cep^{-1/2}, per-lag
  statistics ADDED over axes, one factor per arc (H-15). Re-derive: is the
  separable model Σ = Cep ⊗ C_t (or C_t ⊗ Cep — check the row ordering) what
  the code whitens with?
- **Pooling.** `Pool="all"|"obsType"|"stream"` decides which streams share one
  SHAPE; `estimate="none"|"scale"|"law"` per stream decides what is free.
  Check the identifiability claim: a shape is identifiable only from arcs
  long enough relative to `NumLags`; what does the code do below that
  (skip → `SkipThreshold` refusal at 20 %)?
- **Skips are results.** A non-PD arc is skipped and COUNTED; a pass that
  skips more than `SkipThreshold` refuses. Check a skipped arc keeps its
  previous factor (documented) and that this does not bias the pooled shape.

### 4.3 Prior sigmas by VCE (Koch–Kusche)
`Ω_k = ‖(Lx − x0)/σ‖²`, `r_k = p − tr(L N⁻¹ Lᵀ)/σ²` with `N⁻¹ = Cxx/s0sq`,
`σ² ← σ²·Ω_k/r_k`, from the PREVIOUS pass's solution; groups with `r_k <= 0`
are left alone and flagged; a single-row group cannot estimate its σ.
Re-derive from Helmert's simplified estimator; check the trace uses the
covariance of the CURRENT linearised problem (not the damped one, not a
previous pass's law) and that `s0sq` is the a-posteriori factor the same
pass reports.

### 4.4 IRLS (`adjustRobust`)
Influence Huber / IGG3 / bisquare, scale = MAD of the LAW-whitened residuals
(Yang 1994), factors applied INSIDE each arc's whitener, floored at 1e-12 so
a zeroed row's residual stays recoverable. Check: (i) the fixed point on a
diagonal law equals the textbook IRLS fixed point (pinned to 1e-10); (ii) on
a correlated law the factors act on whitened residuals — is that the
estimator you would write, and what does it do to the correlation structure
of the reweighted arc? (iii) `sigma0` of the reweighted problem is what
`report()` prints — say whether that is the honest scale.

### 4.5 Gauss–Newton with Levenberg damping (`adjustNonlinear`)
Step from the damped normal equations `(N + λI)dx = n`; objective
`F(x) = lᵀPl` of the whitened reduced observations (a sum, not a difference —
the NaN-honest vtv is not usable here, defect 5b′-5). Accept iff `F` fell;
otherwise λ×10 and re-linearise at x. Final Results = UNDAMPED solve at x.
Check: (i) a finite-difference Jacobian test on a fixture (trilateration,
exponential decay) — does the loop converge to the right point at the rate
GN predicts; (ii) the covariance reported is that of the linearised
problem — is it also the covariance of the ESTIMATE when the last accepted
step was damped; (iii) `Lambda0 = 0` (pure GN) disables step rejection
entirely (`opt.Lambda0 > 0 && …`, :102) — decide whether pure GN without a
descent check is acceptable or a missing refusal; (iv) a rejected step
consumes one of `MaxIter` (:93–116): symptom `stopReason = "maxPasses"` on a
problem that would have converged with more rejections.

### 4.6 Diagnostics
Baarda's w-test and Pope's tau on the whitened residual with the redundancy
numbers `r_i` from the hat matrix (exact up to nc ≈ 4000, probe above,
WP-79); Pope's tau goes blind under contamination (DEVELOPMENT). Check the
per-observation diagnostics of a correlated batch live in the WHITENED frame
and only sums are frame-invariant — is that stated where a user reads
`R.diag`? Whiteness: Ljung–Box on lags 1..20 at α = 0.05 per stream — is
K = 20 defensible at the arc lengths in play (2160 epochs, 10 s)?

---

## 5. Checklist (tick each with REPRODUCED / INFERRED / NOT REVIEWED)

- [ ] Units, signs, frames: whitening direction (`W⁻ᵀ` vs `W⁻¹`), residual sign
      `v = Ax̂ − l` vs `l − Ax̂` consistent between accumulator, sweep, and
      `snoop`; PSD scaled by `dt` on the physical axis; DCT normalisation
      (`'ortho'` default vs `'ellmer'`) consistent between forward and inverse.
- [ ] The Helmert redundancy: `tr(R V_n)` estimated, not the trace of H alone
      (WP-32 lesson).
- [ ] Level / shape / arc-factor decomposition: no double counting, no missing
      factor, on a planted 1 / 2 / 0.5 sigma fixture (`demo_gravity`,
      `test_vcePolicy` PLANTED-WEIGHTS REFERENCE, ±25 %).
- [ ] Convergence criterion of each of the three drivers: complete, honest
      `stopReason`, budgets consumed as documented.
- [ ] Priors: pseudo-observation semantics exact (redundancy, residual,
      diagnostics include them); Koch–Kusche step correct; single-row and
      `r_k <= 0` groups handled as documented.
- [ ] Rank deficiency and conditioning: refusal thresholds are the code's
      choice; the condition evidence in `R.diag`/verdict is the number a
      user needs.
- [ ] Multi-axis separable model: Kronecker order matches row order.
- [ ] Hutchinson probes: Z = 100 default; the variance of the trace estimate at
      the redundancies in play; seed in provenance.
- [ ] Plausibility on the physical oracle (`test_physOracle`): max|z| and mean
      z² of coefficient recovery, error degree variances vs formal — "too
      good" is also a finding.
- [ ] What the estimator CANNOT resolve (theory ch. 2 §"Convergence,
      identifiability"): is every stated limitation real, and is any real
      limitation unstated?

---

## 6. Instruments available to you

| instrument | how | what it establishes |
|---|---|---|
| the gate | `matlab -batch "cd('dev/tests'); ci_gate"` (or CI on a draft PR) | all 31 suites loaded, zero failures, warning inventory, per-suite seconds |
| a single suite | `cd dev/tests; runtests('test_ellmer')` — see also `test_vceMath`, `test_vceFunction`, `test_vcePolicy`, `test_physOracle`, `test_precision`, `test_metamorphic`, `test_reliability`, `test_fit` | the pinned claims per area; READ what a point asserts before trusting its green |
| the oracles | O1 closed forms (dense reference solves in tests); O2 the Python mirrors `dev/validation/*.py` (`wp33_loop_mirror.py` is the joint loop); O3 the physical oracle (`examples/demo_gravity_design.m`, `test_physOracle.m`); O7 the GRAVFIRE bundle (the user's, private) | independent authorities; O4 (a GROOPS NEQ file) is OPEN — no external NEQ parity exists |
| watching a run | `lsqvce.fit(data, Estimate="noise", OnPass=@(ctx) …)` — `ctx` carries `.x .Cxx .sigma0 .s0sq .C .Cgroups .sigma2 .arcKeys .streams .sigma .white .whitePval .maxDev .skips .acc`; return `false` to stop | the trajectory of every component per pass; the live accumulator for any engine call |
| the record | `R.history`, `R.convergence()`, `R.summary()`, `R.noise.priorHistory`, `R.plot("psd"|"convergence"|"residuals"|"contribution")` | post-mortem of a run |
| examples | `examples/demo_gravity.m` (VCE on planted sigmas), `demo_fit.m`, `demo_lsqvce.m` (datum prior), `demo_sensors.m` (`Pool="obsType"`), `demo_workflow.m` | end-to-end stories on the v5 surface |
| the theses | `/mnt/project/25000_*.pdf` are PLAIN-TEXT dumps — read with `grep`/`awk`, not a PDF reader | equation-level comparison |

Diagnosis order when a number is wrong: configuration, then criterion, then
code. A cluster of unrelated failures sharing a symptom is session state.

---

## 7. Decisions in force — contradict, do not re-propose

| decision | what it fixed | contradicting it requires |
|---|---|---|
| R5-D4 | the biased raw-residual-ACV "separate" estimator was retired; the joint estimator is THE ONE noise estimator; acceptance is z-score recovery on O3, not equality to v4 | a fixture where the joint estimator is demonstrably worse than a correct separate one |
| R5-D5 | one accumulator; the flat family deleted; no aliases | — |
| H-15 / M3 3a | multi-axis streams are `separable`, whitened per axis, statistics added | a derivation showing the separable model is not what the data has |
| M3 3b | the pooled law lives in the PSD domain, positive by construction, sampled to any arc length | a case where the PSD-domain law is worse than an ACF law |
| M4 4d | "one re-read per pass" by residency (`CacheDays`), diagnostic sweeps on demand | a measurement |
| §12.2 (R5 confirms) | vTv NaN honesty; QR default with visible condition evidence; probe contract + seeds; the physical oracle; six-category taxonomy | — |
| R4 §11.2, R2 §8.4 | earlier "what is right" lists | — |

Recorded, not to be re-opened as findings: the joint estimator's damping is
"damage limitation" (WP-40 r6, mirror-measured); `Parallel.Axis="arcs"` is
gated on a u > 6000 measurement; O4 is open by lack of a file.

---

## 8. Leads noticed while preparing this brief (INFERRED, unverified — yours to decide)

These were read, not run. They are handed over as leads, ranked by the
severity they WOULD have if confirmed; none is a finding until you say so.

| # | lead | where | symptom to look for |
|---|---|---|---|
| L1 | the joint stop rule tests the PSD amplitudes only; levels, arc factors and prior sigmas can still be moving at `stopReason="converged"` | `adjustVCEjoint.m:429–431`; §4.2 | `R.history.sigma` per stream not settled on the pass the loop stopped; `R.noise.priorHistory` still moving |
| L2 | pure Gauss–Newton (`GaussNewton.Lambda0 = 0`) has no descent check — a step that raises F is accepted | `adjustNonlinear.m:102` (`opt.Lambda0 > 0 && …`) | F increasing pass over pass with `accepted = true` in `R.history` |
| L3 | a rejected damped step consumes one iteration of `MaxIter` | `adjustNonlinear.m:93–116` | `"maxPasses"` on a problem with several rejections early |
| L4 | the a-priori `C0 = [1; 0.3; 0.3; …]` is a fixed choice, not derived from the declared law; the fixed point reached may depend on it | `adjustVCEjoint.m:214–226`, help :44–49 | two runs from two non-degenerate starts landing on different `Cgroups` on the same data |
| L5 | the whiteness verdict uses K = 20 lags regardless of arc length and of `NumLags` | `adjustVCEjoint.m:756` (Ljung–Box) | p-values insensitive to correlation beyond lag 20 on 2160-epoch arcs |
| L6 | documentation drift into the theory: help blocks still name the deleted units `adjustVCE`, `adjust`, `robust`, `solveNonlinear` as the default / the reference (`fit.m:115–116, :184, :341`; `adjustVCEjoint.m:14–26, :136`; `Results.m:6, :23, :84, :97, :429–572`; `util/log.m:64, :82`; `util/verbosity.m:82`; `report.m:43`; `test_physOracle.m` See also) — a reader is told the separable mode is "the toolbox DEFAULT" | verified by `ls +lsqvce/+engine` (the files are gone); SWE lane for the fix, YOUR lane for whether the theory text around them is still the claim being tested |

---

## 9. Out of scope for this lane

Architecture, naming, layering, test structure, generated-page drift,
workflow hygiene (SWE reviewer). GRACE/GRACE-FO force models, frames, time
scales, the GRAVFIRE design function (the user's STEP 11; a separate
processing lane). Speed except where a numerical choice is the cause.

---

## 10. Report — format and where it goes

Write `dev/reviews/REVIEW_R6_LSQ_OPTIMIZATION_<yyyy-mm-dd>.md` in exactly
this shape. It is ATTACHED to the record (HANDOVER §13 points at it), never
paraphrased into it.

```markdown
# Review — domain: least-squares adjustment and optimization — main@<sha> — <date>

Scope reviewed: <files / drivers / kernels>, chosen because <reason>.
Tier: <MATLAB reachable? which suites executed? which oracles run?>

## Findings (most severe first)
### F1 — <one-line title> — REPRODUCED | INFERRED
What is wrong: ...
Why (the theory): ...
Observable symptom: ...
Evidence: <command run and its output, or file:lines read>

## Checked and sound
- <area>: <one line on what was verified and how>

## Not reviewed
- <what was out of scope or unreachable, and why>
```

Every finding that the user accepts becomes a work package in HANDOVER §13
with its own branch, oracle and acceptance criterion — that is the delegate's
job, not yours. Your job ends with the report.
