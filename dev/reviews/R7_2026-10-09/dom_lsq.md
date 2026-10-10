# Review — domain: linear model, solvers, robust and nonlinear drivers, diagnostics (dom_lsq) — main@b25bb11 (v6.1.3) — 2026-10-09

Scope reviewed: `+engine/BlockAccumulator.m` (accumulation, `preeliminate`, `solve` incl. the prior branch,
`backSubstitute`, `localCovariance`, `reducedNormalMatrix`), `QRSystem.m`, `NormalEqSystem.m`, `CommonPrior.m`,
`adjustBlock.m`, `adjustRobust.m` + `private/lawWhitenedResiduals.m`, `+kernel/robustWeights.m`,
`adjustNonlinear.m` + `batchSystem.m`, `applyRegularisation.m` (the `'lm'` branch), `residuals.m`,
`redundancy.m`, `diagnose.m`, `diagnoseAccum.m`, `leverage*.m`, `stat.m`, `ellipse.m`, `ftest.m`,
`systematic.m` (Ljung–Box/runs), `snoop.m` (w/tau formation), `fit.m` (dispatch, `priorsOf`, `priorUpdate`,
`runBatch`), `Results.m` (`compare`, `coefficients`, `toTable`, `report` forwarding), `report.m` (global test,
provenance). Chosen because the brief names them; drivers read end to end, kernels re-derived.
Tier: B. No MATLAB. Every REPRODUCED claim below is a line-by-line Python mirror (numpy 2.5 / scipy 1.18) in
`scratchpad/mirrors/dom_lsq/` (m1…m8), compared with a dense closed form, scipy, or Monte Carlo. No suite
executed. HANDOVER §0, §12.2, §12.6, §12.14, §12.29–§12.34 read; DEVELOPMENT §"QR is not optional here",
§"The decision is about cond(N)" read.

**Leads L2 / L3 (brief §8): both STILL HOLD at v6.1.3 — no HANDOVER round touched `adjustNonlinear.m`
after R5 M5 5b′ — and both are symptoms of a deeper defect (F2): the Levenberg damping never reaches the
system it is meant to damp.**

## Findings (most severe first)

### F1 — The QR backend squares the condition number on every accumulator solve: `blocks()` → `R1'R1` → normal-form Schur → `chol` in `merge` — REPRODUCED
What is wrong: the default backend is not a square-root solver on the path users run. Every system is
accumulated into a `QRSystem`, but `preeliminate` reads it back with `loc.blocks(ic, il)`
(`BlockAccumulator.m:1161`), which forms `N = R1.'*R1` (`QRSystem.m:153`); the Schur complement is done in
normal form; the reduced contribution is re-rooted by `[Rc, fail] = chol(red.N)` (`QRSystem.m:178`). This
happens also for `nl = 0` (kinds B/C, `batchSystem` sets `'NumLocal', 0`, `batchSystem.m:74`) and again in
the prior branch of `solve` (`commonBlocks` → `N = Ncc + Nk`, `BlockAccumulator.m:485`, `chol`).
Why: once `N = A'A` is formed, `cond(N) = cond(A)²` and no later factorisation recovers the lost digits —
DEVELOPMENT.md:460–473 says exactly this ("Forming N — even only to Schur-complement it — squares the
condition number … 'qr' eliminates on the observation equations"). The code does not eliminate on the
observation equations. The pinned precision tests measure `QRSystem` in isolation
(`test_precision.m:317` "Measured on the BACKEND IN ISOLATION, because in the full pipeline the effect is
masked") and attribute the masking to regularisation; the mirror shows it is structural (no
regularisation, `nl = 0`). The `eig` range-space branch in `merge` (`QRSystem.m:196`) then drops
eigenvalues below `max(size)*eps(max(dv))` of the SQUARED matrix, i.e. it refuses a full-rank problem at
`cond(A) ≳ 1e8` that the isolated QR solves to 1e-9.
Observable symptom: on the same data, `fit(..., Backend="qr")` loses as many digits as `"normal"` minus
about one order; at `cond(A) = 1e8` the common estimate is wrong by several formal sigmas; at
`cond(A) ≥ 10^8.5` a full-rank problem is refused (`lsqvce:qr:rankDeficient` / `singularR`).
Evidence (m1_qr_roundtrip.py, m1b, m1c; m = 800, nc = 8, singular values `logspace(0,-k,8)`):
```
 k  nl | accumulator qr (default) | QR isolated | normal        (rel. error of x_c)
 6   0 |  2.27e-06                | 1.14e-13    | 5.32e-05
 8   0 |  3.48e-02                | 1.59e-10    | 3.99e-02
 8   3 |  1.02e-01                | 4.77e-08    | 1.57e+00
max|dx|/sigma_x (error in formal-error units, 20 seeds):
 cond 1e7, noise 1e-10: accumulator median 1.38, max 14.9   | QR isolated 6.4e-07
 cond 1e8, noise 1e-4 : accumulator median 2.25, max 33.1   | QR isolated 3.2e-06
full-rank refusals, 200 seeds: cond 1e8.5 -> 42 refused; cond 1e9 -> 102 refused (QR isolated: 0)
```
Contradicts: DEVELOPMENT.md:460–473 and `QRSystem.m` help ("The condition number of R is cond(A), not
cond(A)^2"), as they apply to `BlockAccumulator`; the R5 §12.2 item "QR default with visible condition
evidence" rests on this property. Severity: **S1**.

### F2 — Gauss–Newton: the Levenberg damping is a no-op, so a rejected step is retried unchanged; `R.x` is then the REJECTED point — REPRODUCED
What is wrong: `adjustNonlinear` damps by building the accumulator with `'Reg','lm','RegLambda',lam`
(`adjustNonlinear.m:182`). That option is consumed only by `applyRegularisation` on the LOCAL block
(`BlockAccumulator.m:1200`), and the batch source has `nl = 0`, so `preeliminate` returns from its
`if nl == 0` branch (`BlockAccumulator.m:1163`) before regularisation is ever applied. The common block —
where all Gauss–Newton unknowns live — is never damped. Consequently: (a) after a rejection
(`adjustNonlinear.m:102–116`) the "re-linearise at x, more damping" re-solve returns the identical step,
which is rejected again; λ climbs ×10 per pass with no effect until `MaxIter` or λ > 1e12 (→ `diverged`,
whose message "damping reached … without finding a downhill step" is false — no damping was applied);
(b) the final "undamped" solve (`:156–158`) reports `sol.x = x + dsol.x`, i.e. the last accepted point
plus the full GN step that was just rejected; (c) the descent check is skipped on the first step
(`it > 1`, `:102`), so an uphill first step is accepted even with `Lambda0 > 0`.
Why: Levenberg–Marquardt guarantees descent only if `(N + λD)dx = n` is actually solved; here
`dx` is independent of λ. Leads L2 (no descent check with `Lambda0 = 0`) and L3 (rejections consume
`MaxIter`) hold as stated, but the real consequence is worse: with the DEFAULT `Lambda0 = 1e-6` the driver
stalls on problems pure GN solves (the mirror's `Lambda0 = 0` converges on both starts that stall).
Observable symptom: `R.history` shows identical `vPv` on every rejected pass; `R.x` has a larger
objective than the best accepted iterate; `stopReason "maxPasses"` (warning `adjustNonlinear:maxIter`) or
`"diverged"` after 18 rejections.
Evidence (m3_gauss_newton.py, m3b_stuck.py, m3c_reported.py, m3e_pureGN.py — exponential fit
`y = a·exp(b t)`, 30 points, budget `MaxIter = MaxPasses = 10`):
```
start [10,-5]: vPv 4.9, 36.8R, 36.8R, ... (9 identical rejections)  -> maxIter
   best accepted F = 4.901 ; REPORTED R.x = [2.868 -0.306], F(R.x) = 36.83 ; optimum F = 0.0357
start [0.1,1]: 7 rejections -> F(R.x) = 17.46 (best accepted 5.45)
same starts, Lambda0 = 0 (pure GN): tolF, F = 0.0357 on both
MaxIter 30: reason "diverged" after 21 passes, F(R.x) = 36.83
nl = 0: step with lam = 1e-6 and lam = 1e6 identical: True
Rosenbrock start [-1.2,1]: first step F 24.2 -> 2343 accepted (it > 1 guard)
```
Severity: **S2** (a warning fires, but the returned estimate is a known-worse point and the diagnosis it
prints is wrong).

### F3 — Nonlinear problems with `data.priors`: the prior is imposed on the Gauss–Newton INCREMENT, not on x — REPRODUCED
What is wrong: kind C with `x0` and `priors` is accepted (`fit.m:238`, `runBatch` → `adjustNonlinear(...,
'Prior', pri)`). Each pass solves for the step `dx` with `CommonPrior.normals` adding `K'x0/σ²` to the
right-hand side (`CommonPrior.m:123–130`), i.e. the constraint is `K·dx ~ x0` every pass. The fixed point
`dx = 0` therefore satisfies `J'w + K'x0/σ² = 0`, not the stationarity of the declared objective
`‖w‖² + ‖(Kx − x0)/σ‖²`. With `x0 ≠ 0` the iterate is pushed by x0 on every pass; with `x0 = 0` the prior
is satisfied by the increment and ignored by the estimate.
Why: in a linearised problem a prior on x must enter as `K·dx ~ x0 − K·x_current`. The final solve
(`adjustNonlinear.m:156–158`) repeats the same error.
Observable symptom: `R.x` violates the declared prior (tight σ) or sits away from the regularised optimum
(loose σ), sometimes with `stopReason "converged"`.
Evidence (m3d_prior_on_step.py, the `test_fit` exponential-decay fixture, prior on b):
```
prior b ~ +0.30 (sigma 1e-6): GN reports b = +0.599984 ; declared optimum b = +0.299994
prior b ~ -0.20 (sigma 1e-6): GN reports b = -0.399997 ; declared optimum b = -0.200000
prior b ~ +0.30 (sigma 1e-3): GN reports b = -0.178378 (tolX = "converged") ; optimum b = -0.166254
prior b ~  0.00 (sigma 1e-3): GN reports b = -0.200093 (= the UNCONSTRAINED optimum) ; optimum -0.186610
```
No test combines `x0` with `priors` (`grep "'x0', \[" dev/tests/*.m`: two fits, neither has priors).
Severity: **S1**.

### F4 — A datum prior's covariance depends on the arbitrary datum sigma: "exact for any weight" holds for x, not for Cxx — REPRODUCED
What is wrong: `fit.m:398–399` / Manual 03_data.md:25–26 state the datum `L = B'` is "the exact
minimum-constraint solution for ANY weight". The estimate is (B'x = 0 to 1e-16); the covariance is not:
`BlockAccumulator.solve` returns `Cxx = s0sq * inv(A'A + BB'/σ²)` (`:526`), which equals the
minimum-constraint cofactor plus `σ²·B(B'B)⁻²B'` scaled by `s0sq`. `diagnose` would subtract a datum
correction through `sol.QxxCorr` (`diagnose.m:262`), but nothing in the tree sets that field
(`grep QxxCorr` → only `diagnose.m`). `R.coefficients()`, `toTable`, `D.Qxx`, `ellipse` all read it.
Why: under minimum constraints `B'x̂ ≡ 0`, so its variance is 0; the code assigns it `s0sq·σ²/(B'B)`.
Observable symptom: on `testDatumIsADeclaredPrior`'s levelling network, the standard errors change by
three orders of magnitude when the datum sigma goes 1 → 1e3, while x does not move.
Evidence (m4_common_prior_paths.py, the test_fit fixture):
```
sigma_prior=1   : B'x=2.2e-16  se(x) reported = [0.00495 0.00607 0.00495 0.00607]
sigma_prior=1000: B'x=7.7e-11  se(x) reported = [2.47742 2.47742 2.47742 2.47742]
minimum-constraint se (s0sq*pinv(N))          = [0.00429 0.00554 0.00429 0.00554]
```
Contradicts the recorded claim (R5 M5 5b, §12.12 D8 / fit help) as far as any covariance-derived
number is concerned. Severity: **S1**.

### F5 — Common priors are outside `residuals()` and `redundancy()`: wrong dof, a misfiring `sigma0Mismatch`, a false `identityViolated` — REPRODUCED (arithmetic) / INFERRED (that the warnings fire in MATLAB)
What is wrong: `data.priors` (common) reaches the solve (`solve` prior branch: `vtv` includes the prior
residual, `dof = m + obj.prior.nk - obj.nc`, `BlockAccumulator.m:503, 522`) and `diagnoseAccum` (prior
appended as a batch, `diagnoseAccum.m:130`), but `residuals()` sweeps only system arcs
(`residuals.m:210`, `res.dof = mTot - nParTot`, `:275`) and `redundancy()` checks
`sum(r_i)` over data arcs against `mTot - nParTot` (`redundancy.m:360`). The local-prior streams
(`prior_<name>`) are in the sweep, so this is a COMMON-prior-only second path.
Why: one estimator, two `s0`: the solve counts pseudo-observations in `m` and their residual in `vtv`,
the sweep counts neither. For a datum the sweep's `dof = m − u` also ignores the rank defect
(correct is `m − (u − d)`). The `sigma0Mismatch` text blames the loader ("did NOT return the same data on
the second pass"); the record sets `supersedesSol = true`, so the WRONG dof is presented as the
authority.
Observable symptom: `R.residuals()`, `R.toTable`, `R.plot("residuals")`, and every `adjustRobust` pass
on a fit with `data.priors` warn `lsqvce:residuals:sigma0Mismatch`; `acc.redundancy` /
`vce` / `contribution` warn `lsqvce:redundancy:identityViolated`.
Evidence (m4_common_prior_paths.py):
```
datum (test_fit fixture): sol s0 = 0.009910 (dof 2)  residuals() s0 = 0.014014 (dof 1)  rel 4.1e-01 (guard 1e-6)
ridge prior, m=40,u=6:   sigma=10 rel 4.0e-02 ; sigma=1 rel 6.2e-01 ; sigma=0.1 rel 9.9e-02  -> mismatch every time
redundancy(): sum r(data) = 34.002 / 34.177 / 38.430 vs m-u = 34 -> relChk 5e-5 / 5e-3 / 1.3e-1 (guard 1e-6)
```
Contradicts the D8-2 principle (§12.12: a prior "enters accumulation, elimination, the sweeps,
residuals() and the redundancy count with NO second path") for common priors. Severity: **S2**.

### F6 — A free network without a datum is refused or silently solved depending on rounding inside `chol` — REPRODUCED
What is wrong: the pinned refusal (`test_fit.m:931`, `lsqvce:qr:singularR` for the levelling network
without datum) depends on whether `chol(red.N)` (`QRSystem.m:178`) of the numerically singular `N`
succeeds. When it fails, the `eig` branch drops the null direction and `solve` refuses (`:282`). When
rounding leaves `N` barely positive, `chol` succeeds, the folded factor's smallest pivot is ~1e-8
relative (far above `rtol ≈ 1e-15`), and `x` is returned with an arbitrary level and cofactor sd ~1e7.
The NEQ backend has the same coin flip.
Why: the H-16 lesson ("the refusal sat AT the roundoff boundary") moved from the triangular solve into
`chol`; the rank decision is LAPACK's, not the code's.
Observable symptom: random connected levelling networks without datum: QR refuses 61 %, SOLVES 39 %.
Evidence (m2_rankdef_refusal.py, m2b_detail.py; 2000 random networks):
```
merge eig -> singularR: 1225/2000 ; merge chol -> qr solve ok: 775/2000
u=8: chol path, foldPivotRel 8.8e-09 vs tolRel 1.5e-15, max cofactor sd 4.8e+07, |x - x_minnorm| = 16.6
test_fit fixture: eig path, refused (passes by the fixture's luck)
```
Severity: **S2**.

### F7 — IRLS stop rule is not unit-invariant: `max|dx| < Tol·max(max|x|, 1)` stops after 2 passes for parameters of order 1e-9 — REPRODUCED
What is wrong: `adjustRobust.m:159` uses an absolute floor of 1 on `max|x|`. For dimensionless SH
coefficients (~1e-9…1e-11, the toolbox's own domain) `Tol·1 = 1e-8` exceeds every update, so the loop
stops at pass 2 with `converged = true`.
Why: the criterion must be relative to the parameter scale or to the formal errors (as
`adjustNonlinear`'s `dxRel` already is, `adjustNonlinear.m:119–121`).
Observable symptom: the same IRLS problem expressed in rescaled units stops at pass 14 / 9 / 5 / 2 and
returns a point up to ~1 formal sigma from the fixed point while reporting "converged".
Evidence (m8_irls_units_cov.py, IGG3, m = 400, 10 % blunders):
```
unit 1     : passes 14, |x - fixedpoint|/sigma = 4.1e-08
unit 1e-6  : passes  5, 7.7e-03
unit 1e-9  : passes  2, 1.08e+00
```
Severity: **S2**.

### F8 — The IRLS result's sigma0 and covariance are reported as if least squares; `report()` does not say so — REPRODUCED (numbers) / INFERRED (report wording)
What is wrong: `adjustRobust` help says "its covariance is that of the reweighted problem and is scaled
by the reweighted sigma0 (report() says so)" (`adjustRobust.m:16–18`). `report.m` has no mention of
robust/reweighted/influence anywhere (`grep -i "reweight\|influence\|irls" report.m` → nothing); its
provenance line prints only version/release/date (`report.m:113–117`). The numbers it judges are:
`sigma0 = sqrt(Σ f_i v_i² / (m − u))` with zero-factor rows (floored 1e-12) still counted in `m`, and
`Cxx = s0²(A'FA)⁻¹`.
Why: neither is the scale nor the sampling covariance of an M-estimator (Huber's sandwich). Under Huber
the retained blunders inflate s0; under redescending influences the zero rows deflate it.
Observable symptom: global test, sigma0 CI and coefficient bounds of a robust fit are wrong in both
directions without a word in the verdict.
Evidence (m5_irls_stats.py, m8_irls_units_cov.py; true sigma 1):
```
R.sigma0: clean huber 0.898 / igg3 0.832 ; 10% 25-sigma blunders: huber 2.175, igg3 0.844 (0.890 with dof = alive-u)
reported se / Monte-Carlo sd: clean huber 0.86-0.92, clean igg3 0.73-0.75, 10% huber 1.80-1.86, 10% igg3 0.76-0.79
```
Severity: **S2**.

### F9 — `Results.compare` uses `df2 = n − numel(R.x)`: eliminated locals are not counted — REPRODUCED
What is wrong: `Results.m:432` builds `u = numel(Rk.x)` (commons only) and `rank = u`; `ftest` then takes
`df2 = n − uF` (`ftest.m`). On the streaming kind the locals are parameters too
(DEVELOPMENT.md "The redundancy … counting the eliminated parameters"); `Rk.diag.u` already holds the
full count.
Why: `s0sqF = vPvF/df2` is underestimated, so F is too large.
Observable symptom: false-positive "SIGNIFICANT" verdicts under H0.
Evidence (m5_irls_stats.py (c), n = 1000, nc = 6, nl = 200, 2000 trials, alpha 0.05):
```
rejection rate under H0: code (df2 = n - nc) 0.074 ; correct (df2 = n - nc - nl) 0.049
```
Severity: **S2**.

### F10 — A-posteriori error ellipse uses `chi2inv(conf, 2)` with the estimated `s0²`: under-covers at small redundancy — REPRODUCED
What is wrong: `ellipse.m:102, 121` scale by `D.s0sq` (default `'aposteriori'`) and `k2 = chi2Inv(conf, 2)`.
With an estimated variance factor the exact region is `2·F(2, r)`.
Observable symptom: a "95 %" ellipse covers 86 % at r = 5, 90.5 % at r = 10.
Evidence (m7_ellipse.py, 20 000 MC draws):
```
r =   5: coverage chi2 0.862 ; 2F(2,r) 0.952 ; axis ratio 1.390
r =  10: 0.905 / 0.950 ; r = 30: 0.933 / 0.949 ; r = 200: 0.946 / 0.948
```
The constant itself is right (`sqrt(chi2inv(0.95,2)) = 2.4477`). Matters for kinds B/C networks, not for
streamed GRACE redundancies. Severity: **S3**.

### F11 — `localCovariance` (verbose local summary) scales by `s0sq_hint = ltl_red/dof_red`, which is not a variance factor — REPRODUCED
What is wrong: `BlockAccumulator.m:646` `Cyy = s.s0sq_hint * Qll + G*Cxx*G'`, with
`s0sq_hint = max(ltl_red,0)/max(dof_red,1)` (`:1238, :1280`). `ltl_red` is `l'l` after eliminating the
locals only (commons not solved), so it carries the whole common signal.
Observable symptom: the "locals[...] mean sd" log line of `logLocalSummary` is inflated ~300×. The
user-facing route (`residuals.m`, `Cyy = sol.s0sq*Qll + G*sol.Cxx*G'`) is correct to 5.8e-11.
Evidence (m6_schur_cyy.py): `s0sq_hint = 1.184e+05 vs s0sq 1.116 -> Cyy(log) sd x324.8`.
Severity: **S3**.

## Checked and sound
- Schur elimination, back-substitution and dof with locals (NormalEq arithmetic): x 8.9e-16, s0² and dof (538) identical to the dense block-arrow solve, Cxx 5.8e-11, y 3e-15, Cyy = s0²Nll⁻¹ + G Cxx G' 5.8e-11; the cross-covariance Cxy = −Cxx G' matches the dense solve (not reported by the code). REPRODUCED (m6).
- Accumulated vᵀPv = lᵀPl − nᵀx̂ vs direct sum on a well-conditioned problem: agree to 5e-11 relative; the NaN-honest guards (8 digits NEQ / 15 merge / 24 rotation) only READ, judged reasonable. Partly REPRODUCED (m6), guards INFERRED.
- Prior branch of `solve`: `dof = m + nk − nc` (pseudo-observations counted), `Cxx = s0² N_tot⁻¹`, Koch–Kusche `r_k = p − tr(L N⁻¹ L')/σ²` with `N⁻¹ = Cxx/s0sq` (`fit.m:529`): correct Helmert form. INFERRED.
- `robustWeights`: Huber c = 1.345, bisquare c = 4.685 (95 % efficiency), IGG3 (k0, k1) = (1.5, 3.0) with Yang's weight `(k0/|u|)((k1−|u|)/(k1−k0))²`, MAD·1.4826, std→1 fallback: match the literature. The `vw/s` vs `(vw−med)/s` standardisation is the recorded §12.34 fact; not re-raised. INFERRED + mirrored.
- Pope τ critical value `t·sqrt(r)/sqrt(r−1+t²)`, `t ~ t(r−1)`: 1.91545 vs Monte-Carlo 1.91370 (r = 12). Baarda δ0(0.001, 0.2) = 4.1321, δ0(0.05, 0.2) = 2.8016. `chi2Inv` vs scipy max rel err 2e-15 over r = 1…1e7. `fInv` = scipy to 6 digits. REPRODUCED (m5, m6).
- `diagnose`: w-test uses the a-priori σ0 = 1 (`diagnoseAccum` sets `'sigma0', 1`) with `normCrit`, Pope uses the global a-posteriori s0 with `tauCrit(·, r)`; r_i formula for correlated batches `1 − diag(A Qxx A' P)`; Bonferroni α/n. INFERRED.
- `report` level 1: χ² band `chi2Inv(0.005|0.995, r)/r` on `s0²/σ0²` and the σ0 CI `s0·sqrt(r/χ²)`. INFERRED.
- `systematic` Ljung–Box `Q = n(n+2)Σρ²/(n−k)`, df = h, flagged unreliable for u/n > 0.05; runs test moments. INFERRED.
- `leverage` / `leverageBlocks`: h = a'Q_CC a + b'(Nll⁻¹ + G S G')b − 2a'S G'b with the regularised Nll, the cross term included. INFERRED.
- `adjustNonlinear` TolX in formal-error units (unit-invariant); final covariance is that of the linearised problem at x (true only when no rejection preceded — see F2). INFERRED.

## Not reviewed
- The joint VCE estimator (`adjustVCEjoint`, `psdVCE`, `arcVarianceFactors`, Ljung–Box K = 20 in the joint sweep, lead L5) and its robust pass — other lanes; C1/C2 of the coordinator not extended.
- `snoop` re-solve ladder end to end (only the w/τ formation read); `engine.datum` (no driver caller found, only its own help example); `checkPSD`; `outliers`; the Hutchinson trace path in `redundancy`.
- The physical oracle O3 (`demo_gravity_design`) — needs MATLAB; F1's synthetic fixture stands in for it. A tier-A check: `R.x` of `fit(..., Backend="qr")` vs a dense `[A B] \ l` on the polar-gap design, in formal-error units.
