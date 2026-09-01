# lsqvce — DEVELOPMENT record and reference

*Split out of README.md by L8 / WP-70 (2026-09-02): the README now carries only
the introduction, the quick start and the usage guidelines; everything below is
the development record — status, the research narrative behind the diagnostics,
the module census, correctness history, design decisions, validation, and the
measured chapters. `dev/HANDOVER.md` remains the single status-and-plans
document; this file is the technical companion the README used to be.*

> **Status (2026-09-01, v2.0 — see RELEASE_NOTES.md).** Green everywhere: **23 suite files, 634 test
> points, 0 fail** — verified locally in MATLAB R2026a (including the
> pool suites with real `parfor`/`spmd` workers and the machine-local
> speed budgets) and on every push by GitHub Actions CI (a static `audit`
> job plus the full MATLAB suite on hosted R2026a runners; `main` is
> branch-protected on both). Run the suite with `dev/tests/run_all_lsqvce.m` — not
> `runtests` on the folder. Two review tracks are closed since v1.0: the blind
> five-expert code review (findings F-01…F-15, all resolved — including
> two CRITICAL correctness fixes, see *Correctness history* below) and a
> test-taxonomy audit (six declared test categories, per-unit coverage
> counted mechanically by the audit on every run). `dev/HANDOVER.md` is the
> single status-and-plans document; MIT licensed.
# Diagnostics reference (was lsqdiag)

Solving and **diagnosing** large least-squares problems where `A` and `P` can
only be built in batches.

Streaming: nothing of size `n` is ever formed. `O(u²·nk)` time and `O(u² + nk)`
memory per batch.

---

## The headline finding: **Pope's tau goes blind under contamination**

Discovered while testing v1.0, and the most consequential result in the toolbox.

Pope's studentised residual divides by `s0` — which is estimated **from the data
that contain the blunders**. Several blunders inflate `s0` together, and since
each is then divided by that inflated `s0`, **they hide one another**.

Measured (n=200, u=8, 20σ blunders; fraction of true blunders detected):

| contamination | Pope | **MAD (robust)** |
|---|---|---|
| 2.5 % | 100 % | 100 % |
| **5 %** | **23 %** | **100 %** |
| **10 %** | **0 %** | **100 %** |
| 20 % | 0 % | 99 % |
| 30 % | 0 % | 75 % |

At 10 % contamination **a 100σ blunder is undetectable**: `max|w|` saturates at
~3.3 while the critical value sits at 3.9. Tau is bounded by `√r`; it cannot
grow without limit however gross the error.

**Iterating does not save you.** If nothing trips the threshold on the first
pass, the rejection loop exits immediately having rejected nothing — and reports
*"converged, 0 rejections"*, which reads as clean. Verified: `snoop` with Pope
rejects **zero** of five 25σ blunders.

**The fix**, and it is nearly free:

- `'test','mad'` uses a MAD-based scale (high breakdown point). Recovers all
  five. Costs **zero** extra false alarms on clean data (0 in 300 realisations),
  is unbiased (mean 1.0025), and is only 1.65× less efficient.
- The toolbox **always computes both** and reports `s0/MAD`. That ratio is a
  direct contamination symptom: it fired 100 % of the time Pope had collapsed,
  and 0/300 times on clean data.

```
[ERROR] CONTAMINATION DETECTED. s0 = 7.673 but the robust (MAD) scale = 4.472...
        At this level of contamination the test can be completely blind --
        a clean flag list is NOT evidence of clean data.
```

---

## The cure has its own cliff — and it hides it

The MAD fix (above) is not the end of the story. Pushed further, **every
scale-based indicator switches itself off at exactly the moment the data are
worst.** Measured:

| contamination | 0% | 10% | 30% | **50%** | 60% |
|---|---|---|---|---|---|
| s0/MAD ratio | 1.02 | 4.26 | 3.58 | **0.95** | 0.61 |
| fhat (est. contamination) | 0.02 | 0.11 | 0.30 | **0.00** | 0.00 |
| **s0/σ₀ (a-priori backstop)** | 1.00 | 6.40 | 10.99 | **14.13** | 15.52 |

At 50 % contamination the `s0/MAD` ratio reads **0.95 — below its clean-data
value of 1.02**. An alarm that says *"cleaner than clean"* at total breakdown.
The reason: past ~40 % the MAD adopts the **blunders** as its reference
population and the good data become "the outliers."

So v3 reports **three** indicators, not one:

1. **`fhat`** — the *gauge*. Estimated contamination fraction, accurate to
   ±0.015 below 30 %. Tells you *how much* is bad, which is what determines
   what to do.
2. **`s0/σ₀`** — the *backstop*. **Monotone to 60 % and unbounded**, because σ₀
   is a priori and contamination cannot touch it. The one indicator that cannot
   be turned off.
3. **A counting-based refusal** above 25 %. This depends on *no estimator at
   all*: "you are asking me to delete a third of your data" cannot be falsified.

**Escalation ladder** — `clean → moderate → severe → refuse → breakdown`:

```
[ERROR] TOTAL BREAKDOWN OF THE OUTLIER MACHINERY.
        The robust scale reports almost no outliers (fhat = 0.000) while the
        a-priori backstop says the residuals are 14.1x larger than sigma0
        promised. These are only compatible if MORE THAN HALF the data are bad.
        DO NOT REJECT ANYTHING.
```

`snoop` now **refuses to run** above 25 % contamination (override with
`'force',true`, which still warns). `varComp` refuses too — blunders inflate a
group's variance, which down-weights the group, which hides the blunders
further: the same self-masking pathology one level up.

**The threshold is calibrated, not guessed.** `fhat` has a floor on clean data
(`P(|Z|>2.5) = 0.0124`) and its sd is **1.30× binomial** — the residuals satisfy
`u` constraints and share a MAD estimated from the same sample. A naive 0.02
threshold escalated **87/300 clean datasets**. The calibrated line
(`floor + 4·inflated_sd`, n-dependent) escalates **0/300**, while still detecting
5 % contamination — where Pope starts failing — **100 %** of the time.

---

## |w| and r mean nothing alone

A residual reveals only the **fraction `r_i`** of any error present:

```
E[Δv_i] = −r_i · b_i        (verified exactly)
```

So on a levelling net with a dangling spur, a **40σ blunder at `r = 0` produces
`|w| = 0.000`** — a beautiful small residual on a grossly wrong observation. v3
reports the **four quadrants**:

| | meaning | action |
|---|---|---|
| `\|w\|` hi, `r` hi | clean outlier | reject with confidence |
| `\|w\|` hi, `r` **lo** | **ambiguous** | coin flip — re-measure, don't reject |
| `\|w\|` lo, `r` hi | verified | tested, and passed |
| `\|w\|` lo, `r` **lo** | **UNTESTED** | *not* "fine" — the test had no power |

Plus the **hidden-blunder envelope**: not just the worst MDB, but its
distribution — *"half your observations could harbour an undetected error larger
than 5.5σ"* — and an explicit count of observations with **infinite** MDB, which
drop out of any percentile and make a finite-only summary look reassuring
precisely when it should not.

---

## Modules

Tiers since 2.0 (`dev/docs_src/tiers.csv`, audit check 18). Every name below is the
qualified one you type.

### USER tier — the four names

| Function | Purpose |
|---|---|
| `lsqvce.fit` | the one verb: `data` struct + `Options` → `Results`; dispatches to the engines by the fields present and `Options.Model` |
| `lsqvce.Options` | every knob, typed, defaulted to the engines' own defaults; `disp(opts)` marks the defaults; `engineArgs(target)` is the mapping |
| `lsqvce.Results` | solution + diagnostics + provenance in one object: `R.x`, `R.sigma0`, `R.coefficients()`, `R.report()`, `R.history`, `R.noise`, `R.acc` |
| `lsqvce.report` | the five-level verdict (`Results.report` forwards here) |

### ENGINE tier — what `fit` dispatches to

| Function | Purpose |
|---|---|
| `lsqvce.engine.adjust` | one-call driver: accumulate → solve → diagnose → judge |
| `lsqvce.engine.adjustBlock` | one-call **block-arrow** front door → `Results`; the base of the separable/VCE engine. Pairs with `adjust` (flat) and `adjustVCE` (iterative). _(Was `recover`; that name was REMOVED on 2026-07-18 — all aliases were dropped pre-release.)_ |
| `lsqvce.engine.makeLoader` | turn a `buildDay(day,obsType)` into an accumulator-ready `Loader` (meta / per-arc / one-day cache) |
| `lsqvce.engine.adjustVCE` | iterative separable-GLS engine — **SHAPE** (temporal AR whitening) + **SCALE** (VCE), configured on the constructor (`Shape`/`Scale`/`Outer`/`Joint`/`Axes` groups, constructor-only since 2.0) over `adjustBlock`. Pairs by name with the in-memory shortcut `adjust`. `P.empiricalStochastic()` runs **empirical stochastic modelling** end to end (homogenise → loop → validate); `P.validate(R)` is the formal-vs-empirical **acceptance test**. _(Was `SeparableProblem`; that name was REMOVED on 2026-07-18 — all aliases were dropped pre-release.)_ |
| `lsqvce.engine.adjustVCEjoint` | Ellmer-style JOINT spectral VCE: declared per-stream PSD co-estimated with per-arc variance factors over shared Rademacher probes. **Production guidance (measured, F-15):** run it OPEN LOOP or with at most 3 damped passes (`'StepSchedule','decay','AvgTail',4`) — the undamped closed loop is bistable and never beat open-loop pass-1 quality at any tested volume |
| `lsqvce.engine.Accumulator` | streaming NE/QR; correlated batches; free-network; sparse |
| `lsqvce.engine.diagnose` | per-observation `v`, `r_i`, `w`, MDB, external reliability |
| `lsqvce.engine.systematic` | **Ljung–Box, runs test, covariate binning** — the blind spot |
| `lsqvce.engine.snoop` | iterative rejection via Cholesky **downdate** |
| `lsqvce.engine.robust` | Huber / IGG3 reweighting |
| `lsqvce.engine.varComp` | Helmert variance components, with a **resolution limit** |
| `lsqvce.engine.traceEstimate` | Monte-Carlo (Hutchinson) trace estimator, **Ellmer 2018 §6.5.6** — estimates `tr(N⁻¹Nₘ)` without forming `N⁻¹`. **Opt-in** VCE accelerator (`'TraceMethod','hutchinson'` in `redundancy`/`vce`); exact leverage stays the default |
| `lsqvce.engine.datum` | null-space basis; free-network / minimum-constraint solutions |
| `lsqvce.engine.solveNonlinear` | re-linearisation loop with convergence history |
| `lsqvce.engine.ftest` | nested-model F-test |
| `lsqvce.engine.ellipse` | confidence ellipses |
| `lsqvce.engine.stat` | distribution kernel (no Statistics Toolbox) |
| `lsqvce.engine.selectArAcv` | automatic AR order from a sample ACV (no grid); ACV-domain sibling of `selectArOrder` |
| `lsqvce.engine.resolveSatSpecs`, `lsqvce.engine.satKeyOf` | per-satellite config resolution and the shared satellite-key convention |

### KERNEL tier — numerical building blocks

| Function | Purpose |
|---|---|
| `lsqvce.kernel.*` | the 15 transforms, whitening filters and statistics (`dctI`, `psd2acv`, `levinsonDurbin`, `armaWhiten`, `varWhiten`, `separableWhiten`, `toeplitzCofactor`, `decorrFactor`, `residualAcv`, `whitenessCorr`, …) — listed with one line each in `dev/docs_src/tiers.csv` and `help lsqvce.kernel` |
| `lsqvce.util.log` | toolbox-wide logging: set `logfid` (file id) + `logflag` (screen) once, and **all** status text and warnings route through it |

### Examples and the template (not package code)

| File | Purpose |
|---|---|
| `examples/demo_fit` | *(example)* **START HERE for 2.0**: the user path on all three data kinds (in-memory, robust, batched, nonlinear, streaming gls + vce), the gravity oracle as a `data` struct |
| `template/lsqvce_template.m` | *(template, not package)* ENGINE-tier template for a Loader-based case (procOpt-driven); Returns the Loaders; you fill in two documented stubs (`myBuildDay`, `myEpochCov`). The `observationGroups` loop is real runnable code; the stubs error until written. |
| `examples/demo_grace_dataLoader` | *(example, not package)* the GRACE-FO reference setup — C/D pair merged into one system (range-rate couples both satellites' arc states) plus independent satellite E on the same field. A WORKED case to read alongside the template. |
| `examples/demo_grace` | *(example, not package)* the driver that runs it end to end |
| `examples/demo_gravity_design` | *(example, not package)* fully-normalised Legendre functions and spherical-harmonic potential partials along a sampled orbit — the ONE physically structured design generator; `demo_gravity` and `test_physOracle` both call it, so the recursion and its accuracy claim exist once |
| `examples/demo_gravity_dataLoader` | *(example, not package)* a GRACE-like fixture built on that design: GRACE-C, a trailing GRACE-D, and the inter-satellite potential difference read as an energy-balance KBR stream — with a **different number of local unknowns per system** (2/3/1 per arc), a **ragged arc schedule** (arcs missing by ID, not just by count) and a **different length for every arc** |
| `examples/demo_gravity` | *(example, not package)* the driver: what VCE recovers of the planted stream weights under real column correlation, coefficient recovery in units of the formal error, and error degree variances |

---

## Correctness history — what the reviews changed

The v1.0 claims below survived two adversarial review tracks; two did not
survive unchanged, and honesty requires the pointers. A blind five-expert
review (2026-08) produced findings F-01…F-15; the two CRITICAL ones were
real wrong-number defects, both fixed and regression-locked:

- **F-01** — `vce` with `'TraceMethod','hutchinson'` used the wrong
  denominator (cc-attribution sum instead of the redundancy
  `r_g = n_g − local_g − tr(N⁻¹N_g^cc)`); validated against dense joint
  hat-sums to 1e-8 after the fix.
- **F-03** — the joint spectral paths omitted `B·y` from the residual and
  used a common-only projector; fixed via a shared `jointArcParts`, with
  dense parity in the whitened metric.

F-15 (session finding) established that the undamped closed-loop spectral
alternation is bistable — the production guidance in the `adjustVCEjoint`
row above is the measured consequence. The full evidence trail lives in
`HANDOVER.md` (§8 findings, §9 work packages, §10 test-taxonomy review).

---

## Design decisions, and the evidence

**`vPv = l'Pl − b'x` is a trap.** At σ=1e−6 it is **89 % wrong**; at 1e−9, wrong
by 10⁶. QR is the default and takes `vPv` from the discarded rows; `diagnose`
recomputes it honestly as `sum(p.*v.^2)`.

**Bonferroni is mandatory.** Uncorrected α=0.001 at n=1e7 flags **10 000 good
observations by chance**. Monte-Carlo: 84 false flags → 0.

**`qv` goes negative.** A parameter seen by one observation gives `qv = −2.2e−16`
→ `w = NaN`. Clamped, and the observation is *reported as uncontrolled* rather
than looking clean.

**`chi2Inv` needs log-space Newton.** The χ² density diverges at the origin for
`r<2`, stalling Newton-in-`x` (4e−8 error at r=1). In `log x` it reaches machine
precision.

**Iterating helps via *masking*, not false positives.** Measured: false-positive
rates are *comparable*. The real gain is detection — 8 blunders in n=200:
iterative **100 %**, single-pass **74 %**.

**No fill-reducing permutation.** Measured benefit: **0.0 %** on banded `N`
(already optimal), **3.7 %** on irregular. Against that, a permutation must be
honoured in *three* places, and forgetting one gives plausible numbers on the
**wrong parameters** (measured `qx` error 0.90). Sparse `chol` without
permutation is **bit-identical** to dense. For arrow structures (global params +
per-arc locals) you want Schur elimination, not this toolbox.

**Residuals are datum-invariant.** Verified to 1e-16: `v` and `r_i` are identical
whichever minimum constraint you pick. The datum is a choice; the fit is not.

---

## Validation

**135 checks across nine suites**, against dense formulae, sampling theory, and
brute-force simulation.

**External reliability is exact**, not approximate:

```
predicted shift = 1.765467 sigma      <- from the MDB formula
realised  shift = 1.765467 sigma      <- inject the MDB, re-solve, measure
```

MATLAB unit tests: run `run_all_lsqvce` for the full suite (23 files, 634
test points, ~3 min locally; the speed suite is reported separately by
design). `run_lsqvce('all')` delegates to it, so it is the same run;
`run_lsqvce('<target>')` runs one suite or a name pattern while iterating.
Every test point declares its category (`contract`, `reference`,
`precision`, `speed`, `robustness`, `vectorisation`, plus the cross-cutting
`diagnostics`/`provenance` tags); the audit derives the unit x category
coverage matrix from the test bodies on every run, every tolerance carries
a comment naming its origin, and a public unit without contract coverage
fails the gate.

---

## What no automated test can do for you

**Systematic errors** produce inflated χ² *and* a clean snooping result. That
combination is the signature. `lsqvce.engine.systematic(D, 'covariates', ...)` now
automates the search — Ljung–Box, runs test, and binned means against any
covariate you supply (time, elevation, azimuth, temperature, satellite ID). It
turns *"go look at a plot"* into *"bin 7 of elevation has mean w = 2.9, p = 1e-4"*.

**Uncontrolled observations** (`r_i ≈ 0`) have their residual forced to zero. A
blunder there is invisible. A clean snooping result does not exonerate them —
*it cannot see them*.

At large `n` the χ² band collapses like `1/√r` (±0.4 % at r=1e6), so it rejects
on any imperfection however harmless. **Read the ratio, not the verdict.**

---

## Options

```matlab
lsqvce.engine.adjust(u, src, nb, ...
    'method',     'qr', ...      % 'qr' (default) | 'normal'
    'sparse',     false, ...     % keep N sparse (normal method only)
    'datumB',     B, ...         % free-network / minimum constraint
    'test',       'pope', ...    % 'pope' | 'baarda' | 'mad'  <- see above
    'bonferroni', true, ...
    'robustCheck', true, ...     % always cross-check s0 against MAD
    'store',      'flagged', ... % for huge n
    'progress',   true);         % ETA for long runs
```

### How much the run tells you: `verbose`

`verbose` takes a LEVEL, not a flag. Each level includes the ones below it.

| level | what you get |
|-------|--------------|
| `0` | nothing |
| `1` | outer passes, each stage's start/done with timing, and the convergence report |
| `2` | everything else: per-day accumulation and residuals, per-stream AR order and whiteness, per-group shares and sigmas, system conditioning, inner VCE iterations |


```matlab
R = lsqvce.engine.adjustVCE(L, days, 'verbose', 2, 'GroupBy', 'obsType').run();
% through the user front door: lsqvce.fit(data, lsqvce.Options(Model="vce", Verbose=2))
```

A LOGICAL still works and keeps its old meaning (`false` → 0, `true` → 1), so
existing code is unaffected; levels 2 and 3 are new capability you have to ask
for. At level 1 a pass looks like this — note the convergence report is PER
STREAM, because `max|sigma-1|` alone tells you how far the worst stream is from
unity but not *which* one:

```
===== Iteration 1 of 6 =====================================
      variance components -> adjustment -> stochastic model -> convergence

Variance components  estimating between groups ...
Variance components  done: 4 group(s), max|sigma-1| = 3.104e-02   1.83 s
Adjustment           block-arrow, current weights ...
Adjustment           done   0.94 s
Stochastic model     re-estimating temporal correlation per stream ...

Convergence
  "posC" of "demo-CD-pair"              sigma =  1.0342   white = yes  (p = 0.310)
  "posD" of "demo-CD-pair"              sigma =  0.9871   white = yes  (p = 0.440)
  "rr" of "demo-CD-pair"                sigma =  1.2044   white = NO   (p = 0.002)
  max|sigma-1| = 2.044e-01 (tol 5.0e-03) -> another pass
```

Each row carries BOTH halves of the stopping rule — the variance factor and the
whiteness of the residuals — so a row that fails one but not the other tells you
immediately which half is stuck.

**To see this for yourself rather than take the sample above on trust:**

```matlab
demo_messaging                      % the same run at levels 0, 1 and 2
demo_messaging('levels', [1 2])     % just the default vs full detail
```

**Warnings are not verbosity-controlled.** A run that failed to converge says so
even at level 0, or a silent batch job would hide its own failure.

All output goes through one sink: set `'logFid'` to an `fopen` file id to capture
status text *and* warnings to a file.

---

## Correlated observations (block-diagonal P)

Fully supported, **without whitening**. Pass a full SPD weight **matrix** for a
batch instead of a vector:

```matlab
src = @(k) deal(A(idx{k},:), l(idx{k}), Pk{k});   % Pk{k} is nk x nk SPD
[sol, D, V] = lsqvce.engine.adjust(u, src, nBatch);
```

P is block-diagonal with the **batches as its blocks** — correlation within a
batch, independence between. That is the structure of most real geodetic data
(per-arc, per-pass, per-session noise models).

The diagnostics use **Baarda's generalised statistic**, which works directly on
the *physical* observations:

```
r_i = (Qvv P)_ii                          w_i = (e_i' P v) / sqrt((P Qvv P)_ii)
```

Both are computable batch-by-batch — `G_k = P_k A_k`, then
`(P Qvv P)_ii = P_k(i,i) − g_i' Qxx g_i` — at **O(u²·nk), identical to the
diagonal path**, with nothing of size n ever formed. Verified exact to 2e-15
against a dense reference.

Crucially, this **reduces algebraically to the classical `v_i/σ_vi` when P is
diagonal** (verified to 4.4e-16), so there is one code path, not two — and a
flagged `w_i` names *an observation*, not a linear combination.

### A bug the integrity check caught

`r_i` is a leverage complement bounded in [0,1] **only for diagonal P**. For a
correlated block it can legitimately exceed 1 (measured: **1.0001476**). The
diagonal code clipped it — which silently broke `sum(r_i) = n − u`. The
identity check found it immediately. That is what Level 0 is for.

### What is *not* supported

`lsqvce.engine.snoop` **refuses** correlated batches, and the refusal is principled:
removing one observation from a correlated block takes the **Schur complement**
of the block covariance, changing the covariance of every *surviving*
observation in that block. It is not a row deletion, so the rank-one Cholesky
downdate the function is built on does not apply.

You can still **find** the outliers (`diagnose` works fine). To act on them:
drop the whole affected block, rebuild the block covariance without the rejected
observation and re-run, or use `lsqvce.engine.robust`, which downweights rather than
deletes and leaves the block structure intact.

---

## Combining several data sets (arrow structure)

Several satellite systems observe a **common** parameter block; each has
**local** parameters (arc states, biases, ambiguities) that no other system
sees. The locals are eliminated system by system.

```matlab
cmb = lsqvce.engine.Combine(uC);                      % uC common parameters
for i = 1:nSys
    cmb.addSystem(uL(i), @(k) srcOf(i,k), nBatch(i));
end
sol = cmb.solve();          % xC, plus the locals by back-substitution
D   = cmb.diagnose(sol);    % full per-observation diagnostics
```

`batchFcn` returns `[AC, AL, l, p]` — common partials, local partials, o-minus-c,
weights.

**The definition of "local" that matters:** a parameter is local if *exactly one
data set* can see it. Several **sensors** within one system may share local
parameters — fine, same block. A parameter seen by two **systems** is *common*,
however local it feels; eliminating it as local is not an approximation, it is
wrong.

### QR is not optional here

Relative error in `xC` against a long-double SVD reference, near-collinear common
columns:

| cond(A) | 1e2 | 1e4 | 1e6 | 1e7 | **1e8** |
|---|---|---|---|---|---|
| `'normal'` (Schur) | 9.6e-12 | 5.6e-08 | 1.8e-04 | 2.2e-02 | **4.9e-01** |
| `'qr'` (Householder) | 1.9e-14 | 1.3e-12 | 1.1e-10 | 1.4e-09 | **1.2e-08** |

**Seven to ten orders of magnitude.** At cond 2.4e8 the normal-equation path has
*no correct digits*; Householder still has eight. Forming `N` — even only to
Schur-complement it — squares the condition number, and nothing afterwards undoes
that. `'qr'` eliminates on the *observation equations*.

### Two traps, both measured, both optimistic

**1. The cross-term.** `a'Qxx a` spans *both* blocks:

```
a'Qxx a = aC' Q_CC aC  +  2 aC' Q_Ci aL  +  aL' Q_ii aL
```

Dropping the cross term — tempting, since the reduced system only hands you
`Q_CC` — inflates `sum(r)` by **7.3%** and biases every `r_i` by **+0.066** (up to
+0.27). Observations look *better controlled than they are*, so the MDB is too
optimistic.

**2. The redundancy.** It is `n − u_total`, counting the **eliminated**
parameters. They left the *solve*, not the *model*: they still cost degrees of
freedom. Using `n − uC` deflates `s0²` by **7%**.

Both violate `sum(r_i) = n − u_total`, and the integrity check catches both
immediately. That is what Level 0 is for.

### Cost

| systems × locals | full N | eliminated | ratio |
|---|---|---|---|
| 5 × 10 | 3,025 | 125 | 24× |
| 20 × 20 | 164,025 | 425 | 386× |
| **100 × 50** | **25,050,025** | **2,525** | **9921×** |

Peak memory is `O(uC² + max(uL)²)`, never `O(u_total²)`.

### Sensor and day hierarchies

If your data is organised as **system / sensor / daily batch** with *per-system*
locals, the numerics need nothing new: a `(sensor, day)` batch is just another
batch inside its system's elimination block. Verified — sensor×day batching is
**bit-identical** to a flat structure, and the solution is invariant to batch
order to **3e-16**. Feed days and sensors in any order.

What you *do* need is bookkeeping, or a flagged observation comes back as
*"observation 4193"* instead of *"LRI, day 14"*:

```matlab
bx = lsqvce.engine.BatchIndex();
bx.addSystem('GRACE-C', {'KBR','LRI','GPS'}, days);
bx.setCount('GRACE-C', 'LRI', day, nObs);   % a sensor that was down is
bx.lock();                                  % simply never registered
...
bx.report(D);                               % per-sensor and per-day summary
```

**Per-sensor variance components.** KBR and LRI observe the same arcs, and their
quoted noise levels are rarely both right:

```matlab
VC = cmb.varComp(sol, @(i,k) sensorGroupOf(i,k), nSensor);
```

Verified unbiased: true sigma ratios of 3.0, 0.4 and 2.0 recovered as **3.032,
0.405, 2.021** (mean over 30 realisations).

**A third trap, measured.** The partial redundancy is `r_g = n_g − trace(Qxx N_g)`
with the **full** Qxx. Using only `Qcc` — all the reduced system gives you —
overstates `r_g` by exactly the number of eliminated locals (**6.9 vs 7**), and
since `s_g = q_g / r_g` that *understates* the variance component. `cmb.varComp`
uses the full Qxx.

**Resolution:** with `r_g ≈ 300` per sensor you cannot resolve a sigma ratio
better than ~6%. A 5% "correction" is noise. `VC.significant` says which are real.

**Not supported:** correlated (matrix) weights *inside* an eliminated system —
the generalised algebra there was not validated, and shipping unverified
numbers is the failure mode this toolbox exists to prevent. Decorrelate first
(pass `Lp'*AC, Lp'*AL, Lp'*l` with unit weights), or use the single-dataset path,
which *does* support correlated batches.

---

## Speed: QR vs normal equations

I assert QR's accuracy advantage throughout. Here is what it **costs**.

| u | 5 | 10 | 20 | 50 | 100 | 200 |
|---|---|---|---|---|---|---|
| **QR / normal** | 3.2× | 2.7× | 3.5× | 5.0× | 8.7× | 5.8× |

*(single dataset, n = 2e5, batches of 5000)*

`Combine` is similar — **2.4× to 6.9×**, growing with `uC` and `uL`.

The flop count predicts only **~2×**. The rest is **memory bandwidth**: `A'A` is
a single BLAS-3 SYRK call (cache-optimal), while a Householder sweep has far
worse locality. Same story in MATLAB — same BLAS/LAPACK.

Both paths are **linear in n**. The cost grows with **u, not n**, and it is
**flat in the conditioning** — so you don't pay more for QR precisely when you
need it more.

### The decision is about `cond(N)`, not speed

Measured on the *same* problem:

| cond(A) | normal err | QR err | verdict |
|---|---|---|---|
| 2e1 | 2.8e-14 | 3.9e-14 | normal is fine |
| 2e3 | 4.0e-10 | 3.6e-12 | normal is fine |
| 2e5 | 4.5e-06 | 3.9e-09 | QR recommended |
| 2e7 | 2.5e-02 | 1.8e-07 | **QR mandatory** |
| 2e8 | **2.3e+04** | 6.6e-07 | **QR mandatory** |

**Recipe.** You can't see `cond(A)` without doing the work — but you *can* see
`cond(N)`, and `cond(N) = cond(A)²`:

- `cond(N) < 1e8` → normal loses ~8 of 16 digits. **Fine.** Take the 3–9× speedup.
- `cond(N)` 1e8–1e12 → borderline. Use QR.
- `cond(N) > 1e12` → normal has ~4 digits left. QR is not optional.

Solve **once** with `'normal'` (cheap), read `cond(N)` from the report, and
re-run with `'qr'` only if you must. One cheap solve buys the answer.

**No free lunch:** buffering batches before re-triangularising gains at most
**1.07×**, and gets *worse* beyond a buffer of 5 (cache pressure). Run
`bench_qr_vs_normal` to measure on your own machine and data.

---

## Verified in MATLAB R2026a (v1.0 bring-up record, 2026-07)

*(A dated record of the first MATLAB execution — the current gate is the
current suite in the status block above.)* **125 / 125 tests passed.**
All four demos ran. The suite was finally executed in
MATLAB rather than validated only through a Python port — and it found **four
real bugs the port structurally could not see.**

| Bug | Impact |
|---|---|
| **`kCrit` never assigned on the `'mad'` path** | The robust-scale cure — the toolbox's *headline* fix for Pope's breakdown — was **dead on arrival**. Any user following its own advice hit a hard error. |
| **`extMax > 1.0` as a fixed ERROR threshold** | Fires on **100 % of clean datasets with n ≤ 500** and essentially never for n ≥ 2000. `extMax ~ 4·δ₀/√n`, so a fixed threshold is meaningless. A clean n=500 adjustment was reported `SUSPECT`. |
| **`sprintf` on a sparse `condN`** | The entire sparse path **crashed** in `report`. |
| **`E.xy` returned 2×npts, not npts×2** | `demo_workflow` plotted a *single point* instead of an ellipse — silently, with no error. |

The `extMax` fix is the substantive one. The question is not *"is extMax > 1?"*
but *"is it worse than a healthy network **of this size**?"* — so `D.extRatio`
now compares against the expected value. Verified: clean n=500 → `PASS`
(ratio 1.20); one dominant observation (`min r_i` 0.94→0.096) → `SUSPECT`
(ratio **8.67**).

### The headline claim, live

```
planted blunders : [6 13 15 20 43]        % 5/50 = 10% contamination, 25 sigma
POPE found       : []                     % max|w| = 3.482, crit = 3.897
MAD  found       : [6 13 15 20 43]        % max|w| = 7.618
contamLevel      : severe                 % fhat = 0.100
```

### Speed: I was too pessimistic

Measured **in MATLAB**, QR costs only **1.2–2.4×**, not the 3–9× my NumPy
benchmark suggested — MATLAB's `qr` is far better tuned for tall-skinny
matrices. Both paths agree to **2e-15** when well-conditioned. At
`cond(N) = 3.9e16` the normal path has error **0.806** (no correct digits) while
QR sits at **7.8e-09**.

*Always measure on the platform you ship on: `bench_qr_vs_normal`.*

---

## Documentation

```matlab
lsqvce_setup          % add to path, register the Help browser docs
doc lsqvce            % landing page
doc lsqvce.engine.diagnose   % any function
demo_workflow          % build your OWN case, from scratch  <- start here
run('dev/tests/run_all_lsqvce')   % the FULL unit-test suite
run_lsqvce('adjustvce')           % one suite, while iterating (dev/tests on the temp path)
```

Note: do **not** use `runtests` on a folder. It SKIPS any file that fails to
load and still reports the rest as green, so the suite can look clean while an
entire file never ran. `run_all_lsqvce` iterates the registered list, counts a
load failure as a failure, keeps the warning inventory, and records per-suite
seconds.

**Documentation ships in two forms, both self-contained.** The **master manual**
`lsqvce_manual.pdf` is the single document: theory, every entry point with a
workflow diagram, a step-by-step guide to adapting the template, interpretation,
and the reference tables. It is built with ReportLab — equations are drawn as
vector content, not images, and there is no LaTeX toolchain to install:

```bash
python3 dev/docs_src/build_master_pdf.py        # -> lsqvce_manual.pdf
python3 dev/docs_src/build_theory_pdf.py        # -> lsqvce_theory.pdf  (R4 WP-59)
```

The **theory document** `lsqvce_theory.pdf` (R4, Level L3) is the authoritative
account of the METHOD choices -- backends, conditioning, VCE families,
stochastic model, robust estimation, error handling, limitations -- with every
toolbox claim tied to the work package that pinned it; the manual's theory
chapter remains until decision R4-D3 retires it.

See `dev/docs_src/README_manual.md` for the file map and for how to add an equation
or a diagram. The **HTML manuals** carry the same material in browsable,
per-topic form, with every equation as a **static PNG image**, so the pages
contain **no JavaScript and no external references at all** — open them from
`doc lsqvce` or the Help browser, offline, with no runtime dependency.

- **[User Manual](html/lsqvce_user_manual.html)** — **start here.** Five parts: a
  referenced **Theory** section with typeset equations (Gauss–Markov, BLUE/QR,
  redundancy, data snooping, reliability/MDB, variance components, block-arrow
  Schur reduction, coloured-noise whitening), what you supply, reading the
  diagnosis, the gravity-recovery path, and a step-by-step recipe.
- **[Reference Manual](html/lsqvce_reference_manual.html)** — symbols, equations,
  matrix/return conventions; §4 details the block-arrow QR pre-elimination.
- **[Interpretation Guide](html/lsqvce_interpretation_guide.html)** — what to run
  and how to read the output.
- **[Choosing an engine](html/lsqvce_choosing.html)** — decision flowchart plus
  the exact-vs-Hutchinson cost curves.

The **API reference is the per-function HTML** (one page each, listed by category
on the `doc lsqvce` landing page); every function also has a runnable `EXAMPLE`
block:

```matlab
help lsqvce.engine.adjust      % EXAMPLE + hyperlinked See also
doc lsqvce              % the landing page
```

| File | Purpose |
|---|---|
| `lsqvce_manual.pdf` | **THE master manual** — one document, built with ReportLab (`python3 dev/docs_src/build_master_pdf.py`) |
| `lsqvce_theory.pdf` | **The theory document** (R4 WP-59) — theory, practice, implementation, error handling; `python3 dev/docs_src/build_theory_pdf.py` |
| `html/lsqvce_user_manual.html` | **theory + usage + step-by-step. Start here** |
| `html/lsqvce_reference_manual.html` | symbols, equations, conventions, whitening; §4 block-arrow QR |
| `html/lsqvce_interpretation_guide.html` | how to interpret the results |
| `html/lsqvce_choosing.html` | which engine + exact vs Hutchinson (decision aid) |
| `html/lsqvce.*.html` | **per-function API pages — params, size, type, default, description** |
| `html/` + `helptoc.xml` + `info.xml` | MATLAB `doc` pages |
| `template/lsqvce_template.m` | **your own case: copy, rename, fill the two stubs** |
| `examples/demo_workflow.m` | build-your-own-case walk-through, section by section |
| `examples/demo_lsqvce.m` | six failure modes |
| `dev/tests/test_lsqvce.m` | core science (24) |
| `dev/tests/test_lsqvce_ext.m` | failure modes (34) |
| `dev/tests/test_lsqvce_opts.m` | options, edge cases, correlated P (42) |
| `dev/tests/test_util.m` | `getfielddef`, `mergeStruct` (10) |
| `dev/tests/test_vceFunction.m` | `vce` + `diagnoseAccum` (KeepLocal guard) |
| `dev/tests/test_signalCore.m` | `outliers`, `varWhiten`, `ar2arma` |
| `dev/tests/test_streamUtils.m` | `selectArAcv`, `residualAcv`, `collectResiduals`, `whitenessCorr`, `normalizeSatList`, `resolveSatSpecs`, `satKeyOf` (14) |
| `validation/*.py` | Python mirrors of the numerical properties |

---

Generated with assistance from Claude (Anthropic), Claude Opus 4.8, 2026-07-11.