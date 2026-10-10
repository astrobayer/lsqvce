# Review R7 — four-lane review of lsqvce 6.1.3 (2026-10-09)

*Created 2026-10-09 23:45 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 — review R7, findings only. No code under `+lsqvce/` was changed.*

**Tree:** `main@b25bb11`, version 6.1.3 (`+lsqvce/Contents.m`), 27 suite files.
**Tier:** B — no MATLAB in the session. Evidence is the repo's own audit, grep, and
line-by-line **Python mirrors** of the MATLAB arithmetic (numpy 2.5 / scipy 1.18), kept
in `mirrors/` next to this file. A mirror is not the MATLAB code: every S1/S2 below
must be reproduced at tier A (MCP bridge) as the FIRST step of the work package
that fixes it, before any code is touched.

**Instruments run here**

| instrument | result |
|---|---|
| `python3 dev/docs_src/audit_tree.py` | exit 0, CLEAN, 29 checks (pandoc present, Manual rebuilt) |
| `python3 dev/docs_src/audit_selftest.py` | every injected fault caught (26 cases) |
| skill `audit_help.py` (project sections) | 2 identifier gaps worth reporting (R7-28); rest is template noise |
| skill `count_tests.py` | 0 classes — the script parses classdef only, all 27 suites are function-based; it fails open (R7-46, skill-side) |
| own parse (lane swe_tests) | **576 declared / 564 gate-executed** = HANDOVER §12.34 exactly |

**Method.** Four reviewer lanes run in parallel per the review discipline (findings
only, ranked, REPRODUCED vs INFERRED on every line, "checked and sound" stated):
two software-engineering lanes (code/interfaces; tests/gate/audit) and two domain
lanes (joint VCE estimator; linear model / solvers / robust / Gauss–Newton /
diagnostics). Common brief: `BRIEF_COMMON.md`. The full lane reports are attached
verbatim and are the evidence; this file is the index. The coordinator re-checked the
top findings independently (marked ✔ below): R7-01/02 by a fresh mirror
(`mirrors/coordinator/`), R7-03/06/14/26 and C1 by reading the cited lines.

Severity: **S1** silently wrong numbers · **S2** wrong under reachable conditions, or a
guard/refusal that misfires · **S3** contract/doc/test gap · **S4** hygiene.

---

## Findings (most severe first)

### S1 — silently wrong numbers

| id | finding | lane ref | evidence |
|---|---|---|---|
| R7-01 ✔ | **The declared PSD law becomes a CIRCULANT covariance.** `covModel>psdModel` 'toeplitz' interpolates onto exactly `floor(m/2)+1` frequencies and `psd2acv` integrates by trapz on that grid: the ACV is periodic with period 2(nf−1), so r(m−1) = r(1) (even m) or r(0) (odd m). Default path (Grid 60, whitener "toeplitz") on every arc > 118 epochs. 2111-epoch arc, true law fed in: slope formal sd 90× too small. | dom_vce F2 | REPRODUCED (two independent mirrors: m = 256 → r(255)/r0 = 0.500; m = 255, 2111 → 1.000) |
| R7-02 ✔ | **Estimated law runs away ×10⁷ and the loop says "converged".** `psdVCE` returns C by the ORTHONORMAL inverse DCT-I (lag-0 weight 1/√2), `jointModelOpt>lawSpecs` re-forms S by the plain cosine pair: not inverses of each other, so a positive Ŝ yields negative S bins (35/60 at φ 0.8), floored at eps·max → cond 1e9 → level jumps, σ₀ collapses; level is outside the stop rule. AR(1) φ 0.8: 8/12 seeds "converged" with coefficient z up to 37. Mechanism of the §12.20 "level runs away" item. | dom_vce F1 | REPRODUCED (coordinator mirror: φ 0.8 → 35/60 negative bins, Toeplitz(C) min eig −4.15) |
| R7-03 ✔ | **The default QR backend squares the condition number on the accumulator path.** `QRSystem.blocks` forms `R1'R1`, the Schur complement is done in normal form, `merge` re-roots by `chol` — also for nl = 0. The help's "triangular and well scaled … a very different proposition" is false: cond(R1'R1) = cond(A)². At cond(A) = 1e8: median 2.3 formal σ error; full-rank problems refused from 1e8.5. Contradicts DEVELOPMENT.md:460–473 and the §12.2 "QR default" rationale. | dom_lsq F1 | REPRODUCED (mirror) ✔ code read |
| R7-04 | **Nonlinear fit with `data.priors`: the prior constrains the Gauss–Newton increment, not x** (`K·dx ~ x0` every pass). Tight prior b ≈ 0.30 → fit returns b = 0.60. Untested combination. | dom_lsq F3 | REPRODUCED (mirror) |
| R7-05 | **Datum prior: x exact for any σ, Cxx is not.** `Cxx = s0²·inv(N + BB'/σ²)`; `QxxCorr` is read by `diagnose` but set nowhere. Levelling fixture: se 0.005 (σ 1) vs 2.48 (σ 1e3). | dom_lsq F4 | REPRODUCED (mirror) |
| R7-06 ✔ | **`util.rowWeightKey` formats days with `%g`** → day ids ≥ 1e6 collide (YYYYMMDD: 30 days → 1 key; JDN → 3). With `Robust` on, one day's blunder factors are applied to other days (equal-length days) or a misleading `rowWeightSize` blames the builder. The loader cache keys the same day with `%d` (R1: two definitions disagree). MJD/datenum safe. | swe_code F1 | REPRODUCED (mirror) ✔ code read |
| R7-07 | **`var` law under `Estimate="noise"`**: `oneArc` calls the var Cfcn with an epoch count where it expects rows (wrong Σ: lag-1 0.000 vs 0.667), and the estimated level is applied as a std to a covariance. Nothing refuses the combination; no test runs it. | swe_code F2 | REPRODUCED (arithmetic) / INFERRED (reach) |
| R7-08 | **Joint-loop robust pass on multi-axis streams zeroes clean data.** Epoch statistic = RMS over axes (χ-shaped, median ≈ 1), kernel standardises `vw/s` uncentred → clean 3-axis data: IGG3 zeroes 24 % of epochs (62 % at d = 6). Contradicts §12.29 round-4 "adequate for blunder rejection"; oracle O-105 is scalar only. | swe_code F3 | REPRODUCED (mirror) |
| R7-09 | **Arc factors of one stream normalised by other streams' levels** (γ̄ over all columns, `adjustVCEjoint.m:506`). Stream A's reconstructed σ = 1.01 / 0.79 / 1.15 depending only on stream B's truth. Contradicts the WP-100 Q1 closure for multi-stream runs. | dom_vce F3 | REPRODUCED (mirror) |
| R7-10 | **Koch fold under a correlated law conditions on the rejected epoch's noise being ZERO** (multiplicative D with D_ii → 1e6). φ 0.9: formal/true sd 0.31, estimate 1.38× worse than deleting the epoch. Contradicts WP-105's premise ("the filter adapts around it"). | dom_vce F4 | REPRODUCED (exact linear algebra) |
| R7-11 | **A run with no `"law"` stream stops "converged" on pass 1** under the DECLARED weights; prior sigmas never estimated (demo_lsqvce case 5 returns σ₀ 2.31 under "variance ratio back at ~1"). Generally: returned solution one pass behind the returned model; levels, arc factors, prior sigmas outside the stop rule (brief lead L1 — holds). | dom_vce F6 | REPRODUCED (mirror) |
| R7-12 | **The robust pass biases the stream level low**; "UNBIASED … the weights cancel in α²" is false (s₀ contains tr(Σ_t⁻¹), untouched by weights). Clean AR(1): level 0.83–0.86 law fixed, ≈ 0.93 law free; the oracle's 15 % band cannot see it. | dom_vce F7 | REPRODUCED (mirror) |

### S2 — wrong under reachable conditions, or a guard that misfires

| id | finding | lane ref | evidence |
|---|---|---|---|
| R7-13 | **Per-arc factors never reach the weighting** (only their mean, `jointModelOpt.m:126–134`). Consistent with the WP-100 convention (σ_arc = declared·√factor·σ₀) but the Manual §7.8 says a bad arc "loses weight automatically", and Cxx is 2.6× optimistic when arcs truly differ. **Design question for the user.** | dom_vce F5 | REPRODUCED (consequence) ✔ code read |
| R7-14 ✔ | **Gauss–Newton damping is a no-op**: `'Reg','lm'` damps only the local block, and `preeliminate` returns at `nl == 0` before `applyRegularisation`. Rejected steps are retried unchanged until MaxIter / λ > 1e12; `R.x` is then the REJECTED point (F 36.8 vs best accepted 4.9); first step has no descent check. Leads L2/L3 hold. | dom_lsq F2 | REPRODUCED (mirror) ✔ code read |
| R7-15 | **Common priors are outside `residuals()` / `redundancy()`** → false `sigma0Mismatch` (blames the loader) and `identityViolated`; the wrong dof is marked authoritative. Contradicts D8-2 "no second path". | dom_lsq F5 | REPRODUCED (arithmetic) / INFERRED (warning fires) |
| R7-16 | **Free network without datum: refusal decided by `chol` rounding** (61 % refused / 39 % solved with arbitrary level, sd ~1e7). H-16's lesson moved into `merge`. | dom_lsq F6 | REPRODUCED (2000 networks) |
| R7-17 | **IRLS stop rule not unit-invariant** (`Tol·max(max|x|,1)`): SH-scale coefficients (1e-9) stop at pass 2 "converged", ~1 σ from the fixed point. | dom_lsq F7 | REPRODUCED |
| R7-18 | **Robust σ₀ and Cxx reported as least-squares**; `report()` never says so although the help claims it does. se/MC-sd 0.73–1.86. | dom_lsq F8 | REPRODUCED (numbers) |
| R7-19 | **`Results.compare` df2 omits eliminated locals** → 7.4 % false-positive rate at α 0.05. | dom_lsq F9 | REPRODUCED (2000 trials) |
| R7-20 | **Stop criterion contains α²₀, which ZeroFreqRule "lowest" overwrites** → a settled model never converges (`maxPasses`, `NOISE_UNSETTLED`). | dom_vce F8 | REPRODUCED |
| R7-21 | **`robustWeights` falls back MAD → std → 1 in silence**; near-ties (MAD 3.5e-15) give every ordinary row factor ~7e-15 and the collapse guard cannot fire. | swe_code F4 | REPRODUCED |
| R7-22 | **`Noise.Smooth` / `Noise.FreezeRedundancy` have no consumer**; the GRAVFIRE guide's recipe sets `FreezeRedundancy=100` and runs bit-identically. Audit check 27 sees top-level properties only. | swe_code F5, dom_vce F9 | REPRODUCED (grep) |
| R7-23 | **`priorNotActive` lost from `R.meta.warnings` on kinds B/C** (fires before the inventory reset) and suppressed entirely after a kind-A run fired it. | swe_code F6 | REPRODUCED (mirror of call order) |
| R7-24 | **`diagnoseAccum` keeps every whitened row at full width** (n × uFull) inside `Results`: ≈ 228 GB at GRAVFIRE month scale; README "nothing of size n is ever formed" is false. Extends open WP-86. | swe_code F8 | REPRODUCED (arithmetic) / INFERRED (code) |
| R7-25 ✔ | **The claim WP-110 retired survives in `adjustVCEjoint`**: help l.175–177 (line-wrapped, so check 28's single-line regex misses it) and the runtime `collapsed` message ("beyond the breakdown point of any M-estimator") — check 28 never reads message bodies. The guard compares live EPOCHS with nc + 1 PARAMETERS. | coordinator C1 | REPRODUCED (grep) |

### S3 — contract, documentation and test gaps

| id | finding | lane ref |
|---|---|---|
| R7-26 ✔ | **Seven audit checks (2, 3, 7, 10, 10b, 11, 15) scan `tests/`, which does not exist since L8**; checks 9 and half of 3 read deleted files. Faults injected into a test file leave the full audit CLEAN. (The MATLAB gate still catches load errors.) The learning "the inventory is the tree, never a hand-picked list" recurs. | swe_tests F1 |
| R7-27 | **Category coverage overstated**: check 16 credits comment mentions and the never-gated `test_speed`; honestly counted 18 gaps without exemption, incl. `contract` for `adjustNonlinear`, `diagnose`, `resolveNc`. Six `contract` exemptions cover validating units with "inconvenient" reasons; 24 rows are stale and nothing refuses them. | swe_tests F2, F3 |
| R7-28 | **156 of 301 catalogued ids fired by no test** (3 user-tier: `Results:noAccumulator`, `Results:noVPv`, `fit:lawOnFixedStructure`); `covModel:badCep`, `adjustVCEjoint:badAxisRows` still untested against §12.14 E-10's disposition. `adjustVCEjoint:badSigma2_0` raised, not in its ERRORS block, never fired. *(C2 corrected: `arcMismatch` IS fired, test_ellmer.m:1383.)* | swe_tests F4, coordinator C2 |
| R7-29 | **`testDistributeDaysThroughFitEqualsFused` cannot fail**: the second `Parallel=` resets Mode to "serial" (set.Parallel merges onto defaults); the fused path is compared with itself. | swe_tests F5 |
| R7-30 | **`warning('off','all')` in 16 points** blinds the warning gate there; the onCleanup-in-a-loop in test_reliability alternates the state. | swe_tests F6 |
| R7-31 | **Copied `tol:` comments name the wrong error source** (4 sites); check 16(d) only greps the substring; 146 numeric verify bounds carry no `tol:`. | swe_tests F7 |
| R7-32 | **Fault injection covers 10 of 30 audit checks**; check 14 accepts ids raised only by test files. | swe_tests F11 |
| R7-33 | **HANDOVER §6.1 contradicts the gate**: names `run_lsqvce2`, six categories, "absolute seconds logged, not asserted". | swe_tests F12 |
| R7-34 | **Deleted units / retired behaviour in live text**: runtime remedy `snoop.m:557` ("re-run lsqvce.engine.adjust … use lsqvce.engine.robust"); `Options` EXAMPLE still says Robust × noise is refused (false since 6.1); `adjustRobust` ACCURACY cites `engine.robust` (C3); ~30 help lines; `Results.m:96–98, :760`; `adjustVCEjoint.m:14–17, :1049–1051`; `declaration.m:18–20` / `badParamColumns` undocumented. | swe_code F9, dom_vce F12, coordinator C3 |
| R7-35 | **`diagnose`**: `'verify'` documented and never read; `RMin` default 0.01 (Options) vs 1e-6 (diagnose, reached via `snoop`'s re-diagnosis); many accepted Name-Value pairs undocumented across drivers. | swe_code F10 |
| R7-36 | **`diagnose.m:654` draws `randi` from the GLOBAL stream** past 5000 MDBs: MDB quantiles differ run to run; contradicts WP-36. | swe_code F11 |
| R7-37 | **Unseeded joint runs record no probe seed** (§12.2 "seeds in provenance"). | dom_vce F10 |
| R7-38 | **Ljung–Box K = 20 fixed; pooled per-arc lag sums mis-size it** (64 × 30 epochs: size 0.000 instead of 0.05). | dom_vce F11 |
| R7-39 | **A-posteriori ellipse uses χ²(2) with estimated s0²** → "95 %" covers 86 % at r = 5 (needs 2·F(2, r)). | dom_lsq F10 |
| R7-40 | **`localCovariance` (verbose log) scaled by `s0sq_hint`**, ~300× too large; the user-facing Cyy is correct. | dom_lsq F11 |
| R7-41 | **Parallel log queue never reaches the workers**: worker warnings miss the inventory; `warnOnce` per worker. | swe_code F7 (INFERRED) |
| R7-42 | **Speed budgets**: two absolute-seconds asserts in the gate; one sequential `timeit` pair; `testFftLagSumsBeatTheLoop` times test-local copies, not the kernel. | swe_tests F8 |

### S4 — hygiene

| id | finding | lane ref |
|---|---|---|
| R7-43 | Stream/arc keys composed inline in 8+ places and parsed back by regex (the R1 pattern behind R7-06). | swe_code F12 |
| R7-44 | Five functions > 400 lines without a recorded reason (`report` 826, `diagnose` 1001, `adjustVCEjoint` main 632 …); 74 lines > 120 chars. | swe_code F13 |
| R7-45 | `test_parallel` leaves the pool open; RNG and path not restored by 7+ suites standalone; a mute on an id nothing raises. | swe_tests F9, F10 |
| R7-46 | (skill instrument, not the repo) `count_tests.py` parses classdef only and exits 0 on zero classes. | swe_tests F13 |

---

## Checked and sound (consolidated; details in the lane reports)

- Block-arrow elimination, back-substitution, recovered y and Cyy = s0²N_ll⁻¹ + G·Cxx·G', dof with locals — equal to the dense solve (1e-15 … 6e-11). The **arithmetic** is right; R7-03 is about where it is done.
- Per-lag Ω_n and s_n = tr(R̄V_n) (Gohberg–Semencul path, FFT lag sums) equal the dense forms to 1e-16; the traced matrix is the redundancy, not H (WP-32 holds); Hutchinson via `sharedProbe` unbiased with the theoretical variance at Z = 100.
- DCT-I arithmetic ('ortho' self-inverse, equals scipy; 'ellmer' inverse exact) — R7-02 is the pairing, not the transform.
- Multi-axis Kronecker order: descale + per-axis temporal whitening equals A'[D(R_t⊗I_d)D]⁻¹A to 7e-16; per-axis statistics sum correctly.
- Eq. 6.5.28 after WP-98: exact cancellation; level × shape enters the whitener once on a single stream.
- Koch–Kusche prior update: Helmert form with the prior-including cofactor; E[Ω_k] ≈ r_k.
- Huber 1.345, bisquare 4.685, IGG3 (1.5, 3.0), MAD × 1.4826; Pope τ critical value, Baarda δ₀, χ²/F inverses (to 2e-15 vs scipy).
- Gate logic (`gateVerdict`): red on load error, failures, incompletes, set differences, declared ≠ executed, foreign warning id, blind capture; `ci_gate` errors so CI exits non-zero; assume-filtered points are red.
- Layering (check 21 baseline only), no unqualified sibling calls, no eval/global, no bare catches; probe RNG on private streams; identifier catalogue 301 = 301.
- Test count prediction: 576 declared / 564 executed, reconciled per suite.

## Not reviewed

MATLAB execution of anything (tier B); the physical oracle O3 and the GRAVFIRE bundle O7;
`snoop`'s re-solve ladder end to end; skip-bias and C0 dependence (brief lead L4);
Pool "obsType"/"all" with several law groups; `psdMode` "ar"/"arma"; `release.yml`;
html pages beyond the stale-reference sweep.

---

## Proposed work packages (for the user's decision — none opened)

Order by severity and dependency. Each WP starts by **reproducing its findings at tier A**
from the cited mirror, and ships the mirror's numbers as its oracle.

| WP | scope | findings | oracle / acceptance |
|---|---|---|---|
| 111 | **The PSD → covariance path**: an alias-free ACV (grid ≥ 2m, or the closed-form ACV where the law has one); ONE forward/inverse transform pair so positivity of Ŝ means positivity of the whitener's C; floors at a stated relative level; α²₀ out of the stop rule under "lowest" | R7-01, R7-02, R7-20 | r(m−1) ≈ true ACV on odd/even arcs; GLS under the true AR(1) law vs dense exact Toeplitz in formal-σ units; φ 0.8/0.9 runaway fixture converges with |z| ≤ 3 |
| 112 | **QR on the accumulator path**: eliminate locals on the observation equations (triangular, no `R1'R1`), rank decision on R's pivots by the code's own threshold | R7-03, R7-16 | cond(A) 1e7–1e9 fixture: x_c vs dense `[A B]\l` in formal-σ units; no full-rank refusal; deterministic refusal of the free network |
| 113 | **Gauss–Newton**: damp the common block; prior as `K·dx ~ x0 − K·x`; descent check from pass 1; return the best accepted point | R7-14, R7-04 | exponential/Rosenbrock fixtures reach the optimum; prior fixture hits the declared regularised optimum |
| 114 | **Priors without a second path**: common priors in `residuals()`/`redundancy()`; datum covariance (minimum-constraint cofactor, `QxxCorr` set or the field deleted) | R7-05, R7-15 | σ-independence of datum se; zero `sigma0Mismatch`/`identityViolated` on prior fixtures |
| 115 | **Robust conventions**: struct/`%d` row keys from one helper; a centred / χ-calibrated multi-axis epoch statistic; Koch fold vs nugget/deletion under correlated laws (**user decision**); scale-invariant stop; honest labelling of robust σ₀/Cxx in `report`; fallback warned or refused; `collapsed` wording and epoch/row threshold | R7-06, R7-08, R7-10, R7-12, R7-17, R7-18, R7-21, R7-25, R7-43 | clean 3-axis data: zeroed fraction ≈ the scalar case; YYYYMMDD days pinned; correlated-law fold fixture (formal/true sd ≈ 1) |
| 116 | **Joint-loop policy**: γ̄ per stream; complete stop rule (levels, arc factors, prior sigmas) and the returned solution under the returned model; scale-only runs iterate; `var` + noise refused by name; dead `Noise` fields wired or retired; arc-factor weighting or Manual correction (**user decision**) | R7-07, R7-09, R7-11, R7-13, R7-22 | metamorphic: stream A invariant to stream B's truth; demo_lsqvce case 5 σ₀ ≈ 1 |
| 117 | **Instruments**: audit checks on `dev/tests`; check 16 from code lines and gated suites only, stale-exemption refusal, contract-exemption rule enforced; check 27 on group fields; check 28 on message bodies and wrapped lines; fault cases for every check; the cannot-fail test; warning mutes via `muteWarn`; speed budgets as ratios on the kernel; §6.1 rewritten | R7-26…R7-33, R7-42, R7-45 | each new check fault-injected (fires on broken, silent on healthy) |
| 118 | **Docs truth and hygiene**: deleted-unit references, Options example, undocumented NV pairs, `RMin` one authority, `diagnose` RNG, probe seed in provenance, Ljung–Box K, ellipse 2F, verbose Cyy, worker log queue, README memory claim (with WP-86), long functions | R7-19, R7-23, R7-24, R7-34…R7-41, R7-44 | audit CLEAN with the new checks of WP-117 |

**Decisions the record needs from the user:** (a) R7-13 — do arc factors weight the
solve, or does the Manual change? (b) R7-10 — keep the multiplicative Koch fold (and
document its limit under strong colour) or move to deletion / an additive nugget?
(c) the §12.34 collapse-threshold options, now joined by R7-25's epoch/row unit
question; (d) whether WP-86 (memory) is pulled forward given R7-24.

**Not to be over-read.** Twelve S1 rows is not "the toolbox is wrong everywhere": the
core arithmetic re-derived above is sound. The S1s concentrate in four places — the
PSD → covariance pairing, where the QR elimination happens, the Gauss–Newton/prior
wiring, and the robust × correlated × multi-axis combinations that WP-105/110 added.
