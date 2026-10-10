# Review — domain: variance-component estimator (dom_vce) — main@b25bb11 (6.1.3) — 2026-10-09

Scope reviewed: `+engine/adjustVCEjoint.m`, `psdVCE.m`, `arcVarianceFactors.m`, `jointModelOpt.m`,
`covModel.m` (psd / Cep composition), `sharedProbe.m`, `traceEstimate.m`, `fit.m>priorUpdate/priorFor/streamPolicy`,
`makeLoader.m>batchOf/inflateCep/priorArc`, `+kernel/redundancyLags.m`, `dctI.m`, `psd2acv.m`, `separableWhiten.m`,
`levinsonDurbin.m`, `robustWeights.m`; theory `dev/docs_src/theory_r5/ch_theory_vce.tex`; Ellmer 2018 §6.5 (project copy).
Chosen because the brief's leads L1/L4/L5 and the WP-98/100/105 closures all sit on this path.

Tier: B. No MATLAB. Python 3 (numpy 2.5 / scipy 1.18) mirrors in
`scratchpad/mirrors/dom_vce/` — `kern.py` mirrors `dctI`, `columnLagSums`, `toeplitzInvFirstCol`,
`sigmaInvSubdiagSums`, `redundancyLags`, `jointModelOpt>lawSpecs`, `covModel>psdModel('toeplitz')`
(interp1-linear, floor 1e-10·max, ridge 1e-8, BatchSize 20000) and `psd2acv` (trapz) line by line;
`loop.py` mirrors `adjustVCEjoint` steps 1–10 for one stream (Pool any, `estimate="law"`, Gamma 1,
ZeroFreqRule "lowest", fStream = mean γ²) with EXACT traces (UiZ = U⁻¹) unless stated.
Fidelity check of the loop mirror against the record: it reproduces HANDOVER §12.20 / `demo_noisePolicy`
LIMITATIONS ("law recovers AR(1) up to φ = 0.6; at φ = 0.7 the level runs away") — φ = 0.6, 9×256:
max level/truth 1.49–1.67 over 3 seeds; φ = 0.7: 3.5e6–1.9e8.

Lead status (brief v5.2 → tree 6.1.3):
L1 HOLDS (stop rule = PSD amplitudes of "law" groups only, `adjustVCEjoint.m:476–489, 599`) and is
worse than stated — see F6, F8. L2 still present (`adjustNonlinear.m:102 opt.Lambda0 > 0 && it > 1`;
outside this lane). L3 not reviewed. L4 not separately measured (see F8 for start-at-truth). L5 HOLDS (F11).
L6 PARTLY closed (deleted-unit names gone from fit/report/log; residue in F12).
WP-98 (F-01 recursion) CLOSED — the cancellation it rests on is exact; WP-100 Q2 (one C weights
every arc) not re-opened; WP-100 Q1 is CONTRADICTED for multi-stream runs (F3).

## Findings (most severe first)

### F1 — The estimated law runs away to 10⁷× its level and the loop reports "converged": DCT-'ortho' inverse is not positive-definite and the PSD floors are at machine precision — REPRODUCED
What is wrong: two things compose. (i) `psdVCE.m:248` returns `Cnew = X·Ŝ` with the orthonormal DCT-I,
whose lag-0 row carries k₀ = 1/√2. A positive Ŝ does NOT give a positive-definite Toeplitz(C): every
single-bin positive spectrum gives an indefinite one. `jointModelOpt>lawSpecs` then forms the
cosine series of that C, which goes negative, and floors it at `eps*max(S)` (`jointModelOpt.m:230`);
`covModel` floors at 1e-10·max; `psdVCE` floors Ŝ at an absolute 1e-12 (`MinAmplitude`). (ii) A
spectral bin at ~eps makes the declared Toeplitz near-singular (cond 1e9); the next Helmert step
sees almost no residual in that direction, the level jumps ×10⁷ and σ₀ drops to 0.056 — and since
level and σ₀ are not in the stop rule, α² ≈ 1 on the next pass and the loop stops as "converged".
Why: Wiener–Khinchin positivity holds for the pure cosine pair (row weight 1 at lag 0), not for the
orthonormal DCT-I; "positive by construction" (`jointModelOpt.m:141`, `covModel.m:537`, Manual §7.5)
is true of S, false of the covariance the whitener is built from.
Observable symptom: on a plain AR(1) stream (no locals, 1 stream, exact traces, shipped path, Grid 60,
MaxPasses 10) — φ = 0.8, 20×256, 12 seeds: 8 stop "converged" at pass 4–6 with level 0.7e7–5e7 ×
truth, σ₀ 0.056, and coefficient z-scores (x̂−x)/σ̂ of **11, −6.6, 23, 20** and **−10, 37, −0.8, 17**
where GLS under the true law gives |z| ≤ 2.3. φ = 0.9/0.95: same runaway, `maxPasses`.
Evidence:
- `t9_dctdomain.py`: "60/60 single-bin POSITIVE spectra give an INDEFINITE Toeplitz(C)"; random positive Ŝ:
  4764/5000 indefinite; AR(1) 0.9 spectrum ±30 % per bin: 4929/5000. Pure-cosine columns: min eig −1.4e-14 (PSD).
- trace (φ 0.8 seed 10): pass 3 `lawSpecs raw S min/max -1.9e-02 (neg bins 24)`, `cond T(r) 1.09e+09`,
  level 1.21 → 2.395e7; pass 4 `maxDev 6.31e-03` → converged (`t10_converge_garbage.py`, z-scores in the inline run).
- counterfactuals (`t17_rootcause.py`, `t16_fixpair.py`): flooring the declared series at 1e-3·max instead of eps
  cuts the runaway from 10⁶–10⁸ to 10–50× (still wrong); a Wiener–Khinchin pair alone does not cure it
  (α²_j ≤ 0 bins → absolute 1e-12 floor → same singularity); `StepSchedule="decay"` cures φ ≤ 0.8 but not 0.9 (3.7e6×).
- Grid = arc length (Ellmer's design, 256) runs away too (`t8_grid.py`: level 2.2e7, σ₀ 0.35).
Severity: S1. This is the mechanism of the recorded-not-fixed §12.20 "second, milder" item, and shows it
carries a `converged` verdict, which §12.20 only recorded for the 3-stream variant.

### F2 — The declared covariance is CIRCULANT: epoch 1 is correlated with epoch m like a neighbour; on ODD arcs corr(1, m) = 1.000 — REPRODUCED
What is wrong: `covModel>psdModel` 'toeplitz' interpolates the law onto exactly `nfNeed = floor(m/2)+1`
frequencies (`covModel.m:541–547`) and integrates with trapz on that grid (`psd2acv`). A trapezoid
cosine sum on nf equispaced points is periodic in the lag with period 2(nf−1): for even m, r(m−1) = r(1);
for odd m the period is m−1, so r(m−1) = r(0). The comment "(sharp: at nf = m/2 the matrix is
indefinite, at nf = m/2 + 1 positive definite)" (`:534`) chose the grid for PD, not for correctness.
Reached by every "law" stream with Grid < m/2+1 (the default Grid 60 on any arc > 118 epochs) under
the default `Noise.Whitener="toeplitz"` (`fit.m:308`).
Observable symptom: even m: a spurious end-to-end correlation equal to lag 1 (0.4999 for φ = 0.5);
odd m: Toeplitz cond 3.3e8 (only the 1e-8 ridge keeps it PD). GLS with the law fed in EXACTLY,
arc of 2111 epochs (the GRAVFIRE KBR day length in §12.14): slope true sd **1.46 vs 0.142** without the
wrap, formal sd **0.016** (90× over-confident). The joint loop on 9 arcs × 255 (`t12_odd.py`):
level 1.3e8, |x̂−x| up to 10.7 against 0.69 at m = 256, maxDev stuck at 0.038.
Evidence: `r[255]=r[1]` to 8.8e-15 (m = 256); m = 255/257 `r[m-1]/r[0] = 1.000000`, cond 3.331e+08
(inline run); `t19b.py` isolates the wrap (same S on a 16× finer grid = no-wrap column).
Severity: S1.

### F3 — Arc factors of one stream are normalised by the other streams' levels (cross-stream γ̄) — contradicts the WP-100 closure for multi-stream runs — REPRODUCED
What is wrong: `adjustVCEjoint.m:506` calls `arcVarianceFactors` on ALL columns (every stream, every
policy), so γ̄ (`arcVarianceFactors.m:140`) is a central-half mean over the absolute γ_m of streams that
carry their own level. An `arcFactors=true` stream's factors become γ²_m/γ̄², and γ̄ moves with
UNRELATED streams. WP-100 closed the level/factor split as "a convention, σ₀ carries the rest" — true
on its one-stream fixture only: σ₀ is shared, it cannot absorb a per-stream factor.
Observable symptom (metamorphic; `t15_gammabar.py`, three white streams, Pool "stream", exact traces):
stream A (`scale`, arcFactors, truth 1) reconstructed with the Manual's own rule
declared×√factor×σ₀ = **1.011 / 0.785 / 1.153** when ONLY stream B's truth is 1 / 4 / 0.25;
all runs `converged`. demo_noisePolicy escapes by luck (γ̄² = 1.052 for its mix).
Severity: S1.

### F4 — The Koch fold does not remove a rejected epoch under a correlated law: it conditions on that epoch's noise being zero — contradicts WP-105's premise — REPRODUCED
What is wrong: `makeLoader>inflateCep` (`:1279`) sets σ_i/√w_i with w_i → 1e-12, i.e. C = D R_t D with
D_ii → 10⁶. For correlated R_t, observing l_i with noise d_i·η_i, d_i huge, pins η_i ≈ 0: the
neighbours' noise is then predicted from a value that was never observed. P restricted to the kept
epochs is the inverse of the CONDITIONAL covariance, not of R_keep. WP-105 round 8 says the Koch
weight makes "a blunder one observation of huge variance, the filter adapts around it"; that holds
only for a white or weakly coloured law. The test_vcePolicy control ("blunders declared through Cep,
sigma 1e6") uses the same construction.
Observable symptom (exact linear algebra, `t22b.py`, 256-epoch arc, sinusoid coefficient):
φ = 0.9, 13 epochs zeroed: formal/true sd **0.31**, true/optimal (deletion) **1.38**; φ = 0.95, 5 zeroed:
**0.25** and 1.76; φ ≤ 0.5: ≥ 0.96 (harmless). Lag-0 E[Ω₀]/s₀ = 1.343 at φ 0.9 (`t22_kochfold.py`).
Severity: S1 (formal Cxx 3–6× too optimistic on strongly coloured streams after a robust pass).
Note for the fix: an additive nugget (C + diag(c_i)) or deletion removes information; a multiplicative D does not.

### F5 — Per-arc factors never reach the whitener; the Manual says a bad arc "loses weight automatically" — REPRODUCED (consequence) / mechanism by grep
What is wrong: `model.sigma2` is consumed only by `jointModelOpt.m:131–134` (mean) and the per-stream
level (`:176`); grep finds no other reader in `+engine/makeLoader.m`, `covModel.m`, `BlockAccumulator.m`.
WP-98 fixed the recursion, not the application: with `ArcFactors=true` (the Options default) x̂ and Cxx
are those of equal per-arc weights; the factors are reported, not used. Manual §7.8 item 2 ("The
arc-wise variance factors ARE the down-weighting … loses weight automatically") and the arcSpreadGuard
message (`adjustVCEjoint.m:1049`) disagree with each other and with the code.
Observable symptom (one of 9 arcs 10× in σ, white): true sd of the shipped estimate / Ellmer-weighted
**1.3–3.3**, formal/true sd **2.6** on three of four coefficients (inline run).
Severity: S1 (Cxx of the solve misstates the precision whenever arcs genuinely differ).

### F6 — A run with no "law" stream stops "converged" on pass 1 with maxDev = 0 and returns the solve under the DECLARED weights; prior sigmas are never estimated — REPRODUCED (control flow + arithmetic)
What is wrong: `pvMax = 0` (`adjustVCEjoint.m:476`) and `groups = {}` when no column has policy "law",
so `pv.maxDev = 0 <= Tol` at `:599`. The returned `solF` (`:372`) was accumulated BEFORE the pass's level
update; `priorFor` only updates priors from `lastSol` on pass ≥ 2 (`fit.m:681`), which never comes.
Every per-row-σ stream defaults to "scale" (`fit.m:586`), so demo_lsqvce case 5 is exactly this.
Observable symptom (`t13_scaleonly.py`, demo_lsqvce case-5 geometry, group 2 truly 3× noisier):
`converged = True` at pass 1; returned σ₀ **2.311** (the demo prints it under "variance ratio back at ~1",
`examples/demo_lsqvce.m:153–155`); formal se 1.66× the re-weighted ones; the Manual reconstruction
gives **2.42 / 7.16** for truths 1 / 3; a `data.priors(k).estimate="scale"` group keeps its declared σ.
Also general (L1): the returned solution is one pass behind the returned noise model, and levels, arc
factors and prior sigmas are outside the criterion (theory ch. 2 §"Convergence criteria" promises
"every free component").
Severity: S1.

### F7 — The robust pass biases the level low; "UNBIASED … the weights cancel in α²" is false: s₀ does not see the weights — REPRODUCED
What is wrong: in `sweepArcs>oneArc` the data are descaled by the inflated Cep but Σ_t = `a.Cfcn` is the
unweighted temporal law (`adjustVCEjoint.m:898–901`); s₀ contains tr(Σ_t⁻¹), which no weight touches,
while a down-weighted epoch's residual leaves Ω. Claim at `adjustVCEjoint.m:96–101`, RELEASE_NOTES 6.1,
Manual §7.8, test_vcePolicy:826–837.
Observable symptom: zeroing 20 of 200 epochs changes s₀ **242.93 → 242.91** (inline run). Clean AR(1)
0.5 stream, IGG3, 9×256: level with the law fixed **0.71–0.78**, lvlB = √level·σ₀ **0.83–0.86**
(`t3_robust_level.py`); law free **0.93/0.93/1.01/0.93** (`t11_robust_law.py`). The oracle's 15 % band
(lvlB ≥ 0.85) cannot discriminate; one MATLAB reading of 1.007 is within this spread. Second symptom
of the same per-STREAM MAD: a legitimately 3× noisier arc is 17.6 % zeroed by IGG3 (mean factor 0.63,
`t20_arcs.py`) — the robust pass and the arc factors compete for the same variance.
Severity: S1 (moderate magnitude, documented as unbiased).

### F8 — The stop criterion contains α²₀, which ZeroFreqRule "lowest" overwrites, so a settled run never converges — REPRODUCED
What is wrong: `psdVCE.m:241–242` sets Ŝ₀ = Ŝ₁ but `maxDev` (`:256`) includes |α²₀−1|. Unless the true
S₀ = S₁, α²₀ can never reach 1.
Observable symptom (`t14_dc.py`, started AT the true law, φ 0.5, 9×256): |α²₀−1| creeps 0.009 → 0.025 over
25 passes while max_{j≥1} sits at 0.007; with "keep": 0.001 / 0.001. Shipped default start: maxDev floor
0.023, all of it j = 0 (`t4_closedloop.py`). Result: `stopReason "maxPasses"`, `notConverged` warning and
`NOISE_UNSETTLED` on a settled model.
Severity: S2.

### F9 — `Noise.Smooth` and `Noise.FreezeRedundancy` are documented, defaulted and reach nothing — REPRODUCED (grep)
`grep -rn "Smooth\|FreezeRedundancy" +lsqvce` → only `Options.m:99,100,313`; `engineArgs('adjustVCEjoint')`
(`Options.m:256–261`) forwards neither; `psdVCE` is never passed `'Smooth'`. Manual 04_options.md:115–116
promises "a frequency whose redundancy is below this holds its value"; psdVCE only skips s̃_j = 0 exactly.
Severity: S3.

### F10 — An unseeded joint run records no probe seed — contradicts §12.2 ("probe contract + seeds in provenance") — REPRODUCED (grep)
`sharedProbe.m:115` discards `probeMatrix`'s `seedUsed`; `fit.m:326` records `jointSeed = opts.Seed` (`[]`).
Each pass draws a fresh private seed, none is recoverable, so `Seed=[]` runs are not replayable.
(`vce.m:187` pins one; the joint driver does not.) Severity: S3.

### F11 — Ljung–Box: K = 20 regardless of arc length / Grid, and pooled lag sums mis-size it on short arcs — REPRODUCED
`adjustVCEjoint.m:1100–1101` uses Q = m(m+2)Σρ²/(m−k) with m the TOTAL count, but the lag sums are per
arc (each loses k terms). White null, 2000 trials: 9×256 size 0.040, 8×2160 0.046, **64×30 size 0.000
(E[Q] 12.8 vs 20)**. At 2160-epoch arcs the test is blind beyond 20 lags (0.9 % of the Grid Ellmer uses).
Severity: S3.

### F12 — Stale claims inside the estimator's own text — REPRODUCED (read/grep)
`adjustVCEjoint.m:14–17` "ADJUSTVCEJOINT alternates … the SEPARABLE model, and it is the toolbox DEFAULT"
(self-referential; the separate mode was deleted in M4 4c); `:1049–1051` warning text "never fed back into
the weighting (finding F-01) … Take the ONE-PASS estimate" vs its own help (`:1027–1028`) "WP-98 fixed
that"; `Results.m:98` "use the result of adjustVCE". (L6 residue.) Severity: S4.

## Checked and sound
- redundancyLags s_n, Ω_n (Toeplitz Gohberg–Semencul path, local block K, FFT lag sums): equal dense tr(R̄V_n), ê'Σ⁻ᵀV_nΣ⁻¹ê to 1.3e-16 / 2.2e-16 (`t1_kernels.py`).
- Which matrix is traced: s₀ = tr(R) = 220.145, not tr(H-part) = 12.64; E[Ω_n] (MC 4000) = s_n within MC error (unbiased under the model).
- Hutchinson via sharedProbe normalisation: mean 220.14594 vs exact 220.14478; sd 0.0221 vs closed form 0.0224 at Z = 100; summed s̃_j relative sd ≤ 1.1 % at P = 200 (inline run). Probes shared across arcs per pass (read).
- DCT-I: 'ortho' self-inverse (2e-14), equals scipy `dct(type=1, norm='ortho')` (1.5e-14); 'ellmer' explicit inverse 4.6e-15, cond 1.46 (the positivity defect F1 is the transform's semantics, not its arithmetic).
- psd2acv / fNorm→f placement: S·dt with f = fNorm/dt leaves r(0) invariant; Grid = m reproduces C to 1e-8 (the ridge) (`t5_decl.py`).
- Multi-axis: descale + per-axis temporal whitening equals A'[D(R_t⊗I_d)D]⁻¹A to 7e-16 (axis-major I_d⊗R_t differs 18 %); per-axis s_n, Ω_n sums equal tr(R(V_n⊗I_d)) to 3.5e-16 (`t18_kron.py`).
- Eq. 6.5.28 after WP-98: α²_m invariant to the arc's own s2 (exact cancellation, linear dctI); γ² uses the stream level the whitener applied; level×shape enters the whitener once (`sigmaScale` → `cm.sigma` → S·σ²), C normalised to C(1)=1 in lawSpecs — no double count on a single stream.
- Koch–Kusche (`fit>priorUpdate`): N⁻¹ = Cxx/s0sq is the prior-INCLUDING cofactor (s0sq cancels); E[Ω_k] = 26.63 vs r_k = 26.73 at the true σ; iteration from 0.01 → 0.276 (truth 0.30, one realisation) (`t21_kk.py`).
- Robust standardisation is level-invariant (MAD per stream on Cep-descaled raw residuals; pass-0 sweep; folded factors divided back out) — read, consistent with §12.29.

## Not reviewed
- MATLAB execution of any claim (tier B); mirrors use exact traces, one system, no local parameters, Gamma = 1 unless stated — probe noise and locals would only add variance.
- Skip handling bias and SkipThreshold; C0 dependence (L4) beyond start-at-truth; Pool "obsType"/"all" with several law groups; `psdMode` "ar"/"arma" (levinsonDurbin truncation path read only); streamContribution, vce.m, redundancy.m, applyRegularisation (read; not on the joint user path except applyRegularisation's local `Reg="none"` default).
- adjustNonlinear (L2/L3) and adjustRobust (Yang) — other lanes.
