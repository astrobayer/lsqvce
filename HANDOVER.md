# lsqvce v2 — Handover: Review, Redesign, and Work Packages

**THIS IS THE ONLY STATUS-AND-PLANS DOCUMENT for the v2 toolbox. Do not create
a second one.** The v1 state of this tree is FROZEN for feature work as of
2026-07-19; its own `HANDOVER.md` records that and remains the authoritative
history of v1. Everything about v2 — the design review that motivated it, the
architecture, the task list, and the implementation prompts — lives HERE.

Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Fable 5.
Last substantive update: 2026-08-18 — WP-30 ROUND 1 (container session,
Appendix-D pipeline): full deliverable set on branch `claude/wp-30`
(ci.yml, ci_gate/gateVerdict/test_ci_gate, run_all skip mechanism,
CLAUDE.md, ignore/attributes, references policy); §0 loop + §9.0
preamble amended per D4; ledger row IN PROGRESS; round note in §5 WP-30.
Preceding entry: 2026-08-17, night (v2) — APPENDIX-D ALIGNMENT:
the user's proven pipeline (BEST_PRACTICE_v4-1.md, now in references/)
makes the chat/Cowork CONTAINER the execution vehicle, authenticated by
a fine-grained PAT (Contents/PR/Actions RW, D.3 handling rules) — the
Auth note is corrected accordingly; D.4 session preflight pinned as
binding. Preceding entry: D8 POSTFLIGHT
amendment: order is CI green -> local MATLAB test via MCP (postflight)
-> only then the PR; auth note pinned (no access tokens anywhere — two
one-time browser grants). Preceding entry: D8 CORRECTED per user:
PUBLIC repo (free hosted MATLAB, license choice now an open user
decision), prompting from chat/Cowork (no GitHub triggers, no @claude
action), optional local testing via the official MATLAB MCP Server;
guide rewritten again. Preceding entry (superseded within the hour):
WORKFLOW DECISION D8
(supersedes D6, absorbs D7): GitHub + hosted CI + @claude development
loop, chat/Projects/Cowork as thinking-and-trigger surfaces, transferable
pattern; WP-30 rescoped; `SETUP_GUIDE.md` rewritten for D8.
Previous update 2026-08-17, evening — WORKFLOW DECISION (D6):
university GitLab as plain hosting, local gate + clean-clone verify +
pre-push hook; no CI, no runner, tier 3 off. WP-30 rescoped and renamed;
`SETUP_GUIDE.md` replaces `SETUP_ROUTE_C.md`.
Previous update 2026-08-17 — SINGLE-FOLDER RESTRUCTURE (user
decision): `lsqvce2/` dissolved; THIS document moved to the toolbox root
as `HANDOVER.md`; the frozen v1 handover renamed `HISTORY_v1.md` (banner
added); all paths in this document, the WP-30 YAML and `SETUP_ROUTE_C.md`
swept accordingly; audit + v10review re-run clean after the move.
Previous update 2026-08-16 — REVIEW R2 added (§8 findings, §9
prompts, WP-32…WP-41) and the §1 ledger reconciled to the user's status
statement of the same day; numerical evidence shipped as
`validation/v10review.py`. Session model: Claude Fable 5.
Previous update 2026-07-19, late session — EXTENSION TRACK adopted
(user decision): nine new work packages WP-23…WP-31 plus fold-ins into
WP-7/9/15/16/17/18/20 and §6.1.2. Adoption record: §3.8; ledger rows: §1.
Correction, same session: WP-30 (CI) DEFERRED — no MATLAB-capable CI
runner exists in the user's environment; the local confirming run (§0
loop) remains the only gate until one does.
Also pinned in the same session: the PLOT-HOOK contract (user
requirement: results after EACH step + the processingOptions) —
per-pass firing, ctx payload, containment. See WP-16(c), WP-3, WP-18,
§3.3 / §3.5 / §3.7.

Document map:
- §0  Standing project rules (carried over verbatim in substance) + the loop
- §1  STATUS — the work-package ledger. **This is what is outstanding.**
- §2  Part I — THE REVIEW. Three-expert scrutiny of v1: Ellmer fidelity,
      design flaws, superfluous paths, method limitations.
- §3  Part II — THE REDESIGN. Principles, architecture, data contract,
      configuration, the one outer loop, parallel plan, old→new mapping,
      extension track (§3.8).
- §4  Part III — THE FUNCTION SET. The complete v2 public surface.
- §5  Part IV — WORK PACKAGES & PROMPTS. One copy-paste prompt per package.
- §6  Part V — GUIDELINES. Tests, documentation, coding standards. Every
      prompt in §5 references these; they are stated ONCE, here.
- §7  Appendix — review evidence: file/line pointers into v1.
- §8  REVIEW R2 — the 2026-08-16 blind five-expert scrutiny of the ONE
      tree: verdicts, coverage, findings F-01…F-14 with measured evidence.
- §9  R2 WORK PACKAGES & PROMPTS — one copy-paste prompt per package,
      WP-32…WP-41.

---

## 0. STANDING PROJECT RULES

These do not expire between chats. They are the same rules v1 was built under
(v1 `HANDOVER.md` §0a), restated so this document is self-contained.

1. **Role.** Expert MATLAB programmer and least-squares-adjustment specialist.
2. **Think step by step and MOTIVATE every answer** — why, not just what.
3. **SINGLE POINT OF ACCESS where useful.** One name per thing, one place per
   decision. v2 exists largely to enforce this (see finding F1).
4. **Always prepare unit tests.** New behaviour ships with a test that would
   fail without it. Tests cover CONTRACTS, NUMERICAL PRECISION, and SPEED (§6.1).
5. **Every function gets full help text**: purpose, INPUT and OUTPUT with TYPE,
   DIMENSION and DESCRIPTION, ERRORS it raises, an example, and a
   `Created <date> (CEST, UTC+02:00) with assistance from Claude (Anthropic),
   model <version>` stamp placed after `See also`.
6. **Toolbox documentation is HTML reachable from `doc`** — a page in `html/`
   plus entries in `Contents.m`, `helptoc.xml`, `html/index.html`. **No
   JavaScript/Java anywhere**, equations as static images (v1 build tooling is
   reused, WP-21).
7. **Keep THIS handover current, clean and tidy.** The ledger in §1 is the only
   place a work-package status lives.
8. **No aliases.** A renamed thing has ONE name; the old one is deleted.
9. **Follow the MATLAB coding standards**
   (https://github.com/matlab/rules/blob/main/matlab-coding-standards.md).
10. **Always provide a WRAP after changes.** The user downloads the zip, runs
    the suite in MATLAB, and returns the log AS AN UPLOADED .txt FILE (large
    pastes are lost in transit). *(Amended 2026-08-18, WP-30 round, per the
    D4 mapping: from WP-30's merge on, CI GREEN IS THE GATE — `ci.yml` runs
    `tests/ci_gate.m` on every push — and the local MATLAB run remains (a)
    the postflight second platform and (b) the run for suites CI skips,
    currently `test_parallel`. The zip/log loop stays valid for chat rounds
    without repo access. Amended 2026-08-18 (user decision): a round that
    ends with a READY PULL REQUEST ships NO zip — the PR is the deliverable;
    zips remain for rounds without a merge-ready PR.)*

### The working loop (identical to v1 §4b)
- There is no MATLAB in the sandbox. Validation is STRUCTURAL (function/end
  balance), NUMERICAL (Python vs exact formulae), and STATIC (the audit
  script, ported in WP-0). **The user's MATLAB run is the ground truth.**
- Ship the ENTIRE tree as ONE zip per round via `present_files` — never a
  partial file set.
- Run the full suite with the ONE runner `run_all_lsqvce`; from WP-30
  (branch pushed 2026-08-18) the CI/postflight entry is `tests/ci_gate.m`,
  which CALLS that same runner and adds the coverage-and-log verdict — one
  implementation, one gate on top; read the log,
  not just the pass count: check what LOADED, the warning inventory, and the
  per-suite seconds.

### How to work through this document
Take the next OPEN row in §1 — ORDER OF WORK is the dependency column, not
the numbering (see the ledger note) — open a NEW chat, paste the PROMPT PREAMBLE
(§5.0) followed by that work package's prompt (§5.WP-nn), and attach the
current zip. When the round comes back green, mark the row DONE with the date
and the confirming run, and record anything learned in that WP's section.
Dependencies are listed per package; do not start a package before its
dependencies are DONE.

---

## 1. STATUS — THE WORK-PACKAGE LEDGER

Suite state: the ONE `lsqvce/` tree, 497/498 green (2026-07-19 run, see
`HISTORY_v1.md` §1); R2 packages below are the outstanding work.

**RECONCILIATION (2026-08-16, user decision).** There is ONE toolbox: the
`lsqvce/` tree. The refinement this document planned was realized as
capabilities IN that tree rather than as the renamed §4 surface; the user
confirms all development except WP-30 (CI) and WP-31 (GPU) is implemented.
The rows below are marked DONE on that confirmation, with the tree's own
suite record (497/498 green, 2026-07-19) as the standing MATLAB evidence;
`HISTORY_v1.md` remains the per-round history and §2–§7 stay as the
design record. OUTSTANDING now: WP-30 (DEFERRED), WP-31 (GATED), and the
R2 packages WP-32…WP-41 opened by the second review (§8/§9).

| WP | title | depends on | status |
|----|-------|-----------|--------|
| 0  | Skeleton, setup, test runner, audit script | — | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 1  | `util`: log, verbosity, options | 0 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 2  | `ParamSet` — parameter tables with scopes | 0,1 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 3  | `dataSource` — the user data contract + `checkSource` | 0,1,2 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 4  | Signal kernels: `dctI`, `psd2acv`, `levinsonDurbin`, `armaACV`, `ar2arma` | 0,1 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 5  | `toeplitzWhiten` — exact streaming Toeplitz whitener (NEW) | 4 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 6  | `noiseModel` — spec → whitener dispatch | 3,4,5 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 7  | Backends: `QRSystem`, `NormalEqSystem` (internal) | 1 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 8  | `Assembler` — scoped elimination engine | 2,3,6,7 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 9  | `mapreduce` + prefetch | 1 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 10 | `sweep` — the ONE second data pass with collectors | 8,9 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 11 | `sharedProbe` + Hutchinson trace | 8 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 12 | `vce` — group Helmert components (exact \| hutchinson) | 10,11 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 13 | `ellmerUpdate` — PSD amplitudes + arc factors (workhorse) | 10,11 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 14 | `armaUpdate` — parametric alternative | 10 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 15 | `snoop` — block-native screening + downweighting | 10 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 16 | `iterate` — the ONE outer loop + convergence gate | 12,13,14,15 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 17 | `Results` + `report` + `contribution` — block-native | 10,16 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 18 | `adjust` + `configure` — the single front door | 16,17 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 19 | Template + GRACE-FO synthetic example | 18 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 20 | Validation campaign: Ellmer parity, precision, speed | 18,19 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 21 | Documentation: HTML set, manual, Contents, helptoc | 18,19 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 22 | Final audit, freeze, release notes | 20,21,23–29 (+30/31 if opened) | DONE (2026-08-16, user-confirmed; see reconciliation) |

Extension track (adopted 2026-07-19, user decision — every review-panel
suggestion except the temporal prior between days; adoption record: §3.8):

| WP | title | depends on | status |
|----|-------|-----------|--------|
| 23 | Component uncertainty (Fisher information) | 12,13 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 24 | VCE-determined regularisation (Koch–Kusche) | 8,12,16 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 25 | Robust IRLS — alternative to snooping | 10,15,16 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 26 | `orderSelect` — AIC/BIC/CV model order | 13,14,16 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 27 | `neqExchange` — SINEX + GROOPS NEQ I/O | 2,8 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 28 | Sketched diagnostics (randomised leverage) | 11,17 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 29 | Metamorphic test suite | 18,19 | DONE (2026-08-16, user-confirmed; see reconciliation) |
| 30 | Development workflow: public GitHub + hosted CI + cloud-session loop (D8) | — | DONE (2026-08-18, PR #1, merge 6ddd415; CI green per commit incl. both PR contexts, red test passed, postflight green via MCP bridge on Windows R2026a 513/513; branch protection on `main` requiring `audit`+`matlab` set by the user the same day; evidence + deviations in §5 WP-30 round notes) |
| 31 | GPU backend | 7,20 | GATED — open only on WP-20 evidence |

Rules for the ledger: one row per WP; status is OPEN → IN PROGRESS → DONE
(date, confirming MATLAB run); nothing else is in flight unless it has a row.
ORDER OF WORK is the dependency column, not the numbering. Recommended
slots: WP-27 before WP-20 (GROOPS parity becomes a file diff); WP-22
always last. WP-30 is DEFERRED (user, 2026-07-19: no MATLAB-capable CI
runner — licence or server — is available; only the local loop works);
the moment one exists, open WP-30 right away, whatever else is in flight
— it guards every round thereafter.
Checkpoint/restart, the run journal, and deterministic reductions have NO
rows — they are folded into WP-16 and WP-7/WP-9 respectively (§3.8).
Default slot for every extension package (23–29, 31): AFTER WP-19, so
the standing campaign fixture their tests cite exists; stated
exceptions override (WP-27 anywhere after 8, recommended before WP-20;
WP-30 deferred).


R2 track (opened 2026-08-16 by the second review, §8; prompts in §9;
order of work: 32 → 33 → 34 → 35 → 36 → 40 → 37 → 38 → 41 → 39):

| WP | title | depends on | status |
|----|-------|-----------|--------|
| 32 | Fix `vce` 'hutchinson': redundancy, not trace (F-01) | — | DONE (2026-08-18, PR #2, merge e167cee; CI green both contexts after a setup-matlab cache-restore hang was resolved per the incident runbook, postflight green 515/515 on Windows R2026a; evidence + round note §9 WP-32) |
| 33 | Joint paths: locals in residual + projector (F-03) | — | IN PROGRESS — round 1 MERGED (2026-08-18, PR #3, 7d7c14c; CI + pre-outage postflight green). ROUND 2 (2026-08-19, user decision) on `claude/wp-33-r2`: the deferred coloured-synthetic PSD-recovery test (corr > 0.95 with locals) + the single-run CI trigger scheme (user decision, see §0/§5 amendments); D8 loop pending |
| 34 | Keyed arc identity in `adjustVCEjoint` (F-02) | — | OPEN |
| 35 | `psdVCE` floor policy + negative sigma2 guard (F-04, F-05) | — | OPEN |
| 36 | One probe generator, honest RNG contract (F-06) | — | OPEN |
| 37 | `redundancyLags` without dense W^-1 + hot-loop tolls (F-09) | — | OPEN |
| 38 | Deduplicate keys / arc batching + audit lock (F-08) | 33,36 | OPEN |
| 39 | `arguments` blocks + Name-Value casing, whole surface (F-07) | 32–38,40,41 | OPEN |
| 40 | Method polish: gates, clamps, whiteness (F-10…F-13; + F-15 alternation bistability, registered 2026-08-19) | 32 | OPEN |
| 41 | Documentation truth sweep (F-12, F-14) | 32,33,36,40 | OPEN |

---

## 2. PART I — THE REVIEW

Conducted 2026-07-19 by a three-expert panel over the complete v1 tree
(80 `.m` files, ~22,400 lines, suite green 497/498) with the Ellmer (2018) and
Kvas (2019) theses open beside the code:

- **Expert A — GRACE/GRACE-FO processing.** ITSG-style pipelines, GROOPS,
  the Ellmer §6.5 stochastic model, Kvas near-real-time processing.
- **Expert B — Adjustment theory & statistics.** Helmert VCE, reliability
  theory (Baarda/Pope, redundancy, MDB), robust estimation.
- **Expert C — Software architecture & numerics.** MATLAB API design, HPC and
  the Parallel Computing Toolbox, numerical linear algebra.

Verdict up front, agreed unanimously: **the numerical kernels of v1 are good
and largely verified — the architecture around them is not.** v1 grew three
entry points, four accumulation pipelines, two diagnostic chains and two outer
loops by accretion; the Ellmer machinery is complete but demoted to an opt-in
mode; and the two features the target workload needs most — data snooping on
the block path and arc-scoped elimination — are missing. The right move is a
clean v2 that KEEPS the verified kernels (they are listed per work package)
and rebuilds the structure around one loop, one sweep, one front door.

### 2.1 Method fidelity — deviations from Ellmer (2018), and their verdicts

Each deviation: what v1 does, what Ellmer does, whether the deviation is
justified, and what v2 does about it.

**D1 — The Ellmer model is not the default.**
v1's default `vceMode` is `'separate'`: alternate a per-stream AR/ARMA fit
(shape) with a Helmert group VCE (scale), gated by Ljung-Box whiteness.
Ellmer's joint co-estimation (§6.5.7) exists (`adjustVCEjoint`,
`'vceMode','joint'`) but is opt-in. The project brief states the Ellmer
approach is "the primary goal and workhorse". *Verdict: not justified as the
default.* The alternating scheme is a legitimate ALTERNATIVE (it has a clean
stopping rule and a parametric model that extrapolates beyond the data), but
the roles must be inverted. **v2: `cfg.noise.method = 'ellmer'` is the
default; `'arma'` and `'var'` are the documented alternatives (WP-13/14/16).**

**D2 — The white start.**
Ellmer eq. 6.5.62 prescribes S^j = 1, sigma2_m = 1 (white). v1 refuses the
literal delta autocovariance start and defaults to a mildly correlated C0,
because a Python study (v1 HANDOVER §2 "KEY FACTS") measured that the delta is
a DEGENERATE FIXED POINT of Ellmer's literal update: his per-lag quantities
Omega_n, s_n (eqs. 6.5.22/6.5.23) carry the multiplicative prefactor
sigma2_m * C^n_xx, so with C^n = 0 for n > 0 every lag contribution beyond the
zero lag vanishes, the DCT of a spike is flat, all alpha2_j come out equal,
and a flat spectrum can only ever be RESCALED, never reshaped. *Verdict: the
deviation is real and justified — but it treats a symptom.* The panel's
sharper diagnosis: the degeneracy is a property of defining the variance
components PER LAG (component matrices C^n * V_n) and transforming afterwards.
Define them PER FREQUENCY instead — Sigma = sum_j S^j * T{x_j}, with x_j the
j-th inverse-DCT basis vector, exactly as Kvas (2019, eq. 3.14, after Koch ch.
36) writes it — and the multiplicative Helmert update acts on S^j directly:
every component starts at 1 > 0, none is pinned at zero, and the white start
is legitimate as written. In the spectral formulation the estimator becomes
alpha2_j = (DCT of the unweighted lag sums e~'V_n e~) / (DCT of the unweighted
lag traces tr(Rbar V_n)) — the same time-domain symmetries and the same shared
probe apply; only the C^n prefactor moves out of the lag sums and into the
prior S^j where it belongs. **v2: the spectral-domain component definition is
the default (`'components','spectral'`); `'components','ellmer'` reproduces
the thesis bit-for-bit (including the non-degenerate-C0 guard v1 already
carries) for GROOPS/ITSG parity runs (WP-13).**

**D3 — `numLags` vs `nf` and the rank rule.**
v1 distinguishes the estimator's lag count from the PSD sampling density and
enforces `nf >= m/2 + 1` because its `covModel` precomputes ONE autocovariance
vector per stream sized for the BATCH window and reuses it for every batch.
In Ellmer the two coincide by construction: Nmax is the longest arc and the
DCT-I maps Nmax lags to Nmax frequencies. *Verdict: a sizing artifact of a v1
precompute, correctly documented, but it confuses users and it exists only
because the covariance is keyed to the batch, not the arc.* **v2: the noise
model is keyed to the ARC family (Nmax = longest arc of the stream, lags ==
frequencies == Nmax), which restores Ellmer's structure exactly and deletes
the distinction. The rank rule remains as an internal assertion (WP-6/13).**

**D4 — Batching inside an arc approximates the Toeplitz whitening.**
Ellmer decorrelates WHOLE arcs (chol of the dense Nm x Nm Toeplitz). v1
splits an arc into row batches; for AR/ARMA filter whitening the
overlap-and-discard scheme is EXACT (verified, 1.6e-15), but for the dense
`'toeplitz'` path each batch is whitened by chol of its own submatrix with
`filterMemory = 0` — i.e. a BLOCK-DIAGONAL approximation of the arc
covariance; with a manual overlap it becomes a truncated-innovations
approximation. Neither is flagged as an approximation in v1, and the project
brief explicitly anticipates "arcs split in batches ... limited by the
available RAM and a fully populated Toeplitz-like covariance matrix".
*Verdict: an UNJUSTIFIED silent approximation on exactly the path the
workhorse method uses.* **v2 fixes it properly: a streaming Toeplitz whitener
(WP-5) built on the Levinson/Schur recursion. The k-th whitened row is the
k-th INNOVATION e_k = (l_k - prediction from rows 1..k-1)/sqrt(PEV_k); the
recursion delivers the order-k predictor progressively in O(N) memory, the
predictor converges to the stationary filter (detectable, so the per-row cost
drops to the effective filter length), and the result is EXACT for any batch
size — batching becomes a pure memory decision with no accuracy trade.**
Fallback: dense chol per whole arc when N is small. An optional additive
per-arc covariance term (Ellmer ch. 8 AOC-style Sigma = Sigma_stationary +
Sigma_add) forces the dense path and is supported as `'addCov'` (WP-6).

**D5 — Arc-wise variance factors, gamma normalisation, robust mean.**
v1 implements eqs. 6.5.27–6.5.30 with the interquartile robust mean of
eq. 6.5.57 (`'MeanType','robust'`, default) and the plain mean as an option.
*Verdict: faithful; keep as is (ported in WP-13).*

**D6 — Zero-frequency rule.** v1 `psdVCE` implements Ellmer's guard (set S^0
to the lowest non-zero-frequency amplitude) as `'ZeroFreqRule','lowest'`.
*Verdict: faithful; keep (WP-13).*

**D7 — Hutchinson probes and the shared block.** v1 implements eqs.
6.5.39–6.5.48 including the once-per-iteration U^-1*Zbar reuse and Z = 100
default; the redundancy-seam operator-order bug found in v1 is FIXED and
regression-tested. *Verdict: faithful and verified; port (WP-11).*

**D8 — Reduced-resolution noise iteration.** Ellmer iterates the stochastic
model on a d/o-60 field and computes the final solution at full resolution
(§6.5 intro). v1 iterates at full resolution always. *Verdict: v1 is not
wrong, but it forgoes Ellmer's principal cost control on precisely the loop
that dominates runtime.* **v2: `cfg.iterate.reducedCommon` — a column
selector (or reduced ParamSet) used during the noise iteration, followed by
ONE full-resolution solve under the converged model (WP-16).**

**D9 — Group weighting can double-count under the Ellmer model.**
v1's outer loop ALWAYS runs a Helmert group VCE (scale) beside the shape
update. Under the Ellmer model the PSD amplitudes already carry the absolute
scale of each observation group, and the arc factors are normalised to mean 1
precisely to avoid double counting (eq. 6.5.28–6.5.30). Running an
unconstrained group VCE on top re-introduces the coupling the gamma
normalisation removes — v1 itself documents that joint and alternating
updates "can converge to DIFFERENT models". *Verdict: a genuine
double-counting hazard.* **v2: with `'ellmer'` noise, per-group scale factors
default OFF (the PSD owns the scale); they can be enabled for the a-priori
groups Ellmer himself scales with a single factor (the pseudo-observation /
prior group, eq. 6.5.31), which v2 supports as the prior's variance component
(WP-12/13/16).**

### 2.2 Design flaws and superfluous paths

**F1 — Three "single" points of access with subset semantics.**
`adjust` ⊂ `adjustBlock` ⊂ `adjustVCE`, each with a DIFFERENT data contract
(flat `batchFcn(k)`, loader `(day,obsType,arcID,opt)`, fluent builder). The
flat world is the degenerate case of the block world (one system, zero
locals, one arc), so the hierarchy is capability duplication, not layering.
**v2: ONE front door, `lsqvce.adjust(cfg)` (WP-18). A single GLS pass is
`cfg.iterate.maxOuter = 0`; a flat problem is a source with no locals.**

**F2 — Four accumulation engines.** `Accumulator` (flat), `Combine`
(block-arrow over batch closures), `BlockAccumulator` (block-arrow over
day/arc loaders), plus the `QRSystem`/`NormalEqSystem` backends they share.
`Combine` and `BlockAccumulator` implement the SAME mathematics over
different contracts. **v2: one `Assembler` (WP-8) over one data contract;
the two backends survive as its internals (WP-7). `Accumulator` and
`Combine` are deleted.**

**F3 — Data snooping lives on the wrong path.** `snoop` (iterative rejection
with verified Cholesky downdate and the correct vPv/r_i update) works ONLY on
the flat path; the block path has `outliers` (flags + a mask handle, no
iteration, no downweighting, no wiring into any loop). The brief requires
"data snooping with subsequent downweighting, iterated to convergence" for
GRACE-scale data — i.e. on the block path. **v2: `snoop` (WP-15) is
block-native: it consumes the sweep's whitened residuals and block-exact
redundancies, tests in the whitened metric, maps flags back to
(system, obsType, day, arc, epoch), and applies a persistent action —
`'reject'` or `'downweight'` (Huber/IGG3 variance inflation folded into the
next pass's whitener) — inside the ONE outer loop, one wave per pass, with a
redundancy floor, a rejection cap, and a no-oscillation policy. The flat-only
`snoop`, and `robust` as a separate flat path, are absorbed.**

**F4 — Two diagnostic chains on the block path, and the default one does not
scale.** The block-native chain (`residuals`/`redundancy`/`outliers`) is
exact and streaming. But `Results.report()` routes through `diagnoseAccum`,
which REBUILDS a dense QR factor of the FULL parameter set
(u_full = nc + sum nl). For the target workload (nc ~ 1e4 at d/o ~90–100,
thousands of arc-locals per month) that is a multi-GB dense factor — it
defeats the pre-elimination the engine exists to provide. **v2: ONE
diagnostic chain, block-native. The five-level report (a v1 strength worth
keeping) is fed from the sweep's exact per-observation r_i/w_i and the
REDUCED common system; `diagnoseAccum` and the flat `diagnose` disappear
(WP-10/17).**

**F5 — Elimination granularity.** v1 pre-eliminates locals PER SYSTEM after
the whole campaign is accumulated: Nll is stored DENSE across every arc of
every day of the system (block-diagonality is exploited only inside the
final `preeliminate`). Memory grows with (campaign length)^2 per system, and
nothing can be freed early. Real GRACE processing eliminates arc state
parameters when the arc closes. **v2: parameters carry a SCOPE — `common |
system | day | arc` — declared in the ParamSet (WP-2); the Assembler (WP-8)
eliminates each scope when it CLOSES (arc-locals at end of arc, day-locals at
end of day, ...), keeping the Schur blocks per scope for back-substitution
and residuals. Peak memory becomes O(nc^2 + max_scope nl^2). The v1 C/D-pair
case (range-rate couples both satellites' arc states) is expressed naturally:
those locals are arc-scoped within the pair system.**

**F6 — `adjustVCE` is a 2,479-line god class.** Fluent configuration, the
outer loop, joint mode, plot hooks, Ljung-Box, ACV helpers and half the log
formatting in one file. **v2: the loop is a ~200-line function (`iterate`,
WP-16) that composes named strategies; configuration is data (`configure`,
WP-18); presentation lives in `report`/`Results` (WP-17).**

**F7 — Two outer loops.** `iterateNoiseModel` and `adjustVCE.run` iterate
the same fixed point with different gates. **v2: one loop (WP-16).**

**F8 — Three names for Helmert VCE.** `vce`, `varComp`, `vceHutchinson`
(plus `traceEstimate` and `streamContribution` overlapping underneath).
**v2: one `vce` with `'trace','exact'|'hutchinson'` (WP-12); the streaming
contribution primitive is internal to the sweep (WP-10/11).**

**F9 — Redundant data sweeps.** `residuals`, `redundancy`, `redundancyLags`,
`contribution` and the VCE inner iterations each re-stream the data. On the
target workload the second pass is where the hours go (v1's own estimate:
~4 h/iteration serial at d/o 90). **v2: ONE `sweep` per outer pass (WP-10)
with pluggable COLLECTORS (residual series, per-group vPv and redundancy,
per-lag Omega_n/s_n, per-observation leverage, contribution probes, snooping
statistics) filled in a single pass over the arcs, parallel over arcs via
`mapreduce`, with the shared probe broadcast once.**

**F10 — Stochastic state is threaded as stringly-typed option soup.**
`opt.covOverride` in two shapes (flat and `.bySys.(satKey).(obsType)`),
`satKey` via `makeValidName`, `jointModelOpt` existing solely to encode
"three conventions that otherwise fail SILENTLY". **v2: the estimated noise
state is a first-class per-stream `noiseModel` table owned by the loop and
handed to the Assembler explicitly; loaders never see it (WP-6/8/16).**

**F11 — The 14-column parametrisation cell.** Meaning by convention, shape
checked only, column count hardcoded in the toolbox contract although it is
the USER'S convention. **v2: a MATLAB `table` with a small set of REQUIRED
variables (`name`, `scope`, `system`) and free additional columns; `ParamSet`
(WP-2) validates the required part and passes the whole table through to the
user's assembler untouched. nc/nl per scope are DERIVED from the tables —
the v1 "state nc once" rule survives and generalises.**

**F12 — Surface area.** ~70 public names for what is one method family.
Small glue (`satKeyOf`, `normalizeSatList`, `resolveSatSpecs`,
`collectResiduals`, `whitenessCorr`, ...) is public and documented at the
same level as the engine. **v2: a compact surface — final count 44 (§4: the statistics annex,
the ported kernel set, and the 2026-07-19 extensions push past the early
~30 estimate, still well under half of v1); glue becomes
package-private. General-adjustment extras (`ftest`, `ellipse`, `datum`,
`stat`, `solveNonlinear`, `systematic`) are kept — they are sound and small —
but grouped as a clearly-labelled statistics annex, not the core.**

### 2.3 Limitations of the method (inherent — to document, not to fix)

**L1 — Stationarity per estimation window.** The Ellmer model assumes the
noise is wide-sense stationary within the estimation window (his: one month).
Sensor events break this; the remedy is splitting the window (Kvas §3.2.4
discusses the operational consequence). v2 documents the assumption and makes
the window an explicit configuration, not an implicit one.

**L2 — Arc independence.** No correlation across arcs (Ellmer eq. 6.5.8).
Matches the thesis exactly; v1 already respects it; v2 keeps it and states it.

**L3 — No cross-group correlation.** Sigma is block-diagonal across
observation groups (eq. 6.5.2). Cross-axis correlation WITHIN a stream is
covered by the `'var'`/`'separable'` models (Kvas VAR machinery); cross-TYPE
correlation (e.g. pos↔range-rate) is out of scope for the method.

**L4 — Monte-Carlo jitter.** With Z = 100 probes the estimated redundancies
carry O(1/sqrt(Z)) noise; harmless over iterations (Ellmer's own argument)
but it makes single-pass numbers non-reproducible unless the seed is pinned.
v2 pins a per-run seed by default and documents the trade.

**L5 — Signal leakage into the noise model.** Unmodelled signal and
non-stationary effects alias into the estimated PSD (Ellmer ch. 8 shows the
AOC case and its cure: additive known covariance). v2 supports the additive
per-arc covariance hook (D4) so the cure is expressible; detecting the
leakage remains the analyst's job (the report's whiteness and
factor-vs-geometry diagnostics help).

**L6 — Snooping consumes redundancy and degrades under weak geometry.**
Baarda-type tests are only meaningful where r_i is appreciable; with per-arc
state parameters the first/last epochs of an arc are weakly controlled. v2
keeps v1's refusal floor (`rMinReject`) and reports untestable counts
per arc, and tests in the WHITENED metric — under correlated noise the
per-observation test applies to the innovation, not the raw epoch, and the
report says so explicitly.

---

## 3. PART II — THE REDESIGN

### 3.1 Design principles

1. **One front door.** `R = lsqvce.adjust(cfg)`. Everything the run needs is
   in ONE validated struct produced by `lsqvce.configure`. Power users reach
   the layers below through `R` (the Results), never through parallel entry
   points.
2. **One loop.** A single `iterate` drives: assemble → solve → sweep →
   noise update → snooping → gate. Strategies (`ellmer`/`arma`/`var`, snooping
   on/off) are DATA selected in `cfg`, not separate code paths.
3. **One sweep.** All second-pass quantities come from one pass over the arcs
   per outer iteration, filled by collectors, parallel over arcs.
4. **Exactness by construction.** Whitening is exact for every model at every
   batch size (streaming Toeplitz whitener; filter overlap-discard; dense
   whole-arc fallback). Approximations that remain (Hutchinson trace) are
   explicit, seeded, and reported.
5. **Scoped elimination.** Parameters declare their scope; the Assembler
   frees each scope when it closes. Memory is O(nc^2 + max nl_scope^2)
   regardless of campaign length.
6. **The Ellmer model is the workhorse.** Default noise method, spectral-
   domain components, thesis-parity mode available.
7. **Small surface, deep help.** 44 public names (the §4 LIST is
   authoritative, not the early ~30 estimate); v1's help-text culture
   (motivated ERRORS sections, measured claims) carries over.

### 3.2 Architecture and data flow

    user code                 lsqvce v2 layers
    ---------                 -------------------------------------------
    assembleDay ─┐            L1 DATA      dataSource, ParamSet, checkSource
    epochCov ────┼─► sources
    param tables ┘                │
                                  ▼
                              L2 NOISE     noiseModel (white|psd|arma|var|
                                           separable|custom, +addCov)
                                  │  whiten(A,B,l) exact per arc
                                  ▼
                              L3 ENGINE    Assembler (scoped elimination,
                                           QR|normal, regularisation, prior)
                                  │  reduced common system, Schur blocks
                                  ▼
                              L4 SWEEP     sweep + collectors, sharedProbe,
                                           mapreduce (serial|parfor|spmd)
                                  │  residuals, r_i, Omega_n/s_n, vPv_g, ...
                                  ▼
                              L5 ESTIMATE  ellmerUpdate | armaUpdate, vce,
                                           snoop
                                  │  updated noise state, weights, masks
                                  ▼
                              L6 LOOP      iterate (gate; reducedCommon
                                           stage; final full solve)
                                  ▼
                              L7 RESULTS   Results, report (5 levels),
                                           contribution, convergence views
                                  ▲
    R = lsqvce.adjust(cfg) ───────┘        L8 FRONT DOOR  configure, adjust

Data granularity (unchanged in spirit from v1, formalised): a CAMPAIGN is a
set of SYSTEMS; a system streams DAYS; a day per OBSERVATION TYPE splits into
ARCS (irregular, gap-driven); an arc streams in ROW BATCHES (whole epochs,
epoch-major for multi-axis types). Arcs are the unit of covariance; batches
are a pure memory decision (principle 4).

### 3.3 The user data contract (renamed, simplified)

The user writes at most THREE things. `buildDay` is renamed `assembleDay` —
it assembles the observation equations of one day (verb + object; "build" said
nothing and "Day" alone hid that equations are the product).

1. **`assembleDay(day, obsType, ctx)`** → struct (REQUIRED)
       .A       m x nc   common design for the day
       .B       m x nl   local design, columns ordered per ctx.local
       .l       m x 1    observed minus computed
       .t       m x 1    epochs [s]
       .arcIdx  1 x nArc cell of row indices (optional; default: one arc)
   `ctx` carries `ctx.common` / `ctx.local` (the ParamSet tables), `ctx.proc`
   (free processing options), and injected fields (`day`, `obsType`,
   `system`). Determinism contract as in v1 (the sweep re-reads).
   One signature, always three arguments — v1's five-argument arity dispatch
   is gone (the tables travel inside ctx).

2. **`epochCov(day, obsType, ctx)`** → per-epoch scaling D (OPTIONAL)
   for separable multi-axis streams: `.arc(i).rows/.epochs`, `.sigma` or
   `.Cep` (nEpoch x d x d), exactly the v1 shape (it was right).

3. **Parameter tables** (ParamSet, WP-2): MATLAB tables, one row per unknown,
   required variables `name` (string), `scope` ('common'|'system'|'day'|'arc'),
   `system` (id or missing for common); any further columns are the user's
   own and are passed through untouched.

These are wrapped once:

    src = lsqvce.dataSource(@myAssembleDay, ...
              'system','C', 'observationTypes',{'pos','kbr'}, ...
              'noise', noiseSpecStruct, ...        % a-priori, per obsType
              'common', commonTable, 'local', localTable, ...
              'epochCov', @myEpochCov, 'days', dayList, ...
              'proc', procStruct, ...
              'plotPostSolve', {@myXPlot}, ...   % user hooks; contract:
              'plotPostModel', {@myPsdPlot}, ... % WP-16(c)
              'name','GRACE-C');

`dataSource` (WP-3) replaces `makeLoader` AND the raw-spec struct form — one
accepted shape, not two. It derives nc/nl per scope from the tables (the v1
"state it once" rule), provides the meta/arc protocol, the one-day cache, and
optional prefetch. `checkSource(src)` is a shipped contract validator the
user can run BEFORE any adjustment: determinism (double read, byte compare),
dimensional consistency, arc partition validity, epoch-major layout,
noise-spec resolvability. It doubles as the fixture for the contract tests.

### 3.4 The configuration struct — one source of options

    cfg = lsqvce.configure( ...
        'sources', {srcC, srcE}, ...        % or [src] — cell of dataSources
        'days',    days, ...                % default day list
        'solver',  struct('backend','qr', ...
                          'reg','none', ...       % 'none'|'fixed'|'gcv'|'vce' (WP-24)
                          'lambda',0), ...
        'noise',   struct('method','ellmer', ...     % 'ellmer'|'arma'|'var'|'none'
                          'components','spectral', ...% or 'ellmer' (parity)
                          'numLags','auto', ...       % auto = longest arc
                          'orderSelect','none', ...   % 'aic'|'bic'|'cv' (WP-26)
                          'probes',100, 'seed',1, ...
                          'arcFactors',true, 'meanType','robust', ...
                          'zeroFreqRule','lowest', ...
                          'groupBy','obsType', 'groupScale','auto'), ...
        'outliers',struct('method','snoop', ...      % 'none'|'snoop'|'irls' (WP-25)
                          'action','downweight', ... % snoop: 'reject'|'downweight'
                          'alpha',1e-3, 'bonferroni',true, ...
                          'rMinReject',0.10, 'maxPerPass',[], ...
                          'readmit',false, ...
                          'psi','huber', 'k',1.345), ...  % irls: psi + tuning
        'iterate', struct('maxOuter',8, 'tol',1e-2, ...
                          'reducedCommon',[], ...    % selector for D8 stage
                          'finalFullSolve',true), ...
        'run',     struct('parallel','serial', 'verbose',1, ...
                          'exportDir','', 'showSeconds',0, ...
                          'journal',true, ...        % JSONL per pass (WP-16)
                          'checkpointDir','', ...    % '' = off (WP-16)
                          'resumeFrom','', ...       % checkpoint to resume
                          'gpu',false), ...          % WP-31, gated
        'prior',   [] );                    % CommonPrior spec, optional;
                                            % a VCE group when reg='vce'

`configure` applies defaults, validates types/ranges/consistency ONCE (no
inputParser scattered over 30 functions), and returns a plain struct — the
"interface is a struct" requirement, made canonical. `adjust(cfg)` accepts
only this struct. Fluent builders are gone: they were a second way to say the
same thing (rule 3 in §0). The former `snoop` block is RENAMED `outliers`
(2026-07-19): one policy block hosting both the test-based wave (WP-15) and
the IRLS alternative (WP-25) — motivation in §3.8.

### 3.5 The one outer loop (pseudocode, WP-16)

    function R = iterate(cfg)
        state = initNoiseState(cfg);            % per-stream noiseModel, all prior
        mask  = emptyScreening();
        P     = reducedProblem(cfg);            % D8: reduced common selector
        for it = 1:cfg.iterate.maxOuter
            asm = assemble(P, state, mask);     % L2+L3: whiten, scoped elim
            sol = asm.solve();                  % reduced common solve, chol U
            firePlots('postSolve', ...);        % user hooks, per stream (WP-16c)
            S   = sweep(asm, sol, collectorsFor(cfg));   % L4: ONE pass
            state = updateNoise(cfg, state, S); % ellmer|arma|var (L5)
            firePlots('postModel', ...);        % user hooks, per stream (WP-16c)
            [mask, chg] = screen(cfg, mask, S); % cfg.outliers: snoop|irls (L5)
            g = gate(cfg, state, S, chg);       % see below
            record(history, it, state, S, g);   % + journal + checkpoint (§3.8)
            if g.converged, break, end
        end
        if cfg.iterate.finalFullSolve && wasReduced(P)
            asm = assemble(fullProblem(cfg), state, mask);
            sol = asm.solve();  S = sweep(asm, sol, finalCollectors(cfg));
            firePlots('postSolve', ...);        % both stages once more,
            firePlots('postModel', ...);        % ctx.final = true
        end
        R = packageResults(sol, S, state, history, cfg);

    GATE (all must hold):
      * max|alpha2_j - 1|            < tol      (Ellmer, eq. 6.5.19 sense)
      * max|sigma2_g - 1|            < tol      (only if group scale active)
      * relative PSD / model change  < tol
      * screening stable: no new rejections/downweights this pass
        (outliers 'irls': max weight change < tol instead)
    plus maxOuter cap with the v1-style actionable non-convergence warning.

    maxOuter = 0 (the single-GLS-pass mode, WP-18 — replaces v1
    adjustBlock): the loop body is skipped; ONE assemble + solve + sweep
    at FULL resolution with the a-priori stochastic model, postSolve
    hooks fire once (postModel never — nothing was estimated), and the
    Results carry sweep products but an empty history.

Ordering rationale (Expert B): the snooping wave runs AFTER the noise update
of the same pass so outlier tests always use the current stochastic model;
one wave per pass (not exhaustive inner snooping) keeps the fixed point of
weights, model, and mask coupled — the same argument v1's snoop makes for
one-at-a-time, lifted one level.

### 3.6 Parallelism and performance plan

- `mapreduce` (v1's strategy layer, kept verbatim in spirit, WP-9) remains
  the ONE place parallelism lives. Two parallel axes, chosen by shape:
  day-level accumulation (few systems, many days — GRACE-FO) and arc-level
  sweep (the hot loop; shared probe broadcast once, partial sums reduced in
  task order — v1's validation/v4verify.py proof carries over).
- Prefer a THREAD pool where the loaders allow it (near-zero serialisation;
  R2020b+); process pool for `spmd` multi-node. `cfg.run.parallel` selects;
  `Threshold` guards small runs.
- Prefetch: `dataSource(...,'prefetch',true)` reads day k+1 while k computes
  (v1 PrefetchLoader folded in, WP-9).
- The D8 reduced-common stage is the single biggest lever (Ellmer's own).
- Precision: QR backend default; sd via triangular solves (v1 fix kept);
  rcond guards; streaming Toeplitz whitener removes the block-diagonal error
  (D4); Hutchinson seeded.
- Speed tests with budgets are part of every engine WP (§6.1.3).

### 3.7 Old → new mapping (what carries over, what dies)

| v1 | v2 | note |
|---|---|---|
| adjust / adjustBlock / adjustVCE / adjustVCEjoint | `adjust(cfg)` | one door (F1) |
| Accumulator / Combine / BlockAccumulator | `Assembler` | scoped elim (F2,F5) |
| QRSystem / NormalEqSystem | same, internal | port (WP-7) |
| makeLoader + raw spec form | `dataSource` | one shape (§3.3) |
| buildDay | `assembleDay` | rename + ctx |
| 14-column cells | ParamSet tables | F11 |
| covModel + covOverride plumbing | `noiseModel` | first-class state (F10) |
| — | `toeplitzWhiten` | NEW, exact (D4) |
| residuals / redundancy / redundancyLags / contribution / streamContribution | `sweep` + collectors | one pass (F9) |
| vce / varComp / vceHutchinson / traceEstimate | `vce` (+ internal trace) | F8 |
| psdVCE / arcVarianceFactors / sharedProbe / dctI / covComponent / toeplitzCofactor | `ellmerUpdate` + kernels | workhorse (D1,D2) |
| applyShape / selectArOrder / selectArAcv / iterateNoiseModel | `armaUpdate` | alternative (F7) |
| snoop (flat) / outliers / robust | `snoop` (block-native) | F3 |
| diagnose / diagnoseAccum / report | `report` (block-native) | F4 |
| Results | `Results` (slim) | keep concept |
| mapreduce / PrefetchLoader | `mapreduce` + source prefetch | keep |
| util.log / verbosity / showExport | same | keep, port |
| ftest / ellipse / datum / stat / solveNonlinear / systematic | stats annex | keep |
| levinsonDurbin / armaACV / ar2arma / psd2acv / whittle / varACV / separableWhiten / armaWhiten / varWhiten / decorrFactor / checkPSD / crossCorr | signal kernels | port verified code |
| gcvLambda / applyRegularisation / CommonPrior | same, under Assembler | keep |
| poolGroups | `vce(...,'pool',...)` option | absorbed |
| satKeyOf / normalizeSatList / resolveSatSpecs / collectResiduals / whitenessCorr / getfielddef / mergeStruct / BatchIndex / synthLoader | internal or deleted | F12 |
| plot hooks | kept: `plotPostSolve`/`plotPostModel` on dataSource | ctx.schema 2 — PER-PASS firing now (v1 fired postModel only once); payload = results snapshot + ctx.proc (v1 `procOpt`, renamed at the v1→v2 boundary) + ctx.cfg; contract in WP-16(c) |

---

### 3.8 Extension track (adopted 2026-07-19)

After the redesign was written, a beyond-the-review improvement list was
put to the user. Everything on it was ADOPTED except the temporal prior
between days (the Kvas/ITSG daily-constraint direction) — that one item is
OUT OF SCOPE for v2; do not re-propose it. The scope-aware block structure
of the Assembler (WP-8) is what keeps it possible later. Where each adopted
item lives, and why it earns its place:

| adopted item | lands in | why |
|---|---|---|
| Std + correlation of the variance components | WP-23 (report hook: WP-17) | the Fisher information assembles from traces Hutchinson already delivers — near-free error bars, plus an over-resolution warning for spectral bins |
| VCE-determined regularisation (Koch–Kusche 2002) | WP-24 | the prior becomes one more observation group; lambda stops being a knob and becomes an estimate — the standard cure at higher degree / shorter windows |
| Robust IRLS reweighting | WP-25 | same problem as snooping, different failure modes; both behind cfg.outliers, sharing one weight-composition channel |
| Automatic model-order selection | WP-26 | AR order and bin count stop being silent tuning knobs |
| SINEX + GROOPS NEQ exchange | WP-27 | combination, archival, and GROOPS parity as a FILE DIFF (consumed by WP-20) |
| Checkpoint/restart | folded into WP-16 | a crash costs one pass, not a day; resume reproduces the run BITWISE (possible only because of the determinism row below) |
| Structured run journal (JSONL per pass) | folded into WP-16 | convergence pathologies become diagnosable after the fact, outside MATLAB |
| Deterministic parallel reductions | folded into WP-7 + WP-9; ALWAYS ON | fixed-order pairwise merge makes results bitwise independent of worker count; reproducibility is a property, not a cfg option (one way to do a thing) |
| Sketched diagnostics (randomised leverage) | WP-28 | Hutchinson's philosophy applied to the diagnostics side: redundancy numbers on huge parameter sets without any dense factor |
| Metamorphic test suite | WP-29 + §6.1.2 (new kind) | invariance/equivariance properties catch bug classes fixed-input tests cannot |
| CI pipeline | WP-30 — DEFERRED until a MATLAB-capable runner exists (user, 2026-07-19) | regressions caught between rounds once the infrastructure allows; the local confirming MATLAB run is the ground truth either way |
| GPU backend | WP-31, GATED on WP-20 evidence | only if the speed logs prove dense panel ops dominate; the WP-7 kernel seam is prepared either way |

Two consequences ripple through existing packages and are folded into their
prompts: (i) the cfg `snoop` block is renamed `outliers` (one policy, two
methods); (ii) WP-7/WP-9 owe BITWISE schedule-invariant reductions — which
is also what turns the WP-16 resume test and the WP-29 worker-count property
into equalities instead of tolerance statements.

---

## 4. PART III — THE v2 PUBLIC SURFACE (44 names)

Front door
  configure        - build + validate the ONE cfg struct (defaults authority)
  adjust           - run the adjustment: single pass or full iteration
Data
  dataSource       - wrap assembleDay/epochCov + tables into a source
  ParamSet         - parameter tables with scopes; nc/nl derivation
  checkSource      - user-runnable contract validator for a source
Noise
  noiseModel       - spec -> {whiten, acv, psd, memory}; types white|psd|
                     arma|var|separable|custom; optional addCov
  toeplitzWhiten   - exact streaming Toeplitz whitener (Levinson/Schur)
  dctI, psd2acv, levinsonDurbin, armaACV, ar2arma, whittle, varACV,
  checkPSD, crossCorr                     - signal kernels (ported)
Engine
  Assembler        - scoped elimination, QR|normal, regularisation, prior,
                     solve, back-substitution, local covariance
  sweep            - the one second pass; collector registry
  sharedProbe      - U^-1 * Zbar, once per iteration (Ellmer 6.5.48)
  mapreduce        - serial | parfor | spmd strategy layer
Estimation
  vce              - Helmert group components; 'trace' exact|hutchinson;
                     pooling of under-redundant groups built in
  vceUncertainty   - Fisher std + correlation of the components (WP-23)
  ellmerUpdate     - PSD amplitudes + arc variance factors (spectral|ellmer
                     components; zeroFreqRule; robust mean)
  armaUpdate       - AR/ARMA re-fit per stream + whiteness gate
  orderSelect      - AIC|BIC|CV choice of AR order / bin count (WP-26)
  snoop            - block-native screening: reject | downweight, iterated
  iterate          - the one outer loop (called by adjust; public for experts)
Results
  Results          - x, sigma0, noise state, history, sweep products, acc
  report           - the five-level verdict, block-native
  contribution     - per-group information share (Level 1 + Level 2)
Exchange
  neqRead, neqWrite - normal-equation exchange: SINEX 2.02 | GROOPS (WP-27)
Statistics annex
  ftest, ellipse, datum, stat, solveNonlinear, systematic
Support
  util.log, util.verbosity, util.showExport
Regularisation
  applyRegularisation, gcvLambda, CommonPrior

Everything else from v1 is package-private (`+lsqvce/private/`) or deleted.

---

## 5. PART IV — WORK PACKAGES AND PROMPTS

*Rows 0–29 are DONE (reconciliation, §1); §5.1–§5.29 and the §5.0 preamble
are the HISTORICAL record and keep their original wording (including paths
that no longer exist). Live packages: WP-30 below (design current) and the
R2 set in §9.*

Each package below is a self-contained implementation brief for a fresh chat.
Usage: paste §5.0 (the preamble) first, then the WP prompt, and attach the
current zip. Every WP DELIVERS THREE THINGS TOGETHER: the code, its tests
(contract + precision + speed where applicable, per §6.1), and its
documentation artifacts (help text per §6.2; the HTML page is generated in
WP-21 but the help text must be final in the WP). A WP is DONE only when the
user's MATLAB run of `run_all_lsqvce` is green and the ledger in §1 says so.

### 5.0 PROMPT PREAMBLE — HISTORICAL (rounds 0–29; superseded by §9.0)

> **Do not paste this for new work.** It predates the reconciliation and
> the 2026-08-17 single-folder restructure (`lsqvce2/`, `run_lsqvce2` no
> longer exist). Current work uses the §9.0 R2 preamble. Kept verbatim
> below as the record under which WP-0…WP-29 were executed.

> You are an expert MATLAB programmer and least-squares-adjustment
> specialist working on `lsqvce2/`, the successor of the frozen `lsqvce/`
> toolbox in the attached zip. Think step by step and motivate every answer.
> Read `HANDOVER.md` FIRST: §0 (rules), §1 (ledger), §3 (design),
> §6 (test/doc guidelines) are binding; the WP text below names the sections
> that matter for this task. The Ellmer (2018) and Kvas (2019) theses are in
> the project knowledge; equation numbers below refer to Ellmer unless said
> otherwise. Rules that bind every WP: single point of access, no aliases,
> full help text with typed/dimensioned INPUT/OUTPUT and an ERRORS section
> grouped by cause, the Claude stamp after `See also`, MATLAB coding
> standards, tests per HANDOVER §6.1, no JavaScript/Java anywhere. Where v1
> code is named "port from v1", start from that file — it is verified — and
> adapt it to the v2 contract instead of rewriting from scratch. There is no
> MATLAB in your sandbox: validate structurally, numerically in Python
> against exact formulae, and with the static audit; the user runs
> `run_lsqvce2` and returns the log as an uploaded .txt file. Ship the
> ENTIRE tree as ONE zip via present_files, update the §1 ledger row and the
> WP's own section in `HANDOVER.md`, and end with a wrap.

### WP-0 — Skeleton, setup, runner, audit
**Depends:** — .  **Read:** §0, §3.2, §6.
Create the `lsqvce2/` tree: `+lsqvce/`, `+lsqvce/private/`, `+lsqvce/+util/`,
`tests/` (FLAT, no subfolders), `examples/`, `template/`, `html/`,
`docs_src/`, `info.xml`, `helptoc.xml` stubs, `lsqvce2_setup.m` (path +
Help-browser registration; nothing writes outside the toolbox folder), and
`tests/run_lsqvce2.m` — the ONE runner: runs every `tests/test_*.m`, logs to
`tests/last_run.txt` with per-suite seconds, warning inventory by identifier,
and an explicit list of suites that LOADED (v1 lesson: a suite that fails to
load vanishes from the count). Port `docs_src/audit_tree.py` from v1 and trim
it to checks that apply to an empty tree (function/end balance, help-text
presence, documented-error-ids-exist, H1 correctness, registration
completeness); keep its verified-both-ways discipline. Acceptance: user runs
`lsqvce2_setup` + `run_lsqvce2` on the skeleton (a trivial self-test passes),
`python3 docs_src/audit_tree.py` exits 0.

### WP-1 — util: log, verbosity, options
**Depends:** 0.  **Read:** §3.7, v1 `+util/log.m`, `+util/verbosity.m`.
Port `util.log` (one logfid + one screen flag routing ALL text and warnings)
and `util.verbosity` (levels 0 silent | 1 progress | 2 stage | 3 diagnostic;
banner/heading/labelledRow formatters) unchanged in behaviour. Add
`util.options`: ONE helper implementing defaults-overlay + type/range
validation from a declarative spec (absorbs v1 `getfielddef`/`mergeStruct`;
`configure` in WP-18 is built on it). Tests: port the v1 messaging wiring
test idea — capture the log at each level and assert what appears; options
overlay precedence; unknown-field error. No speed test needed.

### WP-2 — ParamSet
**Depends:** 0,1.  **Read:** §2.2 F5/F11, §3.3.
`lsqvce.ParamSet`: value class (or validated-table functions — decide and
motivate; the panel leans VALUE CLASS wrapping a table, so ordinary table
operations remain available via `P.table`). Required variables: `name`
(string, unique within scope), `scope` (categorical:
common|system|day|arc), `system` (missing for common). Free extra columns
pass through untouched. API: constructor from table; `nc(P)`;
`nl(P, system, scope)`; `columnsOf(P, system, day, arc)` — the column-index
bookkeeping the Assembler uses to place B-columns per scope instance; merge
of several sources' tables with the CONSISTENCY RULE: all sources must agree
on the common table (the v1 resolveNc rule, generalised — disagreement is an
error naming both sources, never a silent pick). Errors grouped: "your table
is malformed" vs "your sources disagree". Tests: derivation, agreement,
column bookkeeping against a hand-built 2-system/2-day/2-arc case;
precision: none; speed: none.

### WP-3 — dataSource + checkSource
**Depends:** 0,1,2.  **Read:** §3.3, v1 `makeLoader.m` (the protocol),
v1 `BlockAccumulator.m` "THE LOADER CONTRACT".
`lsqvce.dataSource(assembleDay, Name,Value)` per §3.3 — ONE accepted form.
Provides: meta discovery (obsTypes, arcs per day without a full read where
declared; else by reading), per-arc slicing, one-day cache keyed
(day,obsType), plot hooks (`plotPostSolve`,`plotPostModel`: cell arrays of
function handles, validated at wrap time; the FIRING discipline and the
ctx payload are WP-16(c)'s contract — checkSource never calls them),
optional prefetch flag (consumed by WP-9), and the
injected ctx (§3.3). It does NOT interpret noise specs (that is WP-6's job);
it carries them. `lsqvce.checkSource(src)`: the user-facing contract
validator per §3.3 — determinism double-read, dimensions vs ParamSet,
arc-partition validity (disjoint, covering, whole epochs, epoch-major),
noise-spec resolvability (calls WP-6 when present, else defers), clear
PASS/FAIL report per check. Tests: the contract validator IS the test
fixture — build a correct synthetic source and 8+ deliberately broken
variants (wrong B width, non-deterministic, split epoch, overlapping arcs,
...) and assert checkSource catches each with the specific identifier.

### WP-4 — Signal kernels (port)
**Depends:** 0,1.  **Read:** §3.7; v1 files of the same names.
Port verbatim-in-behaviour: `dctI` (Norm 'ortho'|'ellmer'), `psd2acv`,
`levinsonDurbin`, `armaACV`, `ar2arma` (the REVERTED v1 version — v1
HANDOVER §2 records that the "robustness" variant regressed; do not
reintroduce it), plus `whittle`, `varACV`, `checkPSD`, `crossCorr`. Keep v1
error identifiers' SEMANTICS (new ids under `lsqvce:` v2 namespace). Tests:
port the v1 numerical fixtures (DCT-I round trip and Ellmer-norm parity;
Levinson vs chol on small Toeplitz; armaACV vs simulated process; ar2arma
resonant round-trip at 1.3e-14 — the regression guard); add speed budgets
for levinsonDurbin (O(p^2)) and psd2acv.

### WP-5 — toeplitzWhiten (NEW)
**Depends:** 4.  **Read:** §2.1 D4 (the whole finding), §3.1 principle 4.
`[Aw,Bw,lw, st] = lsqvce.toeplitzWhiten(acv, A, B, l, st)` — streaming exact
whitening for a symmetric-Toeplitz arc covariance given its autocovariance
vector. Algorithm: Levinson/Schur recursion delivering, per epoch k, the
order-(k-1) prediction coefficients and prediction-error variance; whitened
row k = innovation_k / sqrt(PEV_k) applied to A, B, l alike. State `st`
carries the recursion across successive row batches of the SAME arc, so
batching is exact at any BatchSize (the whole point). Convergence detection:
when the reflection coefficients fall below a tolerance, freeze the filter
(stationary regime) and switch to FIR application — cost drops from O(k) to
O(p_eff) per row. Multi-axis: epoch-major rows, temporal law per epoch,
d columns per epoch handled by applying the scalar recursion per axis (the
separable D is WP-6's job and composes outside). Guards: non-PD acv
(prediction error collapses) → error naming the nf-rule (§2.1 D3). Python
validation REQUIRED before shipping: against dense chol whitening of the
whole arc to <= 1e-12 relative, including batch-boundary crossings and the
frozen-filter regime; document measured errors in the help text like v1
does. Tests: contract (state reuse, batch invariance: results identical for
BatchSize in {7, 64, whole arc}); precision (vs dense chol); speed (budget
vs dense chol at N = 4000: must win by >5x).

### WP-6 — noiseModel
**Depends:** 3,4,5.  **Read:** §2.1 D3/D4, §3.7, v1 `covModel.m`,
`separableWhiten.m`, `armaWhiten.m`, `varWhiten.m`, `decorrFactor.m`.
`nm = lsqvce.noiseModel(spec, arcInfo)` — spec types `white | psd | arma |
var | separable | custom`, optional `spec.addCov` (per-arc additive fixed
covariance → dense whole-arc path, §2.1 D4/L5). Returns an object/struct
with `whiten(idx,A,B,l,st)`, `acv(nLag)`, `psd(f)`, `memory()` (the exact
BatchOverlap contract for filter types; 'streaming' for toeplitz via WP-5;
'wholeArc' for dense). KEY CHANGES vs v1 `covModel`: (a) the PSD/toeplitz
path uses `toeplitzWhiten` — the block-diagonal approximation is DELETED;
(b) the ACV is sized to the ARC family Nmax, not the batch (D3) — lags ==
frequencies, `numLags` option gone; (c) no covOverride merging — the loop
passes the CURRENT model explicitly (F10); the physical-frequency grid rule
(dt from the arcs) is kept. Port the v1 separable (D x R_t kron I_d), arma
FIR, var/Whittle paths behind the same interface. Tests: port v1 covModel
numerics; batch-invariance for EVERY type; separable vs brute-force dense on
a small arc (v1 has the fixture); addCov path vs dense chol; the nf-rule
assertion fires.

### WP-7 — Backends (port, internal)
**Depends:** 1.  Port `QRSystem` and `NormalEqSystem` to
`+lsqvce/private/`. Behaviour identical to v1 (QR default, rcond rank guard
on solve, merge of compressed triangular factors for day-level parallel
reduction — v1 validation/v5verify.py property). Tests: QR vs normal
agreement 1e-13 on well-conditioned; conditioning advantage on a 1e8-cond
fixture; factor-merge invariance across merge schedules, asserted BITWISE.
AMENDED (§3.8): the factor/normal merge is a FIXED-ORDER pairwise tree
keyed by block index, so results are bitwise independent of the worker
schedule — always on, not a cfg option. Keep the panel ops (syrk, qr,
triangular merge) in ONE private kernel layer so the gated WP-31 gpuArray
path can slot in without touching any caller.

### WP-8 — Assembler
**Depends:** 2,3,6,7.  **Read:** §2.2 F2/F5, §3.2, §3.3; v1
`BlockAccumulator.m` (the reference for validation messages, regularisation
hooks, prior, local recovery), `applyRegularisation.m`, `CommonPrior.m`,
`gcvLambda.m`, `datum.m`.
`lsqvce.Assembler(paramSet, opts)` — the ONE engine. `addSource(src, days)`
streams day → obsType → arc → row batches; whitening via the CURRENT
noiseModel handed in explicitly; SCOPED ELIMINATION per §2.2 F5: arc-scope
Schur-eliminated when the arc closes, day-scope at day close, system-scope
at system close; Schur blocks kept per scope instance (KeepLocal) for
back-substitution, local covariance, and the sweep. Regularisation of local
blocks before elimination and of the common block (prior) as in v1;
`solve()` returns the reduced-common solution (x, Cxx or R factor, sigma0,
rank, rcond log); `backSubstitute(scopeKey, x)`; `localCovariance`;
`reducedNormalMatrix`; per-system report table. Keep v1's grouped
error-message discipline ("your declaration is wrong" vs "your assembleDay
does not fit it"). Motivate and document the memory statement O(nc^2 +
max_scope nl^2 + one arc). Tests: block-arrow correctness vs a dense
reference solve on a small campaign (all four scopes populated) to 1e-12;
elimination-order invariance; day/system parallel-reduction equivalence
(with WP-9); regularisation + prior paths; speed budget: assembly of the
synthetic GRACE case (WP-19 fixture, pre-built here in miniature) within a
recorded budget.

### WP-9 — mapreduce + prefetch
**Depends:** 1.  Port v1 `mapreduce.m` (serial|parfor|spmd, Threshold,
Reduce/Combine contracts, purity requirement) unchanged in behaviour; add
`'pool','threads'|'processes'|'auto'` preference (§3.6) with the documented
constraint list for thread pools. Implement source prefetch (v1
PrefetchLoader folded in): with `'prefetch',true` the dataSource loads day
k+1 on a background worker (`parfeval`) while day k is consumed; graceful
serial fallback. Tests: schedule equivalence (serial == parfor == spmd on
the reduction fixtures), threshold fallback, prefetch correctness
(byte-identical days) — mark pool-dependent tests with assumptions as v1
does.
AMENDED (§3.8): the Reduce stage consumes WP-7's fixed-order pairwise tree;
schedule-equivalence tests for accumulation are BITWISE (serial == parfor
== spmd), no longer tolerance-based.

### WP-10 — sweep
**Depends:** 8,9.  **Read:** §2.2 F4/F9, §3.5; v1 `residuals.m`,
`redundancy.m`, `redundancyLags.m`, `contribution.m`,
`streamContribution.m`, `leverage.m`.
`S = lsqvce.sweep(asm, sol, collectors)` — ONE second pass over the arcs
(parallel over arcs via mapreduce; shared probe broadcast once). Collector
registry (each collector: `init`, `arc(ctxArc, whitenedData, resid, ...)`,
`reduce`, `finalize`):
  'residuals'   per-stream whitened residual series (optionally raw), vPv,
                dof, the two-route sigma0 CONSISTENCY CHECK from v1
  'redundancy'  exact per-observation h_i, r_i in the whitened metric from
                the block-native leverage (reduced common + scope Schur
                blocks — NEVER a full-parameter rebuild), per-group r_g
  'lags'        per-arc Omega_n, s_n via the sharedProbe (eqs. 6.5.45–56;
                time-domain symmetries, unrolled sums, epoch-major striding)
  'contribution' Level-1 shares r_g/sum r_g + Level-2 streamed diag(N^-1 N_g)
  'snoopstats'  the w_i, r_i, per-arc maxima the screening consumes
The gate: sum(r_i) == m - u_effective identity (v1 Level 0) computed here.
Tests: every collector vs a dense small-problem reference; identity holds to
1e-10; single-pass equivalence (running all collectors together == running
each alone); parallel == serial; speed: one sweep of the miniature GRACE
fixture within budget, and scaling recorded per arc count.

### WP-11 — sharedProbe + Hutchinson trace
**Depends:** 8.  Port v1 `sharedProbe.m` (U^-1 Zbar once per iteration, eq.
6.5.48, Z default 100, seeded) and the trace-estimation core of
`traceEstimate.m`/`streamContribution.m` as the sweep's internal primitive;
keep the v1 fix `U'\(N*(U\V))` (the operator-order bug is regression-tested).
Tests: trace estimator vs exact on random SPD (mean error ~ 1/sqrt(Z)); the
operator-order regression test; seed reproducibility.

### WP-12 — vce
**Depends:** 10,11.  **Read:** §2.1 D9, §2.2 F8; v1 `vce.m`, `varComp.m`,
`vceHutchinson.m`, `poolGroups.m`.
ONE `lsqvce.vce`: group Helmert components sigma2_g = vPv_g / r_g iterated
to a tolerance, consuming the SWEEP products (never its own data pass);
`'trace','exact'|'hutchinson'`; grouping by obsType|system|arc|handle;
`'minRedundancy'` floor with BUILT-IN pooling: below-floor groups are pooled
(mergeLow|nearest, the v1 poolGroups strategies) and the VCE re-run on the
pooled keys — closing v1 open item A in one place. The rigorous
(H_ij trace matrix) estimator from v1 varComp is kept as
`'estimator','rigorous'|'simplified'` with simplified default (motivate: the
simplified estimator is the one the Ellmer machinery composes with). Tests:
recovery of planted components (v1 fixtures); pooling changes the fixed
point and clears the floor; hutchinson vs exact within Monte-Carlo bounds;
group/stream alignment (the v1 "pin the group order" lesson) asserted.
AMENDED (§3.8): `info` returns the per-group traces and the H matrix so
WP-23 can assemble the component Fisher information — computed products
are kept, not discarded.

### WP-13 — ellmerUpdate (THE WORKHORSE)
**Depends:** 10,11.  **Read:** §2.1 D1–D6, D9 IN FULL; Ellmer §6.5; Kvas
§3.2.4/eq. 3.14; v1 `psdVCE.m`, `arcVarianceFactors.m`, `jointModelOpt.m`
(for the conventions it encoded — they become impossible by construction,
not re-implemented).
`[state, info] = lsqvce.ellmerUpdate(state, S, opts)` — per observation
group: per-frequency amplitude update alpha2_j = Omega~_j / s~_j from the
sweep's lag collectors, with `'components','spectral'` (DEFAULT — no C^n
prefactor in the lag sums; white start legitimate; derive and document per
§2.1 D2) and `'components','ellmer'` (thesis-literal prefactor + the v1
non-degenerate-C0 guard and delta-start warning, for GROOPS parity);
Shat_j = alpha2_j * S_j; zeroFreqRule 'lowest'|'keep'; optional smoothing
(v1 'Smooth', explicitly NOT a rank cure); arc factors per eqs. 6.5.27–30
with the 6.5.57 interquartile mean ('robust' default); prior/pseudo-obs
variance factor per eq. 6.5.31/6.5.58–61 reusing U^-1 Zbar; group-scale
interaction per D9 (PSD owns absolute scale; group VCE off by default under
this method). The updated state feeds WP-6 noiseModels of type 'psd'.
Python validation REQUIRED: a full miniature loop vs an exact dense
implementation of both component conventions, including the fixed-point
check that 'spectral' reshapes a white start and 'ellmer' does not.
Tests: parity fixtures generated by the Python reference (values embedded);
convergence on the synthetic GRACE fixture; the D9 no-double-count test
(enable group scale under ellmer → warning + constrained behaviour);
zero-freq and robust-mean unit tests ported from v1.
AMENDED (§3.8): `info` carries the per-bin trace sums s~_j and Omega~_j so
WP-23 can assemble the amplitude Fisher information without a second sweep.

### WP-14 — armaUpdate
**Depends:** 10.  **Read:** §2.1 D1; v1 `adjustVCE.applyShape`,
`selectArOrder.m`, `selectArAcv.m`, `residualAcv.m`, `iterateNoiseModel.m`.
The parametric alternative: per stream, empirical ACV of the de-scaled
whitened residuals (epoch-major striding for multi-axis — the v1
`lagsPooledOverAxes` lesson), ONE order selector (unify v1's two: input-type
dispatch ACV|series, motivate the criterion — ACV-fidelity with the v1
`acvTol` floor argument), AR fit via Levinson, optional ar2arma, Ljung-Box
whiteness with p-values recorded in history. Produces noiseModel specs of
type 'arma'. Tests: recovery of a planted AR(3) from residuals (v1
fixture); whiteness gate wiring; multi-axis striding regression test.

### WP-15 — snoop (block-native)
**Depends:** 10.  **Read:** §2.2 F3, §2.3 L6; v1 `snoop.m` (KEEP its
measured masking argument and the vPv/r_i identity in the help), 
`outliers.m`, `robust.m`.
`[mask, info] = lsqvce.snoop(mask, S, opts)` — one screening WAVE per outer
pass on the sweep's whitened statistics: w_i = |e~_i|/sqrt(r_i) in the
whitened metric, Bonferroni-corrected alpha, refusal floor rMinReject with
per-arc untestable counts, cap maxPerPass (default: the v1 one-at-a-time
argument lifted to "few per wave, largest first" — motivate against masking
exactly as v1's help does), actions `'reject'` (row mask consumed by the
Assembler next pass) or `'downweight'` (Huber|IGG3 variance-inflation
factors composed into the whitener next pass), `readmit` policy, and a
no-oscillation guard (a row toggling reject/readmit twice is frozen and
reported). MDB per observation in the whitened metric for the report.
Tests: planted blunders on the block fixture found across passes with zero
false positives at the corrected threshold (port the v1 measured scenario);
masking scenario (two equal blunders) resolved over waves; downweight path
converges with sigma0 → 1; refusal floor respected; oscillation guard.
AMENDED (§3.8): opts arrive as cfg.outliers (the block was RENAMED from
`snoop`; method 'snoop' selects this WP). Implement the variance-inflation
composition into the whitener as a GENERAL per-row weight channel — WP-25's
IRLS reuses that channel unchanged.

### WP-16 — iterate
**Depends:** 12,13,14,15.  **Read:** §3.5 IN FULL, §2.1 D8/D9.
Implement the loop exactly as §3.5, ~200 lines: strategies injected from
cfg; the four-part gate; history recording (per pass: alpha2 max-dev, PSD
change, sigma_g, screening changes, per-stream whiteness where 'arma',
timings per stage); reducedCommon stage + final full solve; the v1-style
ACTIONABLE non-convergence warning (name the stream/group and the knob).
Verbosity per §0 rule and the v1 banner discipline (stage labels through
`verbosity.banner` only). Tests: end-to-end convergence on the synthetic
GRACE fixture under all three noise methods; maxOuter cap warning; reduced
stage == full-resolution result within tolerance on a case where the
parametrisations nest; gate unit tests (each clause individually blocks).
AMENDED (§3.8): (a) CHECKPOINT/RESTART — after every outer pass, serialise
the complete loop state (noise state, weights, masks, history, rng state,
cfg hash) to cfg.run.checkpointDir; cfg.run.resumeFrom resumes, refusing a
cfg-hash mismatch with a documented error; a resumed run must reproduce the
uninterrupted run BITWISE (the WP-7/9 determinism makes that assertable).
(b) RUN JOURNAL — with cfg.run.journal, append one jsonencode'd line per
pass to <exportDir>/journal.jsonl: components, sigma0, gate values,
screening counts, chosen orders, per-stage seconds; `convergence()` accepts
a journal file as an alternative source; document the record schema in the
help. Tests add: resume == uninterrupted bitwise; journal parses and equals
history; cfg-hash refusal fires the documented id.
(c) PLOT HOOKS — the user's functions from each dataSource
(`plotPostSolve`, `plotPostModel`) fire EVERY outer pass: postSolve
right after the solve, postModel right after the noise update, and both
once more after the final full solve with ctx.final = true (this is a
deliberate v2 change: v1 fired postModel only once, at the end — the
user needs the results after EACH step). Fired per stream, signature
fcn(ctx); ctx is FIXED and schema-versioned (ctx.schema = 2; fields are
only ever ADDED, never renamed or removed): stage, pass, final,
stream / system / obsType, x + sd (current common estimates and their
std), common (the ParamSet), proc (THE SOURCE'S processingOptions,
round-tripped untouched — v1 ctx.procOpt), cfg (the validated run
configuration, read-only by convention), model + sigmaG (current noise
state; psd / acv / arOrder additionally on postModel), the gate values
OF THE PREVIOUS PASS (both firings precede this pass's gate; NaN on
pass 1), and exportDir / showSeconds from cfg.run (so the v1
util.showExport pattern ports unchanged). This plot context is a
DIFFERENT struct from the assembleDay ctx of §3.3 (different consumer,
different fields — the assembleDay ctx has no schema field); the full
plot-ctx schema is documented ONCE, in the `iterate` help (a PLOT
CONTEXT section), and the `dataSource` help refers to it — single point
of access. CONTAINMENT: a throwing hook produces one
warning with a documented id plus a journal entry and the campaign
CONTINUES — a plotting bug must not kill an hours-long VCE run.
Tests add: spy hooks count 2 x passes + 2 final firings in the right
order; ctx schema fields present and typed; proc round-trip identity;
a throwing hook leaves the numerical result bit-identical and raises
exactly the documented warning id.

### WP-17 — Results + report + contribution
**Depends:** 10,16.  **Read:** §2.2 F4; v1 `Results.m`, `report.m`,
`diagnose.m` (for the five-level CONTENT, not its data path),
`contribution.m`.
`Results`: slim handle carrying x, Cxx accessor, sigma0, converged, history,
noise state, sweep products, the Assembler; methods `report()`,
`residuals()`, `contribution()`, `convergence()` (per-pass view),
`summary()` (per-stream view), `disp`. `report`: the five levels, KEPT —
0 INTEGRITY (sum r_i == m - u, from the sweep), 1 GLOBAL (variance factor vs
chi-square, "read the ratio" wording kept), 2 LOCALISE (screening summary),
3 RELIABILITY (r_i distribution, MDB, untestable counts per arc), 4
PARAMETERS (correlations/conditioning of the REDUCED common block +
per-scope local summaries) — ALL fed from sweep products and the Assembler;
no data re-read, no full-parameter rebuild. `contribution` Levels 1+2 from
the sweep collector. Tests: report levels against hand-computed small
verdicts (port v1 report fixtures); Results self-description (convergence
views work from a saved Results without the builder — the v1 1.1a decision).
AMENDED (§3.8): the redundancy/leverage collector exposes an estimator seam
('exact' now; WP-28 adds 'sketch', auto above a size threshold); report
Level 1 additionally prints sigma2 ± std and the bin-correlation warning
when WP-23's info is present, degrading silently when it is absent.

### WP-18 — configure + adjust
**Depends:** 16,17.  **Read:** §3.4 IN FULL (as amended: outliers block,
orderSelect, journal / checkpoint / resume / gpu keys, solver.reg 'vce').
`configure`: the ONE defaults+validation authority (built on util.options),
returning the plain cfg struct of §3.4; every option documented with type,
dimension, default, and WHICH WP consumes it; cross-field consistency
(e.g. outliers.action requires iterate.maxOuter >= 1; noise 'ellmer' +
groupScale 'on' → the D9 warning path; solver.reg 'vce' requires a nonempty
prior; outliers 'irls' warns on set-but-unused test-only fields). `adjust(cfg)`: validate → iterate →
Results; `maxOuter = 0` is the single-GLS-pass mode (replaces v1
adjustBlock); plot hooks per the WP-16(c) contract — in maxOuter = 0
mode postSolve still fires once, postModel does not (no stochastic
estimation ran; document that). Tests: option matrix (every enum value at least constructed);
defaults snapshot test (changing a default must break a test — the v1
"hardcoded verbosity" lesson); end-to-end smoke at verbose 0/1/3 with the
captured-log assertions.

### WP-19 — Template + GRACE-FO example
**Depends:** 18.  **Read:** §3.3; v1 `template/lsqvce_template.m`,
`examples/demo_grace_dataLoader.m`, `demo_grace.m`, `demo_messaging.m`.
Rebuild the synthetic GRACE-FO case on the v2 contract: the C/D pair
(3-D positions per satellite with separable covariance and fixed per-epoch
D; scalar range-rate coupling both satellites' ARC-scoped states) + the
independent satellite E; parameter tables with all four scopes exercised;
deliberately wrong white a-priori; days with gaps producing irregular arcs.
One template (`template/lsqvce2_template.m`) with two documented stubs
(`myAssembleDay`, `myEpochCov`) that error until written, and one worked
example (`examples/demo_grace.m`) runnable end to end, plus the
every-verbosity-level demo. The worked example REGISTERS example plot
hooks (`plotPostSolve`, `plotPostModel`) demonstrating the WP-16(c) ctx
— per-pass estimate evolution and the estimated PSD — reusing
util.showExport so the headless/export path is exercised too. The example doubles as the standing test
fixture: build it ONCE in setupOnce and share (the v1 task-6 lesson,
including the determinism-vs-solve test-weakening trap — read it in v1
HANDOVER §1.4b and keep both assertions). Tests: the demo runs green and
recovers the planted gravity field within tolerance under 'ellmer'.

### WP-20 — Validation campaign
**Depends:** 18,19.  **Read:** §6.1.
Three suites: (a) ELLMER PARITY — the Python dense reference (WP-13)
extended to a full small campaign; embedded fixture numbers; both component
conventions; GROOPS parity ADDITIONALLY as an NEQ file diff via WP-27
`neqWrite` where a reference NEQ exists; (b) PRECISION — QR vs normal conditioning study, whitener
exactness across batch sizes, elimination vs dense reference, the two-route
sigma0 identity, snoop downdate identity (port the v1 5e-16/200-downdate
fixtures); (c) SPEED — budgets per §6.1.3 for assembly, sweep, whitener,
and the full loop on the WP-19 fixture, logged per suite so regressions are
visible in the runner output. Record all measured numbers in this handover
section when done.

### WP-21 — Documentation
**Depends:** 18,19.  **Read:** §6.2; v1 `docs_src/` (build_missing_pages,
build_master_pdf, pdf_graphics — port; note v1's GUARDED build_manuals.py is
NOT ported), v1 HANDOVER §1.3b–e (the font-coverage and rasterise-and-look
lessons — follow them).
Generate the HTML page set for every public name (§4), `Contents.m` grouped
as §4, `helptoc.xml`, `html/index.html`, and ONE master PDF manual
(ReportLab, drawn equations with the FreeSerif coverage assertion, workflow
diagrams from the reusable helper): Part I orientation + the engine-choice
story (now: one engine, the cfg decisions), Part II theory (BLUE/QR/
redundancy/VCE/block-arrow/Ellmer §6.5 incl. the spectral-vs-ellmer
component discussion of §2.1 D2), Part III the cfg reference, Part IV
writing assembleDay/epochCov with wrong-then-right worked examples, Part V
reading the report (the five levels + the two lies of Level 2), Part VI
installation + first steps. Audit checks: every public name has a page, is
registered, documented error ids exist in code, no JS, zero missing glyphs.

### WP-22 — Final audit + release
**Depends:** 20,21,23–29 (+30/31 only if opened).  Full audit run; the
WP-30 CI pipeline green on the release commit IF WP-30 was opened; help-text completeness sweep (port v1
`review_mfiles.py` with its house-style-aware rules); warning inventory
baselined; version + release notes (`Contents.m` version line); ledger
closed; v2 declared current in this handover; v1 remains frozen.

The nine packages below are the EXTENSION TRACK (§3.8). The numbering is
historical — slot them by the dependency column of the §1 ledger.

### WP-23 — Component uncertainty (Fisher information)
**Depends:** 12,13.  **Read:** §3.8; Ellmer §6.5.39–48 (the traces reused);
the amended `info` contracts in WP-12/WP-13.
`[sig, C, uinfo] = lsqvce.vceUncertainty(vceInfo)`: assemble the Fisher
information of all estimated components, F_kl = tr(Qinv Q_k Qinv Q_l)/2,
from the trace products WP-12/13 already computed — exact where
'trace','exact', Hutchinson otherwise with the probe count and its
4/sqrt(Z) bar propagated into `uinfo`; invert (pinv on rank deficiency,
with a named warning) → per-component std `sig` and correlation matrix `C`.
Surfaced through `Results` and report Level 1 (WP-17 hook): sigma2 ± std;
adjacent spectral bins with |corr| > 0.9 raise the OVER-RESOLVED warning
naming the bins and the numLags knob. No second sweep, no extra solves —
that is the point; say so in the help. Tests: Monte-Carlo — 200 noise
realisations of a small fixture, empirical component std within
4/sqrt(200) of predicted; the over-resolution warning fires on the
over-binned fixture and stays silent on the sane one; pinv-path warning;
exact-vs-hutchinson F within the Monte-Carlo bound.

### WP-24 — VCE-determined regularisation (Koch–Kusche)
**Depends:** 8,12,16.  **Read:** §3.8; Koch & Kusche, J Geod 76 (2002)
259–268; v1 `CommonPrior.m`, `applyRegularisation.m`; Ellmer eq. 6.5.31
(prior-group machinery reused).
`solver.reg = 'vce'`: the prior (CommonPrior/Kaula, COMMON-scope
parameters only — motivate: arc/day locals are data-determined and a prior
there masks assembly defects) becomes one more observation group with its
own component sigma2_p; each outer pass, lambda = sigma0^2 / sigma2_p from
vPv_p and the prior-group redundancy r_p via the existing trace machinery;
the update lives inside `iterate` between noise update and gate; the
(lambda, sigma2_p, r_p) triple is journaled. Guards: reg 'vce' with empty
prior → configure error; sigma2_p at its floor → the ACTIONABLE warning
("prior uninformative — tighten or drop it"). Python reference for one
full lambda iteration REQUIRED. Tests: a truncated-degree fixture that is
rank-deficient WITHOUT the prior converges and recovers truth; lambda
fixed point independent of the start (two starts a decade apart agree
within 1e-3 relative); when the data alone determine x, lambda collapses
and the solution matches the unregularised reference within 1e-10.

### WP-25 — Robust IRLS (alternative to snooping)
**Depends:** 10,15,16.  **Read:** §3.8; the §3.5 gate note; v1 `robust.m`
(the psi functions to port).
`cfg.outliers.method = 'irls'`: every observation, every pass, receives
weight w = psi(t_i)/t_i on the sweep's whitened standardised residual t_i
— no test, no alpha. psi: 'huber' (default k = 1.345; 95% efficiency under
normality — state it in the help) | 'igg3' | 'tukey'. Weights compose
through the SAME per-row channel WP-15 built; gate clause max|dw| < tol
(already in §3.5); method 'irls' deactivates the snoop wave, while MDB /
reliability reporting (Level 3) still comes from the sweep, unchanged.
Tests: clean data → all weights → 1 and the solution equals the unweighted
one within 1e-12; planted-blunder fixture — 'irls' and 'snoop'+downweight
agree on the clean-data parameters within tolerance; a 10%-contamination
fixture is run under both policies and the measured parameter bias of each
is RECORDED in this section (a comparison, not an assertion); psi unit
tests ported from v1 `robust`.

### WP-26 — orderSelect
**Depends:** 13,14,16.  **Read:** §3.8.
`[k, info] = lsqvce.orderSelect(S, spec)` — data-driven choice of (a) the
AR order p for method 'arma' and (b) the spectral bin count for 'ellmer':
criteria 'aic' | 'bic' from the whitened-innovation Gaussian
log-likelihood per stream (document the exact formula including the
parameter counts), and 'cv' = leave-k-arcs-out predictive whitened SSE
(k from spec, default ceil(nArcs/10)). Called by `iterate` at pass 1 when
cfg.noise.orderSelect ~= 'none', re-affirmed once after the reduced stage;
chosen orders journaled and shown in report Level 1. Tests: 'bic' recovers
the true order on synthetic AR(p), p in {2,5,10}, at least 18/20 seeds
(BIC consistency); the bin-count choice lands within one step of the
uncertainty-vs-resolution knee on the over-binned fixture (SHARED with
WP-23: whichever package lands first creates it under tests/); 'cv'
agrees with 'bic' on a large-sample fixture; speed: selection <= 3x one
plain noise update on the WP-19 fixture.

### WP-27 — neqExchange (SINEX + GROOPS)
**Depends:** 2,8.  **Read:** §3.8; the SINEX 2.02 format description
(IERS); the GROOPS normal-equation file layout — VERIFY against the GROOPS
source (github.com/groops-devs/groops, the normals file I/O) IN THE ROUND
and cite file + commit in the help; v1 ParamSet naming.
`neqWrite(file, neq, ps, meta, 'format','sinex'|'groops')` and
`[neq, ps, meta] = neqRead(file)` with autodetect — the argument is
named `neq`, NOT `sys`, because "system" already means a satellite
system in this toolbox (ParamSet scope). The Assembler gains
`getNormals(scope)` returning that struct: the REDUCED common-scope
(N, n, lPl, m_eff)
so exchange never touches local blocks; SH parameter naming mapped to
SINEX CN/SN degree-order codes (document the mapping table); SINEX text at
full printable precision with an honest precision-loss note; the GROOPS
path byte-exact to the verified layout. Motivate in the help: combination
with other centres, archival of expensive accumulations, GROOPS parity as
a file diff (WP-20 consumes it). Tests: round-trip — groops bitwise, sinex
<= 1e-13 relative on neq.(N, n, lPl); a small third-party SINEX fixture
parses
(embed it); getNormals equals the dense reference on the WP-19 fixture;
malformed files raise the documented ids.

### WP-28 — Sketched diagnostics
**Depends:** 11,17.  **Read:** §3.8; §2.2 F4 (what this retires at scale);
Ellmer §6.5.39–48.
Randomised leverage: h_i estimated from the SAME shared probe block (reuse
`sharedProbe`; no additional solves), r_i = 1 - h_i, driving report
Level 3 and `contribution` when the common block exceeds a size threshold
(auto; overridable 'leverage','exact'|'sketch' — the WP-17 seam). Level 0
wording under sketch: the integrity identity is reported WITH its
4/sqrt(Z) Monte-Carlo bar instead of exact equality — say so in the level
text, never silently. Tests: sketch vs exact r_i on a mid-size fixture
within the bound at Z = 100; the sum of sketched r_i within the bound of
m - u; the auto threshold honoured (assert the estimator choice); speed:
sketch <= 0.1x exact on a 5000-parameter fixture.

### WP-29 — Metamorphic suite
**Depends:** 18,19.  **Read:** §6.1.2 (the METAMORPHIC kind, added
2026-07-19); §3.8.
One file `tests/test_metamorphic.m`; every property runs the FULL
`adjust(cfg)` on a small WP-19 fixture variant: (a) arc-permutation
invariance — shuffle arc order within days; x agrees to 1e-12 (the
permutation re-seats values in the fixed reduction tree, so eps-level, not
bitwise — say so in the test comment); (b) unit-rescaling equivariance —
scale one obsType's l and columns by 1e3: x rescales exactly, sigma0 and
the component RATIOS invariant; (c) split-merge — a day split at an arc
boundary equals the whole (fixture variant without day-scope parameters;
motivate in-file); (d) duplication — every row duplicated at half weight
reproduces x and Cxx; (e) worker-count invariance — serial == parfor
BITWISE (the WP-7/9 promise asserted at pipeline level,
pool-assumption-filtered). Each property cites in a comment the bug class
it guards.

### WP-30 — Development workflow: public GitHub + hosted CI + cloud-session loop (D8, adopted and corrected 2026-08-17)
**Depends:** —.  **Status: OPEN. DECISION RECORD D8 — CORRECTED (user,
2026-08-17, night; the first D8 recording of the same night was wrong on
three points and is replaced by this text): GITHUB, PUBLIC REPOSITORY,
HOSTED CI, CLAUDE-SIDE PROMPTING, transferable to other projects.**

**The three corrections:** (1) the repository is PUBLIC — that is what
makes MATLAB on GitHub-hosted runners FREE (Route A; the batch-token
route B and the runner bridge are OBSOLETE for this project). The earlier
privacy requirement is withdrawn by the user; consequences are binding:
measurement data and secrets NEVER enter the repo (already pinned), and
the LICENSE is **MIT** (USER DECISION, pinned 2026-08-18; `LICENSE` file
at the root, (c) 2026 Matthias Weigelt). (2) The work-package
prompt is given ON THE CLAUDE SIDE — in chat/Projects or via Cowork — NOT
via GitHub triggers: no @claude action, no `claude.yml`, no issue
template. The execution vehicle (per the author's PROVEN pipeline,
`references/BEST_PRACTICE_v4-1.md` Appendix D, adopted 2026-08-17): the
CHAT/COWORK CONTAINER ITSELF — an ephemeral Linux box holding a full
clone in /home/claude. It pre-validates numerics in Python, edits, runs
the static gates on a fresh mirror, pushes branch `claude/wp-NN` over
HTTPS, and drives ALL coordination through the GitHub REST API: create
the PR with measured numbers in its body, CHECK PR STATE BEFORE EVERY
PUSH, poll check-runs at suite cadence, cancel/re-run hung workflows.
GitHub is the only durable state — anything unpushed dies with the
session. Claude Code cloud sessions remain the token-free ALTERNATIVE
vehicle (OAuth grant at claude.ai/code), useful from the phone. (3) LOCAL MATLAB MAY BE USED FOR TESTING
when the machine is available, via the OFFICIAL MathWorks MATLAB MCP
SERVER (github.com/matlab/matlab-mcp-server, free/open-source, verified
2026-08-17: installs into Claude Desktop as an .mcpb extension bundle
plus the Filesystem extension, or into Claude Code via `claude mcp add`;
drives the locally licensed MATLAB; tools: start/quit MATLAB, run code,
check code). Amended the same night: the bridge test is the MANDATORY POSTFLIGHT
stage of the loop (below); additional pre-flight use before pushing
stays optional.

**The corrected D8 loop (postflight amendment, same night):** THINK —
discussion and design in chat/Projects (these sessions), decisions pinned
HERE. PROMPT — the WP prompt is handed to Claude in chat or Cowork; the
SESSION'S CONTAINER (Appendix-D pipeline) opens a BRANCH, develops,
PUSHES, and opens a DRAFT pull request immediately after the first push
*(amended 2026-08-19, user decision: `ci.yml` now triggers on
`pull_request` + push-to-`main` only — the previous all-branches `push`
trigger ran every check TWICE once a PR existed, doubling MathWorks
downloads and racing parallel cache saves, the duplicate-key incident's
root exposure; the draft PR restores one merge-preview verdict per dev
push, and `main` keeps its own post-merge run)*.
CI — `ci.yml` triggers on EVERY push to EVERY branch (checks `audit` +
`matlab`; setup-matlab@v3, release R2026a, free on the public repo;
centrepiece `tests/ci_gate.m`); the session POLLS the check-runs API at
suite cadence and fixes reds itself — the user steps in only where
judgment is needed — until green. POSTFLIGHT —
once CI is green: at the PC, Claude checks out the branch and runs the
SAME gate in the LOCAL MATLAB via the MCP bridge. This is a second
platform (Windows + the user's R2026a vs the runner's Linux), not a
repetition; its result is quoted in the PR description. PR — only when
BOTH are green is the pull request opened. REVIEW — the user reads diffs
+ the log artifact, requests changes, MERGES (merging is the user's
alone). RECORD — the same PR updates this HANDOVER's ledger; DONE =
merged to `main`, checks green, date + PR number + short hash.
Enforcement = branch protection requiring `audit` and `matlab`. If the
PC is off, CI-green branches WAIT at postflight — the PR is the promise
that both worlds agreed. Long postflight runs go as detached
`matlab -batch` jobs writing a result `.mat`, polled with `isfile`
(Appendix D.2 step 6), respecting the MCP call window.

**Session preflight (Appendix D.4, binding for every working session):**
`git fetch` FIRST (a stale base is the §6.10 class of release error);
probe the token file and its scopes with one API call; check whether the
MCP bridge is present — an absent bridge means static gates + CI only,
which is a RECORDED VERIFICATION MODE AND A DEBT ENTRY, not an
equivalent; know the MCP call window (~4 min) and the suite duration —
every polling and detaching decision derives from them.

**Auth note (CORRECTED 2026-08-17 after the user produced
`BEST_PRACTICE_v4-1.md`; supersedes the "no tokens anywhere" claim of
the same night, which held only for the cloud-session vehicle):** the
chat/Cowork container has no GitHub identity, so the author provisions
a FINE-GRAINED PERSONAL ACCESS TOKEN, stored as a PROJECT FILE and read
at call time. It serves two needs (Appendix D.3): (1) `git push` over
HTTPS from the container; (2) everything git cannot say — the REST
coordination layer: create the PR, read PR state before every push,
re-target stacked PR bases, poll check-runs, cancel/re-run hung
workflows (the setup-matlab hangs were resolved exactly this way, which
is why ACTIONS RW is provisioned up front — a missing scope is
discovered mid-incident, i.e. too late). SCOPES: Contents RW, Pull
requests RW, Actions RW; least privilege beyond that. HANDLING RULES
(binding, restated in CLAUDE.md): read from its file at call time,
passed in a header; never echoed into command output, never committed,
never baked into a remote URL a log could print. The user's OWN pushes
use Git Credential Manager's browser sign-in; cloud sessions
(alternative vehicle) use the OAuth grant and need no token; public
clones need none; the MCP bridge is purely local.

**DELIVERABLES of this WP (fresh chat, §9.0 preamble):**
`.github/workflows/ci.yml` (hosted variant ONLY: triggers `push` on ALL
branches + `pull_request`; ubuntu-latest,
setup-matlab@v3 `release: R2026a`, LSQVCE_CI honoured, artifact
`tests/last_run_all.txt`); `tests/ci_gate.m` (0 failed, 0 incomplete,
suites run == suite files on disk, log written); `CLAUDE.md` (pointer to
THIS document §0/§1/§6/§8/§9; commands incl. the MCP-bridge note; hard
rules: gate never bare runtests, no aliases, one WP = one branch = one
PR, PR opened only after CI green AND local postflight green (postflight
result cited in the PR description), ledger update in the same PR, Claude stamp in touched help blocks);
`.gitignore`/`.gitattributes`; `references/` (both theses; `BEST_PRACTICE_v4-1.md` placed there
2026-08-17 — Appendices C/D are the operating runbook); the LICENSE
decision surfaced to the user; CLAUDE.md additionally restates the
token-handling rules from the Auth note. In-round validation: YAML parses; first
green run on `main`; the red test (a deliberately failing branch must
show red checks and an unmergeable PR); Cowork/cloud-session hand-off
verified once end to end; record dates HERE.

**TRANSFERABILITY (part of D8):** portable pattern = {one gate entry
point per repo + `ci.yml` + `CLAUDE.md` pointing at that repo's SSOT +
the cloud-session habit + optional MCP bridge for local toolchains}.
New-repo checklist: copy `ci.yml`, swap the gate command and toolchain
step, write the gate and `CLAUDE.md`, grant the repo to Claude's GitHub
access, protect `main`. Lift `ci.yml` into a reusable `workflow_call`
only when a second MATLAB repo exists.

**ROUND NOTE (2026-08-18, WP-30 round 1; container session per Appendix D,
model Claude Fable 5).** Shipped on branch `claude/wp-30`:
`.github/workflows/ci.yml` (D2 reference, corrected to the deliverable
spec: `push` on ALL branches + `pull_request`; audit job = audit_tree +
v10review on python 3.12/numpy, verified exit 0 in-session; matlab job =
setup-matlab@v3 R2026a, `LSQVCE_CI=1`, `lsqvce_setup; ci_gate`, log
artifact `if: always()`); `tests/ci_gate.m` — the ONE CI entry, judgment
factored into the PURE `tests/private/gateVerdict.m` so the verdict is
unit-testable WITHOUT recursively re-running the suite; new suite
`tests/test_ci_gate.m` (14 test points: green case, every red condition —
load error, fail, incomplete, forgotten-on-disk, phantom, missing/empty
log — reason accumulation, input validation; §6.1 fail-on-old-code: the
gate does not exist on the old tree); `run_all_lsqvce` gained the `skip`
MECHANISM (prints + logs every skip on the headline; POLICY stays in
ci_gate, rule 3), and `test_ci_gate` was added to its suite list so the
gate's own coverage identity holds (19 test_*.m on disk = 19 in list,
test_speed excluded as the documented benchmark); `CLAUDE.md` (pointer
per D3 tier 2, token-handling rules restated); `.gitignore` /
`.gitattributes` per D4 (+ `references/local/`). DEVIATION, pinned: the
thesis FULL TEXTS are NOT committed — the repo is PUBLIC and neither
document states a redistribution licence; `references/README.md` carries
the full citations + open-access pointers, and the gitignored
`references/local/` is the drop-in slot for repo-based sessions (the
migration item's need — sessions without chat-project knowledge — is met
without republication). The user may overrule with evidence of a CC
licence. LICENSE for the repo itself: **MIT** (user decision
2026-08-18, `LICENSE` committed in this round). SETUP_GUIDE 4.1 aligned with
what actually ships. Postflight and PR: WAITING at the user's PC per the
loop. CI evidence of this branch: next dated line below.

**CI EVIDENCE (2026-08-18, recorded from the check-runs API + the log
artifact).** Branch `claude/wp-30`, commit 5a1b7c8: `audit` GREEN (9 s;
audit_tree 15 checks clean + v10review T1–T7 pass on the runner) and
`matlab` GREEN (6 m 17 s incl. MATLAB R2026a install; run
32143388026). Log artifact read per §0: **482/482 passed, 0 failed, 0
incomplete**, headline records `1 SUITE(S) SKIPPED ON REQUEST:
test_parallel` (the LSQVCE_CI policy, visible as designed); coverage
identity 19 run + 1 skipped + 1 excluded = 20 test_*.m on disk (the new
`test_ci_gate` ran, 14/14, 0.6 s); warning inventory = only the known
intentional warn-path warnings (datum auto-detect, rejection cap);
wall clock 270.5 s, dominated by test_adjustVCE (232.1 s, expected: VCE
outer loops). The D2 headless question is hereby VERIFIED: all four
figure-opening suites passed on the display-less runner. **RED TEST
(same date):** branch `claude/wp-30-redtest` = a green-running suite
file deliberately absent from the runner list; verdict RED at the
FIRST gate — audit FINDING `[tests] test_redtest_dummy is on disk but
NOT in run_all_lsqvce` (static disk-vs-list check, fail-fast, saves the
MATLAB minutes), `matlab` skipped. Layering, not duplication: the
audit catches this class statically; `gateVerdict` additionally
catches at runtime what statics cannot (load errors, skip accounting,
phantom entries, fails/incompletes, the log artifact) — pinned by
test_ci_gate. Bare `runtests` would have been green on that branch.
Red-test branch deleted after recording. An UNMERGEABLE red PR
additionally requires branch protection = USER step (Guide 4.4; the
session token deliberately has no Administration scope).

**POSTFLIGHT (2026-08-18, run BY CLAUDE via the MATLAB MCP bridge —
user decision of the same day: the bridge session IS the postflight
vehicle when present).** Platform: the user's Windows 11, MATLAB
R2026a Update 4, PCT 16 workers; branch checked out at each step,
gate = detached `matlab -batch postflight_run` (wrapper around
`lsqvce_setup; ci_gate`) polled with `isfile` per Appendix D.2 step 6.
THREE runs, TWO REAL FINDINGS — exactly what the second platform is
for, since CI SKIPS test_parallel: (1) `hasBackgroundPool`'s
exist()-based probe (class+file) is FALSE-NEGATIVE under `matlab
-batch` (measured: license=1, ver=1, class=0, file=0 while
backgroundPool starts with 16 workers) — the prefetch test silently
filtered on every batch run (7bf5c83). (2) The SAME pattern sat in
PRODUCTION code: `PrefetchLoader.poolOK` — prefetch silently degraded
to passthrough in every batch campaign. Fix: ONE probe
`lsqvce.util.hasBackgroundPool` (probe-by-USE, full help, HTML page,
Contents/helptoc/index rows); `poolOK` and the local test guard
DELETED (rules 3+8), all call sites redirected; consistency test
`testBackgroundPoolGuardAgreesWithUse` pins probe == use and FAILS on
the old code under -batch (§6.1) (ce94c2a). INFRA lesson, pinned for
every future detached run: launch as `start "" /min cmd /c "matlab
-batch <script> > log 2>&1"` — WITHOUT the full stdout/stderr
redirection the spawned MATLAB inherits the bridge call's dying
console handle and dies mid-run with iolib:badbit (`-logfile` does not
protect against this; observed on run 1). FINAL VERDICT: ci_gate
GREEN, **513/513 passed, 0 failed, 0 incomplete**, test_parallel
31/31 INCLUDING the previously-filtered prefetch test, 185.4 s wall
clock (2026-08-18 17:25 local). Both worlds now agree: CI green
(Linux, 482/482 + visible skip) AND postflight green (Windows,
513/513, nothing skipped). Also gitignored: `custom_toolbox.json`,
`helpsearch-v*/` (builddocsearchdb artifacts of lsqvce_setup, found
dirtying the postflight tree). Remaining before merge: branch
protection (user, Guide 4.4).

*The sections D1–D5 below predate the D6 decision and are kept as the
evaluation record; D6 is the operative design.*

**D1 — Where MATLAB runs (the question that deferred this package).**
- **Route A — public repo, hosted runners.** `matlab-actions/setup-matlab@v3`
  licenses MATLAB (and requested `products:`) free of charge for PUBLIC
  projects on GitHub-hosted runners. Zero infrastructure; also lets Claude
  read the tree in chat rounds via raw URLs. Requires publishing the code.
- **Route B — private repo, hosted runners.** Needs a MATLAB batch licensing
  token (MathWorks Batch Licensing Pilot form) stored as a secret and mapped
  to `MLM_LICENSE_TOKEN`. Availability depends on the licence type — ask the
  licence administrator first.
- **Route C — private repo, self-hosted runner on the user's MATLAB
  machine.** Uses the existing local licence; the runner automates today's
  "local confirming run". HARD SAFETY RULE: never pair a self-hosted runner
  with a public repo (fork PRs would execute arbitrary code on the machine);
  keep the repo private, run the runner under a dedicated non-admin account.
- Recommendation: A if the code may be public, else C (no new licence).

**Route C prerequisites (pinned 2026-08-16, user inquiry):** (1) PRIVATE
repo, current tree pushed, tag `r2-baseline`. (2) Dedicated NON-ADMIN user
on the MATLAB machine — every CI job runs as this user. (3) Runner from
repo Settings → Actions → Runners → New self-hosted runner (four copy-paste
commands), installed as a SERVICE (Windows: offered by the config script;
Linux/macOS: `./svc.sh install && ./svc.sh start`), label `matlab`;
outbound-only HTTPS, no inbound ports. (4) `matlab` on the PATH the service
sees (Windows system PATH; Linux symlink or the runner's `.env`).
(5) LICENCE GOTCHA — the gate of this route: as the runner's user,
`matlab -batch "disp(version)"` must start with NO interactive prompt;
online licensing needs one desktop sign-in under that account first, a
Designated-Computer licence file usually needs nothing. (6) ci.yml with
`runs-on: [self-hosted, matlab]` and the setup-matlab step removed; the
`audit` job stays on ubuntu-latest (fits the free private-repo hosted
minutes). (7) Hardening: allow-list actions (`actions/*`,
`matlab-actions/*`, `anthropics/*`); review any PR touching
`.github/workflows/` before merge — repo write access equals code execution
on the machine. Machine off => runs queue until the runner reconnects; the
runner self-updates. RESOLVED 2026-08-17: OS = Windows, release = R2026a
(after editing the SYSTEM Path, restart the runner service so it sees it);
tier 3 is decided after the first green run.


**D2 — The pipeline** (`.github/workflows/ci.yml`; pin actions, PIN the
MATLAB `release:` to the user's local version — the local run has been the
ground truth, so CI must match it, not "latest"):

```yaml
name: lsqvce CI
on:
  push: { branches: [main] }
  pull_request:
concurrency:
  group: ${{ github.workflow }}-${{ github.ref }}
  cancel-in-progress: true
jobs:
  audit:                       # < 1 min, no MATLAB: fail fast
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v6
      - uses: actions/setup-python@v5
        with: { python-version: '3.12' }
      - run: pip install numpy
      - run: python docs_src/audit_tree.py
      - run: for f in validation/v10review.py; do python "$f"; done
  matlab:
    needs: audit
    runs-on: ubuntu-latest     # Route C: runs-on: [self-hosted, matlab]
    env:
      MLM_LICENSE_TOKEN: ${{ secrets.MLM_LICENSE_TOKEN }}   # Route B only
      LSQVCE_CI: 1
    steps:
      - uses: actions/checkout@v6
      - uses: matlab-actions/setup-matlab@v3   # Routes A/B; on C omit if the
        with:                                  # local MATLAB is on PATH
          release: R2026a                      # == the user's local release
          # products: Parallel_Computing_Toolbox   # only if pool suites run in CI
      - uses: matlab-actions/run-command@v3
        with:
          command: lsqvce_setup; ci_gate
      - uses: actions/upload-artifact@v4
        if: always()
        with: { name: run_all_log, path: tests/last_run_all.txt }
```

**`tests/ci_gate.m` (new, the ONE CI entry; never bare `runtests`):** calls
`run_all_lsqvce`, then errors unless ALL of: 0 failed, 0 incomplete, suites
RUN == suite files ON DISK (lesson §0a-2: a count of failures is not a
statement about coverage — a suite that fails to LOAD must trip the gate),
and the log was written. Honours `LSQVCE_CI=1`: cap any parallel pool at 2
workers and, until PCT is licensed on the runner, skip the pool suites via
run_all's existing filter, PRINTING what was skipped into the log. Hosted
Linux runners have no display: verify once that the plot-hook/showExport
tests are headless-safe under `-batch` (off-screen figures), else guard them
with the same env flag. Branch protection on `main` requires `audit` +
`matlab` green.

**D3 — The Claude development workflow** (three tiers, they coexist):
1. **Chat rounds (unchanged).** The zip loop keeps working; the repo is the
   archive and the CI the regression guard between rounds.
2. **Claude Code on the MATLAB machine — the real unlock.** With Claude Code
   running where MATLAB lives, the ground-truth gate becomes SELF-SERVE:
   Claude edits, runs `matlab -batch "lsqvce_setup; ci_gate"`
   itself, reads `tests/last_run_all.txt` itself, iterates, commits, pushes;
   CI re-confirms on the PR. The round-trip (zip → user runs → log upload)
   collapses to one session. Repo root gets a `CLAUDE.md` that is a POINTER,
   not a second status document (§0 rule stands): it names
   `HANDOVER.md` §0/§1/§6/§8/§9 as binding, lists the three commands
   (setup, ci_gate, audit), and restates the hard rules (run_all_lsqvce
   never runtests; no aliases; one WP per branch; update ONLY
   HANDOVER.md; Claude stamp in every touched help block).
3. **@claude on GitHub** (`anthropics/claude-code-action@v1`): quick setup
   via `/install-github-app` from Claude Code, or manual per
   code.claude.com/docs/en/github-actions. Auth: `ANTHROPIC_API_KEY` secret
   (per-token billing) or `CLAUDE_CODE_OAUTH_TOKEN` (Pro/Max subscription,
   generated locally with `claude setup-token`). Minimal workflow
   `claude.yml`: triggers `issue_comment`/`pull_request_review_comment`
   (types: created), job permissions `contents: write`,
   `pull-requests: write`, `issues: write`, steps checkout@v6 +
   `anthropics/claude-code-action@v1` with the chosen secret. Usage: one
   GitHub ISSUE per R2 package whose body is the §9.0 preamble + the WP
   prompt; `@claude` implements on a branch and opens the PR; the D2 gates
   judge it; the user reviews and merges. On routes A/B the Claude job may
   itself add setup-matlab so Claude-in-CI runs the suite before pushing.

**D4 — Process mapping (the loop, amended when this WP lands):** one WP =
one branch (`wp-33-joint-locals`) = one PR whose description is the §9
prompt; the §1 ledger row flips DONE only when the PR's `matlab` check is
green ON MAIN — amend §0's loop text and the §9.0 preamble accordingly (the
"user uploads the log" sentence becomes "CI green is the gate; the local run
remains for suites CI skips"). Import: `git init` from the current zip,
first commit tagged `r2-baseline`; `.gitignore`: `*.asv`, `*.autosave`,
`tests/last_run_all.txt`, `exports/`, `.DS_Store`, and a `data/`
exclusion — MEASUREMENT DATA NEVER ENTERS GIT (pinned 2026-08-17): the
gate runs on the synthetic fixtures only, real-data campaigns reference
local paths in their configuration; `.gitattributes`:
`* text=auto`. `html/` IS versioned (the doc set ships). Secrets never in
the tree.

**D5 — GitLab portability (pinned 2026-08-17, user inquiry; CLOSED the
same day: the university GitLab cannot be used — GitHub Route C is
CONFIRMED. Kept as reference for future multi-repo generalisation).** Tiers 1–2 and
the gate are HOST-AGNOSTIC: the git remote, `tests/ci_gate.m` and `CLAUDE.md`
are unchanged. Translation table: `ci.yml` → `.gitlab-ci.yml`; self-hosted
runner + labels → GitLab Runner installed as a Windows SERVICE with the
SHELL executor, runner tag `matlab`, job `tags: [matlab]`; checkout is
implicit; `upload-artifact` → `artifacts: {when: always, paths: [...]}`;
branch protection + required checks → protected branch + the "pipelines must
succeed" merge check; secrets → masked CI/CD variables; PR → MR. Licensing
identical to Route C (there is NO GitLab analogue of GitHub's free
public-repo MATLAB licensing — irrelevant on this route). University
instance: ask the admins (a) whether SHARED runners exist (then the `audit`
job runs there in a `python:3.12` image instead of the user's machine),
(b) whether a campus network licence could later back an institutional
MATLAB runner, (c) whether the instance is VPN-only — the runner service
then needs that connectivity PERMANENTLY. Tier 3 on GitLab EXISTS, BETA,
maintained by GitLab (code.claude.com/docs/en/gitlab-ci-cd, verified
2026-08-17): a `claude` job on merge-request/manual/web triggers with a
masked `ANTHROPIC_API_KEY` (or Bedrock/Vertex OIDC); true @claude COMMENT
mentions need webhook wiring — more DIY than GitHub's app. Generalisation
lever for MANY repos: register the runner at GROUP level (every project in
the group shares the MATLAB machine) and keep ONE `ci-templates` project
that each repo `include:`s with variables (gate command, tags); the GitHub
equivalent is a reusable workflow (`workflow_call`).

**D6 — the "no CI" variant (pinned 2026-08-17; adopted that evening,
SUPERSEDED the same night by D8 — kept as record; `ci_gate.m` carries
over, the hook/verify pieces become optional local extras).**
Question: with Claude Code running the gate locally, is CI needed at all?
Honest answer: the EXECUTION is identical (Route C's runner is this same
machine); CI uniquely adds (a) fresh-clone verification of the exact
commit — catches "passes locally, file never committed", (b) platform
ENFORCEMENT — red cannot merge, (c) server-recorded verdicts and the
tier-3 option (the GitLab @claude integration runs AS a pipeline job, so
no CI = no tier 3 on either host). (a) and (b) can be rebuilt locally:
a clean-clone gate step (clone HEAD into %TEMP%, run `ci_gate` there,
write a verdict file keyed to the commit hash) plus a git PRE-PUSH HOOK
refusing to push without a fresh green verdict for HEAD — enforcement in
the machine, not in prose (§0 lesson 1); CLAUDE.md rule: `--no-verify` is
forbidden. (c) is genuinely lost; acceptable for solo work. If chosen:
remote = university GitLab as PLAIN hosting (SSH/HTTPS; the D5 VPN note
applies; protected `main` + MR-based self-review still work without
pipelines), WP-30 rescopes to "local gate + hooks + repo workflow"
(`ci_gate.m` unchanged and central; new clean-clone verify script +
pre-push hook), the operator guide carries no runner parts — done: `SETUP_GUIDE.md`,
2026-08-17 (~30-40 min total setup), and ledger DONE entries record
date + commit hash + "gate green (clean clone)".

**D7 — fully-remote variant (pinned 2026-08-17, user inquiry; ABSORBED
into the adopted D8 the same night — see the decision record above).** Architecture: Claude
Code CLOUD SESSIONS (Anthropic-managed VMs; started from browser, Desktop
app or phone; clone/push and PR creation are GITHUB-ONLY — GitLab repos
can only be bundled in, no push-back; verified 2026-08-17 against
code.claude.com/docs/en/claude-code-on-the-web) and/or @claude CI jobs;
MATLAB exists ONLY as a hosted-runner CI job judging every PR. GATE ITEM:
a private repo with no self-hosted runner needs a MathWorks BATCH
LICENSING TOKEN (Batch Licensing Pilot / licence administrator) — without
it the only variant is a public repo, which violates the privacy
requirement. If granted: the @claude job runs setup-matlab itself and
executes the gate IN-SESSION (agent not blind); cloud web sessions have
no MATLAB and rely on the PR verdict; hosted minutes are consumed
(2000/month free on private repos). Adoption path if ever chosen: keep
everything from D6 (`ci_gate.m`, hooks, `CLAUDE.md`, repo discipline),
push the repo to GitHub, add `ci.yml` (D2 reference) + `claude.yml` —
an ADDITIVE upgrade, no rebuild. University GitLab cannot host this
variant (no admin-provisioned MATLAB runners available).

**Migration item (pinned 2026-08-17):** Claude Code has NO chat-project
knowledge, so the Ellmer (2018) and Kvas (2019) thesis text files move INTO
the repo under `references/` as a WP-30 deliverable; `CLAUDE.md` points at
them, and the §9.0 preamble sentence "theses are in the project knowledge"
is amended to "theses are in `references/`" when this WP lands. After that,
the implementation loop is Claude Code end to end (repo replaces the zip,
`ci_gate` replaces the log upload); the chat/Desktop apps remain optional
for planning and review only.

**Operator guide:** `SETUP_GUIDE.md` at the toolbox root — REWRITTEN
2026-08-17 (night) for the adopted D8 loop (Parts 0–6 + the Appendix-A
runner bridge; earlier same-day versions covered Route C and D6). HOW-TO
only; design and status stay HERE; the WP-30 round keeps it consistent
with what it ships (check names `audit`/`matlab`, gate command, issue
template, both `runs-on` variants).

**In-round validation for this WP:** the YAML parses (python-yaml), action
versions pinned, one deliberately red commit proving each gate trips (then
reverted); record the first live run's date and link HERE. Coverage floor:
start by RECORDING the real number the first run reports in this section,
then set the floor one point below it and raise it per round.

### WP-31 — GPU backend (GATED)
**Depends:** 7,20 — DO NOT OPEN unless the WP-20 speed logs show dense
panel ops (syrk/qr in accumulation or solve) dominating wall time on the
user's real campaigns; record that evidence here first.
**Read:** §3.8; the WP-7 kernel seam.
`cfg.run.gpu = true`: the WP-7 kernel layer dispatches panel ops to
gpuArray (double precision only — the conditioning story of §2 forbids
single; say so in the help), gathering at the fixed-order merge so the
reduction stays deterministic. The Toeplitz whitener STAYS on CPU (the
recursion is sequential — state why in the help). Graceful: no device →
warning + CPU path; device tests carry assumption filters and count as
incomplete (§6.1.3). Tests: CPU vs GPU <= 1e-10 on the WP-19 fixture;
bitwise equality across gpu on/off is NOT promised (in-panel op order
differs) — document that honestly, the merge order is; speed LOGGED, not
asserted.

---

## 6. PART V — GUIDELINES (referenced by every prompt; stated once)

### 6.1 Test guidelines

**6.1.1 Organisation.** `tests/` is FLAT. One file per public unit:
`test_<name>.m`, function-based tests with `setupOnce` fixtures; expensive
fixtures (the WP-19 campaign) built once and SHARED, with the determinism
assertion kept separate from the solve assertion (v1 §1.4b lesson). The ONE
runner is `tests/run_lsqvce2.m`; never `runtests('tests')`. The runner log
records: pass/fail/incomplete, WHICH FILES LOADED, per-suite seconds, and
the warning inventory by identifier. A green gate means: 0 failed AND every
suite loaded AND no NEW warning identifier AND no speed budget exceeded.

**6.1.2 The test kinds — every WP states which it ships.**
- CONTRACT tests: the promises in the help text. Input validation fires the
  DOCUMENTED identifier (and only that one); dimension/type of every output;
  determinism where promised; state-carrying calls are order/batch
  invariant; error messages name the offending field. Documented-id ↔
  code-id agreement is ALSO enforced statically (audit check, v1 check 14).
- PRECISION tests: against an independent dense/analytic reference, with an
  explicit tolerance and a comment saying where the tolerance comes from
  (machine-eps scaling argument or a measured Python number). House
  tolerances: exact-identity checks <= 1e-10 relative; factorisation
  round-trips <= 1e-12; Monte-Carlo quantities within 4/sqrt(Z). Never
  loosen a tolerance to make a test pass — that is a finding.
- SPEED tests: `timeit`/`tic-toc` medians with a RELATIVE budget (vs a
  recorded baseline op on the same machine, e.g. "sweep <= 25x one dense
  chol of the largest arc"), so budgets survive hardware changes; absolute
  seconds are LOGGED, not asserted. Speed suites are marked so the runner
  reports them separately.
- METAMORPHIC tests (added 2026-07-19; WP-29 hosts the pipeline suite):
  invariance/equivariance properties — permutation, rescaling, split-merge,
  duplication, worker count — asserted through the public API. Each
  property states whether the expectation is BITWISE (fixed-order
  reductions) or eps-level, and why.

**6.1.3 Regression discipline.** Every defect fixed during v2 gets a named
regression test citing the WP section that records it (the v1 habit that
worked). Randomness is seeded per test. Tests that need a pool carry an
assumption filter and are counted as incomplete, not skipped silently.

### 6.2 Documentation guidelines

**6.2.1 Help text template** (binding for every function/class/method):
H1 line `%LSQVCE.NAME  One-sentence purpose.`; call forms; a MOTIVATED
description (why, and any measured numbers backing claims); `INPUT` and
`OUTPUT` blocks with NAME, TYPE, DIMENSION, DESCRIPTION for every argument
and field; `OPTIONS` with defaults; `ERRORS` grouped BY CAUSE ("your table
is malformed" / "your assembleDay does not fit"), every id real; a runnable
`EXAMPLE`; `See also` in caps; the Claude stamp AFTER `See also`:
`Created <date> (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model <version>.` — updates append an `Updated <date> ...` line.

**6.2.2 HTML + manual.** Every public name: a page in `html/` + entries in
`Contents.m`, `helptoc.xml`, `html/index.html`; zero JavaScript/Java; every
equation a static image; the PDF manual is built by the ported ReportLab
pipeline with the glyph-coverage assertion and is verified by RASTERISING
AND LOOKING, never by "no exception" (v1 §1.3b/e lessons).

**6.2.3 The audit is the gate.** `docs_src/audit_tree.py` (WP-0, extended
per WP) must exit 0 before any ship. Checks are verified BOTH WAYS when
added: fire on a broken tree, silent on a healthy one.

### 6.3 Coding standards
The MATLAB coding standards linked in §0 rule 9, plus the house rules v1
paid for: gate behaviour on the structural reason, never a literal false;
public Name/Value options lowerCamelCase through a case-insensitive parser
with PartialMatching OFF; struct field access is case-sensitive — never rely
on parser case-folding for internal structs; no aliases; one name per thing;
`error`/`warn` only through `util.log`; verbosity only through
`util.verbosity`; no function longer than ~400 lines without a written
justification in this handover.

---

## 7. APPENDIX — REVIEW EVIDENCE (v1 file/line pointers)

- F1/F2: `+lsqvce/adjust.m`, `adjustBlock.m`, `adjustVCE.m` (2,479 lines),
  `Accumulator.m`, `Combine.m`, `BlockAccumulator.m`, `QRSystem.m`,
  `NormalEqSystem.m`.
- F3: `snoop.m` (flat-only; call-site audit in v1 HANDOVER §1.3e),
  `outliers.m` (mask handle, no loop), v1 HANDOVER §1.OPEN item B.
- F4: `diagnoseAccum.m` lines ~100–133 (`uFull = nc + sum(nl)`, dense
  `Rfull` via repeated `qr`), routed from `Results.report()`.
- F5: `BlockAccumulator.m` `preeliminate` (lines ~1108–1180): dense Nll per
  system, arc block-diagonality exploited only at elimination time; arc
  locals must span the system-wide B (validateArc).
- D4: `covModel.m` `psdModel` (line ~440: `Cfcn = toeplitz(r(1:numel(idx)))`,
  `info.filterMemory = 0`) + `BlockAccumulator.accumArc` (overlap discard is
  the FILTER correction; dense path batches whiten independently).
- D2: v1 HANDOVER §2 "KEY FACTS" (white start degenerate fixed point,
  measured in Python); Ellmer eqs. 6.5.19–6.5.26; Kvas eq. 3.14.
- D3: v1 HANDOVER §2 (numLags vs nf; the sharp nf >= m/2 + 1 rule);
  `covModel.m` psdUnderResolved warning.
- D7: v1 HANDOVER §2 "Function audit" (the redundancy seam operator-order
  fix) and `sharedProbe.m`/`traceEstimate.m`.
- D9: `adjustVCE.run` (scale VCE every pass) vs `adjustVCEjoint` help
  ("the two can converge to DIFFERENT models").
- F9: separate sweeps in `residuals.m`, `redundancy.m`, `redundancyLags.m`,
  `contribution.m`; the ~4 h/iteration estimate in `BlockAccumulator.m`
  PARALLELISM section.
- F10: `makeLoader.m` RUNTIME COVARIANCE OVERRIDE section; `jointModelOpt.m`
  ("three conventions that otherwise fail SILENTLY").
- F11: `makeLoader.m` 'commonParams'/'localParams' (14 columns, shape-checked
  only); `examples/demo_grace_dataLoader.m` column convention.
- Verified kernels worth porting unchanged: v1 HANDOVER §2 (Ellmer phases
  A–E MATLAB-confirmed; ar2arma reverted version; snoop downdate identities;
  QR/normal benchmarks; mapreduce exactness proofs in `validation/`).

*(end of handover)*

---

## 8. REVIEW R2 — FIVE-EXPERT SCRUTINY OF THE TOOLBOX (2026-08-16)

Blind re-review of the one `lsqvce/` tree, requested by the user on
2026-08-16, conducted INDEPENDENTLY of §2 (findings were not checked against
§2 and any overlap is convergent, not copied). Performed with assistance from
Claude (Anthropic), model Claude Fable 5. Numerical evidence:
`validation/v10review.py` (numpy, seed 20260816, exit 0; test IDs
T1-T7 are cited below). Thesis checks quote Ellmer (2018) directly.

### 8.1 The panel and its verdicts

- **E1 Adjustment theory & VCE.** The exact machinery is right: the block-arrow
  leverage algebra reproduces the dense joint inverse to 2e-17 (T1a), sum(h)=u
  and sum(r)=m-u hold and are CHECKED on every call, and rigorous Helmert in
  `varComp` (H_ij = tr(QN_i QN_j) + diag(n_i - 2 tr(QN_i))) is textbook-correct.
  Two REACHABLE estimator paths are mathematically wrong: F-01 and F-03.
- **E2 GRACE-FO processing.** Ellmer fidelity in the primitives is high
  (6.5.45-6.5.56 transcribed exactly, T7; degenerate-start detection in
  `psdVCE` is a genuine contribution beyond the thesis). The joint drivers
  deviate from the thesis exactly where arc-wise parameters matter — and
  Table 6.2 makes arc-wise parameters (initial states, ACC bias/scale) part of
  Ellmer's estimated vector, so the deviation hits the GRACE case squarely
  (F-03). The block path offers outlier DETECTION (`outliers` -> masks) but no
  automated reject-and-resolve driver; `snoop` lives on the flat path only —
  a capability observation, not a defect.
- **E3 Numerics & precision.** Kernels verified: `levinsonDurbin` recovers an
  AR(2) from its own theoretical ACV to 1e-8 and matches direct Yule-Walker
  (T4a/b); `armaACV`'s FFT autocorrelation equals the direct impulse
  autocorrelation (T4c); `dctI` 'ortho' is self-inverse to 1e-12 and 'ellmer'
  has cond 1.478 (T5). Concerns are policy floors and clamps that silently
  bias (F-04, F-10, F-11) and a dense per-arc inverse the docs deny (F-09).
- **E4 Software engineering & MATLAB.** The self-audit culture (15 static
  checks, clean), error messages in the standard's solution/problem form, no
  TODO debt, and <=100-char lines are exemplary. The one systematic breach is
  input validation: 0 of 79 public functions use `arguments` blocks against
  both the pinned standard and project rule 5 (F-07); duplication (F-08) and
  documentation drift (F-14) are the remaining hygiene items.
- **E5 Performance & robustness.** The square-root engine (`QRSystem`) with
  residual-mass bookkeeping on merge is the strongest piece of the tree. Hot
  loops carry avoidable O(N^3) work (F-09) and small tolls (`condest` on every
  `decorrFactor` success, growing task arrays); the silent skip-and-shift in
  the joint arc sweep (F-02) and unguarded negative variance factors (F-05)
  are the robustness holes.

### 8.2 Scope and coverage (read before trusting the findings)

Deep line-reads (~7 kLOC of 34 kLOC): `vce`, `varComp`, `redundancy`,
`streamContribution`, `leverage`, `leverageBlocks`, `traceEstimate`,
`sharedProbe`, `dctI`, `psd2acv`, `toeplitzCofactor`, `covComponent`,
`psdVCE`, `redundancyLags`, `arcVarianceFactors`, `adjustVCEjoint`,
`iterateNoiseModel`, `adjustVCE` (run loop, applyJoint, sweepStreamLags,
covGridSize, lagsPooledOverAxes, psdIsUsable), `BlockAccumulator`
(solve, preeliminate), `QRSystem`, `robust`, `levinsonDurbin`, `armaACV`,
`armaWhiten`, `decorrFactor`, plus `Contents.m`, README, setup, both runners'
headers and the joint-path tests. The whole tree passed the static battery
(audit_tree.py 15/15; TODO/rng/line-length/inv greps). NOT deep-reviewed:
`Combine` (1051 LOC), `diagnose` (1005), `report` (812), `snoop` internals
(667, help + downdate design read), `makeLoader` (552), `covModel` bodies
(495), `Results`, `BatchIndex`. Residual risk concentrates there; no defect
pattern pointed into them, but §8.3 cannot certify them. There is no MATLAB
in the sandbox: everything below is code-reading plus Python evidence; the
user's MATLAB run remains the ground truth (§0).

### 8.3 Findings

Severity: CRITICAL = wrong numbers on a reachable path; MAJOR = invalid
results can propagate or a binding rule is breached; MINOR = quality/clarity.

**F-01 (CRITICAL) `vce` 'traceMethod','hutchinson' divides by the wrong
quantity.** `+lsqvce/vce.m` lines 167-186 and 222-236 REPLACE the group
redundancy r_g = sum_{i in g}(1-h_i) with `streamContribution`'s
tr(N^-1 N_g^cc). Those are different objects: sum_g r_g = m-u while
sum_g tr(N^-1 N_g^cc) = nc (T2a: 687 vs 6 on the evidence problem). The
Helmert update sigma^2 = vtv/r then iterates against a denominator that is
orders of magnitude too small: measured, the exact path recovers truth within
its own resolution while the substituted path inflates 8.4x for one group and
diverges to 1e23 for the other (T2b/c). Ellmer's own redundancy is always
total-minus-absorbed (eq. 6.5.23: s_n = trace(R_bar V_n) with R_bar =
Sigma^-1 - hat), never the absorbed trace alone. The root cause is a
conceptual conflation that also appears in PROSE: vce.m's help (lines 36-38)
asserts "r_g = tr(N^-1 N_g^cc)", and adjustVCE.m's contribution comment
claims "r_g/sum(r_g) = tr(N^-1 N_g^cc)/(m-u)" — both false identities. The
default 'exact' path is UNAFFECTED. Fix: WP-32.

**F-02 (CRITICAL) `adjustVCEjoint` matches arc factors by position, and arcs
can vanish silently.** sweepArcs (adjustVCEjoint.m ~line 230) indexes
sigma2 by an enumeration counter; oneArc SKIPS empty arcs and non-PD
covariances (lines 264-266) without record. One skip on pass k shifts every
subsequent arc's factor by one from pass k on — factors applied to the wrong
arcs, no error, no warning. `vce` solved the identical problem with keyed
matching and a `groupMismatch` error; the joint driver has no equivalent.
Because the covariance model updates between passes, an arc's Sigma can flip
PD/non-PD MID-RUN, so this is a realistic trigger, not a corner case.
Fix: WP-34.

**F-03 (CRITICAL) Both joint paths ignore local parameters in the residual
and the projector.** adjustVCEjoint.m line 283 and adjustVCE.m
sweepStreamLags compute ebar = lbar - Abar*x; arc.B is never read and the
per-system local solution y is never subtracted, while the toolbox's own
`residuals`/`redundancy` define v = l - A*x - B*y. The redundancyLags
projector likewise uses the common-only design against the REDUCED common
factor. Ellmer's step 7 subtracts A*Delta_x for the FULL Table-6.2 vector —
which includes the daily arc-wise parameters — so this deviates from the
thesis precisely on the GRACE configuration (nl > 0). Measured at truth on
white noise with nl = 4: Omega_0 inflated 3.8x and alpha2 = 4.04 where the
consistent formula gives 1.08 (truth 1.0; T3). The estimated "noise" PSD
therefore contains the arc-parameter signal. WHY THE SUITE IS GREEN: the
joint tests are declared plumbing tests ("the test exercises the DRIVER's
plumbing", tests/test_ellmer.m jointAcc) — no accuracy assertion exists on
this path. Fix: WP-33 (+ the missing accuracy test).

**F-04 (MAJOR) `psdVCE` contradicts its own comment and biases zero-power
bins upward.** Line 169: `Shat = max(abs(Shat), MinAmplitude)` under the
comment "floor it". abs() of a noise-negative amplitude keeps the noise
MAGNITUDE: at a zero-power bin the expected injected power is ~2x what a true
floor admits (T6: ratio 2.01). Decide the policy (floor is the comment's
promise; abs needs a written justification), implement ONE of them, document
it, test it. Fix: WP-35.

**F-05 (MAJOR) `arcVarianceFactors` can return NEGATIVE variance factors.**
alpha2_m can be negative (s_n is a trace of an indefinite product;
redundancyLags says so); gamma2 = alpha2 .* sigma2 keeps that sign, and
sigma2Hat = gamma2/gammaBar^2 (line 113) ships it out unguarded. Any
downstream sqrt turns the arc weight complex; nothing warns. Fix: WP-35.

**F-06 (MAJOR) The RNG contract is broken in two places and the probe
machinery exists three times.** `streamContribution` (lines 162-168) and
`sharedProbe` (lines 91-96) document "default []: private stream, no global
side effect" but call plain randn — the GLOBAL stream advances. Three
Rademacher generators exist (traceEstimate: randi; streamContribution and
sharedProbe: sign(randn) with a +1 patch) — rule 3 says one. And in `vce`
'hutchinson' with the default empty seed, the final-pass probes differ from
the loop's, so vc.r is not the r the loop converged on. Fix: WP-36.

**F-07 (MAJOR) Input validation breaches the pinned standard and project rule
5.** The standard: "Use an arguments block to do validation. Avoid the use of
varargin to handle name-value arguments." The tree: 0 of 79 public functions
have an `arguments` block; 42 files use varargin+inputParser. Rule 5's "typed
signatures" is therefore unmet across the whole surface. Consequence beyond
compliance: inputParser's case-insensitive matching has masked an API-wide
casing split — 201 lowerCamel vs 53 UpperCamel Name-Value names, with the
SAME option spelled both ways in sibling functions ('groupBy'/'GroupBy',
'probes'/'Probes', 'verbose'/'Verbose', 'seed'/'Seed', 'maxIter'/'MaxIter').
The standard says UpperCamelCase. Fix: WP-39 (staged; largest package).

**F-08 (MINOR) Duplication against rule 3.** `makeKeyFcn` exists verbatim in
`redundancy` AND `streamContribution` (grouping presets can drift apart);
the batch-and-whiten arc loop exists three times (arcLeverageAndResidual,
arcCommonProbe, sweepStreamLags); the Rademacher draw three times (F-06).
Fix: WP-38, with an audit check so it cannot re-drift.

**F-09 (MINOR, performance + honesty) `redundancyLags` forms the dense
inverse it says it avoids.** Line 135: X = W\eye(N) — an N x N O(N^3) object
per arc, then O(N^2) per lag — while the help claims "the dense N x N Rbar
is NEVER formed on the main path" (true of Rbar; W^-1 is the same size) and
its own FUTURE-OPT note concedes the structured alternative. For a Toeplitz
Sigma the Sigma^-1 sub-diagonals have Trench/Gohberg-Semencul forms. Related
small tolls: `decorrFactor` runs condest on every successful chol (guard by
nargout), and `redundancy`'s task enumeration grows a struct array
element-wise. Fix: WP-37.

**F-10 (MINOR) `vce` line 188 floors r_g to 1 silently.** A group with true
r_g < 1 gets sigma^2 = vtv/1 — biased LOW — with no warning; contrast
varComp's explicit `estimable` handling. Fix: WP-40.

**F-11 (MINOR) `varComp`'s negative-component clamp destabilises.** s(bad) =
1e-6 multiplies the CUMULATIVE scale by 1e-6, i.e. next pass weights that
group 1e6 stronger — oscillation fuel. Freezing (s = 1) or flooring the
cumulative scale is the stable treatment. Fix: WP-40.

**F-12 (MINOR) Joint-mode convergence gate degenerates; 'joint' is
mislabelled.** In run(), shapeStable compares constant nLag values, so from
pass 2 the gate reduces to `balanced` alone — jinfo.maxDev (the spectral
convergence) is never enforced. Separately, adjustVCE 'vceMode','joint' does
NOT co-estimate the arc factors (applyJoint never calls arcVarianceFactors)
and still ALTERNATES scale (vce) with spectrum (applyJoint) inside a pass —
yet Contents.m tells users adjustVCE(...,'vceMode','joint') "wires all of
this up for you" as the 10-step co-estimation. Given the class's own warning
that alternating and joint schemes can converge to DIFFERENT Sigma, either
the code moves or the claim does. Fix: WP-40 (gate) + WP-41 (claims).

**F-13 (MINOR) Whiteness statistic crosses arc boundaries.**
iterateNoiseModel's maxAcf runs on residuals CONCATENATED across arcs while
its own pooledACV carefully never lets a lagged product cross an arc —
inconsistent, and a small bias for many short arcs. Fix: WP-40.

**F-14 (MINOR) Documentation drift sweep.** (a) redundancyLags help
references a 'Solver' option and a struct-W input that do not exist in code;
(b) the false r_g = tr(N^-1 N_g^cc) prose of F-01 (vce help, adjustVCE
comment, and the Contents `contribution` line calling the redundancy share an
"information share" — the Level-2 diagonal measures the OPPOSITE,
absorbed information); (c) the RNG claims of F-06; (d) the orphan HTML pages
lsqvce.SeparableProblem / lsqvce.recover that the tree's own audit flags as
"redirect stub?" — decide banner-or-delete; (e) vce's accFcn([]) -> unit-
weights convention lives in a code comment, not in the help INPUT block;
(f) two shipped runners (run_all_lsqvce is the gate, run_lsqvce the selector)
— the v1 open item "merge the runners" is restated here, not re-opened.
Fix: WP-41.

**F-15 (MAJOR, NEW 2026-08-19, session finding — not from the blind
review) Ellmer alternation is seed-fragile / bistable.** The closed-loop
spectral VCE (declared covariance fed back each pass; psdVCE update
alpha2.*Sprior with abs+floor+ZeroFreqRule) has an unstable branch: with
EXACT s_n (no probes) the recovered-PSD correlation across data seeds
swings 0.65..0.99 at 8 arcs x 96 rows, and MORE arcs on SHORT arcs make
it WORSE (16x96: corr 0.12; 16x128: -0.28; 16x192: one seed -0.18) — a
bistable fixed point, not estimation noise (SNR arithmetic puts pure
noise at ~0.97). LONG arcs stabilise (16x256: min 0.968 over 8 seeds).
Evidence: validation/wp33_loop_mirror.py (full-pipeline Python mirror).
Mechanism suspects: abs()+floor can lock a sign-flipped branch; the
ZeroFreqRule coupling. Disposition: registered here; candidate for WP-40
(damped/relaxed update or a monotone step).

### 8.4 What is RIGHT (so nobody "fixes" it)

The exact default paths are sound and verified: leverage block algebra
(T1, 2e-17), the sum(r) = m-u identity enforced per call, rigorous Helmert,
Levinson/armaACV/dctI/psd2acv (T4/T5), the per-lag kernel as transcribed
(T7), the QR square-root engine with the rho^2 residual row on merge, the
contamination guard in varComp, psdVCE's degenerate-start detection, the
skip-accounting in sweepStreamLags' warning, and the 15-check self-audit.
None of these may regress while the findings above are fixed; WP-20's parity
fixtures plus v10review.py are the guard.

---

## 9. R2 WORK PACKAGES AND PROMPTS

Same rules as §5: one row per package in the §1 R2 ledger, one NEW chat per
package, paste the §9.0 R2 PROMPT PREAMBLE (NOT §5.0 — the target tree
changed with the §1 reconciliation) then the package's prompt below, attach
the current zip. §6 guidelines bind every package. ORDER OF WORK:
32 -> 33 -> 34 -> 35 -> 36 -> 40 -> 37 -> 38 -> 41 -> 39 (39 last: it
touches every signature, so it must land on an otherwise-quiet tree; re-run
the FULL suite and audit after it).

### 9.0 R2 PROMPT PREAMBLE (paste before every R2 WP prompt)

> You are an expert MATLAB programmer and least-squares-adjustment
> specialist working on the ONE toolbox tree in the attached zip
> (a single folder, `lsqvce/`, since the 2026-08-17 restructure). Think step by step and motivate every answer. Read
> `HANDOVER.md` FIRST: §0 (rules), §1 (R2 ledger), §8 (the finding
> this WP fixes, with its evidence), §6 (test/doc guidelines) are binding.
> The Ellmer (2018) and Kvas (2019) theses are cited in
> `references/README.md`; their texts are in the chat project knowledge
> and, for repo-based sessions, in the gitignored `references/local/`
> (NOT redistributed in the public repo — WP-30 round note, 2026-08-18);
> equation numbers refer to Ellmer unless said otherwise. Rules that bind
> every WP: single point of access, no aliases, full help text with
> typed/dimensioned INPUT/OUTPUT and an ERRORS section, the Claude stamp
> (date, model) after `See also`, the MATLAB coding standards, tests per §6.1
> (every fix ships a test that FAILS on the old code), HTML doc registration
> per §6.2, no JavaScript/Java anywhere. There is no MATLAB in your sandbox:
> validate structurally, numerically in Python against exact formulae
> (extend `validation/v10review.py` where the WP names it), and with
> `python3 docs_src/audit_tree.py` (exit 0 required). The user runs
> `run_all_lsqvce` and returns the log as an uploaded .txt file — that run is
> the only gate. Ship the ENTIRE tree as ONE zip via present_files, update
> ONLY `HANDOVER.md` (the §1 R2 ledger row + a dated round note in
> the WP's own §9 entry; `HISTORY_v1.md` stays frozen as the v1
> history), and end with a wrap.

### WP-32 — Fix `vce` 'hutchinson': estimate the redundancy, not the trace

> Fix finding F-01 (HANDOVER §8.3): `lsqvce.vce` with
> 'traceMethod','hutchinson' substitutes tr(N^-1 N_g^cc) for r_g. Evidence:
> validation/v10review.py T2. Implement the correct streamed redundancy
> r_g = n_g - (local dof absorbed by g) - tr(N^-1 N_g^cc), where the local
> term is exact per system (a system's locals are exclusively its own; for
> sub-system groups derive the split from the block identities used in
> lsqvce.leverage and document the derivation in the help). n_g comes from
> the same sweep. Update vc.r and the help text (kill the false identity in
> lines 36-38). Tests: (1) hutchinson-vs-exact parity on a synthetic
> block-arrow problem within 3x the Monte-Carlo se; (2) sum_g r_g ~ m-u;
> (3) a regression test that FAILS on the old code (sigma inflation > 5x).
> Port v10review.py T2 into the MATLAB suite. No API change beyond the
> corrected numbers.

**ROUND NOTE (2026-08-18, WP-32 round 1, model Claude Fable 5; branch
`claude/wp-32`, riding housekeeping: WP-30 DONE flip + ci.yml
`cache: true`).** DERIVATION (documented in the streamContribution OUTPUT
block, validated in v10review.py T2d to 1e-8 against dense joint-inverse
h-sums, incl. sub-system splits): the streamed cc-attribution is the
symmetrized Schur split N_g^cc = Mcc_g - 1/2(Nlc_g'G + Nlc'G_g); the exact
per-group remainder from the leverage-block identities is
local_g = tr((Nll^-1 + G Sxx G') B_g'B_g) - tr(B_g'A_g Sxx G'), collapsing
to nl_s for whole-system groups; then r_g = n_g - local_g - tr(N^-1 N_g^cc)
equals n_g - sum_{i in g} h_i and sums to m-u. IMPLEMENTED:
streamContribution's ONE sweep additionally streams B_g'B_g + row counts
(arcCommonProbe extended), builds X = Sxx G'/G Sxx G' once per system (one
nc-solve), and outputs .n/.localAbsorbed/.r/.rSE per group (.trace keeps
its meaning: the Level-1 cc-attribution); vce 'hutchinson' divides by .r
in the loop AND reports it as vc.r (false identity in the help lines 36-38
and both comment blocks killed); v10review.py extended with T2d/T2e
(corrected split == exact; corrected iteration recovers truth [1.0396,
1.9766] vs truth [1,2] while the old path gives [8.4, 1.6e23]). TESTS in
test_contribution (suite already owned the hutchinson wiring; no new
files): (1) parity streamed-vs-exact r per obsType group (cuts across all
3 systems) within 3x MC se; (2) exact partitions (sum local = sum nl to
1e-6, sum r = m-u within 3x rSE); (3) regression test_vce_hutchinson_
regression_F01 with NO escape hatch — fails on the old code via
convergence, sum(vc.r) (old: ~nc=6 vs m-u), 20% sigma agreement, and the
explicit >5x-inflation guard; the old wiring test's false-identity prose
corrected. No API change: fields added, names kept.
**Evidence (2026-08-18):** CI green on 9d2d02f (audit incl. extended
v10review; matlab full gate). Postflight via MCP bridge on Windows R2026a:
ci_gate GREEN, 515/515, 0 incomplete, 201.5 s (18:39 local). PR #2 opened
with both verdicts quoted. Housekeeping ridden on this branch: WP-30 DONE
flip + ci.yml `cache: true` (fed1761).
**INCIDENT (2026-08-18, after the record push; pinned for every polling
session):** both runs on 669e3c4 hung >30 min INSIDE setup-matlab — the
first two runs to RESTORE the just-saved `cache: true` entry, started
simultaneously. Resolution per the Auth-note runbook (Actions RW): read
/jobs with step elapsed FIRST (>15 min in setup-matlab = hang, not a slow
suite); cancel + re-run the affected run(s); one re-run went green in
~90 s, PROVING the cache entry intact — so never delete the cache while a
twin restore succeeds; the second run needed one more cancel+re-run.
Merge followed with both contexts green.

### WP-33 — Joint paths: subtract B*y and absorb locals in the projector

> Fix finding F-03 (HANDOVER §8.3): adjustVCEjoint.oneArc and
> adjustVCE.sweepStreamLags compute ebar = lbar - Abar*x, omitting the local
> fit, and probe a common-only projector. Evidence: v10review.py T3; Ellmer
> step 7 + Table 6.2. Per system, recover y once per pass exactly as
> lsqvce.redundancy does (Nll_r \ (nl_r - Nlc*x)) and form
> ebar = W^-T (l - A x - B y). Extend the per-lag projector by the local
> block: R_bar = Sigma^-1 - Sigma^-1 [A B] inv(N_joint) [A B]' Sigma^-1,
> realised with the SAME block identities as lsqvce.leverage (Lc, Ll, G) so
> no joint inverse is formed; the shared probe still covers the common block,
> the local terms are exact and per-system. Both drivers route through ONE
> shared implementation (rule 3). Tests: alpha2 -> 1 at truth on white noise
> WITH nl > 0 (the accuracy test the suite lacks; tolerance from the probe
> se); PSD recovery corr > 0.95 on a coloured synthetic with locals;
> nl = 0 regression: results identical to the old path.

**ROUND NOTE (2026-08-18, WP-33 round 1, model Claude Fable 5; branch
`claude/wp-33`, riding housekeeping: WP-32 DONE flip + the setup-matlab
cache-restore incident runbook).** THE ONE SHARED IMPLEMENTATION (rule 3):
new package-private `+lsqvce/private/jointArcParts.m` — per arc
AtilBar = W^-T(A - B G), Bbar = W^-T B, ebar = W^-T(l - A x - B y); both
drivers call it. KERNEL: `redundancyLags` gains 'LocalBbar'/'LocalLl';
with them Rbar = Sigma^-1 - Rt Rt' - K K', K = W^-1(W^-T B)Ll^-T — the
EXACT local block (nl columns, same row-dot pattern as Rt), while the
shared probe now covers the REDUCED AtilBar, i.e. the operator U actually
factors; the block identity [A B]Njoint^-1[A B]' = Atil Sxx Atil' +
B Nll^-1 B' is validated dense in v10review.py T3c (2e-17) and the fixed
alpha2/s_0 in T3d. nl = 0 branch keeps the pre-fix expressions VERBATIM
(bitwise). DRIVERS: adjustVCEjoint.sweepArcs and adjustVCE.sweepStreamLags
recover y once per system and pass (Nll_r \ (nl_r - Nlc x), exactly as
redundancy), build G and Ll = chol(Nll,'lower') once, and route every arc
through jointArcParts; lagsPooledOverAxes threads the sliced local rows
symmetrically; new error ids lsqvce:adjustVCEjoint:localNotPD /
lsqvce:adjustVCE:localNotPD. TESTS (test_ellmer; no new suite files):
dense joint-projector parity incl. the DISCRIMINATING s_0 = N-(nc+nl)
(0.5 tol vs a 2.0 gap to the common-only N-nc); alpha2 ~ 1 at truth on
white noise WITH nl>0, pieces formed exactly as the drivers form them,
PLUS the pre-fix contrast alpha2Old > 1.5 asserted — the §6.1
fail-on-old-code lever (construction sanity in Python: new 1.004, old
2.14); bitwise nl = 0 regression (empty LocalBbar == no option, NO
tolerance); driver-level flat-PSD check on the white jointCfg fixture
with nl = 2 (leakage of the arc-parameter signal shows up as PSD
structure). The spec's coloured-synthetic PSD-corr test is DEFERRED to
the WP-33 round-2 note pending the postflight verdict on the four tests
above (the white-fixture driver check covers the leakage mechanism; a
coloured fixture needs a new synthetic and its runtime budget — decide
with the user if round 1 lands green).
**Round-1 CI red + fix (2026-08-18, same day):** the dense parity test
itself was wrong twice — it built G/Nll/Nred in the RAW metric while the
fixture's W is coloured (the reduction lives in the WHITENED metric, which
is exactly what the drivers' st.Nll_r/st.Nlc carry), and it asserted
s_0 = N-(nc+nl), which holds only for WHITE arcs (s_0 = trace(Rbar); on
the coloured fixture the exact value is trace(RbarEx) = 27.68, with the
pre-fix common-only projector a structural 3.18 higher — now the
discriminating assertion, 0.5 tol). The KERNEL and DRIVERS needed no
change: the whitened-metric limit agrees with the joint reference to
2.5e-16 (Python, seedable mirror). Note the round-1 CI also proved the
suite-wide safety of the driver rewiring: 487/488, the sole red being
this test's own reference.
**Evidence + incidents (2026-08-18, evening):** POSTFLIGHT via MCP bridge
GREEN — 519/519, 0 incomplete, 186.2 s (20:53 local), test_ellmer 58/58 —
recorded BEFORE the bridge became unavailable (user, same evening: bridge
gone for now; CI GREEN SUFFICES this round, local MATLAB testing resumes
later — a D.4 verification-mode note for FUTURE rounds, no debt for this
WP since its postflight is on record). CI: setup-matlab wedged across
FOUR attempts on 7f16773; escalation extends the incident runbook: the
cache list showed THREE entries with the IDENTICAL key (1.71 GB each —
every wedged run re-saved on its exit path), all deleted via the caches
API; a subsequent NO-cache attempt still wedged (800 s) → the cache is
EXONERATED, upstream transient (MathWorks CDN suspected; GitHub status
operational); attempt 4 ran green. Runbook addition: on repeated wedges,
LIST the caches — duplicate keys are a symptom, not the root cause; a
wedge without any cache present proves upstream. PR #3 opened with both
verdicts; zip-shipping rule amended in §0 rule 10 (user decision: no zip
when a PR is ready).
**ROUND 2 (2026-08-19, model Claude Fable 5; branch `claude/wp-33-r2`).**
User decisions this round: (a) the WP-33 spec's coloured PSD test lands
as round 2 (not WP-40); (b) CI runs ONCE per push — `ci.yml` triggers
rescoped to `pull_request` + push-to-`main`, the loop opens a DRAFT PR
after the first push (D8 loop text and CLAUDE.md amended; motive: the
double-run scheme doubled MathWorks downloads and raced cache saves).
Shipped: test_adjustVCEjoint_psd_recovery_coloured_with_locals on the
new jointCfgColoured fixture (one system, nl = 2, stationary AR(2)
[1 -0.75 0.4] — v10review T4's — one day, 6 arcs x 96 rows, nLag 10,
default non-delta start, MaxIter 6, Probes 150); truth on the
estimator's own grid via dctI(armaACV), corr > 0.95 asserted — fails on
the pre-WP-33 code (arc-parameter leakage + mis-sized s_n pull corr
under the bar). Calibration note: the psdVCE help's own study reached
0.945 from a flat start on ITS fixture; this fixture is sized richer
(arcs per frequency), and a near-miss on CI means MORE ARCS, not a
lower bar.
**Round-2a CI red + lesson (2026-08-19):** corr = -0.31 — not a
near-miss but a FIXTURE-SEMANTICS error: jointAcc's loader whitens with
cfg.cov, so with cfg.cov = ARMA the estimator saw ALREADY-WHITE
residuals and the correct Shat was FLAT (which is exactly why the white
flat-PSD test passes). Fix: generate coloured, DECLARE white —
cfg.userOpt.covOverride.pos = diag(1) threaded through jointAcc
('UserOpt'); synthLoader generates from cfg.cov (line 186) and declares
cm after the covOverride merge (147-150). Pinned as the fixture rule:
a spectral-recovery test needs truth-vs-declared DECOUPLING, or it
tests whitening, not estimation.
**Round-2b (same day): corr 0.65 with FIXED white weighting** — the
psdVCE update is relative to the PRIOR spectrum, so freezing the
declared covariance biases Shat toward the non-delta C0. Fix: the
fixture now runs Ellmer's TRUE alternation — jointAccMapped feeds the
driver's model.C back as the declared covariance each pass (covModel
type 'acf', nearest-lag lookup, zero beyond nLag-1; pass 1 white);
jointAcc reverted to its plain form (one way). MaxIter 8, Probes 200.
**Round-2c (2026-08-19): red = `noArcs`, root cause rank deficiency.**
The r2b 'acf' declaration of the nLag-truncated C is RANK DEFICIENT as a
Toeplitz beyond 2*(nLag-1) rows — covModel's own documented sharp limit
(nf >= m/2 + 1): at nLag = 10 the declared covariance dies at row 19 of
every 96-row arc, chol failed, oneArc skipped ALL arcs, the driver
raised noArcs. Cure per that same in-code doc: jointAccMapped maps
Shat = dctI(model.C) onto a 64-frequency grid (>= 96/2 + 1 = 49) and
declares covModel type 'psd' with psdMode 'toeplitz' — the toolbox's OWN
PD-safe PSD->Cfcn route (spectrum floor + resolution guard); lagInterp
DELETED (rule 8). Lesson for every model-feedback fixture: never declare
a truncated ACV directly; go through 'psd' with nf >= m/2 + 1
(commit 0db4b74).
**Round-2d (2026-08-19): r2c red at corr 0.638 — closed loop retired.**
Signal-to-noise arithmetic (round note is the record): the AR(2) truth
has CV 0.85 across 10 bins and the per-bin dof (~58) put PURE estimation
noise at corr ~0.97 — the deficit was SYSTEMATIC: the feedback fixture's
own declaration error (pchip-upsampled Shat -> psd2acv -> chol) compounds
over the passes; psdVCE's in-help reference reaches corr 0.9999 in the
OPEN-LOOP regime. The WP's spec tests the F-03 fix, not a fixture's
ability to close the feedback loop losslessly — so r2d measures open
loop: declare WHITE (cfg.override -> jointAcc UserOpt), generate AR(2),
C0 = [1; 1e-3...] (psdVCE's documented near-delta start, Sprior ~ flat),
MaxIter 1, Probes 300, 8 arcs x 96 rows: pass-1 Shat = alpha2 .* Sprior
IS the recovery. jointAccMapped and the psd-feedback plumbing DELETED
(rule 8); jointAcc gained the cfg.override pass-through. The closed-loop
mechanics lessons (r2a truth-whitening, r2b prior bias, r2c rank
deficiency nf >= m/2+1) stay pinned above — they are real operating
knowledge for production model-feedback loaders.
**Round-2e (2026-08-19): the measurement record and the honest bar.** r2d
red at corr 0.608 — the near-delta C0 IS the documented degenerate-start
throttle: C ~ delta annihilates every lag n > 0 in the sigma^2*C^n
weighting, alpha2 comes out FLAT in one pass (mirror: alpha2 = 1.69 in
every bin). Full-pipeline Python mirror (exact WP-33 kernel algebra,
Z->inf traces), measured across seeds and sizes: (a) fixed-white
iteration RISES then FALLS — peak corr ~ 0.96 at pass 1-2, FIXED POINT
corr ~ 0.64 without re-whitening, which is exactly where r2b/r2c/r2d all
landed (0.65/0.638/0.608): method behaviour, not a code bug; (b) the
ceiling is structural: pass-1 corr 0.92-0.96 (min 0.921 over 7 seeds)
and NEITHER more arcs (8->24), longer arcs (96->128) NOR more lags
(10->32) lift it above ~0.95 — the truncated-lag DCT estimator's bias
floor on this AR(2); (c) the spec's corr > 0.95 is therefore MEASURED
UNATTAINABLE for this estimator design in every tested regime; (d) corr
does not even DISCRIMINATE F-03: the pre-fix pipeline scores 0.55-0.95
across seeds (overlapping the fixed pipeline) because random B*y leakage
sometimes correlates by accident. CONSEQUENCE (deviation from the WP
spec, surfaced to the user in the round wrap for veto): the coloured
test ships as an ACCURACY FLOOR corr > 0.85 (MaxIter 1, default C0,
fixed white; measured min 0.921, margin 0.07) guarding the chain against
future regression; the §6.1 fail-on-old-code lever for F-03 is and
remains round 1's alpha2-at-truth test (gap 1.0 vs 2.1-4.0, no seed
overlap). OPTION pinned for WP-40: a lossless re-whitening alternation
(ARMA-fit feedback declaration) may reach 0.95+; decide there.
**Round-2e (2026-08-19): r2d red at corr 0.383 — open loop was the wrong
theory; closed loop restored at the STABLE volume.** The sigma^2*C.*omega
weighting (eq. 6.5.22-23) suppresses every lag n>0 under a near-delta
prior (psdVCE's degenerate-start warning says exactly this) — open-loop
flat-prior measurement is self-annihilating, 0.383 explained. A
full-pipeline Python mirror (validation/wp33_loop_mirror.py, exact s_n)
reproduces the r2b/r2c-range corr at 8x96 and EXPOSES finding F-15
(alternation bistability — registered in §8, queued on the WP-40 row).
Test design r2e: jointAccMapped restored (psd route, nf = 160 >=
256/2+1), 16 arcs x 256 rows (mirror-stable: min 0.968/8 seeds), MaxIter
10, Probes 300, seeds pinned (deterministic in MATLAB: the corr is a
sharp regression sentinel), corr over bins 2:end (psdVCE DEFINES
Shat(1)=Shat(2); the AR(2) truth peaks at DC, so including bin 1 only
imports the ZeroFreqRule as an error term).

### WP-34 — Keyed arc identity in `adjustVCEjoint`

> Fix finding F-02 (HANDOVER §8.3): arc factors are matched by enumeration
> position and skipped arcs silently shift the alignment. Key every arc by
> (sysID, day, obsType, arc) exactly as lsqvce.redundancy's keyFcn does;
> carry sigma2 as a keyed map pinned on pass 1; a key set that changes
> between passes raises lsqvce:adjustVCEjoint:arcMismatch (mirror vce's
> groupMismatch text); an arc skipped THIS pass keeps its previous factor and
> is counted in a skip summary. Tests: induce a mid-run skip (make one arc's
> covariance non-PD on pass 2) and verify (a) the other arcs' factors are
> byte-identical to the no-skip run, (b) the skip is reported; plus the
> mismatch error test.

### WP-35 — Amplitude and factor sanity: psdVCE floor policy, negative sigma2

> Fix findings F-04 and F-05 (HANDOVER §8.3). psdVCE: decide floor vs abs for
> negative Shat (the comment promises a floor; v10review.py T6 measures the
> abs bias at 2x on zero-power bins), implement the decision, align comment
> and help, and add a test pinning the behaviour at a zero-power bin.
> arcVarianceFactors: a negative sigma2Hat_m must not leave the function —
> floor at a documented epsilon relative to gammaBar^2, warn with the arc
> indices and the count, and return a .clamped mask; test with a constructed
> negative alpha2. Do not change converged results on healthy inputs (add a
> no-op regression test on the standing fixture).

### WP-36 — One probe generator, honest RNG contract

> Fix finding F-06 (HANDOVER §8.3). Create ONE private Rademacher/Gaussian
> probe generator (rule 3) used by traceEstimate, streamContribution,
> sharedProbe; empty seed draws from a PRIVATE substream derived once (never
> the global stream), matching what the help already promises — update the
> three help blocks to the one true contract. In vce 'hutchinson', the final
> pass reuses the LOOP's probe state so vc.r is the converged r. Tests:
> global-RNG invariance (state before == after with default seed);
> reproducibility (same seed -> identical traces across the three callers);
> vce final-vs-loop consistency.

### WP-37 — `redundancyLags` without the dense W^-1 (+ small hot-loop tolls)

> Fix finding F-09 (HANDOVER §8.3). Replace X = W\eye(N) by structured
> computation of the Sigma^-1 sub-diagonals: for the Toeplitz covariance case
> use the Trench/Gohberg-Semencul representation (O(N^2) total, O(N) memory
> beyond the factor); for a general W fall back to per-lag solves without
> materialising W^-1. Keep the WantRbar testing path. Make the help's memory
> claim true. Also: decorrFactor computes info only when nargout > 1;
> redundancy's task enumeration preallocates. Tests: bit-agreement of s_n
> with the current implementation on the standing fixture (1e-12), plus a
> SPEED test (§6.1) demonstrating the removal of the O(N^3) term on a long
> arc; no API change.

### WP-38 — Deduplicate: keys, arc batching, and an audit lock

> Fix finding F-08 (HANDOVER §8.3). One private makeKeyFcn used by
> redundancy, streamContribution (and any other grouping site found by
> grep); one private batch-and-whiten arc iterator used by
> arcLeverageAndResidual, arcCommonProbe and sweepStreamLags (parameterised
> by the per-batch callback); the probe generator is WP-36's. No public
> surface change, no behaviour change: the suite must be green with byte-
> identical fixture numbers. Add an audit_tree.py check that fails on a
> second definition of the key presets or the batching loop, so the
> duplication cannot return.

### WP-39 — `arguments` blocks and Name-Value casing across the surface

> Fix finding F-07 (HANDOVER §8.3): 0 of 79 public functions validate with
> `arguments`; Name-Value casing is split 201/53 with same-name pairs in both
> casings. Migrate the PUBLIC surface to arguments blocks (typed, sized,
> defaulted — rule 5's "typed signatures"), replacing varargin+inputParser;
> unify every Name-Value name to UpperCamelCase per the standard. Because
> inputParser was case-insensitive, existing user code in either casing keeps
> working under arguments-block name-value matching only for EXACT names:
> state the one canonical casing in every help block and the manual, and note
> the break for previously-sloppy casings in the release notes. Stage it:
> (a) the front doors (adjust, adjustBlock, adjustVCE configurators,
> adjustVCEjoint, iterateNoiseModel), (b) the VCE/diagnostics layer, (c) the
> primitives; full suite + audit after each stage. Tests: signature tests per
> migrated function (wrong type/size/name errors with the documented
> identifier), plus the untouched-numerics regression on the standing
> fixture. This package lands LAST.

### WP-40 — Method polish: gates, clamps, whiteness

> Fix findings F-10, F-11, F-12 (gate half), F-13 (HANDOVER §8.3). vce:
> replace the silent rg = max(rg,1) with varComp-style estimability (freeze
> + warn below a documented threshold). varComp: negative-component
> treatment becomes freeze (s = 1) or a floored CUMULATIVE scale — pick,
> motivate in the help, warn once with the group list. adjustVCE joint gate:
> convergence additionally requires hist.jointMaxDev(end) <= outer Tol;
> history documents it. iterateNoiseModel: maxAcf pools per-arc (no lagged
> product crosses an arc), mirroring pooledACV. Tests: one per change,
> each failing on the old code; fixture regression for the healthy paths.

### WP-41 — Documentation truth sweep

> Fix findings F-12 (claims half) and F-14 (HANDOVER §8.3), touching DOCS
> only (help text, Contents.m, HTML regeneration via the v1 tooling): correct
> the adjustVCE-vs-adjustVCEjoint 'joint' claims (alternating spectrum+scale,
> no arc factors, vs the 10-step co-estimation); kill every r_g =
> tr(N^-1 N_g^cc) assertion and rename the Level-1 share to what it is (a
> redundancy share) with one sentence contrasting it against the Level-2
> absorbed-information diagonal; fix the RNG sentences to WP-36's contract;
> remove the 'Solver'/struct-W ghost from redundancyLags; document
> accFcn([]) in vce's INPUT block; decide banner-or-delete for the two
> orphan HTML pages and make audit check 4 enforce the decision; retarget
> in-code prose references to "HANDOVER" (verbosity, adjustVCE, redundancy,
> vceHutchinson, varComp, poolGroups, tests) to `HISTORY_v1.md` where they
> mean the v1 history. Re-run the
> doc build and the audit; no .m logic changes (help-only edits inside .m
> files are in scope).
