# lsqvce v2 — Handover: Review, Redesign, and Work Packages

**THIS IS THE ONLY STATUS-AND-PLANS DOCUMENT for the v2 toolbox. Do not create
a second one.** The v1 state of this tree is FROZEN for feature work as of
2026-07-19; its own `HANDOVER.md` records that and remains the authoritative
history of v1. Everything about v2 — the design review that motivated it, the
architecture, the task list, and the implementation prompts — lives HERE.

Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Fable 5.
Last substantive update: 2026-08-21 (evening) — WP-52 DONE (PR #31
merged); WP-53 ROUND 1: demo_gravity returns `score`, and the USER's
run settled the §10.13 quality debt — sigma0 0.9998, stream ratios
0.996/1.000/1.019, mean z^2 1.087, anomaly at 1.33 sigma; pinned in the
demo's ACCURACY block. Round note §10.15.
Preceding entry: Last substantive update: 2026-08-21 (later) — WP-52 ROUND 1 on branch
`claude/wp-52`: ragged schedules NATIVE in makeLoader (per-day arc
resolver, .arcId addressing with absent ids as data gaps,
'localsPerArc' deriving nl and owning the offsets); demo rewired and
raggedClosure DELETED. Round note §10.14. Bridge still down: no local
postflight on this WP either.
Preceding entry: Last substantive update: 2026-08-21 — WP-50 DONE (PR #29 merged); WP-51
ROUND 1 on branch `claude/wp-51`: demo_gravity (GRACE-like, energy-balance
KBR stream, per-system local counts 2/3/1, ragged arcs absent by ID,
ragged arc lengths), SH generator moved to examples/ as the ONE
implementation. CI green incl. an end-to-end demo run; LOCAL POSTFLIGHT
MISSING (bridge dropped) and the recovery quality is UNMEASURED — §10.13.
Preceding entry: 2026-08-20 (night, later) — WP-49 DONE (PR #28
merged); WP-50 ROUND 1 on branch `claude/wp-50`: the QR sigma0=0
mechanism found (merge-fold normal-form cancellation), guarded with a
three-anchor measured threshold (good 13.4 @ 8e-5 / destroyed ~19 /
pinned 15), rotation backstop, healthy path locked warning-free;
postflight 593/593. Round note §10.12.
Preceding entry: 2026-08-20 (night) — WP-48 tier (a) DONE (PR
#27 merged), tier (b) BACKLOG per user; WP-49 ROUND 1 on branch
`claude/wp-49` (honest vPv: NaN never a plausible number, ONE
cancellation id per namespace — BREAKING retirement of the negativeVtv
name — residuals supersede logic + flag; postflight 591/591; four pinned
findings incl. the QR sigma0=0 deep-regime observation). Round note
§10.11.
Preceding entry: 2026-08-20 (evening, later) — WP-48 tier (a)
ROUND 1 on branch `claude/wp-48`: shDesign + test_physOracle (8/8 first
pass; full runner 533/533 locally), measured numbers pinned per §6.2.3,
authority rows added; WP-49 queued OPEN. Round note §10.10.
Preceding entry: 2026-08-20 (evening) — WP-47 DONE (PR #25
merged; merge 8325a26); WP-42 reconfirmed BACKLOG by the user; §10.9
fixture-reality review written (closed-loop planted-truth, but Gaussian
designs — no physical structure, no real data) with proposals WP-48
(physically structured oracle) and WP-49 (authoritative vPv), both
awaiting user decision.
Preceding entry: 2026-08-20 (later still) — WP-46 DONE (PR #24
merged) and WP-47 ROUND 1 on branch `claude/wp-47`: gate hardening
(declared-vs-executed rule 6, warning-inventory rule 7 + live probe,
instruments/diagnosis doctrine, §6.1.4 status column); check 17's first
live catch (574→581). Round note §10.8.
Preceding entry: 2026-08-20 (later) — WP-46 ROUND 1 (README
de-drifted; §2.8 of the external geoMap handover adopted: ACCURACY block
+ §6.2.4 sync gate + audit check 17, red-tested); WP-45 REJECTED by the
user (recorded §10.7).
Preceding entry: 2026-08-20 — WP-44 DONE (PR #22 merged) and
POSTFLIGHT SETTLEMENT: bridge back, full local run 574/574 GREEN on
039041c incl. test_parallel (31/31, real pools) and test_speed (11/11)
— ALL D.4 debts cleared (§10.6). R3 track fully closed.
Preceding entry: 2026-08-19 (later still) — WP-44 ROUND 1 on
branch `claude/wp-44` (run provenance: runMeta composer, Results.meta,
producer call sites, report 'Provenance' NV, 3 tests, 2 HTML pages);
round note §10.5; ledger row 44 IN PROGRESS.
Preceding entry: 2026-08-19 (late) — WP-43 DONE (PR #20 merged);
r2 same day: tol: debt 233 -> 0 (all sites annotated, check 16d full
enforcement, red-tested); WP-44 prepared as an executable prompt (§10.3),
ledger row OPEN. Round notes §10.4.
Preceding entry: 2026-08-19 — R3 OPENED (WP-43, test-category
scrutiny per the adapted BEST_PRACTICE §3.3 taxonomy): §6.1.2 rewritten to
the six categories + two cross-cutting tags, §6.1.4 authority table added,
§10 review record (findings T-01…T-08), audit check 16 (countability gate),
571 declared test points, 8 contract gaps closed, 3 contract defects fixed
(T-03a), WP-44 proposed. Ledger: §10.
Preceding entry: 2026-08-18 — WP-30 ROUND 1 (container session,
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
order of work: 32 → 33 → 34 → 35 → 36 → 40 → 37 → 38 → 41 → 39.
**R2 TRACK CLOSED (2026-08-19, WP-39 merge ff1ef9e):** all fifteen
findings F-01…F-15 addressed (F-15 closed UNDERSTOOD+GUIDED, §8
disposition). Every public front door now carries a typed
`arguments` signature with canonical UpperCamelCase names, or a pinned
pass-through deviation (adjust, robust, solveNonlinear, the fluent
configurators) with the same canonical surface. Open beyond the track:
the runner migration (first green full run settles all postflight
debts WP-33-r2…39) and the WP-42 backlog.
**Order amendment (2026-08-19, user decision, option A):** WP-37 was
executed BEFORE WP-40 — a session process error (the order line above was
not consulted at WP start; only the ledger rows were), caught by the USER
after WP-37 had gone green. Technical impact assessed as none (WP-37 was
confined to redundancyLags internals + two tolls, no API change; WP-40's
targets are disjoint). Approved actual order: 36 → 37 → 40 → 38 → 41 → 39.
**Process rule, binding for every session incl. parallel ones: read THIS
order line before opening any WP — the ledger alone does not carry the
order.**):

| WP | title | depends on | status |
|----|-------|-----------|--------|
| 32 | Fix `vce` 'hutchinson': redundancy, not trace (F-01) | — | DONE (2026-08-18, PR #2, merge e167cee; CI green both contexts after a setup-matlab cache-restore hang was resolved per the incident runbook, postflight green 515/515 on Windows R2026a; evidence + round note §9 WP-32) |
| 33 | Joint paths: locals in residual + projector (F-03) | — | DONE (2026-08-19: round 1 PR #3 merge 7d7c14c, round 2 PR #4 merge 24cfe6c; CI green 489/489; user approved the 0.85 accuracy floor + F-15→WP-40 scope by merging; postflight: round 1 green pre-outage, round 2 = D.4 DEBT until the bridge returns → SETTLED §10.6; six r2 lessons + F-15 pinned in §9 WP-33) |
| 34 | Keyed arc identity in `adjustVCEjoint` (F-02) | — | DONE (2026-08-19, PR #5, merge af440dd; CI green 491/491; deviation approved by merge: byte guarantee on the kept factor, RelTol 0.15 renormalisation bound; POSTFLIGHT DEBT → SETTLED §10.6; round note §9 WP-34) |
| 35 | psdVCE floor policy + negative-factor clamp (F-04, F-05) | — | DONE (2026-08-19, PR #6, merge 827ee6b; CI green 493/493 — the run that also validated the self-healing setup live: attempt-1 wedge killed at 720 s, attempt-2 restore 90 s; POSTFLIGHT DEBT → SETTLED §10.6 until the runner/bridge returns; round note §9 WP-35) |
| 36 | One probe generator, honest RNG contract (F-06) | — | DONE (2026-08-19, PR #7, merge 74ab106; CI green 498/498 — the WP's own changes green on the FIRST pass, the sole round-1 red was the test's own traceEstimate-returns-a-scalar mistake; POSTFLIGHT DEBT → SETTLED §10.6; round note §9 WP-36) |
| 37 | `redundancyLags` without dense W^-1 + hot-loop tolls (F-09) | — | DONE (2026-08-19, PR #8, merge 66ffb81; CI green FIRST pass 501/501 incl. the N=4000 speed lever; GS/Levinson validated in Python to 1e-14 pre-port; POSTFLIGHT DEBT → SETTLED §10.6; round note §9 WP-37) |
| 38 | Deduplicate keys / arc batching + audit lock (F-08) | 33,36 | DONE (2026-08-19, PR #11, merge 74d2f64; CI green FIRST pass 507/507 byte-identical; audit lock 10b live, red-tested; GroupBy error-id unification approved by merge; POSTFLIGHT DEBT → SETTLED §10.6; round note §9 WP-38) |
| 39 | `arguments` blocks + Name-Value casing (F-07) | all | DONE — 2026-08-19, staged rollout complete: (a) a1+a1-2+a2 (PRs #13, #14), (b) b1+b2 (PRs #15, #16), (c) c1+c2 (PRs #17, #18), merge ff1ef9e (round notes + pinned deviations §9 WP-39); POSTFLIGHT DEBT → SETTLED §10.6 |
| 40 | Method polish: gates, clamps, whiteness (F-10…F-13; + F-15 alternation bistability) | 32 | DONE (2026-08-19, TWO rounds: r1 F-10…F-13 PR #9 merge 736807d, CI 505/505; r2 F-15 PR #10 merge 8ca16b9, CI 507/507 after six diagnosis rounds with the psd2acv-faithful mirror; ALL THREE r2 decisions ACCEPTED by the user 2026-08-19 — detailed in the DECISION RECORD below the round notes (§9 WP-40); POSTFLIGHT DEBT → SETTLED §10.6 both rounds; round notes §9 WP-40) |
| 41 | Documentation truth sweep (F-12, F-14) | 32,33,36,40 | DONE (2026-08-19, PR #12, merge 739e578; CI green 507/507, zero fixture movement; stub rule + all claim fixes approved by merge; round note §9 WP-41) |

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

**6.1.2 The test categories — every test point declares which it ships.**
*(Rewritten 2026-08-19, WP-43: the adapted taxonomy from
references/BEST_PRACTICE_v4-1.md §3.3, user decision. The former four-kind
list is superseded; METAMORPHIC is no longer a category of its own — its
properties classify as contract (order/batch invariance where promised),
robustness (worker count) or vectorisation (split/merge), which is where
the WP-29 suite's points now sit.)*

| Category | Asserts |
|---|---|
| contract | The promises in the help text: output shapes and types, argument rejection with the DOCUMENTED identifier and only that one, determinism where promised, error messages naming the offending field. Documented-id <-> code-id agreement is ALSO enforced statically (audit check 14). |
| reference | Agreement with an outside authority — a published example, a second toolkit/language, real data. |
| precision | The numerical claim, with an explicit tolerance AND a comment saying where the tolerance came from. House starting points: exact-identity <= 1e-10 relative; factorisation round-trips <= 1e-12; Monte-Carlo quantities within 4/sqrt(Z) (a 2.93-sigma bar, NOT containment). Never loosen a tolerance to make a test pass — that is a finding. |
| speed | A RELATIVE budget against a recorded baseline on the same machine (medians; absolute seconds are LOGGED, not asserted; speed suites marked so the runner reports them separately). |
| robustness | Behaviour at and beyond the edge of validity: degeneracy CONSTRUCTED FROM THE ALGEBRA, conditioning, boundary and empty inputs, non-finite propagation, reproducibility under worker count. |
| vectorisation | A batched call equals the same inputs one at a time — the single cheapest reference test in a vectorised language, and the one most often assumed. Each such property states whether the expectation is BITWISE or eps-level, and why. |

Two CROSS-CUTTING TAGS, not categories: `diagnostics` (does every measured
number reach the report) and `provenance` (does every result carry the
metadata saying where it came from).

MECHANICS (binding, enforced by audit check 16):
- Every test point carries `% categories: <list>` as the FIRST line after
  its function signature, optionally `| tags: <list>`. Only the six
  category names and the two tag names are legal vocabulary.
- The unit x category matrix is DERIVED per run from the `lsqvce.<name>`
  references in each test body — no second bookkeeping file to drift —
  and printed by the audit; gaps show as rows, not as absences nobody
  notices.
- A public unit with zero `contract` coverage is RED unless
  `tests/category_exemptions.csv` carries a reasoned row. An exemption is
  a claim that the test is IMPOSSIBLE, not inconvenient; `contract` may
  never be exempted for a function that validates its arguments.
- TOLERANCE PROVENANCE: every 'AbsTol'/'RelTol' carries a `tol:` comment
  (same line or up to three lines above) naming its origin. Enforced RED
  in the swept suites listed in the audit (start: test_signalCore,
  test_precision); elsewhere the unannotated count is REPORTED as the
  standing debt number, driven to zero suite by suite (WP-39 rollout
  pattern).

**6.1.2b Diagnosis order** *(added 2026-08-20, WP-47; external handover
§2.10)*. When a number is wrong: CONFIGURATION (replay the exact frozen
script, never a recollection of it), then CRITERION (re-derive against the
authority — corrected, never widened until it passes), then CODE —
cheapest repair first. And never loosen a tolerance to make a test pass;
that is a finding (already binding via §6.1.2). Cluster rule: unrelated
failures sharing one symptom are SESSION STATE, not code — check the
state before reading a single diff.

**6.1.3 Regression discipline.** Every defect fixed during v2 gets a named
regression test citing the WP section that records it (the v1 habit that
worked). Randomness is seeded per test. Tests that need a pool carry an
assumption filter and are counted as incomplete, not skipped silently.

**6.1.4 The authority table (added 2026-08-19, WP-43).** One table listing
every external authority the `reference` tests consume, what it certifies,
and where it is consumed — so gaps show up as empty rows rather than as an
absence nobody notices. Amend this table in the same edit that adds or
retires an authority.

*(Status column added 2026-08-20, WP-47 — adopted from the external
geoMap oracle register: dated confirmations AND refutations, ☑ consumed
and confirmed / ☐ open; the refutation discipline is the valuable half.)*

| Authority | Certifies | Consumed by | Status |
|---|---|---|---|
| Ellmer (2018) doctoral thesis | dctI norm convention, PSD-VCE update law, joint-mode default constants | test_ellmer (17 references), §8 review evidence | ☑ consumed since the Ellmer WPs; PARTIAL REFUTATION 2026-08-19 (WP-40 r6): the closed-loop alternation as such does not beat open-loop pass-1 quality — F-15, production guidance pinned |
| Kvas (2019) doctoral thesis | design record for the VCE iteration/combination architecture | design sections §2–§3 | ☐ no direct numeric consumption — empty row, deliberate |
| Baarda / Pope test theory (textbook) | normCrit, delta0/MDB, tau bounds, snooping thresholds | test_lsqvce (test_normCrit_textbook, test_delta0_textbook, tau tests), test_lsqvce_opts (baarda_vs_pope) | ☑ confirmed at the quoted decimals; SCOPE FINDING (v1, 2026-07): Pope's tau goes blind under contamination — the README headline finding |
| Closed-form ARMA(1,1)/AR ACV algebra | armaACV, levinsonDurbin, psd2acv identities | test_precision, test_lsqAccum, test_signalCore | ☑ confirmed 1e-9…1e-12 (tol: comments name the sites) |
| Python mirrors (second language): v8verify.py, v9verify.py, v10review.py, wp33_loop_mirror.py | measured tolerances and loop behaviour pre-port | tol: comments in test_signalCore/test_streamUtils/test_precision; §9 WP-33/40 round notes | ☑ v8: 2e-16 (whitening); WP-37 GS/Levinson 1e-14 pre-port; wp33 mirror: the F-15 instrument — NOTE its own r6 lesson: an unfaithful mirror declaration path produced false accuracy promises, faithfulness is part of the oracle |
| GROOPS / SINEX NEQ files (real data exchange) | neqExchange round-trip | WP-27 suite | ☑ round-trip green since WP-27 |
| Analytic SH algebra (fully-normalised Legendre: addition theorem, closed forms P̄10/P̄20/P̄22) | the shDesign partial generator of the WP-48 physical oracle | test_physOracle (reference tests) | ☑ 2026-08-20: addition theorem to 4.3e-14 at nmax 12, closed forms to 1e-14 (measured, user machine) |
| Published GRACE-FO monthly solution (release to be named) | end-to-end sign/magnitude/pattern of a real month through the SINEX path | WP-48 tier (b) | ☐ user to supply files + name the release |
| Real GRACE-FO L1B-derived NEQ month | the same, input side | WP-48 tier (b) | ☐ user to supply |

### 6.2 Documentation guidelines

**6.2.1 Help text template** (binding for every function/class/method):
H1 line `%LSQVCE.NAME  One-sentence purpose.`; call forms; a MOTIVATED
description (why, and any measured numbers backing claims); `INPUT` and
`OUTPUT` blocks with NAME, TYPE, DIMENSION, DESCRIPTION for every argument
and field; `OPTIONS` with defaults; an `ACCURACY` block *(added 2026-08-20,
WP-46, adopted from the external geoMap handover §2.8)* stating the
function's numerical claim, what it was measured against — naming the
authority by its §6.1.4 row — and the tolerance achieved (where the
measurement could only be taken in MATLAB, say so and give the re-derivation
procedure); `ERRORS` grouped BY CAUSE ("your table
is malformed" / "your assembleDay does not fit"), every id real; a runnable
`EXAMPLE`; `See also` in caps; the Claude stamp AFTER `See also`:
`Created <date> (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model <version>.` — updates append an `Updated <date> ...` line.
ROLLOUT RULE for the ACCURACY block: binding for every NEW function and
every function TOUCHED from 2026-08-20 on; NO bulk migration of the ~80
existing help blocks — the source handover's own record (a 94-file
help-text mutation producing three defects found only by reading the
applied diff) is the pinned rationale, and it matches our WP-43 r2
bulk-annotation lessons.

**6.2.2 HTML + manual.** Every public name: a page in `html/` + entries in
`Contents.m`, `helptoc.xml`, `html/index.html`; zero JavaScript/Java; every
equation a static image; the PDF manual is built by the ported ReportLab
pipeline with the glyph-coverage assertion and is verified by RASTERISING
AND LOOKING, never by "no exception" (v1 §1.3b/e lessons).

**6.2.3 The audit is the gate.** `docs_src/audit_tree.py` (WP-0, extended
per WP) must exit 0 before any ship. Checks are verified BOTH WAYS when
added: fire on a broken tree, silent on a healthy one.
*Extended 2026-08-20 (WP-47, adopted from the external geoMap handover
§2.9), binding:*
- every NEW audit/gate check ships its fault-injection red-test IN THE
  SAME ROUND (plant the defect, see it fire, restore, see green — recorded
  in the round note; checks 16/17 and gateVerdict rules 6/7 are the
  precedents);
- no count, structure or coverage claim from an instrument written in the
  same session as the change it checks, unless the instrument was shown to
  fail on a known defect — otherwise mark the claim PROVISIONAL in
  writing;
- a second gate must read something the first does not (the Python audit
  and the MATLAB runner deliberately count declared points with SEPARATE
  implementations of one definition);
- a bulk mutator is proven on a scratch copy and its APPLIED DIFF is read
  before it touches the tree; one case form is never enough; a one-shot
  mutator is not shipped (our WP-43 r2 lessons are the local evidence).

**6.2.4 Documentation sync gate** *(added 2026-08-20, WP-46; adopted from
the external geoMap handover §2.8.2, adapted)*. Every change ships its
documentation in the SAME change — all surfaces, every time: help text,
the HTML page, `Contents.m`, and (where a claim lives there) `README.md`
and this handover. The audit cross-reads the surfaces and fails on
disagreement; initial scope is check 17 (README title version ==
Contents.m version, README suite/point-count claim == the tree's counts),
extended surface-by-surface as claims are made machine-checkable rather
than by trusting prose. Two adopted principles pinned for those
extensions: (a) COUNT COMPLETENESS IN THE BUILT ARTEFACT, not the pipeline
feeding it — a doc builder must report arguments documented vs arguments
RENDERED, counted in the emitted HTML; (b) EVERY EXAMPLE IS CODE —
extract and lint the help-text EXAMPLE blocks at build time so a stale
example fails the build, not the reader (staged: report the
extraction/lint count first, flip to red after the sweep — the WP-39/43
rollout pattern). A generated manual is verified by rasterising and
looking (already §6.2.2).

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
ZeroFreqRule coupling. **DISPOSITION UPDATE (2026-08-19, WP-40 r2): both
suspects EXONERATED by measurement (the WP-35 floor does not cure it —
some seeds worsen; ZeroFreqRule off and higher nLag only SHIFT the
failures). Mechanism identified: stochastic-approximation instability of
the undamped multiplicative step (per-bin estimation noise re-amplified
through the re-whitening). MITIGATION SHIPPED, opt-in: psdVCE 'Gamma'
(damped update) + adjustVCEjoint 'StepSchedule','decay' (Robbins-Monro
Gamma_it = 2/(2+it)) + 'AvgTail' (Polyak tail averaging). Mirror-measured
across the grid: worst-seed corr from -0.66 to ~0.9 at nearly unchanged
ceiling; defaults are bitwise the old path. Residual risk: the undamped
default remains bistable — RECOMMENDATION for production closed loops:
'StepSchedule','decay','AvgTail',4. Severity after mitigation: user call
(approve by merging PR; suggest MITIGATED).**
**DISPOSITION UPDATE 2 (2026-08-19, WP-40 r6, psd2acv-FAITHFUL mirror):
the r2 mitigation's accuracy promises were an artefact of an unfaithful
mirror DECLARATION path (zero-extended dctI acv vs MATLAB's full m-lag
psd2acv transform). Faithful result: the closed loop does not beat
OPEN-LOOP pass-1 quality at ANY tested volume (open loop 16x256: min
0.96/med 0.97; closed stabilised: P(corr<0.8) = 0.31 at 8x96, 0.75 at
16x256 — the old "stable volume" was the artefact); decay+tail is
DAMAGE LIMITATION (prevents the ~0.2 collapse), not restoration. Also
falsified as cures: lag windows (Lc 10/20/30), pchip-vs-linear
interpolation, log-shrinkage of the declaration (beta 0.5/0.7).
PRODUCTION GUIDANCE (pinned in the adjustVCEjoint help): open loop or
<= 3 damped passes. Severity call and any deeper method work (e.g. a
monotone/trust-region alternation) = user, asked at PR #10.**

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
**RUNBOOK UPGRADE (2026-08-19, user decision after the THIRD wedge, on
WP-35's run):** ci.yml now SELF-HEALS -- setup-matlab attempt 1 has
`timeout-minutes: 12` + `continue-on-error`, a second attempt fires
automatically when the first did not succeed, the job carries a 45-min
backstop, and the check name `matlab` is unchanged (protection
untouched). Manual cancel+re-run remains only for double wedges. Options
recorded for the future if wedges persist: (a) run the matlab job in the
mathworks/matlab Docker image -- removes the MathWorks download entirely
but needs a batch-licensing token as a repo secret (USER decision +
setup); (b) a self-hosted runner on the user's PC is NOT recommended:
on a PUBLIC repo, fork PRs would execute code on that machine.
**RUNNER DECISION (2026-08-19, user):** the self-hosted-runner migration
(private repo + runner on the user's Windows PC = ONE full-suite gate on
the target platform, incl. test_parallel; triggerable from any Claude
session via push / workflow_dispatch / re-run — GitHub is the
dispatcher, no bridge needed) is APPROVED and DEFERRED ("we implement it
later"; the user is away from the machine). Until then: hosted CI with
the self-healing setup stays the gate; postflight debts keep accruing
and are all settled by the first green full run on the runner.

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
**Round-2f (2026-08-19): r2e red at corr 0.182 — the r2a trap resurfaced
in the plumbing.** r2e's test premise is FIXED-WHITE declared weighting
(pass-1 floor, default C0), but the call was switched to jointAcc while
the white cfg.override was dropped — the loader therefore declared
cfg.cov, THE TRUTH, and truth-whitening made Shat track the C0 prior
(alpha2 ~ 1): corr 0.182 on the runner. This also independently
CORROBORATES r2e's measurements: the failure mode is exactly the r2a
lesson, not the floor design. Fix: cfg.override (white diag) restored in
jointCfgColoured (jointAcc threads it as UserOpt); jointAccMapped kept,
annotated as the F-15 closed-loop reproduction fixture (not called by
the floor test). The USER VETO on the lowered 0.85 floor (spec said
0.95, measured ceiling min 0.921 open-loop) remains OPEN and is asked in
PR #4. **r2f VERDICT (same day): CI GREEN — 489/489, 0 failed,
test_ellmer 59/59 incl. the floor test (test_parallel skipped per
LSQVCE_CI; the corr value prints only on failure — the floor held).**

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

**ROUND NOTE (2026-08-19, WP-34 round 1, model Claude Fable 5; branch
`claude/wp-34`, riding the WP-33 DONE flip).** IMPLEMENTED in
adjustVCEjoint only: every enumerated arc gets the key
s<sys>_<obsType>_d<day>_a<arc> (redundancy's 'arc' keyFcn format,
verbatim); the factor store is a containers.Map PINNED on the PASS-1
ENUMERATION (skipped-on-pass-1 arcs included, default 1; Sigma2_0 aligns
to that order, out.arcKeys exposes it); lookups and write-backs are BY
KEY — a skip cannot shift later arcs; oneArc now RETURNS its skip reason
('gap'|'nonPD') and sweepArcs records {key, reason}; a key set that
changes between passes raises lsqvce:adjustVCEjoint:arcMismatch
(vanished-from-enumeration = identity error, covariance failure = kept
skip); skips accumulate in hist.skips and out.skips with the pass
number; verbose prints a per-pass skip summary. New outputs
.arcKeys/.skips; errors arcMismatch/badSigma2_0. TESTS (test_ellmer):
induced non-PD on pass 2 via a covModel 'custom' override closed over
model.iter (the realistic F-02 trigger — the covariance model updating
mid-run); assertions: pinned keys equal across runs, skip reported as
{2, key, 'nonPD'}, the SKIPPED arc's factor BYTE-IDENTICAL to its own
pass-1 value (a MaxIter-1 control run), others RelTol 0.15 against the
no-skip run; plus the arcMismatch error test (meta drops the last arc
from pass 2). DEVIATION, pinned with motivation: the spec's "other arcs
byte-identical to the no-skip run" is unattainable by construction —
arcVarianceFactors normalises every factor by gammaBar (eqs. 6.5.29-30),
a mean over the arcs PRESENT, so removing one column renormalises all
factors; byte identity would require decoupling the estimator. The
byte-level guarantee is therefore placed where it truly holds (the kept
factor of the skipped arc), the alignment guarantee on the keys, and the
tolerance comparison (0.15) bounds the gammaBar renormalisation. §6.1:
both tests ERROR on the pre-fix code (no arcKeys/skips outputs) and the
alignment contract did not exist.
**Round-1 CI red + fix (2026-08-19):** the -eye trigger fired during
ACCUMULATION — decorrFactor hard-errors when NO eigenvalue is positive
(notPD, max = -1), so the pass-2 accumulator died before any sweep.
Corrected trigger: INDEFINITE Sig (one -0.5 eigenvalue, max eig 1 > 0):
decorrFactor floors the negative direction with a warning and
accumulation proceeds (the robust production path), while the sweep's
strict chol fails and oneArc skips — the honest mid-run scenario. Lesson
for skip-inducing fixtures: the accumulator and the sweep gate covariances
DIFFERENTLY (floor-and-warn vs skip); pick the trigger for the stage
under test.

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

**ROUND NOTE (2026-08-19, WP-35 round 1, model Claude Fable 5; branch
`claude/wp-35`).** F-04 DECISION: FLOOR, never abs — abs folds zero-mean
noise at zero-power bins (E|x| = 2 E[x+]; T6 measured the 2x injection),
a floor admits only MinAmplitude, and genuinely positive bins are
identical under both; psdVCE line, comment, MinAmplitude help and stamp
now all say the same thing. F-05: sigma2Hat = gamma2/gammaBar^2 is
clamped at the new documented option ClampEps (default 1e-8, on the
NORMALISED scale — healthy factors are O(1)); warn id
lsqvce:arcVarianceFactors:clamped lists count + first 10 arc indices;
new output .clamped (also on the badMean early return, all-false);
healthy inputs byte-identical (the clamp branch only executes on
any(clamped)). TESTS (test_ellmer, §6.1 each): zero-power-bin pin —
Shat(2) == MinAmplitude exactly where the old abs left 0.5*Sprior(2) =
0.202 (fixture signs verified: Sprior = [0.877 0.404 0.506 0.286], all
positive, deterministic); constructed negative alpha2 — .clamped flag,
exact ClampEps value, no negative factor leaves (old code: negative
escapes, no field); healthy no-op — .clamped all-false and sigma2
matches gamma2/gammaBar^2 within 4 eps (plain mean, closed formula).
**Round-1 CI red + fixes (same day):** (a) the zero-power fixture fed
per-LAG ratios, but psdVCE forms alpha2 IN THE FREQUENCY DOMAIN (Omega
and s are DCT-I transformed first, 6.5.25-26) — the fixture now picks
the target frequency-domain pair and inverse-transforms (ortho DCT-I is
self-inverse); (b) the no-op test's byte-identity claim was too strong:
the function's internal pooling differs from the closed formula by
summation order (CI-measured 2.6e-16 rel) — the honest regression is
mask all-false + 4-eps formula match, and the true byte guarantee is
that the clamp branch never executes on healthy input.

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

**ROUND NOTE (2026-08-19, WP-36 round 1, model Claude Fable 5; branch
`claude/wp-36`, riding the WP-35 DONE flip + the runner
approved-and-deferred decision).** THE ONE GENERATOR:
`+lsqvce/private/probeMatrix.m` — seeded: RandStream('twister', seed),
bit-identical across callers; EMPTY seed: a concrete seed is drawn ONCE
from `+lsqvce/private/probeStream.m` (persistent 'threefry4x64_20',
'shuffle' at first use — per-session, never the global stream) and
RETURNED as seedUsed so callers can pin it. Rademacher recipe unified to
sign(randn), sign(0)->+1 — the recipe TWO of the three callers already
used, so their seeded realisations are unchanged; traceEstimate's seeded
realisations CHANGE (it used 2*randi-1; unbiasedness unaffected —
documented in its help). REWIRED: sharedProbe (local drawProbes DELETED,
rule 8), streamContribution (the false "private stream" comment replaced
by the now-true one), traceEstimate (the global-stream default is gone).
VCE: an empty TraceSeed under 'hutchinson' is pinned ONCE up front via
probeMatrix's seed draw, so every loop pass AND the final vc.r pass probe
with IDENTICAL Z (vc.r == the converged loop redundancy); new output
vc.traceSeedUsed (replay with it reproduces vc bit-identically). TESTS
(no new suite files): test_traceEstimate — global-RNG invariance with
empty seeds (rng state bitwise before==after; FAILS pre-fix) +
same-seed-identical-probes across sharedProbe/traceEstimate (via the
'Probes' bypass, 4-eps); test_contribution — streamContribution's
invariance half (needs the accumulator fixture) + the vce
traceSeedUsed replay test (vc.r and sigma bitwise; FAILS pre-fix: field
absent, final pass drew fresh probes).
**Round-1 red + fix (2026-08-19):** one failure, MY test's API mistake —
traceEstimate returns the SCALAR estimate, the draft dot-indexed
`.trace` on a double (MATLAB:structRefFromNonStruct). Everything the WP
actually changed ran green on the first CI pass (497/498: invariance,
cross-caller identity via Zbar, stream invariance, and the vce bitwise
replay all held).

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

**ROUND NOTE (2026-08-19, WP-37 round 1, model Claude Fable 5; branch
`claude/wp-37`, riding the WP-36 DONE flip).** KERNEL: only the SUM over
each Sigma^-1 sub-diagonal enters s_n, so the dense X = W\eye(N) is gone
from the main path. TOEPLITZ arcs (every stationary covModel: acf, arma,
psd) are DETECTED from three probe columns of Sigma (cols 1, 2, mid via
W'*(W e_j), tol 1e-10 relative) and solved by Levinson-Durbin
(toeplitzInvFirstCol: Durbin + the Yule-Walker identity x = [1;a]/(E c1))
+ Gohberg-Semencul: Sigma^-1 = (1/x1)(L(x)L(x)' - L(z)L(z)'), whose
triangular-Toeplitz sub-diagonal sums are O(N) weighted correlations
S_n(u) = sum_j u(j+n)u(j)(N-n-j+1) — O(N^2) time, O(N) memory. General W
falls back to streaming the columns of W^-T with an nMaxLag ring buffer:
still O(N^3) time (no structure, no shortcut — stated honestly in the
help) but O(N*nLag) memory, which finally makes the help's memory claim
TRUE. WantRbar alone still materialises (testing path, documented). Both
recursions were mirrored EXACTLY in Python first and validated against
dense inverses to 1e-14 (this note is the record; five Toeplitz sizes +
ring fallback). HOT-LOOP TOLLS: decorrFactor computes info (incl. the
per-batch condest) only when nargout > 1; redundancy's task enumeration
caches the meta calls and preallocates (the old tasks(end+1) copied
O(nTask^2)). TESTS (test_ellmer): dense-reference agreement at 1e-12 on
the Toeplitz path (probes AND locals engaged) and on a deliberately
non-stationary covariance (fallback path); the §6.1 SPEED lever at
N = 4000 with a 5 s budget the dense kernel blows by an order of
magnitude. NO API change.

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

**ROUND NOTE (2026-08-19, WP-38 round 1, model Claude Fable 5; branch
`claude/wp-38`, riding the WP-40 DONE flip + the DECISION RECORD for the
three accepted r2 decisions + the WP-42 backlog entry).** ONE
`+lsqvce/private/makeKeyFcn.m` (the redundancy wording verbatim, presets
arc/obsType/system/day + handle pass-through) replaces the two verbatim
copies in redundancy and streamContribution; adjustVCE.vceGroupKeyFcn is
a DIFFERENT class-scoped scheme and deliberately NOT unified (documented
in the new help). ONE `+lsqvce/private/arcBatchWhiten.m` owns the
batch/overlap/Wfcn-Cfcn skeleton; the two callers keep their numeric
payloads and accumulation order as NESTED callbacks — byte-identical
outputs by construction (the shared Cfcn branch also whitens lb, which
the probe callback simply ignores: outputs unchanged, documented). AUDIT
LOCK: audit_tree.py check 10b fails on any second definition of either
name or a definition outside +lsqvce/private/ — red-tested by planting a
copy (exit 1, clean FINDING line) and restoring (exit 0). **PINNED
MINI-DEVIATION (user approves by merge): the GroupBy error id unifies to
lsqvce:makeKeyFcn:badGroupBy** (was lsqvce:redundancy:/
lsqvce:streamContribution:badGroupBy; no test asserted the old ids; both
help ERRORS sections updated). No public surface change; the suite's
byte-identity is carried by the standing fixtures (§6.1 lever = the
audit lock itself, which FAILS on the pre-WP-38 tree).

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

**ROUND NOTE (2026-08-19, WP-39 stage (a) round 1 = a1, model Claude
Fable 5; branch `claude/wp-39`, riding the WP-41 DONE flip).** MIGRATED
to typed `arguments` signatures with EXACT UpperCamelCase names:
(1) adjustVCEjoint — names were already UpperCamel, UNCHANGED; validators
became mustBe*; Norm/MeanType are now string-typed and flow into
psdVCE/dctI's inputParser unchanged (inputParser accepts strings).
(2) adjustBlock — Backend/Reg/RegLambda/RegMatrix/Parallel/Verbose/
ShowSeconds/ExportDir (was lowercase); enum VALUES case-sensitive
lowercase via mustBeMember (strcmpi used to accept any casing); the v1
"ShowSeconds < 0 = show and leave open" contract PRESERVED (an initial
mustBeNonnegative draft would have broken it — caught pre-commit, help
states the v1 contract). Zero call sites needed sweeping.
(3) iterateNoiseModel — UpdateTypes/MaxIter/SigmaTol/WhiteFactor/
WhiteLags/MaxLag/OrderGrid/AcvTol/Detrend/Verbose; the data-dependent
UpdateTypes default (all fields of cov0) kept via the {} sentinel; both
existing call sites already used UpperCamel, no sweep.
TESTS (one per function, the fails-on-old lever = a LEGACY-lowercase
name that the old case-insensitive inputParser ACCEPTED and the
arguments block must reject; ids checked loosely via the 'MATLAB:'
prefix for release robustness): test_ellmer (joint: 'numlags', wrong
type, enum casing), test_adjustBlock ('backend', 'Backend'="QR",
RegLambda<0), test_lsqvce_ext ('maxiter', MaxLag<8).
DOCS: three HTML pages regenerated from the corrected help
(delete+rebuild); the LIVE manuals edited BY HAND because
build_manuals.py is _STALE_GUARDed (rebuilding would revert shipped
fixes, HISTORY 1.4C): user manual gets the Name-Value-casing callout AND
its fluent example switched to the canonical UpperCamel names (mergeNV
is case-insensitive today, so the old lowercase example worked but
breaks at a2 — fixed ahead), reference manual gets the RELEASE NOTE
listing the exact old->new name breaks.
**RELEASE NOTE (verbatim, for the repository notes):** "WP-39 (staged):
lsqvce.adjustVCEjoint, lsqvce.adjustBlock and lsqvce.iterateNoiseModel
now use typed `arguments` signatures. Option names are EXACT
UpperCamelCase (adjustBlock: backend/reg/regLambda/regMatrix/parallel/
verbose/showSeconds/exportDir -> Backend/Reg/RegLambda/RegMatrix/
Parallel/Verbose/ShowSeconds/ExportDir; iterateNoiseModel: updateTypes/
maxIter/sigmaTol/whiteFactor/whiteLags/maxLag/orderGrid/acvTol/detrend/
verbose -> UpperCamel equivalents). Previously-tolerated sloppy casings
and enum-value casings (e.g. 'QR') now error."
**DEFERRED to stage (a) round 2 (a2), deliberately:** lsqvce.adjust
(64 lowercase call sites to sweep) and the adjustVCE fluent
configurators (37 sites; mergeNV is case-insensitive and canonicalises
to the stored UpperCamel field names — the a2 migration replaces mergeNV
semantics with exact-name arguments and sweeps the remaining lowercase
sites). Stages (b) VCE/diagnostics layer and (c) primitives follow.
**Round a1-2 (2026-08-19, same day): PREMISE CORRECTED by CI
measurement.** The a1 red (test_adjustBlock: lowercase 'backend' not
rejected) exposed that MATLAB matches comma-syntax Name-Value NAMES
case-insensitively (and partially) EVEN under `arguments` blocks -- the
"EXACT-name matching" premise of the F-07 casing break was WRONG, and
the two "passing" lowercase tests were FALSE POSITIVES (the
dummy-driven body crashed with a MATLAB:-prefixed id that satisfied the
loose check; lesson pinned: a loose id-prefix check must not share its
prefix with plausible body errors -- pick levers whose old/new
behaviour differs at the SIGNATURE). CORRECTED across help, both live
manuals, master_content_b and the tests: legacy name casings KEEP
WORKING (good for users); the real breaks are (1) enum VALUE casing
(mustBeMember vs strcmpi -- the fails-on-old lever: 'Decay'/'QR'),
(2) typed constraints firing at the signature, (3) unknown option names
rejected instead of tolerated. **RELEASE NOTE (corrected, supersedes
the a1 text):** WP-39 (staged): lsqvce.adjustVCEjoint,
lsqvce.adjustBlock and lsqvce.iterateNoiseModel now use typed
`arguments` signatures with canonical UpperCamelCase option names.
Legacy name casings continue to match (MATLAB name-value semantics).
Breaking: enum values are case-sensitive lowercase ('qr', not 'QR';
'decay', not 'Decay'); wrongly-typed or out-of-range option values and
unknown option names now error at the signature.
**ROUND NOTE (2026-08-19, WP-39 stage (b) = b1, model Claude Fable 5;
branch `claude/wp-39` continuing past the PR #14 merge).** MIGRATED to
typed `arguments` signatures, canonical UpperCamelCase, lowercase
o-struct mappings keeping bodies untouched: vce, varComp (typed
positionals; CheckD left untyped so both [] and the diagnose struct
pass), redundancy (Parallel default "" via char()), redundancyLags
(names were UpperCamel; ebar/nMaxLag stay REQUIRED positionals — an
optional draft would have changed the contract, caught pre-commit),
streamContribution, contribution, psdVCE, sharedProbe (lowercase ->
NumProbes/Seed/Probes/Distribution; 18 legacy call sites keep working
via MATLAB's case-insensitive NV matching — no churn sweep, the a2
decision). PLACEMENT LESSON (caught by the page builder's "no help
block" skip): an `arguments` block inserted BETWEEN the function line
and the help comment silently kills `help`/`doc` for the function —
MATLAB help must directly follow the signature; all eight blocks moved
AFTER the help, and the a1 trio was verified unaffected. INTERNAL
FORWARDING verified: vce->redundancy lowercase pairs ride the same
case-insensitive matching; string enums flow into the remaining
inputParser consumers (dctI) unchanged. TESTS (a1-2 lesson): the
fails-on-old levers assert the SPECIFIC MATLAB:validators:mustBeMember
id, unreachable for the old strcmpi/inputParser code — vce 'Exact',
varComp 'Rigorous', redundancy 'Hutchinson'/'PARFOR', psdVCE 'Lowest',
sharedProbe 'Rademacher'; plus typed-constraint checks (Probes 0, Level
3, Gamma 1.5, nMaxLag 0, non-integer u). DOCS: help OPTIONS casing swept
(help lines only, scripted token map), casing note + stamp per file,
eight HTML pages regenerated, reference-manual release note extended.
REMAINING: stage (b) round 2 = the diagnostics family (diagnose,
residuals, report, outliers, robust, ftest, crossCorr, ellipse,
selectArOrder), then stage (c) primitives.
**b1 round 2 (2026-08-19): field-completeness repair.** CI red (17
errors, all one cause): varComp's b1 arguments block DROPPED three of
the twelve old options (backend/seed/reprobe) — the body read o.backend
into a struct that never carried it. Fixed by restoring Backend/Seed/
Reprobe (typed, canonical) and extending the o-struct mapping. LESSON
PINNED, mechanical and reusable: after any parser->arguments migration,
diff the set of o.<field> READS in the body against the MAPPED fields
(or the declared opt.* names when o = opt) — a five-line Python check
now run for all eight b1 functions: ALL COMPLETE. The old parse block is
the option inventory; migrate against IT, not against the help text.
**b1 round 3 (2026-08-19): the last red = a documented-ERRORS contract
moved.** redundancyLags' nMaxLag <= 0 now fails mustBePositive AT THE
SIGNATURE where the body used to throw lsqvce:redundancyLags:badNMaxLag
(the contract test caught it -- working as designed). RULE adopted and
release-noted: scalar per-argument constraints live in the typed
signature and surface MATLAB validator ids; CROSS-ARGUMENT shape checks
stay in the body with their lsqvce: ids (badW/badUiZ unchanged and
green). Test, ERRORS help and the reference-manual release note updated
accordingly.
**ROUND NOTE b2 (2026-08-19): the remaining VCE helpers + the full
diagnostics layer.** MIGRATED to typed `arguments` (canonical
UpperCamel, lowercase o-struct mappings keep bodies untouched, the
field-completeness check green for all): arcVarianceFactors, poolGroups
(Strategy members are mergeLow/NEAREST -- a guessed mergeAll in the
first spec was caught by the post-hoc grep, fixed pre-commit),
vceHutchinson, traceEstimate, jointModelOpt, diagnose, snoop, outliers,
residuals, ftest, ellipse, crossCorr, systematic, selectArOrder, and
report (whose SECOND option 'width' the early survey missed -- the
verbatim parse block was the inventory, per the b1 lesson).
**PINNED DEVIATION (user approves by merge): robust follows the adjust
pattern** -- inputParser stays (KeepUnmatched forwards unmatched options
to lsqvce.diagnose), canonical UpperCamelCase names + typed validators,
enum VALUES case-sensitive via strcmp, and an r_-remap keeps the body
untouched. TESTS: one grouped b2 signature test (test_lsqvce_ext) with
per-function enum-VALUE levers asserting mustBeMember's specific id
(traceEstimate Distribution, poolGroups Strategy, ellipse Scale,
jointModelOpt PsdMode) and robust's case-sensitive Influence through
the parser's own ArgumentFailedValidation id. DOCS: 16 stamps, 16 pages
regenerated, the reference-manual release note extended. Stage (c)
(infrastructure/classes: BlockAccumulator, Combine, makeLoader,
adjustVCE ctor, Accumulator, mapreduce, datum, solveNonlinear,
toeplitzCofactor, dctI, BatchIndex, PrefetchLoader, Results) remains.
**ROUND NOTE c2 (2026-08-19) -- THE R2-CLOSING ROUND.** The seven class
files: Accumulator + Combine already carried arguments blocks (arrived
in earlier merged commits). MIGRATED: BatchIndex ('By' enum
sensor/day/both; the old lower() call dropped -- the enum is
case-sensitive lowercase now); PrefetchLoader ('Enable' keeps its v1
"auto"|true|false contract -- UNTYPED in the block since a mustBeMember
would kill the logical form, validated in-body with the NEW id
lsqvce:prefetch:badEnable); Results.convergence (Plot/Stream);
BlockAccumulator ctor (nc typed positional; 16 opts, names unchanged
UpperCamel; Backend/Reg/Parallel via mustBeMember; Verbose UNTYPED =
the historical logical-or-level contract; the r remap keeps the ctor
body untouched incl. the CondWarn<=CondAbort cross-check with its
lsqvce:accum:badThresholds id); adjustVCE ctor (GroupBy untyped --
accepts handles; VceMode enum; LogFid/LogFlag empty-sentinels; o remap
replaces o = p.Results FIRST OCCURRENCE ONLY -- two identical lines
exist for other methods, the global-replace assert caught it) and the
three option-taking methods plots / empiricalStochastic / validate,
whose arguments DEFAULTS reference the live object state
(obj.plotOpt.*, obj.shapeOpt.WhiteAlpha, verbosity.progress(obj.verbose)
-- legal: defaults may reference earlier-declared arguments); validate's
narginchk replaced by the block (missing R now = MATLAB:minrhs) while
the isa contract check stays in the body with its documented lsqvce: id
(the b1-r3 rule). TESTS: grouped c2 signature test with the
mustBeMember levers (BlockAccumulator Backend, adjustVCE VceMode), a
typed-constraint check, and BOTH legs of the Enable contract. DOCS:
five pages regenerated, help casing swept + stamps, the release note
declares the WP-39 rollout COMPLETE.
**c2 round 2 (2026-08-19): two standing contract tests moved with the
migration** (the b1-r3 rule's mirror, both caught by CI as designed):
test_ellmer's VceMode-rejection now expects mustBeMember (was the
parser's ArgumentFailedValidation), and test_reliability's
bad-construction test moves ALL FOUR preempted expectations --
nc<=0 twice (expectedPositive -> validators:mustBePositive; the old
validateattributes id) and the two bogus-enum checks
(ArgumentFailedValidation -> validators:mustBeMember). The
CondWarn<=CondAbort cross-check keeps its lsqvce:accum:badThresholds
body id, green. Also this round: the ctor's dead narginchk(2,inf)
deleted (the arguments block enforces arity; missing positionals =
MATLAB:minrhs; rule 8). **With c2 green, WP-39 closes --
and with it the R2 track (F-01..F-15, every finding addressed or
formally dispositioned).** Remaining after R2: the self-hosted runner
migration (settles ALL postflight debts on its first green full run)
and the WP-42 backlog.
**ROUND NOTE c1 (2026-08-19): the infrastructure functions.** MIGRATED
to typed `arguments` (canonical UpperCamel, o-struct mappings, bodies
untouched): makeLoader (15 options; the body's long->short o aliasing
at its top is FED by the mapping -- the completeness checker's
single-struct regex flags it as a false positive, noted), mapreduce,
datum, toeplitzCofactor, dctI (the hot primitive: internal callers
already pass canonical 'Norm'/'Inverse'; enum members ortho/ellmer now
case-sensitive). **PINNED DEVIATION (user approves by merge):
solveNonlinear follows the adjust pattern** (KeepUnmatched forwards to
the solve/diagnostics chain; canonical names + typed validators; r_
remap keeps the body untouched). Standing contract tests checked and
unaffected (makeLoader badDay is a closure-body check; dctI tooShort is
a positional body check). Grouped c1 signature test with enum-VALUE
levers (dctI Norm, toeplitzCofactor Form, mapreduce Mode via
mustBeMember's id; solveNonlinear Method via the parser's id). Six
stamps, six pages regenerated, release note extended. REMAINING = c2,
the class surface: BlockAccumulator (ctor options already canonical),
Combine (FOUR parsers across methods), Accumulator, adjustVCE ctor
(21 options), BatchIndex, PrefetchLoader, Results -- then WP-39 and the
R2 track close.
**c1 round 2 (2026-08-19):** two reds. (1) makeLoader's PROTECTIVE
contract test (removed terse aliases must error) moved to the
signature-level id — arguments blocks reject unknown NV names with
MATLAB:TooManyInputs, not the parser's UnmatchedParameter; the b1-r3
rule extended: unknown-name rejection ids also live at the signature
now. (2) An arity bug in MY c1 test: solveNonlinear takes FOUR
positionals; the draft passed two, 'Method' landed as nBatch and the
parser never saw the pair — fixed with the full positional set. Lesson:
signature tests must satisfy the POSITIONAL contract before probing the
NV surface.

**Round a2 (2026-08-19): adjust + the fluent configurators — two PINNED
DEVIATIONS (user approves by merge), both consequences of the corrected
a1-2 premise.** (1) `adjust` KEEPS inputParser BY DESIGN: its
KeepUnmatched is a documented FEATURE (unmatched options forward to
LSQVCE.DIAGNOSE via stripArgs — strcmpi, so any casing strips), which an
arguments block cannot express. The F-07 substance lands anyway:
canonical UpperCamelCase names on its five own options (inputParser
matches case-insensitively — the 64 lowercase call sites keep working
UNSWEPT, and under the corrected premise a sweep buys nothing: behaviour
is identical to the arguments-migrated trio), TYPED validators on all
five, and the case-sensitive lowercase enum value on 'Method' ('QR' now
errors). The help states the pass-through contract explicitly; forwarded
diagnose options keep diagnose's lowercase surface until its own stage.
(2) The fluent configurators KEEP mergeNV: it already canonicalises
case-insensitively to the stored UpperCamelCase fields — exactly
MATLAB's NV semantics — and an arguments migration would need
observable-change-free sentinel machinery; typed validation stays
downstream. TEST (test_lsqvce): the a1-2 lesson applied — the two
fails-on-old levers assert the SPECIFIC
MATLAB:InputParser:ArgumentFailedValidation id (the old code had no
validators, so that id was unreachable), plus an end-to-end
pass-through+lowercase regression. The 64-site sweep is RETIRED as
unnecessary (premise correction), not deferred.

**ROUND NOTE (2026-08-19, WP-39 stage (a) round a1, model Claude Fable 5;
branch `claude/wp-39`, riding the WP-41 DONE flip). Staging WITHIN (a),
honest to the sweep load:** a1 = the three front doors with 0-2 call
sites — adjustVCEjoint (names already UpperCamel, parser→arguments
only), adjustBlock (Backend/Reg/RegLambda/RegMatrix/Parallel/Verbose/
ShowSeconds/ExportDir; enum VALUES now case-sensitive lowercase via
mustBeMember where strcmpi accepted any casing; 0 name-value call
sites), iterateNoiseModel (UpdateTypes/MaxIter/SigmaTol/WhiteFactor/
WhiteLags/MaxLag/OrderGrid/AcvTol/Detrend/Verbose; both call sites were
already UpperCamel). a2 = adjust (64 lowercase call sites!) + the
adjustVCE fluent configurators (37 sites, mergeNV) — next round. LESSON
(pinned, cost one audit red): the `arguments` block goes AFTER the help
comment block — placed before it, MATLAB help/doc go blind and the page
builder finds nothing. UpdateTypes' data-dependent default
(fieldnames(cov0)) is kept via the documented {} sentinel. SIGNATURE
TESTS (fail on the old code): the old lowercase names now ERROR (the
release-break lever), unknown names error, enum-casing and type
violations error — loose id checks (try/catch) for release robustness.
DOCS per the user's instruction: help OPTIONS sections renamed + stamps;
three HTML pages regenerated; the USER MANUAL gains an "Option names and
typed signatures" paragraph — edited in the LIVE html (build_manuals is
1.4C-stale-guarded and would revert shipped fixes) AND authored in
master_content_b.py for the eventual back-port.
**RELEASE NOTE (collect for the version bump):** WP-39 migrates the
public surface to typed `arguments` blocks. Name-Value names are
UpperCamelCase and EXACT-match; previously-tolerated sloppy casings
('backend', 'maxiter', …) now error. Enum values are case-sensitive
lowercase. Unknown option names error instead of being ignored.
Migrated in a1: adjustBlock, adjustVCEjoint, iterateNoiseModel.

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

**ROUND NOTE (2026-08-19, WP-40 round 1 = F-10…F-13, model Claude Fable 5;
branch `claude/wp-40`, PR #9; rides the WP-37 DONE flip + the ORDER
AMENDMENT (user option A) + the binding read-the-order-line process
rule). F-15 is round 2 by design (WP-33 two-round pattern; repro fixtures
jointAccMapped + validation/wp33_loop_mirror.py stand ready).**
F-10 (vce): the silent rg = max(rg,1) is replaced by varComp-style
estimability — new documented option 'freezeR' (default 1; deliberately
distinct from 'minRedundancy', redundancy's per-observation testability
knob passed through): groups with r_g < freezeR keep their sigma for the
pass (upd = 1), warned ONCE per run (lsqvce:vce:lowRedundancy, keys +
r_g); healthy groups bitwise unchanged. F-11 (varComp): a negative
rigorous component FREEZES its group (s = 1) instead of s = 1e-6 — the
old clamp multiplied the CUMULATIVE scale by 1e-6 and fed a 1e6
reweighting oscillation; warning + help motivate the choice. F-12 gate
half (adjustVCE.run): joint convergence additionally requires
hist.jointMaxDev(end) <= outer Tol at BOTH gate sites — shapeStable
degenerated from pass 2; the claims half stays with WP-41. F-13
(iterateNoiseModel): maxAcf pools PER ARC (per-arc mean removal, pooled
numerators/denominator, no boundary products) via the new
gatherResidualArcs; dead gatherResiduals DELETED (rule 8). TESTS: freezeR
= Inf freeze+warn with sigma bitwise at start (pre-fix: option unknown);
negative-component freeze on a deterministic seed-searched two-group rig,
scale stays O(1) (pre-fix: 1e-6); joint-gate NOT-converged while
jointMaxDev(end) > Tol (pre-fix: converged); F-13 serial smoke on a new
minimal fixture. **PINNED DEVIATION (user approves by merge): the F-13
fails-on-old DISCRIMINATION needs arc-mean-structured whitened residuals
that only runner-tier fixtures produce — deferred to test_parallel on the
self-hosted runner; CI carries the smoke + healthy-path regression.**
**Round-1 red + fix (2026-08-19, same day):** the F-11 rig's group 2 had
m2 = 6 observations — BELOW varComp's minObs default of 10, so it was
excluded from estimation and the negative branch (bad = s<=0 &
estimable) could never fire: 0/25 seeds. Rig fixed to m2 = 12 with true
sigma^2 = 0.03; Python frequency check: negative in 14/40 seeds (first
firing seed 2), so the deterministic seed search over 1:40 is safe.
Lesson: an estimability-gated pathology test must first PASS the
estimability gate.

**ROUND NOTE (2026-08-19, WP-40 round 2 = F-15, model Claude Fable 5;
branch `claude/wp-40-r2`; rides the WP-34 DONE flip). MEASUREMENT SERIES
(all in validation/wp33_loop_mirror.py, permanently extended; 8 seeds per
cell, corr over bins 2:end):** (1) The mirror still modelled the
PRE-WP-35 abs update — re-measured under the shipped FLOOR semantics: the
instability PERSISTS and some seeds WORSEN (8x96 seed 8: 0.94 -> -0.21;
16x128 seed 5: 0.63 -> -0.64) — the abs suspect is dead. (2) ZeroFreqRule
off and nLag 10->16->24 only SHIFT which seeds fail — both suspects
exonerated. (3) Constant damping (gamma 0.3..0.7) kills the collapses but
lowers the ceiling (under-converged). (4) THE CURE: Robbins-Monro decay
gamma_it = 2/(2+it) — full early steps, averaging later — lifts the
worst-seed corr to 0.87-0.91 across ALL grids; Polyak tail averaging of
the returned spectrum adds ~0.02 (and does NOTHING without the decay:
averaging cannot rescue an undamped collapse). SHIPPED, all defaults
bitwise-old: psdVCE 'Gamma' (the primitive, alpha2.^Gamma); adjustVCEjoint
'StepSchedule' ('none'|'decay') computing Gamma_it and threading it, and
'AvgTail' (log-domain Polyak average applied ONLY to the returned
Shat/C, via psdVCE's exact dctI inverse; the loop iterates on raw C).
TESTS: the F-15 repro fixture (jointAccMapped) with decay+tail4 must
clear corr > 0.8 (mirror-measured min ~0.89 leaves probe margin; fails
pre-fix: options unknown); bitwise default regression (one pass, explicit
defaults == omitted).
**Round-2 CI red #2 (2026-08-19, this session):** still
psd2acv:expectedNonnegative, contiguous stack pinned: psd2acv:61 <-
covModel:409 (the 'toeplitz' psd2acv) <- covModel:213 <- synthLoader:152,
i.e. FIVE LINES BELOW the elementwise max() floor -- statically
impossible under my reading of max() semantics, so the floor is replaced
by indexed sanitation (S(~isfinite(S) | S < fl*Smax) = fl*Smax) plus a
diagnostic assert (lsqvce:covModel:psdSanitationEscape names min(S),
#nonfinite, fl, Smax) so any residual escape reports its value class
instead of failing downstream. LATENT BUG found during the chase and
fixed: synthLoader:224, the GENERATION path's psd2acv, was the toolbox's
one call with RAW cm.S -- now floored identically.
**Round-2 CI red #3 -> ROOT CAUSE (2026-08-19): the escaping S was
COMPLEX.** The sanitation assert passed (real-part comparisons: complex
values satisfy S > 0 elementwise on their real parts) while
validateattributes 'nonnegative' rejects complex -- and the source is the
new damped update itself: alpha2 can be negative at near-zero-power bins
(the WP-35 finding), and a FRACTIONAL power of a negative ratio is
complex; alpha2.^Gamma propagated complex Shat -> C -> the declared
spectrum. Fix in the primitive: the WP-35 floor philosophy applies BEFORE
the power -- Gamma < 1 uses max(alpha2,0).^Gamma; Gamma == 1 keeps the
pre-existing expression VERBATIM (bitwise default regression). The
indexed sanitation + diagnostic assert and the synthLoader generation
floor stay (both are real hardening regardless).
**Round-4 (2026-08-19): the round-3 clamp was WRONG-VALUED.** max(alpha2,0)
under Gamma < 1 creates ABSORBING ZERO bins (Shat_j = 0 -> Sprior_j = 0
forever): the crash chain died but corr fell to 0.511. The mirror -- the
instrument the 0.87-0.91 promises were measured on -- HOLDS at
nonpositive ratios (a2 -> 1: an unidentifiable ratio takes no step this
pass; wp33_loop_mirror.py, 'nonpositive: hold'). psdVCE now ports that
semantics exactly; Gamma == 1 stays verbatim. Lesson pinned: a mitigation
port must replicate the MEASUREMENT INSTRUMENT's semantics, not a
plausible variant of them. F-15 §8 disposition updated: mechanism identified,
mitigation shipped opt-in, severity call = user by merge.
**Round-5 DIAGNOSIS SESSION (2026-08-19, measurements only — no code
push; pinned so no session repeats them).** Round-4 red: corr 0.734 vs
bar 0.8 (the mirror had measured 0.87-0.91). Three suspects EXONERATED
by extending the mirror: (1) PROBE NOISE — shared-probe Hutchinson at
Z = 300 is indistinguishable from exact s_n in the mirror (the common
space is only nc = 3 dimensional); 12 runs identical to three decimals.
(2) ARC-FACTOR COUPLING — the driver's s2*C weighting of the pooled
omegaT/sT columns was mirrored (incl. the gammaBar trimmed-robust
renormalisation): NO effect on this homogeneous-arc fixture (factors sit
at ~1). (3) SEED TAIL — 40 data seeds at 12 iters + tail-4: min 0.877,
p5 0.905, P(corr < 0.8) = 0.00 — MATLAB's 0.734 lies OUTSIDE the mirror
distribution, so a genuine mirror<->MATLAB implementation difference
remains. PRIME REMAINING SUSPECT: the DECLARATION path — the mirror
zero-extends a 64-bin dctI acv beyond lag 64, while MATLAB goes covModel
'psd'/psdMode 'toeplitz' through psd2acv (full m-lag cosine transform +
its own floor): a different declared Sigma means a different whitening
trajectory. NEXT STEP: make the mirror's declaration psd2acv-faithful
and re-measure; if the gap closes, the cure is spectral-declaration
resolution/smoothing, not more damping. SIDE FINDING (pinned): LONGER
schedules are WORSE — 20 iters + tail-8 over the same 40 seeds: min
0.773, P(<0.8) = 0.03; the decayed late steps drift and the longer tail
averages over the drift. Keep 12/4 until understood.
**Round-6 (2026-08-19): the mirror made FAITHFUL — the mitigation's
promises fall, the honest claims ship (code commit 8f383f7; this note is
its doc follow-up — the note's anchor missed the disposition wording in
the code commit).** The declaration path was the gap: the mirror
zero-extended a 64-bin dctI acv, MATLAB integrates the full m-lag
psd2acv cosine transform. With the faithful kernel the mirror REPRODUCES
MATLAB's regime (8x96 stabilised: min 0.668, P(corr<0.8) = 0.31 — the
CI 0.734 sits inside), and the r2e "stable volume" 16x256 collapses to
P(<0.8) = 0.75 while OPEN-LOOP pass 1 scores min 0.96/med 0.97 there.
Further falsified as cures (16-seed grids each): lag windows on the
declared acv (Lc 10/20/30), pchip-vs-linear interpolation, log-shrinkage
of the declaration toward white (beta 0.7/0.5). SHIPPED in 8f383f7: the
stabilised test asserts the honest property (damage-limitation floor 0.5
with the undamped ~0.18 contrast in the message; the 0.8 bar rested on
the unfaithful instrument); adjustVCEjoint's StepSchedule/AvgTail help
claims corrected BEFORE they become F-14-class doc drift;
wp33_loop_mirror.py docstring supersedes its old grids. USER DECISIONS
at PR #10: (a) accept the honest bar + production guidance, (b) F-15
severity under the new picture, (c) whether deeper method work
(monotone/trust-region alternation) gets its own WP.

**DECISION RECORD (2026-08-19, USER: "All open decisions accepted" —
detailed here for later pickup, binding):**
**(a) Honest bar + production guidance — ACCEPTED, now BINDING policy.**
The closed-loop spectral VCE through the psd2acv declaration path is
operated OPEN LOOP or with AT MOST 3 damped passes
('StepSchedule','decay' + 'AvgTail' <= 4) in any production
configuration; the stabilised-loop test contract is the
DAMAGE-LIMITATION floor corr > 0.5 (undamped contrast ~0.18) and no
session may raise it back toward an accuracy bar without new
faithful-mirror evidence. The corrected StepSchedule/AvgTail help claims
are the authoritative wording (WP-41's doc sweep must NOT "improve" them
back).
**(b) F-15 severity — ACCEPTED as UNDERSTOOD + GUIDED; F-15 CLOSES with
WP-40.** The pathology is real (bistable/counterproductive closed loop),
the mechanism is measured (declaration-path feedback; six exonerated
suspects listed in the §8 disposition + round notes), the mitigation
limits damage, and the pinned guidance avoids the regime. §8 keeps the
full record; no open finding remains under F-15.
**(c) Deeper method work — ACCEPTED as a BACKLOG entry, NOT scheduled:**
**WP-42 (BACKLOG, user-approved 2026-08-19 for later pickup):
monotone/trust-region alternation.** Scope sketch: replace the
multiplicative alpha2 step by a monotone-descent or trust-region update
of the declared spectrum (candidates: line search on the whitened-domain
negative log-likelihood; a KL/Itakura-Saito trust region per pass);
INSTRUMENT: the psd2acv-FAITHFUL mirror kernel (wp33_loop_mirror.py, r6)
is the mandatory pre-port evidence tool; ENTRY CRITERIA: a production
need for closed-loop operation (e.g. arcs too heterogeneous for one
open-loop pass) or a research interest in beating the open-loop 0.96
reference; SUCCESS CRITERIA: faithful-mirror P(corr < open-loop-quality)
materially reduced across the 8x96 and 16x256 grids, THEN the MATLAB
port. Until picked up, the row lives here — deliberately NOT in the R2
ledger (the R2 track closes with 38/41/39).

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

**ROUND NOTE (2026-08-19, WP-41 round 1, model Claude Fable 5; branch
`claude/wp-41`, riding the WP-38 DONE flip). Docs only — the audit
enforces zero logic drift.** (1) JOINT CLAIMS (F-12 half): adjustVCE's
'joint' help now says what the mode does — pooled per-lag sums + stream
scales per pass, NO per-arc factors — and points to ADJUSTVCEJOINT for
the full ten-step co-estimation. (2) r_g CLAIMS: every
"r_g = tr(N^-1 N_g^cc)" assertion is gone; contribution's .r field is
documented as the cc-attribution trace it is (historical field name
kept), the Level-1 share renamed to the trace share with the one-line
Level-2 contrast; adjustVCE's .contribution doc likewise. (3) RNG
sentences already carried WP-36's contract (checked, no drift found
beyond WP-36's own edits). (4) redundancyLags' 'Solver'/struct-W ghost
removed with an honest tombstone. (5) vce documents accFcn([]) = first
pass, unit weights. (6) STUB DECISION (user approves by merge): the two
removed-name redirect pages (lsqvce.recover, lsqvce.SeparableProblem)
are KEPT — they already carry the "(removed)" title, a visible banner
line and a 0-s refresh; audit check 4 now FORMALLY enforces the rule
(only these two, only with title+refresh; red-tested with a planted
ghost page: exit 1/0). (7) HANDOVER→HISTORY_v1.md retargeting in
adjustVCE, redundancy, vceHutchinson, varComp, poolGroups, verbosity(*),
run_lsqvce, test_adjustVCE — round-note references that genuinely mean
HANDOVER.md stay. ((*) verbosity's line reads "HANDOVER 4d" and means
the v1 doc: retargeted.) (8) Eight HTML pages regenerated from the
corrected help via build_missing_pages (delete+rebuild); audit green.
The r6 StepSchedule/AvgTail claims were NOT touched, per the binding
DECISION RECORD (a).

---

## 10. REVIEW R3 — TEST-CATEGORY SCRUTINY & THE COUNTABILITY GATE

*(Opened 2026-08-19 by user request: scrutinise the toolbox against the
adapted test taxonomy of references/BEST_PRACTICE_v4-1.md §3.3 — six
categories contract / reference / precision / speed / robustness /
vectorisation plus the two cross-cutting tags diagnostics / provenance —
and extend the tests where required. Session model: Claude Fable 5.)*

### 10.1 Method

The count was made MECHANICAL before it was made: every one of the 563
pre-existing test points received a `% categories:` declaration
(body-evidence classification: verifyError/identifier asserts -> contract,
Rel/AbsTol asserts -> precision, timeit -> speed, and so on; try/catch
identifier checks were widened to contract in a second pass), and the new
audit check 16 derives the unit x category matrix from the `lsqvce.<name>`
references in each test body — no second bookkeeping file to drift. The
matrix is printed on every audit run; this section records only the
FINDINGS, because the matrix regenerates itself.

### 10.2 Findings T-01…T-08, dispositions

| # | Finding (measured evidence) | Disposition |
|---|---|---|
| T-01 | Categories not countable: no declarations, §6.1.2 predated the adapted taxonomy (no reference/robustness/vectorisation as named kinds, no cross-cutting tags) | FIXED this round: §6.1.2 rewritten, 563+8 declarations, audit check 16 (red on missing/unknown declarations and on contract gaps without an exemption row) |
| T-02 | Tolerance provenance: 257 Rel/AbsTol assertions, 14 with an origin comment (5%) — the category's defining requirement | CLOSED over two rounds: r1 pinned the `tol:` convention + swept 2 suites; r2 (2026-08-19) annotated ALL remaining 216 sites with per-site derivations, debt 0, staged list retired — check 16d enforces EVERYWHERE (red-tested with a planted site: exit 1/0) |
| T-03 | 8 public units with documented rejection ids but no direct contract test: ReducedContribution (the ONE unit with zero direct tests — both ids unexercised), ar2arma, collectResiduals, residualAcv, varWhiten, whittle, stat, util.hasBackgroundPool. Sub-finding: varWhiten/whittle coverage lived ONLY in the CI-skipped test_parallel | FIXED this round: direct contract tests in the CI-run family suites (combine, signalCore, streamUtils, lsqvce, util); exemption ledger created and EMPTY BY OUTCOME |
| T-03a | Contract DEFECTS found by writing the promised tests: (i) collectResiduals documents badRes for an absent field but threw raw MATLAB:nonExistentField; (ii) varWhiten documents badSigma for "not d x d" but only enforced the PD half; (iii) ar2arma's tooShort ERRORS line described the AR length while the code checks armaN | FIXED this round: (i)/(ii) code now fires the documented id naming the offending field (help text is the contract); (iii) help wording corrected to the code's actual check |
| T-04 | Non-finite propagation asserted only on the accumulator observation path; the signal kernels reject via validateattributes 'finite' but nothing pinned that guard | FIXED this round: testSignalKernelsRejectNonFiniteAtTheDoor (test_reliability) pins psd2acv/levinsonDurbin/armaACV rejection so a refactor cannot drop it silently |
| T-05 | provenance tag: Results carries structural provenance (source, batchFcn, systems, history) but NO run metadata — no toolbox version, no RNG seed record, no options snapshot, no timestamp; nothing in the tree can answer "which run produced this file" a month later | OPEN — PROPOSED WP-44 (API addition, needs user approval; see §10.3). NOT silently implemented: an API change is a decision, not a test extension |
| T-06 | diagnostics tag: measured numbers DO reach the report (report() machine-readable findings; adjustVCE messaging levels carry per-group sigma at level 2) — the relevant tests now carry the tag | NO ACTION beyond tagging; re-review if report surface changes |
| T-07 | armaWhiten batched-equals-column-at-a-time never asserted — the single cheapest vectorisation reference test, assumed since v1 | FIXED this round: testArmaWhitenBatchedEqualsColumnAtATime (test_precision), expectation stated eps-level (10*eps) with the BLAS multi-RHS reason |
| T-08 | reference authorities existed (Ellmer thesis, Baarda/Pope textbook values, closed forms, four Python mirrors) but no standing table of who certifies what | FIXED this round: §6.1.4 authority table; amend it in the same edit that adds/retires an authority |

### 10.3 WP-44 — Run provenance in Results (T-05)

*(Prepared 2026-08-19 per user instruction; status OPEN in the R3 ledger.
Copy-paste prompt in the §9 style:)*

> Close finding T-05: every Results object must carry the metadata saying
> where it came from. Add an IMMUTABLE `meta` property to lsqvce.Results:
> a struct with EXACTLY these fields (add-only across versions, never
> rename — the plot-hook rule): `.toolboxVersion` (char, read from the ONE
> version source, Contents.m line 1), `.matlabRelease` (char,
> version('-release')), `.createdUTC` (char, ISO 8601 with trailing Z),
> `.seeds` (struct of every RNG seed the run resolved, empty struct when
> none — the WP-36 probe contract is the source of truth for which exist),
> `.options` (the RESOLVED option struct the run actually used, not the
> user's raw input). SINGLE POINT OF ACCESS: one private builder
> `+lsqvce/private/runMeta.m` composes the struct; every Results
> constructor call site passes it — no second composer anywhere. The
> Results constructor gains the ONE extra trailing argument (narginchk
> widened accordingly); update every constructor call site in the same
> edit (adjustBlock, adjustVCE, adjustVCEjoint, and any other producer the
> grep finds) — no call site may pass a hand-built struct. Surface: report
> level 0 prints one provenance line (version, release, createdUTC);
> Results help text documents the field list with types and dimensions.
> TESTS (category declarations per §6.1.2): contract | tags: provenance —
> field presence, types, immutability (MATLAB:class:SetProhibited via
> try/catch, NEVER the verifyError output form, §10.4 r1 lesson);
> determinism — two runs with identical inputs/options give identical meta
> EXCEPT `.createdUTC`; robustness — a run with no seeds carries the empty
> struct, not a missing field. DOCS: help text + regenerated HTML page via
> the v1 tooling; Contents.m untouched except if the version line moves
> (it must not). CONSTRAINTS: no behaviour change to any numerical path;
> the audit (incl. checks 14/15/16) and the full suite must stay green;
> byte-identical results for all existing fixtures except the new field.
> Deviations: pin in the round note, approval by merge.

### 10.4 WP-43 round notes

**ROUND 1 (2026-08-19, model Claude Fable 5, branch `claude/wp-43`, PR #20).**
First push: audit GREEN, matlab RED 527/529 — the two failing points were
EXACTLY the two new tests using the output-argument form
`err = tc.verifyError(...)`; both errored (fail + incomplete pairing) on the
function-based TestCase under the hosted R2026a runner. Fix commit 8e52c33:
both sites moved to the repo's try/catch identifier pattern (WP-39).
**LESSON, binding: never take an output from verifyError on the
function-based TestCase — capture via try/catch and assert
`e.identifier`/`e.message` directly.** Counts, verified: static test
points 563 on main -> 571 on the branch (+8, the gap-closure tests);
the CI run executed 529 points with test_parallel skipped.

**ROUND 2 (2026-08-19, same session, branch `claude/wp-43-r2`, riding the
WP-43 DONE flip).** The tol: debt driven 233 -> 0 in one round: three
conservative context-rule passes (82 + 35 + 30 sites) plus a 35-site
correction pass (a `0\b` regex matched the 0 in 0.15 -- every false
"exact identity" comment found and replaced with an honest statistical
one; LESSON: a bulk-annotation rule is code and gets the same review as
code) plus a 30-site hand pass and 5 cluster tails. Every comment states
the origin visible at the site: house bound with the eps/cond argument,
the measuring Python script, the closed form, or the estimator's own SE.
Check 16d flipped from staged to FULL enforcement (every suite), the
report line reworded, and the enforcement red-tested with a planted
unannotated tolerance (exit 1, restored, exit 0).
**r2 CI RED (first push, 4 suites failed to LOAD):** six comments had
been inserted BETWEEN a `...` continuation and its continuation line --
a standalone comment line there is a PARSE ERROR on the R2026a runner
(and two of the six carried the wrong rule text besides). Fix: all six
moved above the statement head with per-site honest wording; audit check
16 gained guard (d0) that REFUSES the pattern tree-wide (a tol: comment
whose preceding line ends in `...`). LESSON: an inserted comment is a
code edit -- it must respect statement boundaries, and the audit now
enforces that mechanically. WP-44 prepared as an
executable prompt (§10.3) per user instruction; ledger row OPEN.

### 10.5 WP-44 round notes

**ROUND 1 (2026-08-19, model Claude Fable 5, branch `claude/wp-44`, riding
the PR #21 merge 7f1856c).** Shipped per the §10.3 prompt: (1) the ONE
composer `+lsqvce/private/runMeta.m` (parses the Contents.m Version line --
the classdef-calls-package-private precedent is jointArcParts, WP-33);
(2) Results: immutable `.meta`, constructor widened narginchk(7,8), the
bare 7-argument path composes empty-seeds/options meta THROUGH runMeta
(never a hand-built struct, never a missing field); (3) both producer call
sites pass composed meta -- adjustBlock: seeds struct() + its resolved `o`;
adjustVCE (blockSolveWith): seeds.jointSeed = jointOpt.Seed (stored as
resolved, [] when unset) + the resolved option groups (vceGroups/vceMode/
jointOpt/shapeOpt/scaleOpt/outerOpt/blockOpt/plotOpt/verbose); (4) report:
additive NV 'Provenance' (default struct(), prints nothing) -- the line
lands as a level-0 NOTE finding code 'PROVENANCE', machine-readable in
V.findings, defensively read (a malformed provenance prints '?', never
breaks the report); R.report() forwards obj.meta, a user-supplied
'Provenance' wins by NV last-match; (5) three tests in test_adjustBlock
(contract+provenance field/type/superset-lock/immutability; robustness
determinism-minus-createdUTC + bare-ctor completeness; diagnostics
report-line + PROVENANCE finding + bare-report absence), all try/catch
identifier pattern; (6) help text Results/report + the two HTML pages
delete+rebuilt via build_missing_pages (WP-41 pattern). DESIGN NOTE: the
provenance line prints through lsqvce.report, not a second printer in
Results.report -- report owns ALL printing (rule 3).

### 10.6 POSTFLIGHT SETTLEMENT (2026-08-20) — ALL D.4 DEBTS CLEARED

The MCP bridge returned 2026-08-20. Full local postflight on the user's
machine (Windows 11, MATLAB R2026a Update 4, PCT present), working copy at
main head 039041c (= the PR #22 merge, so every WP through 44 is in the
run). Executed in four batches to stay under the MCP call timeout
(evaluate_matlab_code was non-functional in this bridge session — every
step ran through run_matlab_file on written driver scripts; the drivers
and their logs were removed after recording, self-deleting cleanup):

| batch | scope | verdict |
|---|---|---|
| A | 10 light suites | 182/182 GREEN, 16.4 s |
| B | 9 heavy suites (adjustVCE alone 45/45 in 140.9 s) | 350/350 GREEN, 152.5 s |
| C | test_parallel ALONE — the CI-skipped suite, real pools (parfor + spmd + threads) | 31/31 GREEN, 11.5 s |
| D | test_speed — reported separately by the runner (§6.1.2); relative budgets on THIS machine | 11/11 GREEN; throughput logged 4.9e6 obs/s |

TOTAL 574/574 — exactly the static test-point count, so the 11-point gap
between CI's 529+31 executed and the static count is fully explained by
the separately-run speed suite; nothing is silently skipped anywhere.

This one green full run on 039041c settles EVERY accumulated D.4 debt
(rows 32–41, 43, 44 and the WP-30-era entries): each debt was "full local
suite green when the bridge returns", main is a superset of every debtor
merge, and the run included the two things CI cannot do — real pools and
machine-local speed budgets. Standing policy unchanged: local postflight
via the bridge before PR-open whenever the bridge is up; CI green remains
the merge gate.

BRIDGE LIMITATION (recorded for the next session):
evaluate_matlab_code failed on every call this session while
detect/check/run_matlab_file worked — drive everything through script
files until proven otherwise.

### 10.7 WP-46 round notes

**ROUND 1 (2026-08-20, model Claude Fable 5, branch `claude/wp-46`).**
User instruction: fix the drifted README, adopt the external geoMap
handover's §2.8 (documentation doctrine), report other adoptables.
Shipped: (1) README de-drifted — status block to the current gate (20
files / 574 points, CI both contexts, branch protection, MIT), Modules
table gains adjustVCEjoint WITH the F-15 production guidance and the
Results.meta provenance, Validation section carries the six-category
taxonomy, the v1 bring-up section retitled as a dated record, and a new
"Correctness history" section pointing at the two CRITICAL review fixes
(F-01/F-03) — a README that quotes v1 claims owes the reader the review
outcome; (2) §6.2.1 gains the ACCURACY block (authority named by its
§6.1.4 row; binding for NEW/TOUCHED functions only — NO bulk migration,
rationale pinned: the source project's 94-file mutation lesson + our own
WP-43 r2 bulk-annotation lessons); (3) new §6.2.4 documentation sync
gate; (4) audit check 17 (README title version == Contents.m; README
"N suite files, M test points" == tree counts) — RED-TESTED three ways
(planted wrong version, wrong count, removed claim; each fired; restored
green). POSTFLIGHT: not applicable — zero .m files changed; the audit is
the relevant gate (ran locally and in CI). Decision also recorded this
round: WP-45 (suite economy) REJECTED by the user 2026-08-20 — the
141 s of test_adjustVCE are test-design redundancy, not toolbox runtime
(the speed suite covers that separately and is green at 4.9e6 obs/s);
rejected, not deferred, so no fresh session re-proposes it.

### 10.8 WP-47 round notes

**ROUND 1 (2026-08-20, model Claude Fable 5, branch `claude/wp-47`).**
User-approved bundle: every WP-46 recommendation except §2.6.3
(verifyAndRecord). Shipped:
(1) **declared-vs-executed per suite** (external 2.6.4): new
`tests/private/declaredPoints.m` (static MATLAB counter, %{ %}-aware —
the block-commented test IS the case a naive grep misses, and its unit
test plants exactly that); run_all results gain `.declared`, the summary
flags a mismatch, gateVerdict rule 6 reds on it naming both numbers. The
Python audit (check 16) counts the same definition with deliberately
SEPARATE code — the second gate reads something the first does not.
(2) **warning inventory + single-allowed-identifier gate** (external 2.5):
measured FIRST — the settlement-run captures and the CI logs show ZERO
MATLAB warnings inside the suites (the only CI 'Warning' lines are runner
environment noise outside the captures), so full enforcement is safe
immediately. run_all turns `warning verbose` on for the run (restored via
onCleanup), parses identifiers from each suite's capture ('<unidentified>'
for a bare Warning line so nothing hides), prints and logs the inventory;
gateVerdict rule 7 allows EXACTLY `lsqvce:internal:testProbe` and — the
adopted "the gate must have something to see" — REQUIRES the probe in
test_ci_gate's inventory when that suite ran (blind capture ≠ warning-free
code). The live probe is testWarningProbeFires, fired un-suppressed.
KNOWN LIMITATION pinned: warnings raised on parallel workers may not
reach the capture.
(3) **§6.2.3 instruments doctrine** (external 2.9): red-test in the same
round, provisional-claim rule, second-gate independence, mutator
discipline — our WP-43 r2 lessons cited as local evidence.
(4) **§6.1.2b diagnosis order** (external 2.10): configuration →
criterion → code; session-state cluster rule.
(5) **§6.1.4 status column** (external oracle register): dated ☑/☐ with
the REFUTATIONS carried (Ellmer closed-loop partial refutation F-15;
Pope-tau scope finding; the wp33 mirror's own faithfulness lesson).
LIVE CATCH: check 17 fired mid-round when the 7 new test points made the
README's 574 stale (→ 581) — the WP-46 gate caught its first real drift
one WP after shipping. Seven new test points, all category-declared,
try/catch identifier pattern.

**Postflight + the r1b/r1c iteration (same day, bridge up, driven via
run_matlab_file — evaluate_matlab_code still dead this session).** The
gate's first full local run REFUTED the round's own pre-measurement: the
inventory held 12 foreign identifiers across 10 suites, not zero — the
"zero warnings" reading had measured the WRONG ARTIFACTS (the settlement
log held only batch C; CI console logs never print the suite captures).
LESSON PINNED (§6.2.3 provisional-claim rule, applied to ourselves): an
inventory claim needs the inventory instrument, not a proxy artifact.
r1b: 33 scoped muteWarn insertions at the provoking tests
(tests/private/muteWarn.m, cleanup-restored) + a suite-scoped mute in
test_adjustVCE's setupOnce (its fixtures run MaxOuter 1–3 to stop BEFORE
convergence by design — breadth justified in the file). r1c after
per-test attribution runs: three placements corrected, two cannot-fire
mutes REMOVED per doctrine, and one real bite — the mute in
testPinvFallbackNeverCrashes blinded that test's own verifyWarning.
LESSON PINNED: NEVER mute an ASSERTED warning — verifyWarning captures
it silently anyway, so an asserted identifier pollutes no inventory.
POSTFLIGHT VERDICT on head 072ce6e (user machine, R2026a U4): batch 1
(all but adjustVCE) 525/525, inventory EXACTLY
{lsqvce:internal:testProbe}, zero declared-mismatches; batch 2
(adjustVCE alone) 45/45, declared 45 == ran 45, inventory empty (the
probe fires only in test_ci_gate — correct); test_speed 11/11 at
4.68e6 obs/s. 570 + 11 = 581 = the static count. All pf47 driver
scripts and .mat artifacts self-removed; working copy restored to main.
RISK CARRIED INTO CI — MATERIALISED AS PREDICTED: the first hosted run
went RED on exactly one runner-specific identifier,
MATLAB:graphics:HardwareUnavailable (headless runner, software OpenGL;
MATLAB fires it ONCE per session at the first graphics use).
Disposition: muteWarn at all three figure-creating tests
(opts/ellipse_confidence_scales_axes, util/ShowExportWritesPng,
util/ShowExportSanitisesPipe) — each can be the session's first graphics
use in a subset run, so all three are genuine provokers; no allow-list
growth. The gate did its job on its first CI outing.
SECOND CI RED (same code, doc-only commit -- run-to-run variance):
lsqvce:normalEq:negativeVtv from
test_precision/testNormalEquationsLoseTheAnswerAtHighCondition -- exact
attribution from the uploaded run_all_log artifact's backtrace. At
cond(A)^2 ~ 1/eps whether v'Pv rounds to -1.2e-16 (warns) or +eps
(silent) is a BLAS coin flip; the local runs and the first CI run landed
silent, so the one-run inventory UNDER-SAMPLED a nondeterministic
provoker. THIRD LESSON PINNED: the muting criterion is CAN-FIRE (the
regime the test drives), not DID-FIRE (one sample) -- sweep by regime,
not by observed inventory alone. r1e disposition: mute at the flaky
provoker AND its regime sibling
(combine/test_qr_beats_normal_when_ill_conditioned, cond 2.4e8, same
coin); the two cancellation tests that ASSERT their warning via
verifyWarning stay unmuted per the r1c rule.

### 10.9 Fixture-reality review and the two proposals (2026-08-20)

*(Written in answer to the user's question: are the adjustVCE tests /
demos closed-loop on real signal values, or random-flow simulation?)*

**What the fixtures ARE.** demo_grace_dataLoader and the test_adjustVCE
suite are CLOSED-LOOP PLANTED-TRUTH tests in the estimation sense: the
observations are generated FROM a known deterministic truth — common
parameters x_c = (1:nc)', deterministic local a-prioris, and a
deterministic planted stochastic law per stream (trueAR: AR(2) for
posC/posD, white rr; trueSigma 0.8/1.3/1.6/1.1; a fixed 3x3 cross-axis
mix) — then the loop must RECOVER the plant: relative stream weights to
the planted ratios, whiteness of the whitened residuals, x against dense
parity, and validate()'s formal-vs-empirical acceptance. Noise
REALISATIONS are seeded-random (reproducible); the planted values are
deterministic. This is NOT random-flow smoke testing.

**What they are NOT.** The DESIGN matrices are random Gaussian
(A = randn, B = 0.3*randn): the parametrisation tables carry SH
degree/order labels, but demoBuildDay reads only the COUNTS — no
spherical-harmonic partials, no orbit geometry, no real GRACE-FO L1B
values anywhere in the tree. Consequence, stated honestly: Gaussian
designs are BENIGN — near-orthogonal, incoherent, well-conditioned —
while a real SH gravity design has strong column correlations
(order resonances), polar-gap conditioning, and orbit-driven leverage
concentration. The VCE loop is therefore proven on the ALGEBRA but
under-tested on the DESIGN STRUCTURE it will meet in production; the
§6.1.4 register carries real data only at the FORMAT level (SINEX
round-trip).

**WP-48 (proposed): the physically structured closed-loop oracle.**
Tier (a), executable with no external data: build A from actual SH
partials (dV/dC_nm, dV/dS_nm) evaluated along a deterministic sampled
near-polar circular orbit — the shAnalysis toolbox on the user's machine
(v3.28) can serve as the partial generator, or self-contained Legendre
recursions — planted coefficient set with a known anomaly, planted
coloured noise; closed loop asserts coefficient recovery, VCE weights
under REALISTIC column correlation, and error-degree-variances against
formal errors. Tier (b), needs user-supplied files: one real GRACE-FO
month through the existing SINEX/NEQ path against a published monthly
solution (register rows added as ☐ 'user to supply'). Tier (a) is the
high-value/no-dependency half.

**WP-49 (proposed, small): authoritative vPv from the residual pass.**
The negativeVtv episode showed l'Pl − b'x cancellation at the regime
edge of the normal backend. The residual pass already recomputes the
same quadratic form stably (regression-locked at 1e-10 agreement); when
a residual pass runs, its vPv should REPLACE the accumulated one as the
quoted number, and the cancellation warning should say so. Considered
AND REJECTED alongside it: iterative refinement on the normal path — a
guard that masks the conditioning instead of the documented cure
(switch to QR); guards-are-not-fixes (§6.3).

### 10.10 WP-48 tier (a) round notes

**ROUND 1 (2026-08-20, model Claude Fable 5, branch `claude/wp-48`).**
Executes the §10.9 proposal, tier (a). Shipped:
(1) `tests/private/shDesign.m` — fully-normalised Legendre recursion +
potential partials along a deterministic circular inclined orbit with
Earth-rotating node; NO random numbers in the design; template help incl.
the new ACCURACY block (first function under the WP-46 rule), authority
"Analytic SH algebra" added to §6.1.4 with measured status.
(2) `tests/test_physOracle.m` (8 points, category-declared): reference
anchors (addition theorem, closed forms), shDesign contract incl. bitwise
determinism, `testDesignIsNotBenign` (the fixture self-guard: if a
refactor ever degenerates the oracle into another Gaussian test, THIS
goes red, not a silent loss of power), and the closed loop — coefficient
recovery inside formal errors (max|z| < 5, mean z² chi² band),
per-degree error-vs-formal ratios, VCE planted-ratio recovery at 0.25
under real column correlation, whiteness/convergence.
(3) Runner + README (21 files / 589 points — check 17 enforced the
update again).
MEASUREMENT DISCIPLINE (§6.2.3, applied): the r0 head carried
PROVISIONAL bands, marked as such in the commit message; the local run
then PINNED the real numbers — addition-theorem residual 4.3e-14 (my
5e-13 guess had a decade of headroom), design max|corr| 0.956 and
single-day cond(N) 4.05e18 vs randn 0.246 / 7.99 (my cond guess of
2.7e5 was FOUR MYTHS too small: one satellite-day is near-singular,
exactly like a real single-arc gravity problem — the 4-day stack
resolves it; this is precisely the structure the Gaussian fixtures
cannot show). Suite 8/8 first pass; full runner 533/533 with the new
suite, inventory exactly the probe, declared==ran==8.
Tier (b) rows sit ☐ in §6.1.4 awaiting the user's files.

### 10.11 WP-49 round notes

**ROUND 1 (2026-08-20, model Claude Fable 5, branch `claude/wp-49`).**
Executes the §10.9 proposal: the authoritative vPv. Shipped:
(1) `NormalEqSystem.solve`: the clamp-to-0 is gone -- a vtv that is
negative OR has lost > 8 digits against l'Pl leaves as NaN, never as a
plausible number (the flat path's measured argument applied to the block
backend), and the warning names the residuals() repair path.
**BREAKING, declared:** the former 'negativeVtv' identifier in the
normalEq namespace is RETIRED into the ONE id
`lsqvce:normalEq:cancellation` -- the sign at eps-garbage level is a BLAS
coin flip (WP-47 r1e evidence), so two ids for one failure mode were an
alias; error identifiers are names (§6.3). The BlockAccumulator
prior-folded path got the same guard under `lsqvce:accum:cancellation`.
(2) `lsqvce.residuals`: NaN-aware supersede logic -- the old mismatch
check was BLIND to NaN (|x - NaN| > tol is false); now a non-finite
accumulated sigma0 raises `lsqvce:residuals:supersedesSol` naming the
numbers to quote, and `res.supersedesSol` flags it mechanically.
(3) Two tests (test_precision): the bare-backend NaN contract, and the
end-to-end block supersede against the QR route.
FOUR ROUND FINDINGS, all pinned:
(a) SESSION-STATE CONFIRMED (§6.1.2b): the first postflight showed 2
fails "no warnings issued" -- my own measurement script had disabled the
two ids and ran the suite before its onCleanup restore. Configuration,
then criterion, then code: it was configuration.
(b) The route-agreement number is 8.0e-5 (this fixture; 3.2e-3 on a
harder 1-sat probe), NOT the 1e-9 I presumed: the residual pass inherits
the normal-path x-hat and pays a second-order penalty at this
conditioning -- still 175x below sigma0's own sampling scatter, so the
superseding number is statistically equivalent. Band pinned at 1e-3.
(c) OBSERVATION (candidate WP-50, not in scope): at extreme dynamic
range (signal 3e4, sigma 1e-5) the QR backend returned sigma0 = 0
EXACTLY -- rotation cancellation produces its own wrong plausible
number. Recorded; the honest-vPv rule may one day need a QR analogue.
(d) Check 14 correctly refused a retired id written in citable
colon-form in help/comments -- retirement notes name the id without the
namespace-colon form.
POSTFLIGHT on head 001af1b (user machine): 535/535 + adjustVCE 45/45
(declared==ran) + speed 11/11 at 4.73e6 obs/s = 591/591; inventories
exactly the probe / empty. Drivers self-removed, working copy restored.

### 10.12 WP-50 round notes

**ROUND 1 (2026-08-20, model Claude Fable 5, branch `claude/wp-50`).**
Executes the 10.11(c) observation as a WP. ROOT CAUSE FOUND on reading
merge(): the block-QR path recovers a pre-eliminated contribution's
residual mass by the NORMAL-form subtraction rho2 = ltl_red - |rhs|^2 --
the very cancellation the other routes guard, inside the QR backend --
and a destroyed (<= 0) rho2 was silently treated as "no residual row",
collapsing vtv toward the EXACT-ZERO sigma0 observed in WP-49; the code's
own comment even named the omission as the catastrophic case. Shipped:
(1) merge(): dof-carrying contributions with destroyed rho2 (<= 0 or
> 15 digits lost) warn ONCE (`lsqvce:qr:cancellation`, at the site with
the numbers) and poison the object (`vtvPoisoned`); dof-free
contributions legitimately carry rho2 ~ 0 and are exempt;
(2) solve(): poisoned -> vtv/s0sq NaN; plus a rotation-round-off
backstop against the augmented column mass (threshold 24);
(3) residuals() supersede closes the loop unchanged (WP-49 machinery);
(4) QRSystem help: WARNINGS + ACCURACY blocks; two tests (the exact
observed deep regime refuse-and-supersede with res.sigma0 ~ 1; the
healthy path locked WARNING-FREE at 1e-10 so a misfiring guard cannot
land silently).
MEASUREMENT TRAIL (three calibration anchors, r0 -> r1 -> r1c):
- bare-QR rotation ramp: rel error only 4e-7 at digitsLost 21.3 -- the
  rotation route is robust far beyond the fold; guard 24 is a backstop.
- the r0 threshold 8 (copied from the normal route) was WRONG here: no
  kappa^2 amplification in this subtraction. r1's 13 then KNIFE-EDGED a
  proven-good fixture: the WP-49 test's QR reference folds at MEASURED
  13.4 digits (ltl_red 7.19e16 / rho2 3.06e3) with sigma0 good to 8e-5
  -- the eps*10^dl bound is ~50x pessimistic there. Destroyed anchor:
  the deep probe at ~19 digits, rho2 NEGATIVE. PINNED 15: 1.6 digits
  above measured-good, ~4 below measured-destroyed.
- r1 lesson (small, pinned): arrayfun over a void method errors on
  output collection -- the warning assertion uses a plain-loop helper,
  and the fold's warning fires at addSystem time, not solve.
POSTFLIGHT on head efec974 (user machine): 537/537 + adjustVCE 45/45
(declared==ran) + speed 11/11 at 5.10e6 obs/s = 593/593; inventories
exactly the probe / empty. Drivers removed, working copy restored.

### 10.13 WP-51 round notes

**ROUND 1 (2026-08-21, model Claude Fable 5, branch `claude/wp-51`).**
User request: turn the WP-48 tier (a) oracle into a GRACE-like DEMO, with a
different local-unknown count per satellite and a different number of arcs
AND arc length per day.
DESIGN DECISION (recorded because it moves a file across the tree): the demo
needs the SH partial generator, and duplicating a numerical kernel breaks
one-way. Rejected: a second copy in examples/, and promotion to a public
`lsqvce.shDesign` (the toolbox is a VCE engine; it should not own gravity
machinery, and §4 is user-governed). Adopted: MOVE it to
`examples/demo_gravity_design.m` and have the oracle suite call it there --
the house pattern already has tests consuming examples/ fixtures
(`test_adjustVCE` -> `demo_grace_dataLoader`). Public API surface unchanged;
ids renamed to `lsqvce:gravityDesign:*`.
RAGGED ARCS -- the mechanism: `makeLoader` resolves `'arcs'` ONCE without
reading a day, and `batchOf` ERRORS on an arc id beyond the day's arcIdx, so
a fixed declaration cannot express a ragged schedule. But the accumulator
probes `'meta'` PER DAY (`dataFcn('meta', d, [], uo)` in accumDay). The
loader therefore wraps `.closure` (the field the engines actually drive,
BlockAccumulator line ~1290) and answers that probe with the day's ACTUAL
arc list; absent arcs are never requested. Arc lengths vary freely -- they
are just row counts.
TOPOLOGY: GRACE-C (potC, 2 locals/arc), trailing GRACE-D (potD, 3), and the
inter-satellite potential DIFFERENCE as an energy-balance KBR stream (kbr, 1
local/arc, design A_D - A_C so it couples both orbits). nl = 32/48/16 --
three different local block sizes in one adjustment. Schedule 3/4/2/4/3
arcs, absent BY ID (day 3 = [1 3], day 5 = [2 3 4]), lengths 140-220.
VERIFICATION STATE -- READ THIS BEFORE TRUSTING THE DEMO'S NUMBERS:
the MCP bridge dropped MID-TASK, so this round has NO local postflight.
CONFIRMED IN A SECOND SESSION (2026-08-21, later): a follow-up session
adopted this branch specifically to settle the debt with a working bridge
and found the bridge (Filesystem AND MATLAB tools) gone MID-SESSION as
well -- the outage is environmental, not a one-off. The debt therefore
stood until 2026-08-21 (later), when the USER ran the demo on his machine
and reported the output. **DEBT SETTLED, numbers pinned into the demo's
ACCURACY block (WP-53):** sigma0 = 0.9998; stream sigmas recovered
0.9956 / 1.9999 / 0.5096 against planted 1.0 / 2.0 / 0.5 (ratios 0.996 /
1.000 / 1.019) -- the weights come back under REAL column correlation;
max|z| = 2.67 over 77 coefficients and mean z^2 = 1.087 (0.5 sd from the
ideal 1.00 at sd = sqrt(2/77) = 0.16) -- the formal errors DESCRIBE the
actual errors; the planted anomaly (5,3,C) truth +0.1764 recovered as
+0.2843 +- 0.0812, i.e. 1.33 sigma; degree ratios 0.23 / 1.29 / 0.79 /
1.00 / 1.07 / 1.15 / 1.00 for n = 2..8 -- the low n=2 is the expected
small-sample scatter of a degree with five coefficients (chi2, 5 dof,
roughly 0.1..2.2 at one sigma), not a defect. The quality claim is now
MEASURED. Original text of the debt, for the record: CI proves the
demo RUNS and that the fixture structure is what it claims; the RECOVERY
QUALITY it prints is still unmeasured. Settle it by running
`[R, score] = demo_gravity;` locally and reading sigma0, score.zMeanSq,
score.zMax, score.degRatio and score.sigmaRatio against the planted
1.0/2.0/0.5 -- then pin or open round 2. That second session also
independently built a variant of this WP (public `lsqvce.shDesign` +
two-system fixture) and DISCARDED it in favour of this branch: one design
per problem, and this one is richer (three systems, arcs absent by id).
Compensation: two smoke tests were added so CI itself executes the demo
end to end (`testDemoGravityRunsEndToEnd`) and asserts the fixture
structure (`testDemoGravityFixtureStructure`) -- both green on the hosted
runner, which proves the demo RUNS and that the ragged/per-system
structure is what it claims. What is NOT verified: the QUALITY of the
recovery the demo prints (recovered-vs-planted weights 1.0/2.0/0.5,
max|z|, degree-variance ratios). NO statistical band is pinned anywhere in
this round -- deliberately, per §6.2.3. DEBT: run `demo_gravity` locally,
read the three sections, and either pin the numbers in the demo help
(COST/ACCURACY) or open a round 2 if the recovery is poor. Prior evidence
that it should be sound: WP-48 measured VCE ratio recovery at 0.25 RelTol
on a comparable structured design (§10.10).

### 10.14 WP-52 round notes

**ROUND 1 (2026-08-21, model Claude Fable 5, branch `claude/wp-52`).**
User question after WP-51: can `makeLoader` support ragged schedules
NATIVELY instead of leaving `raggedClosure` to the user? ANSWER FOUND BY
READING, not by guessing: the gap was never in the protocol. EVERY engine
probes `dataFcn('meta', d, [], uo)` PER DAY (BlockAccumulator.accumDay,
residuals, redundancy, streamContribution, diagnoseAccum, adjustVCE x2,
adjustVCEjoint, and PrefetchLoader passes it through) -- `metaOf` simply
threw the day away. Two of the three asks needed NOTHING: unequal arc
lengths were always native (arcIdx entries are row-index vectors), and
per-SYSTEM local counts were always native (`numLocal` is per loader).
Shipped, all ADDITIVE (no id renames, existing loaders untouched):
(1) `'arcs'` accepts `@(day)` / `@(day,obsType)` or `'fromDay'`;
`resolveArcSpec` + `checkArcIds` validate what a resolver returns.
(2) A day may declare `.arcId`; `batchOf` then addresses BY ID and an
absent id returns an EMPTY batch -- verified first that all eight arc
consumers already guard `isempty(arc)`, which is why this is safe. The
positional out-of-range error SURVIVES for static schedules: a
mis-declared static count is a different mistake from a ragged one.
(3) `'localsPerArc'` gives makeLoader the bookkeeping: it walks the day
list, assigns each (day, arc) a global ordinal and column offset, DERIVES
nl (so it can no longer be stated twice and disagree), and hands each day
its slice as buildDay's 6th argument `arcInfo` (`.arcId .offset .count
.ordinal .nl .offsetOf`). It REFUSES the two unknowable cases: no
`'days'`, or `'fromDay'` (nl must be known without reading data).
(4) The demo was rewired onto the native path and `raggedClosure` plus the
`arcOrd` lookups in the builder were DELETED -- that deletion is the proof
the support is native, and `testDemoGravityUsesTheNativeRaggedPath` locks
it (the loader must carry the plan, and the source must not wrap
`.closure`).
FOUR CI CATCHES (the bridge stayed down, so CI was the only route):
`L.meta()` still called the old 2-arg `metaOf`; `clear()` on a Map is not
a teardown; `adjustBlock` has no `KeepLocal`; and -- the instructive one
-- `L.meta` CAPTURED `L.day` before that field existed (an anonymous
function captures values, not future struct fields), fixed by hoisting the
accessor into a variable. A fifth catch was in my own TEST data, not the
code: per-arc `[1, tau]` locals absorb a common constant and trend
EXACTLY, so the fixture was rank deficient by construction and returned
NaN -- the commons are now harmonics of a global clock. WP-49's
`supersedesSol` surfaced it immediately instead of printing a plausible
number.
TOLERANCE LESSON (§6.2.3 applied to myself): the recovery test first
carried a `RelTol` of 1e-6 that I had reasoned to rather than measured;
CI measured 1.5e-2. The fix was NOT a looser band but a better claim --
the recovery is now asserted against the FORMAL errors (max|z| < 4) plus
sigma0 ~ 1, which is the statement worth making and does not encode one
noise draw.
DEBT: no local postflight -- the MCP bridge was down for this entire WP
(third consecutive session). CI green is the only evidence; the WP-51
recovery-quality debt (§10.13) is untouched and still open.

### 10.15 WP-53 round notes

**ROUND 1 (2026-08-21, model Claude Fable 5, branch `claude/wp-53`).**
Two user findings in one line: "demo_gravity does not know a second
output?" plus the output of a real run. Shipped:
(1) `demo_gravity` now RETURNS the self-assessment it prints (second
output `score`: sigma0, per-stream planted/recovered/ratio, zMax,
zMeanSq, the anomaly with its z, and the degree ratios). Printing a
number the caller cannot get hold of is a demo that can only be READ,
never CHECKED -- the fix is one struct assembled from the values already
computed for the printed lines, so there is ONE source and no recompute.
(2) The measured numbers of the default seeded run are pinned in a new
ACCURACY block in the demo's help, and §10.13's quality debt is settled
with them (see there). They are QUOTED, not asserted: the test locks the
structural contract of `score` (documented fields present, sigma0 equal
to the result's, ratio = recovered/planted, z recomputed) and leaves the
statistics to the formal-error reading, because the test runs MaxOuter 1
and would not reproduce a MaxOuter 4 run's numbers anyway.
READING OF THE RESULT (worth recording): the weights come back under real
order-block correlation (worst ratio 1.019 on the KBR stream, which
carries the fewest locals), and mean z^2 = 1.087 says the formal errors
describe the actual errors on a physically structured design -- the claim
WP-48 set out to make and could only make on Gaussian fixtures before.
DEBT: still no local postflight (bridge down for a fourth session); CI
green remains the only evidence for this WP.

### R3 ledger

| WP | title | depends on | status |
|----|-------|-----------|--------|
| 43 | Test-category audit: adapted taxonomy, countability gate (check 16), gap closure | — | r1 DONE (2026-08-19, PR #20 merged by the user; r1 red 527/529 fixed same day, final head CI GREEN 529/529; POSTFLIGHT DEBT → SETTLED §10.6 per D.4); r2 PR #21 READY (tol sweep 233→0 + full enforcement + guard d0; head efbe8be CI GREEN both contexts; DONE flips on the user's merge; round notes §10.4) |
| 44 | Run provenance in Results (T-05) | 43 | DONE (2026-08-20, PR #22 merged by the user, merge 039041c; CI green FIRST pass 4ea32c0; POSTFLIGHT DEBT → SETTLED §10.6; round note §10.5) |
| 45 | Suite economy (adjustVCE test-run cache) | — | REJECTED (user, 2026-08-20): test-environment cost only, toolbox speed separately gated by the speed suite; see §10.7 |
| 46 | README truth + documentation doctrine (external §2.8 adoption) | — | DONE (2026-08-20, PR #24 merged by the user, merge c11512d; CI green FIRST pass; round note §10.7) |
| 47 | Gate hardening: declared-vs-executed, warning inventory, instruments/diagnosis doctrine, authority status column (external 2.5/2.6.4/2.9/2.10 + oracle-register adoption; §2.6.3 excluded by user) | 46 | DONE (2026-08-20, PR #25 merged by the user, merge 8325a26; two real CI catches on the way — headless-graphics id (r1d) and the nondeterministic negativeVtv BLAS flake, disposed by REGIME (r1e); three pinned lessons; round note §10.8) |
| 42 | Monotone / trust-region spectral alternation (F-15 cure) | — | BACKLOG (user decision reconfirmed 2026-08-20 “for the time being”; entry criteria and the mirror-first success criterion pinned in §9 WP-40/§3222 entry; production guidance — open loop or ≤3 damped passes — covers operation meanwhile) |
| 48 | Physically structured closed-loop oracle: tier (a) deterministic SH-partial design; tier (b) real GRACE-FO month vs published solution | — | tier (a) DONE (2026-08-20, PR #27 merged by the user, merge f1bb951; round note §10.10); tier (b) BACKLOG (user decision 2026-08-20; §6.1.4 rows ☐ remain the entry point) |
| 49 | Authoritative vPv from the residual pass on the normal backend (§10.9) | 48 | DONE (2026-08-20, PR #28 merged by the user, merge a2b6430; BREAKING id retirement declared; round note §10.11) |
| 50 | Honest vtv on the QR backend (10.11 finding c): the merge-fold cancellation was the true sigma0=0 mechanism | 49 | DONE (2026-08-20, PR #29 merged by the user, merge 4cbc280; postflight 593/593; three-anchor threshold calibration; round note §10.12) |
| 51 | demo_gravity: GRACE-like physically structured demo (per-system local counts, ragged arcs, ragged arc lengths) | 48 | IN PROGRESS — PR #30 (2026-08-21, branch `claude/wp-51`; CI GREEN incl. an end-to-end demo run; **NO local postflight — bridge dropped mid-task**; recovery QUALITY unmeasured, see §10.13) |
| 52 | Native ragged schedules in makeLoader (per-day arc resolution, arc ids, owned per-arc local bookkeeping) | 51 | DONE (2026-08-21, PR #31 merged by the user, merge 30e9c31; round note §10.14; postflight DEBT -- bridge down) |
| 53 | demo_gravity returns its self-assessment; measured numbers pinned (settles the §10.13 quality debt) | 51 | IN PROGRESS (2026-08-21, branch `claude/wp-53`; round note §10.15) |
