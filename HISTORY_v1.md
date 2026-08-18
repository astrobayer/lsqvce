> **RENAMED 2026-08-17: this file was `HANDOVER.md` until the single-folder
> restructure (user decision — `lsqvce2/` dissolved, one tree). It is the
> FROZEN v1 HISTORY. The LIVE status-and-plans document is `HANDOVER.md` in
> this folder. In-code prose references like "HANDOVER §4d" or "HANDOVER
> roadmap item 0" refer to THIS file (retargeting them is WP-41).**

# lsqvce — Handover & Roadmap

**THIS IS THE ONLY STATUS-AND-PLANS DOCUMENT. Do not create a second one.**
(User instruction, 2026-07-18.) Side documents drift out of step with each other
and with the code, and then nobody knows which one is current — that is exactly
how this file came to claim a fix was DONE while the bug was still in the code
(§1.0, defect 4). If a plan, a build sketch, or a "start here" brief is worth
keeping, it belongs in a section of THIS file. Two previous side documents,
`ELLMER_BUILD_PLAN.md` and `NEW_CHAT_START_HERE.md`, were folded in and are gone;
their durable content is in §2 (the Ellmer key facts) and §0/§4b (the working
loop). The one legitimate companion is `docs_src/audit_tree.py`, which is not
prose — it is an executable check that re-derives the claims made here.

This top section is the CURRENT STATE and the OPEN TOPICS. The dated blow-by-blow
history is in the APPENDIX at the bottom (newest first) — you rarely need it, but
it is the authoritative record of why each change was made.

Last updated: 2026-07-19. **Suite green: 497/498, 0 failed, 127.6 s** (§1.STATUS).
New here? Read §0a, then §1.OPEN. Nothing else is required to start work.

> **v1 IS FROZEN FOR FEATURE WORK (user decision, 2026-07-19, evening
> session).** A three-expert design review of the whole toolbox concluded
> that the numerical kernels are sound but the architecture is not (three
> entry points with subset semantics, four accumulation engines, snooping on
> the wrong path, a non-scaling default diagnostics chain, the Ellmer model
> demoted to opt-in). A successor is being built from a clean redesign.
> **All new development is tracked in `HANDOVER.md (this folder)`** — that file
> is the single status-and-plans document of the v2 effort, exactly as this
> file is (and remains) the record of v1. Open items A and B below (§1.OPEN)
> are NOT being fixed here; they are subsumed by v2 work packages 12 and 15.
> Only defect fixes that block a v1 user still land in this tree.

---

## 0a. START HERE IN A NEW CHAT

_Read this section, then §1.OPEN. Everything else is background._

### The standing instructions (from the user, verbatim in substance)
These are the PROJECT RULES. They do not expire between chats.

1. **Role.** Expert MATLAB programmer and least-squares-adjustment specialist.
2. **Think step by step, and MOTIVATE every answer.** Say why, not just what.
3. **Develop a SINGLE POINT OF ACCESS where useful.** One name per thing, one
   place per decision. (This rule is why aliases were removed, why every stage
   label goes through `verbosity.banner`, and why the block-spec field list is
   audited across its four definitions.)
4. **Always prepare unit tests.** New behaviour ships with a test that would
   fail without it.
5. **Every function gets full help text**: purpose, INPUT and OUTPUT with TYPE,
   DIMENSION and DESCRIPTION, errors it raises, an example, and a
   `Created <date> (CEST, UTC+02:00) with assistance from Claude (Anthropic),
   model <version>` stamp.
6. **Toolbox documentation is HTML reachable from `doc`.** A new function needs
   a page in `html/` plus entries in `Contents.m`, `helptoc.xml` and
   `html/index.html`. `docs_src/build_missing_pages.py` generates the page;
   `audit_tree.py` check 4 verifies the registration.
7. **Keep THIS handover current, clean and tidy.** One status document (see the
   note at the top of the file), and the test count lives in §1 only.
8. **Avoid aliases when refining functions.** A renamed thing has ONE name; the
   old one is deleted, not deprecated.
9. **Follow the MATLAB coding standards**
   (https://github.com/matlab/rules/blob/main/matlab-coding-standards.md).
10. **Always provide a WRAP after changes.** The user downloads the zip and runs
    the suite in MATLAB, then returns the log.

### How the loop works
The user uploads `lsqvce.zip`; extract and work from it. There is NO MATLAB in
the sandbox and normally no access to the user's filesystem. Validate three
ways: STRUCTURAL (function/end balance against the pristine upload), NUMERICAL
(Python vs exact formulae), and the STATIC AUDIT (`python3
docs_src/audit_tree.py`, exit 0 = clean, 11 checks). **The user's MATLAB run is
the ground truth** and has repeatedly caught what the sandbox cannot see.

Ship the ENTIRE toolbox as ONE zip per round via `present_files` -- never a
partial set, which is what makes the two copies diverge. Full detail in §4b,
including the part about large pastes being lost in transit: logs come back as
an UPLOADED .txt file.

### Hard constraints
- **No JavaScript/Java anywhere** in shipped files, documentation included.
- Run the full suite with **`run_all_lsqvce`**, never `runtests('tests')`.
- `tests/` is FLAT (no subfolders) since 2026-07-18.

### First action in a new chat
Read §1.OPEN for what is outstanding. If the suite has not been run against the
current zip, ask for `run_all_lsqvce` and `tests/last_run_all.txt` before
starting anything new -- and read the log per §4d (warning types, log size,
what did NOT run), not just the pass count.

### The five lessons this project paid for
Each cost a real defect. They are the reason the audit script exists.
1. **A rule recorded only in the handover is not in force.** Put it in the code,
   give it a test, then write it down (§1.0 defect 4).
2. **A count of failures is not a statement about coverage.** `test_util` once
   failed to LOAD and the gate still said "0 failed" while 41 tests vanished.
   Ask what did not run.
3. **A checker that strips comments or strings is right for COUNTING BLOCKS and
   wrong for MATCHING TEXT** -- and one that ignores line CONTINUATIONS misses
   whatever sits on the second line. Three separate defects, §4h.
4. **A hardcoded verbosity is a bug waiting to happen.** Gate on the structural
   reason, never on a literal `false` (three instances in two days).
5. **Verify every match BEFORE mutating anything.** An edit script that writes
   once at the end silently discards the substitutions that already succeeded
   when a later one fails.

## 0. WHAT THIS TOOLBOX IS

`lsqvce` is a MATLAB package (`+lsqvce/`) for least-squares adjustment with
variance-component estimation (VCE), aimed at GRACE-FO gravity recovery but usable
for any large block-structured adjustment. Core capabilities:

- **Three solve entry points**, increasing in capability:
  - `lsqvce.adjust`      — flat, in-memory least squares (small problems).
  - `lsqvce.adjustBlock` — streaming block-arrow QR (common + per-arc locals);
    O(u_c²) memory, independent of the number of arcs.
  - `lsqvce.adjustVCE`   — the iterative separable engine: alternates SHAPE
    (temporal AR whitening per stream) and SCALE (Helmert VCE between groups)
    around the block-arrow solve until both converge.
- **VCE**: exact Helmert (`lsqvce.vce` / `lsqvce.varComp`) and a matrix-free
  Hutchinson variant (`lsqvce.vceHutchinson`, or `varComp('backend','hutchinson')`)
  for large u.
- **Stochastic modelling**: AR/ARMA fitting (`selectArAcv`, `ar2arma`,
  `levinsonDurbin`, `armaACV`), whitening, PSD/ACV conversion.
- **Diagnostics**: redundancy/reliability, residuals, data snooping (`outliers`),
  convergence diagnostics on `R.history`.
- **Loaders**: `lsqvce.makeLoader` wraps a user day-builder into the loader
  interface (per-arc A/B/l, per-epoch covariance, parametrisation labels, plot
  hooks). `demo_separable_data.m` is the worked template for a real pipeline.

### Environment & workflow (IMPORTANT — read before editing)
- **There is no MATLAB/Octave in this sandbox, and normally no access to the
  user's filesystem either.** All validation here is (a) STRUCTURAL —
  function/end balance vs the pristine original — (b) NUMERICAL — Python
  (numpy) against exact formulae — and (c) the STATIC AUDIT
  (`docs_src/audit_tree.py`, §4f). The USER runs the real MATLAB suite and
  returns the log; **those runs are the ground truth** and have repeatedly
  caught things the sandbox cannot see.
- **THE LOOP IS THE ZIP ROUND-TRIP — see §4b.** Ship the ENTIRE toolbox as ONE
  zip every time (currently 344 files, 1.5 MB) via `present_files`, never a
  partial set of changed files. Work tree: unpack the uploaded zip under
  `/home/claude/work/`; keep the pristine upload unpacked separately (e.g.
  `/home/claude/orig/`) so every change can be diffed against it. Output ship:
  `/mnt/user-data/outputs/lsqvce.zip`.
- **Run the full test suite with `run_all_lsqvce`** (18 files; 459 executed
  points as of 2026-07-18). `runtests('tests')` SKIPS the `tests/`
  subfolder — do not use it as "the whole suite".
- Stamp new files: `Created 2026-07-18 (CEST, UTC+02:00) with assistance from
  Claude (Anthropic), model Claude Opus 4.8.`
- **Logs come back as an UPLOADED .txt FILE, not a paste.** Large pastes are
  silently lost in transit (§4b) — the full suite log is ~105 KB. If something
  arrives blank or truncated, ask for the file again rather than inferring the
  result; `view /mnt/user-data/uploads/<name>` reads reliably.

---

## 1. OPEN TOPICS

### 1.STATUS — WHERE THINGS STAND (2026-07-19)

**The suite is GREEN: 497 / 498 passed, 0 failed, 1 incomplete, 127.6 s.**
(MATLAB-confirmed 2026-07-19 16:28; log 107,976 bytes, 54 warnings / 12
categories. Per-suite seconds are logged since §1.1e -- `test_separable` alone
is 108.9 s of the 127.6 s, recorded as open task 6. The incomplete is
`test_parallel/testPrefetchServesFromCacheWhenPrimed`, "Filtered by
assumption" -- the long-known `backgroundPool` filter, not a regression.)
All 18 test files LOADED; check that, not just the failure count (§4d).

Warning inventory, stable across every run this session: **51 lines, 12
categories**, all deliberately provoked --
13 VCE non-convergence and 13 loop non-convergence (both from low `MaxIter` /
`MaxOuter` set to keep tests fast), 6 singular local block (the collinear
fixture), 4 zero redundancy, 3 joint non-settling, 3 datum auto-detect,
2 ARMA non-stationary, 2 rank-deficient (MATLAB's own), 2 `v'Pv` cancellation,
1 plot-hook failure test, 1 rejection cap, 1 maxIter cap. **No new type.**

WHAT IS FINISHED
- The Ellmer (2018) stochastic model, Phases A-E (§2).
- Streaming Hutchinson VCE (§1.1) and the contribution matrix, Levels 1 and 2
  (§1.2) -- both VERIFIED in the code this session, not assumed.
- The interface: one name per thing, `nc` derived from the systems, no aliases.
- The messaging layer: three verbosity levels, plain-words stage labels, shared
  block renderers, and -- as of this session -- a wiring test that reads the
  captured log back and asserts what appears at each level.
- Tests flattened into one folder, release-named files renamed to their subjects.
- `docs_src/audit_tree.py`: 11 static checks, each verified to fire on a broken
  tree AND stay silent on a healthy one.

WHAT IS NOT
- The four tasks below. Nothing else is in flight.

### 1.OPEN — THE AUTHORITATIVE TASK LIST (2026-07-19)
Everything below §1.OPEN is BACKGROUND: closed items kept as the record, plus
the detail behind each task. **This list is what is actually outstanding.**
Order agreed with the user: documentation goes LAST, because it should describe
a settled interface rather than chase one.

**ALL SEVEN ROADMAP TASKS ARE CLOSED as of 2026-07-19** (tasks 3, 5, 6 and 7
awaiting their confirming MATLAB run). What remains is listed under SMALLER /
USER-ONLY below, plus the standing warning that `docs_src/build_manuals.py` is
GUARDED and must be back-ported before it is ever run (§1.3b).

| # | task | why now / blocked by |
|---|------|----------------------|
| 1 | ~~**SIMPLIFY THE USER INTERFACE, then template it**~~ | **DONE 2026-07-19, MATLAB-verified 498/499 + `doc lsqvce` works.** Spec + decisions §1.1a; built over four rounds -- §1.1b template + setup/audit fixes, §1.1c `convergence()` moved to `Results` + `summary()` added, §1.1d registration + demo cleanup + 12 stale HTML examples, §1.1e timing instrumentation. |
| 2 | ~~**Interface guide: preparing `buildDay` and friends**~~ | **DONE 2026-07-19, MATLAB-verified 497/498 (§1.2a, run §1.2b).** Written as user-manual **Parts IV-B and IV-C** (user decision: in the manual, everything a user writes, worked wrong-then-right examples). Also fixed 14 dead fluent-method names in shipped HTML and added audit check 13 to stop the class recurring. |
| 3 | ~~**Merge the two test runners** into one entry point~~ | **DONE 2026-07-19 (§1.3a), awaiting the MATLAB run.** `run_lsqvce('all')` delegates to `run_all_lsqvce`; the old `'all'` ran `runtests` on a FOLDER, i.e. the §1.7b trap. |
| 4 | ~~Purge SHAPE/SCALE from the help text~~ | **DONE 2026-07-19.** ~60 stage-sense occurrences reworded; the seven that remain are ordinary English ("the shape of `systems`", "SCALE each epoch by its Cholesky factor") or `verbosity.m` quoting the retired form to explain what it replaced. Audit check 11 guards user-visible strings from a relapse. |
| 5 | ~~**§1.3 documentation rebuild** — the master PDF~~ | **DONE 2026-07-19 (§1.3b-d), awaiting the MATLAB run.** `lsqvce_manual.pdf`, 18 pages, built with ReportLab (`python3 docs_src/build_master_pdf.py`). The three old PDFs are deleted and their builder retired. |
| 6 | ~~**`test_adjustVCE` runtime: 109 s of a 120 s suite**~~ | **DONE 2026-07-19 (§1.4b), awaiting the MATLAB run.** The dataset is built once in `setupOnce` and reused. Judge it by the PER-SUITE SECONDS, not the pass count. |

| 7 | ~~**Documentation completeness sweep**~~ | **DONE 2026-07-19 (§1.4c), awaiting the MATLAB run.** `review_mfiles.py` 121 -> 1 observation; the remainder is a structural remark about `BlockAccumulator`'s length, not a defect. |

WHAT IS LEFT (user decision 2026-07-19: **b and c are DROPPED**; real-data
processing happens outside this development). Detail in §1.4:

| # | item | verdict |
|---|---|---|
| **A** | **wire `poolGroups` into the `adjustVCE` chain** | **OPEN, recommended.** A missing DATA PATH, not a convenience: `R.history` records no per-group redundancy, which is exactly what pooling needs. |
| **B** | **data snooping is NOT connected to the block-arrow path** | **OPEN.** `lsqvce.snoop` is flat-path only; `adjustBlock`/`adjustVCE` have no outlier rejection at all. |
| C | `build_manuals.py` is GUARDED | **STANDING HAZARD, not obsolete** -- re-verified 2026-07-19, see §1.4C. |
| — | ~~namespace sub-packaging~~ | **DROPPED** by the user. |
| — | ~~split `BlockAccumulator.m`~~ | **DROPPED** by the user. |
| — | ~~run the real GRACE-FO loop~~ | **OUT OF SCOPE** -- happens outside this development. |

### 1.4c TASK 7 (2026-07-19) — documentation completeness, 121 -> 1

**HELP TEXT AND COMMENTS ONLY, no executable line touched. Expect 497/498, and
`test_adjustVCE` still fast (task 6 unaffected).**

`docs_src/review_mfiles.py` went from **121 observations to 1**, and the one
left is `BlockAccumulator.m` flagged as "long with few section banners" -- a
structural remark, not a defect.

DONE, in priority order rather than alphabetically. The ranking was measured:
cross-reference each file against `Contents.m` to see whether it is PUBLIC, then
count its ids. That put `BlockAccumulator` (23), `Combine` (16), `BatchIndex`
(15) and `Accumulator` (14) first -- 68 ids in four classes a user reaches
through `R.acc`.
- **ERRORS sections: 55 files -> 0.** Every id is grouped BY CAUSE, not listed
  alphabetically, e.g. `BlockAccumulator` splits into "your system declaration
  is wrong", "a day from YOUR buildDay does not fit it", "the prior is
  malformed", "asked for something the accumulator does not have". A flat list
  of 23 ids would be technically complete and useless.
- **Provenance: 33 files -> 0**, stamped after `See also` so the help footer
  MATLAB links from is not broken.
- **INPUT/OUTPUT: 13 -> 0**, with sizes and types.

**CHECK 14 EARNED ITS KEEP AGAIN, IMMEDIATELY.** I documented
`lsqvce:vce:badSigma`; the code throws `badSigma0`. My collection regex had
truncated at the digit. The audit caught it in the same round -- and a follow-up
scan found 8 ids toolbox-wide that contain digits, so the same slip could have
been repeated. Without check 14 this would have shipped as a help block promising
an identifier no user could ever catch.

REMAINING (not defects, recorded so they are not rediscovered)
- `BlockAccumulator.m` is 631 code lines with few `% ---` banners. Splitting it
  is a real refactor with real regression risk, and it is well commented
  internally; NOT worth doing for tidiness alone.
- The 12 numerical primitives (`levinsonDurbin`, `whittle`, `armaACV`,
  `varACV`, ...) now have ERRORS sections too, though a user reaches them only
  through `covModel`.

### 1.4b TASK 6 (2026-07-19) — test fixture cached; stale H1s; check 15

**EXPECT 497/498 UNCHANGED, and `test_adjustVCE` MUCH FASTER.** The whole point
is the per-suite seconds, so compare THOSE, not the pass count. Baseline to beat:
`test_adjustVCE 109.1 s` of a `119.7 s` suite.

**THE CAUSE, MEASURED BEFORE TOUCHING ANYTHING.** `setupOnce` already called
`demo_grace_dataLoader()` -- and kept only `cfg`, THROWING THE LOADERS AWAY. So
`makeVCEBuilder` regenerated the entire dataset (4 days x 2 arcs x 300 epochs,
~16,800 rows) on each of its 27 calls, and 7 further tests called the generator
directly. All 27 fixture users also `.run()`, so this was never "wasted data for
tests that do not need it" -- it was the same 4 days rebuilt 27 times.

FIX: keep BOTH returns in `setupOnce` and reuse them.
**WHY SHARING LOADERS IS SAFE, since it is not obvious:** a loader is a STRUCT of
closures over a shared `containers.Map` cache. The struct is a VALUE, so the two
tests that do `L{i}.plotPostSolve = ...` get their own copy and cannot leak hooks
into the next test. The Map is a HANDLE and IS shared -- but it only ever holds
days built by a DETERMINISTIC `buildDay`, so a cache hit is byte-identical to a
rebuild. That is the same property the accumulator already relies on for its
residual pass. Checked: no test asserts on buildDay call counts, and none
depends on a cold cache.
Two read-only direct callers converted; the other five genuinely need a fresh or
specially-configured dataset and were left alone.

**A REAL WEAKENING CAUGHT AND REPAIRED.** `test_deterministic` ran two builders
and compared them. With a shared fixture both would use the SAME cached loaders,
so it would have silently degraded from "the generator is deterministic" to "the
solve is deterministic" -- a cached day trivially equals itself. It now asserts
BOTH: same-loaders/same-solve, AND a freshly built dataset reproducing the shared
one. That costs one extra solve and is worth it, because generator determinism is
exactly what the accumulator's second pass depends on.

**THREE STALE H1 LINES** from the 2026-07-19 rename: `demo_grace_dataLoader.m`
still announced `DEMO_SEPARABLE_DATA`, `demo_grace.m` said `DEMO_SEPARABLE`, and
`test_adjustVCE.m` said `TEST_SEPARABLE`. `help <file>` printed the old name.
Check 2 did not catch them because it guards removed FUNCTION names in code, not
the prose of an H1 line. **NEW AUDIT CHECK 15** closes the class; verified both
ways (fires on a reverted H1, silent on the clean tree).

### 1.4a REVIEW SWEEP (2026-07-19) — what was checked, found and fixed

A full scrutiny pass over all 80 `.m` files (21,599 lines), the HTML and the PDF.
Built `docs_src/review_mfiles.py` for it -- a REVIEW tool, deliberately not a
gate, so it reports and a human judges.

**TWO OF MY OWN CHECKS WERE WRONG BEFORE THE CODE WAS.** Recorded because the
lesson keeps recurring:
- The reviewer first reported **77 "H1 mismatch"** findings. All false: a package
  function's H1 is conventionally `LSQVCE.NAME`, which is correct MATLAB style.
  Fixed to compare the last dotted component -> 0 findings.
- It then reported 13 files with no INPUT section. Several were false: the house
  style also uses `INPUTS`, `OPTIONS`, `OUTPUT V`, `OUTPUT REJ`. **Being stricter
  than the house style reports correct files as defective**, which is worse than
  not checking.
- New audit **check 14** (documented error ids must exist in code) first reported
  **80 findings, nearly all false** -- most WARNING ids are raised through
  `lsqvce.util.log.warn`, which an `error()`-only regex cannot see. A check that
  cries wolf 80 times trains you to ignore it. Corrected -> 4 findings, 3 of them
  real.

FIXED THIS PASS
- **ERRORS sections added** to the functions a user actually meets:
  `makeLoader` (14 ids, grouped into "your declaration is wrong" vs "your
  buildDay is wrong"), `covModel` (19 ids, grouped by cause), `adjustBlock`,
  `adjustVCE`, `snoop`, `separableWhiten`, `varWhiten`.
- **INPUT/OUTPUT blocks** added to `adjust`, `separableWhiten`, `varWhiten`,
  `getfielddef`, `mergeStruct`.
- **`separableWhiten` help cited `lsqvce:separableWhiten:notWholeEpochs`; the
  code throws `lsqvce:separable:notWholeEpochs`.** I introduced that error in
  this same pass and check 14 caught it minutes later.
- **A THIRD STALE ALIAS CLAIM**: `html/lsqvce.adjustBlock.html` promised
  `lsqvce.recover` "remains as a deprecated function alias that warns once
  (`lsqvce:recover:renamed`)". `recover.m` does not exist. Same defect as the
  `SeparableProblem` one in §1.3c -- **that is three; every removed name left a
  false alias promise in the HTML.** Check 14 now guards this class.
- **HANDOVER cleaned**: removed a 390-line session document appended at the end
  that claimed *318/318 tests* and described `SeparableProblem` as current -- a
  stale baseline sitting where a reader looks last. Also fixed a duplicated
  `## 0. WHAT THIS TOOLBOX IS## 0. WHAT THIS TOOLBOX IS` heading.

ROBUSTNESS: audited all 22 `catch` blocks that swallow without rethrowing. **All
legitimate**: `gcvLambda` falls back from the Cholesky to the QZ eigensolver,
`BlockAccumulator.logLocalSummary` skips an unrecoverable system in a
DISPLAY-only path, `adjustVCE` defaults a size probe, and the plot/log wrappers
must never break a solve. No silent failure found in a computational path.

NEW AUDIT CHECK 14 -- documented error ids must exist in the code. An id in help
text is a PROMISE: users grep for it and write `try/catch` against it. Verified
both ways (fires on an invented id, silent on the clean tree).

### 1.3b TASK 5 IN PROGRESS (2026-07-19) — master PDF, ReportLab

**NOT FINISHED. NO `.m` FILE TOUCHED THIS ROUND — expect 497/498 unchanged.**

DECISIONS (user): ONE master PDF **replacing** all three manuals; workflow
diagrams as ReportLab boxes-and-arrows from a REUSABLE helper; equations
**DRAWN**, not text-with-unicode.

BUILT AND VERIFIED SO FAR
- `docs_src/pdf_graphics.py` — the reusable layer. `eq()` draws a display
  equation from a TOKEN LIST (deliberately not a LaTeX parser: ~20 display
  equations, and a half-implemented parser mis-renders silently). `workflow()`
  draws boxes-and-arrows from a declarative row spec with an optional feedback
  loop, so every entry point is drawn by ONE routine and they cannot drift.
  Palette matches the HTML doc set.
- `docs_src/build_master_pdf.py` — styles, page furniture, `code()`,
  `callout()`, `table()`, `caption()`.
- `docs_src/master_content_a.py` — front matter + full TOC, Part 1
  (orientation, engine-choice diagram), Part 2 (theory: eq 2.1-2.12 covering
  BLUE, QR, redundancy/MDB, VCE, block-arrow Schur, AR whitening, the separable
  model).
- **RENDERING VERIFIED BY RASTERISING AND LOOKING**, not by "no exception":
  `pdftoppm` -> `view`. That caught a real defect — the hat in `x̂` and `σ̂₀²`
  was trailing the glyph ("x ˆ") because a naive token list appends the accent
  as a normal character. Added a `('hat','')` token that draws the accent
  CENTRED over the previous glyph and advances zero width; needed `prev_w`
  tracking in every advancing branch. Re-rasterised and confirmed.

### 1.3e FOURTH REVIEW (2026-07-19) — a SILENT font gap, and Part 5/6 expanded

**THE MISSING OPERATOR WAS A FONT GAP, NOT A TYPO — AND IT WAS INVISIBLE TO
EVERY TEST SO FAR.** The user reported eq 2.12 missing the operator between R
and I. The `('g','otimes')` token was present all along: **Liberation Serif does
not contain U+2297**, and a glyph a font LACKS does not draw a box -- it draws
NOTHING and takes the rest of the string with it. Eq 2.12 lost the operator AND
the identity matrix after it, through builds that correctly reported ZERO boxes.
A missing glyph is strictly more dangerous than a .notdef box, because the box
is visible and countable.
- Switched the equation faces to **FreeSerif**, which has full maths coverage.
- Added **`_assert_coverage()`**, run at import: it reads the font's cmap and
  RAISES if any glyph in GREEK/SYM_RAW is absent. VERIFIED BOTH WAYS -- pointed
  back at Liberation it fires and names **five** missing glyphs
  (U+2208 U+221D U+2297 U+2308 U+2309), so the ceiling brackets in eq 2.13 had
  been vanishing too, unnoticed.
- Eq 2.12 also gained the transpose it was missing: C = sigma^2 D (R_t kron I_d) D^T.

**SNOOP'S SCOPE, VERIFIED IN CODE.** The user asked whether data snooping is part
of the 3.1-3.3 chain and whether it is flat-path only. Both answers confirmed by
reading `snoop.m` and every call site: it is **a separate step, flat path
ONLY** -- `snoop(sol, D, batchFcn)` takes what `adjust` returns, no engine calls
it, and `adjustBlock`/`adjustVCE` have NO equivalent. Section 3.4 now says this
in a warn box, gains its own figure (an iteration around `adjust`, drawn with the
`loop=` argument), shows the two-call sequence, and explains what to do on the
block path instead.

**FIGURE 4.1 -- buildDay AND epochCov, not OR.** They were drawn side by side in
one row, which reads as alternatives. Now a single box, with the caption stating
they are supplied TOGETHER (buildDay gives design + observations, epochCov gives
the per-epoch D for the SAME rows) and that a scalar stream is the one case where
epochCov is omitted, because at d = 1 there is no D.

**PART 5 EXPANDED WITH THE THEORY** (user request), now with six equations:
- 5.1 why Level 0 is an IDENTITY not a test (trace of a projector = its rank,
  eq 5.1), and what Level 1 actually tests (chi-square, eq 5.2) -- including
  that a ratio BELOW one is also a failure, the one people ignore.
- 5.2 Lie 1 now derives the mechanism: r_i = 1 - h_i (eq 5.4), so an observation
  with leverage -> 1 determines its own fitted value and cannot disagree with
  itself. Lie 2 shows tau's denominator (eq 5.5) containing the very blunders it
  is meant to catch, which is ALSO why rejection is one at a time.
- 5.3 quantifies the blind spot with the MDB (eq 5.6): halving r_i does not
  halve what you can detect, it multiplies it by sqrt(2).

**PART 6 GAINED INSTALLATION** as new §6.1 (old 6.1-6.4 renumbered 6.2-6.5):
unpack -> `lsqvce_setup` -> `run_all_lsqvce` (expect 497/498) -> optional
`savepath`; a requirements table; the note that nothing writes outside the
toolbox folder; and a first-steps table.

Manual now **20 pages, 33 sections, 0 boxes**.

### 1.3d TASK 5 COMPLETE (2026-07-19) — master manual shipped, old PDFs retired

**`lsqvce_manual.pdf` — 18 pages, 0 boxes, 32 numbered sections, 22 references.**
One `.m` file changed all round (`covModel.m` help text); expect 497/498.

WRITTEN THIS ROUND (`docs_src/master_content_c.py`)
- **Part 5, interpretation** — the five verdict levels and why they are read
  top-down; the two lies (small residual at low r_i; tau blinded by the blunders
  that inflate s_0); the residual-vs-redundancy quadrant; what to do next; and
  the three habits for a gravity solution. CONSOLIDATED from the former
  interpretation guide, not rewritten — that material was already correct and
  re-authoring it would only introduce drift.
- **Part 6, reference** — symbols (from the old reference manual, plus D/R_t/d
  for the separable model), an equation index, the error identifiers a user
  actually meets, and the testing section.
- **22 references**, harvested from `build_manuals.py`'s REFERENCES body and
  extended with Kvas (2019).

RETIRED
- `lsqvce_user_manual.pdf`, `lsqvce_reference_manual.pdf` and
  `lsqvce_interpretation_guide.pdf` **DELETED** from the toolbox root.
- `docs_src/build_pdfs.py` (the HTML->LaTeX->pdflatex chain) **GUARDED and
  refuses to run**, same pattern as `build_manuals.py`. Kept rather than deleted
  because its conversion is the only record of how the old equations were marked
  up.
- The **HTML manuals stay**: only the PDFs were superseded. Their three "(PDF)"
  links were removed from `html/index.html`, which now points at the master.
- `README.md` and `docs_src/README.md` rewritten for the new pipeline. The
  latter now states plainly that there are TWO doc paths and that
  `build_manuals.py` is guarded.
- VERIFIED: no surviving reference anywhere to a deleted PDF; one `.pdf` in the
  root.

DOC-SET STATE, so it is not rediscovered:
| path | builds | status |
|---|---|---|
| `build_master_pdf.py` + `master_content_{a,b,c}.py` | `lsqvce_manual.pdf` | **live, THE manual** |
| `build_manuals.py` | the HTML manuals | **guarded** — behind its output (§1.3b) |
| `build_pdfs.py` | the three old PDFs | **retired** — output deleted |

### 1.3c THIRD REVIEW (2026-07-19) — arrow routing, section 1.3, 'separable'

**`'separable'` IS NOT A LEFTOVER — VALIDATED, NO CODE CHANGE NEEDED.**
The user asked whether Part 4.5's `'separable'` example was stale after the
VCE renaming. It is not: `'separable'` is one of **EIGHT covariance-model TYPE
values** in `covModel.m`'s switch (`diag`, `acf`, `psd`, `arma`, `var`,
`separable`, `custom`, `toeplitz`), with its own validation and error ids. It
names MATHEMATICS -- the covariance separates into a per-epoch spatial factor
and a temporal law, C = D (R_t kron I_d) D -- not an engine. Categorised every
occurrence: 77 prose (the model), 20 `separableWhiten`, 12 type-string/case, 2
historical notes about the class rename. The CLASS `SeparableProblem` became
`adjustVCE`; the TYPE STRING was deliberately left alone.
Since the ambiguity is a fair one, it is now stated explicitly in THREE places:
`covModel.m` help, `html/lsqvce.adjustVCE.html`, and manual section 4.5.

**A REAL STALE CLAIM FOUND WHILE CHECKING** (`html/lsqvce.adjustVCE.html`):
- it promised the old name "remains as a deprecated alias that warns once
  (`lsqvce:SeparableProblem:renamed`) and forwards here" -- **that alias does
  not exist**; it was removed 2026-07-18 under §4a. A user following the page
  gets an undefined function. Corrected to say REMOVED, pointing at the
  tombstone page.
- its purpose line still used **SHAPE / SCALE**, the invented stage vocabulary
  the user rejected on 2026-07-18. Reworded. (Audit check 11 did not catch it
  because it guards `.m` user-visible strings, not generated HTML -- the same
  blind spot that produced check 13.)

**FIGURE 1.1 ARROW OVERLAP FIXED.** The `adjustVCE -> Results` edge spans TWO
rows, and `decision_tree` bent every elbow at the vertical midpoint -- which for
a two-row edge lands ON the intervening row, drawing the line straight through
`adjustBlock`. Multi-row edges now leave the source sideways and route down an
outside gutter (reserved in the layout, else they clip at the frame).
VERIFIED GEOMETRICALLY, both ways: sample each line segment and test whether any
sample falls strictly inside a node rectangle -- **2 crossings with the old
routing, 0 with the new**.

**SECTION 1.3 REMOVED**, per the user: a section holding a single callout does
not earn a heading. The "a residual reveals only the fraction r_i" callout moved
into **2.3 Residuals, redundancy, reliability**, where the redundancy number is
actually defined, so it now lands next to eq 2.6 instead of ahead of it. TOC
updated; manual 14 -> 13 pages.

### HOW THE MASTER MANUAL IS BUILT (read before touching the PDF)

**ONE COMMAND, from the toolbox root:**
```
python3 docs_src/build_master_pdf.py      ->  lsqvce_manual.pdf
```
Build-time only; MATLAB never calls it and nothing here ships as a runtime
dependency. **Full instructions are in `docs_src/README_manual.md`** -- token
tables for equations, node kinds for diagrams, and the two traps below.

| file | role |
|---|---|
| `docs_src/build_master_pdf.py` | entry point: styles, page furniture, `P/B/N/code/callout/table/caption`, `build()`, `main()` |
| `docs_src/pdf_graphics.py` | drawing layer: `eq()`, `workflow()`, `decision_tree()`, `rule()` |
| `docs_src/master_content_a.py` | front matter, Part 1 orientation, Part 2 theory |
| `docs_src/master_content_b.py` | Part 3 entry points, Part 4 template guide |
| `docs_src/master_content_c.py` | Part 5 interpretation, Part 6 reference, references — NOT YET WRITTEN |
| `docs_src/README_manual.md` | the build instructions |

`main()` imports the content modules INSIDE the function: they import names from
`build_master_pdf`, so a top-level import is circular. Parts 5-6 are wrapped in
`try/ImportError` so the build works before `master_content_c.py` exists.

**TRAP 1 -- REPORTLAB'S GRAPHICS LAYER IGNORES THE Symbol FONT.**
`String(..., fontName='Symbol')` inside a `Drawing` emits `.notdef` for EVERY
glyph, while the identical call on a CANVAS renders fine. Equations must be
Drawings (they have to be flowables so they page-break), so Symbol is unusable
here. Every glyph is therefore a TRUE UNICODE CODEPOINT drawn in a registered
Liberation TrueType face -- registered in BOTH `pdf_graphics.py` (equations) and
`build_master_pdf.py` (body text; entities like `&#8968;` in a Paragraph fail
identically). `eq()` RAISES if the faces are missing rather than falling back to
base-14 and emitting boxes.

**TRAP 2 -- VERIFY BY EXTRACTING TEXT, NEVER BY LOOKING AT A RENDER.**
```
pdftotext lsqvce_manual.pdf - | python3 -c "import sys;t=sys.stdin.read();print(sum(t.count(c) for c in '\u25a0\u25a1\ufffd'))"
```
**Must print 0.** A `.notdef` box is a real character in the content stream, so
it is countable. A pixel heuristic on a rasterised page CANNOT distinguish a box
from a table border, a rule or an arrowhead, and it silently measures the wrong
page the moment the document repaginates. Both mistakes were made here and both
reported a clean result on a manual full of boxes.

**PICK THE RIGHT DIAGRAM ROUTINE.** `workflow()` is a PIPELINE: rows stacked,
edges inferred from position. `decision_tree()` is a BRANCHING CHOICE with edges
DECLARED explicitly. Figure 1.1 was drawn three times wrongly with `workflow()`,
because with two nodes above and two below it zips them in order, so the "no"
branch silently acquired an arrow to whatever sat beneath it. **Any diamond with
two outcomes needs `decision_tree()`.**

### FIGURE 1.1, CORRECTED A THIRD TIME (2026-07-19)
The user had to state the logic three times before it was drawn right. Final,
and it is now a `decision_tree` with LABELLED yes/no edges:
- noise model known? **no** -> `adjustVCE` (regardless of size: estimating the
  model is inherently iterative)
- **yes** -> does everything fit in memory? **yes** -> `adjust`; **no** ->
  `adjustBlock`
A known noise model is a PREREQUISITE for both `adjust` and `adjustBlock`; size
is the ONLY thing separating those two.

SECOND USER REVIEW (2026-07-19) -- I HAD DECLARED THE BLACK BOXES FIXED. THEY
WERE NOT, AND MY VERIFICATION WAS THE PROBLEM.

**THE VERIFICATION FAILURE, first, because it matters more than the bug.** I
"proved" the fix with a pixel heuristic counting columns containing a long
unbroken dark run, and reported 4-before / 0-after. Both numbers were worthless:
(a) adding section 2.7 had REPAGINATED the document, so the page I measured
"after" was a nearly blank page-4, not the equations; (b) the heuristic counts
ANY solid dark feature -- table borders, the orange rules, filled arrowheads --
so the title page alone scores 72 with no equation on it. I shipped a green
result from a test that could not see the defect. **THE RELIABLE TEST IS TEXT
EXTRACTION**: `pdftotext | count of U+25A0` -- a .notdef box appears as a real
character in the content stream. 15 before, 0 after, and 52 Greek/math glyphs
now extract correctly. Pixel heuristics on rendered pages are not evidence.

**THE ACTUAL BUG: ReportLab's GRAPHICS layer ignores the Symbol font encoding.**
Proved in isolation: `String(..., fontName='Symbol')` inside a Drawing extracts
as black boxes, while the identical characters drawn on a CANVAS render fine.
Every equation here is a Drawing (they must be flowables so they page-break), so
the whole Symbol-font approach was unusable -- it was never going to work, and
my first "fix" only moved which glyphs went through it.
RESOLUTION: register the **Liberation TrueType faces** (metric-compatible with
Times, full Greek + true minus U+2212 + radical + Kronecker + ceiling brackets)
and write every glyph as its TRUE UNICODE CODEPOINT. GREEK/SYM_RAW rewritten
from Symbol-font byte hacks to real codepoints. The same faces are now
registered for BODY TEXT too -- entities like `&#8968;` and `&#8321;` in a
Paragraph were failing identically, which is where the last 4 boxes lived.
The old "reject non-ASCII" guard is inverted: `eq()` now RAISES if the Unicode
faces failed to register, because the base-14 fallback cannot draw Greek and
silently emitting boxes is the failure being designed out.

FIGURE 1.1, CORRECTED AGAIN. The user pointed out that a known noise model is a
**PREREQUISITE** for `adjust`, not a parallel choice: you hand it the covariance
matrix, so it must exist in memory before the call. The decision order is now
(1) is the noise model known? -> if not, `adjustVCE`; (2) only then, does
everything fit in memory? -> `adjust` if yes, `adjustBlock` if no. Table 1.1
updated to match ("required -- you supply the covariance").

USER REVIEW OF THE FIRST DRAFT (2026-07-19) -- three defects, all fixed
1. **BLACK BOXES throughout Part 2.** Every non-ASCII glyph in a Times/Helvetica
   token rendered as .notdef. MEASURED, not guessed: 'sigma' drawn from Symbol
   gives 4940 dark pixels, the same glyph as U+03C3 from Times gives 538 (the box
   outline). The base-14 text fonts are StandardEncoding. Fixed by routing minus
   (U+2212), radical (U+221A), tilde (U+0303), ceiling brackets and the ellipsis
   through the Symbol font, and **adding a GUARD in `eq()` that RAISES on any
   non-ASCII in an 'r'/'v'/'b' token**. The guard immediately caught a `\u221a`
   I had missed inside eq 2.6, hidden behind dead scaffolding
   (`'sqrt' if False else 'r'`) from the first draft.
   VERIFIED BY PIXEL ANALYSIS, both ways: a .notdef box is a SOLID dark
   rectangle, so count columns with a >14px unbroken dark run -- **4 before the
   fix, 0 after**. That is objective; squinting at a thumbnail is not.
2. **Figure 1.1 asked the wrong question.** It branched on "does it fit in
   memory?", which the user correctly called superfluous: if the normal matrix
   did not fit, no computation would be possible at all. The real distinction is
   whether **A, l AND the full covariance** fit at once (-> `adjust` takes them
   directly) or the problem must be **accumulated in batches** (-> `adjustBlock`
   / `adjustVCE`). Figure and Table 1.1 both rewritten; the noise-model question
   now correctly applies only to the streaming branch.
3. **Section 2.7 did not say how D is stored.** Now documented FROM THE CODE
   (`covModel.m:29`, `epochSlice`), not from assumption -- and the code accepts
   TWO shapes, where the draft implied one:
     `nEpoch x d`      diagonal, per-axis sigmas
     `nEpoch x d x d`  full per-epoch covariance (the realistic 3-D PPP case)
   Explicitly NOT `n x 9` (flattening loses the axis pairing; raises
   `lsqvce:covModel:badCep`) and NOT a large sparse `3n x 3n` (D is
   block-diagonal by epoch; only blocks are stored). Added eq 2.13 giving the
   indexing rule that ties Cep to epoch-major rows,
   `epoch(k) = ceil(k/d)` and `rows(e) = d(e-1)+1 ... de`, which is exactly what
   `epochSlice` computes -- so the row-order hazard is now stated as an
   indexing consequence rather than a bare warning.

PARTS 3 AND 4 WRITTEN (2026-07-19, `master_content_b.py`)
- Part 3: `adjust`, `adjustBlock`, `adjustVCE` (with the feedback loop drawn via
  `workflow(..., loop=)`), `snoop`, and reading a `Results` — each with its own
  diagram and the calls that go with it.
- Part 4, the template guide, seven steps: procOpt -> observationGroups ->
  buildDay -> a-priori covariance -> epochCov -> run/read -> hooks, closing with
  "four ways to get it wrong" (wrong-then-right code per failure) and the
  five-step checklist for a stream that will not whiten. Cross-references the
  theory by equation number (2.10 for merging systems, 2.12/2.13 for row order).

STILL TO WRITE (next round)
- Part 5 (interpretation) and Part 6 (symbols, equation summary, error ids,
  testing), then the 21 references — HARVEST from `docs_src/build_manuals.py`
  (`INTERP`, `REFMAN`, `REFERENCES` bodies) rather than rewriting.
- Retire `lsqvce_user_manual.pdf`, `lsqvce_reference_manual.pdf`,
  `lsqvce_interpretation_guide.pdf` and their `build_pdfs.py` (the pdflatex
  chain) once the master covers them; update README, index.html, helptoc.xml.

**HAZARD FOUND — READ BEFORE TOUCHING THE MANUAL HTML AGAIN.**
`docs_src/build_manuals.py` is the GENERATOR for
`html/lsqvce_{user_manual,reference_manual,interpretation_guide,choosing}.html`,
and it still contains everything fixed downstream this session:
- Parts IV-B/IV-C (task 2, §1.2a) are **ABSENT from the source** — re-running the
  generator would DELETE them.
- The 4 dead fluent names (`vceGroupsBy`/`.shape(`/`.scale(`) and 2 stale `nc`
  examples are **STILL IN THE SOURCE** — re-running would REINTRODUCE them, and
  audit check 13 would then fire on the regenerated HTML.
Those edits were applied to the OUTPUT, not the INPUT. Nobody has re-run the
generator, so the shipped HTML is correct — but it is one command from
reverting. RESOLUTION when the master PDF lands: retire `build_manuals.py`'s
manual bodies (the master becomes the single source) or back-port. Until then,
**do not run `python3 docs_src/build_manuals.py`.**

### 1.3a RENAMES + TASK 3 (2026-07-19) — separable→grace, runners merged

**EXPECT 497/498 UNCHANGED.** Renames and a delegation; no test logic touched.

RENAMES (user: "separable" is a leftover from `SeparableProblem`, now `adjustVCE`)
| was | is now |
|---|---|
| `examples/demo_separable_data.m` | `examples/demo_grace_dataLoader.m` |
| `examples/demo_separable.m` | `examples/demo_grace.m` |
| `tests/test_separable.m` | `tests/test_adjustVCE.m` (it tests that class) |
| `makeProblem` (in that test) | `makeVCEBuilder` |

- 13 files updated, 0 leftover identifiers. **The 61 TECHNICAL uses of
  "separable" were left alone** — the separable covariance model
  `C = sigma^2 D (R_t kron I_d) D` is real vocabulary, not a leftover. Only the
  `demo_separable*` / `test_separable` IDENTIFIERS moved.
- Longest-name-first substitution order, else `demo_separable` eats the prefix of
  `demo_separable_data` and yields `demo_grace_data`.
- HANDOVER.md deliberately NOT rewritten: its history says what was true then.

**WHAT `makeProblem` ACTUALLY DID** (the user asked; it was undocumented): it does
NOT build loaders — `demo_grace_dataLoader` does. It wraps those loaders in a
configured `lsqvce.adjustVCE` with verbosity off, i.e. returns a BUILDER that has
not run. Renamed `makeVCEBuilder` and given help saying so, including its cost.
NOTE it collided with an UNRELATED `makeProblem(n,u,seed,outliers)` in
`test_lsqvce.m` (returns `[A,l,p,xt]`, a flat synthetic set) — same name, two
contracts, two files. The rename removes the collision; the other one stays.

TASK 3 — RUNNERS MERGED
`run_lsqvce('all')` now **delegates to `run_all_lsqvce`** instead of running its
own `runtests(fullfile(root,'tests'))`. That old line was the §1.7b trap itself:
runtests on a FOLDER silently skips a file that fails to load and reports the
rest green. One gate, one implementation.
- `'all'` with a `'Name'` filter is now an ERROR
  (`lsqvce:runLsqvce:allTakesNoFilter`) rather than a filtered pseudo-gate.
- RETURN TYPE DIFFERS BY TARGET and is documented: `'all'` returns
  `run_all_lsqvce`'s per-suite struct array; everything else returns
  `matlab.unittest.TestResult`. They answer different questions (one suite vs
  one test), so they summarise at different granularity.
- Shortcut `'separable'` -> `'adjustvce'`. I first kept `'separable'` as a silent
  alias for muscle memory; **the user rejected it — the no-alias rule (§0a rule
  8) wins.** An unknown name already prints the shortcut list, so the old
  spelling self-documents.
- `run_all_lsqvce` help gained the `.secs` field (added §1.1e, never documented —
  a §4g surface I missed) and a cross-reference to `run_lsqvce`.

**AUDIT CHECK 11 CAUGHT ME MID-EDIT.** I wrote `run_lsqvce('all')  THE GATE` into
a `fprintf`; `GATE` is one of the INVENTED stage names the user rejected on
2026-07-18, and check 11 guards user-visible text against exactly that. Reworded
to "the whole suite". The check was right and I was careless — worth recording as
evidence the guard earns its keep on the person who wrote it.

STALE README FIXED (found while updating runner docs): claimed "16 files, ~420
tests" (it is 18 / 498), and a note recommended
`runtests('tests','IncludeSubfolders',true)` for a `tests/` SUBFOLDER that was
flattened away in §1.7b — i.e. it recommended the unsafe form for a directory
that no longer exists.

### 1.2b RUN CONFIRMED + THE TIMING INSTRUMENTATION PAID OFF IMMEDIATELY

**MATLAB 16:28 — 497/498, 0 failed, 0 load failures, 0 `Verification failed`,
12 warning categories unchanged.** Warnings 56 -> 54 (-2), exactly the two shed
by folding `test_history_carries_whiteAlpha` into `test_summary_is_per_stream`.
Task 2 (§1.2a) and task 1 are both CONFIRMED. Log 107,976 bytes.

**FIRST RUN WITH PER-SUITE TIMING, AND IT ANSWERED THE USER'S QUESTION FROM TWO
ROUNDS AGO IN ONE LINE:**
```
total  : 127.6 s wall clock
  test_separable  : 44/44 pass ... 108.9s     <-- 85% of the entire suite
  every other suite:                0.3s - 2.4s
```
The user noticed a slowdown at §1.1e; the log could not say why, so the runner
was instrumented. The very first instrumented run localised it. **RECORDED AS
TASK 6** — cause measured (`makeProblem` rebuilds the whole 4-day x 2-arc x
300-epoch GRACE demo dataset in each of 44 tests; 39 then run a full VCE loop),
but NOT FIXED THIS ROUND: it touches all 44 tests and belongs in its own round
with a before/after comparison, not bundled into a documentation wrap-up.

BASELINE FOR ALL FUTURE RUNS: **497/498, 127.6 s, 54 warnings / 12 categories,
107,976 bytes.** Compare per-suite seconds from now on, not just counts.

### 1.2a TASK 2 (2026-07-19) — interface guide, as user-manual Parts IV-B/IV-C

**NO `.m` FILE CHANGED. EXPECT 497/498, exactly as §1.1e predicted** (that
prediction stands unaltered by this round: the test-count reduction there was
deliberate).

THE GAP, MEASURED BEFORE WRITING
`buildDay` appeared in exactly ONE place in the whole doc set —
`html/lsqvce.makeLoader.html`. Grep for it in the user manual, the choosing
guide and the interpretation guide: **zero hits in all three.** Worse, the
manual's "Part II — What you supply" documents `batchFcn`, the FLAT seam, so a
gravity user reading the part named for their question was being shown the wrong
function. Part IV covered block-arrow THEORY but not its interface.

BUILT (user decision: in the user manual, everything a user writes, worked
wrong-then-right examples)
- **Part IV-B "The interface you write"** — the four things a user writes
  (buildDay, the loader declaration, epochCov, hooks), each with what the
  reference page deliberately does not say: that you return a WHOLE DAY and the
  accumulator slices it; why the 5-arg arity is preferable (it is the only one
  that receives `procOpt` and lets `nc` be derived); the merge-systems rule; the
  single-char `systemId` constraint; that hooks attach AT CONSTRUCTION and a
  post-hoc bare handle is skipped SILENTLY.
- **Part IV-C "Four ways to get it wrong"** — wrong-then-right code per failure.
  I: non-deterministic builder. II: axis-major row order (the silent one — the
  model strides by `nAxis`, so a 60-epoch correlation length reads as 180 and
  NOTHING errors). III: variances passed as weights. IV: mismatched row counts
  (caught by name, `lsqvce:makeLoader:dayShape`). Closes with a five-step
  checklist for a stream that will not whiten.
- TOC updated; all `href="#..."` anchors verified to resolve.
- CLAIMS VERIFIED IN CODE, not asserted: epoch-major ordering
  (`covModel.m:32`), both error ids exist (`makeLoader.m`), parameter names are
  case-insensitive so only the METHOD names were stale (`mergeNV`, ~line 2144).

**FOURTEEN MORE DEAD NAMES IN SHIPPED DOCS — found while surveying, not sought.**
`P.vceGroupsBy(...).shape(...).scale(...)` — ALL THREE renamed to
`groupComponentsBy`/`stochasticModel`/`varianceComponents` — still appeared in
`lsqvce.adjustVCE.html` (11), `lsqvce_user_manual.html` (6) and
`lsqvce_choosing.html` (3), INCLUDING the manual's headline Part IV example and
the adjustVCE reference page's own syntax line. A user copying any of them gets
`Unrecognized method`. 20 replacements; grep now returns 0.

**AUDIT CHECK 13 BUILT** (`check_html_examples`), the candidate parked in §1.1d.
Scans `<pre>` blocks only — prose may legitimately DISCUSS an old name, but an
example is what gets copied. Keyed old→new so the finding says what to write
instead. Covers the fluent renames, the removed leading `nc`, `P.convergence(`
and `runtests(`. VERIFIED BOTH WAYS (§1.0d): fires with file + term + replacement
on a reintroduced `vceGroupsBy` (exit 1), silent on the clean tree (exit 0).
**IT IMMEDIATELY EARNED ITS KEEP: it found a 15th stale reference I had missed**
— `html/index.html` still recommended `runtests('tests')`, the exact trap §1.7b
removed (it skips files that fail to load and reports the rest as green). Now
replaced with `run_all_lsqvce` plus a warn box explaining why.

WHY THIS KEEPS HAPPENING (worth stating once): checks 3 and 11 guard `.m` files
and user-visible STRINGS. Nothing guarded generated HTML, so every rename since
the doc set was generated has been able to rot silently. That is now closed for
KNOWN removals; a rename made in future must still be added to check 13's table.

### 1.1e TASK 1 ROUND 4 (2026-07-19) — run verified, timing instrumented

**MATLAB RUN 16:02 CONFIRMS ROUNDS 1-3: 498/499, 0 failed, 0 load failures,
0 `Verification failed`, `doc lsqvce` works (user-verified, incl. the new
"Your own case (template)" entry).** Task 1 CLOSED.

**MY PREDICTION OF 499/500 WAS ARITHMETIC ERROR, NOT A MISSING TEST.** All four
new/changed tests ran (confirmed by name in the log). 495 + 3 = 498; round 2's
"expect 498/499" was right and I mis-added when writing round 3. NOTE for
counting: `grep -c '^function test'` OVERCOUNTS BY ONE -- it matches the harness
entry `function tests = test_separable`, which is not a test.

THE USER NOTICED THE RUN TOOK LONGER. It did, and the log could not say why.
- Warnings 51 -> 56 (+5), log 103,684 -> 110,692 bytes. **No new warning
  category -- still 12.** The +5 is entirely inside two existing ones ("VCE did
  not converge" 13->16, "the loop did not converge" 13->15).
- ATTRIBUTED BY PARSING THE LOG: `test_summary_is_per_stream` +3 and
  `test_history_carries_whiteAlpha` +2 = exactly the +5. Every other test's
  warning count is unchanged. Benign -- both use deliberately tight budgets
  (`MaxOuter` 1-2, `MaxIter` 3-4) and assert on SHAPE, so the early-stop warning
  is expected. But each ran a FULL VCE LOOP.
- **FIXED: `test_history_carries_whiteAlpha` DELETED as a standalone test**; its
  two assertions folded into `test_summary_is_per_stream`, which already has a
  completed run in hand. A second full VCE loop to check one field is not worth
  the seconds. `test_separable` 45 -> 44, so **EXPECT 497/498 next run** (a
  DELIBERATE reduction, not a regression).
- **FIXED THE UNDERLYING GAP: `run_all_lsqvce` NOW RECORDS TIME** -- per-suite
  seconds on both the console and the log line, plus a `total : N s wall clock`
  header. A test that quietly starts running a full VCE loop costs seconds and
  shows up NOWHERE in a pass count; that is precisely what happened here and the
  log could not answer it. From the next run on, compare per-suite seconds.

BUG INTRODUCED AND CAUGHT IN THE SAME ROUND (recorded, since it would have been
severe): adding `secs` to the two `results(end+1) = struct(...)` sites without
adding it to the EMPTY INITIALISER at line 54 would raise
`MATLAB:heterogeneousStrucAssignment` on the first suite and break the ENTIRE
RUNNER -- a total failure from a "cosmetic" logging change. Caught by re-deriving
all three field lists and checking they agree. §9's block-spec check (audit 9)
exists for exactly this class in `+lsqvce`; the runner has no such guard.

### 1.1d TASK 1 BUILD, ROUND 3 (2026-07-19) — registration, demo cleanup, stale docs

**AWAITING MATLAB. EXPECT 499/500** (round 2 expected 498/499; this round adds
`testAdjustBlockResultHasNoIterationHistory`, so `test_adjustBlock` 21 -> 22).

**DO adjust / adjustBlock NEED ADAPTING FOR convergence()? NO — verified, and
the answer differs per engine.**
- `lsqvce.adjust` returns `[sol, D, V]` STRUCTS, not a Results. It is the flat
  in-memory shortcut with no iteration, so the readers are not reachable from it
  and should not be.
- `lsqvce.adjustBlock` DOES return a Results (line 128) but never calls
  `setIterative`, so `history` stays `struct([])` -> `requireHistory`'s
  `isempty(h)` fires and both readers raise `lsqvce:Results:noHistory` with a
  message naming the engine that DOES populate a history. Correct by
  construction, no change needed.
- GAP THIS EXPOSED, now closed: round 2 proved that guard only on a HAND-BUILT
  fixture. `testAdjustBlockResultHasNoIterationHistory` now proves it on a REAL
  `adjustBlock` result — the case a user actually hits.

DEMO CLEANUP (§1.1a decisions 10-11)
- **Decision 11 DONE**: the redundant `'AcvTol', 1e-3` override is gone (it
  equalled the default). The three GENUINE overrides now carry a comment saying
  they are DEMO-SPEED choices, naming the defaults they replace.
- **Decision 10 REVERSED — `broadcastHooks` is NOT redundant, and removing it
  would have SILENTLY DISABLED the demo's plots.** The reasoning behind decision
  10 was right about `makeLoader` (a bare handle IS broadcast) and wrong about
  WHERE the demo assigns: it sets `L{i}.plotPostSolve` AFTER construction, and
  `adjustVCE.plotHooks` (~line 1885) reads that field with
  `if ~isa(hookMap,'containers.Map') ... continue` — a bare handle assigned
  post-construction is skipped WITHOUT ERROR. Caught before shipping by reading
  the consumer instead of trusting the decision. The helper stays, now with a
  comment recording exactly why. LESSON: normalisation that happens at
  CONSTRUCTION does not apply to a field assigned afterwards.
- The demo now calls `R.summary()` and `R.convergence()` instead of hand-rolling
  the per-stream print from `R.history` — it should demonstrate the toolbox, not
  work around it.

TEMPLATE REGISTRATION (§0a rule 6) — task 1 deliverable now complete
- `html/lsqvce_template.html` in the house style (warn/danger/note callouts for
  the two-registers point, determinism, epoch-major rows, the immutable-procOpt
  rule, the five fed-not-restated values, merging systems, systemId).
- Registered in `helptoc.xml` (top level, beside the manuals), `html/index.html`
  (a note by the "which front door" box + a row in the docs table), `README.md`
  (two tables), `+lsqvce/Contents.m` (its own "Your own case" block).
- `README.md` line 238 and `Contents.m` both called `demo_separable_data` "THE
  template to adapt" — that claim now belongs to the template; the demo is
  described as the WORKED case to read alongside it.

STALE DOCS FOUND AND FIXED (not sought — hit while registering)
- **TWELVE HTML examples across SEVEN pages still showed the REMOVED leading
  `nc` argument** (`adjustBlock(nc, systems, days)`, `adjustBlock(6, {L}, 1:3)`,
  `adjustVCE(nc, L, days)`). The §4a removal (2026-07-18) updated the `.m` help
  but not the generated HTML — §4g surface 3, the one that gets forgotten. These
  are SHIPPED, CLICKABLE docs telling users to call a signature that now errors.
  Pages: `lsqvce.Results`, `lsqvce.adjustBlock`, `lsqvce.adjustVCE`,
  `lsqvce.makeLoader`, `lsqvce.util.log`, `lsqvce_choosing`,
  `lsqvce_user_manual`. All 12 fixed; grep now returns 0.
- `Contents.m` line 139 still said "the GRACE-FO drivers live in
  examples/gravity" — that folder was removed in §1.3b.

WORTH A FUTURE AUDIT CHECK (not built, to keep this round reviewable): the
entry-point signature check (check 3) covers `.m` call sites but NOT `<pre>`
blocks in `html/`. A check that parsed shipped HTML examples for known-removed
signatures would have caught all twelve. Candidate for check 13.

### 1.1c TASK 1 BUILD, ROUND 2 (2026-07-19) — convergence() moved, summary() added

**AWAITING MATLAB.** Baseline before this round: 495/496, 0 failed
(MATLAB-confirmed 14:15, log 103,684 bytes, 51 warning lines / 12 categories,
0 load failures). **EXPECT 498/499**: `test_separable` 43 -> 46 (+3 net; one
existing test gained a second assertion, two tests are new).

ROUND 1 IS CONFIRMED GREEN. `lsqvce_setup` now runs warning-free (user-verified),
495/496 unchanged, warning inventory byte-identical, log +143 bytes (timestamp
noise). The three dead `addpath` targets and the template landed without effect
on the suite, as intended.

DONE — §1.1a decisions 6 and 7
- **`convergence()` MOVED from `adjustVCE` to `Results`.** A TRUE move: deleted
  from `adjustVCE` (lines 682-756), no alias, per §0a rule 8. The private helper
  `plotConvergence` (was 1029-1082) moved with it and is now
  `methods (Access = private, Static)` -- once `WhiteAlpha` came from the
  history it touched NO object state, and a static declaration is the honest
  one. Signature drops the R argument: `R.convergence(...)`, not
  `P.convergence(R, ...)`.
- **`hist.whiteAlpha` added** (`adjustVCE.m` ~665). This is what breaks the
  dependency on the builder: the plot needs the threshold, and a Results must be
  readable without the object that produced it. Guarded in the plot, so a
  history predating the field just gets no threshold line.
- **`summary()` added to `Results`** -- one row per STREAM (stream, sigma,
  arOrder, white, pValue) from the FINAL pass, via the new file-local `finalRow`
  (returns NaN for an absent/mis-shaped field: a diagnostic must still print
  when one optional field is missing).
- **Error id `lsqvce:adjustVCE:noHistory` -> `lsqvce:Results:noHistory`**, raised
  by the new private `requireHistory`, shared by both readers so they cannot
  drift apart.
- **Tests**: `test_convergence_table_shape` updated to `R.convergence()`;
  `test_convergence_errors_without_history` updated to the new id AND extended to
  assert `summary()` raises the SAME id (if they ever diverge, one grew its own
  copy of the guard); NEW `test_summary_is_per_stream` (pins the shape
  difference between the two readers -- the whole reason both exist) and NEW
  `test_history_carries_whiteAlpha`.
- **Docs**: `Results` class help METHODS block lists all five readers plus the
  complementarity note; `html/lsqvce.Results.html` gains the two rows, a
  "The two iteration readers" table, and the `noHistory` note. `convergence` was
  never in the HTML, so nothing had to be retargeted.
- **Stale comment fixed** in `adjustVCE.m`: the `verbose` property said
  "2 stage | 3 diagnostic"; level 3 was removed 2026-07-19.

VERIFIED
- `grep` for `P.convergence`, `adjustVCE:noHistory`, `plotConvergence` across
  `.m`/`.html`/`.md`: the only hits are the new ones inside `Results.m`.
- `audit_tree.py` exit 0, all 12 checks.
- Balance +0 on all four edited files, using the AUDIT'S OWN `code_lines()`.

A THIRD §4h FALSE POSITIVE, and the lesson landed for real this time.
A throwaway depth counter reported `Results.m` unbalanced at the line
`v = double(M(end, :))` -- `end` as an ARRAY INDEX, the exact class §4e names
(`history(end+1)`). Two successive attempts to patch the regex made it worse
(-1, then +1). The shipped `check_balance` was right all along because it tracks
BRACKET DEPTH and skips `end` inside brackets. **RULE: do not hand-roll a
balance checker; import `docs_src/audit_tree.code_lines`.** Retroactively, last
round's template "+3" was the same class (prose in error strings) -- the
validated checker gives +0 there too.

### 1.1b TASK 1 BUILD, ROUND 1 (2026-07-19) — template + setup/audit fixes

**AWAITING MATLAB.** Baseline before this round: 495/496, 0 failed.
Expect the SAME count -- no test was added or changed this round.

BUILT
- **`template/lsqvce_template.m`** (348 lines, 9 functions), the §1.1a
  deliverable. Returns the Loader cell array; the user calls the toolbox. TWO
  REGISTERS, stated in the header so it cannot mislead: `buildLoaders` is real
  runnable code deriving one Loader per `observationGroups` entry, while
  `myBuildDay`/`myEpochCov` ERROR until written. Carries the commented hook
  stub with the full ctx schema, the epoch-major row-order warning, the
  merge-systems rule, the sysID single-char constraint, and the
  "procOpt is IMMUTABLE INPUT, not shared state" rule.
- **`lsqvce_setup.m`**: added `template/`; **removed THREE DEAD `addpath`
  TARGETS** -- `tests/accum`, `examples/accum`, `examples/gravity`, all deleted
  2026-07-18 (§1.7b, §1.3b) but still added here, so every user running setup
  got three `MATLAB:addpath:nameNonexistentOrNotADirectory` warnings. Same class
  as §1.0 defect 4: a removal recorded as done with one surface left behind.
  Also fixed the help header and the printed hint that still said
  `runtests('tests')` -- the very trap §1.7b removed.
- **`audit_tree.py` CHECK 12** (`check_addpath_targets`): every literal
  `addpath(fullfile(root,...))` target must exist. VERIFIED BOTH WAYS -- reports
  `lsqvce_setup.m:25 addpath target 'template' does not exist` on a broken tree
  (exit 1) and is silent on a healthy one (exit 0). It would have caught all
  three dead entries. Computed paths are skipped, not guessed at.
- **`+util/verbosity.m` help**: `resolve` and `atLeast` said levels `0..3` while
  the code enforces `0..2` (level 3 was merged into 2 on 2026-07-19 and is
  REJECTED with its own message). Help only; no behaviour change. Check 11 does
  not cover this because it guards user-visible STRINGS, not help text.

STRUCTURAL NOTE (a §4h false positive, recorded so it is not chased again)
A crude balance counter reports the template +3. All three are PROSE INSIDE
ERROR-MESSAGE STRINGS -- "Asked **for** day", ".nAxis **if** multi-axis",
"**if** every stream is scalar" -- with MATLAB's `''` escapes defeating a naive
string stripper. The SHIPPED `audit_tree.py` check 7 uses a proper `code_lines()`
and passes. There is no defect; the throwaway checker was the broken one. This is
the third time this lesson has been paid for.

NOT DONE THIS ROUND (deliberately, to keep the round reviewable)
- The `convergence()` move + `summary()` (§1.1a decisions 6-7). Next.
- The demo cleanup (`broadcastHooks`, `AcvTol`) (§1.1a decisions 10-11).
- `template/` is NOT yet registered in `README.md`, `Contents.m` or
  `html/index.html`, and has NO html page. §0a rule 6 requires this before the
  task is closed.

### 1.1a TASK 1 SPECIFICATION (2026-07-19) — audit + decisions

_STATUS: decisions 1-5 BUILT (§1.1b, the template); 6-7 BUILT (§1.1c);
11-12 BUILT and 10 REVERSED (§1.1d — removing `broadcastHooks` would have
silently disabled the demo plots). Task 1 is COMPLETE pending the MATLAB run._

Task 1 asked for an audit of `demo_separable_data.m` as the specification for a
template. This section IS that specification. Every open question was put to the
user and answered; the answers are recorded as DECISIONS below with their
reasoning, so none has to be re-litigated. **No code has been written yet.**

**BASELINE for this work: 495/496, 0 failed** (MATLAB 2026-07-19 11:27, log
103,541 bytes, 51 warning lines, 0 files failed to load, 0 `Verification failed`).
The 155 KB log size quoted elsewhere in this file is STALE — the shrink is
removed failure diagnostics, not suppressed output. Compare against 103 KB and
51 warnings from now on.

#### THE AUDIT — what a user adapts in `demo_separable_data.m` (309 lines)

ONE file, but NINE local functions in four tiers. The tiers are the finding: the
problem is not the file count, it is that ~110 lines of synthetic-data
fabrication are INTERLEAVED with the ~60 lines carrying the real contract.

| tier | functions | what the user does |
|---|---|---|
| MUST REPLACE (your data) | `demoBuildDay`, `demoEpochCov` | write from scratch against the contract |
| MUST EDIT (declarations) | `demoCfg`, `demoLoaders`, `whiteCov` | topology, streams, a-priori covariance |
| MAY EDIT (your convention) | `commonParamTable`, `localParamTable`, `demoAxes` | 14-column meanings are the demo's, not the toolbox's |
| SCAFFOLDING (delete) | `A_design`, `B_design`, `colouredNoise`, `arcIndex` | synthetic data only |

Contracts, verified against the code:
- `demoBuildDay` -> `.A` (m x nc), `.B` (m x nl), `.l` (m x 1), `.t` (m x 1),
  `.arcIdx` (1 x nArc cell), optional `.nAxis`. MUST BE DETERMINISTIC — the
  accumulator re-reads every day for the residual pass.
- `demoEpochCov` -> `.Cep` (nEpochDay x d x d) and `.arc(i).epochs/.rows`. Only
  needed for multi-axis streams.

#### THE `cfg` / `procOpt` FINDING (the duplication flagged in §1.OPEN)

They are NOT duplicates by design, but they ARE duplicated in practice.
- **`cfg` has NO toolbox meaning.** VERIFIED: no file in `+lsqvce/` reads it
  (`synthLoader.m`'s `cfg` argument is an unrelated self-contained loader). It is
  closure state captured by the anonymous builder. `demo_separable.m` touches
  exactly two fields (`cfg.days`, `cfg.trueSigma`).
- **`processingOptions` IS the declared interface**: stored as `L.procOpt`, piped
  as `buildDay`'s 5th argument, surfaced to hooks as `ctx.procOpt`. The toolbox
  never reads a FIELD of it — it is an opaque pass-through.
- **THE DEFECT:** the demo uses both at once. `mEpoch`, `dt`, `arcsPerDay`, `seed`
  are copied into `procOpt`, but `demoBuildDay` reads them from `cfg` and its
  `procOpt` argument is marked `%#ok<INUSD>` — GENUINELY UNUSED. Four values
  stated twice, free to drift. Same bug class as the `cfg.nc` duplication killed
  in §4a. `procOpt` wins because it is the interface.

#### THE USER'S REAL `procOpt` (supplied 2026-07-19)

A full processing configuration, not a bag of builder settings: paths, logging,
forces, constants, ephemeris, earthRotation, metadata, observationGroups, plus
LIVE CURSORS (`iteridx`, `indexToSystem`, `indexToGroups`).

FIVE FIELDS DUPLICATE TOOLBOX STATE: `days` (positional arg to `adjustVCE`),
`loggingLevel` (`verbose`), `loggingOverviewFileID` (`logFid`),
`showAndSaveGraphics` (`plots(...)`), `observationGroups` (`observationTypes` +
`groupComponentsBy`). That is the §4a `nc` defect five times over.

**DECISION (user): the config STAYS the origin and FEEDS the toolbox call.** No
field is deleted from the user's struct; no value is typed twice:
```matlab
procOpt = myPipelineConfig(...);                 % your pipeline, unchanged
L = lsqvce.makeLoader(@myBuildDay, ..., 'processingOptions', procOpt);
R = lsqvce.adjustVCE(L, procOpt.days, ...
        'verbose', procOpt.loggingLevel, ...
        'logFid',  procOpt.loggingOverviewFileID) ...
      .groupComponentsBy('obsType').run();
```
THREE OF THE FIVE ARE NOT ONE-TO-ONE and the template must show the mapping, not
a pass-through: `loggingLevel` (user's scale vs the toolbox's 0/1/2 — a wrong
level silently changes output volume); `observationGroups` (a cell of structs ->
one loader per group, a real derivation); `showAndSaveGraphics` (an `if` around
hook attachment, not an option).

**DECISION (user): DROP the live cursors.** `procOpt` is captured at `makeLoader`
time and piped UNCHANGED thereafter, so `iteridx` would freeze at pass 1 while
the loop iterates around it — a builder reading it would silently see the wrong
pass forever. `ctx.pass` / `ctx.sysID` are computed live and correct by
construction. **RULE for the interface guide: `procOpt` is IMMUTABLE INPUT, not
shared state.**

#### DECISIONS — the template

1. **ONE file.** Loaders + `procOpt` only; NO run/report step (the user calls the
   toolbox). Deliverable name/location still to choose.
2. **`buildDay`/`epochCov` are EMPTY STUBS** with the contract documented; the
   user fills everything. The stub ERRORS until filled.
3. **The `observationGroups` loop IS real runnable code** — one loader per group.
   NOTE the file therefore has TWO REGISTERS: a working loop plus a stub that
   errors. The file must say so plainly or the user will not know which half to
   run and which to write.
4. **`procOpt` is a nested struct-of-structs**; `buildDay` reads ONLY its
   `procOpt` argument (this is what kills the `%#ok<INUSD>`).
5. **A COMMENTED HOOK STUB** is included: `ctx` is a contract (schema 1) and a
   user who never sees it does not know the seam exists.

#### DECISIONS — the toolbox

6. **`convergence()` MOVES from `adjustVCE` to `Results`** — a TRUE move, deleted
   from `adjustVCE`, no alias (§0a rule 8). It is a POST-RUN READER, not part of
   the solve: everything it computes comes from `R.history`; the only thing it
   takes from `P` is `WhiteAlpha`, and only to draw a threshold line. Requiring
   the builder object to read a finished result is the awkwardness being removed.
   BLAST RADIUS, verified: 2 test call sites (`test_separable.m:97,120`), 0 demo
   call sites, 0 html pages referencing it. `WhiteAlpha` reaches `Results` by
   being stashed in `hist` at `adjustVCE.m` ~663, beside `converged`/`nOuter`,
   which also makes the history self-describing.
7. **`summary()` added to `Results`** — the PER-STREAM final state (stream, final
   sigma, final AR order, white yes/no, p-value). It COMPLEMENTS `convergence()`
   rather than duplicating it: `convergence()` is per-PASS reduced across streams;
   `summary()` is per-STREAM reduced across passes. Both read `R.history`, which
   already carries `streams`, `sigma`, `arOrder`, `white`, `whitePval`.
8. **NO `runVCE` wrapper**; the fluent chain and its defaults are UNCHANGED.
9. **Hooks stay USER-WRITTEN in the demo.** Shipping default hooks would nudge
   users into attaching canned plots instead of writing their own, defeating the
   seam's purpose.

#### DECISIONS — the demo

10. **Delete `broadcastHooks`** — REDUNDANT: `makeLoader` already broadcasts a
    bare handle to every obsType. The demo builds a Map by hand for no reason.
    This is a demo bug, not a missing feature.
11. **Drop the `AcvTol` override** — it passes `1e-3`, which IS the default
    (`adjustVCE.m:284`). It does nothing.
12. Label the three GENUINE overrides (`MaxOrder` 40 vs 60, `MaxIter` 8 vs 10,
    `Tol` 1e-3 vs 5e-3) as DEMO-SPEED choices, not recommendations.

#### TWO CORRECTIONS TO CLAIMS MADE EARLIER IN THIS SESSION

Recorded because both are the §1.0 failure mode — asserting from memory instead
of re-deriving from the tree — and the second changed a design decision:
- Claimed the demo's four fluent numbers were undocumented tuning a user would
  inherit blindly. WRONG: the toolbox defaults exist and are documented; three
  are deliberate overrides and the fourth is a no-op. The proposed "move the
  magic numbers into the fluent methods as defaults" would have moved them to
  where they already were. Withdrawn — it is a demo edit, not a toolbox change.
- Claimed nothing in the toolbox showed the sigma trajectory. WRONG:
  `convergence()` already returns a per-pass table and plots two panels. The real
  gap is the PER-STREAM view, which is what decision 7 addresses.

KNOWN WRINKLES, none blocking — recorded so they are not rediscovered as bugs:
- `test_adjustBlock/testRecoverHookFailureNeverBreaksSolve` intends to exercise a
  THROWING hook (`@(ctx) error(...)`) but the identifier that actually reaches
  the catch is `MATLAB:maxlhs` (too many outputs), because an anonymous function
  wrapping `error` has no output to give. The guarantee it checks — the solve
  survives a bad hook — is genuinely verified, so this is a comment/name
  mismatch, not a defect.
- `run_lsqvce('all')` runs every file but MUTES warnings. It is not the gate;
  `run_all_lsqvce` is. Task 2 should remove the chance of confusing them.
- `html/lsqvce.recover.html` and `html/lsqvce.SeparableProblem.html` are redirect
  stubs for removed names. Kept deliberately so old links land somewhere useful;
  the audit reports them as "no backing .m", which is expected.
- `test_speed.m` is a benchmark, not a pass/fail suite, and is deliberately
  excluded from the gate. Run it by hand.


### 1.0 RE-SYNC AUDIT (2026-07-18) — CLOSED, kept as the record

STATUS: **all four defects fixed, suite MATLAB-confirmed green (458/459, 0
failed).** Nothing here is outstanding. It is retained because it documents how
the handover drifted and how that was detected — read it once, then work from
§1.1 onward. The reusable takeaway is the audit script (§4f); the reusable
LESSON is at the end of this section.

This handover had gone OUT OF SYNC with the tree. A file-by-file audit of the
shipped snapshot was run to establish what is actually present. Result: **the
engine code is complete and consistent — every function the handover claims
exists, exists — but three defects were found that the handover did not
record.** All three are now FIXED in the tree; none required an engine change.

**AUDIT METHOD** (repeat this whenever the handover is doubted): check claimed
files exist; cross-check `run_all_lsqvce`'s hardcoded suite list against the test
files on disk; grep for every removed name; compare each public function's
signature against its call sites AND its own help block; check every `+lsqvce/*.m`
has an `html/` page and a `Contents.m`/`helptoc.xml`/`index.html` entry; verify
all doc links resolve. A tokeniser-based function/end balance check was also run.

**WHAT WAS VERIFIED SOUND**
- All 12 newest functions present (`streamContribution`, `contribution`,
  `adjustVCEjoint`, `jointModelOpt`, `psdVCE`, `arcVarianceFactors`,
  `redundancyLags`, `sharedProbe`, `toeplitzCofactor`, `covComponent`, `dctI`,
  `resolveNc`).
- Alias removal (§4a) is GENUINELY complete: `recover.m` and
  `SeparableProblem.m` are gone, `coalesceOpt`/`warnOnce` are gone,
  `PartialMatching` is off, and no removed fluent method (`.shape`, `.scale`,
  `.outer`, `.axes`, `.vceGroupsBy`) is called anywhere in code.
- The `nc` derivation is consistent: `adjustBlock(systems, days)` and
  `adjustVCE(systems, days, ...)` everywhere, and `test_adjustBlock` already builds
  specs carrying `.nc`.
- All 18 test files on disk are registered in `run_all_lsqvce`.
- Structural balance: no genuine mismatch. (`mapreduce.m` and `redundancy.m`
  trip a crude counter because `case 'parfor'` / `case 'spmd'` are STRING
  LITERALS, not block openers — a known false positive, not a defect.)

**DEFECT 1 — six tests written against the REMOVED `nc` signature (REAL BREAK).**
`tests/test_ellmer.m` still called `lsqvce.adjustVCE(4, {}, 1, 'verbose', false)`.
Under the current 2-positional signature that binds `systems = 4`, so
`lsqvce.resolveNc(4)` raises `lsqvce:resolveNc:missing` — the six `vceMode` /
`joint` / `separate` mode tests ERROR before asserting anything. This is a
TEST-side break; the engine is correct. **The suite could NOT have been
463/464 green with this file as shipped**, so that headline count is stale.
FIX: a new local helper `modeSpec()` returns the smallest legal system — a raw
spec struct `struct('sysID',1,'nl',2,'nc',4,'loader',@...)` — which is all
`resolveNc` needs. These six tests only inspect CONFIGURED state (`P.vceMode`,
`P.jointOpt`) and never call `run()`, so the loader handle is never invoked; it
errors if it ever is, so the fixture cannot silently start doing real work.

**DEFECT 2 — twelve functions had NO HTML page (`doc` failed).**
Every other function has a page reachable from `doc lsqvce.<name>`; the 12 newest
shipped without one, and none was in `helptoc.xml` or `html/index.html`
(`Contents.m` covered 11 of 12 — `resolveNc` was missing there too).
FIX: new `docs_src/build_missing_pages.py` parses the standardised help blocks
and emits pages in the existing house style (same CSS, nav, see-also, footer).
It only writes pages that do NOT exist, so it cannot clobber curated ones. All
12 pages generated and registered in `Contents.m`, `helptoc.xml` and
`index.html` (with a new "joint VCE (Ellmer 2018)" group mirroring the
DRIVER -> ADAPTER -> PRIMITIVES layering). Verified: 0 broken links across all
html files, 0 missing helptoc targets, still 0 `<script>` tags.

**DEFECT 3 — stale help/doc text describing deleted things.**
Six help blocks still documented the removed leading `nc` argument
(`adjustVCE`, `adjustBlock`, `Results`, `makeLoader`, `util.log`), `adjustVCE`
and `adjustBlock` still promised the deprecated aliases "for one version" (they
were deleted), and four help EXAMPLES called removed fluent methods. The two
redirect stubs said "renamed" when the names are now REMOVED.
FIX: all corrected. The stubs are KEPT (old links still land somewhere useful)
but now state plainly that the name is gone and calling it is an error; a stray
duplicate `</head>` in both was also repaired.

**DEFECT 4 — `adjustBlock` derived nc BEFORE validating the systems shape
(REAL ENGINE BUG, found by the MATLAB run of the fixed tree).**
§4a of this handover already recorded the required ordering — "in `adjustBlock`,
`resolveNc` must run AFTER `normalizeSystems` validates the shape" — but **the
fix was never applied to the code**. `adjustBlock.m` called
`lsqvce.resolveNc(systems)` on line 96, before `normalizeSystems` on line 109, so
a malformed `systems` reported the confusing generic `lsqvce:resolveNc:missing`
("no system declares a common-parameter count") when the real problem was that
the entry was not a system at all. `testRecoverErrors` caught it:
`verifyError(@() lsqvce.adjustBlock({1},1:3), 'lsqvce:adjustBlock:badCell')` got
`resolveNc:missing` instead.
FIX: the `resolveNc` call moved to AFTER `normalizeSystems`, with the rationale
inline so it is not "optimised" back. `nc` is not used until the accumulator is
constructed, so nothing else moves. This is the ONLY executable engine change of
the re-sync. `adjustVCE` does NOT share the bug: its shape error
(`lsqvce:adjustVCE:badSystems`) fires lazily in `loaderList` at run time, and its
constructor test correctly expects `resolveNc:empty`.

**VERIFIED BASELINE (2026-07-18 20:06, MATLAB, user's box)**
`run_all_lsqvce` → **456 / 458 passed, 1 failed, 1 incomplete.** Per-suite counts
are in `tests/last_run_all.txt`. Reading it per §4d:
- The 1 FAILED was defect 4 above (`test_adjustBlock/testRecoverErrors`), now fixed.
- The 1 INCOMPLETE is `test_parallel/testPrefetchServesFromCacheWhenPrimed`,
  "Filtered by assumption" — the long-known `backgroundPool` filter, not a
  regression.
- `test_ellmer` came back **54/54**, confirming the defect-1 fix.
- Warning inventory: 51 warning lines, 20 distinct messages, ALL of the
  deliberate kinds §2 lists (non-convergence on tiny synthetic runs, datum
  auto-detect, `v'v` cancellation, rank-deficiency probes, the intentional
  hook-failure test). No new type, no warning storm.
- Log size 157,991 bytes — essentially identical to the 157 KB baseline, so no
  hidden warning storm.

**CONFIRMED AFTER THE FIXES (2026-07-18 20:13, MATLAB, user's box):
458 / 459 passed, 0 FAILED, 1 incomplete.** This is the CURRENT GREEN BASELINE.
- 459 executed points, not 458: `test_adjustBlock` went 20 -> 21 with the new
  `testShapeIsValidatedBeforeNcIsDerived` guard, which passes.
- `test_ellmer` 54/54 and `test_adjustBlock` 21/21 — both re-sync fixes confirmed.
- The 1 incomplete is `test_parallel/testPrefetchServesFromCacheWhenPrimed`,
  "Filtered by assumption" (the long-known `backgroundPool` filter). Its
  "Failure Summary" block in the log shows an empty Failed column — it is a
  FILTER, not a failure. Do not read that heading as a regression.
- §4d signals, compared against the 456/458 run: warning inventory
  BYTE-IDENTICAL — same 20 types, same 51 warning lines, same counts, so
  nothing was silently suppressed or newly provoked. Log 155,281 bytes vs
  157,991: the 2,710-byte shrink is exactly the removed failure diagnostic.
  Zero `Verification failed` blocks.

**The audit script now covers both defect classes.** `docs_src/audit_tree.py`
gained a check that flags a `verifyError` expecting a specific shape-validation
id from a `systems` argument declaring no nc — but ONLY while `adjustBlock`
actually derives nc first. Verified both ways: it reports the defect on the
pre-fix tree and is silent on the fixed one, so it cannot cry wolf (§4d).

**THE LESSON (this is the part worth carrying forward).** Defects 1 and 4 are
the same failure mode: *a rule recorded only in the handover is not in force.*
Defect 1 was a test the handover did not know about; defect 4 was a fix the
handover described as DONE while the code still had the bug. Neither is visible
to prose review — re-reading this file would never have found either. So:
- put a rule in the CODE, and give it a TEST, before writing it down here;
- treat any "DONE" claim as a hypothesis until a check re-derives it from the
  tree (`audit_tree.py`) or MATLAB runs it;
- when a count in this file cannot be reproduced, the FILE is wrong, not the run.

### CURRENT POSITION (2026-07-18)
The Ellmer stochastic model (Phases A–E) is COMPLETE. All three interface items
are closed. **The suite is GREEN: 471/472 passed, 0 failed, 1 incomplete**
(MATLAB-confirmed 2026-07-18 20:36; the incomplete is the `backgroundPool`
filter in `test_parallel`). That run covers messaging items 1/3/4/5; ITEM 2
(multi-line stage detail) landed after it and awaits a run. The re-sync found
and fixed four defects to get here (§1.0).

The interface change was verified against a BASELINE: the warning inventory was
IDENTICAL before and after it, with no `renamedMethod` and no `ncConflict`, so
nothing is silently falling back on old behaviour. (The "12 types / 157 KB"
figures once quoted here were superseded by the 2026-07-18 re-sync runs, which
measure **20 distinct warning messages, 51 warning lines, 155,281 bytes** on the
green suite. Use those as the comparison baseline from now on — a change in the
TYPE COUNT is the signal to look at, per §4d.)

**INTERFACE AND MESSAGING — items 1, 3, 4, 5 DONE 2026-07-18** (see §1.0b);
item 2 (multi-line stage detail) is the remaining piece. Then §1.3.

### 1.0b MESSAGING REBUILD (2026-07-18) — items 1, 3, 4, 5

**New file `+lsqvce/+util/verbosity.m`** is the single point of access for
message levels AND banner formatting.

**VERBOSITY WAS NOT PROPAGATED ACROSS INTERNAL BOUNDARIES (2026-07-18, user
caught this).** The user asked why `residuals` still printed with verbose off.
It should not have, and the same root cause produced a SECOND fault in the
opposite direction:

1. **Over-share.** `lsqvce.residuals` defaults its own `verbose` to TRUE —
   correct for a direct call, wrong for one made THROUGH an accumulator that was
   already told to be quiet. `BlockAccumulator.residuals` forwarded `varargin`
   without adding the level, so `adjustVCE` runs with verbose OFF still emitted
   a per-day residual line per system and day: **544 lines in the 478/479 log,
   all from `test_separable`, and the only output those "silent" runs produced.**
   FIX: the forwarder now injects `obj.opt.Verbose` when the caller did not
   name one. An explicit `'Verbose'` still wins, in both directions (tested).
2. **Under-share.** `buildAccumulator` hardcoded `'Verbose', false`, so the
   accumulator inside `adjustVCE` was mute on EVERY pass whatever the user asked
   — meaning the level-2 per-day accumulation detail built for item 2 could
   never appear through `adjustVCE` at all. It was dead code on the main path.

FIX for (2) is not simply "pass `obj.verbose`": `buildAccumulator` is the
`accFcn` VCE rebuilds on every INNER iteration (~4 per outer pass, ~30 per run),
so propagating naively would emit thousands of per-day lines. The level is now an
ARGUMENT, and the two call sites differ deliberately:
- `blockSolveWith` accumulates ONCE per outer pass -> its per-day detail is
  exactly what level 2 is for, so it gets the full level at `stage` and above;
- the `accFcn` handed to VCE repeats ~30x -> held back to level 3.
`iterateNoiseModel` already passed `'Verbose', false` explicitly; that lone
correct call site was the tell that the others had forgotten.

LESSON, and it is the same one as §1.0 defect 4: a default that is right for a
direct call is wrong for an internal one, so EVERY internal call must state the
level rather than inherit a library default. When adding a call to anything that
prints, pass the caller's level explicitly.

**LABELS ARE NOW PLAIN WORDS, NOT BRACKETED TAGS (2026-07-18, user decision).**
The abbreviated form was dropped for two reasons. (a) The codes were INVENTED,
not geodetic: 'scale' and 'shape' are not what anyone calls variance-component
estimation and stochastic modelling, and the user had already rejected those
words once. (b) A tag that abbreviates the message merely repeats it —
`[SCALE   ] variance-component estimation between groups` says it twice.
Labels are now words in a 21-character column (`verbosity.WIDTH`, widen freely),
unbracketed, with the message trimmed to what the label does not already say:

```
Variance components  estimating between groups ...
Variance components  done: 4 group(s), max|sigma-1| = 3.104e-02   1.83 s
Adjustment           block-arrow, current weights ...
Stochastic model     re-estimating temporal correlation per stream ...
Convergence
                     GraceA "kinematicPosition"   sigma =  1.0342   white = yes  (p = 0.310)
                     pair "rangeRate"             sigma =  1.2044   white = NO   (p = 0.002)
                     max|sigma-1| = 2.044e-01 (tol 5.0e-03) -> another pass
```
CASE STILL CARRIES MEANING: Capitalised = outer-loop stage, lowercase = detail
beneath one. `verbosity.labelWidth()` exposes the column so a continuation or a
table indents to it — a hardcoded indent would drift the day WIDTH changes.
THE MAPPING: SCALE->Variance components, SHAPE->Stochastic model,
SOLVE->Adjustment, GATE->Convergence, accum->normal equations,
solve->conditioning, outlier->data snooping (Baarda's term, already used in 9
files), addSys->systems, mapred->parallel execution; resid/redund/vce/pool/
contrib/joint/selectAR/psdIt spelled out in full.

**THE CONVERGENCE TABLE IS LEVEL 1, NOT LEVEL 3 (2026-07-18, user caught this).**
§1 item 3 required per-stream sigmas "in the convergence report, not just
max|sigma-1|". They were built but GATED AT DIAGNOSTIC, so the default level
still showed only the scalar — the requirement was half-met. The per-stream table
now prints at PROGRESS. Level 3 keeps only what it genuinely adds: the per-GROUP
sigmas, which differ from the per-stream view whenever a coarse
`groupComponentsBy` makes several streams share one component.

**BANNER UNIFICATION WAS INCOMPLETE — FOUND IN THE 477/478 LOG (2026-07-18).**
Item 4 was recorded as done when only `adjustVCE` and `BlockAccumulator`'s loop
had been converted; ELEVEN other files still wrote their own literals. The log
proved it: `[resid ]` 499 times and `[resid   ]` 63 times — the SAME tag at two
widths in one run. All 20 remaining sites now route through `verbosity.banner`
(`residuals`, `redundancy`, `vce`, `streamContribution`, `poolGroups`,
`selectArOrder`, `outliers`, `mapreduce`, `iterateNoiseModel`, `contribution`,
`adjustVCEjoint`, `BlockAccumulator`). DELIBERATELY LEFT: `report.m`'s
`[ERROR]`/`[WARN ]` are severity markers in a text table, not stage banners, and
their narrow width is what keeps that column aligned.
`audit_tree.py` gained check 8 for this; verified it reports all 37 sites on the
pre-fix tree and is silent on the fixed one.

**WHY THE LOG SHOWED NO NEW STYLE AT ALL.** `[resid ]` was the ONLY banner in
the entire 477/478 run — no `[SCALE]`, no `[GATE]`, no stage detail. Every test
passes `'verbose', false`, so the whole messaging layer is dark; `residuals`
alone printed because its `verbose` DEFAULTS TO TRUE and internal callers do not
pass the flag. So the suite proves the formatters are correct (they are unit
tested) but exercises almost none of the real output. §1 item 5 said the
messaging was untested; that is now only half fixed — the FORMATTERS are tested,
the WIRING is not. See §1.7.

**MATLAB-CONFIRMED 2026-07-18 20:36 for items 1, 3, 4, 5: 471/472, 0 failed**
(`test_util` 21 -> 34, so all 13 verbosity + 2 render tests pass). §4d signals:
warning inventory IDENTICAL to the previous green run — same 20 types, same 51
lines — and log size 155,295 vs 155,281 bytes. That matters more than the count
here: the `render` fix changed the output path of EVERY message in the toolbox,
and it provoked no new warning and no size drift. ZERO literal backslash-n in
the log. ITEM 2 (below) came after this run and is NOT yet MATLAB-confirmed.

1. **Verbosity levels replace the boolean** (item 1). `0` silent | `1` progress
   | `2` stage | `3` diagnostic, CUMULATIVE (each includes those below), so
   every check is `atLeast`, never an equality test — `v == 2` would silence
   stage detail at level 3, which is how detail gets lost when a level is added.
   BACK-COMPATIBLE: `verbose` still accepts a LOGICAL and `resolve` maps
   `false->0`, `true->1`, exactly the two behaviours the boolean selected. Every
   existing call, test and demo keeps its meaning; levels 2–3 are new capability
   reachable only by passing a number. `obj.verbose` is resolved ONCE in the
   constructor so it is always numeric downstream.
   BOUNDARY NOTE: `lsqvce.vce` and the `report` flags take a LOGICAL and would
   reject a level, so the level is converted at those three call sites
   (`vce` gets `diagnostic(v)` — its inner Helmert iterations are level-3
   material; `report` defaults to `progress(obj.verbose)`).
3. **Per-stream sigmas** (item 3) via the new file-local `sigmaLine`, printed at
   level 3 under `[GATE    ]`. `max|sigma-1|` names how far the worst group is
   from unity but not WHICH group, so a loop stalled on one stubborn stream
   looked identical to one drifting everywhere.
4. **Banner vocabulary unified** (item 4). `verbosity.banner` pads every tag to
   one width, so `[GATE ]`/`[CONTRIB]`/`[SHAPE]` vs `[shape   ]` all align.
   CASE IS MEANINGFUL and preserved: UPPERCASE = outer-loop stage, lowercase =
   per-item detail, and detail lines are indented one level deeper than the
   stage they sit under.
5. **The messaging is now TESTED** (item 5). 13 new tests in `test_util`: the
   level grid, back-compat mapping, banner width/case, and — the point — two
   END-TO-END tests that read the rendered log text back and assert a real
   newline arrived, not the two characters backslash-n.

**LATENT BUG FOUND AND FIXED WHILE DOING THIS.** `log.render` returned a LONE
argument untouched, and `emit` writes with `fprintf(fid,'%s',txt)`, which
expands nothing — so `log.msg('TEXT\n')` put a LITERAL backslash-n in the log.
~20 call sites are written that way (`report`, `snoop`, `robust`, `BatchIndex`,
`Results`, `ftest`, `systematic`, `adjustVCE`) and the suite log shows NONE of
them ever executed, which is exactly why it stayed hidden — and it is the same
bug class the user hit by hand in the joint path. `render` now calls
`sprintf(varargin{1})` on the single-argument path, matching fprintf semantics.
Verified no call site passes a variable or a bare `%` as a lone argument (a
parser over every `log.msg`/`log.warn` call found 0 at risk), and `'%%'` still
renders as one `%`. Two regression tests added.

2. **Multi-line stage detail DONE 2026-07-18** (item 2), to the user's own
   specification, which is recorded here because it is a USER-FACING LAYOUT
   SPEC, not an implementation detail:
```
Accumulate the normal matrix:
 Day 54466
  1. Accumulate observations of type "kinematicPosition":
   -> Group {1,1}: type "kinematicPosition" for "GraceA"   - 00:00:00.519
  2. Accumulate observations of type "rangeRate":
   -> Group {1,2}: type "rangeRate" for "pair"             - 00:00:01.204

Convergence check:
  -> GraceA "kinematicPosition"       sigma =  1.0342   white = yes  (p = 0.310)
  -> GraceB "kinematicPosition"       sigma =  0.9871   white = yes  (p = 0.440)
  -> pair "rangeRate"                 sigma =  1.2044   white = NO   (p = 0.002)
  max|sigma-1| = 2.044e-01 (tol 5.0e-03) -> another pass
```
   Rendered by five helpers on `verbosity` (`heading`, `dayLine`, `stepLine`,
   `groupLine`, `duration`) plus the file-local `convergenceRow` in adjustVCE.
   Indentation encodes the hierarchy (activity > day > step > group); every
   duration is `HH:MM:SS.mmm` in a FIXED COLUMN so a slow group is visible at a
   glance regardless of name length. The convergence table puts BOTH halves of
   the stopping rule on one row -- SCALE (sigma) and SHAPE (whiteness + p) --
   so a row failing one but not the other says which half is stuck, which the
   scalar `max|sigma-1|` cannot.
   PLUMBING: `BlockAccumulator.accumDay` now returns a second output `perType`
   (per-obsType arc count and elapsed time). It is collected UNCONDITIONALLY --
   one `tic` per obsType, not per arc, so the cost is unmeasurable and the
   caller decides whether to print. Both existing callers request one output,
   which is safe in MATLAB. `BlockAccumulator`'s `Verbose` was widened to accept
   a level and is resolved ONCE in the constructor; `mapreduce` takes a LOGICAL,
   so the level is converted at that boundary (`progress(...)`).
   NAME IN THE OUTPUT: the display name is read from `UserOpt.systemName` /
   `.name` if a loader supplies one, else the `sysID` is used. Threading a name
   through the SPEC struct was deliberately NOT done -- that would change the
   accumulator's spec contract for a cosmetic gain.
   LEVEL SPLIT: level 2 gets the indented per-day form; level 1 keeps the
   compact one-line-per-day `[accum   ]` banner, because at progress level the
   point is the shape of the run, not which type took how long.

**A NOTE ON THE DEPRECATED-ALIAS PLAN.** The original sketch proposed keeping
`Verbose` as a deprecated alias mapping `true -> 2`. That was NOT followed, on
purpose: (a) aliases were all deliberately removed on 2026-07-18 (§4a) and
reintroducing one contradicts that decision; (b) mapping `true -> 2` would make
every existing call print MORE than before, so ~30 call sites, all the demos and
the whole suite would silently change behaviour. `true -> 1` keeps the boolean's
exact former meaning, which is what "nothing breaks" has to mean. Levels 2-3 are
new capability, reachable only by asking for them.

STILL OPEN from the messaging list: nothing. Items 1-5 are complete; §1.3
(documentation) is the next agreed step.
2. **Multi-line stage detail**, in the style the user actually uses:
```
Accumulate the normal matrix:
 Day 54466
  1. Preparing partial derivatives for observations of type "satellite":
   -> Group {1,1}: type "kinematicPosition" for "GraceA"   - 00:00:00.519
   -> Group {1,2}: type "kinematicPosition" for "GraceB"   - 00:00:00.368
  2. Accumulate the normal matrix:
   -> Group {1,1}: type "kinematicPosition" for "GraceA"   - 00:00:01.204
```
3. **PER-STREAM SIGMAS** in the convergence report, not just `max|sigma-1|`.
4. The banner vocabulary is INCONSISTENT and should be unified: `[GATE ]` padded
   with a trailing space, `[CONTRIB]` breaking the alignment, `[SHAPE]` vs
   `[shape   ]` mixing case for stage vs per-stream lines.
5. NOTE: the messaging is currently UNTESTED — every test runs `Verbose` off, so
   no test asserts a single line of it. That is how the literal `\n` bug survived
   until the user ran the joint path by hand.

ALSO OPEN: `cfg` in the demo builders may duplicate `processOpt` — the user
flagged this while spotting the `nc` duplication, not yet examined.
**EXAMINED 2026-07-19 — see §1.1a.** It DOES duplicate: four values live in both,
and `demoBuildDay`'s `procOpt` argument is `%#ok<INUSD>`, genuinely unused. The
fix is specified there (not yet built).

RESOLVED (2026-07-18): the two side documents `ELLMER_BUILD_PLAN.md` and
`NEW_CHAT_START_HERE.md` are GONE ON PURPOSE. Per the user's instruction there is
now exactly ONE status-and-plans document — this one (see the note at the very
top). Their durable content already lives here: the Phase A–E measurements and
the "KEY FACTS worth not re-deriving" in §2, the working loop in §0/§4b. Older
sections below still name them as sources; read those as historical citations,
not as files you can open.

`docs_src/build_reportlab.py` (the §1.3 master-PDF builder) also does not exist
yet — §1.3 is scoped but not started, which is correct and expected.

**§1.3 (the six-part documentation) is DEFERRED** until the interface settles.

THE THREE INTERFACE ITEMS, all now closed:
1. ~~naming/surface: joint as a mode flag rather than a separate function~~ DONE.
2. ~~`accFcn(model)` contract for the joint loop~~ DONE 2026-07-18:
   `lsqvce.jointModelOpt(model, streams)` turns the driver's model into loader
   UserOpt, encoding the three conventions that otherwise FAIL SILENTLY (the
   loop diverges rather than erroring): sigma2 is a VARIANCE but a spec's sigma
   is a STANDARD DEVIATION; the override must land at
   covOverride.bySys.(satKey).(obsType); psdMode is read from the UserOpt, not
   the spec. Six tests, one per convention. adjustVCEjoint's help points at it.
3. ~~public vs internal API~~ DECIDED 2026-07-18: **keep them all public**, but
   LAYER them in Contents.m so the reader knows what to call. Reasoning: they
   are not incidental helpers (dctI, covComponent, toeplitzCofactor,
   redundancyLags all have legitimate standalone uses in LSQ/VCE work); they
   already carry full help with equation references; and the stated goal is a
   READABLE REFERENCE IMPLEMENTATION, for which inspectable primitives are a
   feature. Contents.m now presents them as DRIVER (adjustVCEjoint) -> ADAPTER
   (jointModelOpt) -> PRIMITIVES (psdVCE, arcVarianceFactors, redundancyLags,
   sharedProbe, toeplitzCofactor, covComponent, dctI), with the note that most
   users want adjustVCE(...,'vceMode','joint') and never touch the rest.

**ALL THREE INTERFACE ITEMS ARE NOW CLOSED.** §1.3 (documentation) is the next
agreed step.

**§1.3 (the six-part documentation) is DEFERRED at the user's request** until the
interface settles — writing the manual against an API about to change is rework.

The items below are the original roadmap. Rough priority:

### 1.0d MESSAGING REVIEW ROUND 2 (2026-07-19, from demo_messaging)

**THREE LEVELS, NOT FOUR.** Level 3 was merged into level 2 and removed. Once
the inner accumulator stopped repeating its per-day block on every Helmert
iteration, the only thing 3 still added was the per-group sigmas -- too thin to
justify a tier a user then has to choose between. `resolve(3)` now ERRORS with a
message saying so, rather than clamping: someone asking for 3 wants more output,
and a silent clamp leaves them wondering why nothing changed. The two predicates
`stage()`/`diagnostic()` collapsed into one, `detail()` -- two names for one test
is the ambiguity being removed everywhere else.

**THE SYSTEM NAME NOW TRAVELS WITH THE SPEC.** Previously refused as "cosmetic";
reversed, because four separate review points needed it. `oneSpec` and
`adjustVCE.blockSpecs` carry `.name` from the Loader, `processSystem` takes it
and republishes it on `uo.systemName` so the SECOND pass (`lsqvce.residuals`,
which sees only the stored UserOpt) labels its lines identically. Users saw `x1`
and asked what it meant: `matlab.lang.makeValidName('1')` = `'x1'`, correct for a
struct field and meaningless on screen.

**BUGS FOUND IN THIS ROUND**
1. `groupLine([1 t], ...)` HARDCODED the system index, so every system's types
   printed as `{1,k}` -- system 2's lines claimed to be system 1's. Not spotted
   by the user; found while reworking the layout.
2. The residual line passed the observation COUNT into the obsType slot, so it
   read `type "4200 residual(s)" for "1"`.
3. The residual report collapsed all observation types into ONE line per day, so
   a three-stream system looked like a one-stream system. NOT a calculation bug
   -- `otMap` had the per-type numbers all along, the display threw them away.

**FINAL WORDING (2026-07-19).** The pass banner is `===== Iteration 1 of 6 =====`
-- "outer pass" was internal jargon; a user reads an ITERATION counter. The same
word now runs through the non-convergence warning ("did not converge in %d
iteration(s)"), both plot x-labels and the demo output, so one concept has one
name. The VCE's own message still says "VCE did not converge in N iterations",
which is unambiguous because it names VCE.
`v'v` was also spelled out in the cancellation warning ("The weighted sum of
squares v'Pv ..."), matching the residuals block -- it IS v'Pv, the numerator of
sigma0^2, and the user had to ask what the abbreviation meant.

**THE DEPRECATED STAGE NAMES SURVIVED THREE ROUNDS OF SEARCHING (2026-07-19).**
`SCALE / SHAPE / SOLVE / GATE / OUTER PASS` were still in the level-1 and
level-2 output. Two reasons, both instructive:
1. **The words lived on CONTINUATION lines.** `log.msg(['...' ... '...SOLVE ->
   SHAPE temporal AR -> GATE)'])` puts the call on line 1 and the offending text
   on line 2, so every line-by-line search for "a log call containing SCALE"
   found nothing. Audit CHECK 11 now JOINS continuations before matching -- its
   first version did not, and it too reported the tree clean.
2. **One of them was in a WARNING, not a status line.** `notConverged` said
   "did not meet both gates". It surfaced only because the user asked for the
   warning inventory to be read carefully -- the pass/fail count and the status
   output were both clean. That is the case FOR the §4d habit, made concrete.

Also lost to the aborted-edit-batch problem above: the 'each pass: SCALE ...'
sub-line had been rewritten once already and the write never happened.

WHAT CHANGED: the pass banner is now `===== Outer pass 1 of 6 =====` with a
plain sub-line; the warning says "did not converge ... the variance factors and
the residual whiteness were not both satisfied"; ~60 help-text occurrences were
reworded. LEFT ALONE ON PURPOSE: `SHAPE of \`systems\``, `SHAPE of your
problem`, `SHAPE (2-D cell, 14 columns)` and `SCALE each epoch by its Cholesky
factor` are ordinary English, and `verbosity.m` QUOTES the retired bracketed
form in its own help to explain what it replaced.

CHECK 11 guards only USER-VISIBLE strings (log.msg/log.warn/error/warning/
title/xlabel/ylabel/fprintf). Comments stay prose -- "the shape of the noise" is
a legitimate sentence and forbidding it would push the check into being ignored.

**A PROCESS FAILURE WORTH RECORDING (2026-07-19).** An edit script applied
several substitutions and wrote the file ONCE at the end. One `assert` in the
middle failed, the script aborted -- and the substitutions that had already
SUCCEEDED were discarded with it, silently, because nothing had been written yet.
Two of them were lost that way and only surfaced as a test failure two rounds
later (`testVerbosityNumericLevelsPassThrough` still looping `0:3`).
THE RULE: **verify every match BEFORE mutating anything, or write after each
successful substitution.** A partially-applied edit that reports success is worse
than one that fails outright.

**VOCABULARY: `vtv` IS `v'Pv`.** The code name is `vtv` because it accumulates
`w'w` where `w` is the WHITENED residual, so it is the WEIGHTED sum of squares --
the numerator of `sigma0^2`, not a sum of residuals. It now prints spelled out
("weighted sum of squares v'Pv") because the abbreviation meant nothing to a
reader, which the user confirmed by having to ask.

**THE GATE REPORTED 0 FAILED WHILE A WHOLE FILE DID NOT LOAD (2026-07-19).**
The most serious thing in this session. `test_util.m` ended up with FIVE
duplicated function names after an editing slip, MATLAB refuses to load such a
file at all, and `run_all_lsqvce` printed:
```
result : 447 / 448 passed, 0 failed
  test_util            :   0/NaN pass, NaN fail, NaN incomplete
```
**41 tests silently vanished and the headline still said the suite was green.**
The `catch` branch recorded NaN counts and `continue`d, so `grandF` never saw it.
TWO fixes, because either alone would have left the hole half open:
1. `run_all_lsqvce` now COUNTS load errors, names them on the headline
   ('N FILE(S) FAILED TO LOAD - the count above EXCLUDES them'), marks the
   console line '<-- NOTHING IN THIS FILE RAN', and writes a section into the
   log so the failure is visible in the artefact that gets shared.
2. `audit_tree.py` CHECK 10 flags duplicated function names before MATLAB ever
   sees the file. NOTE its first version was broken in an instructive way: the
   output-argument group and its '=' were two separate optionals, so the name
   group matched the LAST LETTER of every function ('testFoo(' -> 'o') and every
   file looked duplicated. A check that reports 114 findings on a healthy tree
   is worse than no check -- it trains you to ignore it. Verified both ways
   after the fix.

THE RULE: **a count of failures is not a statement about coverage.** Ask what
did not run, not just what did not pass.

**A BUG I SHIPPED AND THE USER HIT IMMEDIATELY (2026-07-19).** Adding `.name`
to the block spec broke `demo_messaging` with *"Subscripted assignment between
dissimilar structures"*. The spec struct is written out in FOUR places -- a
`repmat` TEMPLATE and a `struct(...)` literal in each of `adjustBlock` and
`adjustVCE` -- and I updated two of them. MATLAB only rejects the mismatch when
the loop actually RUNS, so nothing caught it: not the parser, not the audit, not
the suite (the messaging tests would have, but the demo ran first).
FIX: both templates now carry `name`, with a comment saying why they must track
`oneSpec`. `audit_tree.py` gained CHECK 9, which compares the block-spec field
list across all four sites and reports the differing ones by file and line;
verified it fires on the broken tree and is silent on the fixed one.
LESSON, and it generalises: **a preallocation template is a duplicate
declaration of a type.** When adding a field to a struct that is assigned into a
preallocated array, the template is the half you will forget.

**LAYOUT NOW SHARED.** Three renderers in `verbosity`: `typeLine` (per-type work,
`Type {sys,type}: "x" of "system"` + optional count + time), `streamLine` (the
one shape behind Stochastic model / Contribution / Convergence) and `fieldBlock`
(titled key-value block, blank line either side). `heading` gained a leading
blank line and lost its colon, so the four or five blocks of a pass separate
without anyone adding blank lines by hand. `stepLine` and `groupLine` were
deleted -- the numbered step restated the type the line already named.

### 1.0c THE TWO MESSAGING-TEST FAILURES (FIXED 2026-07-18)
The wiring test earned its keep on its first run: 484/487, 2 failed.

**FAILURE 1 was the TEST being wrong, not the code.** `testMessagingLevelZeroIsSilent`
asserted "no output at all" at level 0, but the capture contained two WARNINGS.
That is CORRECT behaviour and worth stating as a rule: **warnings are not
verbosity-controlled.** A run that failed to converge must say so even at level
0, or a silent batch job hides its own failure — and `lsqvce.util.log.warn`
deliberately raises a real MATLAB warning so `verifyWarning`/`lastwarn` work.
The test now asserts the ABSENCE OF STATUS MARKERS rather than an empty file
(filtering `Warning:` lines would be fooled by a long warning that wraps).

**FAILURE 2 was a REAL third hardcoded silence.** `addSystems` set
`optL.Verbose = false` before mapping, with the comment "ordered printing
happens below". That is right under `parfor`/`spmd`, where workers' lines
interleave unpredictably — but it applied in SERIAL too, where there is nothing
to interleave, and it is the only path `buildAccumulator` uses. So the level-2
per-day block was produced inside `processSystem` and then thrown away on every
run. FIX: silence only when the requested mode is not serial (`taskSilenced`),
and print the ordered per-system substitute only when the task WAS silenced —
otherwise serial would print both.

THIS IS THE THIRD INSTANCE of the same bug class in two days (the others:
`residuals` inheriting a library default, `buildAccumulator` hardcoding false).
The pattern: **a hardcoded verbosity anywhere is a bug waiting to happen.** If a
component must be quiet for a structural reason, gate it on that reason
(`taskSilenced`), never on a literal `false`.

### 1.1 Streaming / block-arrow Hutchinson VCE — CLOSED (verified 2026-07-18)
VERIFIED IN THE CODE, not assumed. `lsqvce.vce(..., 'traceMethod','hutchinson')`
calls `lsqvce.streamContribution(acc, sol, ...)` at both the per-iteration and
the final-redundancy points. That primitive is exactly what this topic
specified: it estimates `tr(N^-1 N_g^cc)` by Hutchinson probing over a SINGLE
data sweep, never forms `N_g^cc` (which would be O(G*nc^2) and would destroy the
streaming property), and mirrors the per-arc leverage path in
`lsqvce.redundancy`. `vtv` stays exact; only the trace moves.

**A CORRECTION WORTH RECORDING.** This was closed by `streamContribution`, NOT
by the Ellmer work. They are different pieces: Ellmer (§2) estimates the
STOCHASTIC MODEL — per-frequency amplitudes of the temporal covariance;
this topic was about the TRACE ESTIMATOR behind the VCE redundancy numbers. Both
are done, but attributing one to the other would put a wrong justification on
record, and wrong justifications are how this file drifted before (§1.0).

### 1.2 Contribution matrix — CLOSED, with one gap now filled (verified 2026-07-18)
- **Level 1 (scalar share): wired into the loop.** `adjustVCE` computes the
  share every outer pass from `vc.r` (free — VCE already returns it), stores the
  trajectory in `R.history.contribution`, and prints it under the `Contribution`
  label at verbosity level 2. A still-drifting share on the last pass is the
  signal that the loop has not settled.
- **Level 2 (per-parameter `diag(N^-1 N_g^cc)`): implemented** in
  `lsqvce.contribution(..., 'Level', 2)`, over the SAME streaming primitive as
  §1.1, so the full dense `C_g` is never formed. Invariants hold:
  `sum_g t_g = nc` and `sum_g diag(C_g) = ones(nc,1)`.
- **NOT auto-run at convergence, deliberately.** The original plan said "once at
  convergence"; that was not implemented and on reflection should not be. Level 2
  costs an EXTRA streamed sweep with Hutchinson probes, and most runs do not want
  to pay for it. It is opt-in.
- **GAP FOUND AND FILLED:** asking for it meant reaching into `R.acc` and
  `R.sol` by hand, while the analogous residual pass had `R.residuals(...)`.
  Added `R.contribution(...)`, which forwards to `lsqvce.contribution` unchanged
  (two tests: it matches the direct call, and Name/Value pairs reach through so
  the accessor cannot silently cap the caller at Level 1).

### 1.3 Documentation rebuild — ONE master PDF via ReportLab (user-specified)
Replace the current doc set with a single MASTER PDF built with **ReportLab**
(Platypus flowables), which the user prefers over the LaTeX/HTML approach. DECIDED
by the user: (a) ONE master PDF with six parts (single TOC + cross-references),
NOT six separate PDFs; (b) the API part is AUTO-GENERATED from the standardised
help blocks, with a few CURATED side fields per function. Doc USAGE is unchanged —
the PDF sits in the toolbox root, opened as now.

**Tooling (verified available in sandbox):** ReportLab 4.4.10, matplotlib 3.10.8,
pdflatex/pdftoppm/Pillow. New script `docs_src/build_reportlab.py`
(SimpleDocTemplate + Paragraph/Table/Image/TableOfContents). REUSE the existing
`_render_eq` (LaTeX→PNG, in build_manuals.py) so equation images are identical to
today; render diagrams/decision-trees to PNG with matplotlib and drop in as Image
flowables. Keep the authored PROSE where it is (the section bodies in
build_manuals.py) — this changes the RENDERER, not the content, except for the new
sections below. Ship the PDF in the toolbox root; keep the HTML manuals too unless
the user says otherwise (both are JS-free and doc-usable).

**Six parts of the master PDF:**
1. **Theory — extended to cover ALL approaches.** Gauss-Markov + BLUE; flat solve
   (`adjust`); streaming block-arrow QR + Schur elimination (`adjustBlock`); the
   separable SHAPE/SCALE model + its fixed-point iteration (`adjustVCE`); Helmert
   VCE (exact) + matrix-free Hutchinson; redundancy/reliability; and the
   contribution matrix C_s = N⁻¹N_s (§1.2). Each block ends with "which function
   implements this".
2. **Workflow — the LOOPS, one flowchart per approach.** `adjust` (single solve),
   `adjustBlock` (accumulate→eliminate→solve), `adjustVCE` (outer SHAPE↔SCALE with
   the inner VCE loop nested), showing where data is read, where each iteration
   gates, and where plot hooks fire. Map the new `[SCALE]/[SOLVE]/[SHAPE]/[GATE]`
   observability banners onto these diagrams.
3. **User manual — decision trees + build-your-own-case.** A top-level decision
   tree (size of u? need a noise model? exact vs Hutchinson?) routing to the right
   entry point, then a STEP-BY-STEP recipe PER APPROACH: wire data into a
   makeLoader day-builder, choose common vs local parametrisation, set options,
   read the result. `examples/demo_workflow.m` is the front of this.
4. **Interpretation guide — expanded.** Per diagnostic LEVEL: what the report
   shows, what results are POSSIBLE at that level, and what a given value MEANS
   (sigma0 far from 1, redundancy too low, a stream that won't whiten, contribution
   shares, rcond in the solve, etc.). The "read your output" reference.
5. **Failure catalogue — NEW.** Each `lsqvce:...` error/warning id: what it means,
   why it fires, how to fix it. The id set is well-defined and enumerable
   (rankDeficient, priorNotPD, noCep, notHomogenised, notConverged, cancellation,
   nonStationary/unstable, noArcs, badX, ...). `grep -rho "lsqvce:[a-zA-Z]*:[a-zA-Z]*"
   +lsqvce/` to enumerate; make it complete and systematic.
6. **API reference — AUTO-GENERATED (65 functions).** Parse each function's help
   block for: purpose, one-line description, INPUTS/OUTPUTS with
   type-dimension-description. Add CURATED side fields held in a small per-function
   side-table (a dict in the build script): usual POSITION in the processing chain,
   the RELATING THEORY (cross-ref to part 1), and BEST PRACTICES. Auto-gen keeps it
   in sync with code; the curated fields are the only hand-maintained part. The help
   blocks are already standardised (INPUT/OUTPUT with sizes) — lean on that.

**Build-first steps:** (1) enumerate error ids for part 5 and the function list for
part 6; (2) stand up build_reportlab.py with the master-doc skeleton + TOC +
`_render_eq` reuse; (3) matplotlib helpers for the flowcharts/decision-trees;
(4) the auto-gen help parser + the curated side-table; (5) port the existing
theory/interpretation prose in, then extend. Validate: PDF opens, TOC resolves,
every equation/diagram image present, no broken cross-refs. Keep it JS/Java-free
(trivially true — ReportLab emits a plain PDF).

### 1.3b EXAMPLES CONSOLIDATED (DONE 2026-07-18, user decision)
`examples/gravity/` and `examples/accum/` are REMOVED. `gravity/` held
`iterateSeparable.m` (296 lines) which HAND-ROLLED the loop `lsqvce.adjustVCE`
now provides — its own help admitted this — plus `graceLoader.m` (331 lines)
superseded by `lsqvce.makeLoader`. A superseded parallel implementation sitting
beside the template a user is meant to adapt is the worst possible signal about
which one to copy.

`examples/` is now: `demo_separable_data.m` + `demo_separable.m` (THE GRACE-FO
template), `demo_lsqvce.m`, `demo_combine.m`, `demo_sensors.m`,
`demo_workflow.m`, `bench_qr_vs_normal.m`.

Three library help blocks pointed into the removed folders and were repointed;
one of them (`residualAcv`) wrote `lsqvce.graceLoader`, a package prefix that
function never had. `Contents.m`, `README.md` and `html/index.html` updated.

### 1.4 OPEN ITEMS — DETAILED (re-verified against the code, 2026-07-19)

Three items. Nothing is blocking and nothing is a defect. Sub-packaging and the
`BlockAccumulator` split were DROPPED by the user on 2026-07-19; real-data
processing happens outside this development.

---

#### A. Wire `poolGroups` into the `adjustVCE` chain — OPEN, recommended

**This is a missing DATA PATH, not a convenience wrapper.** An earlier note
called it the latter; that was wrong.

WHAT EXISTS: `+lsqvce/poolGroups.m`, documented, with tests in
`test_lsqvce_ext.m`. It reads a finished VCE result, finds groups whose
redundancy is below a floor, and returns a remapping (`.map`, `.remap`,
`.indexRemap`, `.newGroups`, `.changed`). The rationale is sound: a variance
component estimated from r_g observations has a relative standard error of about
0.5*sqrt(2/r_g), so below a floor the estimated sigma_g is mostly noise.

THE GAP, VERIFIED: `poolGroups` consumes an **`lsqvce.vce`** result, but the
manual's recommended path is the **`adjustVCE` fluent chain**, and `R.history`
records only:
```
arOrder contribution converged groups jointMaxDev jointUpdated
sigma streams white whiteAlpha whitePval
```
**No per-group redundancy.** That is exactly the number pooling needs.

WHY IT MATTERS: "coarsen `groupComponentsBy`" is the standard advice when the
variance factors oscillate -- it appears in manual §2.4, §5.4 and the §4.9
checklist, and in `varComp`'s `singular` error text. A user following that advice
today has to work the new grouping out by hand.

SHAPE OF THE FIX: record r_g per group per pass in `hist`, then expose it as
either `R.suggestPooling()` (matching `summary`/`convergence`, which already live
on `Results`) or a fluent `.poolThinGroups(floor)`. **Which one is the design
question; the history change is the work.** Needs its own round with tests.

---

#### B. Data snooping is NOT connected to the block-arrow path — OPEN

VERIFIED 2026-07-19: `lsqvce.snoop` is **flat-path only**, and
`adjustBlock`/`adjustVCE` have **no outlier rejection at all**. Grep of
`adjustBlock.m`, `adjustVCE.m`, `Results.m` and `BlockAccumulator.m` for `snoop`
returns nothing. Manual §3.4 states this plainly, so the DOCUMENTATION is
honest -- but the capability gap is real.

WHY IT IS NOT A SIMPLE FORWARDING JOB, and this is the useful part:
- `snoop(sol, D, batchFcn, ...)` takes the structs `lsqvce.adjust` returns plus
  **the same batch function**, because it RE-SOLVES after every rejection.
- The block path has no `batchFcn`. It has LOADERS, arcs, and a Schur
  complement; observations are streamed and never all resident.
- **The diagnostics half already exists**: `lsqvce.diagnoseAccum` produces the
  same `D` struct `snoop` reads (residuals, r_i, w_i/tau_i, `flagIdx`,
  `contamLevel`), with the common-parameter redundancies CORRECTED for the
  eliminated arc-local parameters. So the inputs are available.
- What is missing is the **re-solve loop**: rejecting observation i means
  rebuilding that arc's contribution and re-forming the Schur complement, not
  re-running a batch function.

DESIGN QUESTIONS TO SETTLE FIRST (do not start coding before these are answered):
1. **Reject or down-weight?** On a streamed path, a rejection changes an arc's
   row set. Down-weighting to zero keeps the geometry fixed and may be far
   cheaper -- and it is what `lsqvce.robust` already does on the flat path.
2. **Whole arcs or single observations?** An arc whose local block is
   near-singular is a MODEL defect, not an outlier (manual §2.5, eq 2.10) --
   the accumulator already flags it. The two cases need different treatment.
3. **One at a time?** The flat path rejects singly because of MASKING (manual
   §5.2, eq 5.5): the blunders inflate s_0, which is the denominator of their
   own test. That argument is unchanged on the block path, so a batch rejection
   would be wrong for the same reason -- but one-at-a-time re-solving a streamed
   system is expensive. This is the real cost question.
4. **Which scale?** `snoop`'s `tooContaminated` guard and the `'mad'` fallback
   assume a single global s_0. With per-stream variance components there is no
   single scale, so the guard needs rethinking rather than porting.

RECOMMENDED SEQUENCING: settle 1-4, then prototype against
`demo_grace_dataLoader` with an injected blunder before touching the engine.

---

#### C. `docs_src/build_manuals.py` is GUARDED — a STANDING HAZARD, re-verified

The user asked whether this is obsolete now that the master PDF exists.
**It is NOT.** Checked 2026-07-19:
- `build_manuals.py` generates FOUR files -- `lsqvce_user_manual.html`,
  `lsqvce_reference_manual.html`, `lsqvce_interpretation_guide.html`,
  `lsqvce_choosing.html` -- and **all four are still shipped and live**
  (30 KB, 22 KB, 13 KB, 8 KB), linked from `helptoc.xml` (4 references) and
  `html/index.html` (3). They are what `doc lsqvce` opens.
- The master PDF replaced the three old **PDFs**, which were deleted. The HTML
  was deliberately KEPT as the browsable per-topic form.
- The guard's premise still holds exactly: the SOURCE has **0** occurrences of
  user-manual Parts IV-B/IV-C and **5** dead fluent method names; the shipped
  OUTPUT has the reverse (Parts IV-B/IV-C present, 0 dead names). Running the
  generator would swap a correct file for a broken one.

So it is neither obsolete nor a task -- it is a trap with a guard on it. Leave
`_STALE_GUARD = True` until someone back-ports Parts IV-B/IV-C and the renames
into the `USERMAN` body.

### 1.7b TESTS FOLDER FLATTENED (DONE 2026-07-18)
`tests/accum/` is gone; `tests/` is FLAT with 19 files. Three problems fixed:
1. **The trap**: `runtests('tests')` does not recurse, so it silently ran 7
   files and skipped 12. Flattening removes it outright — that command now runs
   everything. (`run_all_lsqvce` remains the gate because it also keeps the
   warning inventory, §4d.)
2. **`accum` no longer described its contents** — it held the VCE loop tests,
   contribution, precision, reliability.
3. **Five files were named after RELEASES**, and one after a removed function:
```
   test_v3      -> test_vceMath        leverage, Helmert VCE, snooping, reg, multi-axis
   test_v4      -> test_parallel       mapreduce, addSystems, redundancy, PSD iteration
   test_v7      -> test_vceFunction    lsqvce.vce contract + diagnoseAccum bridge
   test_v8      -> test_signalCore     outliers, varWhiten, ar2arma
   test_v9      -> test_streamUtils    selectArAcv, whitenessCorr, collectResiduals, ...
   test_recover -> test_adjustBlock    'recover' was REMOVED 2026-07-18
```
   `test_vceMath` / `test_vceFunction` follow the authors' own distinction:
   the first checks the VCE *mathematics* by hand, the second the *function*
   contract. Thematically-named files were left alone to keep churn minimal.

Done by `docs_src/flatten_tests.py` (kept as the record of what moved). It also
rewrote each renamed file's `function tests = <name>` entry point, fixed the
`examples/` path in `test_separable` that assumed the extra level, and updated
every reference in the runners, `lsqvce_setup.m`, README, this file and
`html/index.html`. Verified by §4f check 1 (18 registered == 18 on disk).

### 1.7c THE TWO RUNNERS: KEEP BOTH, FOR NOW (2026-07-18, user asked)
They are NOT duplicates:
- `run_all_lsqvce` is the SHIPPING GATE: file-by-file, per-file pass counts,
  writes `tests/last_run_all.txt`, and does NOT mute warnings (the inventory is
  a signal, §4d).
- `run_lsqvce(target)` is for ITERATION: one file or one name pattern, one-line
  summary, per-target log, warnings MUTED by default because the suite
  deliberately provokes them.
The overlap is real (both discover, run, summarise, log) and merging them into
one entry with an `'all'` target would match the "single point of access" rule.
NOT DONE, deliberately: it changes the command the user types while we are still
iterating, and a broken runner costs a whole round-trip. Listed as an open task.

Flattening DID break the two folder-based shortcuts (`'accum'`, `'core'` were
"the two halves"), so they were replaced with thematic single-file targets:
`messaging`, `accum`, `block`, `parallel`, `vce`, `vcemath`, `signal`,
`streams`, plus `all`. NOTE `run_lsqvce('all')` is not the gate — it mutes
warnings.

### 1.7 THE MESSAGING WIRING IS STILL UNTESTED (new, 2026-07-18)
The formatters are unit tested (20 tests in `test_util`), but NO test runs the
engine with messaging ON, so nothing checks that the right line appears at the
right level from the right place. Evidence: the 477/478 log contains exactly one
banner type, `[resid ]`, because every test passes `'verbose', false` and
`residuals` only printed by accident (its `verbose` defaults to TRUE and internal
callers omit the flag).

BUILT 2026-07-18: six tests at the end of `tests/test_separable.m`
(`testMessagingLevel*`, helper `captureRun`). They run a tiny `adjustVCE` at
levels 0/1/2/3 with the log routed to a temp file and assert on the CAPTURED
TEXT —
- level 0 produces NOTHING,
- level 1 has the outer-pass banner and `[GATE    ]` but NO `[shape   ]`,
- level 2 adds the indented per-day block and `[CONTRIB ]`,
- level 3 adds the `Convergence check:` table,
- and each level is a strict superset of the one below.
That last property is the one worth pinning: it is what makes the levels
cumulative, and an equality test (`v == 2`) anywhere would break it silently.

ALSO WORTH FIXING while there: `residuals`' `verbose` default of TRUE is why a
"silent" suite still emitted 562 lines. Every other entry point defaults to
printing too, which is defensible for a user-facing call but means an internal
caller must remember to pass the flag. Consider making internal calls pass the
caller's level explicitly rather than relying on the default.

### 1.6 Verification still pending
The observability/plot-hook change this section used to track was confirmed long
ago. WHAT IS PENDING NOW: **messaging item 2** (multi-line stage detail, §1.0b)
shipped after the 471/472 run, so it is the only unverified change in the tree.
First action next session: `run_all_lsqvce`, expect ~478/479, and check the
warning inventory against 20 types / 51 lines — `accumDay` gained an output, so
a new warning type there is the signal to look at.

---

## 2. CURRENT STATE — what is DONE

### Ellmer (2018) stochastic model — COMPLETE, MATLAB-confirmed (2026-07-18)
(The suite stood at 456/457 when this landed; the CURRENT figure is in §1 —
counts in §2 are historical markers of when a piece was confirmed, not the
present total.) Nine new functions, Phases A–E of the Ellmer build plan:
- `dctI`, `covComponent`, `toeplitzCofactor` — type-I DCT and the lag/Toeplitz
  primitives (Ellmer eqs. 6.5.6, 6.5.12–6.5.17)
- `sharedProbe`, `redundancyLags` — the once-per-iteration `U^-1 Zbar` and the
  per-arc `s_n = tr(Rbar V_n)` / `Omega_n` kernel (eqs. 6.5.45–6.5.56)
- `psdVCE`, `arcVarianceFactors` — per-frequency amplitude VCE and the arc
  factors with the gamma normalisation and robust mean (eqs. 6.5.19–6.5.30, 6.5.57)
- `adjustVCEjoint` — the standalone 10-step driver (§6.5.7)

**`adjustVCE` now has a `vceMode` flag**, the user-facing surface:
```matlab
lsqvce.adjustVCE(nc, systems, days, 'vceMode','joint')
P.joint('numLags',60, 'probes',100, 'psdMode','toeplitz', 'smooth',0)
P.separate()                       % back to the default alternating loop
```
Default is `'separate'` (alternate SHAPE↔SCALE, Ljung-Box gate). `'joint'`
co-estimates every variance component in one Helmert VCE per pass, as Ellmer
does — the two can converge to DIFFERENT models because the components are
coupled, which is exactly what the gamma normalisation guards against.

KEY FACTS worth not re-deriving (all MEASURED in Python against exact formulae;
merged 2026-07-18 from the former ELLMER_BUILD_PLAN.md, which is gone — this is
now the surviving record, so do not go looking for that file):

- **The white start is a DEGENERATE FIXED POINT, not a numerical accident.**
  Ellmer's literal `C=[1 0 0 ...]` is a fixed point of his own weighting: the
  loop reports convergence while the PSD stays wrong. Default C0 is deliberately
  non-degenerate; a delta start warns.

- **A PSD on `nf` frequencies spans EXACTLY `2*(nf-1)` independent lags**, so the
  Toeplitz matrix it generates is rank deficient beyond that BY CONSTRUCTION —
  true even for a flat white spectrum, so it is not a data or noise problem. The
  transition is sharp (nf = 16, flat spectrum): m = 30 -> min eig +5.0e-01;
  m = 31 -> 0.0; m = 32 -> -3.9e-16. THE RULE, verified sharp: **`nf >= m/2 + 1`**
  (m = 300: nf = 150 gives min eig -8.7e-16, nf = 151 gives +6.2e-02).
  Consequences: `NumLags` defaults to `'auto'` (sized `floor(mMax/2)+1` from the
  longest arc), `covModel`'s toeplitz branch WARNS
  `lsqvce:covModel:psdUnderResolved` when the rule is violated, and `psdRidge`
  (default 1e-8 on the zero lag) remains only as a LAST-RESORT guard that MASKS
  the modelling error rather than curing it.
  TESTED AND REJECTED: spectral floors (1e-6..5e-2) and 3-point smoothing do NOT
  fix it; nor does writing the truncated ACV as an `acf` override (truncating an
  autocovariance at lag L is not itself a valid covariance — zeroing the tail
  needs ~30% variance inflation to become PD again).

- **The AR failure was NOT noise — it was asking for more AR parameters than the
  spectrum has independent lags.** Levinson-Durbin on an ACV from a PSD with
  nf = 12: p = 5 gives max|root| 0.00-0.66 (stationary), p = 10 gives 0.02-0.92,
  p = 30 gives EXACTLY 1.000000 — and it fails even for a flat white spectrum.
  The prediction error collapses to ~0 and Levinson divides by it; that
  suspiciously round 1.000000 is the tell. Low-pass filtering would not have
  helped.

- **`psdMode` defaults matter and are read from the UserOpt, not the spec.**
  `'toeplitz'` builds the covariance from the PSD DIRECTLY (Ellmer's own eqs.
  6.5.6/6.5.17): no AR, no Levinson, no stationarity requirement, cost O(b^3) per
  batch. `'ar'`/`'arma'` fit a filter (cheap, O(m*p)) but an estimated PSD can fit
  a NON-STATIONARY AR whose autocovariance does not exist, so `armaACV` throws.
  `applyJoint` once wrote a `'psd'` spec without setting psdMode and silently
  inherited `'ar'` — that, not the estimator, caused
  `lsqvce:armaACV:nonStationary`. Hence `jointOpt.PsdMode` defaults to
  `'toeplitz'`. `psdIsUsable` must probe with the SAME mode AND a realistic
  length (2*numLags, capped 256): an 8-sample probe passes models that then fail
  on 300-row arcs.

- **`numLags` and `nf` are DIFFERENT numbers, and the distinction is a sizing
  artifact of a precompute, not a change to the noise model.** `numLags` is the
  ESTIMATOR's statistical resolution (small is better — each lag comes from
  finite residuals, so more lags means noisier alpha2_j). `nf` is how finely the
  spectrum must be SAMPLED for psd2acv to generate a full-rank Toeplitz matrix
  over the LONGEST window it is asked for (the algebraic rule above). They differ
  only because `covModel` precomputes one autocovariance vector per stream, sized
  for the whole batch window, so it can serve every batch of that stream. IN
  ELLMER THE TWO COINCIDE: his Nmax is the max arc length and the DCT maps Nmax
  lags to Nmax frequencies. Set BatchSize so the precompute window equals the arc
  length and the toolbox collapses to exactly that structure.
  `psdVCE('Smooth', h)` (moving average, reflected edges) is the legitimate
  variance-vs-resolution trade for a NOISY spectrum — explicitly NOT a cure for
  rank deficiency.

- **Rows are EPOCH-MAJOR**, so joint-mode lags must stride by `nAxis`
  (`lagsPooledOverAxes`), or the temporal law comes out stretched by the axis
  count.

- **The covariance frequency grid must use PHYSICAL frequency.** It was once
  built on normalised frequency (0..0.5) while the data has dt = 5 s (Nyquist
  0.1 Hz), claiming power over a band 5x too wide — 184 psd2acv warnings. The
  grid runs `0 .. 1/(2*dt)` with dt taken from the arcs actually swept.

- **ARCS REMAIN INDEPENDENT.** Nothing here models correlation across arcs: a
  batch never spans two arcs (`BlockAccumulator.accumArc` splits ONE arc into row
  batches), each arc gets its own Cfcn/Wfcn, and there is no arc-to-arc term.
  This matches Ellmer exactly — his Q^m_ll is arc-local, Nmax x Nmax, and every
  arc's covariance is a leading principal submatrix of it.

### Documentation is JavaScript-free (hard environment requirement)
The user's environment forbids JS/Java. Final doc state:
- **PDFs** typeset from LaTeX (`docs_src/build_pdfs.py`, converts the authored
  HTML bodies to LaTeX, native math): `lsqvce_user_manual.pdf`,
  `lsqvce_reference_manual.pdf`, `lsqvce_interpretation_guide.pdf` in the root.
- **HTML manuals** with every equation a static PNG rendered from LaTeX at build
  time (`docs_src/build_manuals.py`, images cached by hash in `html/eq/`). Zero
  `<script>`, zero external refs across all 73 html files. Build is Python +
  pdflatex/pdftoppm/Pillow only. The old MathJax runtime AND the Node/JS converter
  are gone.

### Observability + plot-hook plumbing (most recent change; see APPENDIX for detail)
All gated on `Verbose`:
1. Outer loop prints explicit `OUTER PASS k/N` + `[SCALE]`/`[SOLVE]`/`[SHAPE]`/
   `[GATE]` banners with per-stage tic/toc timing.
2. Per-day accumulation and residual progress with elapsed time (fused-serial
   path in `BlockAccumulator.processSystem` and `residuals.m`).
3. Solve step logs `rcond(N)` before → prior/regularisation → `rcond(N_reg)`
   after, plus a COMPACT per-system local-parameter summary (nl, est range, mean
   sd) via `logLocalSummary`.
4. Terse `sigma0`/`dof` lines expanded with inline legends.
5. FIXED, DOCUMENTED plot-hook context `ctx.schema = 1` (fields only ever ADDED,
   never renamed): both stages carry `ctx.x`, `ctx.sd` (sqrt(diag(sol.Cxx))),
   `ctx.common`, `ctx.procOpt`, `ctx.model`, `ctx.sigmaG`, `ctx.R/...`; postModel
   ALSO carries `ctx.psd/psdFreq/acv/arOrder` (via `modelSpectrum`, PSD =
   σ²|B|²/|A|²). Locals are NOT in the hook (they go to the log). New test
   `test_plot_ctx_schema_fields`. PSD and rcond math validated in Python.

### Naming, options, parity — COMPLETE (see APPENDIX + §3 rename table)
- `SeparableProblem`→`adjustVCE`, `recover`→`adjustBlock` (both warn-once aliased).
  Six internal renames (no alias). Public Name/Value options are lowerCamelCase
  (case-insensitive parser). makeLoader terse options expanded to long forms with
  deprecated aliases. See §3.

### Function audit (B1–B5, ~70 functions) — the genuine fixes that STAND
- `decorrFactor` O(n³)→O(n²) broadcast; `robust` reweight O(nB²)→O(nB);
  `Accumulator`/`Combine` `diag(s)*R`→broadcast.
- `solveNonlinear` sd via triangular solve `sqrt(sum((R\eye).^2,2))` (was
  `sqrt(diag(inv(R'R)))`, which squares the condition number).
- `whittle` chol PD-check; `Accumulator.solve` rcond rank guard (chol alone misses
  ~1e-14 singular matrices).
- **redundancy Hutchinson-seam operator bug FIXED**: was `U\(N*(U'\V))`
  (= U⁻¹NU⁻ᵀ, trace ≈ 1.02·nc), corrected to `U'\(N*(U\V))` (= U⁻ᵀNU⁻¹ = I).
- ~60 local-function/method help gaps filled. All validated in Python.

### One regression was made and REVERTED
- Item-5 "ar2arma robustness" (pole-reflection projection + tighter `>=1-1e-12`
  trigger) REGRESSED the resonant round-trip (relErr 0.16 vs 1.3e-14). `ar2arma.m`
  and `test_signalCore.m` are REVERTED to the pristine original (guard `>=1`, warning id
  `:unstable`). The roadmap "ar2arma robustness" item is WITHDRAWN.

### Test-suite health
Six PRE-EXISTING latent bugs (not audit-introduced) that the user's real runs
surfaced, all fixed: `LogFid`/`logFid` case mismatch (25-test cascade),
`firePlots` no-output hook, `datum('B')`→`'datumB'`, and a four-layer
`test_homogenise_...` fix (missing covariance → rank-deficient design →
fragile lastwarn → separable-needs-Cep so the run legitimately can't complete, now
a try/catch that tolerates the expected error and checks lastwarn). That work
brought the suite to 393/393 at the time; see §1 for the current figure.

---

## 3. RENAME QUICK-REFERENCE

| old name | new name | alias kept? |
|---|---|---|
| `recover` | `adjustBlock` | yes (warn-once fn + HTML stub) |
| `SeparableProblem` | `adjustVCE` | yes (warn-once subclass + HTML stub) |
| `covModelFactory` | `covModel` | no |
| `fromAccum` | `diagnoseAccum` | no |
| `psdIterate` | `iterateNoiseModel` | no |
| `descaledResidualAcv` | `residualAcv` | no |
| `maxLagCorr` | `whitenessCorr` | no |
| `selectAROrder` | `selectArOrder` | no |
| option `datum 'B'` | `'datumB'` | sibling-aliased |
| makeLoader `obsTypes`/`cov`/`sysID`/`nl`/`common`/`local`/`procOpt`/`epochCov` | `observationTypes`/`covariance`/`systemId`/`numLocal`/`commonParams`/`localParams`/`processingOptions`/`epochCovariance` | yes (warn-once) |
| all PascalCase public options | lowerCamelCase | n/a (case-insensitive parser) |

Internal PascalCase (`BlockAccumulator`/`PrefetchLoader` `obj.opt.Field`) and the
internal `obsType` token / adjustBlock systems-spec struct fields were left as-is
ON PURPOSE (load-bearing `o.Field` struct access is case-SENSITIVE).

---

## 4. WORKFLOW AND VALIDATION HABITS

### 4a. THE INTERFACE (2026-07-18) — ONE NAME FOR EVERYTHING

**All aliases have been removed.** Pre-release, so backward compatibility was not
a constraint, and two ways to say one thing is precisely the ambiguity that
produced the `nc` duplication. Three families went:

1. **Fluent method aliases** on `adjustVCE` — renamed, old names deleted.
2. **Nine `makeLoader` option aliases** (`obsTypes`/`cov`/`sysID`/`nl`/`nc`/
   `epochCov`/`common`/`local`/`procOpt`) — deleted along with the `coalesceOpt`
   resolver and its `warnOnce` helper. `PartialMatching` is off, so an old name
   is now an `InputParser:UnmatchedParameter` ERROR, not a silent no-op.
   `testMakeLoaderShortAliasesRemoved` guards this: if PartialMatching were ever
   re-enabled, `'sysID'` would silently prefix-match `'systemId'` again.
3. **Two shim FILES** (`recover.m`, `SeparableProblem.m`) — to be deleted.

**There is no `nc` argument.** It is derived from the systems by
`LSQVCE.RESOLVENC`, so the common-parameter count is stated exactly ONCE, where
the system is defined:
- a Loader supplies it via `'commonParams'` (rows) or `'numCommon'` (scalar) —
  symmetric with `'localParams'`/`'numLocal'`, and passing both when they
  disagree is an error (`makeLoader:ncConflict`);
- a raw spec struct supplies it in its `.nc` field, beside the `.nl` it already
  carries.
All systems must agree; disagreement is `resolveNc:conflict`, not a silent pick.
WHY: the demos used to say `cfg.nc = 6` on one line and build a 6-row
`cfg.common` on another — the same number twice, free to drift apart.

```matlab
R = lsqvce.adjustVCE(systems, days) ...      % was (nc, systems, days)
       .groupComponentsBy('obsType') ...
       .stochasticModel('MaxOrder', 40) ...
       .varianceComponents('MaxIter', 8) ...
       .iteration('MaxOuter', 6) ...
       .run();
```

RENAMES (old names DELETED, not aliased — this is development, and two ways to
say one thing is the bug we keep finding):

| was | is | why |
|-----|----|----|
| `.vceGroupsBy()` | `.groupComponentsBy()` | says what is grouped |
| `.shape()` | `.stochasticModel()` | the covariance structure |
| `.scale()` | `.varianceComponents()` | the Helmert weights |
| `.outer()` | `.iteration()` | purpose, not implementation |
| `.axes()` | `.componentsPerEpoch()` | `axes` SHADOWED a built-in |

ORDERING NOTE worth keeping: in `adjustBlock`, `resolveNc` must run AFTER
`normalizeSystems` validates the shape. Deriving first made a malformed
`systems` report a confusing `resolveNc:missing` instead of the specific
`badCell`/`badStruct`.
**This note was written but the code was never changed** — the bug shipped and
was caught by `testRecoverErrors` on 2026-07-18 (defect 4 in §1.0). It is fixed
now, with the rationale inline in `adjustBlock.m` so the order is not "tidied"
back. LESSON: a rule recorded ONLY in the handover is not in force; put it in
the code, and give it a test.

### 4b. THE WORKING LOOP (REVISED 2026-07-18 — the zip round-trip IS the loop)

**THE ZIP ROUND-TRIP IS THE WORKING LOOP. Do not build a loop on direct file
access or on the live MATLAB connection.** An earlier version of this section
described Claude editing `/home/matweigelt/Matlab` in place via the Filesystem
tools and reading logs off the user's disk. That is SUPERSEDED. In practice
there is often no MATLAB and no filesystem access at all, and even when the
tools appear to work they are error-prone: the sandbox copy and the user's
copy drift out of sync silently, which is exactly how a build went unnoticed
and stale once already.

THE LOOP:
1. **The user uploads a zip** of the toolbox (or Claude ships one from the
   previous round).
2. **Claude works in the sandbox** on the unpacked tree, and runs
   `python3 docs_src/audit_tree.py` before shipping (§4f).
3. **Claude ships the ENTIRE toolbox as ONE zip** via `present_files` — never a
   partial set of changed files, because a partial drop is precisely what makes
   the two copies diverge. Every round ends with a shippable zip.
4. **The user unpacks it, runs MATLAB, and returns the result** either as an
   UPLOADED .txt FILE (preferred) or pasted into the chat.

RULES THAT MAKE THIS RELIABLE:
- **Prefer uploading the log as a .txt file over pasting.** `run_all_lsqvce`
  writes `tests/last_run_all.txt` and `run_lsqvce('<target>')` writes
  `tests/last_run_<target>.txt` precisely so there is always a file to attach.
- **LARGE PASTED MESSAGES GET LOST.** A long paste is copied on the user's side
  but does not always arrive — silently, with no error and no partial content.
  The full suite log is ~105 KB, far past that. So: attach the file. If the
  content must be pasted, paste a SMALL EXCERPT (the failure block, the summary
  table) rather than the whole log.
- **If Claude receives nothing, or an obviously truncated log, SAY SO and ask
  for the file again.** Do not infer results from a partial paste, and do not
  proceed as though a run happened. Check the `when :` timestamp on line 2 of
  the log: if it has not moved, the run is the previous one and the user is
  probably still on the old build.
- **The zip Claude ships is the reference.** If a MATLAB result contradicts the
  sandbox, suspect the user is running an older unpacked copy FIRST — confirm
  with `which -all lsqvce.adjustBlock` and a `dbtype` of the changed lines
  before doubting the code.
- Claude cannot run MATLAB. Everything in the sandbox is (a) STRUCTURAL, (b)
  NUMERICAL in Python, or (c) the static audit. **The user's MATLAB run is the
  ground truth** and has repeatedly caught what the sandbox cannot see.

### 4c. TEST RUNNERS AND LOGS
- `run_lsqvce('<target>')` runs a SUBSET and writes
  `tests/last_run_<target>.txt`. Targets: `joint`, `ellmer`, `contrib`,
  `separable`, `accum`, `core`, a file path, or a `'*pattern*'`. Per-target
  names mean successive runs do not overwrite each other.
  `run_lsqvce` with no argument LISTS the shortcuts; `run_lsqvce('ellmer',
  'Name','*psdVCE*')` narrows further, and `run_lsqvce('*degenerate*')` matches a
  test name anywhere. Warnings are MUTED by default there ('Quiet') because the
  suite deliberately provokes them in failure-mode tests and they otherwise bury
  the result — which is why the subset runner is for iterating, never for gating.
- `run_all_lsqvce` is the full gate and writes `tests/last_run_all.txt`. It does
  NOT mute warnings: the inventory is a signal (§4d).
- Both capture the FULL console output — warnings, failure diagnostics, stack
  traces — while the console shows only a summary. Send the FILE, not a paste
  (§4b).

### 4d. READ THE LOG, NOT JUST THE PASS COUNT
**Four times in one session, green tests concealed a real defect.** Every one was
caught by a QUANTITATIVE signal in the log, not by a failing assertion:
- 616 `psdUnderResolved` warnings while 4/4 passed → auto-sizing was not working
- log size 997 KB → 4.8 MB while 4/4 passed → a warning storm had reappeared
- 6 skip warnings where there had been 0 → the axis fix had a bug
- runtime 17 s → 74 s → the usability probe was factorising a dense matrix per
  stream per pass

So after every run, check: **warning counts by type, log size, and runtime** —
not only pass/fail. A test that cannot fail for the reason it names is worse than
no test; if a guard has no counter behind it, a failure there reports as an
all-zero tally (this is exactly how the axis-fix bug hid).

### 4g. RENAMING IS NOT DONE WHEN THE CODE COMPILES (2026-07-18)
A rename has THREE surfaces, and only the first is obvious:
1. **Executable code** — call sites. Caught immediately: it errors.
2. **Help text and `See also` lines** — MATLAB renders `See also LSQVCE.FOO` as a
   CLICKABLE LINK, so a stale one is a dead link in the shipped documentation,
   not a harmless comment. FIFTEEN of these survived the 2026-07-17 rename
   undetected, including in both demo files a user is told to copy.
3. **The generated HTML, `Contents.m`, `helptoc.xml`, `README.md`** — §4f
   checks these.

So: after any rename, grep the OLD name across `+lsqvce/ tests/ examples/
*.md html/` and read every hit. `audit_tree.py` check 2 now covers help text as
well as code, with an allowlist for the deliberate "was REMOVED" history notes.

### 4h. IF A CHECKER CANNOT SEE COMMENTS, IT WILL MISS THINGS (2026-07-18)
This has now caused THREE separate defects, so it is a standing rule rather than
an anecdote:
1. `log.render` returned a lone argument untouched, so `log.msg('TEXT\n')` wrote
   a literal backslash-n. Undetected because no test read the output.
2. The banner check first used `code_lines()`, which BLANKS string literals — so
   it could not see the very `'[resid ]'` text it was looking for, and reported
   a clean tree while 37 hardcoded banners sat in it.
3. The removed-name check had the same flaw for COMMENTS, so 15 stale `See also`
   references passed as clean.

THE RULE: a checker that strips comments or strings is right for BALANCE
COUNTING and wrong for CONTENT MATCHING. When looking for text, read the RAW
source and filter afterwards. When counting blocks, strip first. `audit_tree.py`
does both and says which it is doing at each site — keep that distinction
explicit if you extend it.

### 4e. STRUCTURAL / NUMERICAL CHECKS
- **Structural**: function/end balance vs the pristine original. Crude counters
  are fooled by `history(end+1)`, `Cprior(2:end)`, and single-line
  `if ..., ..., end`; a residual that MATCHES the original means balanced.
  MATLAB's own `checkcode('<file>')` is authoritative — use it when the counter
  disagrees with itself.
- **Casing**: parsed-option access must match its declaration's case. The NV
  parser is case-insensitive, but `o.Field` struct access is case-SENSITIVE.
- **Numerical**: validate any algebra in Python BEFORE trusting MATLAB. This
  session that caught: the DCT-I normalisation, the degenerate white start, the
  `2*(nf-1)` rank rule, `W = inv(Wop)'`, and that slicing `W` breaks `Sigma=W'W`.
- **Skills**: when creating files or running code, read the relevant
  `/mnt/skills/public/<name>/SKILL.md` FIRST.

---

### 4f. THE AUDIT SCRIPT (new 2026-07-18) — CHECK, DON'T TRUST

`python3 docs_src/audit_tree.py` (run from the toolbox root, READ-ONLY, exit 0 =
clean) re-derives from the tree the facts this handover asserts: test-file
registration vs `run_all_lsqvce`'s hardcoded list, removed names in executable
code, entry-point signatures at every call site, html-page + `Contents.m`
coverage for every function, link/helptoc integrity, the JS-free requirement,
and function/end balance. It exists because this document DID drift (§1.0) and
prose cannot be verified by reading it again.

RUN IT after any rename, any new function, and before shipping a zip. It is not
a substitute for MATLAB — it cannot execute anything — but every defect in §1.0
would have been caught by it immediately.

Two notes so its output is read correctly:
- The two redirect stubs (`lsqvce.recover.html`, `lsqvce.SeparableProblem.html`)
  report as "no backing .m". That is EXPECTED and intended — they are kept so
  old links land somewhere useful.
- Its balance check drops `case` LABELS before counting, because `case 'parfor'`
  and `case 'spmd'` embed block-opener keywords as string literals and a naive
  counter reports `mapreduce.m` +2 and `redundancy.m` +1. Those two were
  long-standing FALSE POSITIVES, not defects. If you write your own checker,
  reproduce this or you will chase them again.

## 5. KEY FILE LOCATIONS
- Engine: `+lsqvce/adjustVCE.m` (1400+ lines; run loop ~330–415, firePlots
  ~1240, plotContext + loaderParametrisation + modelSpectrum ~1277, applyShape
  ~1090, checkHomogenised PRIVATE ~801, logFid reads ~233).
- Accumulator: `+lsqvce/BlockAccumulator.m` (solve ~366–490 with rcond +
  logLocalSummary, processSystem fused loop ~706 with per-day timing,
  backSubstitute (sysID,x)->y, localCovariance (sysID,Cxx)->Cyy, commonBlocks/
  reducedNormalMatrix return the TOTAL normal).
- Residuals: `+lsqvce/residuals.m` (per-day loop ~126 with progress timing).
- Others: `+lsqvce/{adjust,adjustBlock,vce,varComp,vceHutchinson,covModel,
  makeLoader,redundancy,residuals,outliers,ar2arma,selectArAcv,armaACV,
  poolGroups}.m`, `+lsqvce/+util/log.m` (printf-style `msg`/`warn`; `warn` always
  raises a real MATLAB warning so `verifyWarning`/`lastwarn` work).
- Tests: `tests/run_all_lsqvce.m` (THE runner), `tests/test_lsqvce*.m`,
  `tests/test_{separable,reliability,precision,recover,v3,v4,v7,v8,v9,...}.m`.
- Messaging: `+lsqvce/+util/verbosity.m` — the LEVELS (`resolve`/`atLeast`/
  `progress`/`stage`/`diagnostic`) AND the layout helpers (`banner` for the flat
  log; `heading`/`dayLine`/`stepLine`/`groupLine`/`duration` for the indented
  stage detail). Formatters that consume it: `convergenceRow`, `sigmaLine`,
  `contribLine`, `groupKeyStrings` (file-locals at the bottom of adjustVCE.m).
- Docs build: `docs_src/build_manuals.py` (HTML+PNG eqs), `docs_src/build_pdfs.py`
  (LaTeX PDFs). Both Python + pdflatex only, no JS.
- Maintenance scripts (NOT shipped docs): `docs_src/audit_tree.py` (§4f, the
  consistency check) and `docs_src/build_missing_pages.py` (generates an html
  page for any function lacking one; only writes pages that do not exist, and it
  scans `+util` too).

---

# APPENDIX — dated session history (newest first)

### 2026-07-19 (evening session): full design review — v1 frozen, v2 scoped
Three-expert scrutiny (GRACE processing / adjustment theory / architecture)
of the complete tree against Ellmer (2018) ch. 2+6 and Kvas (2019). No `.m`
file was changed. Findings, verdicts and the successor design live in
`HANDOVER.md (this folder)` (§2 review, §3 design, §5 work packages WP-0..22
with implementation prompts, §6 test/doc guidelines). Headlines: the Ellmer
machinery (Phases A–E) is complete and verified but not the default; the
dense-Toeplitz batch whitening is a silent block-diagonal approximation
(`covModel` toeplitz path, `filterMemory = 0`) — v2 replaces it with an
exact streaming Levinson/Schur whitener; snooping exists only on the flat
path (open item B); `diagnoseAccum` rebuilds the full-parameter factor and
does not scale (default report route); pre-elimination is per-system, not
per-arc. The spectral-domain component formulation (Kvas eq. 3.14) was
identified as the clean cure for the degenerate white start documented in §2
KEY FACTS. This tree stays green as the reference implementation.

_The section above is the current state. Everything below is the blow-by-blow
history of how it got there, kept for provenance; you rarely need to read it._

### 2026-07-18 (session): streaming contribution primitive + §1.1/§1.2 shipped & confirmed
Implemented the coupled §1.1 (streaming Hutchinson VCE) and §1.2 Level 2
(contribution diagonal) via ONE shared primitive, plus §1.2 Level 1 earlier in the
same arc. New files: `+lsqvce/streamContribution.m` (matrix-free per-group
N^-1 N_g^cc probing over one data sweep) and `+lsqvce/contribution.m` (Level 1 exact
share + Level 2 streamed diagonal). Wired `vce`'s `'TraceMethod','hutchinson'` to a
REAL per-group trace substitution (was a machinery check). New test file
`tests/test_contribution.m` (8 tests) + registered it in `run_all_lsqvce`
(the suite list is hardcoded, not a glob -- a new file is invisible until added).
KEY CORRECTNESS FIX caught offline before shipping: the per-group reduced-common
contribution must use the ORDERING-INDEPENDENT symmetric Schur decomposition
`N_g^cc = Ncc_g - 1/2(Nlc_g' inv(Nll_s) Nlc_s + Nlc_s' inv(Nll_s) Nlc_g)`, NOT a
"first group seen owns the system's Schur term" rule -- the latter is arbitrary and
wrong whenever a system's arcs span several groups (exactly GRACE obsType grouping
of the merged C/D pair). The symmetric form partitions N exactly and is verified in
Python (sum_g tr = nc, sum_g diag = ones, both exact). One MATLAB round-trip found a
TEST-FIXTURE bug (covariance reweighting needed sqrt(variance) threaded through
opt.covOverride, mirroring adjustVCE.covOverrideOpt -- not a direct cfg.cov mutation
by sig); code was correct. Final: 404/405 suite green (lone remainder is the
pre-existing parallel-pool test_parallel filter), test_contribution 8/8. Only §1.3 (docs)
remains on the roadmap. ALSO added a "Relation to Ellmer (2018)" theory subsection
(REFMAN &sect;3a in docs_src/build_manuals.py) + Ellmer as reference [21]; rebuilt the
HTML manuals (equations render to PNG; JS-free). NOTE for future doc builds: this
sandbox's pdftoppm (poppler 24.02) rejects the -transp flag build_manuals.py uses,
so a cache-miss equation won't render here without dropping that flag; the committed
script is correct for the user's environment (transparent PNGs) and was NOT changed
-- the sandbox rebuild used a throwaway monkey-patch.

### 2026-07-17 (session end): documentation rebuild scoped (ReportLab master PDF)

User wants the docs rebuilt with ReportLab (prefers it over LaTeX/HTML), as ONE
master PDF (decided: single doc, six parts, not six PDFs) with an AUTO-GENERATED
API part (decided). Rewrote §1.3 as the full build sketch: new
docs_src/build_reportlab.py (Platypus), reusing _render_eq (LaTeX→PNG) for
equations and matplotlib for flowcharts/decision-trees. Six parts: (1) Theory
covering ALL approaches, (2) Workflow loop-flowcharts per approach, (3) User manual
= decision trees + step-by-step per approach, (4) expanded Interpretation guide,
(5) NEW Failure catalogue (enumerate every lsqvce: id -> meaning+fix), (6) API
auto-gen from the standardised help blocks + a curated per-function side-table
(position-in-chain / theory-link / best-practices). Verified tooling present:
reportlab 4.4.10, matplotlib 3.10.8, pdflatex/pdftoppm/Pillow. Doc usage unchanged;
JS/Java-free trivially. Updated §1.3 + the new-chat brief. Docs only; no code.
This is intended as the FINAL planned addition before the doc build itself.

User: show Level 1 (scalar share) EVERY outer pass, Level 2 (per-parameter
diagonal) only in the FINAL iteration. Folded into §1.2. VERIFIED the per-pass
Level 1 is genuinely free: `vce` already returns `vc.r` (per-group redundancy,
aligned with `vc.groups`), so in adjustVCE.run the share is `vc.r./sum(vc.r)` right
after the vce call — no extra data pass. Design: log it as a per-pass `[CONTRIB]`
line (Verbose-gated), label as CURRENT-pass (it moves, being weight-dependent),
AND store a `contribution` trajectory (nOuter×nGroup) in R.history next to sigma/
arOrder/white. Level 2 runs once at converged weights via the §1.1 streaming
primitive. Level 1 ships independently/immediately; Level 2 needs §1.1. Docs only.

User asked about a "contribution matrix" = information share C_s = N⁻¹ N_s^cc
(Σ_s C_s = I). Captured as new open-topic §1.2, three levels: L1 scalar share is a
QUICK WIN (tr(C_g)=tr(N⁻¹ N_g^cc) is ALREADY red.group(k).r; just normalise to
r_g/(m-u), sums to 1 -> new lsqvce.contribution wrapper); L2 per-parameter diagonal
diag(N⁻¹ N_g^cc) rides on the SAME streaming primitive as §1.1 (build once, two
consumers); L3 full dense C_g is O(G·nc²) — do NOT. Weight-dependent, compute at
converged weights. Cross-linked §1.1<->§1.2 and updated the new-chat brief. Docs
only; no code.

The user flagged that the "expose/hold per-group N_g^cc" plan does NOT scale
(O(G·nc²) memory; a group can span many arcs). Rewrote §1.1 with the memory-honest
design: never materialise N_g^cc — apply it to Z probe vectors by STREAMING each
group's arcs (mirroring the exact per-arc leverage path ALREADY in redundancy.m,
which builds leverage blocks once and re-streams arcs via dataFcn without forming
any N_g^cc), accumulating all groups' N_g^cc·Z in a SINGLE data sweep. Added a
CORRECTION note in the appendix's old item-0 bullets pointing forward to §1.1, and
updated NEW_CHAT_START_HERE.md. Docs only; no code changed.

Rewrote the top current-state section for clarity and put OPEN TOPICS front and
centre (§1): (1) streaming/block-arrow Hutchinson VCE — the one real engineering
gap; (2) additive doc polish; (3) optional small items; (4) user-only GRACE-FO
run; (5) the observability change is not yet MATLAB-confirmed. The dated APPENDIX
below is preserved unchanged. Added `NEW_CHAT_START_HERE.md` (a paste-in brief for
continuing in a fresh chat). No code changed this step.

User asked to make SHAPE<->SCALE explicit, add day/residual progress + timing,
show conditioning in the solve, explain the terse output, and pipe estimates/sds/
common/procOpt to hooks (+PSD at postModel). Decisions: rcond (not cond); locals
to the LOG not the hooks; "processors"=processingOptions; single on/off Verbose;
compact local summary. Five changes, all gated on Verbose:

1. **Outer-loop banners (adjustVCE.run).** Explicit `===== OUTER PASS k/N =====`
   plus `[SCALE]`/`[SOLVE]`/`[SHAPE]`/`[GATE]` stage lines, each with tic/toc
   timing and a one-line "what happens each pass" legend.
2. **Per-day progress (BlockAccumulator.processSystem fused path + residuals.m).**
   `[accum ] system S, day D : k arc(s) [t s]` during accumulation and
   `[resid ] system S, day D : k residual(s) [t s]` during the residual pass.
   Fused-serial path only (worker paths would interleave).
3. **Solve conditioning + local summary (BlockAccumulator.solve).** rcond of the
   reduced common normal BEFORE the prior and, if a prior regularises, AFTER
   (`rcond(N)=... -> common prior -> rcond(N_reg)=...`); validated in Python
   (reg raises rcond ~2.7e-17 -> 1.9e-5). New method logLocalSummary writes a
   COMPACT per-system local summary (nl, est range, mean sd) via backSubstitute
   + localCovariance. NB: backSubstitute is (sysID,x)->y ONE output; localCovar
   is (sysID,Cxx)->Cyy -- do NOT call backSubstitute with Cxx.
4. **Expanded terse output.** sigma0/dof lines now carry inline legends.
5. **Fixed documented ctx schema (plotContext, ctx.schema=1).** Both stages now
   carry ctx.x, ctx.sd (sqrt(diag(sol.Cxx))), ctx.common, ctx.procOpt (via new
   loaderParametrisation, matched by systemId), plus ctx.model/sigmaG/R/... .
   postModel ALSO gets ctx.psd/psdFreq/acv/arOrder via new modelSpectrum
   (S(w)=sigma^2|B|^2/|A|^2; validated == analytic AR(1) to 1e-15; white->flat).
   Locals are NOT in ctx (they go to the log per decision B). Header help +
   PLOTCONTEXT doc updated; fields only ever ADDED, never renamed.

New test: test_plot_ctx_schema_fields (saves one postSolve + one postModel ctx to
MAT, checks all schema-1 fields present, sd length=nc & >=0, procOpt is a struct,
and postModel psd/psdFreq align & psd>=0, arOrder>=0). All balance-checked; PSD
and rcond math validated in Python. Not yet run in MATLAB by the user.

### 2026-07-17 (session): sixth MATLAB run — layer 4 of the homogenise test

Same test, next error: `lsqvce:covModel:noCep` ("'separable' requires cm.Cep")
from buildAccumulator->covModel during run(). ROOT CAUSE (correct engine
behaviour, not a bug): a 3-D position stream is CORRECTLY promoted to the
separable model by the SHAPE step (adjustVCE.m:1108, `type='separable'`), and
separable REQUIRES the per-epoch D (cm.Cep). The test deliberately omits that D
(the whole point -- positions without D must warn), so the run legitimately CANNOT
complete: it warns in the step-0 preflight, then hard-errors downstream on the
very thing it warned about.

The prior fix used verifyWarning, which REQUIRES the call to return -- impossible
here. FINAL FIX: the test now try/catches empiricalStochastic (tolerating the
expected covModel:noCep error) and checks lastwarn holds notHomogenised. Verified
this is reliable: checkHomogenised (empiricalStochastic:587) fires the warning
immediately BEFORE run() (592), and run() errors in the FIRST SCALE step
(buildAccumulator) with NO warning emitted in between (notConverged only fires
after the loop, never reached), so lastwarn is not clobbered. homogenise defaults
true so the preflight runs. The earlier makePosDay sin-column (full-rank) change is
now moot for this test (errors before any solve) but kept as a cleaner design.

This test peeled through FOUR layers: (1) missing covariance arg, (2) rank-
deficient design, (3) fragile lastwarn capture, (4) separable-needs-Cep is a hard
error so the run can't complete. All resolved. Expected: 393/393.

### 2026-07-17 (session): fifth MATLAB run — same test, deeper error

`test_homogenise_warns_when_scalar_positions` got past the covariance fix but then
threw `lsqvce:qr:singularR` ("singular at parameter 2") from the block-arrow solve
inside empiricalStochastic->run->vce. TWO root causes in this pre-existing test,
both fixed test-side:

1. **Rank-deficient toy design.** `makePosDay` built local column B = 0.05*ones,
   which is COLLINEAR with the common constant column A(:,1)=ones. After the
   block-arrow Schur elimination of the local, the reduced common normal was
   numerically singular (Python: min eig 5.7e-14, cond 1.8e14 -> QR R(2,2)~0).
   FIX: B = sin((1:m)'/3) -- independent of the constant. Reduced common normal
   is then well-conditioned (min eig 7.9, cond 19). makePosDay is used ONLY by
   this test, so the change is isolated.
2. **Fragile warning capture.** The test used lastwarn() to grab the warning id
   AFTER empiricalStochastic returns, but the full run can emit a later warning
   (notConverged) that overwrites lastwarn -> false failure. FIX: use
   verifyWarning (detects the specific warning firing ANYWHERE in the call;
   log.warn always raises a real MATLAB warning so this works), and suppress
   notConverged so the tiny non-converging synthetic run stays quiet. verifyWarning
   also requires the call to COMPLETE, which the rank fix now guarantees.

Both validated in Python (Schur-complement conditioning) + structural balance.
Expected: 393/393.

### 2026-07-17 (session): fourth MATLAB run — 1 failure (test_homogenise)

Down to 1/393. `test_homogenise_warns_when_scalar_positions` errored with
`lsqvce:makeLoader:noCov` ("covariance has no field for observationType 'pos'").

Same class as the LogFid cascade: a PRE-EXISTING broken test unmasked once the
LogFid crash stopped aborting test_separable early. The test (identical in the
original) called makeLoader with `'observationTypes',{'pos'}` but NO
`'covariance'`, and makeLoader has always required a covariance entry per obsType
(the noCov check is identical in the pristine original). So the test could never
have completed against this makeLoader; it only surfaced now that the tests before
it pass.

FIX (test-side, matches its intent): supply `cov = struct('pos', struct('type',
'diag','sigma',1))` — a SCALAR (diag) covariance on the 3-axis pos stream, which
is exactly the "scalar positions" case the test targets. makeLoader now builds the
loader, and because the diag cov carries no epochCov hook, the step-0
homogenisation preflight fires `lsqvce:adjustVCE:notHomogenised` as the test
verifies. Confirmed the checkHomogenised logic (nAxesOf('pos')=3 via axes() and the
pos-prefix default; no epochCov key => missing => warn).

PROACTIVE sweeps done (no more of these lurking):
- No other makeLoader test call declares obsTypes without covariance.
- No remaining p.Results.<Case> vs addParameter('<case>') mismatches anywhere in
  +lsqvce (the LogFid one was the only instance).

Expected result: 393/393.

### 2026-07-17 (session): third MATLAB run — 2 failures (both test_separable)

Down from 31 -> 2. Both in test_separable:

1. **test_plot_hooks_fire (pre-existing bug, unmasked by the LogFid fix).**
   `firePlots` called `figH = hookFcn(ctx)`, requiring the hook to RETURN a value.
   A legitimate side-effect-only hook (records to a file, returns nothing) then
   throws MATLAB:TooManyOutputs, so the hook never records and the test's temp
   file is never written. This is in the ORIGINAL too; it only surfaced now
   because the LogFid crash used to abort the test before reaching firePlots.
   FIX: firePlots now requests one output but tolerates a no-output hook (catches
   TooManyOutputs/maxlhs/unassignedOutputs, calls the hook for its side effect,
   treats figH as []). This makes recording/logging-only hooks first-class.
2. **test_convergence_table_shape (my Item-3 test).** Made its column check robust
   (ismember instead of strict VariableNames equality). All the Item-3 code paths
   (applyShape returning pwhite, hist.white/whitePval written together each pass,
   ljungBoxWhiteP == ljungBoxWhite statistic, chi2Cdf present, convergence() table)
   were re-verified sound against the source; the strict-equality assertion was the
   only real fragility.

Both edits balance-checked. Note: the sandbox has no MATLAB, so these were found
by reading the traces + source; the user's real run is the ground truth.

### 2026-07-17 (session): second MATLAB run — regressions found & reverted

User ran `run_all_lsqvce` (395 executed tests, 31 failed). Accounting:

**Already fixed the prior round, but the user ran the PRE-FIX zip** (line numbers
in the log prove it): test_rank_deficiency (rcond guard), test_datum_user_supplied
(datumB), test_redundancy_seam (operator). These pass with the shipped build.

**NEW regressions my audit introduced — fixed this round:**
1. **LogFid cascade (25 test_separable errors).** adjustVCE registered
   `addParameter('logFid'/'logFlag')` but READ `p.Results.LogFid`/`LogFlag`
   (capital L). inputParser preserves the registered case in .Results, so every
   adjustVCE construction threw "Unrecognized field name LogFid". This was a
   PRE-EXISTING latent bug in the original (dormant), surfaced now. Fixed the
   three reads to `p.Results.logFid`/`logFlag`.
2. **ar2arma REVERTED (testAr2armaRoundTrip, testAr2armaProjectionPreservesAcv).**
   My Item-5 "robustness" change (pole-reflection projection + tighter `>=1-1e-12`
   trigger) REGRESSED the resonant round-trip the function exists for: relErr
   0.16 vs the reference 1.3e-14. Root cause: the tighter trigger fired on stable
   near-boundary fits and the reflect/refit path replaced a good ARMA with a worse
   one. Item-5 was a net negative -> ar2arma.m and test_signalCore.m REVERTED to pristine
   original (guard `>=1`, warning id `:unstable`, no reflection). Lesson: my
   Python pole-reflection proof was correct in isolation but broke integrated
   behaviour; only the real MATLAB run caught it.

**NOT bugs:**
- test_vceMath/testBlockSchurIsFaster: a flaky timing MICRObench (0.79 vs 0.60 ms on a
  tiny problem where LAPACK overhead hides the asymptotic block win). Pre-existing,
  untouched by the audit. Made robust: larger problem (nl=1600) + 1.5x margin so it
  measures the real asymptotic advantage, not noise.
- test_parallel/PrefetchServesFromCache: correctly SKIPPED (verifyAssumption) when no
  parallel pool — working as designed, not a failure.

**Count note:** runtests reports ~395 EXECUTED tests (parameterised helpers +
skipped/filtered differ from a raw grep of ~423 function defs). 395 is the real
executed number.

**ITEM-5 IS NOW REVERTED** — the roadmap "ar2arma robustness" item is withdrawn as
a regression. Everything else from the audit (B1-B5 fixes) stands.

### 2026-07-17 (session): first real MATLAB test run — 3 failures fixed

User ran the suite in MATLAB (invaluable — the sandbox has no MATLAB). Two points:

**Test COUNT clarification (not a bug):** `runtests('tests')` reports only ~175
tests because it does NOT recurse into `tests/` (239 more tests). There are
~423 test functions across 16 pass/fail files. The correct entry point is
`run_all_lsqvce` (lists all 16 files) or `runtests('tests','IncludeSubfolders',
true)`. README updated (was wrongly calling `runtests('tests')` "the full suite"
and had stale 24/34 counts).

**3 failures, all fixed:**
1. **REAL CODE BUG (redundancy Hutchinson seam).** `+lsqvce/redundancy.m:272` built
   the trace-check operator as `U\(Ncc*(U.'\V))` = U^-1 N U^-T, which is NOT the
   identity (a similarity transform with trace ~1.02*nc), so tr(N^-1 N) estimated
   ~102 instead of 100. The code comment already said the operator should be
   U^-T N U^-1 = I; the two solves were swapped. Fixed to `U.'\(Ncc*(U\V))`
   (verified in Python: trace exactly nc, per-probe variance ~1e-29). This is the
   experimental opt-in 'TraceMethod','hutchinson' path. test_traceEstimate updated
   to match the corrected operator.
2. **ROBUSTNESS GAP + FLAKY TEST (rank deficiency).** `test_rank_deficiency_errors_
   clearly` built A with col4=col1+col2 and expected solve() to throw. chol on the
   resulting A'A catches the singularity only ~60% of the time (smallest eigenvalue
   ~1e-14 > 0), so the test was seed-dependent. Fixed BOTH: (a) added an rcond
   guard in Accumulator.solve normal path — if chol succeeds but rcond(R)^2 < 1e-12
   the system is numerically rank-deficient and rankFailure() is raised (verified:
   200/200 singular caught, 0/200 false triggers on well-conditioned N incl. the
   cancellation-warning test at rcond^2=5.7e-4); (b) made the test deterministic
   with a fixed RandStream.
3. **TEST BUG (datum param name).** `test_datum_user_supplied_basis` called
   `lsqvce.datum(N,'B',Bu)` but the parameter is `'datumB'` (function correct; the
   bug was in the ORIGINAL upload too). Fixed the test.

All three fixes balance-checked; the two operator/threshold facts validated in
Python.

### 2026-07-17 (session): all JavaScript removed — PDFs from LaTeX, HTML equations as PNG

**Why:** the previous round left one JavaScript file, `docs_src/render_equations.js`
(a build-time Node/MathJax converter). Even build-time-only, JS is not allowed in
this environment. Removed entirely; the whole toolbox now contains ZERO .js/.jar/
.class/.java files and ZERO `<script>` tags.

**PDFs (rebuilt from LaTeX).** New `docs_src/build_pdfs.py` converts the three
authored HTML bodies (INTERP/REFMAN/USERMAN in build_manuals.py) to LaTeX and
compiles with `pdflatex` — headings, tables, lists, code, prose, and the native
`\(...\)`/`\[...\]` math passed through as REAL typeset equations. Outputs
`lsqvce_user_manual.pdf` (4p), `lsqvce_reference_manual.pdf` (3p),
`lsqvce_interpretation_guide.pdf` (2p) in the toolbox root. Verified valid +
visually checked (equations, tables, symbols all typeset correctly).

**HTML (equations as PNG images).** `build_manuals.py` now renders each equation
from LaTeX to a cropped transparent PNG under `html/eq/` (cached by content hash;
101 unique PNGs, ~408 KB) and embeds it as `<img class="eq-inline">` /
`<div class="eq-display"><img>`. Its MathJax `<script>` config was removed and CSS
for the equation images added. Result: no `<script>`, no SVG-as-math, no external
refs in any of the 73 html files.

**Build is Python + LaTeX only** (pdflatex, pdftoppm, Pillow) — no Node, no JS.
docs_src/README.md rewritten to document both build scripts. index.html now links
each manual's PDF next to its HTML; README doc table + equations paragraph updated.

**Verified on the packaged zip:** 0 js/jar/class files, 0 `<script>`, 0
javascript:/onerror, 0 broken equation-image refs, 0 broken links, 3 valid PDFs,
0 LaTeX build artifacts (.aux/.log/.tex) and 0 _pdfbuild/_tmp dirs shipped.

### 2026-07-17 (session): MathJax runtime removed — equations now static SVG

**Why:** the bundled `html/mathjax/tex-svg.js` (2.1 MB of third-party JS) plus its
CDN `onerror` fallback was a supply-chain / execution risk in a locked-down
environment. Removed entirely.

**What was done:** every equation in the three manuals was PRE-RENDERED to a static
inline `<svg>` at build time with `mathjax-full@3.2.2` (run only in the build
sandbox, never shipped), and the MathJax `<script>` config + async loader + CDN
fallback were stripped from the pages. Result, verified across all 73 html files:
ZERO `<script>` tags, ZERO external http refs, ZERO on* handlers. `html/mathjax/`
deleted. Equation counts rendered: user 46+2, reference 64+12, interpretation 6+0.
A real HTML parser confirms every (nested) `<svg>` is balanced; sample equations
rasterised and eyeballed correct. Manuals dropped from a 2.1 MB runtime dep to
652 KB of self-contained static SVG.

**Reproducible build (now TWO steps, see docs_src/README.md):**
1. `python3 docs_src/build_manuals.py` — emits manuals with LaTeX in \\(..\\)/\\[..\\]
   (no loader, no CDN — the generator's MathJax loader was removed too).
2. `cd docs_src && npm install mathjax-full@3.2.2 && node render_equations.js
   ../html/lsqvce_X.html ../html/lsqvce_X.html` for each manual — converts to static
   SVG and strips any leftover MathJax config. `node_modules` is build-time only,
   never shipped.

New file: `docs_src/render_equations.js` (the converter). README.md, docs_src/
README.md, build_manuals.py, and the manual footers updated to describe static SVG.

### 2026-07-17 (session): FULL FUNCTION AUDIT — in progress, batched

**Task:** scrutinise every function for logic/coding errors, memory/time
optimisation, robustness, and complete help text (size/type/description for every
in/out INCLUDING local functions). Batched with handover checkpoints so it can be
resumed if the chat runs out.

**Batch plan** (kernels → engines → support):
- B1 math micro-kernels: levinsonDurbin, armaACV, psd2acv, varACV, decorrFactor,
  whittle, armaWhiten, varWhiten, separableWhiten, residualAcv, whitenessCorr,
  selectArAcv, selectArOrder, ar2arma, checkPSD, gcvLambda
- B2 VCE/trace/redundancy: varComp, vceHutchinson, vce, traceEstimate, redundancy,
  leverage, leverageBlocks, poolGroups, ReducedContribution, applyRegularisation,
  CommonPrior, datum
- B3 diagnostics/stats/robust: diagnose, report, stat, snoop, outliers,
  systematic, robust, ftest, localDiagnostics, crossCorr, ellipse, residuals,
  collectResiduals
- B4 solve engines/accumulators: adjust, adjustBlock, solveNonlinear, Accumulator,
  NormalEqSystem, QRSystem, BlockAccumulator, BatchIndex, Combine, Results,
  diagnoseAccum, mapreduce
- B5 separable engine/loaders/util: adjustVCE, covModel, iterateNoiseModel,
  makeLoader, synthLoader, PrefetchLoader, resolveSatSpecs, normalizeSatList,
  satKeyOf, +util/*, aliases (recover, SeparableProblem)

**B1 — DONE (2026-07-17).** Findings/changes:
- `decorrFactor`: replaced `diag(1./sqrt(d)) * V.'` with column-broadcast
  `(1./sqrt(d)) .* V.'` — O(n^2) not O(n^3), no dense diagonal. Verified in Python
  (identical result; W'W = inv(C) preserved).
- `whittle`: PD check now uses `chol` instead of `min(eig(Sf)) <= 0` (cheaper,
  standard).
- `varWhiten`: added full help to local `blockToeplitz`.
- `selectArAcv`: added full help to local `whiteReturn`.
- Verified numerically: `varACV` Yule-Walker identity to 1e-16; `gcvLambda` GCV
  form; pole/ACV facts. Everything else in B1 (levinsonDurbin, armaACV, psd2acv,
  armaWhiten, separableWhiten, residualAcv, whitenessCorr, selectArOrder,
  checkPSD) reviewed and found correct with complete help incl. local functions.
- NOTED (not changed, deliberate): `selectArAcv` re-runs Levinson per order in its
  sweep (~41x the Levinson work of a single ladder) but the FFT `armaACV` per
  order dominates and can't be laddered without the polynomial; the sweep also
  breaks early at pStar. Left as-is to avoid risking a correct hot path.

**B2 — DONE (2026-07-17).** VCE/trace/redundancy core. Findings:
- Added full help to local functions: `vce/ternary`, `redundancy/emptyPartial`,
  `redundancy/makeKeyFcn`, `poolGroups/ensureCell`.
- Verified `traceEstimate`, `varComp`, `vceHutchinson` are reference-quality
  (matrix-free, unbiased, complete help) — no changes needed.
- Confirmed all `inv()` occurrences in `redundancy`/`leverage`/`ReducedContribution`/
  `CommonPrior`/`datum`/`applyRegularisation` are in COMMENTS/help only; the code
  correctly uses triangular solves (chol + backslash) everywhere. No wasteful
  explicit inverses in any hot path.
- `applyRegularisation` diag() uses (L'*diag(lam)*L, diag(diag(Nll)), V*diag(d)*V')
  are legitimate one-time small-block constructions, not hot-path — left as-is.
- Every B2 function has complete INPUT/OUTPUT (or class properties) documentation.

**B3 — DONE (2026-07-17).** Diagnostics/stats/robust. Findings:
- **Real optimization:** `robust` reweighted-source was O(nBatch^2) — it re-streamed
  from batch 1 to locate batch k on every diagnose call. Fixed by precomputing
  per-batch row offsets ONCE (`batchRowOffsets`) so each `reweighted(k)` is O(1) to
  locate; whole diagnose pass now O(nBatch). Verified in Python: identical fac
  indexing for every batch, complexity O(nB^2)->O(nB).
- Added full help to 11 local functions: `report/{printWrapped,emit}`,
  `stat/{normCdf,normInv,mustBeProb}`, `snoop/{scanWorst,accumB}`,
  `systematic/mark`, `robust/{callB,reweighted,struct2nv}` (+ new `batchRowOffsets`).
- Added explicit INPUT blocks to `ftest` and `ellipse` (sizes/types for solFull/
  solRed and D/pair).
- No `inv()`/`min(eig)`/hot-path diag-multiply issues; `diagnose` uses 'INPUTS'
  header (fully documented). No logic errors found.

**B4 — DONE (2026-07-17).** Solve engines/accumulators. Findings:
- **Numerical fidelity fix:** `solveNonlinear` computed parameter sigmas as
  `sqrt(diag(inv(solTry.R'*solTry.R)))` — forming inv(R'R) SQUARES the condition
  number (defeating the QR backend) and costs two dense products. Replaced with a
  single triangular solve `Rinv = R\eye; sigX = sqrt(sum(Rinv.^2,2))`, which is
  diag((R'R)^-1) exactly but at cond(R) not cond(R'R). Verified in Python (matches
  to 1e-9; cond 4e2 vs 1.6e5).
- **Optimization:** `Accumulator` (row-sign normalise) and `Combine` (reduced QR
  sign normalise) used `diag(s)*R` with s=+/-1 — replaced with broadcast `s.*R`
  (O(n^2), no dense diagonal). Verified identical.
- Added help to ~28 methods/locals: adjust/{callBatch,stripArgs};
  adjustBlock/{oneSpec,isLoaderObj,loaderFcn,loaderDays};
  solveNonlinear/{callBatch,evalVPv,struct2nv,ternaryStr}; and one-line help to the
  private helpers/constructors of BlockAccumulator, BatchIndex, Combine, Results,
  Accumulator, NormalEqSystem, QRSystem (incl. Combine.redundancy's multi-line sig).
- `diagnoseAccum` and `mapreduce` reviewed — fully documented, no changes.
- No `inv()`/`min(eig)` issues remain in code (all such were comments).

**B5 — DONE (2026-07-17).** Separable engine/loaders/util. Findings:
- Added help to file-local functions: `adjustVCE/ternary`; `resolveSatSpecs/
  {emptyTemplate,matchById,override,lookup}`; `makeLoader/{requireHandle,dispatch,
  metaOf,arcVec,dayOf,validateDay,batchOf}`; `PrefetchLoader/{constructor,delete,
  harvest}`.
- No `inv()`/`min(eig)`/hot-path diag issues. covModel, iterateNoiseModel,
  synthLoader, normalizeSatList, satKeyOf, recover, SeparableProblem, and +util/*
  reviewed — fully documented, correct, no changes.

### AUDIT COMPLETE (B1–B5). Whole-toolbox verification:
- **Zero balance mismatches** across every modified file vs the pristine original
  (every edit preserved function/end structure).
- **Zero code-level inv() anywhere** in +lsqvce — every solve uses triangular/
  Cholesky factorisation (the solveNonlinear fix removed the last one).
- Zero broken HTML links.
- Real fixes made: decorrFactor O(n^3)->O(n^2) (broadcast); whittle chol PD check;
  robust O(nBatch^2)->O(nBatch); solveNonlinear cond-squaring inv(R'R) -> triangular
  solve; Accumulator/Combine diag(s)*R -> broadcast. ~60 local-function/method help
  gaps closed across all batches. All numeric claims validated in Python.
- Nature of the audit: the prior code was already high quality, so this pass
  CONFIRMED correctness and filled optimisation/help gaps rather than finding major
  bugs. The genuine wins are the five performance/fidelity fixes above and complete
  help coverage (incl. every local function).


### 2026-07-17 (session): demos, docs→HTML, and roadmap items 2–5

**Item 1a — Hutchinson second scenario in demos.** Real side-by-side exact-vs-
Hutchinson runs added where the FLAT `varComp` path allows a genuine substitution
(`demo_lsqvce`, `demo_workflow`), Z=4000, framed as the large-u method. Honest
pointer-notes (not fake blocks) added to the block-arrow demos (`demo_sensors`,
`demo_separable`, `example_vce_regularisation`) since Hutchinson is not a real
substitution there yet (block-arrow variant deferred by user — a major redesign).
The three demos that run no VCE were left untouched. New test
`test_hutchinson_recovers_planted_ratio_like_demos`.

**Item 1b — build-your-own-case guide.** `demo_workflow.m` fully replaced with a
from-scratch recipe (batchFcn contract → engine decision tree → common vs local
→ exact vs Hutchinson → reading the result), all CURRENT names.

**Item 1c — decision/workload diagrams.** New doc page `html/lsqvce_choosing.html`
with two static SVGs (`lsqvce_cost_vce.svg` cost curves, `lsqvce_decision.svg`
flowchart), linked from index/helptoc/README.

**Docs: PDF → HTML (all four).** The four root PDFs were REMOVED and their content
moved to HTML so the whole doc set travels with the toolbox and survives handover:
`html/lsqvce_user_manual.html`, `_reference_manual.html`,
`_interpretation_guide.html`. The API reference is the per-function HTML pages
(already current) plus the index. Equations were RE-TYPESET fresh in MathJax and
render OFFLINE via a bundled self-contained SVG build at `html/mathjax/tex-svg.js`
(CDN fallback only if that copy is removed). All names updated (adjustBlock/
adjustVCE, diagnoseAccum, etc.). The generator is stashed at
`docs_src/build_manuals.py` (re-runnable). index.html/helptoc.xml/README.md/
lsqvce_setup.m updated to point at the HTML; no shipped file references a .pdf.

**Item 2 — worker-aware logging.** `lsqvce.util.log` gained `startParallel()` /
`stopParallel()` (a `parallel.pool.DataQueue` sink) and a `queue` state key;
`msg`/`warn` route worker output to the client via `sendToQueue`/`drain`. Wired
into `mapreduce`'s parfor pass (safe no-op without the PCT). +4 tests in test_util.

**Item 3 — convergence diagnostics.** `adjustVCE` gained a `convergence(R)` method
returning a per-pass table (maxSigmaOffset, nWhite, minWhitePval, converged) and
an optional 2-panel plot (SCALE gate + SHAPE whiteness p-values). History now also
records `.whitePval` per pass (captured from `ljungBoxWhiteP` in `applyShape`).
+2 tests in test_separable.

**Item 4 — VCE group-pooling helper.** NEW `lsqvce.poolGroups(vc, ...)`: reads a
VCE result, finds groups below a redundancy floor, returns a remap (map/remap/
indexRemap/newGroups) with 'mergeLow' and 'nearest' strategies. It SUGGESTS/builds
a coarser grouping; the caller re-runs VCE. HTML page + Contents/index/helptoc.
+5 tests. Logic validated offline in Python.

**Item 5 — ar2arma robustness.** The stationarity guard now PROJECTS unstable AR
poles inside the unit circle (r → 1/conj(r), monic re-form — validated to preserve
the ACV with s2 UNCHANGED to ~1e-14) and re-fits the MA, accepting the projected
model only if stationary AND no worse than the pure-AR fallback on the ACV round-
trip; it falls back to the input AR only if projection cannot rescue it. Warning
ids: `lsqvce:ar2arma:projected` (rescued) and `lsqvce:ar2arma:nonStationary`
(fallback; replaces the old `:unstable`). +3 tests in test_signalCore.

**Still OPEN / deferred:** block-arrow & streaming Hutchinson substitution
(`vce`/`Combine.varComp`) — a genuine accumulator core change, deferred by user;
real GRACE-FO data run (item 1); namespace sub-packaging (item A leftover).



_Addendum, 2026-07-17 (item 0 tier T2 — varComp 'backend' dispatch). Model: Claude Opus 4.8._

**Roadmap item 0, tier T2 done, in the only form that is cleanly buildable:
`varComp` gained an opt-in `'backend','exact'|'hutchinson'` that delegates the
Hutchinson backend to `lsqvce.vceHutchinson`.** NOT run on R2026a.

### Why NOT the literal "wire into vce" (important scoping result)
The roadmap said "wire T1 into `vce`". On inspection that is NOT cleanly buildable
without the accumulator core change the roadmap itself deferred, because the three
VCE consumers have DIFFERENT architectures:
- flat `varComp` -- assembles per-group Ng{g} by re-streaming -> `vceHutchinson`
  already IS this (T1). A backend swap here is clean.
- block-arrow `vce` -- gets r_g from `lsqvce.redundancy`, which uses per-OBS
  leverage on the TOTAL reduced normal (`acc.reducedNormalMatrix()`); the
  per-group N_g^cc is never formed. Hutchinson needs per-group N_g^cc -> a core
  change. redundancy's own 'hutchinson' path remains only the machinery check.
- `Combine.varComp` -- computes r_g from the FULL Qxx per observation (common +
  local + cross terms via QCi/Qii), not a clean tr(N^-1 N_g); Hutchinson would
  need re-deriving for the full block-arrow Qxx.
So T2 was built where it fits WITHOUT a core change: on the flat `varComp`, whose
signature is identical to `vceHutchinson`. (User agreed to this deviation.)

### What was added (varComp.m, content-only -- no new files)
- `'backend','exact'|'hutchinson'` (default exact). When 'hutchinson' AND no
  'checkD', varComp delegates the WHOLE job to `lsqvce.vceHutchinson(u, batchFcn,
  nBatch, groupFcn, nGroup, ...)` forwarding method/maxIter/tol/verbose/minObs/
  probes/seed/reprobe, and returns its (subset) struct. Dispatch sits AFTER parse,
  BEFORE the contamination guard and the main loop; it `return`s on delegation.
- Guard: if 'checkD' is supplied, it CANNOT be honoured by the estimator (no
  per-obs redundancy), so varComp warns `lsqvce:varComp:hutchNoCheckD` and FALLS
  BACK to the exact backend. This keeps the contamination protection intact.
- New forwarding options `'seed'` (falls back to 'traceSeed') and `'reprobe'`.
- 'backend' is DISTINCT from the existing 'traceMethod': traceMethod only adds a
  machinery check to the EXACT backend; backend chooses the estimator itself. Both
  documented, distinction spelled out in help + HTML.
- Tests (test_lsqvce.m +3): delegates-and-marks (.estimator='hutchinson'), matches
  exact within the band (RelTol 0.05, Z=4000), checkD-falls-back-to-exact (result
  has no .estimator). HTML varComp page: 'backend'/'seed'/'reprobe' documented,
  See also -> vceHutchinson.

### Struct-shape caveat (document for users)
With `backend='hutchinson'`, the returned VC is vceHutchinson's SUBSET struct
(.s/.scale/.sigma/.redund/.redundSE/.estimator, NO .x/.significant/per-obs). A
caller that reads VC.x or VC.significant must stay on the exact backend. This is
intentional (strict Ellmer) and stated in the help.

### Still open under item 0
- The real block-arrow substitution (vce / Combine.varComp) still needs per-group
  N_g^cc (or full-Qxx probe re-derivation) = the accumulator core change. That is
  the streaming variant, still parked.
- Streaming Hutchinson that never forms Ng densely (the big-memory win) likewise
  needs the accumulator to expose per-group contributions.

**[CORRECTION 2026-07-17 — supersedes the two bullets just above.]** The user
flagged that "expose/hold per-group N_g^cc" does NOT scale (O(G·nc²) memory; a
group can span many arcs). The corrected, memory-honest design is in the CURRENT
STATE §1.1 at the TOP of this file: never materialise N_g^cc; apply it to Z probe
vectors by STREAMING each group's arcs (mirroring the EXACT per-arc leverage path
already in redundancy.m), accumulating all groups' N_g^cc·Z in a SINGLE data sweep.
Read §1.1, not these bullets, before starting item 0's streaming variant.

Roadmap §7 item 0: T1 DONE (vceHutchinson), T2 DONE (varComp backend dispatch,
flat). Block-arrow/streaming variants remain open (see §1.1 for the corrected
approach).

---

_Addendum, 2026-07-17 (item 0 tier T1 — lsqvce.vceHutchinson). Model: Claude Opus 4.8._

**Roadmap item 0, tier T1 done: a NEW standalone entry point `lsqvce.vceHutchinson`
implements the full per-group Hutchinson variance-component estimation, strictly
following Ellmer (2018) §6.5.6.** It is a SEPARATE door, so exact `varComp` /
`adjustVCE` are byte-for-byte untouched. NOT run on R2026a, but the trace algebra
and the full end-to-end algorithm were validated numerically offline (see below).

### Why a separate door (the design the user chose)
The estimator produces only PER-GROUP scalar traces; it structurally cannot return
per-observation r_i, so it cannot feed snooping / (w,r)/MDB / the sum(r_i)=m-u
guard. Bolting it into the exact path as a flag would have been a lie. Instead
`vceHutchinson` mirrors `varComp`'s signature and iteration EXACTLY, and deviates
only where Ellmer forces it: it drops the dense inverse and all per-obs outputs,
and documents that loudly.

### Scoping result that made this feasible WITHOUT a core change
The roadmap feared the streaming accumulator couldn't expose per-group N_g^cc.
On inspection, the FLAT path (`varComp`, and hence `vceHutchinson`) already
assembles each group's normal `Ng{g}` itself by re-streaming through
batchFcn/groupFcn -- it never needed the accumulator to hand it N_g^cc. So the
flat T1 path needs NO accumulator change. (The STREAMING block-arrow variant --
never forming Ng densely -- IS where the accumulator blocker bites; that remains
future work, and is the only place the O(u^2) dense-Ng assembly could be removed.)

### The algebra (validated offline in Python vs exact)
With N = U'U and a raw u x Z Rademacher block, Ubar = U^-1 (raw/sqrt(Z)):
- redundancy  tr(N^-1 N_g)          = tr(Ubar' N_g Ubar)      [Ubar on BOTH sides]
- Helmert     tr(N^-1 N_i N^-1 N_j) = sum(sum(W_i .* W_j)),  W_g = U^-T (N_g Ubar)
Ubar and each W_g are formed ONCE per iteration and reused for every group pair,
so no u x u inverse and no O(u^3) product appears. Offline check: H matched exact
to 0.1%; the redundancy trace matched to 0.05% at high Z; the FULL loop (with the
H s = q solve) reproduced exact varComp's sigma/scale/redund to ~1e-7 on a 2-group
problem, and both recovered a known sigma ratio of 2.0 identically.
IMPORTANT BUG CAUGHT DURING VALIDATION: the redundancy trace must use Ubar on both
sides; an earlier draft used Zbar on the left, which was wrong by a factor ~Z.
Fixed and re-validated. (Recorded so it is not reintroduced.)

### Cost / when it helps (be honest with users)
Cost O(u^3 [one Cholesky] + G*Z*u^2) vs varComp's O(G*u^3); crossover u ~ Z, so
the speedup is ~u/Z and only real for u in the THOUSANDS. Below a few thousand
unknowns varComp is faster AND exact AND gives the diagnostics -- the help and the
HTML say so in a warning box. It also trades accuracy for speed (MC error ~1/sqrt(Z),
worse on low-redundancy / near-rank-deficient groups). The dense per-group Ng
assembly (O(G*u^2)) is still paid here; only the inverse and the O(G*u^3) Helmert
products are removed. Removing the dense Ng too is the streaming variant (future).

### Files
- NEW `+lsqvce/vceHutchinson.m` (full help; strict-Ellmer; VC struct is a subset
  of varComp's, marked `.estimator='hutchinson'`, adds `.redundSE`).
- NEW `html/lsqvce.vceHutchinson.html` (with the when-to-use warning and the
  strict-Ellmer "does NOT produce" danger box). Registered in `Contents.m`,
  `helptoc.xml`, `html/index.html`.
- Tests: `tests/test_lsqvce.m` +5 (shape / no-per-obs-fields; matches-varComp
  within the sampling band [tier T3 validation]; simplified matches; deterministic
  with seed). Tolerances verified generous by offline simulation (Z=4000 -> ~1e-7).

### Watch on the box
- Untested on real MATLAB; the offline Python validation used the identical
  algebra, but MATLAB's chol gives UPPER U with N=U'U (the code assumes this;
  `x = U\(U'\b)` and `Ubar = U\Zbar` rely on it). If any assertion trips, first
  confirm the chol convention.
- Z=100 (default) is Ellmer's, tuned for well-conditioned operators averaged over
  many arcs; on a single low-redundancy group it may be too few -- the tests use
  Z=4000 for tight agreement. Real use should profile Z vs the redundSE it reports.
- tier T2 (wire vceHutchinson into `vce`/`adjustVCE` as an opt-in SCALE backend)
  is NOT done -- deliberately, to keep this a clean standalone first. That is the
  natural next step if profiling shows it earns its place at scale.

Roadmap §7 item 0: tier T1 DONE (standalone `vceHutchinson`); T2 (backend wiring)
and the streaming block-arrow variant remain open.

---

_Addendum, 2026-07-17 (item C part 2 — makeLoader option long-forms). Model: Claude Opus 4.8._

**Item C's second half done: makeLoader's terse option names got canonical
long-forms, with the old terse names kept as deprecated aliases (accept-both, warn
once).** This is the verbose-word expansion the roadmap parked. Unlike the case
pass, these are WORD renames — they do NOT match case-insensitively — so they need
real dual-accept mechanics. NOT re-run on R2026a.

### The option map (canonical <- deprecated alias)
| canonical | alias (deprecated) |
|---|---|
| `observationTypes` | `obsTypes` |
| `covariance` | `cov` |
| `systemId` | `sysID` |
| `numLocal` | `nl` |
| `commonParams` | `common` |
| `localParams` | `local` |
| `processingOptions` | `procOpt` |
| `epochCovariance` | `epochCov` |
| `arcs`/`days`/`name`/`cache` | (unchanged — already clear) |

### Mechanics
- makeLoader's inputParser declares BOTH names per pair (same default). After
  parse, a small `coalesceOpt` helper picks whichever the caller supplied
  (canonical wins if both) and writes it back onto the SHORT internal field name,
  so ALL downstream code (`o.obsTypes`, `o.cov`, `L.sysID`, ...) is unchanged.
  Using a deprecated alias warns ONCE per session (`warnOnce`, persistent Map);
  passing both warns about the clash. Verified the four cases
  (canonical-only / alias-only / both / neither) by offline simulation.
- **`p.PartialMatching = false`** was REQUIRED and set: three aliases are PREFIXES
  of their canonical names (`cov`⊂`covariance`, `common`⊂`commonParams`,
  `local`⊂`localParams`), so the default prefix matching would make them ambiguous
  and error. Exact matching keeps both spellings unambiguous. (Side effect: option
  abbreviations like `'obs'` no longer auto-complete — acceptable and stricter.)
  CaseSensitive stays default (off), so casing still doesn't matter.

### Scope boundary (what did NOT change — important)
- **Output/struct fields stay terse**: `L.nc`, `L.obsTypes`, `L.sysID`, and the
  ctx field `ctx.obsType` are the INTERNAL token the roadmap said to keep. They
  are read by 100s of consumers (stream keys `satKey|obsType`, groupBy, ctx). Only
  the makeLoader OPTION was expanded; the field names are unchanged. So a loader
  built with `'observationTypes'` still exposes `L.obsTypes`.
- **The adjustBlock/recover systems-spec struct** `struct('sysID',1,'nl',4,
  'loader',L)` is a SEPARATE interface (its field names are read as `sys.sysID`/
  `sys.nl` by the accumulator). Those field names were NOT renamed — expanding
  them is a distinct change. Tests that build these specs were deliberately left
  as `sysID`/`nl`.
- **`nc` is not an option** (it is derived, `L.nc = size(commonParams,1)`), so
  there is no `numCommon` option; `nc` remains the output field name.

### Docs/tests/examples
- Help block + error messages in makeLoader rewritten to canonical names (alias
  noted in parentheses). HTML page: option table shows canonical names with
  `(alias: ...)` and a deprecation note; example updated.
- Examples: `demo_separable_data.m` (the only example calling makeLoader) migrated
  to canonical names — teaches the new spelling.
- Tests: all makeLoader calls in `test_adjustBlock`/`test_separable` migrated to
  canonical names (span-scoped edit that left the systems-spec structs alone), so
  the suite runs WARNING-FREE. Added 3 alias tests to `test_adjustBlock`
  (`testMakeLoaderCanonicalNames`, `testMakeLoaderDeprecatedAliasWarns`,
  `testMakeLoaderCanonicalWinsOverAlias`).

### Watch on the box
- Old terse option calls in USER code keep working but now print a one-time
  deprecation notice per option per session. To silence during a migration:
  `warning('off','lsqvce:makeLoader:optRenamed_<alias>')`.
- If a user relied on inputParser ABBREVIATIONS for makeLoader options, those no
  longer work (PartialMatching off). Full names (either spelling) are required.

Roadmap §7 item C: **now fully DONE** (case convention + verbose-word expansion).
The whole naming cluster A/B/C is complete; only the OPTIONAL namespace
sub-packaging remains from A.

---

_Addendum, 2026-07-17 (item C — Name/Value option convention: lowerCamelCase). Model: Claude Opus 4.8._

**Roadmap item C done: all public Name/Value options are now lowerCamelCase, and
the duplicate group-key names are reconciled.** Behaviour is unchanged and OLD
SPELLINGS STILL WORK: MATLAB's `inputParser` is case-insensitive (no parser sets
`CaseSensitive=true`), so `'Verbose'` still matches the canonical `'verbose'` —
hence NO deprecation aliases were needed (user decision). NOT re-run on R2026a.

### The convention
lowerCamelCase for every user-facing Name/Value option: `verbose`, `maxIter`,
`tol`, `alpha`, `sigma0`, `backend`, `reg`, `parallel`, `groupBy`, `traceMethod`,
`probes`, `traceSeed`, `numProbes`, `showSeconds`, `exportDir`, `maxOuter`,
`minRedundancy`, `homogenise`, `validate`, `report`, `collect`, `whitened`,
`detrend`, `maxLag`, `orderGrid`, `acvTol`, `threshold`, `plot`, `distribution`,
`size`, `seed`, ... "lowercase the first letter".

### What changed (canonical names)
- **Declarations** in every public entry point lowercased (adjust, adjustBlock,
  adjustVCE, vce, varComp, redundancy, traceEstimate, iterateNoiseModel, outliers,
  residuals, datum, crossCorr, selectArOrder; diagnose/snoop/report/robust/
  systematic were already lowercase). Each matching `o.Xxx`/`p.Results.Xxx` access
  lowercased in lockstep.
- **Help text** and **HTML pages** for those entry points updated (whitelist
  driven — only option tokens, never values like `'qr'`/`'satObs'`).
- **Examples** (all `demo_*`, `examples/accum/example_*`, gravity example) updated
  to MODEL lowerCamel.

### Duplicate-concept reconciliation
- **`vceGroups` -> `groupBy`** on adjustVCE (now matches vce/redundancy/outliers).
  The fluent `.vceGroupsBy(...)` method and the private `vceGroups` property are
  unchanged; no `'vceGroups'` NV literal existed, so nothing breaks.
- **`vce`'s `groups` kept** — a distinct concept (explicit key LIST) vs `groupBy`
  (the STRATEGY); both coexist intentionally.
- **`datum` `'B'` -> `'datumB'`** — strict lowercasing gives a bare `'b'` (unclear;
  breaks the geodetic B-is-the-datum-matrix convention). `'datumB'` matches its
  siblings adjust/robust. Value/variable stays capital B.

### Fluent builder options (.shape/.scale/.outer)
These merge keys into internal PascalCase structs (shapeOpt/scaleOpt/outerOpt)
whose field names are load-bearing (read as `obj.scaleOpt.MaxIter`, asserted by
test_separable). Instead of renaming fields, **`mergeNV` was made CASE-INSENSITIVE**
so `.scale('maxIter',8)` and `.scale('MaxIter',8)` both hit MaxIter. Help documents
lowerCamel; internal struct fields stay PascalCase; the test is unaffected.

### Deliberately left PascalCase (internal)
**BlockAccumulator, PrefetchLoader** options (Backend, Verbose, KeepLocal, ...):
internal classes whose parsed options are stored in `obj.opt` and read as
`obj.opt.Backend` across many methods. Lowercasing the declaration would require
chasing every stored-field read (high risk, no user benefit — public functions
pass these via case-insensitive literals). Left as-is. (I briefly lowercased
BlockAccumulator, caught the `obj.opt.Parallel` mismatch, reverted — noted so the
next session doesn't re-do it.)

### Tests
Tests still pass options in PascalCase and KEEP WORKING via the case-insensitive
parser — not bulk-edited (~350 occurrences, higher risk than reward, and they
exercise the case-insensitive path, which is a feature). `test_separable` reads
`P.scaleOpt.MaxIter` — still valid. Migrating tests to lowerCamel is a clean
optional follow-up.

### Watch on the box
Old-cased option calls in USER code keep working. Genuine word-renames are only
`vceGroups->groupBy` (no literal usage existed) and `datum 'B'->'datumB'`
(sibling-aliased). Do NOT set `p.CaseSensitive=true` anywhere without migrating
every call site first.

Roadmap §7 naming cluster A/B/C: **all three DONE**. Only the OPTIONAL namespace
sub-packaging (deferred half of A) remains.

---

_Addendum, 2026-07-17 (top-level naming sweep incl. recover->adjustBlock — extends roadmap item A). Model: Claude Opus 4.8._

**Six awkward top-level names renamed (HARD renames, no aliases — these are
less-used internals), plus the AR/Ar casing inconsistency fixed.** This extends
the item-A class rename (SeparableProblem -> adjustVCE) to the rest of the public
surface. Behaviour is byte-for-byte unchanged; pure renames. NOT re-run on the
live R2026a box — run `results = run_all_lsqvce();`.

### The renames

| Old | New | Why |
|---|---|---|
| `covModelFactory` | `covModel` | "Factory" is a Java-ism, not MATLAB idiom; it returns a cov model |
| `fromAccum` | `diagnoseAccum` | preposition-first was cryptic; new name says it bridges accumulator -> diagnostics |
| `psdIterate` | `iterateNoiseModel` | verb-first; and it iterates the whole noise model, not just a PSD |
| `descaledResidualAcv` | `residualAcv` | shorter; the D^-1 descaling is the default behaviour, not the identity |
| `maxLagCorr` | `whitenessCorr` | names the PURPOSE (whiteness check), not the mechanism |
| `selectAROrder` | `selectArOrder` | **casing fix** — its ACV-domain sibling was already `selectArAcv` (Ar), so `AROrder` (AR) was inconsistent |
| `recover` | `adjustBlock` | completes the `adjust*` family: `adjust` (flat) / `adjustBlock` (streaming block-arrow) / `adjustVCE` (iterative + VCE). "recover" was gravity-field jargon that did not describe the operation |

**`recover` is the exception to "no aliases":** it is a primary documented DOOR,
not a low-use internal, so (per user) it got the FULL deprecation treatment like
`adjustVCE` — a warn-once function alias (`+lsqvce/recover.m` forwards to
`adjustBlock`, warns `lsqvce:recover:renamed` once per session) AND a redirect
stub at `html/lsqvce.recover.html`. The other six were hard renames. The private
members of `adjustVCE` that referenced recover were also renamed for coherence:
`recoverOpt`->`blockOpt`, `recoverWith`->`blockSolveWith`,
`recoverSpecs`->`blockSpecs`. Two alias tests added to `test_adjustBlock.m`
(`testRecoverAliasForwardsToAdjustBlock`, `testRecoverAliasWarnsRenamed`); the
existing `testRecover*` tests now construct `adjustBlock` and assert
`lsqvce:adjustBlock:*` ids.

### Deliberately NOT renamed (assessment outcome)

Reviewed all 62 public names; kept the rest on purpose:
- **`stat` kept** (user decision): renaming to `dist` would touch ~44
  `lsqvce.stat.<method>` call sites across 7 source files for little gain.
  Revisit under the namespace sub-packaging pass if desired.
- **`Combine` kept**: it is a `classdef ... < handle` (a noun class). A verb like
  `combineDatasets` would mislabel a class as a function. Not awkward enough to
  justify inventing `Combiner`.
- **`selectArAcv` kept** (it already had the correct `Ar` casing; the fix was to
  its sibling).
- **Literature/kernel names kept**: `levinsonDurbin`, `ar2arma`, `armaACV`,
  `armaWhiten`, `whittle`, `psd2acv`, `varACV`, `varWhiten`, `separableWhiten`,
  `decorrFactor`, `crossCorr`, `checkPSD`, `gcvLambda`, etc. These are precise
  domain/literature terms; "friendlier" names would HURT a geodesist who knows
  the term.
- **Verb-first doors already good**: `adjust`, `adjustVCE`, `recover`,
  `makeLoader`, `diagnose`, `report`, `snoop`, `outliers`, `robust`,
  `redundancy`, `vce`, `varComp`, `traceEstimate`, `datum`, `ftest`, `ellipse`,
  `systematic`, `residuals`, `leverage`.
- **Noun classes already consistent** (PascalCase): `Results`, `Accumulator`,
  `BlockAccumulator`, `CommonPrior`, `PrefetchLoader`, `BatchIndex`, `QRSystem`,
  `NormalEqSystem`, `ReducedContribution`.

### Mechanics (what was touched)

- **`.m` files renamed** in `+lsqvce/`; each file's `function` identifier updated
  to match its filename (MATLAB requires this), plus H1 help and self-refs.
- **HTML pages renamed** (`html/lsqvce.<old>.html` -> `<new>.html`), content +
  titles + h1 + every cross-link `href` retargeted; ALL hrefs verified against
  existing files (no broken links). NO redirect stubs (hard rename, per user).
- **Registration**: `Contents.m`, `helptoc.xml`, `html/index.html` entries
  renamed.
- **Call sites** updated everywhere: source (`recover`, `adjustVCE`, `makeLoader`,
  `iterateNoiseModel`, `diagnose`, ...), tests (`test_vceMath/v4/v7/v9`,
  `test_separable`, `test_reliability`, `test_speed`, `example_*`), examples
  (`graceLoader`, `iterateSeparable`, demos), docs (`README.md`, example
  READMEs), and the Python cross-validation scripts (`validation/v4,v7,v9`,
  which reference the MATLAB names in their docstrings).
- **Grep-verified**: zero occurrences of any of the six old names remain
  anywhere in the tree (all file types).

### Convention going into item C

The public surface is now: **verbs for entry points, PascalCase nouns for
classes, literature names for kernels.** Item C (human-readable Name/Value option
names, camelCase) can proceed against these FINAL names. The one convention
still unsettled and shared with C: engine options are PascalCase (`MaxIter`,
`GroupBy`, `ShowSeconds`) while makeLoader options are lowerCamel (`obsTypes`,
`sysID`) — C is where that gets reconciled.

Roadmap §7 item A naming work: class rename + this awkward-name sweep DONE;
namespace sub-packaging still deferred.

---

_Addendum, 2026-07-17 (rename SeparableProblem -> adjustVCE — roadmap item A). Model: Claude Opus 4.8._

**Roadmap item A (the awkward class name) done: `lsqvce.SeparableProblem` is now
`lsqvce.adjustVCE`.** The name was chosen with the user to PAIR with the
in-memory shortcut `lsqvce.adjust`: `adjust` = "just solve it"; `adjustVCE` =
"solve AND discover the stochastic model." Behaviour is byte-for-byte unchanged;
this is a pure rename + a back-compat alias. NOT re-run on the live R2026a box
(no MATLAB here) — run `results = run_all_lsqvce();` (expect `test_separable` +3:
the three alias tests).

What changed:
- **The class** `+lsqvce/SeparableProblem.m` -> `+lsqvce/adjustVCE.m`
  (`classdef adjustVCE < handle`, constructor `adjustVCE(nc,systems,days,...)`).
  All internal error IDs `lsqvce:SeparableProblem:*` -> `lsqvce:adjustVCE:*`
  (badAxes, notConverged, badSystems, badAxes, validateBadResults,
  notHomogenised, plotHookFailed), and every help self-reference updated.
- **Deprecated alias** `+lsqvce/SeparableProblem.m` is now a THIN SUBCLASS
  `classdef SeparableProblem < lsqvce.adjustVCE`. Its constructor forwards all
  args to `adjustVCE` and, the FIRST time it is constructed in a session, warns
  `lsqvce:SeparableProblem:renamed` (persistent flag -> once per session, then
  silent), as agreed with the user. Because it is a subclass, every method
  (vceGroupsBy/shape/scale/outer/axes/plots/run/empiricalStochastic/validate)
  and the returned Results are inherited unchanged. Kept for ONE version.
- **Error-ID note:** the alias inherits `adjustVCE`'s constructor, whose
  inputParser uses `FunctionName = 'lsqvce.adjustVCE'`, so PARSE errors from
  old-name code now name `adjustVCE` — intentional (points users to the new
  name). Only the deprecation warning carries the old id.
- **Docs/registration:** HTML page renamed `html/lsqvce.SeparableProblem.html` ->
  `html/lsqvce.adjustVCE.html` (title/h1/nav/syntax/error-ids updated, plus a
  "renamed from" note linking the stub). The OLD filename is kept as a REDIRECT
  STUB (`<meta http-equiv="refresh">` to the new page + an explanatory
  paragraph), so bookmarks/links don't 404. `helptoc.xml` and `html/index.html`
  main entry now target the new page; every cross-link in other pages
  (Results, makeLoader, recover, util.log, util.showExport) retargeted to
  `lsqvce.adjustVCE.html`. `Contents.m` entry renamed (with a "Was
  SeparableProblem" note); `README.md` row renamed likewise.
- **Examples:** `demo_separable.m`, `demo_separable_data.m`,
  `examples/gravity/iterateSeparable.m` construct/reference `adjustVCE` now.
- **Tests:** `tests/test_separable.m` switched to construct `adjustVCE`
  (class check `'lsqvce.adjustVCE'`, error-id asserts `lsqvce:adjustVCE:*`), and
  gained THREE alias tests (`test_alias_is_subclass_and_constructs`,
  `test_alias_warns_renamed`, `test_alias_matches_new_name`) that prove the old
  name still builds a working `adjustVCE`, warns the rename id, and gives the
  same solution. `test_adjustBlock.m` comment repointed.

**Deferred (part of item A, by design):** grouping the flat `+lsqvce/` namespace
into sub-packages (`+model`/`+accum`/`+cov`/`+diag`/`+util`) was the "bigger
move; may defer" half of item A. NOT done — it touches every qualified call site
and is a separate, larger pass. Only the one awkward class was renamed, as the
roadmap's "at minimum" instruction said.

**Watch on the box:** MATLAB requires the classdef name to match the filename, so
`adjustVCE.m` MUST contain `classdef adjustVCE` (it does) and the alias file
`SeparableProblem.m` MUST contain `classdef SeparableProblem` (it does). After
extracting, rebuild the help search DB (`builddocsearchdb('<root>')` or
`lsqvce_setup`) so the renamed page is indexed. If any user script still calls
`lsqvce.SeparableProblem(...)`, it works but prints the one-time rename notice;
switch it to `lsqvce.adjustVCE`. Item C (human-readable option names) is the
natural next step and shares the camelCase-convention decision with this rename.

Roadmap §7 item A: **CLOSED (class rename; namespace sub-packaging deferred).**

---

_Addendum, 2026-07-17 (entry-point parity audit — roadmap item B). Model: Claude Opus 4.8._

**Roadmap item B done: audited single-point-of-entry parity across every public
front door, produced the parity table, and CLOSED the two genuine gaps it found
(`varComp` TraceMethod; `recover` postSolve plot hooks).** The recent feature
work (trace estimator, parametrisation pipe, plot hooks) had been driven through
`SeparableProblem`/`makeLoader`; this pass verified the siblings and fixed what
was left behind. NOT re-run on the live R2026a box (no MATLAB here) — run
`results = run_all_lsqvce();` (expect the prior grand total +6: `test_lsqvce`
+3, `test_adjustBlock` +3).

### Parity table (entry point × feature)

Legend: ✓ present · ✗ was a gap (now fixed) · n/a not applicable to that door.

| Entry point         | Parametrisation pipe (common/local/procOpt, nc/nl derived) | Plot hooks (postSolve/postModel + ShowSeconds/ExportDir) | Hutchinson TraceMethod (Probes/TraceSeed) |
|---------------------|:---:|:---:|:---:|
| `SeparableProblem`  | ✓ (via loader) | ✓ (both seams) | ✓ (orchestrates vce) |
| `recover`           | ✓ (nc = L.nc)  | ✓ **NEW postSolve** (no postModel: single solve, no model) | n/a (no VCE in recover) |
| `makeLoader`        | ✓ (source of it) | ✓ (carries the handles) | n/a |
| `vce`               | n/a | n/a | ✓ |
| `redundancy`        | n/a | n/a | ✓ |
| `varComp`           | n/a | n/a | ✓ **NEW** (machinery check; exact stays default) |
| `traceEstimate`     | n/a | n/a | ✓ (is the estimator) |
| `adjust`            | n/a (in-memory flat shortcut) | n/a | n/a |

Two gaps were real (exactly as parked): `varComp` had no `TraceMethod` though it
forms the same trace quantities as `vce`; `recover` fired no plot hooks though it
is a legitimate single-solve front door. Everything else was already in sync.

### Fix 1 — `varComp` TraceMethod parity (`+lsqvce/varComp.m`)

`varComp` forms `tr(Qxx N_g)` (the redundancies) and `tr(Qxx N_i Qxx N_j)` (the
Helmert matrix) EXACTLY — precisely the traces `lsqvce.traceEstimate` (Ellmer
§6.5.6) approximates — yet had none of the `'TraceMethod'`/`'Probes'`/`'TraceSeed'`
options its sibling `vce` gained. Added them, mirroring `vce`/`redundancy`
verbatim in surface and semantics:
- **EXACT stays default and authoritative.** The exact per-group traces are
  untouched; `VC.s`/`VC.sigma`/`VC.redund` are bit-for-bit identical with the
  estimator on or off (asserted by a test).
- **`'hutchinson'`** adds a MACHINERY CHECK only — `tr(N⁻¹N)=u` on the final
  normal via the matrix-free estimator — returned in `VC.traceEstimate`
  (.value/.expected/.relErr/.se/.relSE/.probes/.note), same struct shape as
  `red.traceEstimate`. This is deliberately NOT the per-group substitution; that
  remains roadmap item 0 (the per-obs-r_i vs per-group-scalar structural mismatch
  documented in the trace-estimator addendum below is unchanged). Wrapped in
  try/catch → `lsqvce:varComp:hutchFailed`, exact components unaffected on error.
- Docs: `varComp` help (options + output), `html/lsqvce.varComp.html`
  (Name-Value rows + `.traceEstimate` output + see-also to vce/redundancy/
  traceEstimate). Tests: `tests/test_lsqvce.m` +3
  (`test_varComp_exact_has_no_traceEstimate`,
  `test_varComp_hutchinson_check_returns_u`,
  `test_varComp_hutchinson_does_not_perturb_components`) + a `vcSmall` helper.

### Fix 2 — `recover` postSolve plot-hook seam (`+lsqvce/recover.m`)

`recover` now fires each loader's `plotPostSolve` hook once per (system,obsType)
stream after the solve, mirroring `SeparableProblem.firePlots('postSolve')`:
keyed `satKey|obsType`, ONE ctx struct in, a figure out, `showExport` owns
timing/PNG, every call wrapped so a throwing/figure-less hook is logged
(`lsqvce:adjustBlock:plotHookFailed`) and skipped — plots NEVER break the solve.
- **DESIGN DECISION (locked): postSolve ONLY, no postModel.** `recover` does one
  GLS solve with the loaders' covariance and estimates NO stochastic model, so
  the postModel seam (which fires "once after the stochastic model is estimated")
  has nothing to fire on. postModel stays exclusive to `SeparableProblem`. This
  is documented in both the help and the HTML page as the reason, not an omission.
- Engine knobs `'ShowSeconds'`/`'ExportDir'` added to the constructor (same names,
  same meaning as `SeparableProblem`: 0 = export only, <0 = show and leave open).
- ctx fields: `stage`('postSolve'), `stream`, `satKey`, `obsType`, `sysID`,
  `pass`(=1), `R`(the Results), `x`, `exportDir`, `showSeconds`. This is a
  SUBSET of the SeparableProblem ctx (no `model`/`sigmaG`/`nAxes`/`R.history` —
  none exist for a single non-iterative solve); a hook written for recover is
  forward-compatible with the SeparableProblem ctx, not vice-versa. Documented.
- Raw 4-arg closures carry no hooks and are skipped (the loader-object test in
  `normalizeSystems` now also returns the loader objects so their maps are
  reachable).
- Docs: `recover` help (PLOT HOOKS block + two options), `html/lsqvce.adjustBlock.html`
  (rows + section + see-also to SeparableProblem/showExport). Tests:
  `tests/test_adjustBlock.m` +3 (`testRecoverFiresPostSolveHook` — both
  obsTypes fire once with well-formed ctx; `testRecoverNoHookIsNoOp`;
  `testRecoverHookFailureNeverBreaksSolve`) + `recorderHandle`/`recordCtx`/
  `recorderLog` helpers (a containers.Map handle recorder, returns [] so no
  display is needed).

### Loader-consistency follow-up (2026-07-17, same session)

The item-B note below said the hand-written loaders "should be updated for
consistency." On inspection the honest fix was NOT to copy struct fields onto
them — `synthLoader`/`graceLoader` are RAW 4-arg closures, not Loader structs;
you cannot put `.plotPostSolve` on a bare function handle, and the engines
already accept raw closures as first-class (every optional Loader field is read
behind an `isfield` guard). The real inconsistencies were a name collision and
an under-documented contract:

- **Name collision removed.** `examples/gravity/iterateSeparable.m` defined a
  LOCAL function `makeLoader(cfgSat,satId)` that merely binds `graceLoader` into
  a closure — shadowing (by name) the unrelated, higher-level
  `lsqvce.makeLoader`. No runtime bug (bare name resolves to the local in a
  script), but a genuine readability trap. Renamed to `bindGraceLoader` with a
  comment explaining the distinction. Only one call site (line 150) + the def;
  no test depended on it.
- **Loader contract documented.** `lsqvce.makeLoader` help gained a "THE LOADER
  CONTRACT" block: the ONLY hard requirement is the 4-arg `.closure`; every
  other `L.*` field is optional and engine-guarded, so (a) makeLoader can add
  fields without breaking old loaders and (b) a hand-written raw closure is a
  first-class loader. Also documented the previously-undocumented output fields
  `.epochCov`/`.plotPostSolve`/`.plotPostModel`.
- **Template pointer fixed.** `lsqvce.synthLoader` (the "loader TEMPLATE") now
  spells out the TWO supported ways to use a copy — raw closure vs wrapped by
  `lsqvce.makeLoader` — and warns against hand-rolling a `makeLoader` of your
  own (exactly the trap iterateSeparable had fallen into).
- **Test.** `tests/test_adjustBlock.m` +1
  (`testRecoverRawClosureIgnoresPlotOpts`): a raw-closure system with
  ShowSeconds/ExportDir set still solves cleanly and recovers truth — locking in
  that the plot seam is a silent no-op over a bare handle, never an error.

Touched: `examples/gravity/iterateSeparable.m`, `+lsqvce/makeLoader.m` (help
only), `+lsqvce/synthLoader.m` (help only), `tests/test_adjustBlock.m` (+1).
No numeric behaviour changed. Combined with the two fixes above, `test_adjustBlock`
goes +4 this session (was +3) and `test_lsqvce` +3.

### Residual parity notes (checked, no code needed)

- **`makeLoader` vs hand-written loaders** (`graceLoader`, `synthLoader`,
  `PrefetchLoader`): the engine guards every new field with `isfield`, so absence
  is safe (confirmed in `SeparableProblem.plotHooks`/`epochCovHooks` and now in
  `recover.firePostSolvePlots`). Hand loaders that never set
  `nc/common/local/procOpt/plotPostSolve/plotPostModel` simply don't fire hooks —
  correct, not broken. Updating them to carry the fields for consistency is
  cosmetic and left for the naming/ergonomics passes (items A/C).
- **`adjust`** is the in-memory flat-parameter shortcut; it has no loaders, no VCE
  loop, and no block-arrow normal, so none of the three features applies. Its
  parity row is n/a by construction, not a gap.

### Watch on the box

- `varComp`'s Hutchinson check reuses `N` and `R` (the Cholesky) as they stand
  after the iteration loop — i.e. the CONVERGED normal. If `maxIter` is hit
  without convergence they are still the last iteration's normal (fine for a
  machinery check). If a future refactor moves the `chol` out of the loop body,
  re-point the check at the final factor.
- The recover hook test returns `[]` from its hook, so it needs no display; the
  showExport figure-writing path is already covered (and self-skips headless) in
  `test_util`. If `testRecoverFiresPostSolveHook` ever flakes, check that
  `buildDayImpl` still supplies both `pos` and `kbr` in `cfg.cov`.

---

_Addendum, 2026-07-16 (per-obsType plot hooks). Model: Claude Opus 4.8._

**New: per-(system,obsType) plot hooks in the separable loop, with timed show +
PNG export.** Two seams, exactly as requested: postSolve (after each solve,
BEFORE residuals) and postModel (once, after the stochastic model is estimated).

How it is wired:
- **Handles live on the LOADERS.** `makeLoader` gained `'plotPostSolve'` and
  `'plotPostModel'`: each is a single `@(ctx)->figHandle` (broadcast to all the
  loader's obsTypes) OR a containers.Map/struct keyed by obsType. Normalised to
  a containers.Map obsType->handle and carried as `L.plotPostSolve/.plotPostModel`.
  Bad (non-handle) values error `lsqvce:makeLoader:badPlotHandle`.
- **Engine collects + fires.** `SeparableProblem.plotHooks(stage)` mirrors
  `epochCovHooks`, composing the finest key `satKey|obsType`. `firePlots` builds a
  single ctx struct, calls the handle (wrapped in try/catch -- a throwing or
  figure-less hook is logged `lsqvce:SeparableProblem:plotHookFailed` and skipped,
  so plots NEVER break the solve), then hands the figure to showExport.
- **ctx** carries: stage, stream, satKey, obsType, sysID, nAxes, pass, R, x,
  exportDir, showSeconds, model (the fitted AR/PSD+sigma, postModel), sigmaG (the
  stream's VCE weight). It is ONE struct on purpose, so it can be enriched later
  without breaking hook signatures.
- **Engine knobs** (NOT per loader): `'ShowSeconds'` and `'ExportDir'` on the
  constructor, or the fluent `P.plots('ShowSeconds',s,'ExportDir',d)`. ShowSeconds
  0 = export only; <0 = show and LEAVE OPEN. Files:
  `<ExportDir>/<stage>_<makeValidName(satKey|obsType)>_pass<k>.png`.
- **New helper `lsqvce.util.showExport(figH,name,showSeconds,exportDir)`** owns
  the show/pause/export/close so hooks stay boilerplate-free. Defensive: []/invalid
  handle/headless => silent no-op returning ''. exportgraphics (R2020a+) at 150dpi,
  else print -dpng. Sanitises '|' -> '_'.

Design choices locked with the user: handles on loaders (not the engine); keyed
per (system,obsType) via `satKey|obsType`; hooks take one ctx struct and return a
figure; engine owns timing + export.

Demo: `demo_separable('Plots',true,'ExportDir',dir)` attaches EXAMPLE hooks
(postSolve = residual histogram via lsqvce.collectResiduals; postModel = AR order
+ VCE weight bar) and sets ShowSeconds/ExportDir. Off by default (no figures in
the normal demo run).

Docs/registration: makeLoader help + `html/lsqvce.makeLoader.html` (two options);
SeparableProblem help (PLOT HOOKS block + plots() configurator); `showExport`
registered in Contents/helptoc/index + new `html/lsqvce.util.showExport.html`.
Tests: `test_adjustBlock` +3 (broadcast, per-obsType map, bad-handle error);
`test_util` +4 showExport (empty/invalid no-op, writes-png, sanitises-pipe --
the figure-writing ones assumeFail without a display); `test_separable` +1
(`test_plot_hooks_fire` injects recording hooks and asserts both stages fire,
incl. the independent posE). NOT run here (no MATLAB) -- run `run_all_lsqvce()`.
The showExport figure tests self-skip on headless; the hook-firing test returns []
from its hook so it needs no display.

---

_Addendum, 2026-07-16 (demo update + satellite E). Model: Claude Opus 4.8._

**Demo updated for the parametrisation pipe AND given a third, independent
satellite.** `examples/demo_separable_data.m` was rewritten:
- Uses the new `makeLoader` pipe: builds two 14-column parametrisation tables
  (`commonParamTable` -> nc rows, `localParamTable` -> nl rows) and a `procOpt`
  struct, passes them as `'common'`/`'local'`/`'procOpt'`, and lets nc/nl be
  DERIVED (no explicit `'nl'`). `demoBuildDay` is now the 5-arg form and reads
  the counts from the tables (proving the pipe). The 14-column column convention
  is documented in the file header (kind/label/system/sat/n/m/CS/arc/comp/unit/
  x0/sigma0/estFlag/groupTag) — MEANINGS are the demo's convention; makeLoader
  only checks shape.
- **New system: satellite E (`sysID 2`), INDEPENDENT, positions only (`posE`).**
  No inter-satellite link, its own local block `y_E`, same shared common field.
  So the demo now returns `L = {Lpair, Le}` (was a single-cell `{Lcd}`), and VCE
  by obsType gives FOUR components (posC, posD, rr, posE). This contrasts with
  the C/D pair, which is MERGED into one system precisely because `rr` couples
  them — E needs no such merge.

Config changes: `cfg.nl` split into `cfg.nlPair`(=2*nlSat) and `cfg.nlE`(=nlSat);
added `cfg.obsTypesE={'posE'}`, `cfg.allObsTypes` (all 4), `cfg.y_E`,
`trueAR.posE`, `trueSigma.posE`, and the parametrisation tables on cfg. Removed
unused `whichSat` helper.

Driver `examples/demo_separable.m` and `html/lsqvce.SeparableProblem.html` updated
to describe both systems and the pipe. Tests (`tests/test_separable.m`):
count-dependent tests switched from `numel(cfg.obsTypes)` (3) to
`numel(cfg.allObsTypes)` (4); added `test_independent_satellite_E` (asserts
`numel(L)==2`, `Le.sysID==2`, `Le.obsTypes=={'posE'}`, derived `Le.nl==4`/`nc==6`,
and that `posE` appears as a stream end to end). NOT run here (no MATLAB) — run
`run_all_lsqvce()`; test_separable goes 24 -> 25 and several existing counts
change from 3 to 4 streams/groups.

Watch on the box: the loop now solves a 2-system problem; convergence and the
`R.history` shapes are [nOuter x 4]. If anything in `test_separable` still expects
3, it was missed — grep `numel(cfg.obsTypes)` in that file (should be gone).

---

_Addendum, 2026-07-16 (parametrisation pipe). Model: Claude Opus 4.8._

**Interface change: buildDay can now receive the parametrisation.** Previously
buildDay was called `buildDay(day,obsType)` and the parameter DEFINITIONS were
hidden inside the closure; `nc`/`nl` were just column counts. `makeLoader` now
accepts three options and pipes them to buildDay:
- `'common'` R×14 cell (one ROW per common unknown) -> buildDay 3rd arg; `L.nc`
  DERIVED as `size(common,1)`. The common field is SHARED across systems, so
  every loader must pass the same `common`; pass `nc = L.nc` to recover.
- `'local'`  R×14 cell (one ROW per local unknown)  -> buildDay 4th arg; `L.nl`
  DERIVED as `size(local,1)` (errors `lsqvce:makeLoader:nlConflict` if an
  explicit `'nl'` disagrees).
- `'procOpt'` struct (free-form processing options) -> buildDay 5th arg;
  carried as `L.procOpt`.

New signature: `buildDay(day, obsType, common, local, procOpt)`.
**Back-compatible:** the wrapper `callBuildDay` calls buildDay with as many of
the 5 args as it DECLARES (`nargin`), so every existing `@(day,obsType)` builder
(all current demos/tests/loaders use exactly that) is unchanged. A varargin
builder (nargin<0) gets all five. 14-column shape is validated
(`lsqvce:makeLoader:paramCols`).

Touched: `+lsqvce/makeLoader.m` (options, `checkParamArray`, `callBuildDay`,
`dayOf` signature, derived `L.nc`/`L.nl`, carried `L.common/.local/.procOpt`,
help), `+lsqvce/recover.m` (nc-from-L.nc note), `html/lsqvce.makeLoader.html`.
Tests: 5 new points in `tests/test_adjustBlock.m` (back-compat 2-arg, 5-arg
receipt, nc/nl derivation, nl-conflict error, 14-col validation) + `paramRows`/
`captureBuild` helpers. Not run here (no MATLAB) — run `run_all_lsqvce()` on the
box; the makeLoader-area count in test_adjustBlock goes 6 -> 11.

NOTE for next session: the 14 column MEANINGS are the user's convention;
makeLoader only checks shape and pipes through. If a future task needs the
toolbox to INTERPRET the columns (e.g. auto-build A/B from the parametrisation),
that is a new layer on top of this pipe, not a change to it.

---

_Addendum, 2026-07-16 (trace estimator). Model: Claude Opus 4.8._

**Added: `lsqvce.traceEstimate` — the Ellmer §6.5.6 Monte-Carlo (Hutchinson)
trace estimator — plus opt-in plumbing in `redundancy`/`vce`.**

Context: the user asked whether the empirical trace estimator that §6.5.6 of the
Ellmer thesis specifies for the VCE redundancies is implemented. It was not
(the toolbox computes those traces EXACTLY via per-arc leverage, and had
DELIBERATELY rejected a naive Hutchinson — see the note in `leverage.m`). The
user supplied the thesis; §6.5.6 (eqs. 6.5.39–6.5.48, 6.5.60–6.5.61) gives the
exact recipe: Rademacher probes ±1 (6.5.40), Zbar = Z/√Z (6.5.43), estimate
tr(Zbar'XZbar), and — the key optimisation — form U⁻¹Zbar ONCE from the global
Cholesky and reuse it across every arc and the geophysical redundancy. Z = 100.

**Provenance note (read this):** on opening the working copy this session, a
`+lsqvce/traceEstimate.m` AND `tests/test_traceEstimate.m` were already present,
stamped 21:54 today, though NOT in the user's uploaded archive and NOT
registered anywhere. They are a correct, faithful §6.5.6 implementation with a
sound 15-point test suite (unbiasedness, the trace(Zbar'MZbar) identity,
matrix-free equivalence, the eq.-6.5.48 shared-probe reuse checked against exact
trace(N\Ng), 1/√Z error law, Rademacher<Gaussian variance, RNG hygiene). Almost
certainly a stale draft from an interrupted prior attempt. It was AUDITED
against the thesis equations, verified for block-balance, and ADOPTED rather than
rewritten. If you want it re-derived from scratch to be sure of provenance, say so.

**What is wired (scope is deliberately conservative — read carefully):**
- `lsqvce.traceEstimate(X | Xfun, ...)` — standalone estimator. Explicit-matrix
  and matrix-free (`@(V)->X*V` + `'Size'`) paths; `'Probes'` accepts a shared raw
  probe block (the eq.-6.5.48 reuse); `'NumProbes'` (100), `'Seed'` (private
  stream, no global side effect), `'Distribution'` rademacher|gaussian. Returns
  `t` and `info.{se,relse,perProbe,Z,P,distribution}`.
- `redundancy` / `vce` gained `'TraceMethod','exact'(default)|'hutchinson'`,
  `'Probes'`, `'TraceSeed'`. **EXACT STAYS DEFAULT.** With `'hutchinson'`, the
  exact per-arc r_g are STILL computed and remain authoritative; the estimator
  runs a MACHINERY CHECK on the reduced common normal (tr(N⁻¹N) = nc exactly)
  and returns it in `red.traceEstimate` (.value/.expected/.relErr/.se/.relSE/
  .probes/.note). It is NOT yet a per-group redundancy replacement.

**Why not the full per-group substitution (honest limitation):** the exact
`redundancy` produces per-OBSERVATION r_i (needed for snooping, the (w,r)
quadrants, and the sum(r_i)=m−u identity guard). Ellmer's estimator produces a
per-GROUP scalar only, and needs a group-level common-normal assembly the
streaming per-arc path does not currently expose. A faithful full substitution
is a separate, non-trivial code path (a new group-level reduce), and it CANNOT
be validated here (no MATLAB in this environment). Rather than ship an untested
deep rewrite dressed as done, it is left as roadmap. The estimator + its reuse
API are complete and tested; the deep wiring is the remaining work.

Docs/registration: `Contents.m`, `README.md`, `helptoc.xml`, `html/index.html`,
new `html/lsqvce.traceEstimate.html`, and `tests/run_all_lsqvce.m` (suite added).
A 16th test (`test_redundancy_seam_trNinvN_is_nc`) covers the seam's operator
identity deterministically. Expected grand total (from the prior 326): +16
traceEstimate points → ~342, PLUS the estimator suite must be confirmed on the
live box: run `results = run_all_lsqvce();`.

**OPEN TOPIC (roadmap, not yet started): full per-group Hutchinson redundancy.**
Status: OPEN / parked. This is the deep substitution the estimator was meant to
enable; it is deliberately NOT implemented (the shipped `'hutchinson'` path is
only the machinery check). Recorded here in full so the next session can pick it
up cold.

The core issue is a STRUCTURAL MISMATCH, not a coding inconvenience:
- The exact `redundancy` produces PER-OBSERVATION r_i = 1 - h_i, and from those
  derives (a) data snooping's (w_i, r_i) test, (b) the (w,r) quadrant / MDB /
  reliability numbers, and (c) the integrity guard sum(r_i) = m - u (the toolbox's
  strongest self-test of the block-arrow inverse algebra, checked every call).
- Ellmer's estimator produces only a PER-GROUP SCALAR tr(R̄ V_n). That is exactly
  where the O(u^3) saving comes from (it never forms per-obs leverage), but it
  therefore CANNOT return the r_i that (a)/(b)/(c) consume. So it is not a flag
  inside the existing per-arc loop; it is a SECOND, parallel code path with
  different outputs.
- Second wrinkle: Ellmer's speed needs U^-1 Zbar formed ONCE from the GLOBAL
  Cholesky and reused per group (eq. 6.5.48). But GroupBy can be arc/obsType/
  system/day/handle, so each group's contribution to the COMMON normal, N_g^cc,
  must be assembled separately — and the streaming accumulator currently sums all
  contributions into one reduced normal, not separable by group. The group-level
  common-normal assembly has to be added too (check `BlockAccumulator` API first:
  is N_g^cc recoverable without a core change?).

Suggested TIERED plan (do in order, stop when accuracy/benefit no longer justify
the next tier):
  T0. Keep exact as the DEFAULT (already the case). Up to a few thousand unknowns
      it is both faster AND exact and gives the diagnostics. Estimator never
      default.
  T1. A SEPARATE entry point (e.g. `redundancyHutchinson`, or a mode) returning
      r_g and variance components ONLY — explicitly WITHOUT per-obs r_i, snooping,
      or the sum(r_i)=m-u guard, documented as such. Build U^-1 Zbar once,
      assemble N_g^cc per group, estimate each tr(R̄ V) via `lsqvce.traceEstimate`.
  T2. Wire T1 into `vce` as an opt-in SCALE backend (vce only needs r_g and vtv).
      Guard: if Collect/snooping is requested, fall back to exact and warn.
  T3. Validation harness BEFORE trusting it at scale: run both paths on a mid-size
      problem, check estimator r_g vs exact within a sampling band (few·sqrt(Var/Z));
      add a deliberately ill-conditioned case that DEMONSTRATES the degradation the
      `leverage.m` note warns about, so the failure mode is visible, not latent.

Caveats to carry: (i) none of this is verifiable in the no-MATLAB dev env — T1–T2
are untested until run on the R2026a box; (ii) Z=100 works for Ellmer because he
probes the WELL-conditioned operator U^-T N_g U^-1 (redundancies in [0,1]) and
averages over many arcs/iterations — on a low-redundancy group or near-rank-one
common block, 100 probes may be too few, which is what T3 exists to expose.

---

_Addendum, 2026-07-16 (late). Model: Claude Opus 4.8._

**Added: empirical stochastic modelling as a named method on `lsqvce.SeparableProblem`.**
_(Renamed from an earlier `ellmer` on request: the API carries no personal name; the
literature citations are kept.)_
Two new PUBLIC methods (no new numerics; they order + validate the existing loop):

- `R = P.empiricalStochastic(Name,Value)` — empirical stochastic modelling
  (Ellmer 2018, §6.5; Kvas et al. 2019) as ONE call: step-0 HOMOGENISE preflight
  (checks every multi-axis/position stream carries a fixed per-epoch D via the
  loader's `epochCov`/`Cep`; warns `lsqvce:SeparableProblem:notHomogenised` if
  not) → delegates steps 1–6 (SOLVE↔SHAPE↔SCALE) to `run()` unchanged → step-7
  VALIDATE. Options: `'Homogenise'`(true) `'Validate'`(true) `'Alpha'`
  (=shapeOpt.WhiteAlpha) `'Report'`(=verbose). Returns the SAME `lsqvce.Results`
  as `run()`, plus `R.history.empiricalStochastic` (the validation struct) when Validate.
- `E = P.validate(R,Name,Value)` — the acceptance test, standalone on any
  Results from run/empiricalStochastic. Gathers the SCALE evidence (per-group `sigma` from
  the last outer pass), the SHAPE evidence (per-stream Ljung–Box p-value via the
  new `ljungBoxWhiteP`), the formal parameter std `sqrt(diag(R.sol.Cxx))`, and a
  single `.accepted` verdict (sigma0≈1 AND every sigma_g≈1 AND every stream
  white). Options: `'Alpha'` `'SigmaTol'`(=outerOpt.Tol) `'Sigma0Tol'`(0.05)
  `'Report'`. Does NOT re-solve.

New private helpers on the class: `checkHomogenised`, `reportEmpirical`. New
file-locals: `pluckResidual` (thin wrapper over `collectResiduals`, matches
applyShape's accessor), `ljungBoxWhiteP` (p-value sibling of `ljungBoxWhite`,
same statistic). The empirical/formal proxy is `sigma0` (Cxx already carries
sigma0²); a full per-parameter empirical/formal comparison needs an external
truth ensemble and was deliberately NOT invented here.

Docs updated: `+lsqvce/Contents.m` (SeparableProblem entry), `README.md` (row),
`html/lsqvce.SeparableProblem.html` (new "Empirical stochastic modelling" section with the
step→call table and the two `.note`/`.warn` callouts). `helptoc.xml`/`index.html`
already register the page — a method addition needs no new page.
`examples/demo_separable.m` now calls `P.empiricalStochastic('Report',...)` instead of
`P.run()` to showcase it (same Results downstream).

Tests: 8 new points appended to `tests/test_separable.m`
(`test_esm_returns_results`, `_matches_run`, `_attaches_validation`,
`_no_validate`, `test_validate_on_run_result`, `_converged_accepts`,
`_rejects_non_results`, `test_homogenise_warns_when_scalar_positions` + a
`makePosDay` local). Suite goes 16 → 24; expected grand total 318 → 326.
**Not yet re-run on the live R2026a box** — run `results = run_all_lsqvce();`
to confirm (expect 326/326, 0 failed). If `test_homogenise_warns_...` is finicky
about the warning capture, note it uses `lastwarn` under a scoped
`warning('off',id)` + `onCleanup`, matching the toolbox's log/warn routing.

This addendum addresses roadmap item #3 (convergence/validation diagnostics),
realised as the formal-vs-empirical acceptance test. Everything below is the prior session.

---

_Session close, 2026-07-16 (evening). Model: Claude Opus 4.8._

---

_A 2026-07-16 evening session document used to follow here. **REMOVED
2026-07-19**: it claimed 318/318 tests and described `SeparableProblem` as the
current engine, both long superseded (the suite is 497/498 and the class is
`adjustVCE`). Keeping a stale status block at the end of the handover is how a
future reader gets the wrong baseline. Its substance is in the dated entries
above; the full text is in git history._
