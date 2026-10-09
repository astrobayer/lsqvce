# Review — swe-reviewer (lane: TESTS AND THE GATE) — main@b25bb11 (v6.1.3) — 2026-10-09

Scope reviewed: `dev/tests` (27 `test_*.m` suites, `ci_gate.m`, `run_all_lsqvce.m`, `run_lsqvce.m`, `private/`,
`category_exemptions.csv`, fixtures), `dev/docs_src/audit_tree.py` + `audit_selftest.py`, `dev/docs_src/error_ids.csv`,
`.github/workflows/ci.yml`. I reviewed the whole tree rather than a diff, as the review brief asked.
Tier: B. No MATLAB. I ran the repo's own Python audit (exit 0, "CLEAN") and its self-test (26/26 ok). My own parsers,
Python mirrors and fault injections on scratch copies of the tree are in
`scratchpad/mirrors/swe_tests/` (count.py, cov.py, ids.py, tol.py, rng.py, inject*.py, options_merge.py,
oncleanup_loop.py). Nothing under `/home/claude/lsqvce` was modified.
Prediction reconciled: my static count gives 576 declared and 564 gate-executed (576 minus `test_speed` 12), which matches
HANDOVER §12.34 exactly. The per-suite and per-tag figures are under "Checked and sound".

## Findings (most severe first)

### F1 — Seven audit checks scan a `tests/` folder that no longer exists, so no test file has been checked since L8 — REPRODUCED
What is wrong: checks 2, 3, 7, 10, 10b, 11 and 15 iterate `for d in ('+lsqvce', 'tests', 'examples')`
(audit_tree.py:163, 182, 200, 344, 450, 480, 516, 659). The suites moved to `dev/tests` at L8 (WP-69), and `tests/`
does not exist (`ls tests` → "No such file or directory"; `glob('tests/**/*.m')` → 0 files). Check 9 reads
`+lsqvce/adjustBlock.m` / `+lsqvce/adjustVCE.m` (line 416), and so does the second half of check 3 (lines 219-220). Both
files were moved or deleted, so those checks examine nothing, and check 9 still prints "0 definition(s), one field list".
Why: check 10 was written for the 2026-07-19 lesson, where a duplicated local function in `test_util` stopped the file
from loading and 41 tests silently never ran. Its target class is exactly the test files. A check that cannot see its
target is not a check, and its "OK" line is read as evidence.
Observable symptom: I injected each fault into a scratch copy of `dev/tests/test_signalCore.m` and ran the check.
Every check stayed silent, while the same fault in a package file made it fire:
```
10 duplicate local function in a TEST file (load error)   healthy:0 after-fault:0 SILENT
10 duplicate local function in a PACKAGE file (control)   healthy:0 after-fault:1 FIRES
7  unbalanced end in a TEST file                          SILENT   | PACKAGE control: FIRES
2  removed name used in TEST code                         SILENT   | PACKAGE control: FIRES
15 wrong H1 line in a TEST file                           SILENT   | PACKAGE control: FIRES
3  nc as leading arg to adjustBlock in TEST code          SILENT
FULL AUDIT (main()): duplicate fn + unbalanced end in a test file -> 0 findings, "CLEAN"
```
Evidence: `scratchpad/mirrors/swe_tests/inject.py`, `inject3.py`, `inject4.py`. The MATLAB gate's rule 1 (load error)
still catches the load error at run time, so the damage is to the fail-fast tier, not to the verdict.
Severity: S3.

### F2 — Check 16 counts comment mentions and the never-gated `test_speed` as coverage; counted honestly, 18 gaps have no exemption row, including 3 `contract` gaps — REPRODUCED
What is wrong: check 16 assigns a category to a unit from any `lsqvce.<unit>` string in the test body, comments
included (audit_tree.py:739-741, `body = '\n'.join(lines[i+1:end])`). It also scans `test_speed.m`, which `ci_gate`
excludes and which never runs in the gate (`excluded = {'test_speed'}`, ci_gate.m:85).
Why: §6.1.2 says the matrix is "DERIVED per run from the `lsqvce.<name>` references in each test body". A comment
is a reference to a name, not a call. A benchmark that never runs certifies nothing.
Observable symptom: my mirror of check 16 found 17 (unit, category) cells credited only by comments. Examples:
`engine.diagnose` contract (test_metamorphic.m:35, test_vceFunction.m:76 are comments), `engine.resolveNc` contract
(test_adjustBlock.m:236 is a comment), `engine.adjustNonlinear` contract+reference (test_fit.m:869 is a comment),
`engine.ftest` contract+reference, `kernel.varWhiten` reference. It found 9 `speed` cells reached only through
`test_speed`: BlockAccumulator, QRSystem, makeLoader, synthLoader, armaACV, armaWhiten, decorrFactor, levinsonDurbin,
psd2acv. Counting only code, and only gated suites, leaves **18 gaps with no exemption row**. Three of them are
`contract` for units that validate their arguments: `engine.adjustNonlinear` (`arguments` block; `badX0`,
`nonFinite`), `engine.diagnose` (`arguments`; 5 ids) and `engine.resolveNc` (3 ids). All five `diagnose` ids are fired
by no test (see F4).
Evidence: `python3 -I mirrors/swe_tests/cov.py` prints the "credited ONLY by comment mentions" list and
"HONEST MATRIX … total un-exempted gaps: 18". The audit itself prints "category gaps … every one either covered or
argued" and exits 0.
Severity: S3.

### F3 — Six `contract` exemptions break the binding "never for a validating function" rule, their reasons are "inconvenient", and nothing refuses stale rows (24 exist) — REPRODUCED
What is wrong: §6.1.2 says "`contract` may never be exempted for a function that validates its arguments", and the CSV
header repeats it. Six of the nine contract rows exempt units that do validate:
- `engine.runMeta`: narginchk, `runMeta:badArg` (runMeta.m:43-49)
- `engine.ReducedContribution`: `reduced:notSquare`, `:sizeMismatch` (ReducedContribution.m:43-48)
- `engine.residuals`: `arguments`, `badSol`, `noStore`
- `engine.ftest`: `arguments`, 5 ids
- `engine.adjustRobust`: `arguments` block, adjustRobust.m:87-99
- `engine.declaration`: 16 error sites

The reasons are "covered through lsqvce.fit", "its contract points lived in a deleted suite" and "never called by
name" (category_exemptions.csv:17-21, 306). Those are statements of inconvenience, not impossibility. Check 16
enforces neither the "validates" clause nor the doctrine's second mechanical rule ("something must refuse an entry
whose category has since become covered").
Observable symptom: 24 exemption rows exempt a category that check 16's own derivation shows as covered. Examples:
`engine.residuals,contract`, `engine.ftest,contract`, `engine.localColumns,contract`, `kernel.robustWeights,robustness`,
`layout,metamorphic`. Using code-only coverage, 19 rows are still stale. Overall: 314 rows, but only 54 distinct
reasons, and the most common reason is pasted 42 times.
Evidence: `mirrors/swe_tests/cov.py` (the "stale" lists), the grep of `error(`/`arguments` per unit, and audit_tree.py:747-794.
Severity: S3.

### F4 — 156 of 301 catalogued identifiers are fired by no test; 3 are USER-tier; WP-90's E-10 promise ("the surviving ids are pinned") is not met — REPRODUCED
What is wrong: §6.1.2 `contract` says "argument rejection with the DOCUMENTED identifier", and the skill rule is "every
documented identifier fires". 156 of the 301 ids in `error_ids.csv` appear on no code line in `dev/tests`.
- **USER tier (3):** `lsqvce:Results:noAccumulator`, `lsqvce:Results:noVPv`, `lsqvce:fit:lawOnFixedStructure`.
- **ENGINE tier (125):** e.g. `adjustVCEjoint:badSigma2_0`, `:collapsed`, `:notPD`, `:noArcs`, `:needCfcn`; all
  `normalEq:*` and `qr:*` refusals except `singularR`; `covModel:noSigma`, `:badCep`, `:noPSD`; `diagnose:*`;
  `snoop:*`; `vce:badSigma0`, `:groupMismatch`; `makeLoader:ncConflict`, `:priorType`, `:rowWeightSize`.
- **KERNEL tier (21):** `separable:*` (5 errors), `psd2acv:tooFewPoints`, `redundancyLags:badEbar`.
- **Private (7).**
- By kind: 113 errors and 36 warnings.

The full list is printed by `mirrors/swe_tests/ids.py`.
Contradiction of a recorded disposition: HANDOVER §12.14 E-10 lists `covModel:badCep` and `adjustVCEjoint:badAxisRows`
among the ids "no suite fires" and records the disposition "the surviving ids are pinned (WP-90)". Both survive in the
code (covModel.m:195/202/371/380/458; adjustVCEjoint.m:891), and grep finds both in 0 test lines.
Correction to coordinator C2: `lsqvce:adjustVCEjoint:arcMismatch` IS fired, in
test_ellmer.m:1378-1384 (`test_joint_arc_mismatch_raises`, `verifyError(... 'lsqvce:adjustVCEjoint:arcMismatch')`).
`badSigma2_0` is not fired, so C2 holds for that id only.
Also: no audit check enforces "documented ⇒ fired by a test". Check 14 only checks "cited ⇒ raised", and it counts
ids raised by test files as raised (see F11).
Severity: S3.

### F5 — `testDistributeDaysThroughFitEqualsFused` cannot fail: the second `Parallel=` argument resets `Mode` to "serial" — REPRODUCED (Python mirror of Options.m)
What is wrong: test_fit.m:744-745 calls
`lsqvce.Options(Verbose=false, Parallel=struct('Mode',"parfor"), Parallel=struct('DistributeDays',true,'Threshold',2))`.
The constructor assigns each pair in order (Options.m:182-186). `set.Parallel` merges onto the **defaults**, not onto
the current value (Options.m:195, `mergeGroup(struct('Mode',"serial",…), s)`). The second assignment therefore leaves
`Mode = "serial"`, and BlockAccumulator.m:825 (`distrib = opt.DistributeDays && ~strcmpi(opt.Parallel,'serial') …`) is
false. The distributed code path is never entered.
Why: the test then compares the fused path with itself at RelTol 1e-11, and it asserts DistributeDays and Threshold in
the provenance but never `Mode`, so it passes whatever the distributed fold does.
Observable symptom: the mirror prints `resolved Parallel: {'Mode': 'serial', 'DistributeDays': True, 'Threshold': 2}`
and `BlockAccumulator distrib (3 days): False`. Adding `tc.verifyEqual(R2.meta.options.Parallel.Mode, "parfor")` would
turn the test red.
Evidence: `mirrors/swe_tests/options_merge.py`. Assumed MATLAB semantics: Name=Value pairs reach `varargin` in order,
duplicates included. The engine-level twin in test_parallel still pins the distributed fold. What is unproven is the
fit-level claim.
Severity: S3.

### F6 — `warning('off','all')` in 15 test points makes the warning gate blind there, and the loop form in test_reliability alternates the state — INFERRED (state sequence mirrored)
What is wrong: §6.1 / WP-47 make `muteWarn` "THE one way a test silences a warning" (private/muteWarn.m:5) and
forbid silencing ids that cannot fire. Fifteen points still turn **all** warnings off:
- test_contribution.m:170, 224, 269, 302, 360, 383
- test_ellmer.m:377, 402, 1230, 1256, 1286, 1317, 1361, 1378
- test_adjustBlock.m:198
- plus test_reliability.m:374 inside a loop

Why: run_all parses warnings from the console capture (run_all_lsqvce.m:180-186), and gate rule 7 reds on any foreign
id. A new warning raised inside these 15 points is never displayed, so it never reaches the inventory. That is
16 points of the gate's warning rule disabled. At test_adjustBlock.m:198-199 the point named "must warn about nothing"
sets `lastwarn('','')` and never reads it, so the claim in its comment is not asserted at all.
test_reliability.m:371-375 builds `c = onCleanup(@() warning(w))` inside a 5-iteration loop. Reassigning `c` runs the
previous restore, and `w` captures the already-all-off state on every even iteration. The mirror shows the state during
the iterations as ALLOFF, ON, ALLOFF, ON, ALLOFF, and ON at exit. With a 6th branch (an even count), the test would leave
`ALLOFF` for every later suite: the inventory empties, and any later `verifyWarning` breaks. That is the §6.1.2b
"cluster of unrelated failures" mechanism.
Evidence: grep output listed above; `mirrors/swe_tests/oncleanup_loop.py`. The mirror assumes documented MATLAB
semantics: `warning('off','all')` returns the prior state, and an onCleanup is destroyed on reassignment.
Severity: S3.

### F7 — Copied `tol:` comments state the wrong source for loose bounds, and the provenance rule cannot tell — REPRODUCED
What is wrong: §6.1.2 asks for "a comment saying where the tolerance came from". Several bounds carry a comment
pasted from a different kind of comparison:
- test_lsqAccum.m:104: `AbsTol 5e-3` on recovering the truth from noisy data, labelled "two backends … difference
  ~ cond*eps".
- test_lsqAccum.m:223: `RelTol 1e-3` on an AR(20) approximation against a dense Toeplitz, same label.
- test_parallel.m:330: `RelTol 2e-2` on a batched approximation without overlap, labelled "independent dense solve …
  difference ~ cond*eps".
- test_traceEstimate.m:269: the bound on a Hutchinson trace of the identity, labelled "quantile/CDF round-trip
  through base-MATLAB special functions".

The bounds may be right. The stated sources are false, so a future failure there gets diagnosed against the wrong
cause.
Why the instrument misses it: check 16(d) only tests that the substring `tol:` appears within 3 lines
(audit_tree.py:813). In a fault injection, a bare `% tol:` above `'AbsTol', 0.5` passed. `% tol: as above` appears
23 times. The rule also only covers literal `'AbsTol'`/`'RelTol'`. 181 gated `verifyLessThan`/`verifyGreaterThan`
bounds carry numeric literals, and 146 of them have no `tol:` comment. Examples: test_physOracle.m:213 RelTol 0.25,
test_adjustBlock.m:55-56 `< 0.15`.
Evidence: `mirrors/swe_tests/tol.py`; fault case "16 tol: comment saying nothing" in inject.py → SILENT.
Severity: S3.

### F8 — Speed budgets: absolute seconds, a sequential `timeit` pair, and a budget that times only test-local copies of the code — REPRODUCED (grep)
What is wrong:
- **Absolute budgets.** §6.1.2 says absolute seconds are "LOGGED, not asserted", and R5 H-14 retired absolute budgets.
  Two remain in the gate: test_ellmer.m:1160-1182 (`verifyLessThan(tElapsed, 5.0)`, tic/toc) and
  test_speedBudget.m:69-78 (`verifyLessThan(dt, 20)`). The ledger accepted the second (M7 7a), which contradicts §6.1.2.
- **Sequential timing.** test_vceMath.m:534-535 times dense and then block with two separate `timeit` calls. It does
  not use the "ONE timing instrument" `timePair` (private/timePair.m:6-10): no interleaving, no rotation, no
  per-repeat ratio.
- **Budget on a copy.** test_speedBudget.m:36-49 (`testFftLagSumsBeatTheLoop`) times `lagSumsFft`/`lagSumsLoop`,
  which are local functions of the test (lines 80-94). There is no `lsqvce.` call on any code line of the point.
  `kernel.redundancyLags>columnLagSums` can regress and this budget stays green.
- **Repeat counts.** The two real ratio budgets use Repeats 7 and 5 (lines 47, 64). timePair's own help measures a
  31 % spread at 5 repeats against 10 % at 15.

Severity: S3 for the budget on a copy; S4 for the rest.

### F9 — `test_parallel` opens a parallel pool and never closes it; the suites after it run with the pool resident — INFERRED
What is wrong: `mapreduce` defaults to `AutoPool = true` and calls `parpool` when none is open (mapreduce.m:92, 202-208).
The points at test_parallel.m:75, 91, 153, 230, 253, 389, 414, 434 and 568 request parfor/spmd with Threshold 1.
No `delete(gcp)` exists anywhere in `dev/tests` or the runners (grep returns nothing).
Why: validation.md Part 10 says tests must not change the machine they measure. test_parallel is suite 11 of 26.
test_speedBudget (25), test_vcePolicy, test_fit and test_examples then run with workers resident, and the ratio
budgets are measured on a machine that is busier than the one they were pinned on. Any later point that requests
parfor silently uses the warm pool. The CI runner is a persistent, self-hosted PC.
Observable symptom: after `ci_gate`, `gcp('nocreate')` is non-empty. Timing budgets also differ between a full run and
`run_lsqvce('test_speedBudget')` alone.
Severity: S4.

### F10 — Global state that is not restored: RNG, path when suites run standalone, and a mute on an id that does not exist — REPRODUCED (grep)
- **RNG.** `rng(42)` (test_lsqvce.m:18) and `rng(7)` (test_lsqvce_ext.m:22) are set in `setupOnce`, and test_ellmer has
  22 `rng(n)` calls. None restores the RNG. Twelve gated points draw from the global stream without their own seed:
  - test_ellmer.m:228, 357, 438, 462, 527
  - test_lsqAccum.m:73 (a precision identity)
  - test_vcePolicy.m:162, 511
  - test_adjustBlock.m:188, test_fit.m:584, test_traceEstimate.m:210
  - test_streamUtils.m:166, `makeTinyDay`, which is dead: it is never called

  §6.1.3 says "Randomness is seeded per test". Here the inputs depend on suite order, so a failure in a subset run may
  not reproduce in the full run, and the reverse.
- **Path.** Seven suites `addpath` in `setupOnce` with no teardown: test_errorDoctrine:29-30, test_metamorphic:23,
  test_speedBudget:32, test_streamUtils:19, test_vcePolicy:25, test_examples:46 (adds `template/`), and
  test_physOracle:34/41. `run_all_lsqvce` restores the path, but `run_lsqvce` does not (it has no `pathBefore`), and
  neither does a bare `runtests`, which every suite header documents as "Run:". L8 says "tests leave the session's
  path exactly as found".
- **Dead mute.** test_signalCore.m:151-152 turns off `lsqvce:ar2arma:nonStationary`, an id no code raises (it is not
  in `error_ids.csv`), and restores it with `warning('on', …)` rather than to the prior state. That is a
  cannot-fire mute, which muteWarn's own help forbids.

Severity: S4.

### F11 — Fault injection covers 10 of 30 checks; check 16's self-test covers only its vocabulary; check 14 accepts ids raised only by tests — REPRODUCED
What is wrong: audit_selftest.py has cases for checks 16, 19, 20, 21, 24, 25, 26, 27, 28 and 29 only. Checks 1-15, 10b,
17, 18, 22 and 23 have no fault fixture. That is 20 of the 30 checks, and they include the blind ones in F1. The two
check-16 cases test category vocabulary only (an unknown category fires, `metamorphic` is legal). Nothing exercises
the matrix derivation (F2), the contract-gap finding, the exemption logic (F3) or the `tol:` rule (F7).
Check 14 builds its "thrown" set from `dev/tests/*.m` as well as the package (audit_tree.py:620-621). In a fault
injection, a help block citing `lsqvce:zz:ghost` plus a test-only `error('lsqvce:zz:ghost')` left both check 14 and the
full audit silent.
Two more points:
- **Two definitions of "declared".** Checks 16 and 17 count a test signature inside a `%{ … %}` block as declared.
  `declaredPoints.m` does not, although its help says check 16 "counts the same definition". With a block-commented
  test, the audit and the runtime counter disagree by one.
- **C1, extended.** Check 28 reads only help blocks before "See also" plus md/template/examples. It never reads error
  message bodies, `dev/tests`, or `html/*.html`.

Evidence: `inject.py` / `inject4.py` cases 14 and 17; the self-test list printed by grep.
Severity: S3.

### F12 — The "stated once" test guideline (HANDOVER §6.1) contradicts the gate it governs — REPRODUCED (grep)
What is wrong: §6.1.1 names `tests/run_lsqvce2.m` as "THE ONE runner" (line 2203). The runner is
`dev/tests/run_all_lsqvce.m`, wrapped by `ci_gate`. §6.1.2 still defines **six** categories and says METAMORPHIC "is no
longer a category of its own". Its "MECHANICS (binding, enforced by audit check 16)" paragraph lists the six names and a
TOL_SWEPT start list. Check 16 has enforced **seven** categories since R5 M0 (audit_tree.py:701-703, WP-78), and
TOL_SWEPT now covers every suite (line 710).
Why: §0 makes §6.1 the guideline every prompt references. A session that follows it writes six-category declarations
and points readers to a runner that does not exist. The speed row, "absolute seconds LOGGED, not asserted", is also
contradicted by test_speedBudget (F8).
Severity: S3.

### F13 — The skill's `count_tests.py` reports 0 classes on this tree and exits 0 — REPRODUCED (note on the instrument, not the repo)
Why it fails: `count_file` returns `None` unless a line matches `^\s*classdef` (count_tests.py `CLASSDEF`). All 27
lsqvce suites are **function-based** (`tests = functiontests(localfunctions)`). There is no classdef, no
`methods (Test)` block and no `TestParameter` anywhere (grep), so the script sees nothing. It then prints "0 class(es),
0 declared point(s)" with exit code 0: it fails open. As validation.md Part 3 says, "a guard that prints nothing when
nothing is wrong is evidence only when an independent count agrees". The script should refuse when it finds zero
classes, or parse function-based files. My own parse (`mirrors/swe_tests/count.py`) uses `declaredPoints.m`'s
definition and mirrors MATLAB's rule for local functions (`test*`/`*test`, minus setup/teardown).
Severity: S4.

## Checked and sound
- **Prediction.** The 27 files declare 576 points in total. The gate executes 564 (everything except `test_speed`, 12).
  That equals HANDOVER §12.34 (576/564) and DEVELOPMENT.md's "27 suite files, 576 test points". Points declared per
  suite:
  - 30: adjustBlock, 21: ci_gate, 16: contribution, 71: ellmer, 17: errorDoctrine, 2: examples, 39: fit
  - 7: instruments, 27: lsqAccum, 16: lsqvce, 26: lsqvce_ext, 27: lsqvce_opts, 12: metamorphic, 7: options
  - 25: parallel, 11: physOracle, 31: precision, 37: reliability, 6: signalCore, 12: speed, 3: speedBudget
  - 15: streamUtils, 18: traceEstimate, 43: util, 7: vceFunction, 19: vceMath, 31: vcePolicy

  For every suite, MATLAB's local-function test rule and `declaredPoints.m` agree. There are no `TestParameter`
  multipliers.
- **Categories.** Every one of the 576 points has a legal `% categories:` line. Category totals: contract 304,
  precision 184, robustness 162, reference 54, metamorphic 21, speed 18 (6 of them gated), vectorisation 14. Tags:
  diagnostics 68, provenance 13. No point lacks a verify/assert/assume call.
- **Gate logic** (gateVerdict.m, read): it is red on a load error (NaN total), on failures, on incompletes, on
  disk⇄run set differences in both directions, on a missing or empty log, on declared ≠ executed, on a foreign warning
  id, and on a blind capture (the probe is missing). `ci_gate` raises `ci_gate:red` (ci_gate.m:91-94), so `-batch` and
  CI exit non-zero. test_ci_gate (21 points) feeds synthetic runs into every rule.
- **Skips cannot pass silently.** A point filtered by `assume*` becomes Incomplete, which is red under rule 3; this
  covers test_parallel's 8 PCT points, test_physOracle:408 and test_vcePolicy:172. A requested skip is printed and
  logged.
- **State the runner restores.** `run_all_lsqvce` restores the path and the 'verbose' warning state through onCleanup.
  test_examples restores `DefaultFigureVisible` and closes figures. `muteWarn` restores through `addTeardown`. The log
  tests restore the log state through cleanups.
- **Statistical tests** use seeded private `RandStream`s. The unseeded draws in F10 sit in contract and robustness
  points, plus one precision identity that holds for any input.
- **Audit and CI wiring.** `audit_tree.py` exits 0 and `audit_selftest.py` passes 26/26, reproduced here with pandoc
  present. ci.yml runs audit + selftest + v10review, then the `matlab` job (`needs: audit`) runs ci_gate on the
  self-hosted runner and uploads `last_run_all.txt` with `if: always()`.
- **Timing helper.** `timePair` follows the doctrine as written: warm-up, interleaving, rotation, the median of
  per-repeat ratios, and a band in the record.

## Not reviewed
- MATLAB runtime behaviour: whether `runtests` throws or returns on a file that fails to parse, worker-side warnings,
  and how warnOnce's persistent state interacts with suite order. These are not reachable at tier B; F6 and F9 are
  inferred.
- The numerical correctness of the oracles behind reference/precision values (domain lane), examples, html pages, and
  `release.yml` beyond its header.
- Whether the 11 "raised but undocumented in help" ids (check 20's reported M1 debt) are right per unit (documentation
  lane).
