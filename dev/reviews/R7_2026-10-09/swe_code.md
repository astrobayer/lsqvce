# Review — swe (lane swe_code) — lsqvce main@b25bb11 (v6.1.3) — 2026-10-09

Scope reviewed: the code under `+lsqvce` (user tier fit / Options / Results / report / layout,
`+engine`, `+kernel`, `+util`, `private/`) and its interfaces, chosen because the brief asks for a
whole-package code review rather than a diff. Estimation theory is out of scope (other lane); where
a finding touches it (F3) the code fact is stated and the remedy is left to that lane.
Tier: B. No MATLAB. Instruments run: `dev/docs_src/audit_tree.py` (CLEAN, 29 checks, rc 0; output
kept at scratchpad/audit_swe.txt), the audit's own `raised_ids()` re-used in a mirror, grep, and
Python mirrors of the arithmetic in `scratchpad/mirrors/swe_code/` (rowkey.py, rmsweights.py,
vardense.py, madfallback.py, madfallback2.py, inventory.py, groupfields.py, nvdefaults.py,
ids.py, funlen.py). REPRODUCED means a mirror or grep output shows it; INFERRED means read only.
Already reported by the coordinator and not repeated: C1, C2, C3.

## Findings (most severe first)

### F1 — `util.rowWeightKey` formats the day with `%g`: day ids of 7 or more digits collide, and robust factors from one day are applied to another — REPRODUCED
What is wrong: `+lsqvce/+util/rowWeightKey.m:26` `s = sprintf('%g', double(v));` is used for the
sysID, the day and the arcId. `%g` keeps 6 significant digits, so every day id ≥ 1e6 that differs
only after the 6th digit gets the same key. The same day is keyed with `%d` in the loader cache
(`makeLoader.m:1062` `key = sprintf('%d|%s', day, obsType);`), so the cache is right and the weight
map is wrong. That is two definitions of one fact (R1), and they disagree.
Why: the map is written per arc in sweep order (`adjustRobust.m:168`, `adjustVCEjoint.m:659`) and
read per arc (`makeLoader.m:1261-1263` epochWeights, `:1314-1315` rowWeights). It is also used to
undo the factors (`adjustVCEjoint.m:732`, `private/lawWhitenedResiduals.m:59`). When keys collide,
the last day written wins, and every colliding day reads that day's factors back.
Observable symptom: `data.days` in YYYYMMDD form (20260901..20260930) or as Julian Day Numbers
(2460580...) with `Robust ~= "none"`. If arcs have equal length (fixed-length days, the usual case),
day k's blunder weights are applied silently to day 30, and the factors are divided back out by the
wrong vector. If the lengths differ, `lsqvce:makeLoader:rowWeightSize` is raised with a message
that blames the builder ("must return the same arc on every call"). The contract allows these ids:
Manual 03_data.md:192 says "day identifiers ... integers", and `declaration.m:75` requires only
finite integers. MJD (59396) and datenum (739000) are below 1e6, so they are safe. A numeric sysID
≥ 1e6 collides in the same way.
Evidence: rowkey.py mirror:
`YYYYMMDD Sep 2026  days=30 rowWeightKeys unique=1  dayKey unique=30  e.g. 1|2.02609e+07|kbr|1`;
`day-of-campaign 1e6+k days=30 rowWeightKeys unique=4`; JDN 2460580..2460599 gives 3 unique keys.
The only pin is `test_errorDoctrine.m:267` (day 59396), which is below the threshold.
Severity: S1 (silently wrong numbers on equal-length days). A struct key, or `%d` / `%.17g` taken
from one helper, closes it.

### F2 — a `var` law under `Estimate="noise"` uses the wrong covariance in the joint sweep, and its level is applied as a std to a covariance — REPRODUCED (arithmetic) / INFERRED (reachability)
What is wrong: (a) A `var` stream's loader sets `out.nAxis = d` (`makeLoader.m:1361`).
`adjustVCEjoint>oneArc` then takes the multi-axis branch (`:876`) and calls
`Sig = a.Cfcn(1:mE, uo); % TEMPORAL, over epochs` (`:899`). For `var`, though, `Cfcn` is
`varDense(Av, Sg, d, numel(idx))` (`covModel.m:318`), whose argument is a ROW count (`:588-596`).
It returns the d-interleaved block-Toeplitz matrix over mE/d epochs, not the per-axis temporal
covariance over mE epochs. (b) No policy refuses or special-cases `var`. Its d×d `sigma` makes
`fixedStructure` true (`fit.m:585`), so the policy becomes "scale". `applyCovOverride.m:51`
`cm.sigma = getfielddef(cm,'sigma',1) * cm.sigmaScale` then multiplies the innovation COVARIANCE by
a STANDARD-DEVIATION ratio (`jointModelOpt.m:174-181`, `lvl = sqrt(...)`), so the level applied is
the square root of the one estimated.
Why: `oneArc`'s multi-axis branch assumes a temporal law over epochs composed with Cep, which is the
WP-90 design. `var` is the one law that "carries its own cross-axis structure" (covModel help
:26), and nothing excludes it from the joint loop.
Observable symptom: `fit(data, Estimate="noise")` with any `covariance.<t>.type = 'var'`. With an
even number of epochs per arc, the Ω/s statistics are formed against a wrong Σ (VAR(1), a=0.5, d=2:
the lag-1 entry used is 0.000 against a true 0.667) and the run returns numbers without a warning.
With an odd count it dies with `lsqvce:covModel:notWholeEpochs`, an internal id that does not say
var + noise is unsupported. No test runs `var` under `Estimate="noise"` (grep: `'var'` appears only
in test_parallel / test_streamUtils, all Estimate "none").
Evidence: vardense.py: `nE=6: Sig used by oneArc (6x6) == per-axis temporal covariance? False;
max|diff| = 0.667`; `nE=7: oneArc -> lsqvce:covModel:notWholeEpochs (mrows=7, d=2)`.
Severity: S1 where reached (wrong numbers, silent). The minimal honest fix is a named refusal in
`fit>streamPolicy`.

### F3 — the joint loop's multi-axis robust statistic is strictly positive, and the kernel standardises `vw/s` without centring: on CLEAN 3-axis data IGG3 zeroes 24 % of the epochs — REPRODUCED
What is wrong: `adjustVCEjoint>epochResiduals` reduces a d-axis epoch to its RMS,
`u = sqrt(mean(reshape(v, d, []).^2, 1))` (`:741`), and `perStreamWeights` passes that straight to
`kernel.robustWeights` (`:704`). The kernel computes `s = 1.4826*MAD` about the median, but then
uses `w = vw / s` (`robustWeights.m:68`), not `(vw - med)/s`. A chi/√d variable has a median near 1
and a small MAD, so `w` sits at about 2.3 for a typical clean epoch, past IGG3's k0 = 1.5.
This is the "bodily offset" route that WP-110 (§12.34) measured as the way the guard fires. Here it
is built into the statistic.
Why: §12.29 round 4 recorded the decision "the epoch statistic is the RMS, standardised by the
stream's MAD like a scalar residual -- a chi-shaped variable treated as symmetric; adequate for
blunder rejection". The measurement below CONTRADICTS "adequate". It is reported as a contradiction
of a recorded decision, not as a re-proposal. The scalar (d = 1) path is unaffected, and that is the
only path the oracle O-105 exercises (two scalar AR(1) streams).
Observable symptom: `Robust="igg3"|"bisquare"|"huber"` with `Estimate="noise"` on a multi-axis
stream (the GRAVFIRE 3-axis positions, demo_grace). Clean epochs are down-weighted or zeroed every
pass, and the levels and law come from a fraction of the data. The collapse guard does not fire,
because nAlive stays large.
Evidence: rmsweights.py (200 000 clean N(0,1) epochs, unit level):
`d=1 igg3 zeroed=0.3%`, `d=3 igg3 zeroed=24.1% mean factor=0.364`,
`d=3 huber factor<0.5=33.7%`, `d=6 igg3 zeroed=62.6%`, `d=12 igg3 zeroed=97.4%`.
Severity: S1 (silently wrong weights on reachable input). The remedy belongs to the estimation lane:
a centred or chi-calibrated scale.

### F4 — `kernel.robustWeights` falls back MAD → std → 1 in silence, and the fallback changes the answer by two orders of magnitude (R4) — REPRODUCED
What is wrong: `robustWeights.m:62-67`:
`if s <= 0 || ~isfinite(s), s = std(vw); if s <= 0 || ~isfinite(s), s = 1; end, end`. There is no
warning, no flag in the outputs, and nothing in `R.history.scale` that says which scale was used.
The help documents the fallback (`:11-12`), but a documented silent degrade is still one (the E-02
precedent, §12.14). A NEAR-tie case also has no guard: a floating-point MAD of 1e-15 is accepted.
Why: the MAD is exactly 0 whenever at least half the law-whitened residuals are equal. That happens
with quantised observations, or with rows of zero redundancy (r_i = 0, e.g. arcs saturated by local
parameters), whose residual is 0. The std is the non-robust scale the MAD exists to replace: it is
inflated by the blunders it should expose.
Observable symptom: crossing 50 % ties flips the scale from 0.035 to 3.48 (std), and Huber then
keeps 0.235 of a 20-sigma blunder's weight instead of 0.002. With near-ties (1e-15 instead of 0)
the MAD branch is kept at s = 3.5e-15, and every ordinary row gets a factor of about 7e-15. Huber
never reaches 0, so the collapse guard (`fac > 0`) cannot fire.
Evidence: madfallback2.py: `zero-redundancy rows 49%: scale=MAD s=0.035 ... blunders (median)=0.002`
/ `51%: scale=std s=3.480 ... =0.235`; `near-ties 51%: scale=MAD s=3.48e-15 median Huber factor of
the ordinary rows=6.94e-15`. madfallback.py (quantised data, 75 % at the median):
`scale path=std`. The fallback has no test that fires it. `test_metamorphic.m:268` builds its
fixture specifically to AVOID it.
Severity: S2 (reachable on degenerate input, wrong weights without a message). Return the path
taken and warn once, or refuse.

### F5 — `Options.Noise.Smooth` and `Options.Noise.FreezeRedundancy` are documented knobs with no consumer; the user's GRAVFIRE recipe sets one — REPRODUCED
What is wrong: neither field is forwarded by `engineArgs` (`Options.m:255-261` forwards Grid,
SkipThreshold, MeanType, Pool and Tuning) and neither is read anywhere in the package.
`adjustVCEjoint` calls `psdVCE(..., 'Norm', o.Norm, 'Gamma', gammaIt)` (`:481`) without `'Smooth'`,
and it has no redundancy-freeze option at all.
Why: this is the WP-101 class ("a documented option with no consumer cannot ship again"). Audit
check 27 is silent because it enumerates only TOP-LEVEL properties (`audit_tree.py:1094`
`props = re.findall(r'^        ([A-Z]\w*)\s', ...)`), and the group name `Noise` is consumed. The
instrument cannot see group fields.
Observable symptom: `lsqvce.Options(Noise=struct('FreezeRedundancy', 100))`, exactly the command in
`dev/GRAVFIRE_INTERFACE.md:346` and the mapping row `:332`
(`minimumNumberObservations 100 → Noise.FreezeRedundancy`), runs bit-identically to the default.
The user believes low-redundancy components are frozen, and nothing tells them they are not.
`doc lsqvce.Options` (html:141-142) and Manual 04_options.md:115-116 describe both behaviours.
Evidence: groupfields.py: `Noise Smooth ... <-- NO CONSUMER`,
`Noise FreezeRedundancy ... <-- NO CONSUMER`. Every other one of the 23 group fields has a reader.
`grep -rn "FreezeRedundancy\|Smooth" +lsqvce` returns only Options.m and psdVCE's own option.
Severity: S2 (the user's choice is silently ignored, R5). Extend check 27 to group fields.

### F6 — `Results.meta.warnings` and the once-per-run warnings: wrong order on kinds B/C, and state that leaks between runs — REPRODUCED (mirror of the call order)
What is wrong: in `fit>runBatch`, `warnOnce('lsqvce:fit:priorNotActive', ...)` (`fit.m:274`) runs
BEFORE `runBatchOnBlockPath` resets the inventory (`fit.m:541`) or the x0 path does (`fit.m:289`).
The id is wiped from `R.meta.warnings`. Separately, `warnOnce` (`log.m:165`) consults the
session-persistent `FIRED` list, which is cleared only at the next reset. If a previous kind-A run
fired the id, a kind-B run's `warnOnce` before its reset is SUPPRESSED, and no warning is shown at
all.
Why: the doctrine is "warnings once per run + meta inventory" (R5 M1). The inventory is the
provenance record, and a "proceed with a warning" outcome turning silent is the R4/R5 failure.
Observable symptom: kind B/C with a `data.priors` group `from = 3`: the warning appears, and
`R.meta.warnings` is `{}`. After any kind-A fit that raised it, the same kind-B fit warns nothing.
The only test (`test_vcePolicy.m:472-477`) is on kind A, which resets first (`fit.m:172`), so it
passes.
Evidence: inventory.py: `run 1: warning raised=['lsqvce:fit:priorNotActive'] R.meta.warnings=[]`;
`kind A then kind B: kind-B warning raised=[] kind-B R.meta.warnings=[]`.
Severity: S2.

### F7 — the parallel log sink never reaches the workers: worker warnings miss `R.meta.warnings`, and `warnOnce` fires once per WORKER — INFERRED
What is wrong: `mapreduce.m:160` `logQ = lsqvce.util.log.startParallel();` stores the DataQueue in
the CLIENT's persistent `log.state('queue')` (`log.m:139`). `logQ` is never referenced inside the
`parfor` body (`:166-168` `raw{j} = feval(taskFcn, j)`). On a worker, `state('queue')` is that
process's own empty persistent, so `sendToQueue` returns false (`log.m:253-254`) and nothing is
forwarded. The documented "WORKER-AWARE LOGGING" (`log.m:49-57`) is therefore inert. Worker
warnings (`makeLoader:noCep` warnOnce, `covModel:negativePSD`, `separable:startupNotPD`,
`decorrFactor:notPD`, ...) never enter the client's `FIRED` inventory.
Observable symptom: `Parallel.Mode="parfor"` with ≥ Threshold systems or days: `R.meta.warnings`
differs from the serial run of the same data, and `noCep` prints once per worker.
`test_parallel` checks numbers, not the inventory.
Evidence: lines read as cited. `grep -rn "logQ" +lsqvce` gives the single assignment, unused.
Severity: S3 (provenance incomplete; the numbers are unaffected).

### F8 — `diagnoseAccum` materialises every whitened row at full width (n × (nc+Σnl)) and `Results` keeps it: "nothing of size n is ever formed" is false on every fit — REPRODUCED (arithmetic) / INFERRED (code)
What is wrong: `diagnoseAccum>whitenArcIntoBatches` builds `fullA = zeros(numel(keepRows), uFull)`
for every batch (`:217-221`) and appends it to `batches`. The array lives on in
`batchFcn = @(k) batchOf(batches, k)` (`:162`), which `lsqvce.Results` stores (`Results.m:180`).
Every driver routes through it: `fit.m:322`, `adjustBlock.m:176`, `adjustRobust.m:182`,
`adjustNonlinear.m:159`, `snoop.m:464`.
Why: handover F4 (§2.2, line 707) recorded only the dense Rfull (uFull²) and planned removal in M4
WP-86, which is still OPEN (§12.5). The n × uFull batch store is a larger, unrecorded part of the
same defect. §12.19 already carried "the peak-memory statement ... must be corrected ... in ONE
pass", and README.md:20 still says "nothing of size `n` is ever formed".
Observable symptom: GRAVFIRE-scale month (32 axis-series × 2160 epochs × 30 days; nc = 97²,
nl = 144/day): the batch store is about 228 GB, plus a 1.5 GB Rfull. Every `fit` ends in
out-of-memory after the solve has succeeded.
Evidence: `n=2073600 rows, uFull=13729, batches n x uFull x 8 B = 227.7 GB` (computed).
Severity: S2 at target scale (fails, not wrong); S3 for the false README claim.

### F9 — help text, an EXAMPLE and a runtime error message still name deleted units or retired behaviour (extends C3) — REPRODUCED
What is wrong: the sweep below matched every `lsqvce.<tier>.<unit>` / `LSQVCE.<TIER>.<UNIT>`
against the tree.
- Runtime ERROR MESSAGE (worse than help, because it is the remedy a user is given):
  `snoop.m:557` "re-run lsqvce.engine.adjust, or (c) use lsqvce.engine.robust".
- Live claim on `doc lsqvce.Options` and the Manual html: `Options.m:136-138` EXAMPLE says
  "Robust ~= "none" together with Estimate = "noise" is refused by fit (lsqvce:fit:notYet)". This
  has been false since 6.1 (WP-105); see `html/lsqvce.Options.html:177` and
  `html/lsqvce_manual.html:2279`. Check 28's seeded pattern `Not combinable with .Estimate="noise"`
  does not match this wording, so the ratchet is silent.
- `Results.m:96-98` ERRORS "use the result of adjustVCE"; `Results.m:760` "adjustVCE stores
  hist.whiteAlpha", although nothing writes `whiteAlpha` (grep: read-only), so the threshold line
  never draws.
- robust/adjust/Accumulator/solveNonlinear: `robustWeights.m:6,38,42`; `adjustRobust.m:60` (C3);
  `snoop.m:115,129,149,451`; `ellipse.m:20`; `ftest.m:32,77`; `diagnoseAccum.m:5,64`;
  `diagnose.m:20,95`; `datum.m:85`; `adjustNonlinear.m:14,49`.
- Units deleted in M2/WP-90/WP-92: `collectResiduals.m:33` (KERNEL.RESIDUALACV);
  `whitenessCorr.m:24` (ENGINE.SELECTARACV); `whittle.m:12`, `crossCorr.m:42,70`
  (ENGINE.SELECTARORDER); `satKeyOf.m:8,19` (RESOLVESATSPECS, NORMALIZESATLIST, and drivers
  ITERATESEPARABLE / GRACELOADER that exist nowhere).
- `declaration.m:18-20` "badParamColumns / ... (the fit help lists them)": fit's ERRORS block does
  not list `lsqvce:fit:badParamColumns` (grep count 0). It is raised at `declaration.m:138,143,147,315`
  and documented in no help block (ids.py).
Why: audit check 2's removed-name list (`audit_tree.py:160-162`, and SeparableProblem / recover for
help) predates the M2–WP-92 deletions, and check 24 inspects EXAMPLE blocks only.
Evidence: the sweep printed `MISSING 6 lsqvce.engine.adjust`, `MISSING 4 lsqvce.engine.robust`,
`MISSING 2 lsqvce.engine.adjustVCE`, `MISSING 2 lsqvce.engine.Accumulator`, plus 9 upper-case unit
names. The per-site list is above.
Severity: S3.

### F10 — documented NV set ≠ accepted set, a dead documented option, and two defaults for one fact (R1) — REPRODUCED
- `diagnose` documents `'verify' re-check sum(r) == n-u and vPv consistency (default true)` (`:43`).
  `o.verify` is assigned (`:123`) and never read; the integrity numbers are computed unconditionally
  (`:739-743`). `Verify=false` changes nothing. The same function accepts the undocumented
  `RobustCheck` and `Progress` (`:116,120`), although `RobustCheck` is a documented Options field.
- Defaults that disagree: `Options.Diagnostics.RMin = 0.01` (`Options.m:316`) against
  `diagnose opt.RMin = 1e-6` (`diagnose.m:113`, help `:31`). `fit` passes 0.01, but `snoop`'s
  internal re-diagnosis (`snoop.m:465-467`) passes neither and gets 1e-6. Within one workflow,
  `R.diag` and `snoop`'s D therefore classify "uncontrolled" observations differently.
- Accepted but undocumented in the help block: `adjustVCEjoint` Pool / Estimate / ArcFactors /
  OnPass / OnPassContext (`:286-303`); `adjustRobust` Prior / OnPass / OnPassContext;
  `adjustNonlinear` Prior / Diagnostics / OnPass / OnPassContext; `adjustBlock` Prior;
  `snoop` NBatch / MaxRejectFrac; `BlockAccumulator` LocalBlockSize; `Results.plot` Visible.
Evidence: nvdefaults.py output (17 rows); `grep -n "verify" diagnose.m` gives lines 23, 43 and 123
only.
Severity: S3.

### F11 — `diagnose` draws from the GLOBAL RNG for the MDB reservoir — REPRODUCED (grep)
What is wrong: `diagnose.m:654` `j = randi(nMdbSeen);` once more than 5000 finite MDBs have been
seen (`mdbCap = 5000`, `:207`). Every fit with n > 5000 advances the user's global stream, and
`R.diag.mdbQuantiles` (p50/p90/p99, reported by `report.m:538-548`) changes from run to run on
identical input. WP-36 (§F-06) established "the global RNG is never touched or advanced" for the
probes. This site was outside that sweep.
Observable symptom: `rng(1); a = rand; rng(1); fit(data); b = rand;` gives `a ~= b` for n > 5000.
Two identical fits report different MDB quantiles.
Severity: S3.

### F12 — R1 hygiene: the stream key and the arc key are composed inline and parsed back by regex — REPRODUCED (grep)
`matlab.lang.makeValidName(char(string(sysID)))` + `'|'` + obsType is built in 8 places (fit.m ×3,
adjustVCEjoint ×2, jointModelOpt, applyCovOverride, Results), while `engine.satKeyOf` (a third, `s_`
-prefixed convention) has no caller in the package. The arc key `'s%s_%s_d%s_a%s'` is formatted in
`adjustVCEjoint.m:820` and in `private/makeKeyFcn.m:39`, and recovered by regex
`'^s(.+?)_(.+)_d[^_]+_a'` in `fit.m:626` and `adjustVCEjoint.m:1000,1016`. coding.md R1 says
"Cache keys are structs, not concatenated strings parsed back apart". F1 is the failure this
pattern produces.
Severity: S4.

### F13 — length and line-length rules not met, without a record — REPRODUCED
Per FUNCTION (not per file): `report` 826 lines (592 code), `diagnose` 1001 (496 code),
`adjustVCEjoint` main 632, `makeLoader` main 489, `snoop` 486. HANDOVER §6.3 (`:2358`) requires "no
function longer than ~400 lines without a written justification in this handover", and grep finds
none. 74 lines exceed the MathWorks 120-character limit (188 exceed 100). The worst files are
makeLoader (16) and Options (7). `makeLoader.m:308` is a single 600-character struct literal.
Severity: S4.

## Checked and sound
- Layering: audit check 21 gives 4 upward edges, all in the recorded baseline. A grep of
  `+kernel` / `+util` for `lsqvce.engine.` finds only `toeplitzCofactor → engine.covComponent`
  (baselined).
- Unqualified sibling calls in packages: none (a grep for bare `<unit>(` of every package unit; the
  only hits are MATLAB's `log()`).
- `eval` / `evalin` / `assignin` / `global`: none in `+lsqvce`. Empty or bare catches: 0 (check 19).
- Identifier catalogue: 301 raised = 301 catalogued (check 20). The 11 ids raised by
  `declaration` / `normalizeSystems` are documented on their public fronts (fit / adjustBlock),
  except `badParamColumns` (F9).
- Probe RNG: `probeMatrix` / `probeStream` use private RandStreams and the global stream is not
  touched (read). The seed is recorded as given; an unset seed is recorded as [] by design.
- Warning-state restore: `log.warn` uses `onCleanup(@() warning(st))` held in a function-scope
  variable; mapreduce's queue cleanup also uses onCleanup. Both are sound idiom.
- containers.Map state across passes: `adjustRobust` keeps one W written per pass and keyed per arc.
  `adjustVCEjoint>robustSweep` builds a fresh Map each sweep. `epochResiduals` undoes the folded
  factors from the same W it wrote. Sound, apart from the key format (F1).
- Prior-stream exclusion is consistent in both robust drivers (`lawWhitenedResiduals.m:55`,
  `adjustVCEjoint.m:695`). A data stream named `prior_<x>` is refused only when it equals a
  declared group's name (`makeLoader.m:437-441`). Recorded behaviour, not reopened.
- parfor/spmd: `mapreduce` folds parfor results in input order, and spmd only with an associative
  Combine. The accumulator forces `DistributeDays=false` inside system-level parallelism (no nested
  pools). Read; test_parallel covers the numbers.
- `Options`: the arguments contract matches the help for top-level properties. Retired knobs are
  refused at assignment (WP-101). Group merge semantics are as documented (Manual 04:105-106).
- `epochWeights` / `rowWeights` length guards exist on read (`makeLoader.m:1264,1317`).

## Not reviewed
- Runtime behaviour on MATLAB (tier B): no suite run, no codeIssues. F2, F7 and F8 are inferred
  from code paths and carry no execution evidence.
- Estimation-theory correctness of the Ellmer update, Koch/Yang conventions and the Helmert prior
  update (other lane). F3 states only the measured code fact.
- `report.m` message content and the Manual / theory PDFs beyond the stale-reference sweep.
- `BlockAccumulator` Schur/QR numerics (H-16 closed and recorded) and `QRSystem` / `NormalEqSystem`
  internals. Only interface-level reads were done.
