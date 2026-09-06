# CLAUDE.md — session entry point for lsqvce

Stable orientation for a Claude session: purpose, architecture, decisions,
conventions, pitfalls. STATUS lives in ONE place only: `dev/HANDOVER.md`
(the SSOT — §12.11 the R5 closing summary, §12.12 the parameter-table
interface D8, §12.13 the docs level; the ledger rows and round notes carry
every decision with its reason). Do not create a third document. The single
documentation authority for users is the Manual (`doc lsqvce`).

## Purpose

`lsqvce` is a MATLAB toolbox (R2026a; PCT optional) for least-squares
adjustment with variance component estimation, built for satellite gravity
recovery (GRACE / GRACE-FO, the user's GRAVFIRE pipeline) and equally usable
on a levelling network. One user verb: `R = lsqvce.fit(data, Name=Value)`.
Version 5.2 (Sept 2026) is the result of the R5 review track (PRs #48–#84):
one accumulator, one stochastic-model estimator, one loop with three
policies, one user surface, regularisation declared on the data.

Repository `matweigelt/lsqvce` (PUBLIC — measurement data and secrets never
enter it). Matthias is the sole developer and merger; Claude is the
implementation delegate (branch → commit → PR → poll CI → verdict).

## Architecture (tier rule: user → engine → kernel/util, never upward; `dev/docs_src/tiers.csv`)

- **User tier** `+lsqvce/`: `fit`, `Options` (14 knobs + 5 groups), `Results`
  (`report`, `summary`, `convergence`, `residuals`, `contribution`, `compare`,
  `plot`, `toTable`), `report`, `layout`.
- **Engine** `+engine/`: `declaration` (validate the data struct, build the
  loaders — the ONE declaration step, reused by `layout` / `toTable`);
  `makeLoader` (per-system day loader: scoped parameter plan, validity,
  priors as the observation type `prior_<name>`); `BlockAccumulator` with
  `QRSystem` / `NormalEqSystem` (block-arrow accumulation, local blocks
  eliminated per system and day); the drivers `adjustBlock` (one GLS solve),
  `adjustRobust` (IRLS), `adjustNonlinear` (damped Gauss–Newton),
  `adjustVCEjoint` (the joint Ellmer estimator: PSD-domain law per pool group,
  level per stream, per-arc factors, whiteness collector, `lastSol` to the
  accumulator factory for prior updates); `jointModelOpt`, `covModel`,
  `batchSystem` (kinds B/C as one local-free system); `diagnoseAccum`,
  `diagnose`, `residuals`, `redundancy`, `outliers`, `snoop`, `ftest`;
  `CommonPrior`; `localColumns` (the one column unroll of a local plan).
- **Kernel** `+kernel/`: `redundancyLags` (FFT lag sums), `armaWhiten`,
  `armaACV`, `psd2acv`, `dctI`, `levinsonDurbin`, `robustWeights`, `whittle`.
- **Util** `+util/`: `callUser` (wraps every user function; original error as
  cause), `log` (`warnOnce`, inventory → `R.meta.warnings`), `onPass`,
  `toolboxVersion` (the one version reader), …
- **Around it**: `examples/` (demos incl. `demo_gravfire`), `template/`
  (`gravfire_data.m`: the GRAVFIRE bundle → `data`, STEP 1–7),
  `html/` (the Manual + one page per public unit), `lsqvce_theory.pdf`
  (R5 theory, 54 pp.), `lsqvce_manual.pdf`, `RELEASE_NOTES.md`, `dev/`.

## Key design decisions (R5; reasons in HANDOVER §12)

- **Three data kinds on one path**: (A) streaming — `commonParams`, `days`,
  `systems(k).{buildDay, covariance, arcs, localParams, priors,
  processingOptions}`, optional `priors`, `paramColumns`; (B) in-memory `A, l,
  p` (p diagonal or full SPD); (C) batched `batchFcn`, `numUnknowns`,
  `numBatches`, optional `x0` (Gauss–Newton). All ride the block accumulator.
- **Parameter tables are n×15 cells**: column 1 = scope token (`'common'|'day'`
  for commons; `'arc'|'day'|'system'` for locals), columns 2–15 the user's (his
  14 GRAVFIRE columns). The toolbox derives counts and the block layout ONLY;
  the user builds A and B in `buildDay` with `info` (offsets, rows present).
  `data.paramColumns = struct('start','end','name')` names the validity and
  name columns; a row is instantiated only inside its validity.
  `lsqvce.layout(data)` announces the columns before fit; `R.toTable(data)`
  writes estimate / σ back into the expanded tables.
- **Regularisation is declared, never an option** (`Options.Reg*`, `DatumB`,
  `Sigma0`, `Sparse` are gone). `data.priors(j)` (common, L in the columns of
  A) and `data.systems(k).priors(l)` (local, L in the columns of B) are ONE
  mechanism: one struct per group, declared by `group` (a regexp on the table's
  name column; the toolbox builds L, one row per matched row per instance) OR
  by `L`, never both; `x0`, ONE `sigma` per group (relative weights inside a
  group go into L), `from` (activation pass), `name`, `estimate="none"|"scale"`
  (Helmert / Koch–Kusche per pass; several rows make σ estimable; local groups
  ride as the stream `prior_<name>`, so their σ is the stream level — identical
  to common). A datum is the prior on the null space (exact, weight-free).
  THE RULE: parameters and their columns are the user's; regularisation is the
  toolbox's; the interface is an operator in the user's columns, checked by
  width.
- **`Estimate="none"|"noise"`, `Robust`, `data.x0` are policies of the same
  loop**; the v4 Models map to `Pool` / `ArcFactors` (RELEASE_NOTES table).
- **`Options.OnPass`** fires after every pass with the pass state, the live
  accumulator, `data` and `processingOptions`; `false` stops
  (`stopReason="userStop"`). An opening-brief requirement — dropped once in
  M4, restored in M5; a review-plan item never outranks the brief.
- **Error doctrine**: every catch refuses / reports / counts; identifiers
  catalogued (`dev/docs_src/error_ids.csv`, closed both ways); refusals are
  thresholds the code chooses (QR rank by `rank()` tolerance), never the
  factoriser's behaviour at the boundary; NaN-honest v'Pv; warnings once/run.
- **Oracles over parity**: `demo_gravity` / `test_physOracle` (planted truth,
  z-scores) are counted tests; they found the defects parity tests could not.
- **Docs chain enforced by the audit**: help block → API page (check 22) →
  Manual (check 23: HTML = fresh build, every Options knob / group field
  described) → identifier catalogue (check 20). Manual chapters 5 (API), 10
  (tiers), 12 (migration), A (identifiers) are GENERATED.

## The three commands

```matlab
lsqvce_setup                    % install (root + examples, Help index, savepath); 'check' | 'uninstall'
cd('dev/tests'); ci_gate        % THE gate (never bare runtests); manages its own path
```
```bash
python3 dev/docs_src/audit_tree.py          # static audit, exit 0 before every commit
python3 dev/docs_src/audit_selftest.py      # when the audit itself changed
python3 dev/docs_src/build_missing_pages.py all   # regenerate the API pages from the help blocks
python3 dev/docs_src/build_manual.py        # rebuild the Manual (pandoc + xelatex; CI installs pandoc)
```

## Conventions agreed

- **Dev loop** (HANDOVER §11.5): one branch per level/round from `main`; audit
  green before every commit; CI green on the head before the ready-flip
  (draft PR → poll CI in 240–295 s chunks → GraphQL
  `markPullRequestReadyForReview` → verdict comment); Matthias merges;
  `git checkout main` after. Commit identity `-c user.name="Claude"
  -c user.email="noreply@anthropic.com"`. Every touched file gets a dated
  stamp naming the model and work package. A wrap after every change.
- **Verification tiers**: A = MATLAB MCP available (run `ci_gate`, the
  examples and the GRAVFIRE bundle `C:\Users\matth\Documents\MATLAB\testdata.mat`
  locally); **B (current)** = no MCP: CI runs the MATLAB gate
  (`test_parallel` excluded, it needs the user's PCT); the user runs `ci_gate`
  + `demo_gravfire` locally before merging.
- **Token**: `/mnt/project/GitHub_Access_Token.txt` (name may change; `ls
  /mnt/project`); read at call time with `tr -d '[:space:]'`, never echoed,
  never committed, never in a logged URL.
- **Help text**: H1, usage lines, INPUT / OUTPUT with type and dimensions,
  ERRORS with identifiers, See also, `Created <date> (CEST) with assistance
  from Claude (Anthropic), model <name> -- <level/WP>`; the help block ABOVE
  any `arguments` block; no line starting with an all-caps word matching a
  page title. MATLAB coding standards
  (github.com/matlab/rules/blob/main/matlab-coding-standards.md); no aliases;
  one name per thing.
- **Tests**: `functiontests`; every point declares `% categories:` (contract,
  precision, reference, robustness, vectorisation, speed, metamorphic — every
  public unit has each category tested or argued in
  `dev/tests/category_exemptions.csv`); every tolerance has a `% tol:`
  provenance line; recovery asserts z-scores against formal errors; speed
  budgets are ratios by `timePair` pinned to their own measurement; provoked
  warnings muted with `muteWarn` (gate rule 7); a filtered point is a red
  flag, not a pass; the inventory of a rewrite is the tree (glob), never a
  hand-picked list.
- **Docs**: RELEASE_NOTES entry per version; HANDOVER ledger row + dated
  round note per round (the same PR); Manual sources in
  `dev/docs_src/manual/*.md` (numbered chapters); references (Ellmer 2018,
  Kvas 2019) in `dev/references/README.md`, full texts not redistributed.

## Status snapshot (2026-09-06 — the authority is HANDOVER §12)

- `main` at 5.2 + docs level; PR #84 (setup-check fix) open, ready → tag
  `v5.2.1` after it. Gate: 26 suites, ~527 points, green; CI ~4 min.
- Open: GRAVFIRE STEP 11 (the user's design function behind `buildDay`);
  `Parallel.Axis="arcs"` gated on a measurement at u > 6000 (PCT runner,
  WP-30); self-hosted Windows runner for `test_parallel` / `test_speed` on
  CI; a possible second review track on 5.2.

## Known pitfalls

- Arc keys are `s<sys>_<obsType>_d<day>_a<arc>`; observation types may
  contain underscores (`prior_<name>`) — parsers anchor on `_d[^_]+_a`.
- A local prior row must lie in one day's local block (system-block rows
  count as day 1); couple days through common parameters.
- `exist('pkg.fn','file')` is 0 for package functions — use `which`.
- `builddocsearchdb` creates `helpsearch-v4_<locale>` (git-ignored).
- The loader caches one built day (`Options.CacheDays`); a re-linearising
  source (`State`) is never cached.
- Correlated batches: per-observation diagnostics are in the whitened frame;
  only sums are frame-invariant.
- The QR fold's v'Pv is NaN-honest (an exact fit with ltl > 0 warns
  `qr:cancellation`); Gauss–Newton uses l'Pl for step control.
- `Options.Noise.Whitener` default is `"toeplitz"` (the `"ar"` reduction
  diverges on short PSD grids).
- `struct()` with a cell value expands into a struct array — braces; never
  name a local `keys`; no comment lines inside `...`-continued expressions;
  `containers.Map.Count` is uint64; struct arrays need one field set.
- Regex migrations across tests overreach onto engine-tier argument lists —
  check the `lsqvce.Options(` / `fit(` context.
- Do not resurrect the deleted v4 manual generator (`master_content_*`).
