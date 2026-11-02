# The examples

The folder `examples/` holds runnable demonstrations of the user path:
declare a `lsqvce.TaskDefinition`, choose `lsqvce.Options`, call `lsqvce.fit`, and read
the `lsqvce.Results`. The folder `template/` holds starting points that you
copy into your own project. Every example generates its own synthetic data
from a seeded random stream, except `demo_gravfire`, which reads a private
data bundle that is not shipped. All of them print their findings to the
Command Window. A few of them also return their results so that you can
inspect them afterwards.

`lsqvce_setup` puts `examples/` on the path. It does not add `template/`:
copy those files into your project. `demo_gravfire` adds `template/` to the
path by itself.

## Which example for which question

| Question | Example | Kind of data |
|---|---|---|
| I am new: what does a complete streaming run look like, and where do `Prepare` and `Build` go? | `demo_tutorial` | streaming, one system, two streams |
| How do I wire my own data into a task and read the verdict? | `demo_workflow` (script) | batched; streaming recipe |
| What do the three data kinds look like in one call each? | `demo_fit` | in-memory, batched, streaming |
| What does the report say when an adjustment is wrong in a specific way? | `demo_lsqvce` (script) | in-memory, one streaming case |
| How do I tell the toolbox what I know about each stream's noise? | `demo_noisePolicy` | streaming, one system |
| How do I regularise, or give a rank-deficient problem a datum? | `demo_priors` | streaming, one system |
| How are several systems with their own local parameters combined? | `demo_combine` | streaming, three systems |
| How do I declare system / sensor / day, with data gaps, and find a bad day? | `demo_sensors` | streaming, two systems |
| What does a full noise-estimation run on multi-axis streams look like? | `demo_grace` (+ `demo_grace_dataLoader`) | streaming, two systems |
| Does the noise estimation still work on a physically structured gravity design? | `demo_gravity` (+ `demo_gravity_dataLoader`, `demo_gravity_design`) | streaming, three systems |
| How do I map a GRAVFIRE preprocessing bundle onto a task? | `demo_gravfire` + `template/gravfire_data.m` | streaming, real bundle |
| How do I build GRACE-FO range-rate observation equations from Level-1B files? | `template/kbr_variational.m` | streaming, real data |
| What does each verbosity level print? | `demo_messaging` | streaming |
| Should I use `Backend="normal"` or `Backend="qr"`? | `bench_qr_vs_normal` | batched |

## Suggested reading order

0. **`demo_tutorial`**: the Manual's chapter *Your first adjustment*, as a file.
1. **`demo_workflow`**, then **`demo_fit`**: the three constructors of a task, `lsqvce.Options`, and how to read `R.diag` and `R.report()`.
2. **`demo_lsqvce`**: what the diagnostics say when something is wrong.
3. **`demo_noisePolicy`**: the per-stream noise policy (`Estimate`,
   `ArcFactors` on `addStream`).
4. **`demo_priors`**: `task.layout()` and declared regularisation.
5. **`demo_combine`**, then **`demo_sensors`**: several systems, local
   parameters, data gaps, attribution by system, sensor and day.
6. **`demo_grace`**, then **`demo_gravity`**: complete noise-estimation runs
   on GRACE-like problems.
7. **`demo_gravfire`** / **`gravfire_data`** and **`kbr_variational`**: the
   templates for real data.
8. **`demo_messaging`** and **`bench_qr_vs_normal`**: what a run prints and
   which backend to use, on your own machine.

## Examples, helpers and templates

Not every file in `examples/` is a demonstration. Three of them are
**helper files** that the demos call:

| Helper | Called by | Returns |
|---|---|---|
| `demo_grace_dataLoader` | `demo_grace`, `demo_messaging` | `[cfg, data, task]`: a GRACE-FO-shaped task with a random design (`data` is its compiled declaration, for the engine-tier tests) |
| `demo_gravity_dataLoader` | `demo_gravity`, `demo_fit` | `[cfg, data, task]`: a gravity task with a physical design and planted truth (`data` as above) |
| `demo_gravity_design` | `demo_gravity_dataLoader` | `[A, meta]`: spherical-harmonic potential partials along an orbit |

Each helper also runs on its own and returns a ready task that you can pass to
`lsqvce.fit`.

`demo_workflow` and `demo_lsqvce` are **scripts**, not functions. They
begin with `clear`, so they wipe the base workspace. Run them in a fresh
session, or section by section in the Editor. All other examples are
functions.

## How to run them

The call forms below are the ones the test suite uses to run every file in
`examples/` once (with figures hidden). Each call is the cheapest form that
still runs the whole example:

```matlab
bench_qr_vs_normal('u', [10 20], 'n', 2e4, 'nk', 2000)
demo_combine()
demo_fit()
demo_grace('verbose', false)
demo_grace_dataLoader()
demo_gravfire(bundleFile)             % the suite writes a generated bundle to bundleFile
demo_gravity()
demo_gravity_dataLoader()
orbit = struct('u', linspace(0, 4*pi, 400).', 'lam0', linspace(0, 0.3, 400).', ...
               'incl', 89*pi/180, 'r', 1.07);
demo_gravity_design(4, orbit)
demo_lsqvce                           % script
demo_messaging()
demo_noisePolicy()
demo_priors()
demo_sensors()
demo_tutorial()
demo_workflow                         % script
```

`demo_gravfire` runs on a bundle the suite generates
(`dev/tests/fixture_gravfire_bundle.m`, random numbers in the real structure), so
`template/gravfire_data.m` runs with it. `demo_gravfire` with a path that does not
exist prints what it would need and returns. `kbr_variational` needs Level-1B files
and is not part of the run.

---

## `demo_tutorial`: your first streaming run, with `Prepare` and `Build`

```matlab
[R, s] = demo_tutorial();
[R, s] = demo_tutorial('verbose', true, 'plot', true);   % defaults: false, false
```

The file of the chapter *Your first adjustment*, line for line. One system watches a
periodic signal with two instruments (`"gauge"`, white; `"radar"`, AR(1) but declared
white) on two arcs a day for five days; every arc has an unknown offset, the `"arc"` local
`bias`. The system's `Prepare` (`prepareTutorialDay`) computes the epochs and the signal
basis of a day once, and both streams' `Build` read them from `ctx.Prepared` and give the
bias partial by name (`Local={"bias", ones(m, 1)}`).

It runs `fit` twice — once under the declared noise, once with `Estimate = "noise"` — and
prints the two estimates of the two common parameters against the planted truth. Both
are close to the truth; the declared run's formal errors are about 1.6 times too small,
because the radar's colour was missing. The noise run converges in four passes, σ₀ = 1.000,
and the report's verdict is `PASS`.

`s` returns the truth, both runs' estimates, standard errors and σ₀, each stream's final
level relative to the declared one, and the verdict. The test suite checks these claims
(`test_fit`, `testDemoTutorialIsTheManualsTutorial`).

## `demo_workflow`: build your own case

**Purpose.** This is the from-scratch recipe for wiring your own data into a
`lsqvce.TaskDefinition`, choosing what to estimate, and reading the result in the right
order.

**Run.** `demo_workflow` is a script with no arguments and no outputs. It
runs `clear; clc; close all`. The help text recommends running it section by
section (Ctrl+Enter).

**Outline and logic.** The script is organised as Sections 0 to 5:

- *Section 0* lists the three ways to make a task:
  - in memory: `lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p)`;
  - batched: `lsqvce.TaskDefinition.fromBatches(@batchFcn, NumBatches=K,
    NumUnknowns=u)`, or `X0=x0` for the nonlinear form;
  - streaming: `lsqvce.TaskDefinition(Days=days)`, filled with `addCommon`,
    `addSystem` (its `Prepare`), `addLocal` and `addStream` (its `Build` and `Noise`).

  It also states the batched callback contract
  `[A, l, p] = batchFcn(k)`, with `p` the weights `1./sigma.^2`. The
  callback is called twice per batch and must return identical data both
  times.
- *Section 1* synthesises 3000 observations with 8 unknowns and cuts them
  into 30 batches with a closure over precomputed batch edges. It
  sanity-checks `batchFcn(1)` with an `assert` and makes
  `task = lsqvce.TaskDefinition.fromBatches(batchFcn, NumBatches=nb, NumUnknowns=u)`.
- *Section 2* is a decision tree in comments: which kind, and when to use
  `Estimate="noise"`. It then runs `R = lsqvce.fit(task)` with default
  options.
- *Section 2b* reads the verdict in a fixed order:
  1. integrity, `R.diag.sumRRelErr`;
  2. variance, `R.diag.varRatio`;
  3. blunders, `R.diag.nFlag`;
  4. reliability.

  The one-line verdict comes from `R.report('Verbose', false).verdict`.
- *Section 3* is a recipe only, with no runnable code. It shows a streaming task with
  `addCommon`, one `addSystem` with a `Prepare`, scoped `addLocal` calls and one
  `addStream` per observation type with its `Build` and `Noise`, and points to
  `demo_tutorial`, `demo_combine`, `demo_sensors` and `demo_grace`.
- *Section 4* plants a second instrument group that is three times noisier
  than claimed. It declares the two groups as two streams `g1`, `g2` of one system
  `"groups"` on one day, with commons `x1` … `x8`. Each stream's `Build` (the local
  function `groupArc`) adds one arc and passes the claimed per-row sigmas as `Sigma`;
  both streams declare `Noise=lsqvce.Noise.white(1), Estimate="scale"`. The run
  uses `Options(Estimate="noise", Pool="stream", ArcFactors=false)`. The
  script prints the recovered sigma ratio from `Rv.history.sigma(end,:)`
  (the truth is 3) and `Rv.sigma0`.
- *Section 5* prints the condensed decision tree and the rule for reading
  a result.

## `demo_fit`: the three data kinds, one call each

**Purpose.** This example shows every kind of task that `lsqvce.fit`
accepts, end to end on the user tier: `lsqvce.fit`, `lsqvce.Options`,
`lsqvce.Results` and `R.report`. It also shows robust estimation, a datum
declared as a prior, Gauss–Newton through `X0`, and watching a run with
`OnPass`.

**Run.**

```matlab
out = demo_fit();                  % prints the verdicts
out = demo_fit(Verbose=false);     % silent
```

`out` has the fields `.matrix`, `.robust`, `.datum`, `.batch`,
`.nonlinear`, `.gls` and `.vce`. Each of them is an `lsqvce.Results`.
`out.truth` holds the planted values: `x = [1; -2; 0.5]`,
`[a; b] = [0.7; -0.2]`, and `xGravity` of the gravity loader.

**Logic.**

1. **In-memory.** The data are 30 observations with
   `A = [1, randn, randn]`, seeded with 17, and the task is
   `lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p)`. The same problem is then run with a
   25-sigma gross error at row 7 and `Options(Robust="huber")`. The demo
   prints the IRLS weight factor of that row from
   `out.robust.history(end).factor(7)`.
2. **Free network.** The problem has height differences between four
   points and is rank-deficient by a common level. The datum is declared
   as `task.addPrior("datum", L=ones(1,4), Value=0, Sigma=1)`. The demo
   prints `L*x`. It then fits the same matrices without the prior and prints the
   error identifier with which `fit` refuses the problem.
3. **Batched, linear.** The same 30 rows are split into two batches with
   `lsqvce.TaskDefinition.fromBatches(bf, NumBatches=2, NumUnknowns=3)`.
4. **Batched, nonlinear.** The model is `l = exp(a t) + b`. The callback
   `bfn(k, x)` returns the Jacobian at `x`, and the task is
   `fromBatches(bfn, NumBatches=2, X0=[0.5; 0])` (the number of unknowns is `numel(X0)`). The options are
   `Options(MaxPasses=15, Tolerance=1e-9)`, and the demo prints `[a; b]`.
5. **Streaming.** The task of `demo_gravity_dataLoader()` is solved once
   with defaults (`out.gls`). It is then solved with
   `Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=2, OnPass=watch)`,
   where `watch` prints `ctx.pass`, `ctx.maxDev` and `ctx.sigma` after
   every pass. Two passes are not enough to converge, so the demo
   suppresses the warning `lsqvce:adjustVCEjoint:notConverged`. It prints
   `out.vce.history.sigma`.

**Implementation.** The local function `show` calls
`R.report('Verbose', true)` and prints `R.sigma0`. According to the help
block, the in-memory estimates land within 3 formal sigma, and the
nonlinear form recovers `[a; b]` to 5e-3.

## `demo_lsqvce`: six failure modes and their verdicts

**Purpose.** This example teaches you to read the diagnostics. Each of six
scenarios isolates one failure, so that you see what `R.diag` and
`R.report()` say when one specific thing is wrong, and what they say when
nothing is.

**Run.** `demo_lsqvce` is a script with no arguments. It runs `clear; clc`.
It draws one figure in case 3.

**Logic.** All cases share a setup of 2000 observations, 8 unknowns and
per-row sigma in [0.5, 1.5], seeded with 2026. Each case fits an in-memory task,
`lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p)`:

1. **Clean:** the baseline.
2. **Outliers:** three gross errors at rows 137, 892 and 1455. The case
   prints `R2.diag.flagIdx` and a note on rejecting only the largest `|w|`
   at a time.
3. **Systematic error:** an unmodelled `1.2 sin(2 pi t)`. The case prints
   `nFlag` and plots `R3.diag.w` against time with the `±kCrit` lines,
   above the clean case for comparison.
4. **Uncontrolled observation:** unknown 8 is seen only by observation 500.
   The case prints `nFlag` and the redundancy number `rNum(500)`.
5. **Wrong weights:** group 2 is three times noisier than claimed. The case
   prints `varRatio` and `nFlag`, then shows the cure: a streaming
   task with two streams `g1` and `g2`, each declaring its group's per-row sigmas as
   `Noise=lsqvce.Noise.white(sigma)` (its `Estimate` default is then `"scale"`; local
   helper `groupDay`), and `Estimate="noise", Pool="stream", ArcFactors=false`. It
   prints the recovered sigma ratio (truth 3) and `sigma0`.
6. **Collinear parameters:** column 2 is almost column 1. The case prints
   `R6.diag.corrMax`.

**Output.** A closing summary prints one line per case with
`R.report('Verbose', false).verdict`.

## `demo_noisePolicy`: what you know about each stream

**Purpose.** This example shows the per-stream policy of `addStream`.
`Estimate` is `"none"`, `"scale"` or `"law"`, and `ArcFactors` is `true` or `false`. Each choice recovers a planted truth.
For the policy semantics, see the chapter *Choosing the stochastic model*.

**Run.**

```matlab
[R, s] = demo_noisePolicy();
[R, s] = demo_noisePolicy('maxPasses', 20, 'verbose', true);   % defaults: 10, false
```

`R` is the `lsqvce.Results`. `s` summarises what was measured:

- `.stream`, `.declared`, `.level` and `.truth`, per stream;
- `.noneFixed`;
- `.rho1`: the estimated and planted lag-1 correlation of `rr`;
- `.arcFactor` and `.arcRatio`;
- `.passes`;
- `.badEstimate`: the identifier with which `addStream` refuses a misspelt policy.

**The case.** The data have one system, three days and three arcs a day of
256 epochs (5 s), four common harmonic coefficients, and one `'arc'`-scoped
bias. There are four streams:

| Stream | Truth | Declared | Policy |
|---|---|---|---|
| `gps` | white, sigma 0.2 | `Noise.white(0.2)` (correct) | `"none"` |
| `acc` | AR(1), phi 0.85, level 0.5 | `Noise.psd(f, S)` of the same filter, level 1 | `"scale"` |
| `rr` | AR(1), phi 0.5 | `Noise.white(1)` | `"law"` |
| `mag` | white 0.3, arc 2 three times worse | `Noise.white(0.3)` | `"scale"`, `ArcFactors=true` |

**Implementation.**

- `policyTask` declares the task: one system `"policy"` with the `"arc"` local `bias`
  and one `addStream` per stream, with its `Noise`, `Estimate` and `ArcFactors`. For
  `acc` it builds a `psd` law from `firstOrderPsd`, which gives the one-sided PSD on a
  frequency grid.
- Every stream has its own `Build` (`buildGps`, `buildAcc`, `buildRr`, `buildMag`); all
  four call `policyDay`, which places the bias column of each arc with
  `ctx.cols("bias", Arc=id)` and adds the arc with `B=`.
- `noiseSource` draws the planted noise from a stream seeded on day and
  stream, so that every `Build` is deterministic.
- All four policies run in one call with
  `Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=o.maxPasses, Noise=struct('Smooth', 1))`.
  `ArcFactors=false` is the run default. `mag` overrides it on its own
  `addStream`. The loop converges in 3 passes. Without the smoothing, the 60-bin
  law of `rr` needed 14 passes and settled on the same values to three
  digits.

**Output.**

- A table of the level per stream (declared, estimated, planted) from
  `R.noise.model{k}.sigma`, with streams looked up by name.
- The lag-1 correlation of `rr` from `R.noise.Cgroups`.
- The mean arc factor per arc of `mag` from `R.noise.model{k}.sigma2`, and
  the ratio of the bad arc to the good arcs (planted: 9).
- The identifier with which `addStream` refuses `Estimate="shape"`
  (`MATLAB:validators:mustBeMember`), before any data is read.
- `R.summary()`.

The help block records the measured values of one run. For example, `gps`
keeps its declared 0.2 exactly, and the bad-arc variance ratio is 9.45
against the planted 9. Only the ratio between arc factors is identifiable,
because the level absorbs their mean.

## `demo_priors`: regularisation in the announced layout

**Purpose.** In a scoped declaration you do not know in advance which
column of `A` or `B` an unknown ends up in. `task.layout()` announces
the layout, with every column's name, before `fit` runs, so you can write a prior operator `L` in
those columns. For the full reference, see the chapter *Priors and regularisation*.

**Run.**

```matlab
[R, s] = demo_priors();
[R, s] = demo_priors('verbose', true);
```

`R` is the regularised solve. `s` holds:

- `.nc` and `.nl`;
- `.nullDim`;
- `.unregularised`: the error identifier;
- `.sigma0` and `.zSignal`;
- `.nPseudo`, `.formGap`, `.datumGap` and `.priorLevel`.

**The case.** The data have four harmonic signal coefficients (`'common'`)
and one offset per day (`task.addCommon("offset", Scope="day")`). The system has
two `"arc"` locals (`bias`, `drift`) and one `"system"` local (`tempScale`, a scale on a
known auxiliary series). There are three days with two arcs of 150 epochs each.
The day offset and the arc biases cannot be separated, so the normal
matrix has one null direction per day.

**Logic.**

1. Call `lay = task.layout()`. Print `lay.nc`, `lay.systems(1).nl`
   and the column map from the names in `lay.common` and `lay.systems(1).local`, and
   find the day-offset columns with `lay.colCommon(row, day)`.
2. Count the null directions on a densely assembled `[A B]` (`rank`). Fit
   without priors and record the refusal identifier.
3. Declare two common groups with `task.addPrior`: a tight `"datum"` (the
   offsets average to zero) and a loose `"smooth"` (first differences),
   both `L=` operators with `Estimate="none"`.
4. Declare a local group `"drift"` with `sys.addPrior`, written
   once as `L=` through `lay.systems(1).col(row, day, arc)` and once as the
   pattern `On="^drift$"`. Both solutions are compared.
5. Move the datum (`Value = 0.5`) and measure that the signal stays put while
   the offsets move.
6. Run `Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=4)`
   with the drift group at `Estimate="scale"`. The group appears as its
   own stream `x1|prior_drift`, and its level is printed.

**Implementation.**

- `calibTask` declares the instrument without priors; every variant (no priors, the
  datum, the `L` form, the `On` form, the moved datum, the noise run) is a task of its
  own, and `addDatum` adds the two common groups.
- `buildCalib` is the stream's `Build`. It takes every column from `ctx.cols` — the
  commons, this day's offset, `tempScale`, and each arc's `bias` and `drift` with
  `Arc=id` — and hands them to `priorArcs`, which builds the arcs with `B=`.
- `nullDimension` assembles the dense rank check, calling the same `priorArcs` with the
  columns of the layout. That check is for demonstration only.
- `fitWithoutGroups` gives the pseudo-observation count.
- `R.coefficients().se` gives the z-scores of the signal.

## `demo_combine`: several systems by elimination

**Purpose.** This example shows how to declare several systems that share
common parameters, where each system also has local parameters that only
it observes. `lsqvce.fit` eliminates the locals system by system. The
example also shows how to check the degrees of freedom. A parameter is
local only if exactly one data set observes it.

**Run.** `demo_combine()` takes no arguments and returns nothing. It runs
`clc`.

**Logic and implementation.**

- The data have six common unknowns and three systems with 12, 8 and 15
  `'system'`-scoped locals. The systems have 40, 30 and 50 days of 25
  observations each, with unit noise, seeded with 2026.
- The task has no calendar of its own; each `task.addSystem(name, Days=1:nDay(i))`
  carries its own, with its `"system"` locals (`addLocal(names, Scope="system")`) and one
  stream `"obs"`, `Noise=lsqvce.Noise.white(1)`.
- Each stream's `Build` comes from `makeBuild`, which reads the pre-generated rows, one
  arc per day, and places the locals with `ctx.cols(localNames)`.
- One call `lsqvce.fit(task, lsqvce.Options(Verbose=false))`.

**Output.**

- An integrity block that compares `n - u_total` (the eliminated locals
  included) with `R.diag.r`.
- A table of the common parameters: estimate, truth, error, and sigma from
  `R.sol.Cxx`.
- `R.sigma0`, `R.diag.varRatio` and `R.report`.
- A cost block that compares the full and the reduced normal-matrix sizes.

## `demo_sensors`: system / sensor / day, gaps, and a bad day

**Purpose.** This example declares the full hierarchy:

- a system is one `task.addSystem`;
- a sensor is one stream of its system, `sys.addStream`;
- a day is an entry of the task's `Days`.

It shows a data gap, finds a bad day through one variance factor per
sensor-day, and then estimates the per-sensor noise levels without that day.
The order is the point: a variance component is an average over its stream,
so an event inside the stream corrupts the level.

**Run.** `demo_sensors()` takes no arguments and returns nothing.

**The case.** There are six common unknowns and two systems:

- `GRACE-C` (`id 'C'`) carries `KBR`, `LRI` and `GPS` and has 4 locals;
- `GRACE-D` (`id 'D'`) carries `KBR` and `LRI` and has 3 locals.

The data cover 10 days with 60 observations per sensor and day. The true
sigmas are KBR 1.0, LRI 0.4 and GPS 2.5, but every sensor is declared with
sigma 1. GRACE-D's LRI is missing on days 4–5, and GRACE-C's KBR carries a
6-sigma offset on day 7.

**Logic.**

1. Generate the data into a `containers.Map` per system.
2. Declare the task. Each system's `Arcs` is a resolver from `makeSchedule` that
   returns no arcs on an outage day, which is how a data gap is declared. Every sensor
   is a stream with `Noise=lsqvce.Noise.white(1)`, the claim.
3. Solve once with `Options(Verbose=false)` and print `R.sigma0` (the
   symptom).
4. Find the bad day. Run
   `Options(Estimate="noise", Pool="obsType", ArcFactors=true, MaxPasses=20, Noise=struct('Grid', 8))`.
   Every arc, here one system, sensor and day, gets its own variance
   factor. The demo looks each factor up by `lsqvce.util.arcKey` in
   `R.noise.arcKeys` and prints it divided by the median of its stream,
   because the factors are absolute against the claimed sigma 1. An arc
   above 5 times its stream's median is marked `<<< EVENT`.
5. Remove each event's arc from the store, so the arc resolver returns no
   arc for it, as for an outage. Then run
   `Options(Estimate="noise", Pool="obsType", ArcFactors=false, Noise=struct('Grid', 8))`.
   The KBR streams of both systems share one spectral shape, and every
   stream keeps its own level. The demo prints `R.summary().sigma` per
   stream against its truth.

`Noise.Grid 8` is used because an arc holds only 60 samples. With the
default 60-bin law, the bins that the data barely determine decay toward
zero, and the loop needed 168 passes. The help block records the measured
values: the event is 15.9 times its stream's median, every other sensor-day
lies between 0.5 and 1.8, and the levels come back within 4.0 % of the
truth. Estimated with the event still in the data, C|KBR came back 46 %
high.

## `demo_grace` and `demo_grace_dataLoader`: the noise loop on multi-axis streams

**Purpose.** This is a complete `Estimate="noise"` run on a GRACE-FO-shaped
synthetic problem with three-axis position streams and a scalar
inter-satellite stream. Every stream starts from a deliberately wrong
a-priori model: white noise with sigma 1. The help text calls it the
copy-paste starting point. To use it on your data, replace the loader with
your own task and its `Build` functions.

**Run.**

```matlab
R = demo_grace();                                  % maxPasses 10, verbose true
R = demo_grace('maxPasses', 8, 'verbose', false);
[cfg, ~, task] = demo_grace_dataLoader();          % the task alone
```

**Outline.** `demo_grace` holds the run and the reading.
`demo_grace_dataLoader` holds the declaration, with the local functions:

- `graceTask` (the task), `demoCfg` (the planted truth) and `catalogueInfo`;
- one `Build` per stream, `buildPosC`, `buildPosD`, `buildRr`, `buildPosE`, all calling
  `graceDay`, which draws the day and places each satellite's locals with
  `ctx.cols(satelliteLocals(ctx, sat))`;
- `commonParamTable` and `localParamTable`, the demo's 15-column catalogues;
- the generators `A_design`, `B_design` and `colouredNoise`.

**The declaration.** There are two systems:

- **System 1** merges satellites C and D. It carries the streams `posC`,
  `posD` and `rr`. The `rr` design couples the arc states of both
  satellites, and one system has only one local block, so C and D are
  merged with a stacked local block `[y_C | y_D]`.
- **System 2** is satellite E with `posE` only.

Both systems share six common unknowns. The data cover four days with two
arcs a day and 300 epochs per arc (5 s). The planted sigmas are 0.8, 1.3,
1.6 and 1.1.

The demo keeps its own 15-column parameter catalogue: column 3 names each parameter
(`addCommon` / `addLocal` with `Scope="system"`), and the other columns become the
template's `Info`; a `Build` finds a satellite's locals through
`ctx.Local.Info{k}.satellite`. The run settings (epochs, step, arcs per day, seed) are
the systems' `Settings`, read as `ctx.Settings`. Every stream declares
`Noise=lsqvce.Noise.white(1)`, the position streams with `Axes=3`.

For the position streams, each arc carries the per-epoch 3×3 covariance as `Sigma=`,
and the toolbox combines it with the temporal law. Rows are in epoch-major order.

**Run options and output.** The run uses
`Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=o.maxPasses, Tolerance=1e-3, Noise=struct('Smooth', 1))`.
`Tolerance=1e-3` shows the settled values in more digits. With
`Noise.Smooth 1` the loop converges in 4 passes. The plain 60-bin law
needed 27 passes and settled on the same levels. The demo prints:

- the true sigmas;
- `R.history.sigma`, one row per iteration and one column per stream;
- `R.summary()`, one row per stream;
- `R.convergence()`, one row per iteration;
- `R.report()`;
- `R.converged` and `R.history.nIter`.

The estimated shape is available through `R.plot("psd")` and `R.noise.S`.
In `R.summary()` the `arOrder` column stays `NaN`, because the estimator
fits a spectrum and not an AR order.

## `demo_gravity`, `demo_gravity_dataLoader` and `demo_gravity_design`: a physically structured gravity recovery

**Purpose.** This example runs the same noise loop on a design whose
columns are real spherical-harmonic potential partials along a near-polar
orbit, so the normal matrix has order-block correlation, polar-gap
conditioning and degree attenuation. It also shows how to read a gravity
solution in units of its own formal error.

**Run.**

```matlab
R = demo_gravity();
[R, score] = demo_gravity('MaxPasses', 6, 'Verbose', true);
R = demo_gravity('Plots', true, 'ExportDir', '/tmp/figs');   % two figures, saved as PNG
```

The name-value options are:

- `MaxPasses` (default 10; the run converges in 4);
- `Verbose` (default `false`);
- `Plots` (default `false`);
- `ExportDir` (default `''`).

`score` contains:

- `.sigma0`;
- `.stream`, `.sigmaPlanted`, `.sigmaRecovered` and `.sigmaRatio`;
- `.zMax` and `.zMeanSq`;
- `.anomaly` (`.true`, `.est`, `.sigma`, `.z`);
- `.degree` and `.degRatio`.

An `ExportDir` that does not exist is an error. The help states a cost of about 1.6 s:
5 days, 16 arcs, 3 streams, 77 common and 96 local unknowns.

**Outline.** There are three files:

- `demo_gravity` runs the problem and reads it. Its local functions are
  `section`, `streamType`, `degreeVariances`, `plotDegreeVariances`,
  `plotCoefficientErrors` and `exportIf`.
- `demo_gravity_dataLoader()` returns `[cfg, data, task]`. Its local functions
  are `gravityTask`, `gravityCfg`, one `Build` per stream (`buildPotC`, `buildPotD`,
  `buildKbr`, all calling `gravityArcs`) and `arcOrbit`.
- `demo_gravity_design(nmax, orbit)` returns `[A, meta]`. `A` has
  `(nmax+1)^2 - 4` columns (degrees 2..nmax, C and S terms). `meta.n`,
  `meta.m` and `meta.cs` label the columns. `orbit` needs the fields `u`,
  `lam0`, `incl` and `r` (with `r > 1`); bad input raises
  an error that names the faulty field.

**The declaration.** The field goes to degree `nmax = 8`, which gives 77
coefficients. The planted field follows a closed formula plus one anomaly at
(5, 3, C). The commons are named by term, degree and order (`"C2_0"`, `"S2_1"`, …) and
carry `degree`, `order` and `term` as `Info`.

There are three systems:

| System | Stream | Locals per arc |
|---|---|---|
| `GRACE-C` | `potC` | 2 |
| `GRACE-D` | `potD`, trailing C along-track | 3 |
| `KBR` | `kbr`: the potential difference D − C, whose design couples both orbits | 1 |

The locals (`bias`, `drift`, `quadratic`) are `"arc"`-scoped. The arc schedule is
ragged: 3, 4, 2, 4 and 3 arcs per day. Arcs are missing by ID on days 3 and 5, and every
arc has its own length between 140 and 220 epochs. Because of the `"arc"` locals, the
schedule is declared on each system as a resolver,
`Arcs=@(day, ~) cfg.arcsOfDay{...}`.

The declared noise has the right AR shape and unit sigma (`lsqvce.Noise.arma` for the
potential streams, `lsqvce.Noise.white` for `kbr`). The planted sigmas are 1.0, 2.0 and
0.5.

`gravityArcs` adds one arc per scheduled id of `ctx.ArcIds`, with the locals' columns
from `ctx.cols(names, Arc=a)`.

**Run options and output.** The run uses
`Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=opt.MaxPasses, Noise=struct('Smooth', 1))`.
It converges in 4 passes. The per-arc locals absorb the zero-frequency bin
of the 60-bin law, so without the smoothing that bin settled slowly and the
loop needed 72 passes. The demo prints three sections:

1. **The problem:** sizes, schedule and local counts.
2. **The weights:** planted and recovered sigma per stream from
   `R.noise.model{k}.sigma`, with their ratio. Only the ratios matter,
   because a common factor goes into `sigma0`.
3. **The field:** `z = (R.x - xTrue) ./ sqrt(diag(R.sol.Cxx))` (`R.sol.Cxx` is scaled by σ₀² already),
   max `|z|`, mean `z^2`, the anomaly with its sigma, and the ratio of the
   error degree variance to its formal value per degree.

The help block records the numbers of one measured run.

## `demo_gravfire` and `template/gravfire_data.m`: a GRAVFIRE bundle

**Purpose.** This example maps a GRAVFIRE preprocessing bundle
(`commonData`, `myOptions`, `systemsData`) onto a task and runs it through
`lsqvce.fit`. The bundle contains:

- Stokes coefficients as the common block;
- per-day initial conditions, accelerometer parameters and empirical
  parameters as locals;
- kinematic positions of two satellites and the K-band range-rate as
  observation groups.

The design matrices are not part of the bundle. You supply them.

**Run.**

```matlab
R = demo_gravfire(bundleFile);                                 % Estimate "none"
R = demo_gravfire(bundleFile, Estimate="noise", Passes=3, Verbose=1);
```

- `bundleFile` defaults to `""`, which means `testdata.mat` in the folder
  above the toolbox root.
- `Estimate` is `"none"` (default, one solve) or `"noise"`.
- `Passes` is the `MaxPasses` of the noise run (default 3).
- `Verbose` is 0..2 (default 1).

The bundle is private and not shipped. If the file does not exist, the demo
prints what it would do and returns `[]`.

**Logic.** With a bundle, the demo:

1. adds `template/` to the path and loads the file;
2. calls `task = gravfire_data(..., 'DryRun', true)`;
3. runs `lsqvce.fit(task, Verbose=...)`, or for `"noise"`
   `lsqvce.fit(task, Estimate="noise", Pool="obsType", Noise=struct('Grid', 2160), MaxPasses=Passes, CacheDays=Inf, OnPass=watch)`,
   with the not-converged warning suppressed;
4. prints `R.report` when `Verbose > 0`.

`fit` accepts the options as name-value pairs here, without an
`lsqvce.Options` object. `DryRun` puts random design matrices in place of
the real design step, so the demo tests the pipeline, not gravity recovery.

**`gravfire_data(commonData, myOptions, systemsData, Name=Value)`.** The
name-value options are:

- `'DesignFcn'`: your function `[A, B] = f(o, day, ctx)`, with `o` the bundle of that
  (day, group) and `ctx` the stream's `lsqvce.DayContext`, returning
  `m × ctx.NumCommon` and `m × ctx.NumLocal` (every local column from `ctx.cols`) with
  rows in epoch-major order;
- `'DryRun'`: deterministic random `A` and `B`, seeded per day and
  observation type.

`[task, data] = gravfire_data(...)` returns the `lsqvce.TaskDefinition` and, second,
its compiled declaration (`task.declaration()`) for engine-tier work. The task has the
calendar of `myOptions.days`, one common per Stokes row named `type_degree_order`
(`c_2_0`, `s_2_1`, …) with the 14 GRAVFIRE columns as `Info`, and one system `"GRACE"`
(`Id=1`, the bundle options as `Settings`) with one day's local templates as
`Scope="day"` and three streams. It declares no arc schedule: the pieces each `Build`
adds are the arcs.

The comments are numbered as steps (the steps of `dev/GRAVFIRE_INTERFACE.md`):

1. Name the common rows uniquely and add them with their 14 columns as `Info`.
2. Take the calendar from `myOptions.days`.
3. Name the streams: `posA`/`posB` for `kinematicPosition`, `kbr`
   for `range-rate`. Each gets `lsqvce.Noise.white(sqrt(variance(g)))`.
4. Take one day's rows of the expanded local table as `"day"`-scoped
   templates, with times in seconds of the day.
5. (No arc schedule.)
6. Bind everything `Build` needs into the system's `Settings`.
7. Add the system, its locals and its streams.

The local function `buildGravfireDay` is every stream's `Build`. It:

- builds epoch-major `l` and `t` for the whole day;
- converts `iSpp` into a per-epoch 3×3 `Sigma`;
- calls the design step;
- adds one arc per piece with `out.addArc(..., B=, Sigma=, Excluded=)`, the flagged
  epochs kept on the grid with weight 0.

Regularisation comes from column 11 of the tables:
`{'Tikhonov', fromPass, sigma[, label]}`. On common rows, `priorsFromColumn11`
builds the groups that `task.addPrior` declares in the `L` form. On local rows,
`gravfireLocalPriors` builds the groups that `sys.addPrior` declares in the `On`
(regular-expression) form. Labelled groups get `Estimate="scale"`; unlabelled rows get
`"none"`. The rows of one
label must share sigma and `fromPass`, and on local rows also the a-priori value
(column 10): the `On` form expands its rows in the toolbox's order, so one value per
member cannot be written there. The test suite runs the whole template on a generated
bundle (`dev/tests/fixture_gravfire_bundle.m`).

The function raises these errors:

| Identifier | Cause |
|---|---|
| `gravfire:designMissing` | no `DesignFcn` and no `DryRun`, raised when a day is read |
| `gravfire:unknownQuantity` | an observation group has an unsupported quantity |
| `gravfire:designShape` | `DesignFcn` returned `A` or `B` with the wrong size |
| `gravfire:badConstraint` | a malformed column-11 entry, or rows of one label that differ in sigma, `fromPass` or (local rows) the a-priori value |

The help also refers to `dev/GRAVFIRE_INTERFACE.md`, the step-by-step guide.

## `template/kbr_variational.m`: GRACE-FO range-rates by the variational equations

**Purpose.** This is a real-data template. It turns GRACE-FO Level-1B data
into a `lsqvce.TaskDefinition` for a gravity field in fully normalised spherical
harmonics:

- range-rates (KBR1B) of the C/D pair;
- reduced-dynamic positions (GNV1B);
- accelerometer (ACT1B) and attitude (SCA1B) data.

For each arc and satellite, the reference orbit and the state transition
`Psi` are integrated together with an 8th-order Adams–Bashforth–Moulton
integrator, started with RK4. The acceleration, its gradient and the
coefficient partials are analytic. The range-rate rows are formed on the
fly, so `Psi` is never stored.

**Run.**

```matlab
cfg = struct('days', 60000:60029, 'dataDir', '/data/gracefo/L1B', ...
             'gravityField', '/data/models/GOCO06s.gfc', ...
             'trf2crf', @myIers2010, 'degree', [2 60], 'arcsPerDay', 4);
[task, kit] = kbr_variational(cfg);
R  = lsqvce.fit(task, lsqvce.Options(Estimate="noise", Noise=struct('Grid', 1080)));
dC = kit.unitCommon .* R.x;                            % physical corrections
cfg.fieldUpdate = dC;  task = kbr_variational(cfg);    % second linearisation
```

**Configuration.** Four fields of `cfg` are required:

- `days` (MJD);
- `dataDir`;
- `gravityField`: an ICGEM `.gfc` file or a struct `.GM .R .C .S`;
- `trf2crf`: a function `[R, Rdot] = f(tGps)`.

The other fields default as follows:

| Field | Default |
|---|---|
| `degree` | `[2 60]` |
| `fieldDegree` | `degree(2)` |
| `arcsPerDay` | 1 |
| `arcs` | none; give a resolver to schedule arcs |
| `step` | 5 s |
| `posSampling` | 30 s |
| `sigma.kbr` | 2e-7 |
| `sigma.pos` | 0.02 |
| `obsTypes` | `{'kbr','posC','posD'}` |
| `files`, `formats` | GRACE-FO Level-1B v04 ASCII |
| `quatConvention` | `'srf2crf'` |
| `accToSrf` | `eye(3)` |
| `accBias0` | 0 |
| `accScale0` | 1 |
| `fieldUpdate` | none |
| `background` | no other forces |
| `maxAccGap` | 10 s |

**Outputs.** `[task, kit, data] = kbr_variational(cfg)`. The `task` has:

- the commons `C_n_m` and `S_n_m`, each with its unit (1e-10) as `Info`;
- the calendar `cfg.days`;
- one system `"GRACE-FO"` (`Id=1`) with:
  - `"arc"` locals: initial position and velocity per satellite (`C_x0`, …), and the KBR
    bias, drift and once-per-revolution cos/sin terms;
  - `"day"` locals: accelerometer bias and scale per satellite;
  - the arc schedule (`cfg.arcsPerDay` or `cfg.arcs`), the `Prepare` `prepareKbrDay`,
    and one stream per entry of `cfg.obsTypes`, each with the `Build` `buildKbrDay`.

Every parameter carries a declared unit as `Info`, so the physical value is
`unit .* estimate`. Position streams carry their sigmas in each arc's `Sigma`. `data` is
the compiled declaration (`task.declaration()`), for engine-tier work.

`kit` exposes the physics as function handles:

- `unitCommon` and `unitLocal`;
- `coeffIndex`, `harmonics`, `makeField` and `fieldEval`;
- `propagate` and `accModel`;
- `eraRotation`, `readGfc`, `readL1B`, `prepareDay` and `buildDay`.

`prepareKbrDay`, the system's `Prepare`, integrates every arc of the day once each time
the day is built and forms the rows of every stream, taking the columns of each arc's
locals from `ctx.cols(names, Arc=id)` (helper `arcColumns`); `buildKbrDay` picks its
stream out of `ctx.Prepared` and adds each arc with `Arc=id`. Only the current day's
integration is kept; `CacheDays` keeps the built rows. `kit.prepareDay` and
`kit.buildDay` take `(day, ctx)`, so only the toolbox calls them.

**Limitations.** The help lists these:

- `kit.eraRotation` is a stand-in for synthetic data only.
- L1B layouts and frame conventions must be checked against the product
  handbook.
- Forces other than the static field and the accelerometer enter only
  through `background`.
- KBR quality flags are not screened.
- Memory per stream is `8 m nc` bytes per day, about 0.5 GB for KBR at
  5 s and degree 60.

Errors carry `kbrvar:` identifiers, for example `kbrvar:badConfig`,
`kbrvar:noFile`, `kbrvar:accGap` and `kbrvar:emptyArc`.

## `template/lsqvce_template.m`: the engine-tier starting point

**Purpose.** This template declares a `lsqvce.TaskDefinition` from a `procOpt`
configuration struct and builds from it the cell array of `Loader` objects that the
engine-tier functions take (`lsqvce.engine.adjustBlock`, `lsqvce.engine.adjustRobust`).
Use it when your pipeline is driven by `procOpt` and you want the `Loader` objects in
hand. On the user path you do not need them: `lsqvce.fit(task)` runs the task.

**Run.** Copy the file, rename it, and write the stub `myBuild`:

```matlab
[L, task] = myCase(procOpt);
R = lsqvce.engine.adjustBlock(L, procOpt.days);
R.report();
```

**Structure.**

- **`declareTask`** is real code. It makes one `task.addSystem` per entry of
  `procOpt.observationGroups`, with its `"day"` / `"arc"` / `"system"` locals and one
  `addStream` per observation type, each with a white, unit a-priori `Noise`
  (`aprioriNoise`) and the `Build` `myBuild`. The tables of `procOpt` carry the scope in
  column 1 and a unique name in column 2.
- The Loaders come from `lsqvce.engine.declaration(task.declaration())`, the step `fit`
  uses.
- **`myBuild(day, ctx)`** is a stub. It raises `lsqvce_template:notImplemented`
  until you replace it. It must be deterministic and add the day's arcs with
  `out.addArc(A, l, t, ...)`, the per-epoch accuracy as `Sigma`, with multi-axis rows in
  epoch-major order.
- A `procOpt` without `observationGroups` or `days` raises
  `lsqvce_template:badProcOpt`.

## `demo_messaging`: what each verbosity level prints

**Purpose.** This example lets you compare the verbosity levels on real
output before you adopt a level for long runs. See also the chapter *Watching a run*.

**Run.**

```matlab
demo_messaging()                          % levels [0 1 2]
demo_messaging('levels', [0 1])
demo_messaging('toFile', '/tmp/levels.txt')
```

There are no outputs; everything is printed. A file that cannot be opened
raises `demo_messaging:badFile`.

**Logic.** For each level, `runAtLevel` redirects the toolbox log with
`lsqvce.util.log.set('fid', fid, 'screen', false)` and builds
`demo_grace_dataLoader()`. It runs
`Options(Estimate="noise", Pool="stream", ArcFactors=false, Verbose=lvl, MaxPasses=1)`,
then reads the captured text and prints it under a header. The header gives
the level and `lsqvce.util.verbosity.describe(lvl)`.

**What you will see.**

- **Level 0** prints nothing but warnings. Warnings are not controlled by
  the verbosity level, and a single pass does not converge.
- **Level 1** prints the stage label for each swept system and a
  convergence line for each pass.
- **Level 2** prints the same as level 1 when you go through `lsqvce.fit`.
  The help block explains why.

The per-stream view is not printed. Read it afterwards with `R.summary()`
or `R.convergence('Plot', true)`.

## `bench_qr_vs_normal`: which backend, measured on your machine

**Purpose.** This benchmark measures the cost of `Backend="qr"` against
`Backend="normal"` on your hardware, and shows when the accuracy of QR is
needed.

**Run.**

```matlab
bench_qr_vs_normal()                                       % u = [10 30 50 100], n = 1e5, nk = 5000, reps = 3
bench_qr_vs_normal('u', [10 50 100], 'n', 2e5, 'nk', 5000)
```

There are no outputs; the benchmark prints two tables.

**Logic.** Both runs use the same batched task,
`lsqvce.TaskDefinition.fromBatches(src, NumBatches=nb, NumUnknowns=u)`, with
`Options(Backend="normal")` and `Options(Backend="qr")`.

1. **Timing.** For each `u`, the benchmark times both runs, keeping the
   best of `reps` repetitions. It prints the times, their ratio, and the
   relative agreement of `R.x`. The timed call is the whole `fit`,
   including diagnostics, so the difference between the two times is the
   cost of the backend.
2. **Accuracy against conditioning.** With `u = 10` and `n = 20000`, the
   benchmark makes column 2 nearly collinear with column 1, for
   `eps = 1e-1 … 1e-8`. It compares both backends with an SVD reference
   and prints `cond(R.sol.Cxx)`, which equals `cond(N)`, the two errors,
   and a verdict. When `N` is singular in double precision, the normal
   path refuses with `lsqvce:normalEq:notSPD`. The benchmark catches that
   refusal and reports QR as mandatory.

The recipe the benchmark prints: solve once with `"normal"`, read
`cond(R.sol.Cxx)`, and switch to `"qr"` only when that number calls for it.
For the backend option itself, see the chapter *`lsqvce.Options`* of the API reference.
