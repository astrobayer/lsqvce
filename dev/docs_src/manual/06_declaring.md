# Declaring your problem

This chapter shows how to map a real problem onto a `lsqvce.TaskDefinition`. Every call
is listed with its settings, their values, their defaults and the practice that works;
the full reference of every class is in the API reference and in `help lsqvce.<name>`.

## The picture

You describe one adjustment with four kinds of call: *what is estimated for everybody*
(`task.addCommon`), *which systems deliver data* (`task.addSystem`, one per satellite,
receiver or instrument group), *what each system observes* (`sys.addStream`, one per
observation type) and *what each system must estimate besides* (`sys.addLocal`). Priors
are declared where their parameters live (`task.addPrior`, `sys.addPrior`).

You write two kinds of function, and each has one fixed home:

| function | given to | runs | returns |
|---|---|---|---|
| `Prepare(day, ctx)` | `task.addSystem(name, Prepare=@p)` | once per system and day, before its streams | anything; it reaches every `Build` of the day as `ctx.Prepared` |
| `Build(day, ctx)` | `sys.addStream(name, Build=@b, ...)` | once per stream and day | `ctx.newDay()` filled with `addArc` |

Everything else — reading the days in order, whitening every arc with its noise law,
eliminating the local parameters, iterating the variance components — is done by the
toolbox.

The design is block-arrow shaped: every arc sees the common parameters through its `A`
and its own local parameters through its local design. The local blocks are eliminated
arc by arc, so the memory need is set by the number of common parameters, not by the
number of observations.

## The task: `lsqvce.TaskDefinition`

```matlab
task = lsqvce.TaskDefinition(Days=59396:59426, Name="March 2021");
```

| setting | values | default | best practice |
|---|---|---|---|
| `Days` | integer row vector, each day once (e.g. MJD) | required unless every system sets its own `Days` | the full campaign here; narrow per system only where a satellite is missing |
| `Name` | string | `""` | a run label; it appears in `disp(task)` |

`disp(task)` prints what was declared: the parameter templates per scope, the systems and
their streams, the priors. Nothing is read before `lsqvce.fit(task)`.

## Parameters: `addCommon` and `addLocal`

Ask of every parameter: *over what span is it constant?* The answer is its scope.

| constant over … | call | `Scope` | examples | best practice |
|---|---|---|---|---|
| the whole adjustment, seen by all systems | `task.addCommon` | `"common"` (default) | gravity-field coefficients, station coordinates | keep nc as small as the science allows: memory and time grow with nc² |
| one calendar day, seen by all systems | `task.addCommon` | `"day"` | daily common terms | rare; every day must carry data |
| one arc of one system | `sys.addLocal` | `"arc"` (default) | initial position and velocity, per-arc biases | needs an arc schedule (`Arcs`) on the system or stream |
| one day of one system | `sys.addLocal` | `"day"` | accelerometer scale and bias | many days make a prior's σ estimable |
| the whole adjustment, one system | `sys.addLocal` | `"system"` | antenna offsets, instrument constants | not with the bootstrap over days (`solveContributions`) |

```matlab
task.addCommon(compose("C%d_%d", n, m), Info=table(n, m));   % names + your own columns
sys.addLocal(["x0" "y0" "z0" "vx0" "vy0" "vz0"], Scope="arc");
sys.addLocal(["accScale" "accBias"], Scope="day");
sys.addLocal("manoeuvre", Scope="day", Valid=[59410.25 59410.75]);   % one day only
```

| setting | values | default | best practice |
|---|---|---|---|
| `names` | string array, unique within the task (commons) or the system (locals) | required | readable, regex-friendly names; a shared prefix lets one prior cover a family |
| `Scope` | `"common"`, `"day"` (commons); `"arc"`, `"day"`, `"system"` (locals) | `"common"` / `"arc"` | see the table above |
| `Valid` | `[start end]` in fractional MJD, or `[]` | `[]` (always) | for one-off parameters on `"day"` / `"arc"` rows; such a row exists only on the days its interval overlaps |
| `Info` | a table with one row per name | none | extra columns your `Build` needs (degree, order, axis); each template's row reaches `Build` as a struct, `ctx.Common.Info{k}` / `ctx.Local.Info{k}` |

A name is a **template**: an `"arc"` local becomes one parameter per arc, a `"day"` local
one per day. The toolbox counts and places them; your code always refers to them by name.
A local may not share a name with a common parameter.

## Systems: `addSystem` — where the daily preparation goes

```matlab
grace = task.addSystem("GRACE-FO", Prepare=@prepareGrace, Arcs=4, ...
                       Settings=struct('dataDir', '/data/gracefo'));
```

| setting | values | default | best practice |
|---|---|---|---|
| `name` | string, unique in the task | required | a short label; it prefixes every stream in the report |
| `Prepare` | function handle `@(day, ctx)` | none | every computation two or more streams of a day need; deterministic |
| `Settings` | struct (paths, models, thresholds) | `struct()` | plain values only; one reader for several setups = one `Settings` per system, never a global |
| `Days` | integer row vector | the task's `Days` | narrow only for a system absent on some days |
| `Arcs` | a count, an id vector, or `@(day, stream) ids` | none: the arcs are those `Build` adds | required by `"arc"` locals; a resolver for ragged schedules and manoeuvres |
| `Id` | number or single character | its position | give every system an id or none |

A system is one satellite, one receiver, or a **pair whose observation depends on the
local parameters of both**: the range-rate between two GRACE spacecraft belongs to one
system that holds both satellites' parameters. Systems share only the common parameters;
data that share no common parameters are separate tasks.

### What `Prepare` is for

| question | answer |
|---|---|
| when does it pay? | when two or more streams of one system share an expensive step: an orbit integration, the variational partials, a large file, Earth-rotation matrices. With one stream it gains nothing |
| what does it get? | the day and a `lsqvce.DayContext`: `ctx.Settings`, `ctx.ArcIds` (the system's arcs of the day), the parameter tables |
| what may it return? | anything: a struct, a cell, an object. The toolbox does not look inside |
| how often does it run? | once each time the day is built: twice per pass with the default `CacheDays = 1`, once per run with `CacheDays = Inf`, and on every parallel worker that builds the day |
| how much memory? | the current day's result only; it is dropped when the next day is prepared |
| what must it be? | deterministic, like `Build` |
| what if it fails? | `lsqvce:makeLoader:prepareDay`, naming the day, with your error as the cause |

## Streams: `addStream` — where the daily build goes

```matlab
grace.addStream("kbr",  Build=@buildKbr,  Noise=lsqvce.Noise.white(1e-7));
grace.addStream("posC", Build=@buildPosC, Noise=lsqvce.Noise.white(1), Estimate="scale", Axes=3);
```

| setting | values | default | best practice |
|---|---|---|---|
| `name` | a MATLAB identifier, unique in the system | required | the observation type (`"kbr"`, `"posC"`, `"acc"`); it names the stream in the report |
| `Build` | function handle `@(day, ctx)` | required | one function per stream; take shared results from `ctx.Prepared`, never recompute them |
| `Noise` | a `lsqvce.Noise` | required | start from `lsqvce.Noise.white(sigma)` at a realistic level (*Choosing the stochastic model*) |
| `Estimate` | `"law"`, `"scale"`, `"none"` | `"law"` (`"scale"` for a per-row `sigma`) | `"law"` when the colour is unknown; `"scale"` for a trusted shape; `"none"` for reference streams |
| `ArcFactors` | `true`, `false` | `Options.ArcFactors` (`true`) | `false` for few or short arcs |
| `Axes` | integer ≥ 1 | 1 | 3 for positions; rows epoch-major (x₁ y₁ z₁ x₂ …) |
| `Arcs` | as `addSystem` | the system's | only when this stream has its own schedule |

### Writing `Build`

```matlab
function out = buildKbr(day, ctx)
    P = ctx.Prepared;                                  % prepareGrace's result for today
    out = ctx.newDay();
    for a = 1:numel(ctx.ArcIds)
        arc = P.kbr(a);                                % your arc, from the shared work
        m = numel(arc.t);
        out.addArc(arc.A, arc.observed - arc.computed, arc.t, ...      % A, l, t [s]
                   Local={stateNames, arc.Pstate; "kbrBias", ones(m, 1)});
    end
end
```

| `addArc` setting | values | default | best practice |
|---|---|---|---|
| `A` (1st) | m × nc, columns in `addCommon` order | required | the partials of the common parameters; `"day"` commons at `ctx.cols(names)` |
| `l` (2nd) | m × 1 | required | observed minus computed, raw: never whitened or weighted |
| `t` (3rd) | m × 1, seconds | required | the epochs; the whitener needs the lags |
| `Local` | n × 2 cell `{names, partials; ...}` | none | **the usual way**: name every local partial; the scope decides whether it is this arc's, this day's or the system's column |
| `B` | m × `ctx.NumLocal` (may be sparse) | none | when the partials already form one matrix; every column from `ctx.cols(names, Arc=id)`. One of `Local` or `B` per arc |
| `Sigma` | e × 1, e × d or e × d × d per epoch | none | what the preprocessing knows per epoch; then declare the law with σ = 1 |
| `Excluded` | m × 1 logical | none | keep rejected epochs on the grid with weight 0 (zero in `A`, local partials and `l`) |
| `Arc` | a scheduled arc id | the next | only to add arcs out of order |

When the local partials already exist as one matrix, build `B` from columns by name:

```matlab
B = sparse(m, ctx.NumLocal);
B(:, ctx.cols(stateNames, Arc=ctx.ArcIds(a))) = arc.Pstate;
B(:, ctx.cols("accScale")) = arc.Pacc;                 % a "day" local: no arc needed
out.addArc(arc.A, l, arc.t, B=B);
```

## Arcs and schedules

The arcs `Build` adds **are** the arcs of the day. Without `Arcs` the toolbox takes them
as they come, and days may differ in their number of arcs. Declare a schedule when `"arc"`
locals must be counted before any data is read:

```matlab
sys = task.addSystem("S", Arcs=4);                          % arcs 1..4 on every day
sys = task.addSystem("S", Arcs=@(day, stream) arcTable{day});   % a ragged schedule
```

With a schedule, `ctx.ArcIds` lists the day's arcs; `addArc` fills them in order, `Arc=id`
places one out of order, and a scheduled arc never added is a data gap. A day without any
data in a system belongs off that system's calendar (`Days`).

## Multi-axis observations and their accuracy

A stream with several axes per epoch — a kinematic orbit, a gradiometer — declares `Axes`
and returns its rows epoch-major (x₁ y₁ z₁ x₂ …). What the preprocessing knows about the
accuracy of every epoch travels with the arc as `Sigma`:

```matlab
out.addArc(A, l, t, Sigma=sigmaPerEpoch);                          % n x 1   scalar stream
out.addArc(A, l, t, Sigma=[sx sy sz]);                             % n x 3   per-axis sigmas
out.addArc(A, l, t, Sigma=permute(reshape(Q.', 3, 3, []), [3 1 2]));   % n x 3 x 3 covariance
```

The vector forms are standard deviations and the matrix form is a covariance — check one
epoch by hand before a campaign. `Sigma` is relative to the law's `sigma`: declare
`lsqvce.Noise.white(1)` when `Sigma` already carries the absolute accuracy.

## A complete declaration: a GRACE-FO pair

The range-rate between the two spacecraft depends on both orbits, so the pair is **one**
system with three streams: the kinematic positions of each satellite and the
inter-satellite range-rate. All three rest on the same orbit integration, which therefore
goes into the system's `Prepare`:

```matlab
task = lsqvce.TaskDefinition(Days=59396:59426);
task.addCommon(shNames(60));                                   % SH coefficients to degree 60
task.addPrior("kaula", On="^[CS]", Sigma=1e-10, Estimate="scale");

state = ["x0" "y0" "z0" "vx0" "vy0" "vz0"];
grace = task.addSystem("GRACE-FO", Id='G', Prepare=@prepareGracePair, ...   % both orbits, once a day
                       Arcs=@(day, ~) arcSchedule(day), Settings=struct('dataDir', '/data/gracefo'));
grace.addLocal([state + "C", state + "D"], Scope="arc");       % initial states of C and D
grace.addLocal(["kbrBias" "kbrCos1cpr"], Scope="arc");         % KBR empiricals
grace.addLocal(["accScale" "accBias" "accDrift"] + ["C"; "D"], Scope="day");   % ACC per day
grace.addStream("posC", Build=@buildPosC, Noise=lsqvce.Noise.white(0.02), Axes=3);
grace.addStream("posD", Build=@buildPosD, Noise=lsqvce.Noise.white(0.02), Axes=3);
grace.addStream("kbr",  Build=@buildKbr,  Noise=lsqvce.Noise.white(1e-7));  % colour: estimated

R = lsqvce.fit(task, Estimate="noise", Seed=1);
R.report
```

with the functions of the system:

```matlab
function P = prepareGracePair(day, ctx)
    P.C = integrateWithPartials('C', day, ctx.Settings);   % r, v and d[r;v]/dp
    P.D = integrateWithPartials('D', day, ctx.Settings);
end

function out = buildPosC(day, ctx), out = kinematicArcs(ctx.Prepared.C, "C", ctx); end
function out = buildPosD(day, ctx), out = kinematicArcs(ctx.Prepared.D, "D", ctx); end
function out = buildKbr(day, ctx),  out = rangeRateArcs(ctx.Prepared.C, ctx.Prepared.D, ctx); end
```

A runnable builder of this shape for real Level-1B data is `template/kbr_variational.m`
(chapter *The examples*).

## Checklist before the first run

- Every parameter has a name and the scope over which it is constant.
- Every stream has its own `Build`; what several streams share is in the system's
  `Prepare`.
- `Build` adds one arc per `addArc`, with `A` (m × nc, `addCommon` order), `l` observed
  minus computed, `t` in seconds; nothing is whitened.
- Local partials go by name (`Local=`), or as `B` with columns from `ctx.cols` — never by
  offset arithmetic.
- Calling `Build` twice for the same day gives identical output.
- `Sigma`, if given, matches the rows (epochs × axes = m).
- `disp(task)` shows what you meant before the first `fit`.
