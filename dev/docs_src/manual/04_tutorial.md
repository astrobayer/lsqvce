# Your first adjustment

This chapter has two parts. Part 1 fits a straight line in three lines of code. Part 2
builds a small streaming problem from scratch: two instruments, ten arcs, an unknown
offset per arc, and noise whose colour is unknown. That is the shape of every real case
the toolbox is built for. Part 2 is the file `examples/demo_tutorial.m`; run it, then
read on.

## Part 1 — A straight line

Copy this into the Command Window:

```matlab
rng(1);
t     = (0:0.1:10).';
A     = [ones(size(t)), t];                 % offset and slope
xTrue = [2; 0.5];
sigma = 0.1;                                % the noise of every observation
l     = A*xTrue + sigma*randn(size(t));     % simulated measurements

task = lsqvce.TaskDefinition.fromMatrices(A, l, Weights=ones(size(t)) / sigma^2);   % WEIGHTS: 1/sigma^2
R = lsqvce.fit(task);
R.report
R.coefficients()
```

`R.report` ends in `VERDICT: PASS`, and `R.coefficients()` prints

```
    index    estimate       se        lower      upper
      1       1.9845     0.019632     1.9461      2.023
      2       0.5022     0.003392    0.49555    0.50885
```

the estimate of offset and slope, its standard error, and a 95 % confidence interval.
Both intervals contain the truth (2 and 0.5). `R.sigma0` is 0.994: the weights you gave
match the scatter of the data.

Two things to remember from part 1:

- `Weights` holds **weights**, `1/σ²`, not standard deviations.
- A σ₀ far from 1 means the weights are wrong in scale. Fix that before anything else.

If your problem fits in memory and the noise is known, this is all you need. The rest of
the chapter is for data that do not fit this mould.

## Part 2 — A streaming problem

### The case

A platform observes a periodic signal of one hour period,

> y(t) = a·cos(2πt / 3600) + b·sin(2πt / 3600),

with two instruments on one clock, a **gauge** and a **radar**. It records two arcs a day
of two hours each (240 epochs, 30 s apart) on five days. Every arc starts with an unknown
offset that both instruments see. The gauge has white noise with σ = 0.05. The radar has
correlated noise (AR(1), φ = 0.5) with σ = 0.02, but we pretend we do not know that and
declare it white.

| what | here | in general |
|---|---|---|
| common parameters | `a`, `b` | what you want to know |
| local parameters | one offset per arc (10 in all) | nuisance parameters, eliminated |
| system | the platform | one satellite, receiver, instrument group |
| observation types | `'gauge'`, `'radar'` | the kinds of measurement |
| what both streams share per day | the epochs and the signal basis | an orbit, a file, a rotation |

### Step 1 — Declare the task: its parameters and calendar

```matlab
task = lsqvce.TaskDefinition(Days=1:5);
names = ["cosAmp" "sinAmp"];
task.addCommon(names);                             % Scope="common" (default)
```

A task is one object you fill with `add` calls. `addCommon` declares the parameters every
system sees, by name; `Scope="common"` means one value for the whole run.

### Step 2 — Declare the system: its Prepare, locals and streams

```matlab
truth = [0.8; -0.3];                               % the simulator's truth
sys = task.addSystem("platform", Prepare=@prepareTutorialDay, Arcs=2, ...   % arcs 1 and 2
                     Settings=struct('epochs', 240, 'step', 30, 'truth', truth));
sys.addLocal("bias", Scope="arc");                 % one offset per arc
sys.addStream("gauge", Build=@buildGauge, Noise=lsqvce.Noise.white(0.05));
sys.addStream("radar", Build=@buildRadar, Noise=lsqvce.Noise.white(0.02));   % colour unknown
```

Line by line:

- `addSystem` declares the platform. Its `Prepare` function runs **once per day**, before
  the streams of that day; whatever it returns reaches every stream as `ctx.Prepared`.
- `Arcs=2` is the arc schedule: arcs 1 and 2 on every day. A system without `"arc"`
  parameters may leave it out (its arcs are then simply the arcs `Build` adds). With
  `"arc"` parameters it is required, because their number must be known before any data
  is read; leaving it out is an error that says so.
- `Settings` is a struct of your own. The toolbox passes it unchanged to your functions
  as `ctx.Settings`. Here it carries the simulator's settings. In a real case it carries a
  data directory, a model file, a reference frame.
- `addLocal("bias", Scope="arc")` asks for one parameter per arc. The toolbox counts them
  (5 days × 2 arcs = 10) and places them; your code names them.
- Each `addStream` is one observation type with its own `Build` function, which runs
  **once per day** of that stream, and its own noise law. `lsqvce.Noise.white(0.05)` is
  white noise of standard deviation 0.05.

### Step 3 — Write `Prepare`: what both instruments share

Both streams need the epochs of every arc and the signal basis at those epochs. Computing
this in each `Build` would do it twice per day. `Prepare` does it once:

```matlab
function P = prepareTutorialDay(day, ctx)
    cfg  = ctx.Settings;
    nArc = numel(ctx.ArcIds);
    P = struct('t', {cell(1, nArc)}, 'G', {cell(1, nArc)}, 'y', {cell(1, nArc)});
    for a = 1:nArc
        start  = (day - 1) * 86400 + (ctx.ArcIds(a) - 1) * 43200;  % arcs at 0 h and 12 h
        t      = start + (0:cfg.epochs - 1).' * cfg.step;
        w      = 2 * pi / 3600;
        G      = [cos(w * t), sin(w * t)];
        rs     = RandStream('mt19937ar', 'Seed', 100 * day + ctx.ArcIds(a));
        bias   = 0.5 * randn(rs);                      % this arc's true offset
        P.t{a} = t;
        P.G{a} = G;
        P.y{a} = G * cfg.truth + bias;
    end
end
```

It receives the day and a `lsqvce.DayContext`, `ctx`: here its `Settings` and `ArcIds`,
the arcs of this day. Whatever it returns, `P`, is handed to every `Build` of the same
day.

### Step 4 — Write `Build`: the rows of one instrument

```matlab
function out = buildGauge(day, ctx)
    out = instrumentDay(day, ctx, @(rs, m) 0.05 * randn(rs, m, 1));    % white
end

function out = buildRadar(day, ctx)
    phi = 0.5;                                                         % AR(1)
    out = instrumentDay(day, ctx, @(rs, m) arNoise(rs, m, phi, 0.02));
end

function out = instrumentDay(day, ctx, noise)
    P   = ctx.Prepared;                                % prepareTutorialDay's result
    out = ctx.newDay();
    c   = char(ctx.Stream);
    for a = 1:numel(ctx.ArcIds)
        m  = numel(P.t{a});
        rs = RandStream('mt19937ar', 'Seed', 1000 * day + 10 * ctx.ArcIds(a) + double(c(1)));
        out.addArc(P.G{a}, P.y{a} + noise(rs, m), P.t{a}, ...    % A, observed - computed, t [s]
                   Local={"bias", ones(m, 1)});                  % this arc's offset
    end
end
```

(`arNoise` draws AR(1) noise of standard deviation 0.02; it is in the example file.) The
four rules every `Build` follows:

1. **One `addArc` per arc**, with the design `A` of the common parameters (one column per
   `addCommon` name, in that order), `l` and the epochs `t`, all with the same number of
   rows.
2. **Local parameters by name, never by column.** `Local={"bias", ones(m, 1)}` gives the
   partial of the arc's `bias`; the toolbox knows which column that is on this day and arc.
3. **Raw rows.** `l` is observed minus computed, unweighted and unwhitened. The toolbox
   applies the noise law itself.
4. **Deterministic.** The toolbox reads every day more than once. Here the random noise
   comes from a stream seeded on (day, arc, instrument), so every call returns the same
   rows. Your file reader is deterministic anyway.

### Step 5 — A first run under the declared noise

```matlab
R1 = lsqvce.fit(task);
R1.coefficients()
```

`Estimate = "none"` is the default: one solve with the weights you declared. Always run
this first. It shows that the declaration is accepted (a declaration mistake is refused by
the call that makes it, with its name) and that σ₀ is plausible. Here σ₀ = 0.993.

### Step 6 — Estimate the noise

```matlab
opts = lsqvce.Options(Estimate="noise", Seed=1);
R = lsqvce.fit(task, opts);
R.report
```

The loop solves, looks at the residuals of every stream, updates each stream's noise law
and arc factors, and solves again. It stops when nothing changes by more than
`Tolerance`; here after four passes (`R.stopReason` is `"converged"`). `Seed` makes the
random probes that estimate the redundancies reproducible.

The report ends in `VERDICT: PASS`. Its Level 1 says *"The stochastic model settled in 4
pass(es)"* and *"Global test PASSES"*.

### Step 7 — Compare

`demo_tutorial` prints:

```
param       truth       declared noise      estimated noise
cosAmp     0.8000     0.8011 +- 0.0005     0.8012 +- 0.0008
sinAmp    -0.3000    -0.3006 +- 0.0005    -0.3004 +- 0.0008
sigma0                           0.993                1.000   (converged)
```

The estimates hardly move. What moves is **how well you know them**. In the first run,
σ₀ was fine, because the radar's *level* was declared right. But its *colour* was not:
correlated errors carry less information than white ones, so the first run's standard
errors are about 1.6 times too small. After the noise estimation they are honest. This
is the most common way an adjustment goes wrong without any error message, and the
reason `Estimate = "noise"` exists.

### Step 8 — Look

```matlab
R.summary()              % one row per stream: its final level relative to the declared one
R.plot("psd")            % the estimated noise spectrum of each stream
R.plot("residuals")      % the whitened residuals per stream
```

`R.summary()` lists the streams as `x1|gauge` and `x1|radar`: `x1` is the system's id,
given automatically because we declared none (`addSystem(..., Id='P')` chooses one). The radar's
spectrum is no longer flat: the estimator found the colour we left out.

## Where to go from here

- Replace the simulator in `prepareTutorialDay` and the two `Build` functions with your
  own reader, and the two `addCommon` names with your parameters. *The workflow, step by
  step* goes through the decisions in order.
- *Declaring your problem* covers several systems, day and system scopes, data gaps,
  multi-axis observations and per-epoch accuracies.
- *Choosing the stochastic model* says which noise law to start from and how to read
  the estimate.
- *Reading the results and the report* explains every line of the report.
