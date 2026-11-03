# Tutorials

Three tutorials, one for each kind of task. Each is a file in `examples/` that you can run,
and each section below follows its file step by step. Pick the one that matches your data;
the matrix tutorial is the shortest way in.

| tutorial | file | kind of task | use it when |
|---|---|---|---|
| 1 — a levelling network | `demo_tutorialMatrix` | `TaskDefinition.fromMatrices(A, l, Weights=p)` | the design matrix fits in memory and the observations have no local parameters |
| 2 — a long record in batches; a nonlinear model | `demo_tutorialBatch` | `TaskDefinition.fromBatches(@batchFcn, ...)` | the design is too large for one matrix, the weights are correlated in blocks, or the model is nonlinear |
| 3 — two instruments on a platform | `demo_tutorialStreaming` | `TaskDefinition(Days=...)` and its add-calls | arc-wise data read day by day, parameters local to an arc or a day, noise to be estimated |

All three end the same way: `R = lsqvce.fit(task)` returns a `lsqvce.Results`, and
`R.report` says whether to believe it.

## Tutorial 1 — A levelling network: `demo_tutorialMatrix`

### The case

Six new points P1 … P6 and a benchmark BM of known height 100 m are joined by 15 levelled
height differences. Each line has a length L between 1 and 4 km, and its standard deviation
is 1 mm per √km:

> h_ij = H_j − H_i + e_ij,   σ_ij = 0.001 m · √L_ij

The unknowns are the heights H1 … H6: u = 6 parameters, n = 15 observations, redundancy 9.

### Step 1 — The design matrix, the observations, the weights

```matlab
A = zeros(n, u);
for i = 1:n
    if to(i)   > 0, A(i, to(i))   =  1; end        % + H_to
    if from(i) > 0, A(i, from(i)) = -1; end        % - H_from
end
l = h + Hbm * (from == 0) - Hbm * (to == 0);       % observed minus computed: BM moved to l
w = 1 ./ sigma.^2;                                 % WEIGHTS, not standard deviations
task = lsqvce.TaskDefinition.fromMatrices(A, l, Weights=w);
```

Each row of `A` is one observation, each column one unknown. The benchmark is not an
unknown: its known height moves to the observation side, so `l` is *observed minus
computed*. `Weights` holds **weights**, `1/σ²`, one per observation — not standard
deviations. With these weights the a-priori σ₀ is 1: the declared sigmas are absolute.

### Step 2 — Fit

```matlab
R = lsqvce.fit(task);
T = R.coefficients();
V = R.report();
```

`demo_tutorialMatrix` prints

```
point       truth   estimate    se [mm]
P1       101.2340   101.2348       0.93
P2       102.8710   102.8721       1.05
...
P6       101.9050   101.9066       0.92
sigma0 = 0.991, verdict PASS
```

`R.coefficients()` is a table of estimate, standard error and 95 % bounds. σ₀ = 0.991: the
scatter of the data matches the weights. A σ₀ far from 1 says the weights are wrong in
scale — fix that before anything else.

### Step 3 — How well is every observation controlled?

```matlab
Q = R.reliability();               % one row per observation: v, r, w, mdb, quadrant
```

The redundancy numbers r sum to n − u = 9. An observation with r near 0 is hardly checked
by the others: an error in it passes into the heights unseen. Here the weakest line is
line 8 (P6 to P2), r = 0.37; its minimal detectable bias is 8.0 σ = 8.3 mm. `mdb` is in
units of the observation's own sigma: multiply by its σ for millimetres.

### Step 4 — A gross error

The tutorial adds 20 mm to line 7 (P5 to P6) and fits again:

```matlab
Rp = lsqvce.fit(taskBad);                                    % the default test
baarda = lsqvce.Options(Diagnostics=struct('Test', "baarda"));
Rb = lsqvce.fit(taskBad, baarda);                            % Baarda's w-test
[~, worst] = max(abs(Rb.diag.w));                            % reject ONE: the worst
Rk = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(A(keep, :), lb(keep), Weights=w(keep)), baarda);
Rr = lsqvce.fit(taskBad, Robust="huber");                    % or down-weight instead
```

```
default (Pope) test: sigma0 3.87, flagged 7, verdict PASS WITH WARNINGS
Baarda w-test:       flagged [7 8 9 10], worst |w| 11.3 in line 7
line 7 removed:     sigma0 0.99, flagged zeros(1,0)
Robust="huber":      weight factor of line 7 = 0.127
```

- The report's default test (Pope's τ) scales the residuals by s0, the σ₀ estimated from
  the data — including the blunder. It finds line 7 here, but in a small network a blunder
  can inflate s0 enough to hide itself.
- Our weights are absolute, so Baarda's w-test, which scales by the a-priori σ₀ = 1, is the
  right test. It flags line 7 **and** lines 8, 9, 10: a blunder smears into the residuals of
  the lines that share its points.
- So reject **one** observation at a time, the one with the largest |w|, and fit again.
  Without line 7, σ₀ is back to 0.99 and nothing is flagged.
- `Robust="huber"` instead down-weights the blunder (weight factor 0.127) and keeps all
  observations.

`R.snoop()` automates the reject-and-refit loop. In a network this small it warns that it
cannot tell a blunder from a σ₀ declared too small (`lsqvce:snoop:indeterminate`); with
a few hundred observations or more it can.

### Step 5 — A free network and its datum

Drop the benchmark: BM becomes an unknown too, and the network can float up and down. The
normal equations are singular and `fit` says so by name (`lsqvce:qr:singularR`). A datum is
a declared prior, not an option:

```matlab
free = lsqvce.TaskDefinition.fromMatrices(Af, h, Weights=w);
free.addPrior("datum", L=ones(1, u + 1) / (u + 1), Value=100.9, Sigma=1);   % the mean height
Rf = lsqvce.fit(free);
```

The prior fixes the mean height at 100.9 m. Any σ gives the same solution: a prior on the
null space is a minimum constraint, it removes the defect without straining the data. The
height differences — the quantities the network does determine — agree with step 2 to
1.7e-13 m.

## Tutorial 2 — A long record in batches; a nonlinear model: `demo_tutorialBatch`

### The case

Ten years of hourly sea level at a tide gauge, 87 600 epochs:

> y(t) = c + d·t + a1·cos(ωt) + b1·sin(ωt) + a2·cos(2ωt) + b2·sin(2ωt),   ω = 2π / year

an offset, a trend and the annual and semi-annual tides: six unknowns. The noise is
correlated from hour to hour, AR(1) with φ = 0.6 and σ = 2 cm. You would not hold 87 600
rows in memory; in a real case they are spread over files. So the rows reach the toolbox one
month at a time, 120 batches, through a function you write.

### Step 1 — The batch function

```matlab
batchFcn = @(k) monthRows(cfg, k, "white");
[A1, l1, p1] = batchFcn(1);                % look at one batch first: 730 x 6, 730 x 1, 730 x 1
```

`[A, l, p] = batchFcn(k)` returns the rows of batch k: the design `A` (m × u), the
observations `l` (m × 1, observed minus computed) and the weights `p` (m × 1, `1/σ²`). The
toolbox calls it **twice per batch** — to solve, and again for the diagnostics — so it must
return the same rows every time. `monthRows` draws its noise from a stream seeded on k.

### Step 2 — Declare and fit

```matlab
task = lsqvce.TaskDefinition.fromBatches(batchFcn, NumBatches=120, NumUnknowns=6);
Rw = lsqvce.fit(task);
```

Nothing of size 87 600 × 6 is formed: the normal equations are accumulated batch by batch.

### Step 3 — The batches change nothing

The tutorial stacks the 120 batches into one matrix and fits that with `fromMatrices`. The
two solutions agree to 1.2e-15, relative. Cutting the rows into batches is bookkeeping.

### Step 4 — Correlated weights within a batch

The weight vector of step 2 says the hours are independent. They are not. A batch may
return its **full m × m weight matrix** instead, the inverse of its covariance — for AR(1)
noise a tridiagonal matrix:

```matlab
taskFull = lsqvce.TaskDefinition.fromBatches(@(k) monthRows(cfg, k, "ar1"), ...
                                             NumBatches=120, NumUnknowns=6);
R = lsqvce.fit(taskFull);
```

```
param           truth                white                AR(1) se ratio
offset        0.12000    0.12021 +- 0.00014    0.12021 +- 0.00027     2.00
trend         0.00320    0.00313 +- 0.00002    0.00313 +- 0.00005     2.00
...
sigma0                               0.999                0.999   verdicts PASS / PASS
```

The estimates hardly move, and σ₀ is 0.999 in **both** runs: σ₀ checks the level of the
noise, not its colour. What changes are the formal errors. They double, as theory says for
AR(1) noise and slowly varying parameters, √((1 + φ)/(1 − φ)) = 2. The white run puts the
trend 3.5 of its formal sigma from the truth; the AR(1) run puts it 1.4 sigma away. The
batches themselves are taken as independent; correlation across batches is the streaming
kind's job (tutorial 3).

### Step 5 — A nonlinear model

A GNSS station relaxes after an earthquake, u(t) = a·(1 − exp(−t/τ)) + c, observed daily
for two years (two batches of a year). The batch function now receives the current
parameters x and returns the Jacobian `J` at x and `w = l − f(x)`:

```matlab
task5 = lsqvce.TaskDefinition.fromBatches(@(k, x) relaxRows(nl, k, x), NumBatches=2, X0=nl.x0);
R5 = lsqvce.fit(task5, MaxPasses=20, Tolerance=1e-10);
```

`X0` is the start, here [0.03 0.2 0] against the truth [0.042 0.35 0.003]. `fit` runs damped
Gauss–Newton: it converges in 5 iterations to τ = 0.3416 ± 0.0059 (planted 0.35).
`R5.history` holds every iteration (`vPv`, `dxRel`, the damping `lambda`). A start too far
from the solution may converge to a different minimum: start near what you know.

## Tutorial 3 — A streaming problem: `demo_tutorialStreaming`

### The case

The file is `examples/demo_tutorialStreaming.m`. This is the shape of every real case the
toolbox is built for: days, arcs, parameters local to an arc, and noise whose colour is
unknown.


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

`demo_tutorialStreaming` prints:

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

- Tutorial 1 or 2 fits your case: replace the simulated rows with your own, keep the call.
- Tutorial 3 fits your case: replace the simulator in `prepareTutorialDay` and the two
  `Build` functions with your own reader, and the two `addCommon` names with your parameters. *The workflow, step by
  step* goes through the decisions in order.
- *Declaring your problem* covers several systems, day and system scopes, data gaps,
  multi-axis observations and per-epoch accuracies.
- *Choosing the stochastic model* says which noise law to start from and how to read
  the estimate.
- *Reading the results and the report* explains every line of the report.
