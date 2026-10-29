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

data = struct('A', A, 'l', l, 'p', ones(size(t)) / sigma^2);   % p are WEIGHTS: 1/sigma^2
R = lsqvce.fit(data);
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

- `p` holds **weights**, `1/σ²`, not standard deviations.
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

### Step 1 — Declare the parameters and the calendar

```matlab
data = struct();
data.commonParams = {'common', 'cosAmp'; 'common', 'sinAmp'};
data.days = 1:5;
```

One row per parameter. Column 1 is the scope; column 2 is a name of your choice.

### Step 2 — Declare the system

```matlab
sys = struct();
sys.name        = 'platform';
sys.localParams = {'arc', 'bias'};             % one offset per arc
sys.arcs        = 2;                           % arcs 1 and 2 on every day
sys.processingOptions = struct('epochs', 240, 'step', 30, ...
                               'truth', [0.8; -0.3]);
sys.prepareDay  = @prepareTutorialDay;         % once per day
sys.buildDay    = @buildTutorialDay;           % once per day and observation type
sys.covariance.gauge = struct('type', 'diag', 'sigma', 0.05);
sys.covariance.radar = struct('type', 'diag', 'sigma', 0.02);   % colour unknown
data.systems = sys;
```

Line by line:

- `localParams = {'arc', 'bias'}` asks for one parameter per arc. The toolbox counts
  them (5 days × 2 arcs = 10) and tells `buildDay` where each one sits.
- `arcs = 2` is the arc schedule: arcs 1 and 2 on every day. A system without `'arc'`
  parameters may leave it out (its arcs are then simply the cells `buildDay` returns).
  With `'arc'` parameters it is required, because their number must be known before
  any data is read; leaving it out is an error that says so.
- `processingOptions` is a struct of your own. The toolbox passes it unchanged to your
  functions as `procOpt`. Here it carries the simulator's settings. In a real case it
  carries a data directory, a model file, a reference frame.
- The **field names of `covariance` are the observation types**. There is no separate
  list. Each gets a noise law; `'diag'` is white noise of level `sigma`.

### Step 3 — Write `prepareDay`: what both instruments share

Both streams need the epochs of every arc and the signal basis at those epochs. Computing
this in `buildDay` would do it twice per day. `prepareDay` does it once:

```matlab
function P = prepareTutorialDay(day, ~, ~, procOpt, info)
    nArc = numel(info.arcId);
    P = struct('t', {cell(1, nArc)}, 'G', {cell(1, nArc)}, 'y', {cell(1, nArc)});
    for a = 1:nArc
        start  = (day - 1) * 86400 + (info.arcId(a) - 1) * 43200;   % arcs at 0 h and 12 h
        t      = start + (0:procOpt.epochs - 1).' * procOpt.step;
        w      = 2 * pi / 3600;
        G      = [cos(w * t), sin(w * t)];
        rs     = RandStream('mt19937ar', 'Seed', 100 * day + info.arcId(a));
        bias   = 0.5 * randn(rs);                  % this arc's true offset
        P.t{a} = t;
        P.G{a} = G;
        P.y{a} = G * procOpt.truth + bias;
    end
end
```

It receives the day, the two parameter tables (unused here, hence `~`), your
`processingOptions` and `info`. `info.arcId` lists the arcs of this day. Whatever it
returns, `P`, is handed to every `buildDay` of the same day.

### Step 4 — Write `buildDay`: the rows of one instrument

```matlab
function D = buildTutorialDay(day, obsType, ~, ~, procOpt, info)
    P    = procOpt.prepared;                       % what prepareTutorialDay returned
    nArc = numel(info.arcId);
    D = struct('A', {cell(1, nArc)}, 'B', {cell(1, nArc)}, ...
               'l', {cell(1, nArc)}, 't', {cell(1, nArc)});
    for a = 1:nArc
        m  = numel(P.t{a});
        rs = RandStream('mt19937ar', 'Seed', 1000 * day + 10 * info.arcId(a) + double(obsType(1)));
        switch obsType
            case 'gauge'
                e = 0.05 * randn(rs, m, 1);                            % white
            case 'radar'
                phi = 0.5;                                             % AR(1)
                e = filter(1, [1 -phi], 0.02 * sqrt(1 - phi^2) * randn(rs, m + 200, 1));
                e = e(201:end);
        end
        B = zeros(m, info.nl);
        B(:, info.offsetOf(info.arcId(a)) + (1:info.count(a))) = 1;   % this arc's offset
        D.A{a} = P.G{a};                           % m x 2: the common parameters
        D.B{a} = B;                                % m x 10: the local parameters
        D.l{a} = P.y{a} + e;                       % observed minus computed
        D.t{a} = P.t{a};                           % epochs in seconds
    end
end
```

The four rules every `buildDay` follows:

1. **One cell per arc** in each of `A`, `B`, `l`, `t`, all with the same number of rows.
2. **Columns by `info`, never by hand.** `info.nl` is the total number of local
   parameters; `info.offsetOf(arcId)` is where this arc's block starts; `info.count(a)`
   is its width.
3. **Raw rows.** `l` is observed minus computed, unweighted and unwhitened. The toolbox
   applies the noise law itself.
4. **Deterministic.** The toolbox reads every day more than once. Here the random
   noise comes from a stream seeded on (day, arc, instrument), so every call returns the
   same rows. Your file reader is deterministic anyway.

Because `buildDay` uses `procOpt` (the 5th argument), it must declare at least five
arguments; the toolbox refuses a `prepareDay` whose result could never arrive.

### Step 5 — A first run under the declared noise

```matlab
R1 = lsqvce.fit(data);
R1.coefficients()
```

`Estimate = "none"` is the default: one solve with the weights you declared. Always run
this first. It shows that the declaration is accepted (every declaration error names the
field at fault) and that σ₀ is plausible. Here σ₀ = 0.993.

### Step 6 — Estimate the noise

```matlab
opts = lsqvce.Options(Estimate="noise", Seed=1);
R = lsqvce.fit(data, opts);
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
given automatically because we declared none (set `sys.id` to choose one). The radar's
spectrum is no longer flat: the estimator found the colour we left out.

## Where to go from here

- Replace the simulator in `prepareTutorialDay` and `buildTutorialDay` with your own
  reader, and the two `commonParams` rows with your parameters. *The workflow, step by
  step* goes through the decisions in order.
- *Declaring your problem* covers several systems, day and system scopes, data gaps,
  multi-axis observations and per-epoch accuracies.
- *Choosing the stochastic model* says which noise law to start from and how to read
  the estimate.
- *Reading the results and the report* explains every line of the report.
