# The workflow, step by step

This chapter walks once through a complete adjustment. Each step names the chapter that
treats it in depth. `examples/demo_workflow.m` is the same path as a script to run
section by section.

## Step 1 — Decide which kind of problem you have

| your problem | kind | you supply |
|---|---|---|
| the design matrix fits in memory, one block of unknowns | **matrices** | `data.A`, `data.l`, `data.p` |
| too large for memory, or nonlinear, one block of unknowns | **batched** | `data.batchFcn`, `data.numBatches`, `data.numUnknowns` (or `data.x0`) |
| data arrive per day and arc, several instruments, nuisance parameters per arc or day, correlated noise | **streaming** | `data.commonParams`, `data.days`, `data.systems` with your `buildDay` |

The first two are complete after one line:

```matlab
data = struct('A', A, 'l', l, 'p', 1 ./ sigma.^2);     % weights, not sigmas
R = lsqvce.fit(data);
R.report
```

The rest of this chapter is the streaming kind, the one the toolbox is built for.

## Step 2 — List the parameters

Write two tables. `data.commonParams` holds the parameters every observation can see (for
gravity recovery: the spherical-harmonic coefficients); `localParams` of a system holds
the nuisance parameters that live on one arc, one day or the whole system (initial
states, accelerometer calibration, empirical terms). Column 1 of every row says how often
it is instantiated; the remaining columns are yours:

```matlab
data.commonParams = [repmat({'common'}, nc, 1), myCoefficientNames];   % nc rows
s.localParams = [ {'arc', 'bias'}; {'arc', 'drift'}; {'day', 'accScale'} ];
```

The toolbox derives every count and every column position from these tables. See
*Declaring your problem*.

## Step 3 — Write `buildDay`

One function returns the data of one day and one observation type, one cell per arc:

```matlab
function D = myBuildDay(day, obsType, common, ~, procOpt)
    arcs = readMyArcs(procOpt.dataDir, day, obsType);   % your reader, deterministic
    n = numel(arcs);
    D = struct('A', {cell(1, n)}, 'B', {cell(1, n)}, 'l', {cell(1, n)}, 't', {cell(1, n)});
    for a = 1:n
        m = numel(arcs(a).omc);
        D.A{a} = myPartials(arcs(a), common);   % m x nc
        D.B{a} = zeros(m, 0);                   % no local parameters in this sketch
        D.l{a} = arcs(a).omc;                   % observed minus computed
        D.t{a} = arcs(a).tSeconds;              % epochs in seconds
    end
end
```

Return raw, unwhitened rows; never scale them by your weights. With local parameters, add
a 6th argument `info` and place every column by its offsets (*Declaring your problem*).

## Step 4 — Declare the system and its noise

```matlab
s = struct();
s.buildDay = @myBuildDay;
s.covariance.range = struct('type', 'diag', 'sigma', 2e-7);   % white noise, level 2e-7
s.processingOptions = struct('dataDir', '/data/mission');
data.days = 1:30;
data.systems = s;
```

The fieldnames of `covariance` are your observation types. Start with a `diag` law and a
realistic level; the estimator can find the colour (*Choosing the stochastic model*).

## Step 5 — First run: one solve under the declared model

```matlab
R = lsqvce.fit(data);          % Estimate = "none"
R.report
```

Check two things before anything else: that the run completes without a declaration error
(every error names the field at fault), and that `R.sigma0` is plausible. A σ₀ far from 1
means the declared accuracies are wrong in scale — fix that before estimating anything on
top of it.

## Step 6 — Estimate the stochastic model

```matlab
opts = lsqvce.Options(Estimate="noise", Pool="stream", Seed=1);
R = lsqvce.fit(data, opts);
R.report
R.plot("psd"); R.plot("convergence");
```

The loop re-solves until the noise laws settle. `Pool` decides which streams share one
spectral shape; `ArcFactors` (default on) gives every arc its own variance factor; `Seed`
makes the probe-based redundancies reproducible.

## Step 7 — Add regularisation if the problem needs it

A free network needs a datum; a high-degree gravity field needs a Kaula rule. Both are
declared as priors in your columns:

```matlab
data.priors(1) = struct('group', '^C|^S', 'x0', 0, 'sigma', 1e-10, 'estimate', "scale");
```

See *Priors and regularisation*.

## Step 8 — Read the report and act on it

`R.report` prints a verdict and the evidence behind it: the global test, flagged
observations, the reliability, the whiteness of every stream. Each adverse verdict has a
standard next step — snooping, a richer noise law, a functional model correction. See
*Reading the results and the report*.

## Step 9 — Take the results out

```matlab
T = R.coefficients();                  % estimates, formal errors, confidence bounds
[C, L] = R.toTable(data);              % your own tables with estimates written in
save('run.mat', 'R');                  % R.meta records version, seeds and options
```
