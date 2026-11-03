# Introduction

**lsqvce** is a MATLAB toolbox for least-squares adjustment with variance component
estimation. You describe your measurements and your unknowns; the toolbox estimates the
unknowns, works out how noisy each kind of measurement really is, and tells you whether
the result can be trusted.

It was built for satellite gravity recovery (GRACE / GRACE-FO). There, millions of
correlated observations from several instruments determine a few thousand common
parameters, together with many nuisance parameters per arc or per day. It works just as
well on a small problem: a levelling network, a calibration, or a curve fit.

## What it does for you

- **Solves the adjustment.** A linear or linearised Gauss–Markov model. The *common*
  parameters are shared by all data. The *local* parameters belong to one arc, one day
  or one instrument. They are eliminated on the fly, so the full design matrix is never
  formed and memory stays small.
- **Learns the noise from the data.** With `Estimate = "noise"` the toolbox estimates,
  for every observation stream, the spectral shape of its noise, its level, and a
  variance factor per arc. You only have to give a starting guess.
- **Regularises by declaration.** A prior (a Kaula rule, minimum constraints, a
  background model) is a group of pseudo-observations in your own columns. It is treated
  like data and can have its weight estimated.
- **Judges the result.** Global test, outlier tests, redundancy numbers, minimal
  detectable biases and the whiteness of every stream are reduced to one verdict by
  `R.report`.
- **Is robust on request.** M-estimation with Huber, IGG3 or bisquare weights, in the
  same loop as the noise estimation. Nonlinear models are solved by damped Gauss–Newton
  (the batched kind with `x0`); robust weights and Gauss–Newton do not combine yet.
- **Lets you watch.** A hook runs after every pass, and one plot answers each question
  afterwards.

## Three lines

```matlab
lsqvce_setup                                   % once, to install
R = lsqvce.fit(task, Estimate="noise");        % the one verb
R.report                                       % the verdict
```

`task` is a `lsqvce.TaskDefinition` that describes your problem, filled with flat `add`
calls. The options say how to run it. `R` holds everything that comes back. The toolbox
has these user-facing names:

| name | what it is |
|---|---|
| `lsqvce.TaskDefinition` | the declaration: parameters, systems, streams, priors (`addCommon`, `addSystem`, `addPrior`) |
| `lsqvce.SystemDefinition` | one system of a task, made by `addSystem`: `addLocal`, `addStream`, `addPrior` |
| `lsqvce.Noise` | the noise law of a stream: `white`, `arma`, `psd`, `acf`, `var`, `custom` |
| `lsqvce.DayContext`, `lsqvce.DayData` | what your daily functions receive and return |
| `lsqvce.fit` | runs the adjustment: `R = lsqvce.fit(task, opts)` |
| `lsqvce.Options` | every setting of a run, with safe defaults |
| `lsqvce.Results` | what `fit` returns: estimates, errors, diagnostics, plots |
| `lsqvce.report` | prints the verdict (`R.report` calls it for you) |
| `lsqvce.layout` | shows the column layout `fit` will use, to write priors in |
| `lsqvce.solveContributions` | re-solves a fit for any weighting of its days (bootstrap) |

## Is it for you?

| your situation | where to start |
|---|---|
| a design matrix `A`, observations `l` and weights that fit in memory | chapter *Tutorials*, tutorial 1 (the *matrix* kind) |
| a problem too large for memory, weights correlated in blocks, or a nonlinear model | chapter *Tutorials*, tutorial 2 (the *batched* kind) |
| data that arrive per day and per arc, from several instruments, with nuisance parameters and correlated noise | chapters *Key concepts* and *Tutorials*, tutorial 3 (the *streaming* kind) |

## The three documents

| document | answers | where |
|---|---|---|
| **this Manual** | *How do I use it?* Installation, concepts, a tutorial, the workflow, declaring your problem, the noise model, priors, reading the report, the examples | `doc lsqvce`, `lsqvce_manual.pdf` |
| **the API reference** | *What exactly does each interface accept and return?* The task definition, `Prepare` and `Build`, the noise laws, `Options`, `Results`, every error identifier | `html/lsqvce_api.html`, `lsqvce_api.pdf` |
| **the theory** | *What is computed, and why?* The mathematical framework | `lsqvce_theory.pdf` |

Every function also has its own page: `doc lsqvce.<name>` or `help lsqvce.<name>`.

## How to read this Manual

**If you are new**, read in this order:

1. *Installation* — five minutes.
2. *Key concepts* — the vocabulary every later chapter uses: common and local
   parameters, systems, days, arcs, streams, and the two functions you write
   (a stream's `Build` and, optionally, a system's `Prepare`).
3. *Tutorials* — one runnable example per kind of task (`examples/demo_tutorialMatrix.m`,
   `demo_tutorialBatch.m`, `demo_tutorialStreaming.m`). Run one, then change it.
4. *The workflow, step by step* — the same path for your own data, with every decision
   explained.

**When you write your own case**, keep *Declaring your problem* and *Choosing the
stochastic model* open. Come back to *Reading the results and the report* after every
run. *Troubleshooting* lists the common error messages and what they mean.

**For large campaigns**, read *Large problems* before the first long run.
