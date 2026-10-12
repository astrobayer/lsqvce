# Purpose

**lsqvce** is a MATLAB toolbox for least-squares adjustment with variance component
estimation. It was built for satellite gravity recovery (GRACE / GRACE-FO), where millions
of correlated observations from several instruments determine a few thousand common
parameters together with many arc-wise and day-wise local parameters, and it works just as
well on a levelling network or a curve fit.

## What it does

- **Adjusts** a linear or linearised Gauss–Markov model. The common parameters are shared
  by all data; local parameters belong to one arc, one day or one system and are
  eliminated on the fly, so the full design matrix is never formed.
- **Estimates the stochastic model** of every observation stream from the data: the
  spectral shape of its noise, its level, and a variance factor per arc
  (`Estimate = "noise"`).
- **Regularises by declaration.** A prior — a Kaula rule, minimum constraints, a
  background model — is a group of pseudo-observations in your columns; it carries a
  redundancy and can have its weight estimated.
- **Judges the result.** Global test, outlier tests, redundancy numbers, minimal
  detectable biases, whiteness of every stream, reduced to a verdict by `R.report`.
- **Is robust on request**: M-estimation with Huber, IGG3 or bisquare weights, and damped
  Gauss–Newton for nonlinear models, on the same loop.
- **Lets you watch**: a hook after every pass, and one plot per question afterwards.

## The three lines

```matlab
lsqvce_setup                                   % once
R = lsqvce.fit(data, Estimate="noise");        % the one verb
R.report                                       % the verdict
```

`data` is a plain struct that describes your problem; the options say how to run it; `R`
is everything that comes back.

## The three documents

| document | answers | |
|---|---|---|
| **this Manual** | *how do I use it?* — install, the workflow, declaring your problem, choosing the stochastic model, priors, reading the report, the examples | `doc lsqvce`, `lsqvce_manual.pdf` |
| **the API reference** | *what exactly does each interface accept and return?* — the `data` struct, `buildDay`, covariance specifications, `Options`, `Results`, every identifier | `html/lsqvce_api.html`, `lsqvce_api.pdf` |
| **the theory** | *what is computed, and why?* — the mathematical framework | `lsqvce_theory.pdf` |

Every unit also has its own page: `doc lsqvce.<name>` or `help lsqvce.<name>`.

## Reading order

New users: chapter *Installation*, then *The workflow, step by step* with
`examples/demo_workflow.m` open beside it, then *The examples* to find the one closest to
your problem. Come back to *Declaring your problem* and *Choosing the stochastic model*
when you write your own `buildDay`, and to *Reading the results and the report* every time
you run.
