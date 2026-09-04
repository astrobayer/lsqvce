# Overview

**lsqvce** is a MATLAB toolbox for least-squares adjustment with variance component
estimation, built for satellite gravity recovery (GRACE / GRACE-FO) and equally at home
on a levelling network. Version 5 (September 2026) is the result of the R5 review track:
one accumulator, one stochastic-model estimator, one loop with three policies, one user
surface.

## What it does

* **Adjusts** a linear or linearised Gauss--Markov model with a common parameter block
  and per-system, per-day local blocks eliminated on the fly (the block-arrow structure of
  satellite gravimetry), on one accumulator with a QR or normal-equation backend.
* **Estimates the stochastic model** of every observation stream jointly (Ellmer 2018):
  a spectral shape per pool group, a level per stream, optional per-arc variance factors,
  by Hutchinson-probed redundancies -- `Estimate="noise"`.
* **Regularises by declaration**: priors on common and on local parameters are
  pseudo-observation groups in your columns, counted in the redundancy, with their
  sigma fixed or estimated by Helmert / Koch--Kusche.
* **Diagnoses**: Baarda / Pope data snooping, redundancy numbers, reliability, the
  whiteness of every stream, a five-level verdict (`R.report`).
* **Robust** M-estimation (IRLS) and damped Gauss--Newton on the same loop.
* **Watches**: `Options.OnPass` fires after every pass with the state and the live
  accumulator; `R.plot` afterwards.

## The three sentences you need

```matlab
lsqvce_setup                                  % once
R = lsqvce.fit(data, Estimate="noise", Pool="obsType");   % the one verb
R.report                                      % the verdict
```

`data` is a plain struct (chapter 3); everything about *how* to run is a
`lsqvce.Options` (chapter 4); everything that comes back is a `lsqvce.Results`
(chapter 5). The theory behind every number is `lsqvce_theory.pdf` (chapter 11).

## Reading order

New here: chapter 2 (install), then `demo_workflow` next to chapter 3, then chapter 4
when you want to change a knob. Migrating from 4.x: chapter 12 first.
