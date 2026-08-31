# lsqvce v3.0 — least-squares adjustment with variance component estimation

Declare your problem as a plain `data` struct, `fit` it, and get a judged
result: estimates with formal errors, per-observation diagnostics, and a
five-level verdict that says whether to believe them.

**The user surface is four names** at `lsqvce.` — `fit`, `Options`, `Results`,
`report` — plus the `data` struct specified in `DATA_CONTRACT.md`. Everything
that does the work lives one level down (`lsqvce.engine.*` drivers,
accumulators, loaders, VCE, diagnostics; `lsqvce.kernel.*` transforms, filters,
statistics) and stays reachable for experts.

What it is built for: streaming generalised least squares at satellite-gravimetry
scale — block-arrow systems (common parameters shared across systems, arc-local
nuisance parameters eliminated per system), coloured noise whitened per arc from
a declared covariance model, Helmert variance component estimation with
Monte-Carlo redundancy, robust estimation and data snooping, and diagnostics
that refuse to invent numbers (a `sigma0` of `NaN` is a result, not a bug).
Base MATLAB only — no Statistics or Optimization Toolbox; nothing of size `n`
is ever formed. MIT licensed. Version and release notes: `RELEASE_NOTES.md`;
the development record (status, correctness history, design decisions,
measurements) lives in `dev/DEVELOPMENT.md`.

## Installation

```matlab
cd <toolbox root>
lsqvce_setup           % adds the root (lsqvce.* packages) + examples/,
                       % registers the Help browser docs -> doc lsqvce
```

The setup adds **only** the root and
`examples/` to the path; `dev/tests` and `template/` are not path entries —
tests manage their path themselves, and the template is meant to be copied
into your own project.

## Quick start

Three lines, whatever kind of data you have:

```matlab
opts = lsqvce.Options(Model="vce");        % what to do (defaults = the engines' own)
R    = lsqvce.fit(data, opts);             % the one verb
R.report                                   % the five-level verdict
```

`data` is a plain struct; `fit` tells the three kinds apart by the fields present:

```matlab
% (B) in memory, diagonal weights, no arc-wise unknowns
data = struct('A', A, 'l', l, 'p', p);                         R = lsqvce.fit(data);
R = lsqvce.fit(data, lsqvce.Options(Robust="irls"));           % with M-estimation

% (C) batched callback [A,l,p] = batchFcn(k); with data.x0 the nonlinear form
data = struct('batchFcn', @(k) src(k), 'numUnknowns', u, 'numBatches', nBatch);

% (A) streaming: days, arcs, arc-wise unknowns, coloured noise (GRACE-FO)
data.common.numCommon = nc;   data.days = 1:30;
data.systems = struct('buildDay', @buildDayC, 'observationTypes', {{'kbr'}}, ...
                      'covariance', cov, 'systemId', 'C', 'localsPerArc', 6);
R = lsqvce.fit(data, lsqvce.Options(Model="vce", GroupBy="satObs"));
```

**New here?** In this order:

1. `doc lsqvce` → **Getting started** (`html/lsqvce_getting_started.html`): which
   kind of data you have, `buildDay` in six rules, reading the answer.
2. **`DATA_CONTRACT.md`** — the binding specification of the `data` struct and of the
   one function you write, `buildDay`, with a two-satellite worked example.
3. **`examples/demo_fit.m`** — the user path end to end on all three data kinds,
   including the GRACE-like physical oracle as a `data` struct; copy it.
4. `lsqvce_theory.pdf` — why the toolbox computes what it computes, and where it
   refuses to (NaN is a result).

## Usage guidelines

**Choosing `Options.Model`.** `"gls"` is one generalised solve under the
declared covariance. `"vce"` iterates: re-estimate the temporal noise model per
stream and the Helmert variance components, reweight, repeat until white and
balanced. `"pooled"` estimates one spectrum over all streams per pass.
`"ellmer"` is the ten-step joint co-estimation of PSD and keyed per-arc factors
(ITSG scheme); for it, set `NumLags` at least as large as the longest arc.
On the in-memory kinds, `Robust="irls"` switches to M-estimation
(Huber / IGG3 / bisquare). Engine tier for experts: `html/lsqvce_choosing.html`
is the decision flowchart between the drivers `fit` dispatches to.

**Declaring the stochastic model.** Each observation stream carries a
covariance spec: `diag` (weights), `arma` / `ar` (a whitening filter),
`acf` / `psd` (an exact Toeplitz law), plus an optional fixed per-epoch scaling
(e.g. a kinematic-orbit epoch covariance). Batches are overlapped by the filter
memory automatically, so batching never changes the result. Details and the six
`buildDay` rules: `DATA_CONTRACT.md`.

**Options.** `lsqvce.Options` is a typed value class: name-value on
construction, `disp(opts)` marks which knobs deviate from the engines' own
defaults. The groups that reach the VCE loop (`Shape`, `Scale`, `Outer`,
`Joint`, `Axes`) are documented in `doc lsqvce.Options`; you rarely need more
than `Model`, `GroupBy`, `MaxOuter`, `Verbose`, `Seed`.

**Reading the answer.** `R.report` prints the verdict; `R.coefficients()`
returns estimates, formal errors and 95 % bounds as a table; `R.x`, `R.sigma0`,
`R.history`, `R.noise`, `R.meta` (full provenance: version, release, UTC time,
every seed, the resolved options) carry the rest. The five levels of the
verdict:

| Level | Question | Key output |
|---|---|---|
| **0 Integrity** | Does the arithmetic hold? | `sum(r_i) == n−u` |
| **1 Global** | Is the model consistent? | `s0²/σ0²` vs χ² |
| **2 Localise** | Where is it broken? | snooping, Bonferroni, scale cross-check |
| **3 Reliability** | *What can you not check?* | redundancy numbers, MDB |
| **4 Parameters** | What do you actually know? | correlations, conditioning |

Level 3 is the one people skip and the one that decides whether the result can
be trusted. A `sigma0` of `NaN` means the quadratic form cancelled in double
precision — the toolbox reports that instead of clamping
(`lsqvce_theory.pdf`, section 7).

**Combining several data sets.** Add one system per satellite / instrument via
`data.systems(k)`; common parameters are stated once, each system's arc-local
parameters are eliminated by a Schur complement, and the variance components
weight the systems against each other. Regularisation is always declared, never
silent (`Reg`, `RegLambda`; recorded in the diagnostics).

**Parallelism.** Days within a system (`DistributeDays=true`), arcs in the
redundancy pass, and systems parallelise; nothing inside an arc (filter state)
or inside the dense solve (BLAS) does. Serial and distributed accumulation are
pinned equal to 1e-11 by test.

**Documentation map.**

```matlab
doc lsqvce                    % landing page: Getting started, guides, API
doc lsqvce.Options            % every knob
open lsqvce_manual.pdf        % the master manual
open lsqvce_theory.pdf        % theory, practice, error handling
edit(fullfile('template','lsqvce_template.m'))   % copy to start your own case
```

Developers: the suite, the audit, the doc builders and the full development
record live under `dev/` (`dev/DEVELOPMENT.md`; `dev/tests/run_all_lsqvce.m`
manages its path temporarily and restores it on return).

## License

MIT — see `LICENSE`.
