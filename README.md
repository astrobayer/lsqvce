# lsqvce v6.23.0 — least-squares adjustment with variance component estimation

Declare your problem as a plain `data` struct, `fit` it, and get a judged
result: estimates with formal errors, per-observation diagnostics, and a
five-level verdict that says whether to believe them.

**The user surface is five names** at `lsqvce.` — `fit`, `Options`, `Results`,
`report` and `layout` (the column layout to declare priors in) — plus the `data`
struct specified in the API reference (`doc lsqvce`). Everything
that does the work lives one level down (`lsqvce.engine.*` drivers,
accumulators, loaders, VCE, diagnostics; `lsqvce.kernel.*` transforms, filters,
statistics) and stays reachable for experts.

What it is built for: streaming generalised least squares at satellite-gravimetry
scale — block-arrow systems (common parameters shared across systems, arc-local
nuisance parameters eliminated per system), coloured noise whitened per arc from
a declared covariance model, Helmert variance component estimation with
Monte-Carlo redundancy, robust estimation and data snooping, and diagnostics
that refuse to invent numbers (a `sigma0` of `NaN` is a result, not a bug).
Base MATLAB only — no Statistics or Optimization Toolbox. The **solve** streams: the
normal equations are accumulated batch by batch and nothing of size `n` is formed
there. The **diagnostics** stage streams too: the whitened rows are re-read arc by arc,
not stored (6.2.0, WP-86 part 1), and the cofactor of the full parameter set (commons plus
all locals) is kept as blocks — the commons factor and each system's local blocks, about
1 GB for one GRAVFIRE-scale month (9 409 commons, 4 320 locals), where 6.2 built dense
13 729² matrices of 1.5 GB each (6.3.0, WP-86 part 2). The work per observation stays
O(nc²) because the commons dominate: about 1.7·10⁸ flops per row for the main diagnostic
pass at that scale, ~4·10¹⁴ for 2.07 M rows (1.7× less than the dense factor). The VCE
redundancy pass and `snoop`'s clean re-solve use the same blocks (6.4.0, WP-86 part 3).
MIT licensed. Version and release notes:
`RELEASE_NOTES.md`; the development record (status, correctness history, design
decisions, measurements) lives in `dev/DEVELOPMENT.md`.

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
opts = lsqvce.Options(Estimate="noise");   % what to do (defaults = the engines' own)
R    = lsqvce.fit(data, opts);             % the one verb
R    = lsqvce.fit(data, Estimate="noise"); % the same, Name=Value straight on fit
R.report                                   % the five-level verdict
```

`data` is a plain struct; `fit` tells the three kinds apart by the fields present:

```matlab
% (B) in memory, diagonal weights p = 1./sigma.^2, no arc-wise unknowns
data = struct('A', A, 'l', l, 'p', p);                          R = lsqvce.fit(data);
%     (a full SPD weight matrix is kind C's: return it from batchFcn -- R6 WP-101)
R = lsqvce.fit(data, Robust="huber");                           % M-estimation (IRLS)
data.priors = struct('L', B.', 'sigma', 1, 'name', 'datum');    % a free network's datum

% (C) batched callback [A,l,p] = batchFcn(k); with data.x0 the nonlinear form
data = struct('batchFcn', @(k) src(k), 'numUnknowns', u, 'numBatches', nBatch);

% (A) streaming: days, arcs, scoped parameter tables, coloured noise (GRACE-FO)
data.commonParams = [repmat({'common'}, nc, 1), yourColumns];   data.days = 59396:59425;
data.systems = struct('buildDay', @buildDayC, 'covariance', cov, 'id', 'C', ...
                      'localParams', localTemplates);       % buildDayC returns one CELL
R = lsqvce.fit(data, Estimate="noise", Pool="obsType", OnPass=@myInspect);   % per arc (v6)
```

Every kind runs on the same accumulator and the same loop; `Robust`, the
nonlinear form and `OnPass` (your inspection hook after every pass) work on
all of them. The v5 → v6 migration is one change (`RELEASE_NOTES.md`).

**New here?** `doc lsqvce` opens the landing page (`html/index.html`). Three documents, each with one job:

1. **The Manual** (`html/lsqvce_manual.html`, `lsqvce_manual.pdf`) — how to use the toolbox:
   installation, the workflow step by step, declaring your problem, choosing the stochastic
   model, priors, reading the report, and every example explained.
2. **The API reference** (`html/lsqvce_api.html`, `lsqvce_api.pdf`) — what every interface
   accepts and returns: the `data` struct, `buildDay`, covariance specifications, `Options`,
   `Results`, every error identifier.
3. **The theory** (`lsqvce_theory.pdf`) — the mathematical framework.

Then run `examples/demo_workflow.m` section by section, and copy the example closest to your
problem.

## Usage guidelines

**Choosing what to estimate.** `Estimate="none"` (default) is one generalised
solve under the declared laws. `Estimate="noise"` is the joint (Ellmer)
variance-component estimator on every stream whose spec allows it: the
spectral shape per pool group (`Pool="all" | "obsType" | "stream"`), one level
per stream, and — with `ArcFactors=true` — a keyed factor per arc. The v4
words: `"vce"` = `Estimate="noise", Pool="stream", ArcFactors=false`;
`"pooled"` = `Pool="all", ArcFactors=false`; `"ellmer"` = `Estimate="noise"`.
Set `Noise.Grid` at least as large as the longest arc. `Robust="huber" |
"igg3" | "bisquare"` is IRLS on the same loop, on every kind.

**Declaring the stochastic model.** Each observation stream carries a
covariance spec — its temporal law: `diag` (white), `arma` / `ar` (a whitening
filter), `acf` / `psd` (an exact Toeplitz law) — plus the per-stream policy
`estimate = "none" | "scale" | "law"` and `arcFactors`
(`examples/demo_noisePolicy.m` runs all four against a planted truth). The per-epoch accuracy
of the observations (a sigma per epoch, or a 3 × 3 covariance per epoch for
positions) comes back from `buildDay` with the day as `D.Cep`; the toolbox
composes the two (v5.3). Since v6 the day comes back as ONE CELL PER ARC
(`A`, `B`, `l`, `t`, `Cep`) and the arc schedule is optional — the cells are the
arcs. Batches are
overlapped by the filter memory automatically, so batching never changes the
result. Priors on the common parameters — a datum, a Kaula rule, a
background model — are `data.priors`, written in the columns
`lsqvce.layout` announces (`examples/demo_priors.m`). Details and the six `buildDay` rules:
the API reference (`doc lsqvce`).

**Options.** `lsqvce.Options` is a typed value class: 14 knobs and five
groups (`Noise`, `RobustOpts`, `Diagnostics`, `Parallel`, `GaussNewton`);
name-value on construction or straight on `fit`; `disp(opts)` marks which
knobs deviate from the defaults. You rarely need more than `Estimate`,
`Pool`, `MaxPasses`, `Verbose`, `Seed`, `OnPass`.

**Watching a run.** `OnPass=@(ctx) ...` is called at the end of every pass
with the current solution, the noise state, the live accumulator, your
declaration and your `processingOptions`; plot or save what you like, return
`false` to stop. After the run, `R.plot("residuals" | "contribution" | "psd"
| "convergence")`.

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
silent, and never an option: `data.priors(k)` on the common parameters and
`data.systems(k).priors(l)` on the local ones, each a pseudo-observation group
that carries its own redundancy (the Manual, chapter Priors and regularisation).

**Parallelism.** Days within a system (`DistributeDays=true`), arcs in the
redundancy pass, and systems parallelise; nothing inside an arc (filter state)
or inside the dense solve (BLAS) does. Serial and distributed accumulation are
pinned equal to 1e-11 by test.

**Documentation map.**

```matlab
doc lsqvce                    % landing page: Getting started, guides, API
doc lsqvce.Options            % every knob
open lsqvce_manual.pdf        % the Manual
open lsqvce_api.pdf           % the API reference
open lsqvce_theory.pdf        % the mathematical framework
edit(fullfile('template','lsqvce_template.m'))   % copy to start your own case
```

Developers: the suite, the audit, the doc builders and the full development
record live under `dev/` (`dev/DEVELOPMENT.md`; `dev/tests/run_all_lsqvce.m`
manages its path temporarily and restores it on return).

## License

MIT — see `LICENSE`.
