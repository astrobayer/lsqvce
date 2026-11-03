# lsqvce v7.0.0 — least-squares adjustment with variance component estimation

Declare your problem as a `lsqvce.TaskDefinition`, `fit` it, and get a judged
result: estimates with formal errors, per-observation diagnostics, and a
five-level verdict that says whether to believe them.

**The user surface** at `lsqvce.` is the declaration — `TaskDefinition` (with the
`SystemDefinition` its `addSystem` returns), `Noise` (the temporal laws), and the
`DayContext` / `DayData` your `Prepare` and `Build` functions receive and fill — and the
run: `fit`, `Options`, `Results`, `report`, `layout` (the column layout to declare priors
in) and `solveContributions` (a fit's days re-weighted without re-reading them: a
bootstrap). The API reference (`doc lsqvce`) specifies every one. Everything
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

Declare the problem once, then run it:

```matlab
task = lsqvce.TaskDefinition(Days=59396:59425);          % the calendar (MJD)
task.addCommon(coeffNames);                               % parameters shared by all data
grace = task.addSystem("graceC", Prepare=@prepareOrbit, Arcs=4);   % one instrument
grace.addLocal(["x0" "y0" "z0" "vx0" "vy0" "vz0"], Scope="arc");   % eliminated per arc
grace.addLocal(["accScale" "accBias"], Scope="day");
grace.addStream("kbrr", Build=@buildKbrr, Noise=lsqvce.Noise.psd(f, S));
grace.addStream("pos",  Build=@buildPos,  Noise=lsqvce.Noise.white(1), Axes=3);
task.addPrior("kaula", On="^[CS]", Sigma=1e-10, Estimate="scale");

R = lsqvce.fit(task, Estimate="noise", Pool="obsType");  % the one verb
R.report                                                  % the five-level verdict
T = R.toTable(task);                                      % every parameter, estimate, std
```

The two functions you write have one home each and take `(day, ctx)`: the system's
`Prepare` does the work all streams of a day share (`ctx.Prepared` in every `Build`), and
each stream's `Build` returns its arcs of the day:

```matlab
function out = buildKbrr(day, ctx)
    orb = ctx.Prepared;                      % from prepareOrbit(day, ctx)
    out = ctx.newDay();
    for a = ctx.ArcIds
        [A, l, t, P] = kbrrArc(orb, a);      % your model: A m x nc, l = O - C, t in seconds
        out.addArc(A, l, t, Local={["x0" "y0" "z0" "vx0" "vy0" "vz0"], P}, Arc=a);
    end
end
```

Partials of local parameters are given by name; the toolbox places them in the right
columns. Data that fit in memory, or a design read in batches, take one line each:

```matlab
R = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(A, l, Weights=1./sigma.^2));
R = lsqvce.fit(lsqvce.TaskDefinition.fromBatches(@batchFcn, NumBatches=K, NumUnknowns=u));
R = lsqvce.fit(lsqvce.TaskDefinition.fromBatches(@jacFcn, NumBatches=K, X0=x0));  % nonlinear
```

Every kind runs on the same accumulator and the same loop; `Robust` and `OnPass`
(your inspection hook after every pass) work on all of them, and the nonlinear form
(`X0`) on the batched kind — but not together with `Robust` (`lsqvce:fit:notYet`). The
v6 → v7 migration table (`data` struct → `TaskDefinition`) is in `RELEASE_NOTES.md`.

**New here?** `doc lsqvce` opens the landing page (`html/index.html`). Three documents, each with one job:

1. **The Manual** (`html/lsqvce_manual.html`, `lsqvce_manual.pdf`) — how to use the toolbox:
   installation, the key concepts, a tutorial you can run (`examples/demo_tutorial.m`), the
   workflow step by step, declaring your problem (the task, its systems and streams, and
   where `Prepare` and `Build` go), choosing the stochastic model, priors, reading the report,
   every example explained, and troubleshooting.
2. **The API reference** (`html/lsqvce_api.html`, `lsqvce_api.pdf`) — what every interface
   accepts and returns: `TaskDefinition`, `Prepare` and `Build`, the `Noise` laws, `Options`,
   `Results`, every error identifier, and the compiled declaration the engines read.
3. **The theory** (`lsqvce_theory.pdf`) — the mathematical framework.

Then run `demo_tutorial` and read the Manual's chapter *Your first adjustment* beside it;
`examples/demo_workflow.m` walks the general path section by section. Copy the example
closest to your problem.

## Usage guidelines

**Choosing what to estimate.** `Estimate="none"` (default) is one generalised
solve under the declared laws. `Estimate="noise"` is the joint (Ellmer)
variance-component estimator on every stream whose spec allows it: the
spectral shape per pool group (`Pool="all" | "obsType" | "stream"`), one level
per stream, and — with `ArcFactors=true` — a keyed factor per arc. The v4
words: `"vce"` = `Estimate="noise", Pool="stream", ArcFactors=false`;
`"pooled"` = `Estimate="noise", Pool="all", ArcFactors=false`; `"ellmer"` =
`Estimate="noise"`. `Noise.Grid` is the number of spectral bins of the estimated law
(default 60; keep it with the default `Noise.Smooth = 1`). `Robust="huber" |
"igg3" | "bisquare"` is IRLS on the same loop, on every kind except the nonlinear form.

**Declaring the stochastic model.** Each stream carries its temporal law with its level,
`Noise=lsqvce.Noise.white(sigma)`, `.arma(ar, sigma)`, `.psd(f, S)`, `.acf(@rho, sigma)`,
`.var(ar, Sigma)` (axes correlated across time; not under `Estimate="noise"`) or
`.custom(@covFcn)`, and on `addStream` what the estimator may change about it:
`Estimate="law" | "scale" | "none"` and `ArcFactors` (`examples/demo_noisePolicy.m` runs
all four against a planted truth). The per-epoch accuracy of the observations (a sigma
per epoch, or a 3 × 3 covariance per epoch for positions) travels with each arc,
`out.addArc(..., Sigma=s)`; the toolbox composes the two. Batches are overlapped by the
filter memory automatically, so batching never changes the result. Priors — a datum, a
Kaula rule, a background model, calibration constraints — are `task.addPrior` on the
common parameters and `sys.addPrior` on a system's locals, by name pattern (`On=`) or as
an operator in the columns `task.layout()` announces (`examples/demo_priors.m`). Details
and the `Build` rules: the API reference (`doc lsqvce`).

**Options.** `lsqvce.Options` is a typed value class: 12 knobs and five
groups (`Noise`, `RobustOpts`, `Diagnostics`, `Parallel`, `GaussNewton`);
name-value on construction or straight on `fit`; `disp(opts)` marks which
knobs deviate from the defaults. You rarely need more than `Estimate`,
`Pool`, `MaxPasses`, `Verbose`, `Seed`, `OnPass`.

**Watching a run.** `OnPass=@(ctx) ...` is called at the end of every pass
with the current solution, the noise state, the live accumulator, your
compiled declaration and every system's `Settings`; plot or save what you like, return
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
(`lsqvce_theory.pdf`, section 1.3).

**Combining several data sets.** Add one system per satellite / instrument with
`task.addSystem`; common parameters are stated once, each system's arc-local
parameters are eliminated by a Schur complement, and the variance components
weight the systems against each other. Regularisation is always declared, never
silent, and never an option: `task.addPrior` on the common parameters and
`sys.addPrior` on the local ones, each a pseudo-observation group
that carries its own redundancy (the Manual, chapter Priors and regularisation).

**Parallelism.** Days within a system (`Parallel.DistributeDays = true`), arcs in the
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
record live under `dev/` (`dev/DEVELOPMENT.md`). The gate is
`cd('dev/tests'); ci_gate`, which manages its path temporarily and restores it on return.

## License

MIT — see `LICENSE`.
