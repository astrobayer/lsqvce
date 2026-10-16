# Large problems

A GRACE-FO month at degree 60 has millions of observations, about 3700 common and tens of
thousands of local parameters. This chapter lists the settings that decide whether such a
run fits in time and memory. The defaults are safe; change them when you have measured.

## Memory

The accumulator holds the reduced common system, nc × nc. Local parameters are
eliminated arc by arc and never enter it. The largest transient object is one built day:
m × nc doubles per observation type — at degree 60 about 190 MB for a six-hour three-axis
arc.

`CacheDays` (default 1) is the number of built days kept resident per system. An
iterative run reads every day twice per pass; `CacheDays = Inf` builds each day once per
run if the calendar fits in memory. A system's `prepareDay` result is kept for the
current day only, whatever `CacheDays`: it adds one day's worth, and runs again whenever
a day is rebuilt. If reading the raw files is the expensive part, cache those inside your
`buildDay` closure instead.

## Time

| setting | effect |
|---|---|
| `Noise.Whitener = "ar"` | whitens a long arc in O(m·p) instead of the O(m³) of the dense Toeplitz factor, which is formed over the whole arc (m² memory, refused above `maxDenseArc` = 20 000 epochs); check the AR fit once with `lsqvce.engine.checkPSD`. It speeds up the SOLVE: under `Estimate = "noise"` the estimator's sweep still forms every arc's covariance densely (m² memory, O(m³) per arc and pass), and the same `maxDenseArc` refuses a longer arc there — cut long arcs, or estimate the law on a subset and solve the whole with `Estimate = "none"` |
| `Probes` | the cost of the redundancy sweep grows linearly with it; 100 is the standard, fewer gives noisier redundancies |
| `Parallel = struct('Mode', "parfor", 'DistributeDays', true)` | distributes the day loop of the accumulation over workers; the result is bitwise identical to the serial run |
| `MaxPasses`, `Tolerance` | the number of passes of an iterative run |

## Backend

`Backend = "qr"` (default) keeps the square-root information form and is the safer choice
on ill-conditioned common blocks; `"normal"` forms the normal equations and is faster on
well-conditioned ones. Both give the same numbers to rounding where both work.
`examples/bench_qr_vs_normal.m` measures the trade-off on your machine.

## Before a long run

1. Run one day with `Estimate = "none"` and read the report.
2. Run a few days with `Estimate = "noise"` and `MaxPasses = 3` and look at
   `R.plot("psd")` and `R.summary()`.
3. Fix `Seed` so the long run is reproducible, and save `R` with its `meta`.
