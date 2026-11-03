# Noise laws

`sys.addStream(name, Noise=n, Estimate=e, ArcFactors=f)` declares a stream's noise in two
parts: the **temporal law** with its **level**, a `lsqvce.Noise` object, and what the
estimator may change about it, the stream's `Estimate` and `ArcFactors`. The per-epoch
accuracy is not declared here: it travels with each arc as `Sigma` (chapter *Prepare and
Build*). The toolbox composes the covariance of an arc as C = D (R_t ⊗ I_d) D. Full
reference: `doc lsqvce.Noise`, `doc lsqvce.engine.covModel`.

## Laws

| constructor | arguments | describes | whitening cost (b rows per batch, m per arc) |
|---|---|---|---|
| `lsqvce.Noise.white(sigma)` | `sigma` scalar, or a column of per-row standard deviations | white noise | O(b) |
| `lsqvce.Noise.arma(ar, sigma, MA=ma)` | `ar = [1 a₁ … a_p]`, `MA` (default 1), `sigma` (innovation std) | a known whitening filter | O(b·p) |
| `lsqvce.Noise.psd(f, S)` | `f` (n_f×1, Hz, one-sided, **from 0**), `S` (n_f×1, unit²/Hz) | a measured noise spectrum | O(m³) per arc (Toeplitz) or O(b·p) (AR filter) |
| `lsqvce.Noise.acf(@rho, sigma)` | `rho(dt, opt)`, `sigma` | correlation as a function of the lag | O(m³) per arc |
| `lsqvce.Noise.var(ar, Sigma)` | `ar` (d×d×p), `Sigma` (d×d) | vector AR: cross-axis temporal correlation (no per-arc `Sigma`) | O(b·p·d²) |
| `lsqvce.Noise.custom(@covFcn)` | `covFcn(t, opt)` | any dense covariance | O(m³) per arc |

Every constructor takes `MaxDenseArc=m` (default 20 000 epochs), the longest arc
factorised densely. `sigma` is always a **standard deviation**, never a variance.

A `psd` grid starts at f = 0. The autocovariance integrates the piecewise-linear spectrum
over the grid only, so a grid that starts higher leaves out the band below its first point.
On 512 bins starting at the first nonzero bin, an AR(0.99) law lost 34 % of its variance
that way. A grid whose missing band, estimated as f(1)·S(1), exceeds 10⁻⁶ of the variance is
refused (`lsqvce:psd2acv:noZeroFrequency`). The toolbox does not guess the spectrum there:
prepend f = 0 with its value, for example `[0; f]` and `[S(1); S]` for a flat extension.

## What the estimator may change — on `addStream`

| setting | values | default | meaning under `Estimate = "noise"` |
|---|---|---|---|
| `Estimate` | `"law"`, `"scale"`, `"none"` | `"law"`; `"scale"` for a per-row `sigma` | `"law"`: spectral shape and level are estimated; `"scale"`: the declared shape is kept, only the level is estimated; `"none"`: the stream is used exactly as declared |
| `ArcFactors` | logical | `Options.ArcFactors` | per-arc variance factors for this stream |

A per-row `sigma` column makes the law a fixed structure: it can carry a level but not a
spectral shape: its default is `"scale"`, and `Estimate = "law"` on it is refused (`lsqvce:fit:lawOnFixedStructure`).
A per-arc `Sigma` does not trigger this — per-epoch accuracy together with an estimated
temporal law is the supported combination.

`Options.Pool` decides which streams share one estimated spectral shape; every stream
keeps its own level.

## How a law is whitened

- `Options.Noise.Whitener` selects how every `psd` law whitens — the estimated law and
  every declared one: `"toeplitz"` (default, exact), `"ar"` or `"arma"` (a fitted
  filter). Check an AR fit with `lsqvce.engine.checkPSD` (the relative error of the
  reproduced autocovariance should be below 1e-3).
- Filter laws are applied in batches of up to 20 000 rows that overlap by the filter
  memory (for an MA part, until the impulse response of 1/B(z) has decayed to rounding),
  so the result does not depend on the batching. Dense laws (`acf`, `custom`, Toeplitz
  `psd`) are factorised over the whole arc (exact; O(m³), m² memory: 3.2 GB at 20 000
  epochs). An arc longer than the law's `MaxDenseArc` is refused
  (`lsqvce:covModel:denseArcTooLong`): give a `psd` law a filter (`Whitener = "ar"`), cut
  the arcs shorter, or raise `MaxDenseArc` when the memory is there.
- An `arma` law whose MA polynomial has a root on or outside the unit circle is refused
  (`lsqvce:covModel:maNotInvertible`). Its exact start-up is a dense block as long as the
  impulse response of 1/B(z), so `MaxDenseArc` bounds that too: an MA root close to the
  unit circle on a long arc is refused (`lsqvce:covModel:denseArcTooLong`). The same limit
  bounds the whitener of an arc with many robust nuggets or excluded rows (m·k epoch
  products for k corrected epochs).
- `Whitener = "arma"` reduces the AR fit of the spectrum to an ARMA (Durbin's method,
  `lsqvce.kernel.ar2arma`). A reduction that comes back non-stationary or with a
  non-invertible MA part falls back to the AR fit, with the warning
  `lsqvce:ar2arma:unstable` or `lsqvce:ar2arma:nonInvertibleMA`.
- Under `Estimate = "law"` the declared law contributes only its `sigma` (the starting
  level): the shape is estimated from a white start.
