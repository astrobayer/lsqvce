# Covariance specifications

`data.systems(k).covariance.(obsType)` is one struct per observation type. Its `type`
selects the **temporal law** — how the errors of one stream correlate down the epochs —
and `sigma` its **level**. The per-epoch accuracy is not declared here: it travels with
the day as `Cep` (chapter `buildDay`). The toolbox composes the covariance of an arc
as C = D (R_t ⊗ I_d) D. Full reference: `doc lsqvce.engine.covModel`.

## Laws

| `type` | fields | describes | whitening cost (b rows per batch, m per arc) |
|---|---|---|---|
| `'diag'` | `sigma` (scalar, or a vector of per-row standard deviations) | white noise | O(b) |
| `'arma'` | `ar` = `[1 a₁ … a_p]`, `ma` (default 1), `sigma` (innovation std) | a known whitening filter | O(b·p) |
| `'psd'` | `f` (n_f×1, Hz, one-sided), `S` (n_f×1, unit²/Hz) | a measured noise spectrum | O(m³) per arc (Toeplitz) or O(b·p) (AR filter) |
| `'acf'` | `acf = @(dt, opt)`, `sigma` | correlation as a function of the lag | O(m³) per arc |
| `'var'` | `ar` (d×d×p), `sigma` (d×d) | vector AR: cross-axis temporal correlation (no `Cep`) | O(b·p·d²) |
| `'custom'` | `fcn = @(t, opt)` | any dense covariance | O(m³) per arc |

`sigma` is always a **standard deviation**, never a variance.

## Policy fields

| field | values | default | meaning under `Estimate = "noise"` |
|---|---|---|---|
| `estimate` | `"law"`, `"scale"`, `"none"` | `"law"` | `"law"`: spectral shape and level are estimated; `"scale"`: the declared shape is kept, only the level is estimated; `"none"`: the stream is used exactly as declared |
| `arcFactors` | logical | `Options.ArcFactors` | per-arc variance factors for this stream |

A per-row `sigma` vector makes the spec a fixed structure: it can carry a level but not
a spectral shape, so `estimate = "law"` on it is refused
(`lsqvce:fit:lawOnFixedStructure`). A `Cep` does not trigger this — per-epoch accuracy
together with an estimated temporal law is the supported combination.

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
  `psd`) are factorised over the whole arc (exact; O(m³), m² memory).
- An `arma` law whose MA polynomial has a root on or outside the unit circle is refused
  (`lsqvce:covModel:maNotInvertible`).
- Under `estimate = "law"` the declared law contributes only its `sigma` (the starting
  level): the shape is estimated from a white start.
