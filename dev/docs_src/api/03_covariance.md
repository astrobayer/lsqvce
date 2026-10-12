# Covariance specifications

`data.systems(k).covariance.(obsType)` is one struct per observation type. Its `type`
selects the **temporal law** — how the errors of one stream correlate down the epochs —
and `sigma` its **level**. The per-epoch accuracy is not declared here: it travels with
the day as `Cep` (chapter `buildDay`). The toolbox composes the covariance of an arc
as C = D (R_t ⊗ I_d) D. Full reference: `doc lsqvce.engine.covModel`.

## Laws

| `type` | fields | describes | whitening cost per batch of b rows |
|---|---|---|---|
| `'diag'` | `sigma` (scalar, or a vector of per-row standard deviations) | white noise | O(b) |
| `'arma'` | `ar` = `[1 a₁ … a_p]`, `ma` (default 1), `sigma` (innovation std) | a known whitening filter | O(b·p) |
| `'psd'` | `f` (n_f×1, Hz, one-sided), `S` (n_f×1, unit²/Hz) | a measured noise spectrum | O(b·p) (AR filter) |
| `'acf'` | `acf = @(dt, opt)`, `sigma` | correlation as a function of the lag | O(b³) |
| `'var'` | `ar` (d×d×p), `sigma` (d×d) | vector AR: cross-axis temporal correlation (no `Cep`) | O(b·p·d²) |
| `'custom'` | `fcn = @(t, opt)` | any dense covariance | O(b³) |

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

- A declared `psd` law is whitened by an AR(30) filter fitted to its autocovariance;
  `Noise.Whitener` selects the route of the law the joint estimator estimates
  (`"toeplitz"`, `"ar"`, `"arma"`). Check an AR fit with `lsqvce.engine.checkPSD` (the
  relative error of the reproduced autocovariance should be below 1e-3).
- An arc is whitened in batches of up to 20 000 rows. Filter laws overlap consecutive
  batches by their memory, so the result does not depend on the batching; dense laws
  (`acf`, `custom`, Toeplitz) are factorised per batch, and the correlation between
  batches of a longer arc is dropped.
- Under `estimate = "law"` the declared law contributes only its `sigma` (the starting
  level): the shape is estimated from a white start.
