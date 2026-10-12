# Priors

A prior is a group of pseudo-observations `L·x = x₀ + v` with weight `1/σ²`, declared in
your columns:

```matlab
data.priors(j)            = struct(...)   % on the common parameters: L in the columns of A
data.systems(k).priors(l) = struct(...)   % on system k's local parameters: L in the columns of B
```

## Fields

| field | meaning |
|---|---|
| `group` | a regular expression on the name column of the table (`data.paramColumns.name`, default column 2). The toolbox builds L: one row per matched template row per instance (per day for a `'day'` row, per arc for an `'arc'` row, once for `'common'` / `'system'`). |
| `L` | the explicit operator (q × nc, or q × nl) in the columns `lsqvce.layout(data)` announces. Give `group` **or** `L`, never both. |
| `x0` | a-priori values (q × 1, or a scalar); default 0 |
| `sigma` | ONE standard deviation for all rows of the group; relative weights inside a group go into L |
| `from` | the pass at which the group becomes active (default 1). A single solve reports a group with `from > 1` (`lsqvce:fit:priorNotActive`) and skips it. |
| `name` | label (default: the pattern); a local group becomes the stream `prior_<name>` of its system |
| `estimate` | `"none"` (default): σ fixed; `"scale"`: σ is estimated in every pass of `Estimate = "noise"` |

Malformed fields, or an `L` of the wrong width, are refused by name
(`lsqvce:fit:badPrior`).

## Rules

- **Several rows make σ estimable.** A group with one row has redundancy ≤ 1; its σ is
  left unchanged and flagged in `R.noise.priorHistory`.
- **Local rows live in one day.** Every row of a local group has its non-zero columns
  inside one day's local block or the system block; days are coupled through common
  parameters only.
- **A datum is a group of its own.** Minimum constraints `L = Bᵀ` on the null space,
  `x0 = 0`, give the minimum-constraint solution and covariance for any σ. A group that
  only partly lies in the null space, or overlaps another group, is information.

## Where the result is reported

| | common groups | local groups |
|---|---|---|
| σ per pass | `R.noise.priorHistory` | `R.noise.model{k}.sigma` of stream `satKey|prior_<name>` |
| final | `R.noise.priors` | the same stream entry |

The rows of every prior enter `R.residuals()` and the redundancy count.

## From a parameter table

A row of a GRAVFIRE-style parameter table carrying `{'Tikhonov', fromPass, sigma[, label]}`
in its 11th column means exactly one such struct; `template/gravfire_data.m` emits them:
rows sharing a label become one `group` with `estimate = "scale"`, rows without a label
one fixed-σ group each.
