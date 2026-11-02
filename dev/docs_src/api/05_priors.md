# Priors

A prior is a group of pseudo-observations `L·x = x₀ + v` with weight `1/σ²`, declared on
the object whose parameters it constrains:

```matlab
task.addPrior(name, Name=Value)   % on the common parameters: L in the columns of A
sys.addPrior(name, Name=Value)    % on system sys's local parameters: L in the columns of B
```

## Settings

| name | meaning |
|---|---|
| `name` | the group's name, unique per owner; a local group becomes the stream `prior_<name>` of its system |
| `On` | a regular expression on parameter names. The toolbox builds L: one row per matched name per instance (per day for a `"day"` template, per arc for an `"arc"` template, once for `"common"` / `"system"`) |
| `L` | the explicit operator (q × nc, or q × nl) in the columns `task.layout()` announces. Give `On` **or** `L`, never both (`lsqvce:TaskDefinition:priorTarget`); in-memory tasks take `L` only |
| `Value` | a-priori values (q × 1, or a scalar); default 0 |
| `Sigma` | ONE standard deviation for all rows of the group (required); relative weights inside a group go into L |
| `From` | the pass at which the group becomes active (default 1). A single solve reports a group with `From > 1` (`lsqvce:fit:priorNotActive`) and skips it |
| `Estimate` | `"none"` (default): σ fixed; `"scale"`: σ is estimated in every pass of `Estimate = "noise"` (in a run that cannot estimate it — `Estimate = "none"`, or the in-memory kinds — the group keeps its declared σ, with the warning `lsqvce:fit:priorScaleHeld`) |

A pattern that matches no name, or an `L` of the wrong width, is refused when `fit` starts
(`lsqvce:fit:badPrior`).

## Rules

- **Several rows make σ estimable.** A group with one row has redundancy ≤ 1: its σ is
  re-estimated every pass from that one residual (only a group with redundancy r ≤ 0 is
  left unchanged, flagged in `R.noise.priorHistory`). Give a `"scale"` group several rows,
  or declare it `"none"`.
- **Local rows live in one day.** Every row of a local group has its non-zero columns
  inside one day's local block or the system block; days are coupled through common
  parameters only.
- **A datum is a group of its own, or several.** Minimum constraints `L = Bᵀ` on the null
  space, value 0, give the minimum-constraint solution and covariance for any σ. The datum
  may be split over groups (translation and rotation declared separately): every group the
  data cannot see belongs to it, as long as it does not overlap a group the data does see.
  A group that only partly lies in the null space, or overlaps an information group, is
  information.
- **Datum rows are not observations.** Their redundancy is 0 by construction; the
  diagnostics keep them in Σr and vPv but do not test or count them (`D.nTests`,
  `D.datumRows`; quadrant `"datum"` in `R.reliability()`). A local group's row counts as a
  datum row when the data see less than a millionth of it (r ≤ 10⁻⁶).
- **Prior rows are declarations.** `R.snoop()` never rejects a prior row, common or local,
  and `R.systematic()` leaves the priors' rows out of its sequences (`D.priorRows`).

## Where the result is reported

| | common groups | local groups |
|---|---|---|
| σ per pass | `R.noise.priorHistory` | `R.history.sigma`, column `satKey|prior_<name>` |
| final σ | `R.noise.priors` (absolute) | `R.noise.model{k}.sigma` of stream `satKey|prior_<name>`: the factor **relative** to the group's declared σ (σ_final = declared σ × this factor) |

The rows of every prior enter `R.residuals()` and the redundancy count.

## From a parameter table

A row of a GRAVFIRE-style parameter table carrying `{'Tikhonov', fromPass, sigma[, label]}`
in its 11th column means exactly one prior group; `template/gravfire_data.m` emits them:
rows sharing a label become one group with `Estimate = "scale"`, rows without a label one
fixed-σ group each.
