# Regularisation: priors and groups

## The rule

*Parameters and their columns are yours* — the tables, `buildDay`, plots, `toTable`.
*Regularisation is the toolbox's* — the mechanism, the redundancy, the activation pass,
the estimation of its strength. The interface between them is an operator in your
columns, announced by `lsqvce.layout`, checked by width.

## What a prior is

Tikhonov regularisation of a parameter is an additional observation equation,
"parameter j is observed to be x₀ with standard deviation σ":

    L · x = x₀ + v,   weight 1/σ²

where L is the coefficient matrix of the pseudo-observations (one row per
pseudo-observation) exactly as A is the coefficient matrix of the observations. The
toolbox therefore treats a prior as data: it enters the normal equations, it has a
redundancy (how much of it the data overrule), and estimating its σ is variance-component
estimation on that group (Koch–Kusche).

## The declaration — one mechanism for common and local parameters

```matlab
data.priors(j)            = struct(...)    % on the common parameters, L in the columns of A
data.systems(k).priors(l) = struct(...)    % on system k's local parameters, L in the columns of B
```

One struct per **group**, with the fields

| field | meaning |
|---|---|
| `group` | a regular expression matched against the *name column* of the table (column 2; `data.paramColumns.name` to say otherwise). The toolbox builds L: one row per matched template row per instance (per day for a `'day'` row, per arc for an `'arc'` row, once for `'system'`/`'common'`). |
| `L` | the explicit operator (q × nc or q × nl) in the columns that `lsqvce.layout(data)` announces. Give `group` **or** `L`, never both. |
| `x0` | the a-priori values (q × 1, or a scalar; default 0) |
| `sigma` | ONE standard deviation for the rows of the group. Relative weights inside a group go into L (divide a row by its own relative σ); the struct's σ is the common factor. |
| `from` | the pass at which the group becomes active (default 1). In a single solve a `from > 1` group is reported (`lsqvce:fit:priorNotActive`) and skipped. |
| `name` | the label (default the pattern); local groups become the stream `prior_<name>`. |
| `estimate` | `"none"` (default; σ fixed) or `"scale"`: σ is estimated every pass of `Estimate="noise"` by the Helmert rule σ² ← σ²·Ω/r — the fixed point Ω = r ("the prior's weighted residual sum equals its redundancy"), so λ = s₀²/σ² is the regularisation the data support. |

**Several rows make σ estimable.** A group with one row has redundancy ≤ 1 and cannot
estimate its own σ (it is left unchanged and flagged in `R.noise.priorHistory`). The
x/y/z accelerometer scales of both satellites over a month are one group with 180 rows.

**Local rows live in one day.** A row of a local group has its non-zero columns inside
one day's local block (or the system block); coupling two days is done through the
common parameters. Different days of the same parameter are different rows — that is how
the estimate pools information across the calendar.

## How they are treated

Common groups are stacked into one prior on the common block, folded into the solve and
appended to the diagnostic sweep (they carry a redundancy). Local groups ride as their own
observation stream `prior_<name>` of the system — pseudo-observation arcs with A = 0,
B = L/σ, l = x₀/σ and a unit law — so accumulation, elimination, the sweeps, `residuals()`
and the redundancy count see them as data, and under `Estimate="noise"` the group's σ is
the stream's level: identical behaviour to the common groups, no second estimator.
`R.noise.priors`, `R.noise.priorHistory` (common) and `R.noise.model{k}.sigma` for stream
`satKey|prior_<name>` (local) report the result.

## The datum of a free network

Minimum constraints Bᵀx = 0 are the group `L = B'`, `x0 = 0`: when B spans the null space
the solution is the exact minimum-constraint one for **any** weight (the data's range and
the null space are orthogonal). Without it the accumulator refuses a rank-deficient
system by name (`lsqvce:qr:singularR`); it never guesses.

## The layout

```matlab
lay = lsqvce.layout(data);
lay.nc, lay.common(j)                  % the template row (and day) of every column of A
lay.systems(k).nl, lay.systems(k).local(c)   % (row, day, arc) of every column of B
c = lay.systems(k).col(row, day, arc); c = lay.colCommon(row, day);   % lookups (0 = not instantiated)
```

## Your table's column 11

`{'Tikhonov', fromPass, sigma[, label]}` on a row means exactly one of the structs above;
`template/gravfire_data.m` (STEP 6/7) emits them: rows sharing a label become one group in
the `group` form with `estimate="scale"`, rows without a label one fixed-σ group each.
