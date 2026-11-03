# Priors and regularisation

## The idea

A prior says "these parameters are observed to be x₀, with standard deviation σ". It is
written as pseudo-observations `L·x = x₀ + v` on named parameters, and the toolbox treats
it as data: it enters the normal equations, it has a redundancy (how much of it the data
overrule), and its σ can be estimated like the level of any other stream. Regularisation
is therefore never an option on the run: it is a declared part of the problem.

One call declares a prior group, on the object whose parameters it constrains:
`task.addPrior` for common parameters, `sys.addPrior` for the local parameters of that
system.

| setting | values | default | best practice |
|---|---|---|---|
| `name` | string, unique per owner | required | names the group in `R.noise` and the report |
| `On` | a regular expression on parameter names | — (give `On` or `L`) | the normal way; one row per matched parameter and instance (per day for a `"day"` local, per arc for an `"arc"` local) |
| `L` | q × nc (task) or q × nl (system), in the columns `task.layout()` announces | — | a datum (`L = Nᵀ` of the null space), degree-dependent weights, sums and differences |
| `Value` | scalar or q × 1 | 0 | the a-priori value (1 for scale factors) |
| `Sigma` | positive scalar | required | one standard deviation for the whole group |
| `Estimate` | `"none"`, `"scale"` | `"none"` | `"scale"` only for a group of many rows: a one-row group cannot estimate its σ |
| `From` | integer pass ≥ 1 | 1 | 2–3 lets the first passes see the data alone; a single solve skips such a group |

## Typical priors

**A datum for a free network.** Minimum constraints on the null space N of the design
(`L = Nᵀ`, value 0) fix the datum without distorting the solution. Declare the datum as a
group of its own:

```matlab
task.addPrior("datum", L=Nnull.', Sigma=1);
```

The solution and its covariance are the minimum-constraint ones for any σ, also when the
datum is declared as several groups (for instance translation and rotation separately).
The datum rows are constraints, not observations: the report does not count them as
uncontrolled. Without a datum a rank-deficient system is refused by name; the toolbox never
guesses one.

**A Kaula rule on a gravity field.** One group per rule, σ estimated:

```matlab
task.addPrior("kaula", On="^[CS]", Sigma=1e-10, Estimate="scale");
```

`On` is a regular expression on the names given to `addCommon`; the toolbox builds one
row per matched coefficient. Relative weights inside a group (degree-dependent Kaula
factors) go into an explicit `L` instead: divide each row by its own relative σ.

**Calibration parameters of an instrument.** Priors on local parameters sit on the system:

```matlab
grace.addPrior("accScale", On="^accScale", Value=1, Sigma=0.02, Estimate="scale");
```

The accelerometer scales of both satellites over a month form one group with many rows —
exactly what makes their σ estimable. A group with one row cannot estimate its own σ.

## Writing `L` yourself

`task.layout()` announces which column of A and of each system's local design holds which
parameter on which day and arc, with its name, so an operator can be built for any
combination a pattern cannot express:

```matlab
lay = task.layout();
c = find([lay.common.name] == "offset");             % the columns of a "day" common, one per day
L = [ones(1, numel(c)) / numel(c); diff(eye(numel(c)))];   % mean zero, smooth day to day
Lc = zeros(size(L, 1), lay.nc);  Lc(:, c) = L;
task.addPrior("offsetDatum", L=Lc, Sigma=1e-3);
```

`lay.common(j)` names the parameter of column j of A (`.name`, `.row`, `.day`);
`lay.systems(k).local(c)` does the same for column c of system k (`.name`, `.row`, `.day`,
`.arc`). `examples/demo_priors.m` builds priors both ways and shows that they agree.

## Reading the result

`R.noise.priors` and `R.noise.priorHistory` report the σ of every common group per pass;
a local group appears as the stream `prior_<name>` of its system in `R.noise.model`.
λ = σ₀²/σ² is the regularisation the data support. `R.report` counts the redundancy of
the prior rows like that of any observation.

A group with `Estimate = "scale"` is estimated only under `Estimate = "noise"` on the
streaming kind. In any other run it keeps its declared σ, and `fit` says so once
(`lsqvce:fit:priorScaleHeld`).

## Activating a prior later

`From = 3` activates a group in the third pass of an iterative run: the first passes then
see the data alone. A single solve skips such a group with a warning.
