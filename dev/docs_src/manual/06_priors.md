# Priors and regularisation

## The idea

A prior says "these parameters are observed to be x₀, with standard deviation σ". It is
written as pseudo-observations `L·x = x₀ + v` in your own columns, and the toolbox treats
it as data: it enters the normal equations, it has a redundancy (how much of it the data
overrule), and its σ can be estimated like the level of any other stream. Regularisation
is therefore never an option on the run: it is a declared part of the problem. The
fields are listed in the API reference (*Priors*).

## Typical priors

**A datum for a free network.** Minimum constraints on the null space N of the design
(`L = Nᵀ`, `x0 = 0`) fix the datum without distorting the solution. Declare the datum as
a group of its own:

```matlab
data.priors(1) = struct('L', Nnull.', 'x0', 0, 'sigma', 1, 'name', 'datum');
```

The solution and its covariance are the minimum-constraint ones for any σ. Without a datum
a rank-deficient system is refused by name; the toolbox never guesses one.

**A Kaula rule on a gravity field.** One group per rule, σ estimated:

```matlab
data.priors(1) = struct('group', '^[CS]', 'x0', 0, 'sigma', 1e-10, 'estimate', "scale");
```

`group` is a regular expression on the name column of `commonParams`; the toolbox builds
one row per matched coefficient. Relative weights inside a group (degree-dependent Kaula
factors) go into an explicit `L` instead: divide each row by its own relative σ.

**Calibration parameters of an instrument.** Priors on local parameters sit on the system:

```matlab
data.systems(k).priors(1) = struct('group', '^accScale', 'x0', 1, 'sigma', 0.02, ...
                                   'estimate', "scale");
```

The accelerometer scales of both satellites over a month form one group with many rows —
exactly what makes their σ estimable. A group with one row cannot estimate its own σ.

## Writing `L` yourself

`lsqvce.layout(data)` announces which column of A and B holds which template row on which
day and arc, so an operator can be built by name:

```matlab
lay = lsqvce.layout(data);
L = zeros(2, lay.nc);
L(1, lay.colCommon(5, [])) = 1;            % template row 5 of commonParams
L(2, lay.colCommon(6, [])) = 1;
```

`examples/demo_priors.m` builds priors both ways.

## Reading the result

`R.noise.priors` and `R.noise.priorHistory` report the σ of every common group per pass;
a local group appears as the stream `prior_<name>` of its system in `R.noise.model`.
λ = σ₀²/σ² is the regularisation the data support. `R.report` counts the redundancy of
the prior rows like that of any observation.

## Activating a prior later

`from = 3` activates a group in the third pass of an iterative run: the first passes then
see the data alone. A single solve skips such a group with a warning.
