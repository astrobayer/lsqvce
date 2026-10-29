# Declaring your problem

This chapter shows how to map a real problem onto the `data` struct. The field-by-field
contract — every field, its type, its default and every error — is in the API reference
(chapters *The `data` struct* and *`buildDay`*).

## The picture

You describe one adjustment with three answers: *what is estimated for everybody*
(`data.commonParams`), *which days are processed* (`data.days`), and *which systems
deliver data* (`data.systems`, one entry per satellite, receiver or instrument). Each
system points to one function you write, `buildDay`, that returns the design rows and
observations of one day for one observation type. Everything else — reading the days in
order, whitening every arc with its noise law, eliminating the local parameters,
iterating the variance components — is done by the toolbox.

The design is block-arrow shaped: every arc sees the common parameters through its
`A` and its own local parameters through its `B`. The local blocks are eliminated arc by
arc, so the memory need is set by the number of common parameters, not by the number of
observations.

## Choosing the scope of each parameter

Ask of every parameter: *over what span is it constant?*

| constant over … | table | token | examples |
|---|---|---|---|
| the whole adjustment, seen by all systems | `commonParams` | `'common'` | gravity-field coefficients, station coordinates |
| one calendar day, seen by all systems | `commonParams` | `'day'` | daily common terms |
| one arc of one system | `localParams` | `'arc'` | initial position and velocity, per-arc biases |
| one day of one system | `localParams` | `'day'` | accelerometer scale and bias |
| the whole adjustment, one system | `localParams` | `'system'` | antenna offsets, instrument constants |

Rows are **templates**: the toolbox instantiates every row once per scope and derives the
column counts. With `data.paramColumns` a row can also carry a validity interval and is
then instantiated only on the days it covers — a one-day manoeuvre parameter is just a
row with a one-day validity.

## Writing `buildDay` with local parameters

Declare the 6th argument and never compute a column offset yourself:

```matlab
function D = buildDayGrace(day, obsType, common, ~, procOpt, info)
    raw  = readGraceDay(procOpt.dataDir, day, obsType);  % your reader, deterministic
    nArc = numel(info.arcId);                            % the schedule's arcs, in order
    D = struct('A', {cell(1, nArc)}, 'B', {cell(1, nArc)}, ...
               'l', {cell(1, nArc)}, 't', {cell(1, nArc)});
    for a = 1:nArc
        arc = raw.arc(a);                                % empty on a data gap
        D.A{a} = shPartials(arc, common);                % m_a x info.nc
        D.l{a} = arc.observed - arc.computed;            % observed minus computed
        D.t{a} = arc.tSeconds;
        B = sparse(numel(D.l{a}), info.nl);
        B(:, info.offsetOf(info.arcId(a)) + (1:info.count(a))) = arcPartials(arc);
        if info.dayCount > 0                             % this day's 'day' block
            B(:, info.dayOffset + (1:info.dayCount)) = accPartials(arc);
        end
        D.B{a} = B;
    end
end
```

`B` may be sparse; it usually is, because an arc touches only its own block.

## Arcs and schedules

The cells `buildDay` returns **are** the arcs of the day. Without an `arcs` field the
toolbox takes them as they come, and days may differ in their number of arcs. Declare a
schedule only when `'arc'` parameters need to be counted before any data is read:

```matlab
s.arcs = 4;                                  % arcs 1..4 on every day
s.arcs = @(day, obsType) arcTable{day};      % a resolver for ragged schedules
```

A data gap is an empty cell; a day without any data in a system belongs off that
system's calendar (`s.days`).

## Several systems

Each satellite or instrument is one entry of `data.systems`, with its own `buildDay`,
covariance specs, local parameters and calendar. Systems share only the common
parameters. An observation that depends on the local parameters of two satellites — the
range-rate between two GRACE spacecraft — belongs to **one** system whose local table
holds both satellites' parameters. One reader for several setups takes the setup as data:
plain values in `processingOptions`, its 5th argument, or a closure that binds a
configuration (API reference, chapter `buildDay`, *Where the handle goes*).

## Sharing the work of a day: `prepareDay`

`buildDay` is called once per day **and observation type**. When several streams of a
system rest on one expensive computation — the orbits of both satellites and their
partials for the whole day, used by the kinematic orbits *and* by the range-rates — that
computation would run once per stream. Give the system a `prepareDay` instead. It runs
**once per build of a day**, before the first stream of the day is built, and every
`buildDay` of that day finds its result in `procOpt.prepared`.

```matlab
s.prepareDay = @myPrepareDay;    % P = myPrepareDay(day, common, local, procOpt, info)
s.buildDay   = @myBuildDay;      % D = myBuildDay(day, obsType, common, local, procOpt, info)
```

| question | answer |
|---|---|
| when does it pay? | when two or more observation types of one system share an expensive step: an orbit integration, the variational partials, a large file, Earth-rotation matrices. With one observation type it gains nothing |
| what does it get? | as many of `(day, common, local, procOpt, info)` as it declares; `info` is the same day slice `buildDay` gets |
| what may it return? | anything: a struct, a cell, an object. The toolbox does not look inside |
| what must `buildDay` do? | declare at least five arguments, so that it receives `procOpt` (otherwise `lsqvce:makeLoader:prepareUnused`, reported by `fit` under its own prefix), and read `procOpt.prepared` |
| how often does it run? | once each time the day is built: twice per pass with the default `CacheDays = 1`, once per run with `CacheDays = Inf`, and on every parallel worker that builds the day |
| how much memory? | the current day's result only; it is dropped when the next day is prepared |
| what must it be? | deterministic, like `buildDay` |
| what if it fails? | `lsqvce:makeLoader:prepareDay`, naming the day, with your error as the cause |

The tutorial (*Your first adjustment*, `examples/demo_tutorial.m`) is the smallest
working example; the GRACE-FO pair below and `template/kbr_variational.m` are the real
case. The full contract is in the API reference, chapter `buildDay`, *What all streams of
a day share*.

## Multi-axis observations and their accuracy

A stream with several axes per epoch — a kinematic orbit, a gradiometer — is returned
epoch-major (x₁ y₁ z₁ x₂ …). What the preprocessing knows about the accuracy of every
epoch travels with the day as `D.Cep`, one cell per arc:

```matlab
D.Cep{a} = sigmaPerEpoch;                                % n x 1   scalar stream
D.Cep{a} = [sx sy sz];                                   % n x 3   per-axis sigmas
D.Cep{a} = permute(reshape(Q.', 3, 3, []), [3 1 2]);     % n x 3 x 3 full covariance
```

The vector forms are standard deviations and the matrix form is a covariance — check one
epoch by hand before a campaign. `Cep` is relative to the `sigma` of the covariance spec:
declare `sigma = 1` when `Cep` already carries the absolute accuracy.

## A complete declaration: a GRACE-FO pair

The range-rate between the two spacecraft depends on both orbits, so the pair is **one**
system with three observation types: the kinematic positions of each satellite and the
inter-satellite range-rate. All three rest on the same orbit integration, which therefore
goes into `prepareDay`:

```matlab
data.commonParams = [repmat({'common'}, nSh, 1), shParameterTable(60)];
data.days = 1:31;

state = {'x0'; 'y0'; 'z0'; 'vx0'; 'vy0'; 'vz0'};
s = struct();
s.id = 'G';  s.name = 'GRACE-FO';
s.prepareDay = @prepareGracePair;          % both orbits and their partials, once per day
s.buildDay   = @buildDayGracePair;         % picks the rows of one observation type
s.covariance.posC = struct('type', 'diag', 'sigma', 0.02);
s.covariance.posD = struct('type', 'diag', 'sigma', 0.02);
s.covariance.kbr  = struct('type', 'diag', 'sigma', 1e-7);   % colour: let the estimator find it
s.arcs = @(day, obsType) arcSchedule(day);                    % 'arc' rows need a schedule
s.localParams = [ ...
    [repmat({'arc'}, 6, 1), strcat(state, 'C')];                       % initial state C
    [repmat({'arc'}, 6, 1), strcat(state, 'D')];                       % initial state D
    {'arc', 'kbrBias'; 'arc', 'kbrCos1cpr'};                           % KBR empiricals
    [repmat({'day'}, 6, 1), {'accScaleC'; 'accBiasC'; 'accDriftC'; ...
                             'accScaleD'; 'accBiasD'; 'accDriftD'}]];  % ACC per day
s.processingOptions = struct('dataDir', '/data/gracefo');
data.systems = s;

R = lsqvce.fit(data, Estimate="noise");
R.report
```

with the two functions of the system:

```matlab
function P = prepareGracePair(day, ~, ~, procOpt, info)
    P.C = integrateWithPartials('C', day, procOpt, info);   % r, v and d[r;v]/dp
    P.D = integrateWithPartials('D', day, procOpt, info);
end

function D = buildDayGracePair(day, obsType, ~, ~, procOpt, info)
    P = procOpt.prepared;                  % prepareGracePair's result for this day
    switch obsType
        case 'posC', D = kinematicRows(P.C, day, info);
        case 'posD', D = kinematicRows(P.D, day, info);
        case 'kbr',  D = rangeRateRows(P.C, P.D, day, info);
    end
end
```

A runnable builder of this shape for real Level-1B data is `template/kbr_variational.m`
(chapter *The examples*).

## Checklist before the first run

- `buildDay` returns `A`, `B`, `l`, `t` as one cell per arc, each arc with one row count;
  `l` is observed minus computed; nothing is whitened.
- Calling `buildDay` twice with the same arguments gives identical output.
- `t` is in seconds.
- `size(A, 2)` equals the derived nc, `size(B, 2)` the derived nl — build them from
  `info.nc` / `info.nl` and place blocks by the `info` offsets.
- Every observation type is a field of `covariance`.
- `Cep`, if given, has one cell per arc and matches the rows (n·d = m).
- No unknown field names in `data`: a typo is an error on purpose.
