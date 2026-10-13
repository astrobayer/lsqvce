#!/usr/bin/env python3
"""WP-132 mirror: lsqvce.kernel.stat.chi2Inv in the far lower tail (review low 30).

The old route seeded Newton (in y = log x, step damped to |dy| <= 2) from
Wilson-Hilferty, and from max(0.1 r, 1e-3) where WH is non-positive. Below
p ~ 1e-40 the root lies more than 200 units of y to the left, so the 100
damped steps ran out and the iteration stopped wherever it was: r = 1,
p = 1e-60 came back 9e31 times too large, r = 0.1 failed from p = 1e-6.

The new route (a line-for-line port of stat.chi2Inv below) adds:
  - a second seed, the root of the series' leading term
        F(x) <= (x/2)^a / Gamma(a+1),  a = r/2  (equality as x -> 0),
    taken when it lands closer than WH in log F;
  - below the median, Newton on log F (log F(e^y) is concave: the log-gamma
    density is log-concave), so the iterates approach the root monotonically;
  - a bracket [yLo, yHi] with bisection when a step leaves it;
  - a quantile below realmin returned as the leading term (subnormal or 0).

Reference: scipy.stats.chi2.ppf. Printed: worst relative error and most
iterations over r in [0.1, 1e6], p in [1e-307, 1 - 1e-6], for both routes,
and the largest change in the central range p in [1e-3, 0.995] (r >= 1).

Run:  python3 dev/validation/wp132_chi2inv_mirror.py
Created 2026-10-13 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-132.
"""
import warnings
import numpy as np
from scipy.special import gammainc, gammaln, ndtri
from scipy.stats import chi2

warnings.filterwarnings('ignore')            # log(0) and 0/0 are handled below, as in MATLAB
REALMIN = np.finfo(float).tiny


def dfdy(xc, a):
    """x * pdf(x), formed in logs."""
    return np.exp((a - 1) * np.log(xc) - xc / 2 - gammaln(a) - a * np.log(2) + np.log(xc))


def old_route(p, r):
    p, r = np.atleast_1d(p).astype(float), np.atleast_1d(r).astype(float)
    z = ndtri(p)
    x0 = r * (1 - 2 / (9 * r) + z * np.sqrt(2 / (9 * r))) ** 3
    bad = ~(x0 > 0) | ~np.isfinite(x0)
    x0[bad] = np.maximum(r[bad] * 0.1, 1e-3)
    y = np.log(x0)
    for it in range(100):
        xc = np.exp(y)
        step = (gammainc(r / 2, xc / 2) - p) / dfdy(xc, r / 2)
        step[~np.isfinite(step)] = 0
        step = np.clip(step, -2, 2)
        yn = y - step
        if np.all(np.abs(yn - y) <= 1e-15 * np.maximum(np.abs(y), 1)):
            return np.exp(yn), it + 1
        y = yn
    return np.exp(y), 100


def new_route(p, r):
    p, r = np.atleast_1d(p).astype(float), np.atleast_1d(r).astype(float)
    a = r / 2
    z = ndtri(p)
    x0 = r * (1 - 2 / (9 * r) + z * np.sqrt(2 / (9 * r))) ** 3
    y = np.log(2) + (np.log(p) + gammaln(a + 1)) / a          # the small-x root
    wh = (x0 > 0) & np.isfinite(x0)
    yW = np.log(x0[wh])
    mW = np.abs(np.log(gammainc(a[wh], np.exp(yW) / 2)) - np.log(p[wh]))
    mS = np.abs(np.log(gammainc(a[wh], np.exp(y[wh]) / 2)) - np.log(p[wh]))
    yWh = y[wh].copy()
    yWh[mW <= mS] = yW[mW <= mS]
    y[wh] = yWh
    tiny = np.log(2) + (np.log(p) + gammaln(a + 1)) / a < np.log(REALMIN)
    low = p < 0.5
    yLo, yHi = np.full(y.shape, -np.inf), np.full(y.shape, np.inf)
    it = 0
    for it in range(1, 101):
        xc = np.exp(y)
        Fc = gammainc(a, xc / 2)
        below = Fc < p
        yLo[below] = y[below]
        yHi[~below] = y[~below]
        d = dfdy(xc, a)
        step = (Fc - p) / d
        step[low] = (np.log(Fc[low]) - np.log(p[low])) * Fc[low] / d[low]
        step[low & (Fc == 0)] = -2
        step[~np.isfinite(step)] = 0
        step = np.clip(step, -2, 2)
        step[tiny] = 0
        yn = y - step
        out = ~tiny & np.isfinite(yLo) & np.isfinite(yHi) & ~((yn > yLo) & (yn < yHi))
        yn[out] = (yLo[out] + yHi[out]) / 2
        if np.all(np.abs(yn - y) <= 1e-15 * np.maximum(np.abs(y), 1)):
            y = yn
            break
        y = yn
    return np.exp(y), it


def err(x, ref):
    if ref < REALMIN:                        # subnormal / zero quantile: absolute
        return 0.0 if abs(x - ref) <= 1e-310 else 1.0
    return abs(x - ref) / ref


def main():
    rs = [0.1, 0.5, 1, 1.5, 2, 3, 5, 10, 30, 100, 300, 1e3, 1e4, 1e5, 1e6]
    ps = [1e-307, 1e-300, 1e-200, 1e-100, 1e-60, 1e-45, 1e-40, 1e-20, 1e-10, 1e-6,
          1e-3, 0.005, 0.025, 0.1, 0.5, 0.9, 0.975, 0.995, 0.999, 1 - 1e-6]
    for name, f in (('old', old_route), ('new', new_route)):
        worst, most, off = 0.0, 0, 0
        for r in rs:
            for p in ps:
                x, it = f(p, r)
                e = err(x[0], chi2.ppf(p, r))
                worst, most = max(worst, e), max(most, it)
                off += e > 1e-11
        print(f'{name}: worst relative error {worst:.2e}, most iterations {most}, {off} cases off by > 1e-11')
    x, _ = old_route(1e-60, 1)
    print(f'r = 1, p = 1e-60: old {x[0]:.3e}, reference {chi2.ppf(1e-60, 1):.3e}')
    # the central range the toolbox uses must not move
    rc = np.repeat([1, 2, 5, 10, 100, 1e3, 1e4, 1e6], 6)
    pc = np.tile([1e-3, 0.005, 0.025, 0.5, 0.975, 0.995], 8)
    d = np.max(np.abs(new_route(pc, rc)[0] / old_route(pc, rc)[0] - 1))
    print(f'central range: largest change {d:.1e}')
    x, it = new_route(1e-60, 1)
    assert abs(x[0] / chi2.ppf(1e-60, 1) - 1) < 1e-12
    assert d < 1e-13
    worst = max(err(new_route(p, r)[0][0], chi2.ppf(p, r)) for r in rs for p in ps)
    assert worst < 1e-11


if __name__ == '__main__':
    main()
