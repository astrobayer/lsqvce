#!/usr/bin/env python3
"""WP-128 mirror: the critical value of the MAD-scaled outlier test
(diagnose / snoop with Test = "mad"; lsqvce.kernel.stat.madCrit).

The statistic is w_i = (v_i / sqrt(qv_i)) / s_MAD, s_MAD the centred MAD of
those raw scores. Pope's tau (the critical value the code used) holds only
when the scale CONTAINS the tested residual -- that is what bounds tau by
sqrt(r). A MAD scale does not, and carries its own sampling noise.

Part 1, null fits (no blunders), alpha = 1e-3 Bonferroni over n tests: the
family-wise false-alarm rate (any |w| above the threshold) per candidate.
Part 2, 10 % contamination at 25 sigma: the share of draws in which at least
all but one blunder are flagged, and Pope's (blind) count.

User decision 2026-10-13: t(r/2) -- the candidate at the nominal level for
every n. The price is power on short contaminated networks (n = 50); from
n ~ 200 every candidate finds them all, Pope none (test_lsqvce_ext>
test_pope_breaks_down_and_mad_does_not runs there).

Run:  python3 dev/validation/wp128_mad_crit_mirror.py   (~2 min)
Created 2026-10-13 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-128.
"""
import numpy as np
from scipy import stats

ALPHA = 1e-3


def tau_crit(a, r):
    t = stats.t.ppf(1 - a / 2, r - 1)
    return t * np.sqrt(r) / np.sqrt(r - 1 + t * t)


def crits(n, r):
    a = ALPHA / n
    return {'tau': tau_crit(a, r), 'normal': stats.norm.ppf(1 - a / 2),
            't(r)': stats.t.ppf(1 - a / 2, r), 't(r/2)': stats.t.ppf(1 - a / 2, r / 2),
            't(0.367r)': stats.t.ppf(1 - a / 2, 0.367 * r)}


def mad_w(A, l):
    x = np.linalg.lstsq(A, l, rcond=None)[0]
    v = A @ x - l
    qv = 1 - np.diag(A @ np.linalg.solve(A.T @ A, A.T))
    raw = v / np.sqrt(qv)
    s = 1.4826 * np.median(np.abs(raw - np.median(raw)))
    r = A.shape[0] - A.shape[1]
    pope = np.abs(v) / (np.sqrt(v @ v / r) * np.sqrt(qv))
    return np.abs(raw) / s, pope


def null_fwer(n, u, trials, seed):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, u))
    k = crits(n, n - u)
    hits = dict.fromkeys(k, 0)
    for _ in range(trials):
        w, _ = mad_w(A, rng.standard_normal(n))
        for name, kc in k.items():
            hits[name] += w.max() > kc
    return {name: h / trials for name, h in hits.items()}


def power(n, u, n_bad, draws):
    k = crits(n, n - u)
    found = dict.fromkeys(k, 0)
    pope_max = 0
    for sd in range(draws):
        rng = np.random.default_rng(1000 + sd)
        bad = rng.choice(n, n_bad, replace=False)
        A = rng.standard_normal((n, u))
        l = A @ rng.standard_normal(u) + rng.standard_normal(n)
        l[bad] += 25
        w, wp = mad_w(A, l)
        for name, kc in k.items():
            found[name] += np.sum(w[bad] > kc) >= n_bad - 1
        pope_max = max(pope_max, int(np.sum(wp[bad] > k['tau'])))
    return {name: f / draws for name, f in found.items()}, pope_max


def main():
    print('null family-wise false alarms (nominal 1e-3) | power: >= all-but-one of 10 % blunders')
    for n, u in ((30, 6), (50, 8), (200, 8)):
        f = null_fwer(n, u, 20000, 7)
        p, pope_max = power(n, u, n // 10, 300)
        print(f'n = {n:3d}: ' + '  '.join(f'{c} {f[c]:.4f}/{p[c]:.2f}' for c in f)
              + f'  | Pope flags at most {pope_max} blunders')
        assert f['tau'] > 4 * ALPHA or n >= 200
        assert ALPHA / 4 < f['t(r/2)'] < 2 * ALPHA
        if n == 200:
            assert p['t(r/2)'] == 1.0 and pope_max == 0


if __name__ == '__main__':
    main()
