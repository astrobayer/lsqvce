#!/usr/bin/env python3
"""WP-128 mirror: the critical value of the MAD-scaled outlier test
(diagnose / snoop with Test = "mad"; lsqvce.kernel.stat.madCrit).

The statistic is w_i = (v_i / sqrt(qv_i)) / s_MAD, with s_MAD the centred MAD
of those raw scores. Pope's tau (the critical value the code used) holds only
when the scale CONTAINS the tested residual -- that is what bounds tau by
sqrt(r). A MAD scale does not, and its own sampling error is that of a
variance estimate with the MAD's asymptotic efficiency 0.367: the statistic
is Student t with nu = 0.367 r.

Null simulation (no blunders): least squares with n observations and u
parameters, alpha = 1e-3 Bonferroni over n tests; the family-wise false-alarm
rate (any |w| above the threshold) per critical value. Nominal: 1e-3.

Run:  python3 dev/validation/wp128_mad_crit_mirror.py   (~1 min)
Created 2026-10-13 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-128.
"""
import numpy as np
from scipy import stats


def tau_crit(a, r):
    t = stats.t.ppf(1 - a / 2, r - 1)
    return t * np.sqrt(r) / np.sqrt(r - 1 + t * t)


def main(trials=40000):
    rng = np.random.default_rng(13)
    alpha = 1e-3
    print(f'{"n":>5} {"r":>5} | FWER: tau     normal   t(0.367 r)   (nominal {alpha})')
    for n, u in ((30, 6), (50, 8), (200, 8)):
        r = n - u
        a_eff = alpha / n
        k = {'tau': tau_crit(a_eff, r), 'normal': stats.norm.ppf(1 - a_eff / 2),
             'mad': stats.t.ppf(1 - a_eff / 2, 0.367 * r)}
        A = rng.standard_normal((n, u))
        H = A @ np.linalg.solve(A.T @ A, A.T)
        qv = 1 - np.diag(H)
        hits = dict.fromkeys(k, 0)
        for _ in range(trials):
            e = rng.standard_normal(n)
            raw = (e - H @ e) / np.sqrt(qv)
            s = 1.4826 * np.median(np.abs(raw - np.median(raw)))
            wmax = np.abs(raw).max() / s
            for name, kc in k.items():
                hits[name] += wmax > kc
        f = {name: h / trials for name, h in hits.items()}
        print(f'{n:5d} {r:5d} |       {f["tau"]:.4f}   {f["normal"]:.4f}   {f["mad"]:.5f}')
        assert f['tau'] > 4 * alpha and f['mad'] < 2 * alpha


if __name__ == '__main__':
    main()
