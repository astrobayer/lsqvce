#!/usr/bin/env python3
"""WP-125 mirror: the level of a stream WITHOUT arc factors -- the arithmetic
mean of the arcs' Helmert ratios against the pooled Helmert ratio.

Each arc a has m_a white observations of variance sigma^2 and p local
parameters (a polynomial bias/drift/...), so its redundancy is r_a = m_a - p
and its weighted residual sum Omega_a = v'v ~ sigma^2 chi^2(r_a). With the
declared level 1 the arc's Helmert factor is alpha_a^2 = Omega_a / r_a.

  mean of ratios   s2_mean   = mean_a(Omega_a / r_a)
  pooled ratio     s2_pooled = sum_a Omega_a / sum_a r_a   (the Helmert
                                estimator of ONE component)

Both are unbiased for exact redundancies; the pooled one has the smaller
variance whenever the arcs differ in length (it weights an arc by its
redundancy, the mean weights a 30-epoch arc like a 3000-epoch one). With the
probe-estimated redundancies of the joint estimator (r_a + noise) the mean of
ratios is also BIASED (E[1/r] > 1/E[r], Jensen). Both effects are measured.

Run:  python3 dev/validation/wp125_level_mirror.py
Created 2026-10-12 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-125.
"""
import numpy as np


def one(rng, lens, p, sigma2, probe_rel):
    om, r = [], []
    for m in lens:
        t = np.linspace(-1, 1, m)
        B = np.vander(t, p, increasing=True)
        e = rng.standard_normal(m) * np.sqrt(sigma2)
        v = e - B @ np.linalg.lstsq(B, e, rcond=None)[0]
        ra = m - p
        if probe_rel > 0:                        # Hutchinson-like noise on r_a
            ra = ra * (1 + probe_rel * rng.standard_normal())
        om.append(v @ v)
        r.append(ra)
    om, r = np.array(om), np.array(r)
    return np.mean(om / r), om.sum() / r.sum()


def study(label, lens, p=6, probe_rel=0.0, n=4000, seed=125):
    rng = np.random.default_rng(seed)
    res = np.array([one(rng, lens, p, 1.0, probe_rel) for _ in range(n)])
    mean_, pool = res[:, 0], res[:, 1]
    print('%-46s  mean-of-ratios: bias %+.4f sd %.4f | pooled: bias %+.4f sd %.4f'
          % (label, mean_.mean() - 1, mean_.std(), pool.mean() - 1, pool.std()))
    return mean_, pool


if __name__ == '__main__':
    print('truth sigma^2 = 1; %d arcs per stream; p = 6 locals per arc' % 20)
    study('equal arcs (20 x 400), exact r', [400] * 20)
    study('mixed arcs (10 x 30, 10 x 2000), exact r', [30] * 10 + [2000] * 10)
    study('mixed arcs, probe noise 10% on r', [30] * 10 + [2000] * 10, probe_rel=0.10)
    study('equal arcs, probe noise 10% on r', [400] * 20, probe_rel=0.10)
