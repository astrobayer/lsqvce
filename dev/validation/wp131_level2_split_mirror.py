#!/usr/bin/env python3
"""WP-131 mirror: the Level-2 contribution split (streamContribution) against
the Level-1 share (redundancy / BlockCofactor.leverage).

One system, two observation types sharing the system's local block (bias,
drift, t^2), grouped by obsType. The reduced common normal N = Ncc - Nlc' G
(G = inv(Nll) Nlc) can be split over groups in more than one way:

  half-Schur (the old Level 2):  N_g = Ncc_g - (Nlc_g' G + G' Nlc_g) / 2
  leverage   (Level 1, new L2):   N_g = Z_g' Z_g,  Z_g = A_g - B_g G
                                      = Ncc_g - Nlc_g' G - G' Nlc_g + G' Mll_g G

Both sum to N; they differ per group by tr(Sxx G' B_g'(A_g - B_g G)), zero
only when the group is the whole system. Level 1 reports t_g = sum_{i in g}
z_i' Sxx z_i -- the leverage split -- so the old Level-2 diagonal did not sum to
it, and a group's share could be negative. The leverage split is positive
semidefinite per group, its diagonal sums to t_g exactly, and with the exact
local term tr(Qll Mll_g) the streamed redundancy n_g - local_g - t_g is the
exact sum of (1 - h_i).

Run:  python3 dev/validation/wp131_level2_split_mirror.py
Created 2026-10-13 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-131.
"""
import numpy as np


def main():
    rng = np.random.default_rng(31)
    nc = 4
    groups = []
    for g, (m, t0, t1, w) in enumerate(((60, 0, 40, 1.0), (40, 25, 100, 3.0))):
        t = np.linspace(t0, t1, m) / 100
        A = rng.standard_normal((m, nc))
        A[:, 0] = 3 * t + 0.05 * rng.standard_normal(m)   # a common column close to the local drift
        B = np.c_[np.ones(m), t, t ** 2]
        groups.append((w * A, w * B))
    Aall = np.vstack([a for a, _ in groups])
    Ball = np.vstack([b for _, b in groups])
    Nll, Nlc = Ball.T @ Ball, Ball.T @ Aall
    G = np.linalg.solve(Nll, Nlc)
    N = Aall.T @ Aall - Nlc.T @ G
    Sxx = np.linalg.inv(N)
    Qll = np.linalg.inv(Nll)
    H = np.c_[Aall, Ball] @ np.linalg.inv(np.c_[Aall, Ball].T @ np.c_[Aall, Ball]) @ np.c_[Aall, Ball].T
    h = np.diag(H)
    row = 0
    tot_half = tot_lev = np.zeros((nc, nc))
    for g, (A, B) in enumerate(groups):
        m = A.shape[0]
        Z = A - B @ G
        t_l1 = np.sum(Z * (Z @ Sxx))                                    # Level 1 (exact)
        Ncc_g, Nlc_g, Mll_g = A.T @ A, B.T @ A, B.T @ B
        Nhalf = Ncc_g - 0.5 * (Nlc_g.T @ G + G.T @ Nlc_g)
        Nlev = Ncc_g - Nlc_g.T @ G - G.T @ Nlc_g + G.T @ Mll_g @ G
        dh, dl = np.diag(Sxx @ Nhalf), np.diag(Sxx @ Nlev)
        locA = np.trace(Qll @ Mll_g)
        r_exact = np.sum(1 - h[row:row + m])
        print(f'group {g}: Level-1 t_g {t_l1:.4f} | half-Schur sum diag {dh.sum():+.4f} '
              f'(min {dh.min():+.3f}) | leverage sum diag {dl.sum():.4f} (min {dl.min():+.3f}) | '
              f'r: exact {r_exact:.4f}, n - tr(Qll Mll_g) - t {m - locA - dl.sum():.4f}')
        assert abs(dl.sum() - t_l1) < 1e-10 and abs(m - locA - dl.sum() - r_exact) < 1e-9
        assert abs(dh.sum() - t_l1) > 0.05
        tot_half, tot_lev = tot_half + Nhalf, tot_lev + Nlev
        row += m
    assert np.allclose(tot_half, N) and np.allclose(tot_lev, N)
    print('both splits sum to N; only the leverage split matches Level 1 per group')


if __name__ == '__main__':
    main()
