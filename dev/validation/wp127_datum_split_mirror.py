#!/usr/bin/env python3
"""WP-127 mirror: a free network's datum declared as SEVERAL prior groups
(BlockAccumulator.datumRows / datumPart).

A 2-D distance network, not centred, with its datum split into a translation
group Kt and a rotation group Kr. Off-centre Kt Kr' = -sum y_i (nonzero), so
the WP-119 rule -- a group is a datum when the data AND every OTHER group are
blind to it -- refused both halves, and the covariance fell back to
inv(N + K'K/sigma^2), which grows with the datum sigma.

The WP-127 rule: a group is BLIND when the data cannot see it; a blind group
is a datum when it is orthogonal to the rows of the NON-blind groups only.
The cofactor is then pinv(Ni), Ni = N + Ki'Ki, computed as
    inv(Ni + c^2 Kd'Kd) - Dm Dm',  Dm = Kd' pinv(Kd Kd') / c.
A ridge group the data sees stays information.

Printed: the datum rows found and max |C - pinv(N)| under both rules.
Run:  python3 dev/validation/wp127_datum_split_mirror.py
Created 2026-10-12 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-127.
"""
import itertools

import numpy as np


def network():
    """The test_fit fixture: five points, every distance, linearised."""
    P = np.array([[20, 5], [28, 6], [30, 13], [24, 16], [19, 11]], float)
    rows = []
    for i, j in itertools.combinations(range(5), 2):
        e = P[j] - P[i]
        e /= np.linalg.norm(e)
        r = np.zeros(10)
        r[2 * i:2 * i + 2] = -e
        r[2 * j:2 * j + 2] = e
        rows.append(r)
    A = np.array(rows)
    Kt = np.kron(np.ones((1, 5)), np.eye(2))
    Kr = np.reshape(np.vstack([-P[:, 1], P[:, 0]]).T, (1, -1))
    return A, Kt, Kr


def datum_rows(K, grp, S, tol_s, old):
    n_s = np.linalg.norm(S)
    blind = {g: np.linalg.norm(S @ K[grp == g].T) <= tol_s * n_s * np.linalg.norm(K[grp == g])
             for g in np.unique(grp)}
    is_d = np.zeros(len(grp), bool)
    for g in np.unique(grp):
        if not blind[g]:
            continue
        Kg = K[grp == g]
        other = (grp != g) if old else np.array([not blind[h] for h in grp])
        Ko = K[other]
        if Ko.size and np.linalg.norm(Ko @ Kg.T) > 1e-10 * np.linalg.norm(Ko) * np.linalg.norm(Kg):
            continue
        is_d[grp == g] = True
    return is_d


def cofactor(K, is_d, sig, N):
    Ks = K / sig
    if not is_d.any():
        return np.linalg.inv(N + Ks.T @ Ks)
    Kd, Ki = Ks[is_d], Ks[~is_d]
    Ni = N + Ki.T @ Ki
    c = np.sqrt(np.linalg.norm(Ni)) / np.linalg.norm(Kd)
    Dm = Kd.T @ np.linalg.pinv(Kd @ Kd.T) / c
    return np.linalg.inv(Ni + c ** 2 * Kd.T @ Kd) - Dm @ Dm.T


def main():
    A, Kt, Kr = network()
    N = A.T @ A
    K = np.vstack([Kt, Kr])
    grp = np.array([1, 1, 2])
    print("Kt Kr' =", (Kt @ Kr.T).ravel())
    for sig in (1.0, 1e3):
        for old in (True, False):
            is_d = datum_rows(K, grp, N, 1e-10, old)
            err = np.abs(cofactor(K, is_d, sig, N) - np.linalg.pinv(N)).max()
            print(f"sigma {sig:g}  {'WP-119' if old else 'WP-127'}: datum rows {is_d}, "
                  f"max|C - pinv(N)| = {err:.2e}")
            assert (err > 0.1) if old else (err < 1e-8)
    ridge = np.zeros((2, 10))
    ridge[0, 0] = ridge[1, 3] = 1
    is_d = datum_rows(np.vstack([K, ridge]), np.array([1, 1, 2, 3, 3]), N, 1e-10, False)
    print("with a ridge group the data sees:", is_d)
    assert not is_d.any()


if __name__ == "__main__":
    main()
