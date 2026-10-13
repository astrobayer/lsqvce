#!/usr/bin/env python3
"""WP-130 mirror: the datum test of BlockAccumulator.datumRows under a change
of parameter units.

A prior group is a datum candidate ("blind") when the data cannot see its
rows: S Kg' = 0 (S the data factor R, or N). The test compared ||S Kg'|| with
tol * ||S||_F * ||Kg||_F -- the norm of the WHOLE system -- so a parameter in
small units (its column 1e-7 of the typical one, nc = 400) passed as blind and
got the minimum-constraint covariance: sd 0 against 1.9e3. The rounding of
the product S Kg' is bounded by eps * |S| |Kg'| (elementwise); the new test
compares with THAT norm, which is invariant to column scaling: an exact null
vector stays at ~eps, a scaled but seen column at ~1.

Cases: (a) the WP-127 split datum (translation + rotation of a 2-D distance
net): both groups blind under both tests; (b) the scaled parameter with a
prior of its own: blind under the old test, seen under the new; (c) a ridge
group the data see: seen under both.

Run:  python3 dev/validation/wp130_datum_scale_mirror.py
Created 2026-10-13 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-130.
"""
import itertools

import numpy as np

TOL_QR = 1e-8


def old_blind(S, Kg):
    return np.linalg.norm(S @ Kg.T) <= TOL_QR * np.linalg.norm(S) * np.linalg.norm(Kg)


def new_blind(S, Kg):
    return np.linalg.norm(S @ Kg.T) <= TOL_QR * np.linalg.norm(np.abs(S) @ np.abs(Kg.T))


def net():
    P = np.array([[20, 5], [28, 6], [30, 13], [24, 16], [19, 11]], float)
    rows = []
    for i, j in itertools.combinations(range(5), 2):
        e = (P[j] - P[i]) / np.linalg.norm(P[j] - P[i])
        r = np.zeros(10)
        r[2 * i:2 * i + 2], r[2 * j:2 * j + 2] = -e, e
        rows.append(r)
    Kt = np.kron(np.ones((1, 5)), np.eye(2))
    Kr = np.reshape(np.vstack([-P[:, 1], P[:, 0]]).T, (1, -1))
    return np.array(rows), Kt, Kr


def main():
    A, Kt, Kr = net()
    R = np.linalg.qr(A, mode='r')
    for name, K in (('translation', Kt), ('rotation', Kr)):
        o, n = old_blind(R, K), new_blind(R, K)
        print(f'(a) {name:11s} old blind {o}, new blind {n}')
        assert o and n
    rng = np.random.default_rng(3)
    m, nc, sc = 3000, 400, 1e-7
    A = rng.standard_normal((m, nc))
    A[:, -1] *= sc
    R = np.linalg.qr(A, mode='r')
    e = np.zeros((1, nc))
    e[0, -1] = 1
    o, n = old_blind(R, e), new_blind(R, e)
    sd_true = np.sqrt(np.linalg.inv(R.T @ R + (e.T @ e) / (1e3 / sc) ** 2)[-1, -1])
    print(f'(b) scaled column 1e-7: old blind {o} (sd reported 0), new blind {n} (sd {sd_true:.4g})')
    assert o and not n
    ridge = np.zeros((2, nc))
    ridge[0, 0] = ridge[1, 3] = 1
    o, n = old_blind(R, ridge), new_blind(R, ridge)
    print(f'(c) ridge group:        old blind {o}, new blind {n}')
    assert not o and not n


if __name__ == '__main__':
    main()
