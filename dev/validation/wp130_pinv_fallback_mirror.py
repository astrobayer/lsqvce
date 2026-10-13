#!/usr/bin/env python3
"""WP-130 mirror: a rank-deficient local block (BlockAccumulator's pinv fallback).

One accidentally duplicated local column (B = [1, t, t^2, 1], t over 1e4 s)
makes R11's pivots fail QRSystem.RankTol. The old fallback squared the block
and took pinv(Nll, 1e-10 * norm): eig(Nll)/max = [2e-31, 1e-16, 1e-9, 1],
so the weak but DETERMINED t^2 direction was cut with the null one; dof
counted mObs - nl. The new route (BlockAccumulator.rankReduced) decides the
rank ON the factor: R11 = U S V', k = #{s > RankTol max s}; the locals live on
V1, the rows U2'[R12 r1] join the commons, dof = mObs - k.

Reference: lstsq on the design without the duplicate (rank-equivalent).
Printed: x error in formal sd, dof, sigma0, sum of the leverages (= rank).

Run:  python3 dev/validation/wp130_pinv_fallback_mirror.py
Created 2026-10-13 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-130.
"""
import numpy as np

RANKTOL = 1e3 * np.finfo(float).eps


def old_route(A, B, l, trunc=1e-10):
    nc, nl = A.shape[1], B.shape[1]
    Nll, Nlc, Ncc = B.T @ B, B.T @ A, A.T @ A
    Z = np.linalg.pinv(Nll, rcond=trunc) @ np.c_[Nlc, B.T @ l]
    Nred, nred = Ncc - Nlc.T @ Z[:, :nc], A.T @ l - Nlc.T @ Z[:, -1]
    x = np.linalg.solve(Nred, nred)
    ltl = l @ l - (B.T @ l) @ Z[:, -1]
    vtv = ltl - nred @ x
    dof = len(l) - nl - nc
    return x, dof, vtv / dof


def new_route(A, B, l):
    m, nc = A.shape
    nl = B.shape[1]
    Rp = np.linalg.qr(np.c_[B, A, l], mode='r')          # locals first, as reducedRows
    R11, R12, r1 = Rp[:nl, :nl], Rp[:nl, nl:nl + nc], Rp[:nl, -1]
    rows = Rp[nl:, nl:]
    U, sv, Vt = np.linalg.svd(R11)
    kp = sv > RANKTOL * sv.max()
    U1, U2, V1 = U[:, kp], U[:, ~kp], Vt.T[:, kp]
    k = kp.sum()
    rowsK = np.r_[U2.T @ np.c_[R12, r1], rows]             # the commons' reduced rows
    Rc = np.linalg.qr(rowsK, mode='r')
    x = np.linalg.solve(Rc[:nc, :nc], Rc[:nc, -1])
    vtv = Rc[nc, -1] ** 2 if Rc.shape[0] > nc else 0.0
    dof = m - k - nc
    # local part on the kept basis: y = V1 S1^-1 (U1'r1 - U1'R12 x), H = S1^-1 V1'
    S1 = sv[kp]
    G = V1 @ ((U1.T @ R12) / S1[:, None])
    H = (V1 / S1).T
    Rci = np.linalg.inv(Rc[:nc, :nc])
    Z = A.T - G.T @ B.T
    h = np.sum((Rci.T @ Z) ** 2, 0) + np.sum((H @ B.T) ** 2, 0)
    return x, dof, vtv / dof, h.sum(), k


def main():
    rng = np.random.default_rng(1)
    m, nc = 200, 4
    A = rng.standard_normal((m, nc))
    t = np.linspace(0, 1e4, m)
    B = np.c_[np.ones(m), t, t ** 2, np.ones(m)]               # the duplicate is the last column
    l = A @ rng.standard_normal(nc) + B[:, :3] @ [1, 1e-3, 1e-7] + 0.01 * rng.standard_normal(m)
    D = np.c_[A, B[:, :3]]
    ref = np.linalg.lstsq(D, l, rcond=None)[0][:nc]
    vr = l - D @ np.linalg.lstsq(D, l, rcond=None)[0]
    dofR = m - D.shape[1]
    sd = np.sqrt(vr @ vr / dofR * np.diag(np.linalg.inv(D.T @ D))[:nc])
    xo, dofo, s2o = old_route(A, B, l)
    xn, dofn, s2n, hs, k = new_route(A, B, l)
    print(f'reference: dof {dofR}, sigma0 {np.sqrt(vr @ vr / dofR):.6f}')
    print(f'old pinv(Nll): x err/sd {np.abs((xo - ref) / sd).max():7.2f}, dof {dofo}, sigma0 {np.sqrt(s2o):.6f}')
    print(f'new on factor: x err/sd {np.abs((xn - ref) / sd).max():7.2e}, dof {dofn}, sigma0 {np.sqrt(s2n):.6f}, '
          f'rank {k}, sum h {hs:.6f} (u = {nc + k})')
    assert np.abs((xo - ref) / sd).max() > 10
    assert np.abs((xn - ref) / sd).max() < 1e-4 and dofn == dofR
    assert abs(np.sqrt(s2n) - np.sqrt(vr @ vr / dofR)) < 1e-9 and abs(hs - (nc + k)) < 1e-8


if __name__ == '__main__':
    main()
