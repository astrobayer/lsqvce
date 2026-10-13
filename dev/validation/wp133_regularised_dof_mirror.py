#!/usr/bin/env python3
"""WP-133 mirror: the redundancy of an adjustment whose LOCAL block was
regularised (review 3 finding 6; BlockAccumulator 'Reg').

A ridge on a local block, N = N_data + blkdiag(0, Lambda), makes the hat
matrix's trace u - tr(Q Lambda), so sum_i r_i = m - u + tr(Q Lambda) with
Q = inv(N) -- not m - u. The verdict counted D.r = m - u, so its own
integrity check (sum r against D.r) read SUM_R ERROR and the verdict INVALID
for every regularised solve. tr(Q Lambda) is per system tr(Q_ll Lambda_s),
Q_ll the LOCAL block of the full cofactor = inv(Nll_r) + G Qcc G' (what
BlockAccumulator.localCovariance(id, Qcc, 1) returns, Qcc the commons'
cofactor). tsvd floored the null directions instead of dropping them; on
the kept eigenbasis (the pinv the help promised) the rank counts exactly.

Checked: two systems with locals, one ridge (tikhonov), one lm block, and a
common prior; sum r_i + r_K against m + nk - u + defect; the normal-form
vtv of the ridge solve is |v|^2 + |v_K|^2 + y'Lambda y (the solve subtracts
the penalty so that its sigma0 is the residual pass's); and tsvd's
truncation against pinv.

Run:  python3 dev/validation/wp133_regularised_dof_mirror.py
Created 2026-10-13 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-133.
"""
import numpy as np


def main():
    rng = np.random.default_rng(133)
    nc, nl1, nl2, m1, m2 = 4, 3, 2, 60, 50
    A1, B1 = rng.standard_normal((m1, nc)), rng.standard_normal((m1, nl1))
    B1[:, 2] = B1[:, 1] + 1e-4 * rng.standard_normal(m1)          # weakly determined
    A2, B2 = rng.standard_normal((m2, nc)), rng.standard_normal((m2, nl2))
    u = nc + nl1 + nl2
    Afull = np.zeros((m1 + m2, u))
    Afull[:m1, :nc] = A1; Afull[:m1, nc:nc + nl1] = B1
    Afull[m1:, :nc] = A2; Afull[m1:, nc + nl1:] = B2
    K = rng.standard_normal((2, nc)); sK = 0.5                     # a common prior, 2 rows
    lam1 = 1e-2 * np.eye(nl1)                                      # tikhonov on block 1
    Nll2 = B2.T @ B2
    lam2 = 1e-3 * np.diag(np.diag(Nll2))                           # lm on block 2
    Lam = np.zeros((u, u))
    Lam[nc:nc + nl1, nc:nc + nl1] = lam1
    Lam[nc + nl1:, nc + nl1:] = lam2
    Kf = np.zeros((2, u)); Kf[:, :nc] = K / sK
    N = Afull.T @ Afull + Kf.T @ Kf + Lam
    Q = np.linalg.inv(N)
    r = 1 - np.einsum('ij,jk,ik->i', Afull, Q, Afull)              # data rows
    rK = 2 - np.trace(Kf @ Q @ Kf.T)                               # prior rows
    lhs = r.sum() + rK
    # the defect from the BLOCKS, as the accumulator computes it
    Qcc_full = Q[:nc, :nc]
    defect = 0.0
    for (B, A, L, sl) in ((B1, A1, lam1, slice(nc, nc + nl1)), (B2, A2, lam2, slice(nc + nl1, u))):
        Nll_r = B.T @ B + L
        G = np.linalg.solve(Nll_r, B.T @ A)
        Qll = np.linalg.inv(Nll_r) + G @ Qcc_full @ G.T             # localCovariance(id, Qcc, 1)
        assert np.allclose(Qll, Q[sl, sl], rtol=1e-10, atol=1e-12)
        defect += np.sum(Qll * L)
    rhs = (m1 + m2) + 2 - u + defect
    print(f'ridge + lm + prior: sum r + r_K = {lhs:.12f}; m + nk - u + defect = {rhs:.12f} '
          f'(defect {defect:.6f}); without the defect the check is off by {defect / (m1 + m2 + 2 - u):.2e}')
    assert abs(lhs - rhs) < 1e-9
    # the residual MASS: the normal-form vtv of a ridge solve is the PENALISED
    # objective l'l - b'x = |v|^2 + |v_K|^2 + y'Lambda y; the residual pass and
    # the verdict sum the data (and prior) rows only, so the solve subtracts
    # y'Lambda y -- else residuals() reads a sigma0 mismatch of y'Lambda y / vtv
    l = Afull @ rng.standard_normal(u) + rng.standard_normal(m1 + m2)
    x0K = rng.standard_normal(2)
    xh = np.linalg.solve(N, Afull.T @ l + Kf.T @ (x0K / sK))
    v = Afull @ xh - l; vK = Kf @ xh - x0K / sK
    pen = xh @ Lam @ xh
    vtv_nf = l @ l + (x0K / sK) @ (x0K / sK) - (Afull.T @ l + Kf.T @ (x0K / sK)) @ xh
    print(f'normal-form vtv {vtv_nf:.12f} = |v|^2 + |v_K|^2 + y\'Lambda y '
          f'{v @ v + vK @ vK + pen:.12f} (penalty {pen:.3e})')
    assert abs(vtv_nf - (v @ v + vK @ vK + pen)) < 1e-9 * vtv_nf
    # tsvd: floor vs truncation on an exactly rank-deficient block
    B = rng.standard_normal((40, 3)); B[:, 2] = B[:, 0] + B[:, 1]
    A = rng.standard_normal((40, nc)); l = rng.standard_normal(40)
    Nll = B.T @ B
    d, V = np.linalg.eigh(Nll)
    fl = 1e-10 * d.max()
    keep = d >= fl
    V1, d1 = V[:, keep], d[keep]
    # the reduced commons with the truncated (pinv) block
    G = V1 @ ((V1.T @ (B.T @ A)) / d1[:, None])
    Ncc = A.T @ A - (B.T @ A).T @ G
    ncc = A.T @ l - (B.T @ A).T @ (V1 @ ((V1.T @ (B.T @ l)) / d1))
    x = np.linalg.solve(Ncc, ncc)
    xref = np.linalg.lstsq(np.c_[A, B], l, rcond=None)[0][:nc]
    print(f'tsvd on the kept basis: rank {keep.sum()} of 3, |x - lstsq| = {np.abs(x - xref).max():.1e}')
    assert keep.sum() == 2 and np.abs(x - xref).max() < 1e-10


if __name__ == '__main__':
    main()
