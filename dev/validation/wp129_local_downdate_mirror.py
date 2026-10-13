#!/usr/bin/env python3
"""WP-129 mirror: snoop's rank-one downdate of a local block the solve
eliminated EXACTLY on its QR factor (BlockCofactor.downdate).

WP-128 built the cofactor's local half H = R11^-T from the QR factor, so an
ill-conditioned local block (a quadratic polynomial in raw seconds over
1e4 s, cond(B) ~ 1e8) keeps every direction. The downdate, though, still
REFACTORED from the squared block: H = localBlock(Nll - p b b'). From
cond(B) ~ 7e7 its rcond is below eps, the eigen branch runs, and the
smallest eigenvalue of Nll (~ eps * max) comes back with an O(1) relative
error or is cut -- after ONE rejection the cofactor of that block is wrong
again, and the leverage sum over the survivors misses u.

The fix downdates the factor itself: R11~ = cholupdate(R11, sqrt(p) b, '-')
(positive pivots), H = R11~^-T. The commons update and G's rank-one update
are unchanged. Reference: a fresh QR of the design without the row.

Printed per arc length T: cond(B), the leverage sum over the survivors
(must equal u = nc + nl) and max |Qll - Qll_ref| / max |Qll_ref| for the old
refactor and for the downdated factor.

Run:  python3 dev/validation/wp129_local_downdate_mirror.py
Created 2026-10-13 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-129.
"""
import numpy as np

EPS = np.finfo(float).eps


def chol_downdate(R, x):
    """Upper R (positive pivots), returns R~ with R~'R~ = R'R - x x' (LINPACK dchdd idea,
    hyperbolic rotations row by row); None when the result is not positive definite."""
    R = R.copy()
    x = x.astype(float).copy()
    n = R.shape[0]
    for k in range(n):
        r2 = R[k, k] ** 2 - x[k] ** 2
        if r2 <= 0:
            return None
        r = np.sqrt(r2)
        c, s = r / R[k, k], x[k] / R[k, k]
        R[k, k] = r
        if k + 1 < n:
            R[k, k + 1:] = (R[k, k + 1:] - s * x[k + 1:]) / c
            x[k + 1:] = c * x[k + 1:] - s * R[k, k + 1:]
    return R


def local_block_old(Nll):
    """BlockCofactor.localBlock with tol = [] (what downdate called): H'H = pinv(Nll)."""
    S = (Nll + Nll.T) / 2
    nl = S.shape[0]
    if 1 / np.linalg.cond(S) >= EPS:
        try:
            L = np.linalg.cholesky(S)
            return np.linalg.inv(L), 'chol'
        except np.linalg.LinAlgError:
            pass
    dv, V = np.linalg.eigh(S)
    keep = (dv > nl * np.spacing(dv.max())) & (dv > 0)
    H = np.zeros((nl, nl))
    H[:keep.sum()] = (V[:, keep] / np.sqrt(dv[keep])).T
    return H, f'eig, kept {keep.sum()}/{nl}'


def lev_sum(A, B, H, Rc, G):
    """sum_i a_i' Qxx a_i over the rows [A B], from the blocks (BlockCofactor.leverage)."""
    Z = A.T - G.T @ B.T
    hc = np.sum(np.linalg.solve(Rc.T, Z) ** 2)
    return hc + np.sum((H @ B.T) ** 2)


def main():
    rng = np.random.default_rng(129)
    m, nc, p = 200, 4, 1.0
    A = rng.standard_normal((m, nc))
    out = 57                                         # the row snoop rejects
    keep = np.arange(m) != out
    for T in (3e2, 2e3, 1e4, 5e4):
        t = np.linspace(0, T, m)
        B = np.c_[np.ones(m), t, t ** 2]
        nl = B.shape[1]
        # the cofactor before the rejection (WP-128: from the factor)
        R11 = np.linalg.qr(B, mode='r')
        R11 = np.sign(np.diag(R11))[:, None] * R11
        Ri = np.linalg.inv(R11)
        G = Ri @ (Ri.T @ (B.T @ A))
        Rc = np.linalg.qr(A - B @ G, mode='r')
        Rc = np.sign(np.diag(Rc))[:, None] * Rc
        Nll = B.T @ B
        # BlockCofactor.downdate: the commons and G (unchanged by WP-129)
        b, ac = B[out], A[out]
        c = Ri @ (Ri.T @ b)
        d = 1 - p * b @ c
        z = ac - G.T @ b
        pt = p / d
        Rc2 = chol_downdate(Rc, np.sqrt(pt) * z)
        G2 = G - pt * np.outer(c, z)
        # the local half: old refactor vs downdated factor
        H_old, branch = local_block_old(Nll - p * np.outer(b, b))
        R2 = chol_downdate(R11, np.sqrt(p) * b)
        H_new = np.linalg.inv(R2).T
        # reference: a fresh QR without the row
        Rref = np.linalg.qr(B[keep], mode='r')
        Qref = np.linalg.inv(Rref) @ np.linalg.inv(Rref).T
        u = nc + nl
        s_old = lev_sum(A[keep], B[keep], H_old, Rc2, G2)
        s_new = lev_sum(A[keep], B[keep], H_new, Rc2, G2)
        e_old = np.abs(H_old.T @ H_old - Qref).max() / np.abs(Qref).max()
        e_new = np.abs(H_new.T @ H_new - Qref).max() / np.abs(Qref).max()
        print(f'T = {T:7g} s  cond(B) = {np.linalg.cond(B):.1e}  old: {branch:15s} '
              f'sum h {s_old:8.4f}, Qll err {e_old:.1e} | factor: sum h {s_new:.4f}, '
              f'Qll err {e_new:.1e}  (u = {u})')
        assert abs(s_new - u) < 1e-6 and e_new < 1e-6
        if T >= 1e4:
            assert abs(s_old - u) > 0.5 or e_old > 1e-2


if __name__ == '__main__':
    main()
