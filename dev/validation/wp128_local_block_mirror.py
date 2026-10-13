#!/usr/bin/env python3
"""WP-128 mirror: the diagnostics cofactor of a local block the solve
eliminated EXACTLY on its QR factor (BlockCofactor, BlockAccumulator).

The QR backend eliminates the locals on the factor whenever the pivots of
R11 pass RankTol (cond(B) up to ~4.5e12). The cofactor was rebuilt from the
SQUARED block Nll = B'B: from cond(B) ~ 7e7 its rcond is below eps, the
eigen branch runs, and the eigenvalues below 1e-10 max are cut -- one
redundancy unit lost per cut direction while the solve kept them all.
A quadratic local polynomial in raw seconds is enough: t^2 over 1e4 s.

The fix forms H = R11^-T and G = R11^-1 R12 from the factor (H'H = inv(Nll)).
Printed per arc length T: cond(B), the branch the old code took, and the
leverage sum (must equal u = nc + nl) for the factor and the old cofactor.

Run:  python3 dev/validation/wp128_local_block_mirror.py
Created 2026-10-13 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-128.
"""
import numpy as np


def main():
    rng = np.random.default_rng(2)
    m, nc = 200, 4
    A = rng.standard_normal((m, nc))
    eps = np.finfo(float).eps
    for T in (3e2, 2e3, 1e4, 5e4):
        t = np.linspace(0, T, m)
        B = np.c_[np.ones(m), t, t ** 2]
        nl = B.shape[1]
        Nll, Nlc, Ncc = B.T @ B, B.T @ A, A.T @ A

        def lev(Qll):
            G = Qll @ Nlc
            Qcc = np.linalg.inv(Ncc - Nlc.T @ G)
            Z = A.T - G.T @ B.T
            return np.sum(Z * (Qcc @ Z)) + np.sum(B.T * (Qll @ B.T))

        R11 = np.linalg.qr(B, mode='r')
        Ri = np.linalg.inv(R11)
        exact = lev(Ri @ Ri.T)
        S = (Nll + Nll.T) / 2
        if 1 / np.linalg.cond(S) >= eps:
            branch, Qold = 'chol', np.linalg.inv(S)
        else:
            dv, V = np.linalg.eigh(S)
            keep = dv > 1e-10 * dv.max()
            branch = f'eig, kept {keep.sum()}/{nl}'
            Qold = (V[:, keep] / dv[keep]) @ V[:, keep].T
        old = lev(Qold)
        print(f'T = {T:7g} s  cond(B) = {np.linalg.cond(B):.1e}  old branch: {branch:15s} '
              f'sum h: factor {exact:.4f}, old {old:.4f}  (u = {nc + nl})')
        assert abs(exact - (nc + nl)) < 1e-6
        if T >= 1e4:
            assert old < nc + nl - 0.5


if __name__ == '__main__':
    main()
