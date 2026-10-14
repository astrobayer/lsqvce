#!/usr/bin/env python3
"""WP-135 mirror: gcvLambda (third whole-code review, finding 16).

(a) The functional. The code minimised ||f.*b||^2 / (sum f)^2 with b = V'n --
    "GCV" of N y = n read as if N were the design matrix. The generalised
    cross-validation of the local problem  min ||B y - l||^2 + lam ||L y||^2
    (Golub, Heath & Wahba 1979) is

        G(lam) = ||B y_lam - l||^2 / (m - tr H_lam)^2,

    and in the generalised eigenbasis of (N, L'L) (V'L'LV = I, V'NV = diag d,
    b = V'n, f = lam / (d + lam)):

        ||B y - l||^2 = rss_perp + sum (b^2/d) f^2,   rss_perp = l'l - sum b^2/d
        tr H          = nl - sum f.

    It needs m (rows) and l'l, which the pre-elimination holds. The commons
    enter at their a-priori values (x = 0 of the linearised problem): the
    lambda is chosen per system, before x is known.
(b) The grid. L'L + eps*I gave a direction L does not penalise (a first
    difference, a subset) the eigenvalue ~|N|/eps, and logspace(dmax*1e-12,
    dmax) put every candidate >= 1e3 |N|; worse, the Cholesky form of a pencil
    with cond(L'L) ~ 1/eps returned wrong FINITE eigenvalues too (one negative
    for a PSD N). Now the standard form: null(L) (eigenvalues of L'L below
    1e3 nl eps |L'L|) is fitted freely and eliminated like a common (a Schur
    step on N, n, l'l), and on the penalised complement L'L = diag(e) turns
    the pencil into one symmetric eig; the grid spans its eigenvalues.

Run:  python3 dev/validation/wp135_gcv_mirror.py
Created 2026-10-14 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-135.
"""
import numpy as np
from scipy.linalg import eigh

EPS = np.finfo(float).eps


def gcv_old(N, n, L):
    nlp = N.shape[0]
    LtL = L.T @ L + EPS * np.eye(nlp)
    d, V = eigh((N + N.T) / 2, (LtL + LtL.T) / 2)
    b = V.T @ n
    dmax = d.max()
    grid = np.logspace(np.log10(dmax * 1e-12), np.log10(dmax), 80)
    g = np.array([np.linalg.norm(lam / (d + lam) * b) ** 2 / np.sum(lam / (d + lam)) ** 2
                  for lam in grid])
    return grid[np.argmin(g)], grid


def gcv_new(N, n, L, m, ltl):
    """The MATLAB port (gcvLambda after WP-135): standard form, null(L) eliminated."""
    nlp = N.shape[0]
    N = (N + N.T) / 2
    LtL = L.T @ L
    e, U = np.linalg.eigh((LtL + LtL.T) / 2)
    emax = max(e.max(), 0.0)
    # tol: eig of L'L is accurate to ~nl*eps*|L'L|; 1e3 above it is null(L)
    pen = e > 1e3 * nlp * EPS * emax if emax > 0 else np.zeros(nlp, bool)
    if not pen.any():
        return 0.0, None, None
    W, Z, eW = U[:, pen], U[:, ~pen], e[pen]
    nw, nz = W.T @ n, Z.T @ n
    Nww, Nwz, Nzz = W.T @ N @ W, W.T @ N @ Z, Z.T @ N @ Z
    if Z.shape[1] > 0:                               # null(L): fitted freely -- a Schur step
        cz = np.linalg.cholesky(Nzz)                 # refused when the data miss null(L)
        G = np.linalg.solve(cz, Nwz.T)               # cz^-1 Nzw
        g = np.linalg.solve(cz, nz)
        Nww = Nww - G.T @ G
        nw = nw - G.T @ g
        ltl = ltl - g @ g
    sc = 1 / np.sqrt(eW)                             # L'L = diag(eW) on W: the standard form
    d, Q = np.linalg.eigh((sc[:, None] * Nww) * sc[None, :])
    b = Q.T @ (sc * nw)
    dmax = d.max()
    if not dmax > 0:
        return 0.0, None, None
    ok = d > nlp * EPS * dmax                        # identified by the data
    rss_perp = max(ltl - np.sum(b[ok] ** 2 / d[ok]), 0.0)
    grid = np.logspace(np.log10(dmax * 1e-12), np.log10(dmax), 80)
    gv = np.empty(grid.size)
    for k, lam in enumerate(grid):
        f = np.where(ok, lam / (np.where(ok, d, 0.0) + lam), 1.0)
        rss = rss_perp + np.sum(np.where(ok, b ** 2 / np.where(ok, d, 1.0), 0.0) * f ** 2)
        gv[k] = rss / max((m - nlp + np.sum(f)) ** 2, np.finfo(float).tiny)
    k = int(np.argmin(gv))
    return grid[k], grid, gv


def gcv_dense(B, l, L, grid):
    """Golub-Heath-Wahba on the same grid, by direct solves."""
    m = B.shape[0]
    N = B.T @ B
    n = B.T @ l
    out = []
    for lam in grid:
        M = N + lam * (L.T @ L)
        y = np.linalg.solve(M, n)
        r = B @ y - l
        out.append(r @ r / (m - np.trace(np.linalg.solve(M, N))) ** 2)
    return np.array(out)


def main():
    ok = True
    rng = np.random.default_rng(135)
    print('== (a) the functional: polynomial local blocks, white noise ==')
    for m, deg in ((500, 5), (500, 7), (2000, 6)):
        t = np.linspace(0, 1, m)
        B = np.vander(t, deg + 1, increasing=True)
        L = np.eye(B.shape[1])
        ytrue = rng.standard_normal(B.shape[1])
        l = B @ ytrue + rng.standard_normal(m)
        N, n = B.T @ B, B.T @ l
        lo, _ = gcv_old(N, n, L)
        ln, grid, g = gcv_new(N, n, L, m, l @ l)
        gd = gcv_dense(B, l, L, grid)
        rel = np.max(np.abs(g / gd - 1))
        kd = int(np.argmin(gd))
        err = lambda lam: np.linalg.norm(np.linalg.solve(N + lam * L.T @ L, n) - ytrue)
        print(f'  m {m} deg {deg} cond {np.linalg.cond(N):.1e}: old lambda {lo:.2e} '
              f'(|y-y_true| {err(lo):.3g}), new {ln:.2e} ({err(ln):.3g}), dense GCV {grid[kd]:.2e}; '
              f'curve vs dense max rel {rel:.1e}')
        # tol: the eigenbasis curve against direct solves of cond up to ~1e10
        # normals: the dense trace / solve loses ~cond*eps relative, 1e-5 bounds it
        ok &= rel < 1e-5 and ln == grid[kd]
    print('== (b) RegMatrix with a null space ==')
    m, nl = 400, 6
    t = np.linspace(0, 1, m)
    B = np.c_[np.ones(m), t, np.cos(2 * np.pi * t), np.sin(2 * np.pi * t), t ** 2,
              np.cos(4 * np.pi * t)]
    ytrue = np.array([1.0, 0.5, -0.3, 0.8, -0.2, 0.4])
    l = B @ ytrue + 0.05 * rng.standard_normal(m)
    N, n = B.T @ B, B.T @ l
    normN = np.linalg.norm(N, 2)
    for name, L in (('identity', np.eye(nl)), ('first difference', np.diff(np.eye(nl), axis=0)),
                    ('subset (last 2)', np.eye(nl)[-2:])):
        _, go = gcv_old(N, n, L)
        ln, grid, g = gcv_new(N, n, L, m, l @ l)
        gd = gcv_dense(B, l, L, grid)
        rel = np.max(np.abs(g / gd - 1))
        print(f'  {name:17s}: old grid min {go[0] / normN:.1e} |N|, new {grid[0] / normN:.1e} |N|; '
              f'new lambda {ln:.2e} = dense {grid[int(np.argmin(gd))]:.2e}; curve rel {rel:.1e}')
        ok &= grid[0] < 1e-6 * normN and ln == grid[int(np.argmin(gd))] and rel < 1e-6
    print('ALL CHECKS PASS' if ok else 'A CHECK FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
