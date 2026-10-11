"""WP-112 mirror: the accumulator path in SQUARE-ROOT form (R7-03, R7-16).

The new path, line by line as it will be built in MATLAB:
  loc = R of [A B l] (QRSystem.add, commons first, as accumulated)
  preeliminate (Backend "qr"): re-triangularise the (d+1) x (d+1) factor with
      the LOCAL columns first, Rp = qr(R(:, [il ic end])) -- an orthogonal
      transform of R, so Rp' Rp = P' [A B l]'[A B l] P without ever forming it;
      then [R22 r2; 0 rho] is the square root of the reduced contribution
      (R22' R22 = Ncc - Nlc' Nll^-1 Nlc, rho^2 + r2' r2 = ltl_red).
  merge: QR of the stacked rows -- no chol, no eig.
  solve: x = R \\ r, refused when min|diag R| <= TOL * max|diag R| (the rank
      decision on the pivots, by the code's threshold).
Compared with the shipped path (normal-form Schur + chol in merge, the review's
m1 mirror) and with the dense reference lstsq([A B], l).
"""
import numpy as np
import scipy.linalg as sla

eps = np.finfo(float).eps


def qr_r(M, rows):
    R = np.triu(np.linalg.qr(M, mode='r'))
    return R[:min(R.shape[0], rows), :]


def accumulate(A, B, l, batch=100):
    d = A.shape[1] + B.shape[1]
    R = np.zeros((0, d + 1))
    M = np.column_stack([A, B, l])
    for i in range(0, len(l), batch):
        R = qr_r(np.vstack([R, M[i:i + batch]]), d + 1)
    return R


def reduced_rows(R, nc, nl):
    """[R22 r2; 0 rho] from the commons-first factor R of [A B l]."""
    d = nc + nl
    perm = list(range(nc, d)) + list(range(nc)) + [d]
    Rp = qr_r(R[:, perm], d + 1)
    rows = Rp[nl:, nl:]                       # (nc+1 or fewer) x (nc+1)
    return rows, Rp[:nl, :nl]


def new_path(A, B, l, tol):
    nc, nl = A.shape[1], B.shape[1]
    R = accumulate(A, B, l)
    rows, R11 = reduced_rows(R, nc, nl) if nl else (R, None)
    G = qr_r(rows, nc + 1)                    # merge into an empty common system
    if G.shape[0] < nc:
        return None, 'rankDeficient', np.nan
    dR = np.abs(np.diag(G[:nc, :nc]))
    piv = dR.min() / dR.max()
    if piv <= tol:
        return None, 'singularR', piv
    x = sla.solve_triangular(G[:nc, :nc], G[:nc, -1])
    return x, 'ok', piv


def shipped_path(A, B, l):
    nc, nl = A.shape[1], B.shape[1]
    R = accumulate(A, B, l)
    d = nc + nl
    R1 = R[:d, :d]; r = R[:d, -1]
    N = R1.T @ R1; nn = R1.T @ r
    Ncc, nc_ = N[:nc, :nc], nn[:nc]
    if nl:
        Z = np.linalg.solve(N[nc:, nc:], np.column_stack([N[nc:, :nc], nn[nc:]]))
        Ncc = Ncc - N[nc:, :nc].T @ Z[:, :nc]; nc_ = nc_ - N[nc:, :nc].T @ Z[:, -1]
    Ncc = (Ncc + Ncc.T) / 2
    try:
        Rc = np.linalg.cholesky(Ncc).T; rhs = sla.solve_triangular(Rc.T, nc_, lower=True)
    except np.linalg.LinAlgError:
        dv, V = np.linalg.eigh(Ncc)
        keep = (dv > nc * np.spacing(dv.max())) & (dv > 0)
        Rc = np.diag(np.sqrt(dv[keep])) @ V[:, keep].T; rhs = (V[:, keep].T @ nc_) / np.sqrt(dv[keep])
    G = qr_r(np.column_stack([Rc, rhs]), nc + 1)
    if G.shape[0] < nc:
        return None, 'rankDeficient'
    dR = np.abs(np.diag(G[:nc, :nc]))
    if dR.min() <= nc * np.spacing(dR.max()):
        return None, 'singularR'
    return sla.solve_triangular(G[:nc, :nc], G[:nc, -1]), 'ok'


def fixture(rs, k, nl, m=800, nc=8, noise=1e-4):
    U, _ = np.linalg.qr(rs.randn(m, nc)); V, _ = np.linalg.qr(rs.randn(nc, nc))
    A = U @ np.diag(np.logspace(0, -k, nc)) @ V.T
    B = rs.randn(m, nl)
    l = np.column_stack([A, B]) @ rs.randn(nc + nl) + noise * rs.randn(m)
    return A, B, l


def formal_err(x, A, B, l):
    """max |x - x_ref| / sd_ref, the dense reference and its formal sd."""
    nc = A.shape[1]
    M = np.column_stack([A, B])
    Q, Rd = np.linalg.qr(M)
    xr = sla.solve_triangular(Rd, Q.T @ l)
    v = l - M @ xr
    s0 = np.sqrt(v @ v / (len(l) - M.shape[1]))
    Ri = sla.solve_triangular(Rd, np.eye(M.shape[1]))
    sd = s0 * np.sqrt(np.sum(Ri ** 2, axis=1))
    return np.max(np.abs(x - xr[:nc]) / sd[:nc])


def network(rs, u, datum):
    pairs = [(i, (i + 1) % u) for i in range(u)]
    for _ in range(rs.randint(1, u)):
        i, j = rs.choice(u, 2, replace=False); pairs.append((i, j))
    A = np.zeros((len(pairs) + (1 if datum else 0), u))
    for k, (i, j) in enumerate(pairs):
        A[k, i] = -1; A[k, j] = 1
    if datum:
        A[-1, 0] = 1                         # one height fixed as an observation
    A *= rs.uniform(0.5, 3.0)
    l = A @ (rs.randn(u) * 10) + 0.01 * rs.randn(A.shape[0])
    return A, np.zeros((A.shape[0], 0)), l


def main():
    TOL = None
    print('== pivot ratios min|diag R|/max on the NEW path')
    rs = np.random.RandomState(1)
    free = [new_path(*network(rs, rs.randint(4, 12), False), tol=0)[2] for _ in range(2000)]
    rs = np.random.RandomState(2)
    fixed = [new_path(*network(rs, rs.randint(4, 12), True), tol=0)[2] for _ in range(2000)]
    print('free networks  : max %.2e  median %.2e' % (np.nanmax(free), np.nanmedian(free)))
    print('with datum     : min %.2e' % np.nanmin(fixed))
    for k in (8, 9, 10, 11, 12):
        rs = np.random.RandomState(100 + k)
        p = [new_path(*fixture(rs, k, nl), tol=0)[2] for nl in (0, 3) for _ in range(100)]
        print('cond 1e%-2d      : min %.2e' % (k, np.nanmin(p)))
    TOL = 1e3 * eps
    print('\n== the decision at TOL = 1e3 eps = %.2e' % TOL)
    for name, seed, datum in (('free', 1, False), ('datum', 2, True)):
        rs = np.random.RandomState(seed)
        st = [new_path(*network(rs, rs.randint(4, 12), datum), tol=TOL)[1] for _ in range(2000)]
        print('%-6s networks: ok %4d  refused %4d  (2000)' % (name, st.count('ok'), 2000 - st.count('ok')))
    rs = np.random.RandomState(1)
    st = [shipped_path(*network(rs, rs.randint(4, 12), False))[1] for _ in range(2000)]
    print('free networks, SHIPPED path: ok %d  refused %d' % (st.count('ok'), 2000 - st.count('ok')))
    print('\n== x_c in formal-sd units vs the dense [A B] \\\\ l (20 seeds), and refusals of full-rank problems')
    print('%5s %3s | %-28s | %-28s' % ('cond', 'nl', 'NEW: median / max / refused', 'SHIPPED: median / max / refused'))
    for k, noise in ((7, 1e-10), (8, 1e-4), (8.5, 1e-4), (9, 1e-4)):
        for nl in (0, 3):
            rs = np.random.RandomState(int(10 * k) + nl)
            en, es, rn, rsh = [], [], 0, 0
            for s in range(20):
                A, B, l = fixture(rs, k, nl, noise=noise)
                xn, sn, _ = new_path(A, B, l, TOL)
                xs, ss = shipped_path(A, B, l)
                if sn == 'ok': en.append(formal_err(xn, A, B, l))
                else: rn += 1
                if ss == 'ok': es.append(formal_err(xs, A, B, l))
                else: rsh += 1
            f = lambda e: ('%.1e / %.1e' % (np.median(e), np.max(e))) if e else '- / -'
            print('1e%-3g %3d | %-20s %2d     | %-20s %2d' % (k, nl, f(en), rn, f(es), rsh))


if __name__ == '__main__':
    main()
