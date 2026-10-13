"""WP-133 mirror: a deleted epoch's row in the nugget whitener (covModel>nuggetRows).

The nugget whitener X = G F A (R7 WP-115 part 3b) is not causal: the direction a
zero-weight epoch j deletes is q = F e_j / |F e_j|, spread over the rows j, j+1,
... of the whitened arc (an AR(1) law: q_j^2 = 1/(1+phi^2), 0.55 at phi = 0.9).
A row-wise consumer (redundancy, the verdict's replay) that marks ROW j as the
excluded one then misses 1 - r_j of the identity sum(r) = m_live - u.

The fix: an orthogonal T (k Householder reflections) that maps the deleted
directions span(F U_del) onto the unit rows e_del. X'X is unchanged (T'T = I),
so the normals, the solve and the whitened residual sum are those of the
unrotated whitener; the deleted rows of T X are ~0 and every other row is
untouched in aggregate. Checked here:

  1. (TX)'(TX) == X'X == A' Sigma^-1 A                       (rel 1e-12: an
     orthogonal T changes the normals by rounding only)
  2. the deleted rows of TX are ~0 (each |row| <= 1e-5 |X|)  (the 1e-12 floor:
     the row keeps sqrt(1e-12) ~ 1e-6 of its weight)
  3. sum(r_i) over the LIVE rows of TX == m - k - u          (abs 1e-6)
     while the unrotated X misses it by 0.05-0.17 here (3e-4..9e-4 relative,
     against the verdict's 1e-6 integrity line: R.diag would read INVALID)
  4. the same with partial factors mixed in (the loop's passes): the deletion
     directions (nugget q > 1e11, adjustVCEjoint's convention) are aligned,
     the partial ones stay where they are; 1. and 3. hold
"""
import numpy as np
from scipy.linalg import toeplitz, cholesky, solve_triangular, qr, eigh

def nugget_whiten(T, q, A, rotate):
    m = T.shape[0]
    L = cholesky(T, lower=True)
    FA = solve_triangular(L, A, lower=True)
    c = q * np.diag(T)
    z = np.flatnonzero(q > 0)
    if z.size == 0:
        return FA
    E = np.zeros((m, z.size)); E[z, np.arange(z.size)] = 1
    FU = solve_triangular(L, E, lower=True)
    Q, Rq = qr(FU, mode='economic')
    K = np.diag(1 / c[z]) + Rq.T @ Rq
    B = Rq @ np.linalg.solve(K, Rq.T)
    beta, V = eigh((B + B.T) / 2)
    g = 1 - np.sqrt(np.clip(1 - beta, 0, None))
    X = FA - Q @ (V @ (g[:, None] * (V.T @ (Q.T @ FA))))
    if not rotate:
        return X
    dl = q[z] > 1e11                                   # the deletions
    if not dl.any():
        return X
    Qd, _ = qr(FU[:, dl], mode='economic')            # span(F U_del), orthonormal
    zd = z[dl]
    Vh = []
    for i in range(zd.size):                           # Householder: Qd(:,i) -> -s e_zd(i)
        x = Qd[:, i].copy()
        s = 1.0 if x[zd[i]] >= 0 else -1.0
        v = x.copy(); v[zd[i]] += s                    # v = x + s e (no cancellation)
        v /= np.linalg.norm(v)
        Qd = Qd - 2 * np.outer(v, v @ Qd)
        Vh.append(v)
    for v in Vh:
        X = X - 2 * np.outer(v, v @ X)
    return X

def live_sum_r(X, live):
    N = X.T @ X
    h = np.einsum('ij,ji->i', X, np.linalg.solve(N, X.T))
    return np.sum(1 - h[live])

rng = np.random.default_rng(133)
ok = True
m = 200
t = np.arange(m) * 5.0
A = np.c_[np.ones(m), (t - t.mean()) / np.ptp(t), np.sin(2*np.pi*t/1500.), np.cos(2*np.pi*t/700.)]
u = A.shape[1]
for phi in (0.5, 0.9):
    T = 1.7 * toeplitz(phi ** np.arange(m))
    for mixed in (False, True):
        f = np.ones(m)
        z = rng.choice(m, 30, replace=False)
        f[z[:10]] = 1e-12                              # rejected (the loader's floor)
        if mixed:
            f[z[10:]] = rng.uniform(0.05, 0.95, 20)    # partial
        q = 1 / f - 1
        Sig = T + np.diag(q * np.diag(T))
        Nd = A.T @ np.linalg.solve(Sig, A)
        X0 = nugget_whiten(T, q, A, False)
        X1 = nugget_whiten(T, q, A, True)
        e1 = np.linalg.norm(X1.T @ X1 - Nd) / np.linalg.norm(Nd)
        e0 = np.linalg.norm(X0.T @ X0 - Nd) / np.linalg.norm(Nd)
        dele = z[:10]
        rowmax = np.max(np.linalg.norm(X1[dele], axis=1)) / np.max(np.linalg.norm(X1, axis=1))
        live = np.ones(m, bool); live[dele] = False
        s1 = live_sum_r(X1, live); s0 = live_sum_r(X0, live)
        want = m - dele.size - u
        good = e1 < 1e-12 and rowmax < 1e-5 and abs(s1 - want) < 1e-6
        ok &= good
        print(f'phi {phi} mixed {mixed}: normals rel {e1:.1e} (unrotated {e0:.1e}); '
              f'deleted rows {rowmax:.1e}; live sum(r) {s1:.9f} want {want} '
              f'(unrotated {s0:.4f}) {"OK" if good else "FAIL"}')
print('ALL OK' if ok else 'FAILED')
raise SystemExit(0 if ok else 1)
