"""WP-119 mirror (R7 residue, §12.48): a DATUM GROUP inside a stacked prior.

WP-114 recognised the prior as a datum only when the WHOLE stacked prior lay in the
data's null space. A datum stacked with a regularising group kept the solve's
covariance -- sigma-dependent on the datum directions although x is not.

The unit of meaning is the user's GROUP (data.priors(k)), not the row: a ridge on
every parameter also has a part on the null space, and that part is a Bayesian
prior (se = its sigma), not a datum. So:

  group g is a DATUM  <=>  the data and every OTHER group are blind to its rows:
      N Kg' = 0   and   Kh Kg' = 0  for all h ~= g
  (v'(N + sum_h Kh'Kh)v = 0 with every term PSD: the same as Kg in the null space of
  N + the other groups.)  Tolerances as WP-114: |N Kg'| <= 1e-10 |N||Kg| (|R Kg'| <=
  1e-8 |R||Kg| on the QR factor), |Ko Kg'| <= 1e-10 |Ko||Kg|.

Kd = the datum groups' rows, Ki = the rest, Ni = N + Ki'Ki. Then N Kd' = 0 and
Ki Kd' = 0, so span(Kd') is the null space of Ni, x is the minimum-constraint
solution of Ni for ANY datum sigma, and its cofactor is pinv(Ni), computed stably
(WP-114's form with Ni for N):
  inv(Ni + c^2 Kd'Kd) - Dm Dm',  Dm = Kd'(Kd Kd')^-1 / c,  c = sqrt(|Ni|_F) / |Kd|_F
on the QR backend from the factor of [R; Ki; c Kd], c = |R|_F / |Kd|_F.
A group only PARTLY in the null space stays information (one declaration, one
meaning): its x depends on its sigma.
"""
import numpy as np

def fixture(rs):
    # levelling loop on 4 heights (null: a common level) + 2 extra parameters
    # (a scale and a tilt) seen by every row and by 3 absolute rows of their own
    u = 4; pairs = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)]
    A = np.zeros((len(pairs) + 3, u + 2))
    for k, (i, j) in enumerate(pairs):
        A[k, i] = -1; A[k, j] = 1; A[k, u] = 0.3 * (k + 1); A[k, u + 1] = (-1) ** k
    A[len(pairs):, u:] = rs.randn(3, 2)
    l = A @ np.r_[0, 1.5, -0.7, 2.2, 0.2, -0.1] + 0.01 * rs.randn(A.shape[0])
    return A, l

def solve(A, l, Kw):
    M = np.vstack([np.column_stack([A, l]), np.column_stack([Kw, np.zeros(Kw.shape[0])])])
    G = np.linalg.qr(M, mode='r'); n = A.shape[1]
    x = np.linalg.solve(G[:n, :n], G[:n, -1])
    Ri = np.linalg.inv(G[:n, :n])
    return x, G[n, -1] ** 2 / (A.shape[0] + Kw.shape[0] - n), Ri @ Ri.T

def datum_mask(Kw, grp, seen):
    """seen(Kg) -> data see the rows; returns the row mask of the datum groups."""
    isd = np.zeros(Kw.shape[0], bool)
    for g in np.unique(grp):
        Kg, Ko = Kw[grp == g], Kw[grp != g]
        if seen(Kg): continue
        if Ko.size and np.linalg.norm(Ko @ Kg.T) > 1e-10 * np.linalg.norm(Ko) * np.linalg.norm(Kg): continue
        isd[grp == g] = True
    return isd

def datum_part_normal(N, Kw, grp):
    isd = datum_mask(Kw, grp, lambda Kg: np.linalg.norm(N @ Kg.T) > 1e-10 * np.linalg.norm(N) * np.linalg.norm(Kg))
    if not isd.any(): return None
    Kd, Ki = Kw[isd], Kw[~isd]
    Ni = N + Ki.T @ Ki
    c = np.sqrt(np.linalg.norm(Ni)) / np.linalg.norm(Kd)
    Rm = np.linalg.cholesky(Ni + c ** 2 * Kd.T @ Kd).T
    return Rm, Kd.T @ np.linalg.inv(Kd @ Kd.T) / c

def datum_part_qr(R, Kw, grp):
    isd = datum_mask(Kw, grp, lambda Kg: np.linalg.norm(R @ Kg.T) > 1e-8 * np.linalg.norm(R) * np.linalg.norm(Kg))
    if not isd.any(): return None
    Kd, Ki = Kw[isd], Kw[~isd]
    c = np.linalg.norm(R) / np.linalg.norm(Kd)
    Rm = np.linalg.qr(np.vstack([R, Ki, c * Kd]), mode='r')[:R.shape[1]]
    return Rm, Kd.T @ np.linalg.inv(Kd @ Kd.T) / c

def cofactor(Rm, Dm):
    Ri = np.linalg.inv(Rm)
    return Ri @ Ri.T - Dm @ Dm.T

rs = np.random.RandomState(0)
A, l = fixture(rs)
N = A.T @ A
R = np.linalg.qr(A, mode='r')
datum = np.r_[np.ones(4), 0, 0][None]
ridgeX = np.array([[0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1.]])
print('nullity of N:', N.shape[0] - np.linalg.matrix_rank(N))

sep = {   # name: (groups as (L, sigma or None = sigma_d), the info rows' weighted L)
    'pure datum (WP-114)':             [(datum, None)],
    'datum + ridge on the extras':     [(datum, None), (ridgeX, 0.5)],
    'ridge on the extras + datum, reversed order': [(ridgeX, 0.5), (datum, None)],
}
for name, groups in sep.items():
    print('\n' + name)
    for sd in (1.0, 1e3, 1e6):
        Ls = [L / (sd if s is None else s) for L, s in groups]
        Kw = np.vstack(Ls); grp = np.concatenate([np.full(L.shape[0], k) for k, L in enumerate(Ls)])
        Ki = np.vstack([L / s for L, s in groups if s is not None] or [np.zeros((0, 6))])
        x, s0sq, Qs = solve(A, l, Kw)
        Qref = np.linalg.pinv(N + Ki.T @ Ki)
        Qn = cofactor(*datum_part_normal(N, Kw, grp)); Qq = cofactor(*datum_part_qr(R, Kw, grp))
        print('  sigma_d %-6g x[0] %.6f  se0 solve %-12.4g datum part %.4f  rel |Q - pinv(N+Ki\'Ki)|: normal %.1e  QR %.1e'
              % (sd, x[0], np.sqrt(s0sq * Qs[0, 0]), np.sqrt(s0sq * Qn[0, 0]),
                 np.abs(Qn - Qref).max() / np.abs(Qref).max(), np.abs(Qq - Qref).max() / np.abs(Qref).max()))

print('\nno datum group -- the solve keeps its covariance (None on both backends)')
info = {
    'ridge on everything (Bayesian on the null part)': [np.eye(6)],
    'datum + ridge on everything':                     [datum, np.eye(6)],
    'datum + ridge on heights 1, 2':                   [datum, np.eye(6)[:2]],
    'ONE group: datum row + extra rows (partial)':     [np.vstack([datum, ridgeX])],
    'one row mixing level and scale':                  [np.r_[np.ones(4), 1, 0][None]],
}
for name, Ls in info.items():
    Kw = np.vstack(Ls); grp = np.concatenate([np.full(L.shape[0], k) for k, L in enumerate(Ls)])
    print('  %-50s normal: %s  QR: %s' % (name, datum_part_normal(N, Kw, grp) is None, datum_part_qr(R, Kw, grp) is None))
# partial group: x depends on its sigma, so the solve's covariance is the honest one
for sd in (1.0, 1e3):
    x, _, _ = solve(A, l, np.vstack([datum, ridgeX]) / sd)
    print('  partial group, sigma %-5g: x[4] = %.6f (moves with sigma: information, not a datum)' % (sd, x[4]))
