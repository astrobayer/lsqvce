"""engine.ellipse with Scale 'aposteriori' (default): C = s0^2 Q, k2 = chi2inv(conf, 2).
Coverage of the true (x_i, x_j) by that ellipse, Monte Carlo, vs 2*F(2, r) scaling."""
import numpy as np
from scipy import stats
rs = np.random.RandomState(4)
for m, u in ((8, 3), (13, 3), (33, 3), (203, 3)):
    A = rs.randn(m, u); Q = np.linalg.inv(A.T @ A)[:2, :2]; Qi = np.linalg.inv(Q); r = m - u
    k2 = stats.chi2.ppf(.95, 2); kF = 2 * stats.f.ppf(.95, 2, r); hit = hitF = 0; T = 20000
    for _ in range(T):
        e = rs.randn(m); x = np.linalg.lstsq(A, e, rcond=None)[0]; v = e - A @ x; s2 = v @ v / r
        d = x[:2] @ Qi @ x[:2] / s2; hit += d <= k2; hitF += d <= kF
    print('r = %3d: coverage with chi2 (code) %.3f ; with 2F(2,r) %.3f ; axis ratio sqrt(2F/chi2) = %.3f' % (r, hit / T, hitF / T, np.sqrt(kF / k2)))
