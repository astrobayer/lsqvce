"""(a) adjustRobust: the reported sigma0 of the reweighted problem, dof = m - u
    (zero-factor rows, floored at 1e-12, still counted as observations).
(b) stat kernels vs scipy: tauCrit (Pope), delta0 (Baarda), chi2Inv, fInv, ellipse k.
(c) Results.compare: df2 = diag.n - numel(R.x) (commons only) vs the true n - u_full."""
import numpy as np
from scipy import stats, special

def robust_weights(vw, infl, tune=None):
    med = np.median(vw); s = 1.4826 * np.median(np.abs(vw - med))
    if not (s > 0 and np.isfinite(s)):
        s = np.std(vw, ddof=1)
        if not (s > 0 and np.isfinite(s)): s = 1.0
    w = vw / s; aw = np.abs(w); n = len(w)
    if infl == 'huber':
        c = tune or 1.345; f = np.ones(n); f[aw > c] = c / aw[aw > c]
    elif infl == 'igg3':
        k0, k1 = tune or (1.5, 3.0); f = np.ones(n); mid = (aw > k0) & (aw <= k1)
        f[mid] = (k0 / aw[mid]) * ((k1 - aw[mid]) / (k1 - k0))**2; f[aw > k1] = 0
    else:
        c = tune or 4.685; f = np.zeros(n); mk = aw < c; f[mk] = (1 - (aw[mk] / c)**2)**2
    return f, s

def irls(A, l, infl, maxit=30, tol=1e-8):
    m, u = A.shape; fac = np.ones(m); x = None
    for it in range(maxit):
        p = np.maximum(fac, 1e-12); sw = np.sqrt(p)
        xn = np.linalg.lstsq(sw[:, None] * A, sw * l, rcond=None)[0]
        wres = sw * (l - A @ xn)                       # sweep residual (factor inside the whitener)
        vw = wres / np.sqrt(np.maximum(fac, 1e-12))     # lawWhitenedResiduals: factor divided out
        f, s = robust_weights(vw, infl)
        if x is not None and np.max(np.abs(xn - x)) < tol * max(np.abs(xn).max(), 1):
            x = xn; break
        x = xn; fac = f
    vtv = np.sum(wres**2); s0 = np.sqrt(vtv / (m - u))
    alive = np.sum(fac > 0)
    return x, s0, np.sqrt(vtv / (alive - u)), fac

def _main():
    print('(a) sigma0 of the reweighted problem; true sigma = 1, u = 6, m = 2000, blunders +-25')
    rs = np.random.RandomState(3); m, u = 2000, 6
    A = rs.randn(m, u); xt = rs.randn(u)
    for frac in (0.0, 0.05, 0.10, 0.20):
        for infl in ('huber', 'igg3', 'bisquare'):
            e = rs.randn(m); nb = int(frac * m); idx = rs.choice(m, nb, replace=False)
            e[idx] += 25 * np.sign(rs.randn(nb))
            l = A @ xt + e
            x, s0, s0alive, fac = irls(A, l, infl)
            cl = np.setdiff1d(np.arange(m), idx)
            print('  contam %3.0f%% %-8s: R.sigma0 = %.3f  (dof m-u); with dof = alive-u: %.3f; zero rows %4d; clean-row rms %.3f'
                  % (100 * frac, infl, s0, s0alive, np.sum(fac == 0), np.std(e[cl])))

    print('\n(b) stat kernels')
    def tauCrit(alpha, r):
        t = stats.t.ppf(1 - alpha / 2, r - 1); return t * np.sqrt(r) / np.sqrt(r - 1 + t**2)
    # independent check of Pope's tau: |tau| = sqrt(r) |T|/sqrt(r-1+T^2), T ~ t(r-1);
    # Monte-Carlo tau from studentised residuals of a real LS problem
    rs = np.random.RandomState(7); r = 12; mm, uu = 20, 8
    Am = rs.randn(mm, uu); H = Am @ np.linalg.pinv(Am); Qv = np.eye(mm) - H
    taus = []
    for _ in range(40000):
        v = Qv @ rs.randn(mm); s = np.sqrt(v @ v / (mm - uu)); taus.append(v[0] / (s * np.sqrt(Qv[0, 0])))
    print('  tauCrit(0.05, r=12) code = %.5f ; Monte-Carlo 97.5%% quantile of tau = %.5f; bound sqrt(r)=%.4f'
          % (tauCrit(0.05, r), np.quantile(np.abs(taus), 0.95), np.sqrt(r)))
    def delta0(a, b):
        k = stats.norm.ppf(1 - a / 2); from scipy.optimize import brentq
        return brentq(lambda d: 1 - stats.norm.cdf(k - d) + stats.norm.cdf(-k - d) - (1 - b), 0, 50)
    print('  delta0(0.001, 0.20) = %.4f (Baarda textbook 4.13); delta0(0.05,0.20)=%.4f (2.80)' % (delta0(1e-3, .2), delta0(.05, .2)))
    print('  chi2inv(0.95,2) = %.4f -> ellipse factor sqrt = %.4f (2.4477)' % (stats.chi2.ppf(.95, 2), np.sqrt(stats.chi2.ppf(.95, 2))))
    # fInv as coded: betaincinv(p, d1/2, d2/2) -> x = d2/d1 * xb/(1-xb)
    for p, d1, d2 in ((0.95, 3, 40), (0.99, 1, 1000)):
        xb = special.betaincinv(d1 / 2, d2 / 2, p); print('  fInv(%.2f,%d,%d) code %.5f scipy %.5f' % (p, d1, d2, d2 / d1 * xb / (1 - xb), stats.f.ppf(p, d1, d2)))

    print('\n(c) Results.compare dof: n = 1000 rows, nc = 6 commons, nl = 200 locals; the full model has one extra COMMON')
    rs = np.random.RandomState(11); n, nc, nl = 1000, 6, 200
    Ac = rs.randn(n, nc + 1); Bl = np.zeros((n, nl))
    for i in range(n): Bl[i, i % nl] = 1.0
    null_p = []; null_true = []
    for trial in range(2000):
        l = Ac[:, :nc] @ rs.randn(nc) + Bl @ rs.randn(nl) + rs.randn(n)      # extra common = 0 (H0 true)
        Mf = np.column_stack([Ac, Bl]); Mr = np.column_stack([Ac[:, :nc], Bl])
        vf = l - Mf @ np.linalg.lstsq(Mf, l, rcond=None)[0]; vr = l - Mr @ np.linalg.lstsq(Mr, l, rcond=None)[0]
        vPvF, vPvR = vf @ vf, vr @ vr
        uF, uR = nc + 1, nc                          # numel(R.x): commons only
        df1, df2c = uF - uR, n - uF; Fc = (vPvR - vPvF) / df1 / (vPvF / df2c)
        df2t = n - (nc + 1 + nl); Ft = (vPvR - vPvF) / df1 / (vPvF / df2t)
        null_p.append(Fc > stats.f.ppf(0.95, df1, df2c)); null_true.append(Ft > stats.f.ppf(0.95, df1, df2t))
    print('  rejection rate under H0 at alpha = 0.05: code (df2 = n - nc) %.3f ; correct (df2 = n - nc - nl) %.3f'
          % (np.mean(null_p), np.mean(null_true)))

if __name__=='__main__': _main()
