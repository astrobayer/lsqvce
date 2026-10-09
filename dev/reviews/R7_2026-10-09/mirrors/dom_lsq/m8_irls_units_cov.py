"""(a) adjustRobust stop rule  max|dx| < Tol * max(max|x|, 1)  (Tol 1e-8) is not unit-invariant:
    the same problem with x expressed in 1e-9 units (A*1e9) stops at pass 2.
(b) the reported covariance of the IRLS result  s0^2 (A'FA)^-1  vs the Monte-Carlo sampling
    covariance of the estimator, clean and contaminated."""
import numpy as np, contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    from m5_irls_stats import robust_weights

def irls(A, l, infl, maxit=30, tol=1e-8):
    m, u = A.shape; fac = np.ones(m); x = None
    for it in range(1, maxit + 1):
        sw = np.sqrt(np.maximum(fac, 1e-12))
        xn = np.linalg.lstsq(sw[:, None] * A, sw * l, rcond=None)[0]
        wres = sw * (l - A @ xn)
        vw = wres / np.sqrt(np.maximum(fac, 1e-12))
        f, s = robust_weights(vw, infl)
        dx = np.inf if x is None else np.max(np.abs(xn - x))
        x = xn; facUsed = fac
        if it > 1 and dx < tol * max(np.abs(x).max(), 1):
            break
        fac = f
    p = np.maximum(facUsed, 1e-12); N = (A * p[:, None]).T @ A
    s0sq = np.sum(p * (l - A @ x)**2) / (m - u)
    return x, it, s0sq * np.linalg.inv(N)

rs = np.random.RandomState(8); m, u = 400, 4
A = rs.randn(m, u); xt = rs.randn(u)
e = rs.randn(m); idx = rs.choice(m, 40, replace=False); e[idx] += 8 * np.sign(rs.randn(40))
l = A @ xt + e
xa, ita, _ = irls(A, l, 'igg3', maxit=200, tol=1e-14)        # the fixed point
for sc in (1.0, 1e-3, 1e-6, 1e-9):
    x, it, _ = irls(A / sc, l, 'igg3')                         # x in units of sc  -> x_true*sc
    print('(a) parameter unit %-6g: passes %2d, |x - fixedpoint| / sigma-ish = %.2e' % (sc, it, np.max(np.abs(x / sc - xa)) / 0.05))

print('\n(b) reported se vs Monte-Carlo sd of the IRLS estimate (m = %d, u = %d, N(0,1) noise)' % (m, u))
for frac, infl in ((0.0, 'huber'), (0.0, 'igg3'), (0.10, 'huber'), (0.10, 'igg3'), (0.10, 'bisquare')):
    xs = []; ses = []
    for t in range(1500):
        e = rs.randn(m); nb = int(frac * m)
        if nb: ii = rs.choice(m, nb, replace=False); e[ii] += 25 * np.sign(rs.randn(nb))
        x, it, C = irls(A, A @ xt + e, infl)
        xs.append(x); ses.append(np.sqrt(np.diag(C)))
    sd = np.std(np.array(xs), axis=0); se = np.mean(np.array(ses), axis=0)
    print('   contam %3.0f%% %-8s: mean reported se / MC sd = %s' % (100 * frac, infl, np.round(se / sd, 3)))
