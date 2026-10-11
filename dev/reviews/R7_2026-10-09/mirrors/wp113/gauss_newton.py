"""WP-113 mirror: Gauss-Newton with the damping and the prior where they belong
(R7-14, R7-04). The new adjustNonlinear, step for step:

  accumulate(x)   rows [J(x), w = l - f(x)] (whitened); F(x) = w'w + |(K x - x0)/s|^2
                  (the DECLARED objective, prior included); d = diag(J'J)
  step(x, lam)    one least-squares solve of the stacked rows
                    [J ; K/s ; sqrt(lam) diag(sqrt d)] dx ~ [w ; (x0 - K x)/s ; 0]
                  -- the prior on x as K dx ~ x0 - K x (R7-04), Marquardt damping on
                  the COMMON block as pseudo-observations dx ~ 0 (R7-14)
  loop            descent check from pass 1 (F(x+dx) > F(x) -> rejected, lam x10, the
                  step re-solved at the SAME accumulation); accepted -> lam / 10
                  (below 1e-12 -> 0); lam = 0 rejected -> 1e-3
  result          the best ACCEPTED point x (never x + a rejected step)
Compared with the shipped loop (the review's m3 mirror) and with scipy's optimum.
"""
import numpy as np
from scipy.optimize import least_squares


def new_gn(model, x0, MaxIter=10, TolX=1e-8, TolF=1e-10, Lambda0=1e-6, prior=None):
    def acc(x):
        J, w = model(x)
        F = w @ w
        if prior is not None:
            K, t, s = prior
            F += np.sum(((K @ x - t) / s) ** 2)
        return J, w, F, np.sum(J ** 2, axis=0)

    def step(J, w, x, lam, d):
        rows, rhs = [J], [w]
        if prior is not None:
            K, t, s = prior
            rows.append(K / s); rhs.append((t - K @ x) / s)
        if lam > 0:
            rows.append(np.sqrt(lam) * np.diag(np.sqrt(np.maximum(d, 1e-300)))); rhs.append(np.zeros(len(x)))
        M = np.vstack(rows); r = np.concatenate(rhs)
        dx = np.linalg.lstsq(M, r, rcond=None)[0]
        Q = np.linalg.inv(M.T @ M)
        return dx, np.sqrt(np.diag(Q))

    x = np.array(x0, float); lam = Lambda0
    J, w, F, d = acc(x)
    dx, sig = step(J, w, x, lam, d)
    hist = []; reason = 'maxIter'
    for it in range(1, MaxIter + 1):
        xN = x + dx
        JN, wN, FN, dN = acc(xN)
        if FN > F:
            lam = 1e-3 if lam == 0 else lam * 10
            hist.append((it, FN, False))
            if lam > 1e12: reason = 'diverged'; break
            dx, sig = step(J, w, x, lam, d)
            continue
        dxRel = np.max(np.abs(dx) / sig)
        Fold = F
        x, J, w, F, d = xN, JN, wN, FN, dN
        lam = lam / 10 if lam / 10 >= 1e-12 else 0.0
        hist.append((it, F, True))
        dx, sig = step(J, w, x, lam, d)
        if dxRel < TolX: reason = 'tolX'; break
        if it > 1 and abs(F - Fold) / max(F, 1e-300) < TolF: reason = 'tolF'; break
    return x, reason, hist


def optimum(model, x0, prior=None):
    def res(x):
        r = model(x)[1]
        if prior is not None:
            K, t, s = prior
            r = np.concatenate([r, (K @ x - t) / s])
        return r
    return least_squares(res, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15).x


def Fobj(model, x, prior=None):
    F = model(x)[1] @ model(x)[1]
    if prior is not None:
        K, t, s = prior
        F += np.sum(((K @ x - t) / s) ** 2)
    return F


# exponential, non-zero residual (review m3b): y = a exp(b t)
t = np.linspace(0, 4, 30); y = 3 * np.exp(-1.3 * t) + 0.05 * np.sin(7 * t)
def expfit(x):
    f = x[0] * np.exp(x[1] * t)
    return np.column_stack([np.exp(x[1] * t), x[0] * t * np.exp(x[1] * t)]), y - f

def rosen(x):
    return np.array([[20 * x[0], -10], [1, 0]]), np.array([10 * (x[1] - x[0] ** 2), 1 - x[0]])

# the test_fit exponential-decay fixture (review m3d)
rs = np.random.RandomState(5)
tt = np.linspace(0, 1, 40); xt = np.array([0.7, -0.2])
lobs = np.exp(xt[0] * tt) + xt[1] + 1e-3 * rs.randn(40)
sp = 1 / 1e-3
def decay(x):
    J = np.column_stack([tt * np.exp(x[0] * tt), np.ones(40)])
    return sp * J, sp * (lobs - (np.exp(x[0] * tt) + x[1]))


if __name__ == '__main__':
    print('== R7-14: exponential fit, budget MaxIter 10 / 30')
    xs = optimum(expfit, np.array([1.0, -1.0]))
    for x0 in ([1.0, 0.5], [10.0, -5.0], [0.1, 1.0]):
        for mi in (10, 30):
            x, why, h = new_gn(expfit, x0, MaxIter=mi)
            print('start %-12s MaxIter %2d: %-8s passes %2d rejected %2d  F %.6g (optimum %.6g)  |x-x*| %.1e'
                  % (x0, mi, why, len(h), sum(1 for e in h if not e[2]), Fobj(expfit, x), Fobj(expfit, xs), np.linalg.norm(x - xs)))
    print('\n== R7-14: Rosenbrock')
    for x0 in ([-1.2, 1.0], [2.0, -1.0], [0.5, 2.0], [-0.5, -0.5], [3.0, 3.0]):
        x, why, h = new_gn(rosen, x0, MaxIter=50)
        ups = sum(1 for e in h if e[2] and e[1] > Fobj(rosen, np.array(x0)) and e is h[0])
        print('start %-12s: %-8s passes %2d rejected %2d  |x-(1,1)| %.1e  first accepted step uphill: %s'
              % (x0, why, len(h), sum(1 for e in h if not e[2]), np.linalg.norm(x - 1), bool(ups)))
    print('\n== R7-04: prior on b of the decay fixture')
    K = np.array([[0., 1.]])
    for tp, sg in ((0.30, 1e-6), (0.30, 1e-3), (0.0, 1e-6), (0.0, 1e-3), (-0.2, 1e-6)):
        pr = (K, np.array([tp]), sg)
        x, why, h = new_gn(decay, [0.5, 0.0], MaxIter=20, prior=pr)
        xo = optimum(decay, np.array([0.5, 0.0]), pr)
        print('b ~ %+.2f (sigma %g): GN b = %+.6f a = %.6f (%s, %d)  optimum b = %+.6f a = %.6f  |dx| %.1e'
              % (tp, sg, x[1], x[0], why, len(h), xo[1], xo[0], np.linalg.norm(x - xo)))
