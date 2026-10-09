"""Mirror of lsqvce.engine.adjustNonlinear (step control) on top of the accumulator
arithmetic for a local-free system (nl = 0) with a CommonPrior:

  accumulate(x, lam):  rows [J(x), w = l - f(x)] whitened by sqrt(p);
                       nl = 0 -> 'lm' Reg acts on the EMPTY local block (no damping of commons!)
                       solve(): no prior -> plain LS step;  prior -> N = JtJ + K'K/s^2,
                       n = Jtw + K'x0/s^2  (CommonPrior.normals), dx = N \\ n
                       F = ltl = w'w (data only, commonBlocks)
  loop exactly as adjustNonlinear.m:93-142.
"""
import numpy as np

def accumulate(model, x, lam, prior):
    J, w = model(x)
    N = J.T @ J; n = J.T @ w; F = w @ w
    # BlockAccumulator 'Reg','lm' -> applyRegularisation on Nll: nl = 0 -> preeliminate
    # returns BEFORE applyRegularisation (G-07 branch). So lam has NO effect.
    if prior is not None:
        K, x0, s = prior
        N = N + K.T @ K / s**2; n = n + K.T @ x0 / s**2
    dx = np.linalg.solve(N, n)
    return dx, F, N

def adjust_nonlinear(model, x0, MaxIter=20, TolX=1e-8, TolF=1e-10, Lambda0=1e-6, prior=None, trace=False):
    x = x0.copy(); lam = Lambda0
    dx, vPv, N = accumulate(model, x, lam, prior)
    hist = []; reason = 'maxIter'
    for it in range(1, MaxIter + 1):
        xNew = x + dx
        dxNew, vPvNew, NNew = accumulate(model, xNew, lam, prior)
        if Lambda0 > 0 and it > 1 and vPvNew > vPv:
            lam *= 10; hist.append((it, vPvNew, False, lam))
            if lam > 1e12: reason = 'diverged'; break
            dx, _, N = accumulate(model, x, lam, prior); continue
        sig = np.sqrt(np.diag(np.linalg.inv(N)))
        dxRel = np.max(np.abs(dx) / sig)
        x = xNew; dx, N = dxNew, NNew
        lam = max(lam / 10, 1e-12); hist.append((it, vPvNew, True, lam))
        relF = abs(vPvNew - vPv) / max(abs(vPvNew), 1e-300)
        if dxRel < TolX: reason = 'tolX'; break
        if it > 1 and relF < TolF: reason = 'tolF'; break
        vPv = vPvNew
    # final: undamped solve at x, reported estimate x + step
    dxf, _, _ = accumulate(model, x, 0, prior)
    return x + dxf, reason, hist

# ---- 1. the test_fit exponential-decay fixture, with and without a common prior ----
rs = np.random.RandomState(5)
t = np.linspace(0, 1, 40); xt = np.array([0.7, -0.2])
lobs = np.exp(xt[0]*t) + xt[1] + 1e-3*rs.randn(40)
def expo(x):
    J = np.column_stack([t*np.exp(x[0]*t), np.ones(40)])
    return J, lobs - (np.exp(x[0]*t) + x[1])
p = 1/1e-3**2; sp = np.sqrt(p)
def expoW(x):
    J, w = expo(x); return sp*J, sp*w

xh, why, h = adjust_nonlinear(expoW, np.array([0.5, 0.0]))
print('no prior          : x =', xh, why, 'passes', len(h))

# a declared prior on x(2): b ~ 0.30 with sigma 1e-3  (L = [0 1], x0 = 0.30)
K = np.array([[0., 1.]]); x0p = np.array([0.30]); sg = 1e-3
xh2, why2, h2 = adjust_nonlinear(expoW, np.array([0.5, 0.0]), prior=(K, x0p, sg))
# the MAP / regularised LS solution, by dense minimisation of the declared objective
from scipy.optimize import least_squares
res = least_squares(lambda x: np.concatenate([expoW(x)[1], (K@x - x0p)/sg]), np.array([0.5, 0.0]), xtol=1e-15, ftol=1e-15, gtol=1e-15)
print('prior b~0.30(1e-3): x =', xh2, why2, 'passes', len(h2), '| accepted:', [a for (_,_,a,_) in h2].count(True))
print('   declared-objective optimum (scipy) x =', res.x)

# a datum-like prior with x0 = 0 on b: b ~ 0 sigma 1e-3
xh3, why3, h3 = adjust_nonlinear(expoW, np.array([0.5, 0.0]), prior=(K, np.array([0.0]), sg))
res3 = least_squares(lambda x: np.concatenate([expoW(x)[1], (K@x - 0)/sg]), np.array([0.5, 0.0]), xtol=1e-15, ftol=1e-15, gtol=1e-15)
print('prior b~0   (1e-3): x =', xh3, why3, 'passes', len(h3))
print('   declared-objective optimum (scipy) x =', res3.x)

# ---- 2. descent check on the FIRST step (it > 1 guard) and Lambda0 = 0 ----
# Rosenbrock-like residuals where the pure GN first step overshoots
def model2(x):
    r = np.array([10*(x[1]-x[0]**2), 1-x[0], 0.0])
    J = -np.array([[-20*x[0], 10], [-1, 0], [0, 0]])
    return np.column_stack([J[:2].T]).T[:2] if False else (J[:2], r[:2])
# F increase on the first step: pick a start where the GN step overshoots
def F(x): return np.sum(model2(x)[1]**2)
x0 = np.array([-1.2, 1.0])
J, w = model2(x0); dx1 = np.linalg.lstsq(J, w, rcond=None)[0]
print('\nRosenbrock start', x0, 'F(x0)=%.4g  F(x0+dx_GN)=%.4g  (first step raises F: %s)' % (F(x0), F(x0+dx1), F(x0+dx1) > F(x0)))
xr, whyr, hr = adjust_nonlinear(model2, x0, MaxIter=50)
print('  default Lambda0=1e-6: first-pass record', hr[0], '-> accepted although F rose' if hr[0][2] and hr[0][1] > F(x0) else '')
xr0, whyr0, hr0 = adjust_nonlinear(model2, x0, MaxIter=50, Lambda0=0)
ups = sum(1 for k in range(1, len(hr0)) if hr0[k][1] > hr0[k-1][1])
print('  Lambda0=0: %d passes, %d accepted F-increases, reason %s' % (len(hr0), ups, whyr0))

# ---- 3. lam has no effect on a local-free system ----
J, w = expoW(np.array([0.5, 0.0]))
print('\nnl = 0: step with lam=1e-6 vs lam=1e6 identical? ', np.allclose(accumulate(expoW, np.array([0.5, 0.]), 1e-6, None)[0],
                                                             accumulate(expoW, np.array([0.5, 0.]), 1e6, None)[0]),
      '(BlockAccumulator.preeliminate nl==0 returns before applyRegularisation)')
