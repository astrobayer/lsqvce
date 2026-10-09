"""adjustNonlinear at a rejection loop: what Results.x is (x_last_accepted + undamped step)."""
import numpy as np, contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    from m3_gauss_newton import adjust_nonlinear
from scipy.optimize import least_squares
t = np.linspace(0, 4, 30); y = 3*np.exp(-1.3*t) + 0.05*np.sin(7*t)
def m(x):
    f = x[0]*np.exp(x[1]*t); J = np.column_stack([np.exp(x[1]*t), x[0]*t*np.exp(x[1]*t)]); return J, y - f
F = lambda x: np.sum(m(x)[1]**2)
for x0 in (np.array([10.0, -5.0]), np.array([0.1, 1.0])):
    xr, why, h = adjust_nonlinear(m, x0, MaxIter=10)
    acc = [e for e in h if e[2]]
    print('start %s: reason %s, %d rejections of %d passes; best accepted F=%.3f; REPORTED x=%s F(reported)=%.3f'
          % (x0, why, sum(not e[2] for e in h), len(h), acc[-1][1], np.round(xr, 4), F(xr)))
r = least_squares(lambda x: m(x)[1], np.array([1.0, 0.5]), method='lm')
print('optimum', np.round(r.x, 4), 'F=%.4f' % F(r.x))
# with MaxIter = 30: 1e12 needs log10(1e12/1e-6)=18 rejections -> 'diverged'
xr, why, h = adjust_nonlinear(m, np.array([10.0, -5.0]), MaxIter=30)
print('MaxIter 30: reason', why, 'passes', len(h), 'reported F=%.3f' % F(xr))
