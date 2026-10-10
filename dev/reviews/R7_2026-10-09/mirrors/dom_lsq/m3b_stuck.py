"""After a rejection the 'damped' re-solve is identical (nl = 0: Reg 'lm' never reaches
the common block), so the same step is rejected again until MaxIter / 1e12."""
import numpy as np
from m3_gauss_newton import adjust_nonlinear
# Rosenbrock from several starts; fit's default budget MaxIter = MaxPasses = 10
def model2(x):
    r = np.array([10*(x[1]-x[0]**2), 1-x[0]])
    J = np.array([[20*x[0], -10], [1, 0]])     # J = -dr/dx  (w = l - f, J = df/dx)
    return J, r
for x0 in ([-1.2, 1.0], [2.0, -1.0], [0.5, 2.0], [-0.5, -0.5], [3.0, 3.0]):
    for mi in (10, 20, 50):
        x, why, h = adjust_nonlinear(model2, np.array(x0), MaxIter=mi)
        rej = sum(1 for e in h if not e[2]); 
        print('start %-12s MaxIter %2d: reason %-8s passes %2d rejected %2d  x=%s  |x-x*|=%.2e' % (x0, mi, why, len(h), rej, np.round(x,6), np.linalg.norm(x-np.array([1,1]))))
# a non-zero-residual exponential fit where GN overshoots from a far start
t = np.linspace(0, 4, 30); y = 3*np.exp(-1.3*t) + 0.05*np.sin(7*t)
def m(x):
    f = x[0]*np.exp(x[1]*t); J = np.column_stack([np.exp(x[1]*t), x[0]*t*np.exp(x[1]*t)]); return J, y - f
for x0 in ([1.0, 0.5], [10.0, -5.0], [0.1, 1.0]):
    x, why, h = adjust_nonlinear(m, np.array(x0), MaxIter=10)
    print('expfit start %-12s: reason %-8s passes %2d rejected %2d x=%s' % (x0, why, len(h), sum(1 for e in h if not e[2]), np.round(x,4)))
    print('   vPv trajectory', ['%.3g%s' % (e[1], '' if e[2] else 'R') for e in h])
