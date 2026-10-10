import numpy as np, contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    from m3_gauss_newton import adjust_nonlinear
t = np.linspace(0, 4, 30); y = 3*np.exp(-1.3*t) + 0.05*np.sin(7*t)
def m(x):
    f = x[0]*np.exp(x[1]*t); J = np.column_stack([np.exp(x[1]*t), x[0]*t*np.exp(x[1]*t)]); return J, y - f
F = lambda x: np.sum(m(x)[1]**2)
for L0 in (1e-6, 0.0):
    for x0 in ([10.0, -5.0], [0.1, 1.0], [1.0, 0.5]):
        try:
            x, why, h = adjust_nonlinear(m, np.array(x0), MaxIter=10, Lambda0=L0)
            ups = sum(1 for k in range(1, len(h)) if h[k][2] and h[k][1] > h[k-1][1])
            print('Lambda0=%g start %-12s: %-8s passes %2d rejected %d accepted-uphill %d  F(reported)=%.4f'
                  % (L0, x0, why, len(h), sum(not e[2] for e in h), ups, F(x)))
        except Exception as e:
            print('Lambda0=%g start %s: %s' % (L0, x0, e))
