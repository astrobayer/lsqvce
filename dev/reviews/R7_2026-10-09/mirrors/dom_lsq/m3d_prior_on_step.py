"""Nonlinear + data.priors: the CommonPrior K x ~ x0 is applied to the Gauss-Newton
INCREMENT dx (the accumulator solves for dx; CommonPrior.normals adds K'x0/s^2 to n),
not to x. Fixed point and reported value vs the declared objective."""
import numpy as np, contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    from m3_gauss_newton import adjust_nonlinear, expoW
from scipy.optimize import least_squares
K = np.array([[0., 1.]])
for x0p, sg in ((0.30, 1e-6), (0.30, 1e-3), (0.0, 1e-6), (-0.2, 1e-6)):
    xh, why, h = adjust_nonlinear(expoW, np.array([0.5, 0.0]), prior=(K, np.array([x0p]), sg))
    r = least_squares(lambda x: np.concatenate([expoW(x)[1], (K @ x - x0p) / sg]), np.array([0.5, 0.0]),
                      xtol=1e-15, ftol=1e-15, gtol=1e-15)
    print('prior b ~ %+.2f (sigma %g): GN reports b = %+.6f (%s, %d passes); declared optimum b = %+.6f; a: %.5f vs %.5f'
          % (x0p, sg, xh[1], why, len(h), r.x[1], xh[0], r.x[0]))
