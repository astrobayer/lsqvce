import numpy as np
from kern import *
np.set_printoptions(linewidth=150, precision=4, suppress=True)
def ar1(phi, n): return phi**np.arange(n)
for (L, m, phi) in [(60, 256, 0.5), (60, 256, 0.9), (60, 100, 0.9), (60, 2160, 0.95), (2160, 2160, 0.95), (16, 256, 0.9)]:
    C = ar1(phi, L)
    fN, S = lawSpecs(C)
    r = covModelToeplitzACV(fN, S, 5.0, m)
    lags = [0, 1, 2, L//4, L//2, L-2, L-1]
    print(f'L={L:4d} m={m:4d} phi={phi}:  r_n/C_n at lags {lags}:', np.round(r[lags]/C[lags], 4))
    if m > L:
        k = np.arange(L, min(m, 2*L+5))
        print(f'      lags >= L (C says 0 beyond L-1; true AR {phi}^n): r[L..L+4] =', np.round(r[L:L+5], 4),
              ' true =', np.round(ar1(phi, L+5)[L:], 4), ' r[2L-3]=%.4f vs r[1]=%.4f' % (r[2*L-3], r[1]) if m > 2*L-3 else '')
