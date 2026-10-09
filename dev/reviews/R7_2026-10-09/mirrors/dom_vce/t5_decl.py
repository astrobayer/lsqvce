import numpy as np
from kern import *
np.set_printoptions(linewidth=160, precision=4, suppress=True)
print('Declaration path C -> lawSpecs -> covModel(toeplitz) -> r : absolute error on lags 0..L-1 and leakage beyond')
for (L, m, phi) in [(60,256,0.5),(60,256,0.9),(60,256,0.95),(16,256,0.9),(60,2160,0.95),(2160,2160,0.95),(256,256,0.9)]:
    C = phi**np.arange(L)
    fN, S = lawSpecs(C)
    r = covModelToeplitzACV(fN, S, 5.0, m)
    k = min(L, m)
    err = np.abs(r[:k]-C[:k])
    beyond = np.abs(r[L:]).max() if m > L else 0
    print(f'L={L:5d} m={m:5d} phi={phi}: max|r_n-C_n| (n<L) = {err.max():.3e} at n={err.argmax()}, |r_0-1|={abs(r[0]-1):.2e}, max|r_n| n>=L = {beyond:.3e} (truth {phi**L:.2e})')
