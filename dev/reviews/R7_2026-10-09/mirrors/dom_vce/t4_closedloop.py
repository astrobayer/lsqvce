import numpy as np
from loop import *
np.set_printoptions(linewidth=160, precision=4, suppress=True)
for phi, L, m, nArc in [(0.5, 60, 256, 9), (0.9, 60, 256, 40), (0.95, 60, 256, 40), (0.9, 16, 256, 40)]:
    acv = phi**np.arange(4*m)/(1-phi**2)
    arcs = make_data(acv, nArc, m, seed=11)
    H = run(arcs, L, passes=14)
    print(f'--- AR(1) phi={phi}, Grid L={L}, {nArc} arcs x {m}; exact traces, shipped path')
    for h in H:
        print('  pass %2d  maxDev %.3e  |a2_0-1| %.3e  max_j>=1|a2_j-1| %.3e   level %.4f -> %.4f  sigma0 %.4f  rho1(C) %.4f' %
              (h['it'], h['maxDev'], h['dev0'], h['dev1'], h['level'], h['levelNew'], h['s0'], h['C'][1]/h['C'][0]))
    h = H[-1]; Cn = h['C']/h['C'][0]; rn = h['r']/h['r'][0]; tn = phi**np.arange(L)
    lg = [1, 2, 5, 10, 20, 30, L-2, L-1]
    lg = [g for g in lg if g < L]
    print('  lags', lg); print('  reported C/C0 ', Cn[lg]); print('  whitener r/r0 ', rn[lg]); print('  truth         ', tn[lg])
