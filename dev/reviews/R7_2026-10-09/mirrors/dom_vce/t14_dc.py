import numpy as np
from loop import *
from kern import dctI
# Start AT the truth (C0 = true ACV shape), phi = 0.5 (declaration path accurate there), 9 x 256, Grid 60
phi = 0.5; acv = phi**np.arange(1024)/(1-phi**2)
arcs = make_data(acv, 9, 256, seed=11)
C0 = acv[:60]/acv[0]
for zr in (True, False):
    import loop
    H = run(arcs, 60, passes=25, C0=C0, zero_rule=zr)
    print('ZeroFreqRule %s: |a2_0-1| per pass:' % ('lowest' if zr else 'keep  '), ' '.join('%.3f' % h['dev0'] for h in H[::3]),
          '| max_j>=1:', ' '.join('%.3f' % h['dev1'] for h in H[::3]))
St = dctI(C0); print('true DCT amplitudes S_0, S_1 = %.4f %.4f (ratio %.3f) -> the lowest rule forces S_0 = S_1' % (St[0], St[1], St[0]/St[1]))
