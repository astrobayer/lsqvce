import numpy as np
from loop2 import *
for phi, L, m, nArc, seed in [(0.9, 60, 256, 40, 11), (0.9, 60, 256, 40, 12), (0.7, 60, 256, 9, 11), (0.95, 60, 256, 40, 11), (0.5, 60, 256, 9, 11)]:
    acv = phi**np.arange(4*m)/(1-phi**2)
    arcs = make_data(acv, nArc, m, seed=seed)
    for ex in (False, True):
        try:
            H = run(arcs, L, passes=12, exact_whitener=ex)
        except Exception as e:
            print('  FAILED', e); continue
        print(f'phi={phi} seed={seed} {nArc}x{m} whitener={"fine-grid (r=C, no wrap)" if ex else "SHIPPED interp+psd2acv"}')
        print('   maxDev :', ' '.join('%.1e' % h['maxDev'] for h in H))
        print('   level  :', ' '.join('%.3g' % h['levelNew'] for h in H), '  (truth %.3f)' % (1/(1-phi**2)))
        print('   rho1   :', ' '.join('%.3f' % (h['C'][1]/h['C'][0]) for h in H))
