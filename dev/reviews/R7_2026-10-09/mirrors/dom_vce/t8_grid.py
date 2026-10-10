import numpy as np
from loop import *
print('Same data, shipped code path; only Noise.Grid differs (60 = default, 256 = arc length, Ellmer design)')
for phi, nArc, seed in [(0.9, 40, 11), (0.9, 40, 12), (0.7, 9, 11), (0.95, 40, 11)]:
    acv = phi**np.arange(1024)/(1-phi**2)
    arcs = make_data(acv, nArc, 256, seed=seed)
    for L in (60, 256):
        H = run(arcs, L, passes=10)
        print(f'phi={phi} {nArc}x256 seed={seed} Grid={L}:')
        print('   maxDev :', ' '.join('%.1e' % h['maxDev'] for h in H))
        print('   level  :', ' '.join('%.4g' % h['levelNew'] for h in H), '  (truth %.3f)' % (1/(1-phi**2)))
        print('   rho1   :', ' '.join('%.3f' % (h['C'][1]/h['C'][0]) for h in H), ' sigma0 last %.3f' % H[-1]['s0'])
