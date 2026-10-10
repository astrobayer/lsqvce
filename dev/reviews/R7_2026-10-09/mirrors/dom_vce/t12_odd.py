import numpy as np
from loop import *
# colour-fixture geometry (9 arcs, AR(1) 0.5, Grid 60, shipped path), EVEN vs ODD arc length
for m in (256, 255):
    acv = 0.5**np.arange(1024)/(1-0.25)
    arcs = make_data(acv, 9, m, seed=11)
    try:
        H = run(arcs, 60, passes=8)
    except Exception as e:
        print(m, 'failed', e); continue
    xs = np.array([h['x'] for h in H]); xt = np.array([.8, -.45, .3, .15])
    print(f'm={m}: maxDev', ' '.join('%.1e' % h['maxDev'] for h in H))
    print('       level ', ' '.join('%.3g' % h['levelNew'] for h in H), ' truth 1.333; sigma0 last %.3f' % H[-1]['s0'])
    print('       rho1  ', ' '.join('%.3f' % (h['C'][1]/h['C'][0]) for h in H))
    print('       |x - xtrue| last pass', np.round(np.abs(xs[-1]-xt), 4))
