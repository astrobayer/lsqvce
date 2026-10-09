"""Root-cause counterfactuals for the closed-loop runaway (shipped: lawSpecs eps floor,
covModel psdFloor 1e-10).  Variants change ONE thing in the declaration:
  shipped : S = lawSpecs(C) floored at eps*max, then covModel (floor 1e-10*max)
  floor1e3: the same series floored at 1e-3*max instead
Data/loop otherwise identical (Grid 60, 256-epoch arcs, exact traces)."""
import numpy as np, loop, kern
from scipy.linalg import toeplitz
orig = kern.lawSpecs
def lawSpecs_floor(C, fl):
    fN, S = orig(C)
    L = len(C); Cn = C/max(C[0], np.finfo(float).eps)
    Sraw = 2*(Cn[0] + 2*(np.cos(2*np.pi*np.outer(fN, np.arange(1, L))) @ Cn[1:]))
    return fN, np.maximum(Sraw, fl*Sraw.max())
for name, fl in [('shipped', None), ('floor1e-3', 1e-3)]:
    loop.lawSpecs = orig if fl is None else (lambda C, fl=fl: lawSpecs_floor(C, fl))
    for phi, nArc in [(0.8, 20), (0.9, 40), (0.7, 9)]:
        acv = phi**np.arange(1024)/(1-phi**2)
        res = []
        for seed in (10, 11, 12):
            H = loop.run(loop.make_data(acv, nArc, 256, seed=seed), 60, passes=10)
            res.append(max(h['levelNew'] for h in H)/(1/(1-phi**2)))
        print(f'{name:9s} phi={phi}: max over passes of level/truth, 3 seeds: {np.round(res, 2)}')
