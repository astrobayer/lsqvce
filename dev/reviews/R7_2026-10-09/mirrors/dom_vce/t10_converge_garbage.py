import numpy as np
from loop import *
from scipy.linalg import toeplitz
found = 0
for phi, nArc in [(0.9, 40), (0.7, 9), (0.8, 20), (0.95, 40)]:
    acv = phi**np.arange(1024)/(1-phi**2)
    for seed in range(10, 22):
        arcs = make_data(acv, nArc, 256, seed=seed)
        try:
            H = run(arcs, 60, passes=10)          # the shipped default budget MaxPasses = 10
        except Exception as e:
            print(f'phi={phi} seed={seed}: mirror chol failed ({e})'); continue
        stop = next((h for h in H if h['maxDev'] <= 1e-2), None)
        h = stop or H[-1]
        C = h['C']; ev = np.linalg.eigvalsh(toeplitz(C/C[0])).min()
        tag = 'CONVERGED' if stop else 'maxPasses'
        print(f'phi={phi} {nArc}x256 seed={seed}: stop={tag:9s} at pass {h["it"]:2d}  maxDev {h["maxDev"]:.2e}  level {h["level"]:.4g} (truth {1/(1-phi**2):.3f})  '
              f'rho1 {C[1]/C[0]:.3f} (truth {phi})  sigma0 {h["s0"]:.3f}  min eig Toeplitz(C/C0) {ev:.3f}')
