"""Counterfactual: the same closed loop with a WIENER-KHINCHIN-consistent DCT-I pair
(forward F: S_j = C_0 + 2 sum_{n=1}^{N-2} C_n cos(pi n j/(N-1)) + C_{N-1}(-1)^j;
inverse F^-1: C_n = [S_0 + 2 sum S_j cos + S_{N-1}(-1)^n] / (2(N-1)) -- positive S =>
C is the ACV of a periodic process => Toeplitz(C) PSD) in place of the 'ortho' X."""
import numpy as np, loop, kern
from scipy.linalg import toeplitz
N = 60
n = np.arange(N); F = 2*np.cos(np.pi*np.outer(n, n)/(N-1)); F[:, 0] = 1; F[:, -1] = (-1.0)**n
Fi = np.linalg.inv(F)
rng = np.random.default_rng(0)
print('check: positive S -> min eig Toeplitz(F^-1 S) over 2000 random spectra: %.2e' %
      min(np.linalg.eigvalsh(toeplitz(Fi@np.exp(rng.normal(size=N)))).min() for _ in range(2000)))
def psdVCE_wk(omT, sT, C, zero_rule=True, floor=1e-12):
    Om = (F@omT).sum(1); s = (F@sT).sum(1)
    a2 = np.where(np.abs(s) > 0, Om/s, 1.0)
    Sh = np.maximum(a2*(F@C), floor)
    if zero_rule: Sh[0] = Sh[1]
    return Fi@Sh, a2, Sh
loop.psdVCE = psdVCE_wk
for phi, nArc in [(0.8, 20), (0.9, 40), (0.7, 9), (0.95, 40)]:
    acv = phi**np.arange(1024)/(1-phi**2)
    for seed in (10, 11, 12):
        arcs = loop.make_data(acv, nArc, 256, seed=seed)
        H = loop.run(arcs, 60, passes=10)
        stop = next((h for h in H if h['maxDev'] <= 1e-2), None); h = stop or H[-1]
        print(f'phi={phi} {nArc}x256 seed={seed}: {"CONVERGED" if stop else "maxPasses"} pass {h["it"]:2d} maxDev {h["maxDev"]:.1e} '
              f'level {h["level"]:.3f} (truth {1/(1-phi**2):.3f}) rho1 {h["C"][1]/h["C"][0]:.3f} sigma0 {h["s0"]:.3f}  '
              f'max level over passes {max(hh["levelNew"] for hh in H):.3g}')
