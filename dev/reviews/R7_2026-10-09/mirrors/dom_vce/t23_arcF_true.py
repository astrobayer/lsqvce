"""Same closed loop with ArcFactors = TRUE (the Options default): per-arc sigma2hat =
gamma2/gammaBar^2 (robust central-half mean), stream level = mean of the arc factors."""
import numpy as np, loop
from loop import *
from kern import dctI
def gb(g):
    M = len(g); s = np.sort(g); lo = int(np.floor(M*.25)); hi = M-int(np.floor(M*.25)); return s[lo:hi].mean()
def run_af(arcs, L=60, passes=10):
    m = arcs[0][0].shape[0]; nc = 4; C = np.r_[1, 0.3*np.ones(L-1)]; s2 = np.ones(len(arcs)); first = True; H = []
    for it in range(passes):
        level = 1.0 if first else s2.mean()
        r = level*covModelToeplitzACV(*lawSpecs(C), 5.0, m); W = cholesky(toeplitz(r))
        Ab = [np.linalg.solve(W.T, A) for A, l in arcs]; lb = [np.linalg.solve(W.T, l) for A, l in arcs]
        N = sum(a.T@a for a in Ab); n = sum(a.T@b for a, b in zip(Ab, lb)); x = np.linalg.solve(N, n)
        s0 = np.sqrt((sum(b@b for b in lb)-n@x)/(m*len(arcs)-nc)); UiZ = np.linalg.solve(cholesky(N), np.eye(nc))
        omT = np.zeros((L, len(arcs))); sT = np.zeros((L, len(arcs)))
        for k, (a, b) in enumerate(zip(Ab, lb)):
            s, om = redundancyLags(a, W, UiZ, b-a@x, L); omT[:, k] = s2[k]*C*om; sT[:, k] = s2[k]*C*s
        Cn, a2, _ = psdVCE(omT, sT, C)
        a2m = dctI(omT).sum(0)/dctI(sT).sum(0); g2 = a2m*(np.full(len(arcs), s2.mean()))
        s2 = np.maximum(g2/gb(np.sqrt(np.abs(g2)))**2, 1e-8); C = Cn; first = False
        H.append((it+1, np.max(np.abs(a2-1)), level, s0, np.sqrt(level)*s0, C[1]/C[0], np.linalg.cond(toeplitz(r))))
    return H
xt = np.array([.8, -.45, .3, .15])
for phi, nArc in [(0.8, 20), (0.9, 40)]:
    acv = phi**np.arange(1024)/(1-phi**2)
    for seed in (10, 11, 12):
        H = run_af(make_data(acv, nArc, 256, seed), passes=10)
        stop = next((h for h in H if h[1] <= 1e-2), H[-1])
        print(f'ArcFactors=true phi={phi} seed={seed}: stop pass {stop[0]} maxDev {stop[1]:.1e} level {stop[2]:.3f} sigma0 {stop[3]:.3f} '
              f'sqrt(level)*sigma0 {stop[4]:.3f} (truth {np.sqrt(1/(1-phi**2)):.3f}) rho1 {stop[5]:.3f}  max cond(declared) over passes {max(h[6] for h in H):.1e}')
