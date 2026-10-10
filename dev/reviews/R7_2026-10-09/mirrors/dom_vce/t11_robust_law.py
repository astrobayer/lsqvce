"""Koch-convention robust pass with the LAW re-estimated (Pool stream, arcFactors
false), shipped declaration path, exact traces: the clean-stream oracle of
test_vcePolicy (AR(1) 0.5, 9 arcs x 256, level 1, IGG3) -- what is
lvlB = sqrt(level) * sigma0 ?  Also: the per-row contribution to s_0 vs Omega_0
of a zero-weight row."""
import numpy as np
from scipy.linalg import toeplitz, cholesky
from kern import dctI, redundancyLags, lawSpecs, covModelToeplitzACV
from t3_robust_level import igg3
from loop import psdVCE
def run(seed, robust, passes=16, phi=0.5, nArc=9, m=256, L=60, dt=5.0):
    rng = np.random.default_rng(seed)
    Wt = cholesky(toeplitz(phi**np.arange(m)/(1-phi**2)*(1-phi**2)))  # innovation 1 -> process var 1/(1-phi^2)... use unit-innovation:
    Wt = cholesky(toeplitz(phi**np.arange(m)/(1-phi**2)))
    t = np.arange(nArc*m)*dt; w0 = 2*np.pi/(0.5*86400)
    A = np.c_[np.sin(w0*t), np.cos(w0*t), np.sin(2*w0*t), np.cos(2*w0*t)]
    l = A@np.array([.8, -.45, .3, .15]) + np.concatenate([Wt.T@rng.normal(size=m) for _ in range(nArc)])
    C = np.r_[1, 0.3*np.ones(L-1)]; level = 1.0; w = np.ones(len(l))
    if robust:                                         # PASS 0
        r = covModelToeplitzACV(*lawSpecs(C), dt, m); W = cholesky(toeplitz(r))
        N = np.zeros((4, 4)); n = np.zeros(4)
        for k in range(nArc):
            rows = slice(k*m, (k+1)*m); Ab = np.linalg.solve(W.T, A[rows]); lb = np.linalg.solve(W.T, l[rows]); N += Ab.T@Ab; n += Ab.T@lb
        v = l - A@np.linalg.solve(N, n); s = 1.4826*np.median(np.abs(v-np.median(v))); w = igg3(v/s)
    for it in range(passes):
        r = level*covModelToeplitzACV(*lawSpecs(C), dt, m); W = cholesky(toeplitz(r))
        N = np.zeros((4, 4)); n = np.zeros(4); ltl = 0; parts = []
        for k in range(nArc):
            rows = slice(k*m, (k+1)*m); d = np.sqrt(w[rows])
            Ab = np.linalg.solve(W.T, A[rows]*d[:, None]); lb = np.linalg.solve(W.T, l[rows]*d)
            N += Ab.T@Ab; n += Ab.T@lb; ltl += lb@lb; parts.append((Ab, lb))
        x = np.linalg.solve(N, n); s0 = np.sqrt((ltl - n@x)/(len(l)-4))
        UiZ = np.linalg.solve(cholesky(N), np.eye(4))
        omT = np.zeros((L, nArc)); sT = np.zeros((L, nArc))
        for k, (Ab, lb) in enumerate(parts):
            s, om = redundancyLags(Ab, W, UiZ, lb - Ab@x, L)
            sc = level if it > 0 or robust else 1.0
            omT[:, k] = level*C*om; sT[:, k] = level*C*s
        Cn, a2, _ = psdVCE(omT, sT, C)
        g2 = (dctI(omT).sum(0)/dctI(sT).sum(0))*level
        C = Cn; level = g2.mean()
        if robust:
            v = l - A@x; sc = 1.4826*np.median(np.abs(v-np.median(v))); w = igg3(v/sc)
    return np.sqrt(level)*s0, level, s0, C[1]/C[0], w
for robust in (False, True):
    res = [run(sd, robust) for sd in range(5)]
    print(f'robust={robust}: lvlB = sqrt(level)*sigma0 per seed', np.round([r[0] for r in res], 3),
          ' level', np.round([r[1] for r in res], 3), ' sigma0', np.round([r[2] for r in res], 3), ' rho1', np.round([r[3] for r in res], 3))
