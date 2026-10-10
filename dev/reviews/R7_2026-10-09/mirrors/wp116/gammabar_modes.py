"""WP-116 (R7-09): t15b.py with the gamma-bar normalisation as a CHOICE (argv[2]):
  global  = shipped (gamma-bar over every column)
  stream  = gamma-bar per stream only
  stream+level = per-stream gamma-bar times the stream's absolute level mean(gamma2)
Original header: Mirror of adjustVCEjoint step 9 (lines 471-541) + arcVarianceFactors (robust
gammaBar over ALL columns) + streamMeans + jointModelOpt level, on three WHITE
streams sharing 4 common parameters, Pool "stream", exact traces (UiZ = U^-1):
  A: estimate "scale", arcFactors TRUE,  declared sigma 1, truth 1
  B: estimate "scale", arcFactors FALSE, declared sigma 1, truth 4   (mis-declared)
  C: estimate "law",   arcFactors FALSE, declared sigma 1, truth 1   (white)
Question: does stream A's level converge to its truth (1), or to 1/gammaBar^2?"""
import numpy as np
from scipy.linalg import toeplitz, cholesky
from kern import dctI, redundancyLags, lawSpecs, covModelToeplitzACV
from loop import psdVCE
L, m, nArc, nc = 16, 256, 9, 4
import sys; truth = {'A': float(sys.argv[3]) if len(sys.argv) > 3 else 1.0, 'B': float(sys.argv[1]), 'C': float(sys.argv[5]) if len(sys.argv) > 5 else 1.0}; MODE = sys.argv[2]
arcF = {'A': True, 'B': False, 'C': len(sys.argv) > 4 and sys.argv[4] == 'Caf'}
pol = {'A': 'scale', 'B': 'scale', 'C': 'law'}
def gammabar_robust(g):
    M = len(g); s = np.sort(g); lo = int(np.floor(M*0.25)); hi = M - int(np.floor(M*0.25))
    return s[lo:hi].mean() if hi > lo else g.mean()
def run(seed, passes=8, Tol=1e-2):
    rng = np.random.default_rng(seed)
    t = np.arange(m)*5.0; w0 = 2*np.pi/(0.5*86400)
    arcs = []
    for s in 'ABC':
        for k in range(nArc):
            tt = t + k*m*5.0
            A = np.c_[np.sin(w0*tt), np.cos(w0*tt), np.sin(2*w0*tt), np.cos(2*w0*tt)]
            arcs.append((s, A, A@np.array([.8, -.45, .3, .15]) + truth[s]*rng.normal(size=m)))
    keys = list(range(len(arcs))); s2Map = {k: 1.0 for k in keys}; Claw = np.r_[1, 0.3*np.ones(L-1)]
    first = True
    for it in range(1, passes+1):
        lvl = {s: np.mean([s2Map[k] for k in keys if arcs[k][0] == s]) for s in 'ABC'}   # streamMeans
        if first: lvl = {s: 1.0 for s in 'ABC'}
        Ws = {}
        for s in 'ABC':
            if pol[s] == 'law':
                r = lvl[s]*covModelToeplitzACV(*lawSpecs(Claw), 5.0, m)
            else:
                r = np.zeros(m); r[0] = lvl[s]                       # diag law, sigma = 1*sqrt(level)
            Ws[s] = cholesky(toeplitz(r))
        Ab = [np.linalg.solve(Ws[s].T, A) for s, A, l in arcs]; lb = [np.linalg.solve(Ws[s].T, l) for s, A, l in arcs]
        N = sum(a.T@a for a in Ab); n = sum(a.T@b for a, b in zip(Ab, lb)); x = np.linalg.solve(N, n)
        s0 = np.sqrt((sum(b@b for b in lb) - n@x)/(m*len(arcs) - nc))
        UiZ = np.linalg.solve(cholesky(N), np.eye(nc))
        omT = np.zeros((L, len(arcs))); sT = np.zeros((L, len(arcs)))
        for k, (s, A, l) in enumerate(arcs):
            sn, om = redundancyLags(Ab[k], Ws[s], UiZ, lb[k] - Ab[k]@x, L)
            s2 = 1.0 if first else s2Map[k]
            omT[:, k] = s2*Claw*om; sT[:, k] = s2*Claw*sn                 # ONE C weights every arc (shipped)
        lawc = [k for k in keys if pol[arcs[k][0]] == 'law']
        Cn, a2, _ = psdVCE(omT[:, lawc], sT[:, lawc], Claw); maxDev = np.max(np.abs(a2-1)); Claw = Cn
        s2Cols = np.array([s2Map[k] for k in keys])
        s2App = s2Cols.copy()
        for s in 'ABC':
            idx = [k for k in keys if arcs[k][0] == s]; s2App[idx] = s2Cols[idx].mean()
        alpha2m = dctI(omT).sum(0)/dctI(sT).sum(0)
        gamma2 = alpha2m*s2App
        if MODE == 'global':
            gB = gammabar_robust(np.sqrt(np.abs(gamma2))); sig2hat = gamma2/gB**2
        else:
            sig2hat = np.empty_like(gamma2); gB = np.nan
            for s in 'ABC':
                idx = [k for k in keys if arcs[k][0] == s]
                gs = gammabar_robust(np.sqrt(np.abs(gamma2[idx])))
                sig2hat[idx] = gamma2[idx]/gs**2
                if MODE == 'stream+level':
                    sig2hat[idx] *= gamma2[idx].mean()
        fS = {s: gamma2[[k for k in keys if arcs[k][0] == s and not arcF[s]]].mean() for s in 'ABC' if not arcF[s]}
        for k in keys:
            s = arcs[k][0]
            s2Map[k] = sig2hat[k] if arcF[s] else fS[s]
        first = False
        lvlNew = {s: np.mean([s2Map[k] for k in keys if arcs[k][0] == s]) for s in 'ABC'}
        print(f'  pass {it}: maxDev(C law) {maxDev:.2e}  gammaBar^2 {gB**2:.3f}  sigma0 {s0:.3f}  level A/B/C (truth 1/16/1): '
              f'{lvlNew["A"]:.3f} / {lvlNew["B"]:.3f} / {lvlNew["C"]:.3f}   reconstructed sigma A/B/C (truth 1/4/1): '
              f'{np.sqrt(lvlNew["A"])*s0:.3f} / {np.sqrt(lvlNew["B"])*s0:.3f} / {np.sqrt(lvlNew["C"])*s0:.3f}')
        if maxDev <= Tol:
            print(f'  -> CONVERGED (stopReason "converged") at pass {it}'); break
for sd in (1,):
    print('seed', sd); run(sd)
