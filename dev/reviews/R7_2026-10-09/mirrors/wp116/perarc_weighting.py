"""WP-116 part 2 (R7-13): do the per-arc factors WEIGHT the solve?  Mirror of
adjustVCEjoint step 9 + jointModelOpt on three WHITE streams sharing 4 parameters
(the gammabar_modes.py geometry), Monte Carlo over seeds.
  A: "scale", arcFactors TRUE, 9 arcs, arc 0 has sigma OUT (default 10), others 1
  B: "scale", arcFactors FALSE, truth sigma 2
  C: "law",   arcFactors FALSE, truth sigma 1 (white)
MODE (argv[1]):
  shipped615  weighting: stream MEAN level; eq. 6.5.28 x stream mean;
              sigma2_m = gamma2_m / gammaBar_s^2 * mean(gamma2_s)   (v6.1.5)
  meanabs     weighting: stream MEAN; sigma2_m = gamma2_m (absolute Helmert)
  perarc      weighting: PER-ARC factor; eq. 6.5.28 x the per-arc factor;
              sigma2_m = gamma2_m
Reports: A's factors (truth 100 / 1), the stream levels, and the honesty of Cxx:
true sd of x over seeds vs the mean formal sd (sqrt diag s0^2 N^-1)."""
import sys, numpy as np
from scipy.linalg import toeplitz, cholesky
from kern import dctI, redundancyLags, lawSpecs, covModelToeplitzACV
from loop import psdVCE
MODE = sys.argv[1]; OUT = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
NSEED = int(sys.argv[3]) if len(sys.argv) > 3 else 40
TRACE = NSEED == 1   # one seed: print the pass trajectory
L, m, nArc, nc = 16, 256, 9, 4
truth = {'B': 2.0, 'C': 1.0}
arcF = {'A': True, 'B': False, 'C': False}; pol = {'A': 'scale', 'B': 'scale', 'C': 'law'}
xt = np.array([.8, -.45, .3, .15])
def gbr(g):
    M = len(g); s = np.sort(g); lo = int(np.floor(M*.25)); hi = M-int(np.floor(M*.25))
    return s[lo:hi].mean() if hi > lo else g.mean()
def run(seed, passes=15, Tol=1e-3):
    rng = np.random.default_rng(seed); t = np.arange(m)*5.0; w0 = 2*np.pi/(0.5*86400)
    arcs = []
    for s in 'ABC':
        for k in range(nArc):
            tt = t + k*m*5.0
            A = np.c_[np.sin(w0*tt), np.cos(w0*tt), np.sin(2*w0*tt), np.cos(2*w0*tt)]
            sd = (OUT if k == 0 else 1.0) if s == 'A' else truth[s]
            arcs.append((s, A, A@xt + sd*rng.normal(size=m)))
    keys = range(len(arcs)); s2 = np.ones(len(arcs)); Claw = np.r_[1, 0.3*np.ones(L-1)]
    for it in range(1, passes+1):
        lvl = {s: s2[[k for k in keys if arcs[k][0] == s]].mean() for s in 'ABC'}
        Ws = []
        for k, (s, A, l) in enumerate(arcs):
            w = s2[k] if (MODE == 'perarc' and arcF[s]) else lvl[s]
            r = w*covModelToeplitzACV(*lawSpecs(Claw), 5.0, m) if pol[s] == 'law' else np.r_[w, np.zeros(m-1)]
            Ws.append(cholesky(toeplitz(r)))
        Ab = [np.linalg.solve(W.T, a[1]) for W, a in zip(Ws, arcs)]; lb = [np.linalg.solve(W.T, a[2]) for W, a in zip(Ws, arcs)]
        N = sum(a.T@a for a in Ab); n = sum(a.T@b for a, b in zip(Ab, lb)); x = np.linalg.solve(N, n)
        s0sq = (sum(b@b for b in lb) - n@x)/(m*len(arcs) - nc)
        UiZ = np.linalg.solve(cholesky(N), np.eye(nc))
        omT = np.zeros((L, len(arcs))); sT = np.zeros((L, len(arcs)))
        for k in keys:
            sn, om = redundancyLags(Ab[k], Ws[k], UiZ, lb[k]-Ab[k]@x, L)
            omT[:, k] = s2[k]*Claw*om; sT[:, k] = s2[k]*Claw*sn
        lawc = [k for k in keys if pol[arcs[k][0]] == 'law']
        Cn, a2, _ = psdVCE(omT[:, lawc], sT[:, lawc], Claw); maxDev = np.max(np.abs(a2-1)); Claw = Cn
        s2App = s2.copy()
        if MODE != 'perarc':
            for s in 'ABC':
                idx = [k for k in keys if arcs[k][0] == s]; s2App[idx] = s2[idx].mean()
        gamma2 = dctI(omT).sum(0)/dctI(sT).sum(0)*s2App
        new = s2.copy()
        for s in 'ABC':
            idx = np.array([k for k in keys if arcs[k][0] == s]); g = gamma2[idx]
            if not arcF[s]:
                new[idx] = g.mean()
            elif MODE == 'shipped615':
                new[idx] = g/gbr(np.sqrt(np.abs(g)))**2*g.mean()
            else:
                new[idx] = g
        dev = np.max(np.abs(new/s2-1)); s2 = new
        if TRACE: print(f'    pass {it}: maxDev {maxDev:.1e} levelDev {dev:.1e}  bad/good {s2[0]:.2f}/{np.median(s2[1:nArc]):.3f}')
        if maxDev <= Tol and dev <= Tol: break
    # the returned solution under the returned model (R7-11)
    lvl = {s: s2[[k for k in keys if arcs[k][0] == s]].mean() for s in 'ABC'}
    Ab = []; lb = []
    for k, (s, A, l) in enumerate(arcs):
        w = s2[k] if (MODE == 'perarc' and arcF[s]) else lvl[s]
        r = w*covModelToeplitzACV(*lawSpecs(Claw), 5.0, m) if pol[s] == 'law' else np.r_[w, np.zeros(m-1)]
        W = cholesky(toeplitz(r)); Ab.append(np.linalg.solve(W.T, A)); lb.append(np.linalg.solve(W.T, l))
    N = sum(a.T@a for a in Ab); n = sum(a.T@b for a, b in zip(Ab, lb)); x = np.linalg.solve(N, n)
    s0sq = (sum(b@b for b in lb) - n@x)/(m*len(arcs) - nc)
    return it, x, np.sqrt(s0sq*np.diag(np.linalg.inv(N))), s2, lvl, s0sq
R = [run(sd) for sd in range(NSEED)]
X = np.array([r[1] for r in R]); F = np.array([r[2] for r in R])
s2A = np.array([r[3][:nArc] for r in R])
print(f'MODE {MODE}  OUT {OUT}  seeds {NSEED}  passes median {int(np.median([r[0] for r in R]))} max {max(r[0] for r in R)}')
print(f'  A factors: bad arc {np.median(s2A[:,0]):.2f} (truth {OUT**2:.0f})  good arcs {np.median(s2A[:,1:]):.3f} (truth 1)')
print(f'  levels A/B/C x s0^2 (truth {(OUT**2+8)/9:.2f} / 4 / 1): ' + ' / '.join(f'{np.median([r[4][s]*r[5] for r in R]):.3f}' for s in 'ABC'))
print(f'  true sd x  {np.std(X-xt, axis=0, ddof=1).round(5)}')
print(f'  formal sd  {F.mean(0).round(5)}   ratio formal/true {(F.mean(0)/np.std(X-xt, axis=0, ddof=1)).round(2)}')
