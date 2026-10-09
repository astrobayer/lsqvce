"""Mirror: the Koch-convention robust pass inside the joint loop, on a CLEAN AR(1)
stream whose law is declared exactly. Question: what level does the loop's
Helmert step return when epochs are folded with IGG3 factors (inflateCep), and
what is level*sigma0 (the oracle's lvlB)?  Exact traces (no probes)."""
import numpy as np
from scipy.linalg import toeplitz, cholesky
from kern import dctI, redundancyLags
rng = np.random.default_rng(5)
def igg3(u, k0=1.5, k1=3.0):
    a = np.abs(u); f = np.ones_like(a)
    mid = (a > k0) & (a <= k1); f[mid] = (k0/a[mid])*((k1-a[mid])/(k1-k0))**2; f[a > k1] = 0
    return f
phi, mA, nArc, nc, nLag = 0.5, 256, 9, 4, 60
acv = phi**np.arange(mA)/(1-phi**2)          # true = declared law (process variance 1.333)
Rt = toeplitz(acv); Wt = cholesky(Rt)
t = np.arange(nArc*mA)*5.0; w0 = 2*np.pi/(0.5*86400)
Afull = np.c_[np.sin(w0*t), np.cos(w0*t), np.sin(2*w0*t), np.cos(2*w0*t)]
def run(seed, robust, passes=6):
    r = np.random.default_rng(seed)
    e = np.concatenate([Wt.T@r.normal(size=mA) for _ in range(nArc)])
    l = Afull@np.array([.8, -.45, .3, .15]) + e
    level = 1.0; w = np.ones(len(l))
    out = []
    for it in range(passes):
        # accumulate under D = diag(sqrt(level)/sqrt(w)) (R_t kron 1) D   (inflateCep + sigmaScale)
        N = np.zeros((nc, nc)); n = np.zeros(nc); ltl = 0; arcs = []
        for k in range(nArc):
            rows = slice(k*mA, (k+1)*mA)
            d = np.sqrt(w[rows])/np.sqrt(level)        # descale factors D^-1
            Ad = Afull[rows]*d[:, None]; ld = l[rows]*d
            Ab = np.linalg.solve(Wt.T, Ad); lb = np.linalg.solve(Wt.T, ld)
            N += Ab.T@Ab; n += Ab.T@lb; ltl += lb@lb; arcs.append((Ab, lb, rows))
        x = np.linalg.solve(N, n); vtv = ltl - n@x; s0sq = vtv/(len(l)-nc)
        U = cholesky(N); UiZ = np.linalg.solve(U, np.eye(nc))
        # sweep: per arc, descale by Cep (with the folded factors), temporal W from Cfcn = level*Rt
        W = cholesky(level*Rt)
        a2 = []
        for (Ab, lb, rows) in arcs:
            d = np.sqrt(w[rows])                       # descale by Cep only (level is in Cfcn)
            Ad = Afull[rows]*d[:, None]; ld = l[rows]*d
            AbarS = np.linalg.solve(W.T, Ad); ebar = np.linalg.solve(W.T, ld - Ad@x)
            s, om = redundancyLags(AbarS, W, UiZ, ebar, nLag)
            C = acv[:nLag]/acv[0]
            OmT = dctI(level*C*om); sT = dctI(level*C*s)
            a2.append(OmT.sum()/sT.sum())
        gamma2 = np.array(a2)*level
        level = gamma2.mean()                          # fStream (arcFactors = false)
        if robust:                                      # Koch: raw residual / declared Cep, MAD per stream, IGG3
            v = l - Afull@x
            s = 1.4826*np.median(np.abs(v - np.median(v)))
            w = igg3(v/s)
        out.append((level, np.sqrt(s0sq)))
    return out, w
for robust in (False, True):
    L = []; P = []
    for sd in range(6):
        o, w = run(sd, robust)
        L.append(o[-1][0]); P.append(np.sqrt(o[-1][0])*o[-1][1])
    print(f'robust={robust}: level (variance) per seed', np.round(L, 3), ' lvlB = sqrt(level)*sigma0:', np.round(P, 3),
          ' mean w %.3f, zeros %d' % (w.mean(), (w == 0).sum()))
