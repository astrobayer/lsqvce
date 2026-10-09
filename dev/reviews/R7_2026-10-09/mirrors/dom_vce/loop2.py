"""Closed-loop mirror of adjustVCEjoint (one stream, Pool any, arcFactors=false,
estimate "law", no locals), through the SHIPPED declaration path
lawSpecs -> covModel('toeplitz': interp + psd2acv + ridge) -> chol.
Exact traces (UiZ = U^-1) unless Z given."""
import numpy as np
from scipy.linalg import toeplitz, cholesky
from kern import dctI, redundancyLags, lawSpecs, covModelToeplitzACV

def psdVCE(omT, sT, C, zero_rule=True, floor=1e-12):
    Om = dctI(omT).sum(1); s = dctI(sT).sum(1)
    a2 = np.where(np.abs(s) > 0, Om/s, 1.0)
    Sp = dctI(C); Sh = np.maximum(a2*Sp, floor)
    if zero_rule: Sh[0] = Sh[1]
    return dctI(Sh), a2, Sh

def make_data(acv_true, nArc, m, seed, dt=5.0):
    rng = np.random.default_rng(seed)
    Wt = cholesky(toeplitz(acv_true[:m]))
    t = np.arange(m)*dt; w0 = 2*np.pi/(0.5*86400)
    arcs = []
    for k in range(nArc):
        tt = t + k*m*dt
        A = np.c_[np.sin(w0*tt), np.cos(w0*tt), np.sin(2*w0*tt), np.cos(2*w0*tt)]
        l = A@np.array([.8, -.45, .3, .15]) + Wt.T@rng.normal(size=m)
        arcs.append((A, l))
    return arcs

def run(arcs, L, passes=12, C0=None, Z=None, seed=0, zero_rule=True, dt=5.0, verbose=False, exact_whitener=False):
    m = arcs[0][0].shape[0]; nc = arcs[0][0].shape[1]
    C = np.r_[1, 0.3*np.ones(L-1)] if C0 is None else C0.copy()
    level = 1.0; hist = []
    rng = np.random.default_rng(seed)
    for it in range(passes):
        if exact_whitener:
            r = fine_acv(C, m, level)
        else:
            fN, S = lawSpecs(C)
            r = level*covModelToeplitzACV(fN, S, dt, m)
        Sig = toeplitz(r); W = cholesky(Sig)
        Ab = [np.linalg.solve(W.T, A) for A, l in arcs]; lb = [np.linalg.solve(W.T, l) for A, l in arcs]
        N = sum(a.T@a for a in Ab); n = sum(a.T@b for a, b in zip(Ab, lb))
        x = np.linalg.solve(N, n); vtv = sum(b@b for b in lb) - n@x
        s0 = np.sqrt(vtv/(m*len(arcs)-nc))
        U = cholesky(N)
        if Z is None: UiZ = np.linalg.solve(U, np.eye(nc))
        else:
            raw = np.sign(rng.normal(size=(nc, Z))); raw[raw == 0] = 1
            UiZ = np.linalg.solve(U, raw/np.sqrt(Z))
        omT = np.zeros((L, len(arcs))); sT = np.zeros((L, len(arcs)))
        for k, (a, b) in enumerate(zip(Ab, lb)):
            s, om = redundancyLags(a, W, UiZ, b - a@x, L)
            sP = np.zeros(L); oP = np.zeros(L); sP[:len(s)] = s; oP[:len(om)] = om
            s2 = level if it > 0 else 1.0
            omT[:, k] = s2*C*oP; sT[:, k] = s2*C*sP
        Cnew, a2, Sh = psdVCE(omT, sT, C, zero_rule)
        # arcVarianceFactors with s2Applied = stream level; fStream = mean(gamma2)
        Om = dctI(omT).sum(0); ss = dctI(sT).sum(0); a2m = Om/ss
        gamma2 = a2m*(level if it > 0 else 1.0)
        g2 = gamma2[np.isfinite(gamma2) & (gamma2 > 0)]
        levelNew = g2.mean()
        hist.append(dict(it=it+1, maxDev=np.max(np.abs(a2-1)), dev0=abs(a2[0]-1), dev1=np.max(np.abs(a2[1:]-1)),
                         level=level, levelNew=levelNew, s0=s0, C=C.copy(), r=r[:L].copy(), x=x))
        C = Cnew; level = levelNew
    return hist

def fine_acv(C, m, level=1.0, ovs=16):
    """Continuous inverse of lawSpecs' cosine series (fine grid, no wrap):
    r_n = int_0^0.5 S(f) cos(2 pi f n) df  -> = C_n/C_0 for n < L when the floor is inactive."""
    fN, S = lawSpecs(C)
    L = len(C); nf = max(ovs*m, 4*L) + 1
    f = np.linspace(0, 0.5, nf)
    Sf = 2*(1 + 2*(np.cos(2*np.pi*np.outer(f, np.arange(1, L))) @ (C[1:]/C[0])))
    Sf = np.maximum(Sf, np.finfo(float).eps*Sf.max())
    n = np.arange(m)
    return level*np.trapezoid(Sf[None, :]*np.cos(2*np.pi*np.outer(n, f)), f, axis=1)
