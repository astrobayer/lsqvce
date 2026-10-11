"""Candidate D: the law is a POSITIVE spectrum on L bins, piecewise-linear in f;
its covariance r = K S is the exact (alias-free) integral, and every bin's
component V_j = toeplitz(K[:, j]) is the ACV of a hat spectrum -> PSD.
Helmert per component: alpha2_j = (K' Omega)_j / (K' s)_j >= 0 automatically."""
import sys; import os; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'dom_vce'))
import numpy as np
from scipy.linalg import toeplitz, cholesky, solve_triangular
from kern import redundancyLags

def Kmat(L, m, ovs=64):
    fj = np.arange(L)/(2*(L-1)); d = fj[1]
    f = np.linspace(0, 0.5, ovs*(L-1)+1)
    H = np.maximum(0, 1 - np.abs(f[:, None] - fj[None, :])/d)       # hats
    tau = np.arange(m)
    Cs = np.cos(2*np.pi*np.outer(tau, f))
    return 2*np.trapezoid(Cs[:, :, None]*H[None, :, :], f, axis=1)    # m x L

def data(phi, nArc, m, seed):
    rng = np.random.default_rng(seed)
    acv = phi**np.arange(m)/(1-phi**2); Wt = cholesky(toeplitz(acv))
    t = np.arange(m)*5.0; w0 = 2*np.pi/(0.5*86400); xt = np.array([.8, -.45, .3, .15]); arcs = []
    for k in range(nArc):
        tt = t + k*m*5.0
        A = np.c_[np.sin(w0*tt), np.cos(w0*tt), np.sin(2*w0*tt), np.cos(2*w0*tt)]
        arcs.append((A, A@xt + Wt.T@rng.normal(size=m)))
    return arcs, xt, acv

def run(arcs, K, passes=30, tol=1e-2):
    m, L = K.shape; S = np.ones(L)
    for it in range(passes):
        r = K@S; W = cholesky(toeplitz(r))
        Ab = [solve_triangular(W, A, trans='T') for A, l in arcs]; lb = [solve_triangular(W, l, trans='T') for A, l in arcs]
        N = sum(a.T@a for a in Ab); x = np.linalg.solve(N, sum(a.T@b for a, b in zip(Ab, lb)))
        UiZ = np.linalg.inv(cholesky(N))
        Om = np.zeros(m); sS = np.zeros(m)
        for a, b in zip(Ab, lb):
            s, om = redundancyLags(a, W, UiZ, b - a@x, m)
            Om += om; sS += s
        a2 = (K.T@Om)/(K.T@sS)
        S = a2*S
        if np.max(np.abs(a2-1)) <= tol: return dict(it=it+1, conv=True, r=K@S, x=x, Q=np.linalg.inv(N), a2=a2)
    return dict(it=passes, conv=False, r=K@S, x=x, Q=np.linalg.inv(N), a2=a2)

m, L = 256, 60
K = Kmat(L, m)
for phi, nArc in [(0.5, 9), (0.8, 20), (0.9, 40)]:
    for seed in (10, 11, 12):
        arcs, xt, acv = data(phi, nArc, m, seed)
        h = run(arcs, K)
        z = (h['x'] - xt)/np.sqrt(np.diag(h['Q']))
        print(f'phi {phi} seed {seed}: {"CONV" if h["conv"] else "maxP"} @{h["it"]:2d}  r0 {h["r"][0]:.3f} (truth {acv[0]:.3f})  '
              f'rho1 {h["r"][1]/h["r"][0]:.3f}  rho5 {h["r"][5]/h["r"][0]:.3f} (truth {phi**5:.3f})  |z|max {np.abs(z).max():.2f}  min a2 {h["a2"].min():.3f}', flush=True)
