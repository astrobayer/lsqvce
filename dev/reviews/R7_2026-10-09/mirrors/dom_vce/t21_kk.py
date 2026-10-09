"""Koch-Kusche prior update of fit>priorUpdate: E[Omega_k] = r_k at the true sigma,
with N = the prior-INCLUDING normal matrix (Cxx/s0sq), and the iteration's fixed point."""
import numpy as np
rng = np.random.default_rng(2)
m, u, p = 120, 30, 30
A = rng.normal(size=(m, u)); L = np.eye(u); sig_true = 0.3
def solve(l, x0, sig):
    N = A.T@A + L.T@L/sig**2; n = A.T@l + L.T@x0/sig**2
    x = np.linalg.solve(N, n); v = np.r_[l - A@x, (x0 - L@x)/sig]
    s0sq = v@v/(m + p - u); Cxx = s0sq*np.linalg.inv(N)
    return x, Cxx, s0sq
om, rr = [], []
for t in range(4000):
    xt = sig_true*rng.normal(size=u); l = A@xt + rng.normal(size=m)
    x, Cxx, s0sq = solve(l, np.zeros(u), sig_true)
    Ninv = Cxx/max(s0sq, 1e-300)
    om.append(np.sum((L@x)**2)/sig_true**2); rr.append(p - np.trace(L@Ninv@L.T)/sig_true**2)
print('at the true sigma: E[Omega_k] = %.3f, r_k = %.3f (unbiased if equal)' % (np.mean(om), np.mean(rr)))
xt = sig_true*rng.normal(size=u); l = A@xt + rng.normal(size=m); s = 0.01
for it in range(30):
    x, Cxx, s0sq = solve(l, np.zeros(u), s); Ninv = Cxx/s0sq
    o = np.sum((L@x)**2)/s**2; r = p - np.trace(L@Ninv@L.T)/s**2; s = s*np.sqrt(o/r)
print('iterated from 0.01: sigma -> %.4f (truth %.2f, one realisation)' % (s, sig_true))
