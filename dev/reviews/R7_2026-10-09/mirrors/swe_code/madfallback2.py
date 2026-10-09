# robustWeights scale chain s = 1.4826*MAD -> std(vw) -> 1 (robustWeights.m:62-67), no warning.
# Reachable when >= half the law-whitened residuals are exactly equal: rows with zero redundancy
# (r_i = 0: e.g. arcs carrying as many local parameters as rows) have v_i = 0 exactly.
import numpy as np
def scale(vw):
    med=np.median(vw); s=1.4826*np.median(np.abs(vw-med))
    if s>0 and np.isfinite(s): return s,'MAD'
    s=np.std(vw,ddof=1)
    return (s,'std') if (s>0 and np.isfinite(s)) else (1.0,'1')
rng=np.random.default_rng(3); n=10000
for frac0 in (0.49, 0.51):
    v=rng.standard_normal(n); z=rng.choice(n,int(frac0*n),replace=False); v[z]=0.0
    rest=np.setdiff1d(np.arange(n),z); bl=rng.choice(rest,int(0.03*n),replace=False); v[bl]+=20
    s,p=scale(v); fac=np.minimum(1,1.345/np.abs(v/s+1e-300))
    print('zero-redundancy rows %.0f%%: scale=%s s=%.3f  Huber factor of the 20-sigma blunders (median)=%.3f'
          % (100*frac0,p,s,np.median(fac[bl])))
# floating-point near-ties (zero-redundancy rows carry ~1e-15, not 0): the fallback is NOT taken
for frac0 in (0.49, 0.51):
    v=rng.standard_normal(n); z=rng.choice(n,int(frac0*n),replace=False); v[z]=1e-15*rng.standard_normal(z.size)
    s,p=scale(v); w=np.abs(v/s); fac=np.minimum(1,1.345/np.maximum(w,1e-300))
    clean=np.setdiff1d(np.arange(n),z)
    print('near-ties %.0f%%: scale=%s s=%.3g  median Huber factor of the ordinary (redundant) rows=%.3g'
          % (100*frac0,p,s,np.median(fac[clean])))
