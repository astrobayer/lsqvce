# Mirror of test_streamUtils.m:105-127: can var(ew)=1 at RelTol 0.25 tell a correct
# separable whitener from defective ones? e = sig * L * AR(ar) per axis, epoch-major.
import numpy as np
from scipy.signal import lfilter
d,nE=3,400; ar=[1,-0.5,0.1]; sig=0.9
Ld=np.array([[1,0,0],[0.3,1,0],[0.2,0.4,2.5]]); Dep=Ld@Ld.T; L=np.linalg.cholesky(Dep)
res={k:[] for k in ['correct','diagD_only','D_ignored(AR only)','AR_ignored(D only)','L_transposed','sigma_ignored']}
for seed in range(200):
    rs=np.random.default_rng(seed)
    W=np.vstack([lfilter([1],ar,rs.standard_normal(nE)) for _ in range(d)])
    W=sig*(L@W)            # d x nE
    def arw(X): return np.vstack([lfilter(ar,[1],x) for x in X])   # exact AR inverse filter (ignoring startup)
    Li=np.linalg.inv(L)
    res['correct'].append(np.var(arw(Li@W/sig)[:,5:],ddof=1))
    Dd=np.diag(1/np.sqrt(np.diag(Dep)))
    res['diagD_only'].append(np.var(arw(Dd@W/sig)[:,5:],ddof=1))
    res['D_ignored(AR only)'].append(np.var(arw(W/sig)[:,5:],ddof=1))
    res['AR_ignored(D only)'].append(np.var((Li@W/sig)[:,5:],ddof=1))
    res['L_transposed'].append(np.var(arw(np.linalg.inv(L.T)@W/sig)[:,5:],ddof=1))
    res['sigma_ignored'].append(np.var(arw(Li@W)[:,5:],ddof=1))
for k,v in res.items():
    v=np.array(v); print('%-22s median var(ew)=%.3f  min %.3f max %.3f  passes RelTol0.25: %3d%%'%(k,np.median(v),v.min(),v.max(),100*np.mean(np.abs(v-1)<=0.25)))
