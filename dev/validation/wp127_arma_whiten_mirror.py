#!/usr/bin/env python3
"""WP-127 mirror: the ARMA whitener with an MA part (lsqvce.kernel.armaWhiten).

The whitener of an ARMA(p,q) process is A(z)/B(z): recursive, so its state
matters. The code before WP-127 ran the zero-state filter and replaced only
the leading p rows by the exact startup chol(toeplitz(r(1:p)))^-1; the state
of 1/B(z) was then built from the transient, and rows p+1.. were not
chol(C)^-1 (|F'F - inv(C)| / |inv(C)| = 1.6 % at a = [1 -0.5], b = [1 0.4]).

The fix: whiten rows 1..K exactly (K = p + N, N the decay of 1/B(z) to eps,
covModel's filterMemoryOf), then continue the recursion
    u_k = (a * x)_k - sum_j b_j u_(k-j)
from the EXACT innovations u_k = L(k,k) w_k of those rows, divided by s.
Beyond K the innovations-algorithm coefficients equal b to rounding, so the
operator IS chol(C)^-1. Batched with overlap K (later batches zero-state) it
equals the whole-arc pass.

Printed per case: |F - chol(C)^-1| / |chol(C)^-1| and batched vs whole.
Run:  python3 dev/validation/wp127_arma_whiten_mirror.py
Created 2026-10-12 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-127.
"""
import numpy as np
from scipy.signal import lfilter
from scipy.linalg import toeplitz, cholesky, solve_triangular
def acv(a,b,n,N=60000):
    h=lfilter(b,a,np.r_[1.0,np.zeros(N-1)]); return np.array([h[:N-k]@h[k:] for k in range(n)])
def fmem(p,b):
    if len(b)<2 or not np.any(b[1:]): return p
    rho=max(abs(np.roots(b))); N=max(16,int(np.ceil(np.log(np.finfo(float).eps)/np.log(rho))))
    for _ in range(12):
        h=lfilter([1.0],b,np.r_[1.0,np.zeros(N-1)])
        if np.all(abs(h[-(len(b)-1):])<=np.finfo(float).eps*abs(h).max()): break
        N*=2
    return p+N
def whiten(a,b,s,r,X,first):
    p=len(a)-1; m=X.shape[0]
    Y=lfilter(a,b,X,axis=0)/s
    if not first: return Y
    if len(b)==1 or not np.any(b[1:]):
        nS=min(p,m)
        if nS>0:
            L=cholesky(toeplitz(r[:nS]),lower=True); Y[:nS]=solve_triangular(L,X[:nS],lower=True)
        return Y
    K=min(len(r),m)
    L=cholesky(toeplitz(r[:K]),lower=True)
    W=solve_triangular(L,X[:K],lower=True)          # exact rows 1..K
    if K==m: return W
    u=np.diag(L)[:,None]*W                          # exact unnormalised innovations
    v=lfilter(a,[1.0],X,axis=0)                     # AR part, exact for k>p
    y=np.vstack([lfilter(b,[1.0],u,axis=0), v[K:]])
    z=lfilter([1.0],b,y,axis=0)
    return np.vstack([W, z[K:]/s])
def main():
    rng=np.random.default_rng(1)
    for a,b,s,m in (([1,-0.5],[1,0.4],0.7,300),([1,-0.9],[1,0.8],1.3,600),([1,-1.6,0.9,-0.25],[1,0.5,0.3],0.5,400),([1,-0.5],[1,0.95],1.0,900),([1,-0.5],[1,0.4],0.7,20)):
        a=np.array(a,float); b=np.array(b,float); p=len(a)-1
        K=fmem(p,b); r=s**2*acv(a,b,max(K,2)); C=s**2*toeplitz(acv(a,b,m))
        F=whiten(a,b,s,r,np.eye(m),True)
        Lc=cholesky(C,lower=True); Fr=np.linalg.inv(Lc)
        e1=np.linalg.norm(F-Fr)/np.linalg.norm(Fr)
        # batched: batch bs with overlap K
        bs=max(K+5, m//3); X=rng.standard_normal((m,3)); Wf=whiten(a,b,s,r,X,True); Wb=np.zeros_like(X)
        for i0 in range(0,m,bs):
            i1=min(i0+bs,m); j0=max(0,i0-K)
            Yb=whiten(a,b,s,r,X[j0:i1],i0==0); Wb[i0:i1]=Yb[i0-j0:]
        e2=np.abs(Wb-Wf).max()/np.abs(Wf).max()
        print(a,b,'m',m,'K',K,'|F-chol^-1|rel',f'{e1:.1e}','batched vs whole',f'{e2:.1e}')

        assert e1 < 1e-12 and e2 < 1e-12


if __name__ == '__main__':
    main()
