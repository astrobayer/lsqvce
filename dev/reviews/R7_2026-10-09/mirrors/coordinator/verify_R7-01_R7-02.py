import numpy as np
from scipy.fft import dct
# --- F1: positive Shat -> ortho inverse DCT-I -> lawSpecs cosine series
L=60
f=np.arange(L)/(2*(L-1))
for phi in (0.5,0.8,0.9):
    Shat=(1-phi**2)/np.abs(1-phi*np.exp(-2j*np.pi*f))**2*2   # AR(1) one-sided shape, positive
    C=dct(Shat,type=1,norm='ortho')                          # psdVCE.m:248 (ortho self-inverse)
    Cn=C/max(C[0],np.finfo(float).eps); n=np.arange(1,L)
    S=2*(Cn[0]+2*np.cos(2*np.pi*np.outer(f,n))@Cn[1:])         # jointModelOpt>lawSpecs
    T=np.array([[C[abs(i-j)] for j in range(L)] for i in range(L)])
    print(f"F1 phi={phi}: min Shat {Shat.min():.3g}  lawSpecs S<0 bins {np.sum(S<0)}/{L}  minS/maxS {S.min()/S.max():.3g}  minEig Toeplitz(C) {np.linalg.eigvalsh(T).min():.3g}")
# --- F2: covModel toeplitz grid nf=floor(m/2)+1, trapz psd2acv
def acv(m,phi,dt=1.0):
    nf=m//2+1; fN=1/(2*dt); fg=np.linspace(0,fN,nf)
    S=2*dt*(1-phi**2)/np.abs(1-phi*np.exp(-2j*np.pi*fg*dt))**2
    tau=np.arange(m)*dt
    return np.trapezoid(S[None,:]*np.cos(2*np.pi*np.outer(tau,fg)),fg,axis=1)
for m in (256,255,2111):
    r=acv(m,0.5)
    print(f"F2 m={m}: r1/r0={r[1]/r[0]:.4f} r[m-1]/r0={r[-1]/r[0]:.4f} (true {0.5**(m-1):.1e})")
