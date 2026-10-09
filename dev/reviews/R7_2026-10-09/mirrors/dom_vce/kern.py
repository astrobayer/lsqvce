"""Line-by-line mirrors of lsqvce kernels (dom_vce lane, R7)."""
import numpy as np
from scipy.linalg import toeplitz, cholesky, solve_triangular

def dctI_matrix(N, norm='ortho'):
    m = np.arange(N)[:, None]; n = np.arange(N)[None, :]
    k = np.ones(N); k[0] = 1/np.sqrt(2)
    if norm == 'ortho': k[N-1] = 1/np.sqrt(2)
    X = np.sqrt(2/(N-1)) * np.outer(k, k) * np.cos(m*n*np.pi/(N-1))
    return (X + X.T)/2

def dctI(x, norm='ortho', inverse=False):
    x = np.asarray(x, float)
    X = dctI_matrix(x.shape[0], norm)
    A = X if (not inverse or norm == 'ortho') else np.linalg.solve(X, np.eye(len(X)))
    return A @ x

def columnLagSums(X, nL):
    X = np.atleast_2d(X.T).T if X.ndim == 1 else X
    N = X.shape[0]
    nfft = int(2**np.ceil(np.log2(2*N-1)))
    F = np.fft.fft(X, nfft, axis=0)
    P = np.sum((F*np.conj(F)).real, axis=1)
    r = np.fft.ifft(P).real
    return r[:nL]

def toeplitzInvFirstCol(c):
    c = np.asarray(c, float); N = len(c)
    if N == 1: return np.array([1/c[0]])
    r = c/c[0]; a = np.array([-r[1]]); E = 1 - a[0]**2
    for n in range(2, N):              # MATLAB n = 2:N-1 -> r(n+1), r(n:-1:2)
        k = -(r[n] + a @ r[n-1:0:-1])/E
        a = np.r_[a + k*a[::-1], k]
        E = E*(1 - k**2)
    return np.r_[1, a]/(E*c[0])

def sigmaInvSubdiagSums(W, nL):
    N = W.shape[0]; dS = np.zeros(nL)
    c1 = W.T @ W[:, 0]
    x = toeplitzInvFirstCol(c1)
    z = np.r_[0, x[N-1:0:-1]]
    x1 = x[0]
    for n in range(nL):
        j = np.arange(1, N-n+1)                     # MATLAB 1-based
        wgt = (N - n - j + 1)
        dS[n] = (np.sum(x[j+n-1]*x[j-1]*wgt) - np.sum(z[j+n-1]*z[j-1]*wgt))/x1
    return dS

def redundancyLags(Abar, W, UiZ, ebar, nMax, LocalBbar=None, LocalLl=None):
    N = Abar.shape[0]; nMax = min(nMax, N)
    Winv = np.linalg.inv(W)
    Rt = Winv @ (Abar @ UiZ)
    dS = sigmaInvSubdiagSums(W, nMax)
    dR = columnLagSums(Rt, nMax)
    s = np.zeros(nMax)
    if LocalBbar is not None and LocalBbar.shape[1] > 0:
        K = Winv @ np.linalg.solve(LocalLl, LocalBbar.T).T
        dK = columnLagSums(K, nMax)
        s[0] = dS[0]-dR[0]-dK[0]; s[1:] = 2*(dS[1:]-dR[1:]-dK[1:])
    else:
        s[0] = dS[0]-dR[0]; s[1:] = 2*(dS[1:]-dR[1:])
    om = np.full(nMax, np.nan)
    if ebar is not None:
        et = Winv @ ebar
        acf = columnLagSums(et, nMax)
        om[0] = acf[0]; om[1:] = 2*acf[1:]
    return s, om

def lawSpecs(C):
    C = np.asarray(C, float); L = len(C)
    Cn = C/max(C[0], np.finfo(float).eps)
    fN = np.arange(L)/(2*(L-1))
    n = np.arange(1, L)
    S = 2*(Cn[0] + 2*(np.cos(2*np.pi*np.outer(fN, n)) @ Cn[1:]))
    S = np.maximum(S, np.finfo(float).eps*S.max())
    return fN, S

def psd2acv(f, S, dt, n):
    tau = np.arange(n)[:, None]*dt
    return np.trapezoid(S[None, :]*np.cos(2*np.pi*tau*f[None, :]), f, axis=1)

def covModelToeplitzACV(fNorm, S, dt, m, sigma=1.0, fl=1e-10, BatchSize=20000, ridge=1e-8):
    f = fNorm/dt; S = S*dt*sigma**2
    Smax = S.max(); S = np.where(S < fl*Smax, fl*Smax, S)
    nmax = min(BatchSize, m)
    nfNeed = nmax//2 + 1
    if len(f) < nfNeed:
        fF = np.linspace(f[0], f[-1], nfNeed)
        S = np.interp(fF, f, S); f = fF
        S = np.where(S < fl*Smax, fl*Smax, S)
    r = psd2acv(f, S, dt, nmax)
    r[0] *= (1+ridge)
    return r
