import numpy as np
from scipy.linalg import toeplitz, cholesky
from kern import lawSpecs, covModelToeplitzACV
rng = np.random.default_rng(7)
for phi in (0.5, 0.95):
  for m in (256, 2160, 2111):
    C = phi**np.arange(60)
    r = covModelToeplitzACV(*lawSpecs(C), 5.0, m); acvT = phi**np.arange(m)
    Ws, Wt = cholesky(toeplitz(r)), cholesky(toeplitz(acvT))
    t = np.arange(m)*5.0
    A = np.c_[np.ones(m), (t-t.mean())/np.ptp(t), np.sin(2*np.pi*t/5400.)]
    def gls(W, l):
        Ab = np.linalg.solve(W.T, A); lb = np.linalg.solve(W.T, l); N = Ab.T@Ab
        return np.linalg.solve(N, Ab.T@lb), np.sqrt(np.diag(np.linalg.inv(N)))
    Es, Ei = [], []
    for k in range(40):
        l = Wt.T@rng.normal(size=m); xs, fs = gls(Ws, l); xi, fi = gls(Wt, l); Es.append(xs); Ei.append(xi)
    Es, Ei = np.array(Es), np.array(Ei)
    print(f'phi={phi} m={m:5d}: declared corr(1st,last)={r[m-1]/r[0]:.3f}; empirical sd [bias,slope,sin] shipped {np.round(Es.std(0),3)} '
          f'vs true-law GLS {np.round(Ei.std(0),3)}; shipped FORMAL {np.round(fs,3)}')
