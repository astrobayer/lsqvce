"""Attribute the GLS damage to the WRAP alone: same 60-lag law C, declared (i) by the
shipped covModel 'toeplitz' path, (ii) by the same S integrated on a 16x finer grid
(no wrap within the arc).  Truth = AR(1).  Sinusoid + slope + bias per arc."""
import numpy as np
from scipy.linalg import toeplitz, cholesky
from kern import lawSpecs, covModelToeplitzACV
from loop import fine_acv
rng = np.random.default_rng(7)
for phi in (0.5, 0.9):
  for m in (256, 255, 2160, 2111):
    C = phi**np.arange(60)
    rs = covModelToeplitzACV(*lawSpecs(C), 5.0, m); rf = fine_acv(C, m); acvT = phi**np.arange(m)
    t = np.arange(m)*5.0; A = np.c_[np.ones(m), (t-t.mean())/np.ptp(t), np.sin(2*np.pi*t/5400.)]
    Rt = toeplitz(acvT)
    def ev(r):
        Pi = np.linalg.inv(toeplitz(r)); Q = np.linalg.inv(A.T@Pi@A); G = Q@A.T@Pi; E = G@Rt@G.T
        return np.sqrt(np.diag(Q)), np.sqrt(np.diag(E))
    fs, ts = ev(rs); ff, tf = ev(rf)
    print(f'phi={phi} m={m:5d}: TRUE sd [bias slope sin]  shipped {np.round(ts,3)}  no-wrap {np.round(tf,3)} | FORMAL shipped {np.round(fs,3)} no-wrap {np.round(ff,3)}')
