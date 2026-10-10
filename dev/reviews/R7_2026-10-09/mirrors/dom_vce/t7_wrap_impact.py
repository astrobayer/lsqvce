"""Impact of the circulant declaration (covModel 'toeplitz' at nf = floor(m/2)+1)
on a GLS solve, with the law fed in EXACTLY (C = true AR(1) ACV on L = 60 lags,
so the only error is the declaration path).  Even vs odd arc length."""
import numpy as np
from scipy.linalg import toeplitz, cholesky
from kern import lawSpecs, covModelToeplitzACV
rng = np.random.default_rng(3)
phi, L, dt, nArc = 0.5, 60, 5.0, 30
C = phi**np.arange(L)                      # unit-variance AR(1) shape (lawSpecs normalises)
for m in (256, 255, 2160, 2111):
    acvT = phi**np.arange(m)
    r = covModelToeplitzACV(*lawSpecs(C), dt, m)
    Ssh, Sid = toeplitz(r), toeplitz(acvT)
    print(f'm={m:5d}: corr(epoch 1, epoch m) declared {r[m-1]/r[0]:.4f} (truth {acvT[m-1]:.1e}); cond declared {np.linalg.cond(Ssh):.2e} vs truth {np.linalg.cond(Sid):.2e}')
    Wt = cholesky(Sid); Ws = cholesky(Ssh)
    t = np.arange(m)*dt
    # design: per-arc bias + slope (what an arc-wise local model looks like) + one common sinusoid
    A = np.c_[np.ones(m), (t - t.mean())/np.ptp(t), np.sin(2*np.pi*t/5400.)]
    def gls(Sw, l):
        Ab = np.linalg.solve(Sw.T, A); lb = np.linalg.solve(Sw.T, l)
        N = Ab.T@Ab; x = np.linalg.solve(N, Ab.T@lb); v = lb - Ab@x
        return x, np.linalg.inv(N), v@v/(m-3)
    E_sh, E_id, s0s = [], [], []
    for k in range(nArc):
        l = Wt.T@rng.normal(size=m)             # truth x = 0
        x1, Q1, s1 = gls(Ws, l); x2, Q2, s2 = gls(Wt, l)
        E_sh.append(x1); E_id.append(x2); s0s.append(s1)
        if k == 0: Qsh, Qid = Q1, Q2
    E_sh, E_id = np.array(E_sh), np.array(E_id)
    print('    empirical sd of x (slope, sinusoid): shipped %s  ideal %s ; formal shipped %s ; mean sigma0^2 shipped %.3f'
          % (np.round(E_sh.std(0)[1:], 4), np.round(E_id.std(0)[1:], 4), np.round(np.sqrt(np.diag(Qsh))[1:], 4), np.mean(s0s)))
