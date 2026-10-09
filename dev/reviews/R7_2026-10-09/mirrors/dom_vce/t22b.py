import numpy as np
from scipy.linalg import toeplitz
rng = np.random.default_rng(9)
m = 256; t = np.arange(m)*5.0
A = np.c_[np.ones(m), (t-t.mean())/np.ptp(t), np.sin(2*np.pi*t/1500.)]
print('Koch fold (C = D R_t D, D_ii = 1/sqrt(w_i), w = 1e-12 on k epochs) vs deleting the k epochs; formal = (A^T P A)^-1, true = G R G^T')
for phi in (0.2, 0.5, 0.8, 0.9, 0.95):
    R = toeplitz(phi**np.arange(m)); out = []
    for k in (1, 5, 13):
        z = np.sort(rng.choice(np.arange(5, m-5), k, replace=False)); keep = np.setdiff1d(np.arange(m), z)
        d = np.ones(m); d[z] = 1e-6
        P = (d[:, None]*np.linalg.inv(R))*d[None, :]
        Q = np.linalg.inv(A.T@P@A); G = Q@A.T@P; E = G@R@G.T
        Qd = np.linalg.inv(A[keep].T@np.linalg.inv(R[np.ix_(keep, keep)])@A[keep])
        out.append('k=%2d: formal/true %.2f, true/optimal %.2f' % (k, np.sqrt(Q[2, 2]/E[2, 2]), np.sqrt(E[2, 2]/Qd[2, 2])))
    print(f'phi={phi}:  ' + ' | '.join(out) + '   (sinusoid coefficient)')
