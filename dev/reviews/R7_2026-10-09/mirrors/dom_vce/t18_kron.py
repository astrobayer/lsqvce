import numpy as np
from scipy.linalg import toeplitz, cholesky
from kern import redundancyLags
rng = np.random.default_rng(4)
nE, d, P = 40, 3, 5; m = nE*d
Rt = toeplitz(0.7**np.arange(nE)); Wt = cholesky(Rt)
Cep = np.zeros((nE, d, d))
for i in range(nE):
    M = rng.normal(size=(d, d)); Cep[i] = M@M.T + d*np.eye(d)
Lk = [np.linalg.cholesky(Cep[i]) for i in range(nE)]
D = np.zeros((m, m))
for i in range(nE): D[i*d:(i+1)*d, i*d:(i+1)*d] = Lk[i]
C_epoch_major = D@np.kron(Rt, np.eye(d))@D.T     # rows (i-1)d+k  -> R_t kron I_d
A = rng.normal(size=(m, P)); l = rng.normal(size=m)
# shipped path: descale per epoch by L_i^-1, then per axis rows k::d whiten with Wt
Dinv = np.linalg.inv(D); Ad = Dinv@A; ld = Dinv@l
Ab = np.zeros_like(Ad); lb = np.zeros_like(ld)
for k in range(d):
    Ab[k::d] = np.linalg.solve(Wt.T, Ad[k::d]); lb[k::d] = np.linalg.solve(Wt.T, ld[k::d])
N1 = Ab.T@Ab; Cinv = np.linalg.inv(C_epoch_major); N2 = A.T@Cinv@A
print('N (descale + per-axis temporal) vs A^T [D(R_t kron I_d)D]^-1 A: rel err %.2e' % (np.abs(N1-N2).max()/np.abs(N2).max()))
Calt = D@np.kron(np.eye(d), Rt)@D.T
print('... vs the axis-major ordering I_d kron R_t: rel err %.2e (must be large)' % (np.abs(N1 - A.T@np.linalg.inv(Calt)@A).max()/np.abs(N2).max()))
# per-axis s_n added == tr(R (V_n kron I_d)) in the descaled frame
U = cholesky(N1); UiZ = np.linalg.solve(U, np.eye(P))
x = np.linalg.solve(N1, Ab.T@lb)
s_sum = 0; om_sum = 0
for k in range(d):
    s, om = redundancyLags(Ab[k::d], Wt, UiZ, lb[k::d]-Ab[k::d]@x, 6); s_sum = s_sum + s; om_sum = om_sum + om
Sig = np.kron(Rt, np.eye(d)); Si = np.linalg.inv(Sig); At = Ad
R = Si - Si@At@np.linalg.solve(At.T@Si@At, At.T)@Si
ed = ld - Ad@x
ref_s = []; ref_o = []
for n in range(6):
    Vn = (np.eye(nE, k=n)+np.eye(nE, k=-n)) if n else np.eye(nE)
    V = np.kron(Vn, np.eye(d)); ref_s.append(np.trace(R@V)); ref_o.append(ed@Si@V@Si@ed)
print('sum over axes s_n vs tr(R (V_n kron I_d)): rel err %.2e;  Omega_n: %.2e' %
      (np.abs(s_sum-ref_s).max()/np.abs(ref_s).max(), np.abs(om_sum-ref_o).max()/np.abs(ref_o).max()))
