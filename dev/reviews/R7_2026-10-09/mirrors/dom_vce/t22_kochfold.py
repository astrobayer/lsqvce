"""Koch fold (makeLoader>inflateCep: sigma_i / sqrt(w_i), C = D R_t D) with w_i -> 1e-12
under a CORRELATED temporal law: is a zero-weight epoch REMOVED, or does it become a
phantom 'noise = 0' observation for its neighbours?  Exact linear algebra."""
import numpy as np
from scipy.linalg import toeplitz
rng = np.random.default_rng(9)
m = 256
t = np.arange(m)*5.0
A = np.c_[np.ones(m), (t-t.mean())/np.ptp(t), np.sin(2*np.pi*t/1500.)]
for phi in (0.5, 0.9):
    R = toeplitz(phi**np.arange(m))
    zero = np.sort(rng.choice(m, 13, replace=False))           # 5 % of epochs at weight 0
    keep = np.setdiff1d(np.arange(m), zero)
    w = np.ones(m); w[zero] = 1e-12
    d = np.sqrt(w)                                              # D^-1
    Ri = np.linalg.inv(R)
    P_fold = (d[:, None]*Ri)*d[None, :]                         # (D R D)^-1
    N_fold = A.T@P_fold@A
    N_del = A[keep].T@np.linalg.inv(R[np.ix_(keep, keep)])@A[keep]
    Q_fold = np.linalg.inv(N_fold); Q_del = np.linalg.inv(N_del)
    # empirical: noise ~ R (blunders irrelevant: folded rows are multiplied by 1e-6)
    X = []
    G_fold = Q_fold@A.T@P_fold; G_del = Q_del@A[keep].T@np.linalg.inv(R[np.ix_(keep, keep)])
    E_fold = G_fold@R@G_fold.T; E_del = G_del@R[np.ix_(keep, keep)]@G_del.T
    print(f'phi={phi}: 13/256 epochs at w=1e-12 via the Koch fold vs deleting them')
    print('   formal sd  fold  ', np.round(np.sqrt(np.diag(Q_fold)), 4), '   formal sd delete', np.round(np.sqrt(np.diag(Q_del)), 4))
    print('   TRUE sd of the fold estimator', np.round(np.sqrt(np.diag(E_fold)), 4), '  (formal/true =', np.round(np.sqrt(np.diag(Q_fold)/np.diag(E_fold)), 3), ')')
    # s_0 per Ellmer for the fold: trace of R_fold = P - P A N^-1 A' P in the descaled frame vs the deleted problem
    Rf = P_fold - P_fold@A@Q_fold@A.T@P_fold
    print('   E[Omega_0] (true noise R through the fold) %.2f  vs  s_0 = tr(R_fold-model) %.2f  -> level estimate %.3f'
          % (np.trace(Rf@R@Rf.T@np.linalg.pinv(P_fold) if False else (Rf@R@Rf).trace()/1), np.trace(Rf@R) if False else 0, 0) if False else '', end='')
    # Omega_0 in redundancyLags = etilde' etilde with etilde = Sigma_t^-1 (d o v), Sigma_t = R  -> e' D^-1 R^-2 D^-1 e; model expectation tr(R_ell V0) in the descaled frame
    S = Ri; Hd = (d[:, None]*A)                                   # descaled design
    Nd = Hd.T@S@Hd; Rbar = S - S@Hd@np.linalg.inv(Nd)@Hd.T@S       # the R redundancyLags traces (Sigma_t = R, descaled frame)
    Pd = np.eye(m) - Hd@np.linalg.inv(Nd)@Hd.T@S                   # residual operator in the descaled frame
    Cov_e = Pd@(d[:, None]*R*d[None, :])@Pd.T                      # covariance of descaled residual under TRUE noise R
    EOm = np.trace(S@S@Cov_e); s0 = np.trace(Rbar)
    print(f'   lag-0: E[Omega_0] under the true noise {EOm:.2f} vs s_0 {s0:.2f}  -> alpha2_0 expectation {EOm/s0:.3f} (1 = unbiased)')
