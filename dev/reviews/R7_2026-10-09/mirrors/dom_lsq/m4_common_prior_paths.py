"""Common prior (CommonPrior / data.priors) through the routes:
 (1) BlockAccumulator.solve prior branch: vtv = ltl + x0'x0/s^2 - n'x ; dof = m + nk - nc ; Cxx = s0sq*N_tot^-1
 (2) engine.residuals:   vtv = sum over DATA arcs of w'w ; dof = mTot - (nc + sum nl)   (no prior term)
 (3) engine.redundancy:  r_i = 1 - h_i with h from N_tot (prior included); check vs mTot - nParTot
 (4) Cxx for a datum prior vs the minimum-constraint covariance s0sq * pinv(N)."""
import numpy as np
pairs = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)]; A = np.zeros((5, 4))
for k, (i, j) in enumerate(pairs): A[k, i] = -1; A[k, j] = 1
rs = np.random.RandomState(5); l = A @ np.array([0, 1.5, -0.7, 2.2]) + 0.01 * rs.randn(5)
B = np.ones((4, 1)); u = 4; m = 5
print('--- datum prior L = B\' on a free levelling network (test_fit fixture) ---')
for sig in (1.0, 1e3):
    K = B.T; x0 = np.zeros(1); nk = 1
    N = A.T @ A + K.T @ K / sig**2; n = A.T @ l + K.T @ x0 / sig**2
    x = np.linalg.solve(N, n); ltl = l @ l
    vtv_sol = ltl + x0 @ x0 / sig**2 - n @ x; dof_sol = m + nk - u; s0_sol = np.sqrt(vtv_sol / dof_sol)
    v = l - A @ x; vtv_res = v @ v; dof_res = m - u; s0_res = np.sqrt(vtv_res / max(dof_res, 1))
    h = np.einsum('ij,jk,ik->i', A, np.linalg.inv(N), A); rsum = np.sum(1 - h)
    Cxx = vtv_sol / dof_sol * np.linalg.inv(N); Cmc = vtv_sol / dof_sol * np.linalg.pinv(A.T @ A)
    print('sigma_prior=%g: B\'x=%.1e  sol: vtv=%.6e dof=%d s0=%.6f | residuals(): vtv=%.6e dof=%d s0=%.6f  rel %.2e (guard 1e-6)'
          % (sig, (B.T @ x)[0], vtv_sol, dof_sol, s0_sol, vtv_res, dof_res, s0_res, abs(s0_res - s0_sol) / s0_sol))
    print('   sum r_i(data) = %.6f vs mTot-nParTot = %d' % (rsum, m - u))
    print('   se(x) reported = %s   minimum-constraint se = %s' % (np.round(np.sqrt(np.diag(Cxx)), 5), np.round(np.sqrt(np.diag(Cmc)), 5)))

print('\n--- an informative prior (regularisation) on a full-rank problem ---')
rs = np.random.RandomState(1); m, u = 40, 6
A = rs.randn(m, u); xt = rs.randn(u); l = A @ xt + 0.1 * rs.randn(m)
K = np.eye(u); x0 = np.zeros(u)
for sig in (10.0, 1.0, 0.1):
    N = A.T @ A + K.T @ K / sig**2; n = A.T @ l + K.T @ x0 / sig**2
    x = np.linalg.solve(N, n)
    vtv_sol = l @ l + x0 @ x0 / sig**2 - n @ x; dof_sol = m + u - u; s0_sol = np.sqrt(vtv_sol / dof_sol)
    v = l - A @ x; s0_res = np.sqrt(v @ v / (m - u))
    # diagnoseAccum route: batches = data + prior rows; r = n - uFull with n incl prior rows
    nD = m + u; rD = nD - u
    h = np.einsum('ij,jk,ik->i', A, np.linalg.inv(N), A); rdata = np.sum(1 - h)
    print('sigma=%5g  sol.s0=%.5f (dof %d)  residuals().s0=%.5f (dof %d)  rel diff %.2e  -> sigma0Mismatch %s;'
          '  redundancy(): sum r(data)=%.3f vs m-u=%d relChk %.1e -> identityViolated %s'
          % (sig, s0_sol, dof_sol, s0_res, m - u, abs(s0_res - s0_sol) / s0_sol, abs(s0_res - s0_sol) / s0_sol > 1e-6,
             rdata, m - u, abs(rdata - (m - u)) / (m - u), abs(rdata - (m - u)) / (m - u) > 1e-6))
