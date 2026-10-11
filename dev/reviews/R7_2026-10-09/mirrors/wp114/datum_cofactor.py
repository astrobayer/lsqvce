"""WP-114 mirror (R7-05, R7-15): the datum prior's covariance and the prior in
the residual / redundancy identities.

A prior K x ~ x0 (sigma) is a DATUM when the data do not see K' at all
(N K' = 0). Then x is the minimum-constraint solution for any sigma, and its
data cofactor is N^+ -- independent of sigma. The solve's Q = (N + K'K/s^2)^-1
equals N^+ + s^2 K'(KK')^-2 K; subtracting that term from Q cancels catastrophically
when s is large. The stable form re-weights the datum rows to the data's scale:
  Qmc = (N + c^2 K'K)^-1 - K'(KK')^-2 K / c^2,   c = sqrt(|N|_F) / |K|_F
both terms O(1/|N|), from the factor of [R; cK] (square root, no N formed).
"""
import numpy as np

def network(rs, u):
    pairs = [(i, (i + 1) % u) for i in range(u)] + [tuple(rs.choice(u, 2, replace=False)) for _ in range(u // 2)]
    A = np.zeros((len(pairs), u))
    for k, (i, j) in enumerate(pairs): A[k, i] = -1; A[k, j] = 1
    l = A @ (rs.randn(u) * 10) + 0.01 * rs.randn(len(pairs))
    return A, l

def solve_with_prior(A, l, K, x0, s):
    M = np.vstack([np.column_stack([A, l]), np.column_stack([K / s, x0 / s])])
    G = np.linalg.qr(M, mode='r'); n = A.shape[1]
    R1 = G[:n, :n]; x = np.linalg.solve(R1, G[:n, -1])
    vtv = G[n, -1] ** 2 if G.shape[0] > n else 0.0
    dof = A.shape[0] + K.shape[0] - n
    Ri = np.linalg.inv(R1)
    return x, vtv / dof, Ri @ Ri.T, dof

def datum_cofactor(A, K):
    R = np.linalg.qr(A, mode='r')
    N = R.T @ R
    c = np.sqrt(np.linalg.norm(N)) / np.linalg.norm(K)
    Rm = np.linalg.qr(np.vstack([R, c * K]), mode='r')
    Rmi = np.linalg.inv(Rm[:A.shape[1], :A.shape[1]])
    Dm = K.T @ np.linalg.inv(K @ K.T) / c
    is_datum = np.linalg.norm(N @ K.T) <= 1e-10 * np.linalg.norm(N) * np.linalg.norm(K)
    return Rmi @ Rmi.T - Dm @ Dm.T, is_datum

rs = np.random.RandomState(0)
A, l = network(rs, 6)
K = np.ones((1, 6))
Np = np.linalg.pinv(A.T @ A)
print('levelling network u = 6, datum K = ones (sum of heights ~ 0)')
for s in (1.0, 1e3, 1e6):
    x, s0sq, Q, dof = solve_with_prior(A, l, K, np.zeros(1), s)
    Qmc, dat = datum_cofactor(A, K)
    print('sigma %-6g: datum %s  se(solve Q) %s  se(min-constraint) %s  |Qmc - pinv(N)| %.1e  dof %d'
          % (s, dat, np.round(np.sqrt(s0sq * np.diag(Q)), 5), np.round(np.sqrt(s0sq * np.diag(Qmc)), 5),
             np.abs(Qmc - Np).max(), dof))
# naive subtraction at sigma 1e6
x, s0sq, Q, dof = solve_with_prior(A, l, K, np.zeros(1), 1e6)
Qn = Q - 1e12 * K.T @ np.linalg.inv(K @ K.T) @ np.linalg.inv(K @ K.T) @ K
print('naive Q - s^2 K\'(KK\')^-2 K at sigma 1e6: |err| %.1e (vs stable %.1e)' % (np.abs(Qn - Np).max(), np.abs(datum_cofactor(A, K)[0] - Np).max()))
# a regularising prior is NOT a datum
Kr = np.eye(6)[:2]
print('ridge rows on heights 1, 2: datum detected?', datum_cofactor(A, Kr)[1])

# R7-15: the prior in the residual / redundancy identities (ridge prior, full-rank data)
rs = np.random.RandomState(1)
m, u = 40, 6
A2 = rs.randn(m, u); l2 = A2 @ rs.randn(u) + rs.randn(m)
for s in (10.0, 1.0, 0.1):
    K2 = np.eye(u); t2 = np.zeros(u)
    x, s0sq, Q, dof = solve_with_prior(A2, l2, K2, t2, s)
    v = l2 - A2 @ x; vp = (K2 @ x - t2) / s
    s0_res_old = v @ v / (m - u)
    s0_res_new = (v @ v + vp @ vp) / (m + u - u)
    h = np.einsum('ij,jk,ik->i', A2, Q, A2)
    rK = u - np.trace(K2 @ Q @ K2.T) / s ** 2
    print('ridge sigma %-4g: sol s0^2 %.6f | residuals old %.6f (rel %.1e) new %.6f (rel %.1e) | sum r_data %.4f + r_K %.4f = %.4f vs m+nk-u = %d'
          % (s, s0sq, s0_res_old, abs(s0_res_old / s0sq - 1), s0_res_new, abs(s0_res_new / s0sq - 1),
             np.sum(1 - h), rK, np.sum(1 - h) + rK, m))
