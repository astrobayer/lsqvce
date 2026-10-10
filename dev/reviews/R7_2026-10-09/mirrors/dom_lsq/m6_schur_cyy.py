"""Schur elimination (NormalEqSystem path, two systems with locals), back-substitution,
Cyy = s0sq*Nll^-1 + G Cxx G', dof, vtv = l'l - n'x vs direct sum, against the dense
block-arrow solve. Plus: localCovariance's s0sq_hint (verbose log) vs the true s0sq.
Plus stat.chi2Inv (the coded Wilson-Hilferty + Newton) vs scipy."""
import numpy as np
from scipy import stats, special
rs = np.random.RandomState(2)
nc, sys = 5, [(300, 4), (250, 3)]
blocks = []; Ncc = np.zeros((nc, nc)); ncv = np.zeros(nc); ltl = 0; mred = 0
dense_rows = []; dense_l = []; off = nc; xt = rs.randn(nc) * 100
for (m, nl) in sys:
    A = rs.randn(m, nc); B = rs.randn(m, nl); y = rs.randn(nl) * 50
    l = A @ xt + B @ y + rs.randn(m)
    M = np.column_stack([A, B]); N = M.T @ M; n = M.T @ l
    Nll = N[nc:, nc:]; Nlc = N[nc:, :nc]; Z = np.linalg.solve(Nll, np.column_stack([Nlc, n[nc:]]))
    Ncc += N[:nc, :nc] - Nlc.T @ Z[:, :nc]; ncv += n[:nc] - Nlc.T @ Z[:, -1]
    ltl_red = l @ l - n[nc:] @ Z[:, -1]; ltl += ltl_red; mred += m - nl
    blocks.append((Nll, Nlc, n[nc:], A, B, l, ltl_red / (m - nl)))
    dense_rows.append((A, B, off)); dense_l.append(l); off += nl
x = np.linalg.solve(Ncc, ncv); vtv = ltl - ncv @ x; dof = mred - nc; s0sq = vtv / dof
Cxx = s0sq * np.linalg.inv(Ncc)
U = off; Ad = np.zeros((sum(m for m, _ in sys), U)); r0 = 0
for (A, B, o) in dense_rows:
    Ad[r0:r0 + len(A), :nc] = A; Ad[r0:r0 + len(A), o:o + B.shape[1]] = B; r0 += len(A)
ld = np.concatenate(dense_l); xd = np.linalg.lstsq(Ad, ld, rcond=None)[0]; vd = ld - Ad @ xd
dofd = len(ld) - U; s0d = vd @ vd / dofd; Cd = s0d * np.linalg.inv(Ad.T @ Ad)
print('x  rel err %.1e ; s0sq %.10f vs dense %.10f (dof %d vs %d) ; vtv accum %.10e direct %.10e'
      % (np.linalg.norm(x - xd[:nc]) / np.linalg.norm(xd[:nc]), s0sq, s0d, dof, dofd, vtv, vd @ vd))
print('Cxx rel err %.1e' % (np.linalg.norm(Cxx - Cd[:nc, :nc]) / np.linalg.norm(Cd[:nc, :nc])))
o = nc
for (Nll, Nlc, nl_, A, B, l, hint) in blocks:
    k = Nll.shape[0]; y = np.linalg.solve(Nll, nl_ - Nlc @ x); G = np.linalg.solve(Nll, Nlc)
    Cyy = s0sq * np.linalg.inv(Nll) + G @ Cxx @ G.T
    Cxy_dense = Cd[:nc, o:o + k]
    print('  sys nl=%d: y rel err %.1e  Cyy (residuals route) rel err %.1e  Cxy = -Cxx G\' rel err %.1e ;'
          ' localCovariance s0sq_hint = %.4g vs s0sq %.4g -> Cyy(log) sd x%.1f'
          % (k, np.linalg.norm(y - xd[o:o + k]) / np.linalg.norm(xd[o:o + k]),
             np.linalg.norm(Cyy - Cd[o:o + k, o:o + k]) / np.linalg.norm(Cd[o:o + k, o:o + k]),
             np.linalg.norm(-Cxx @ G.T - Cxy_dense) / np.linalg.norm(Cxy_dense), hint, s0sq,
             np.sqrt(np.mean(np.diag(hint * np.linalg.inv(Nll) + G @ Cxx @ G.T)) / np.mean(np.diag(Cyy)))))
    o += k

def chi2Inv(p, r):
    z = -np.sqrt(2) * special.erfcinv(2 * p); x0 = r * (1 - 2 / (9 * r) + z * np.sqrt(2 / (9 * r)))**3
    if not (x0 > 0 and np.isfinite(x0)): x0 = max(r * 0.1, 1e-3)
    y = np.log(x0)
    for _ in range(100):
        xc = np.exp(y); F = special.gammainc(r / 2, xc / 2) - p
        dF = np.exp((r / 2 - 1) * np.log(xc) - xc / 2 - special.gammaln(r / 2) - (r / 2) * np.log(2) + np.log(xc))
        st = F / dF if np.isfinite(F / dF) else 0; st = max(min(st, 2), -2); yn = y - st
        if abs(yn - y) <= 1e-15 * max(abs(y), 1): y = yn; break
        y = yn
    return np.exp(y)
worst = 0
for p in (0.005, 0.025, 0.95, 0.975, 0.995):
    for r in (1, 2, 5, 50, 1e4, 1e7):
        worst = max(worst, abs(chi2Inv(p, r) / stats.chi2.ppf(p, r) - 1))
print('chi2Inv max rel err vs scipy over p in {.005..995}, r in {1..1e7}: %.1e' % worst)
