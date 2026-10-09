"""Mirror of lsqvce BlockAccumulator with Backend 'qr' vs QRSystem in isolation.

Line-by-line mirror of:
  QRSystem.add            R = triu(qr([R; [A l]],0)), truncated to dim+1 rows
  QRSystem.blocks         N = R1'R1, nn = R1'r, ltl = r'r + R(d+1,end)^2
  BlockAccumulator.preeliminate  nl==0: red = (Ncc, nc_rhs, ltl, mObs)
                                 nl>0 : Schur in NORMAL form (Nll \\ [Nlc nl])
  QRSystem.merge          chol(red.N) -> Rc, rhs = Rc'\\n  (else eig range-space)
  QRSystem.solve          rank test min|diag R1| <= max(size)*eps(max), x = R1\\r
"""
import numpy as np
import scipy.linalg as sla

eps = np.finfo(float).eps


class QRSystem:
    def __init__(s, dim):
        s.dim = dim; s.R = np.zeros((0, dim + 1)); s.m = 0
    def add(s, A, l):
        M = np.vstack([s.R, np.column_stack([A, l])])
        R = np.triu(np.linalg.qr(M, mode='r'))
        s.R = R[:min(R.shape[0], s.dim + 1), :]; s.m += len(l)
    def blocks(s, ic, il):
        d = s.dim; R1 = np.zeros((d, d)); r = np.zeros(d); k = min(s.R.shape[0], d)
        R1[:k] = s.R[:k, :d]; r[:k] = s.R[:k, -1]
        N = R1.T @ R1; nn = R1.T @ r
        ltl = r @ r + (s.R[d, -1] ** 2 if s.R.shape[0] > d else 0.0)
        return N[np.ix_(il, il)], N[np.ix_(il, ic)], N[np.ix_(ic, ic)], nn[il], nn[ic], ltl, s.m
    def merge(s, N, n, ltl, dof):
        try:
            Rc = np.linalg.cholesky(N).T            # upper, MATLAB chol
            rhs = sla.solve_triangular(Rc.T, n, lower=True)
            path = 'chol'
        except np.linalg.LinAlgError:
            dv, V = np.linalg.eigh((N + N.T) / 2)
            keep = (dv > max(N.shape) * np.spacing(dv.max())) & (dv > 0)
            Rc = np.diag(np.sqrt(dv[keep])) @ V[:, keep].T
            rhs = (V[:, keep].T @ n) / np.sqrt(dv[keep])
            path = 'eig(drop %d)' % (~keep).sum()
        rho2 = ltl - rhs @ rhs
        blk = np.column_stack([Rc, rhs])
        if rho2 > 0:
            row = np.zeros(s.dim + 1); row[-1] = np.sqrt(rho2); blk = np.vstack([blk, row])
        M = np.vstack([s.R, blk]); R = np.triu(np.linalg.qr(M, mode='r'))
        s.R = R[:min(R.shape[0], s.dim + 1), :]; s.m += dof
        return path
    def solve(s):
        d = s.dim
        if s.R.shape[0] < d: return None, 'rankDeficient'
        R1 = s.R[:d, :d]; r = s.R[:d, -1]
        dR = np.abs(np.diag(R1))
        if dR.min() <= max(R1.shape) * np.spacing(dR.max()): return None, 'singularR'
        return sla.solve_triangular(R1, r), 'ok'


def accumulator_qr(A, B, l):
    """BlockAccumulator(nc,'Backend','qr') + one system with nl locals."""
    nc, nl = A.shape[1], B.shape[1]
    loc = QRSystem(nc + nl); loc.add(np.column_stack([A, B]), l)
    ic = np.arange(nc); il = nc + np.arange(nl)
    Nll, Nlc, Ncc, nlr, ncr, ltl, m = loc.blocks(ic, il)
    if nl == 0:
        red = (Ncc, ncr, ltl, m)
    else:
        Z = np.linalg.solve(Nll, np.column_stack([Nlc, nlr]))
        red = (Ncc - Nlc.T @ Z[:, :nc], ncr - Nlc.T @ Z[:, -1], ltl - nlr @ Z[:, -1], m - nl)
    g = QRSystem(nc); path = g.merge(*red)
    x, st = g.solve()
    return x, st, path


def normal_eq(A, B, l):
    """NormalEqSystem path (Schur in normal form, chol solve)."""
    nc, nl = A.shape[1], B.shape[1]
    M = np.column_stack([A, B]); N = M.T @ M; n = M.T @ l
    Ncc = N[:nc, :nc]; nn = n[:nc]
    if nl:
        Nll = N[nc:, nc:]; Nlc = N[nc:, :nc]; Z = np.linalg.solve(Nll, np.column_stack([Nlc, n[nc:]]))
        Ncc = Ncc - Nlc.T @ Z[:, :nc]; nn = nn - Nlc.T @ Z[:, -1]
    try:
        L = np.linalg.cholesky(Ncc); return sla.cho_solve((L, True), nn), 'ok'
    except np.linalg.LinAlgError:
        return None, 'notSPD'


def qr_isolated(A, B, l):
    s = QRSystem(A.shape[1] + B.shape[1]); s.add(np.column_stack([A, B]), l)
    x, st = s.solve(); return (None if x is None else x[:A.shape[1]]), st


def main():
  rs = np.random.RandomState(5)
  m = 800
  print('nc=8 commons, cond(A_c) = 10^k ; relative error of x_c vs float128-free LAPACK lstsq reference')
  print('%4s %4s | %-22s | %-12s | %-12s' % ('k', 'nl', 'accumulator qr (default)', 'QR isolated', 'normal'))
  for nl in (0, 3):
      for k in (4, 6, 7, 8):
          rs = np.random.RandomState(5 + k)
          U, _ = np.linalg.qr(rs.randn(m, 8)); V, _ = np.linalg.qr(rs.randn(8, 8))
          A = U @ np.diag(np.logspace(0, -k, 8)) @ V.T
          B = rs.randn(m, nl)
          xt = rs.randn(8 + nl)
          l = np.column_stack([A, B]) @ xt + 1e-10 * rs.randn(m)
          xr = np.linalg.lstsq(np.column_stack([A, B]), l, rcond=None)[0][:8]
          def err(x): return np.nan if x is None else np.linalg.norm(x - xr) / np.linalg.norm(xr)
          xa, sa, pa = accumulator_qr(A, B, l)
          xq, sq = qr_isolated(A, B, l)
          xn, sn = normal_eq(A, B, l)
          print('%4d %4d | %9.2e %-12s | %9.2e %s | %9.2e %s' % (k, nl, err(xa), sa + '/' + pa, err(xq), sq, err(xn), sn))

if __name__=="__main__": main()
