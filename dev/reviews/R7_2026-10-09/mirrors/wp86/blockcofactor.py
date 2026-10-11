"""WP-86 part 2: mirror of lsqvce.engine.BlockCofactor against a DENSE reference.

The diagnostics used to factor the FULL block-arrow normal matrix (u = nc + sum nl,
13 729 at a GRAVFIRE month: 1.5 GB per dense u x u copy, and diagnose kept four).
BlockCofactor answers every question diagnose / snoop / ellipse / report ask of
Qxx = inv(N_full) from the blocks the accumulator already has:

    Rc     nc x nc upper, Rc'Rc = N_red  (reduced commons, prior folded in)
    per system s:  H_s (nl x nl, H'H = pinv(Nll_s)),  G_s = pinv(Nll_s) Nlc_s,  Nll_s

    half(A)    W with W'W = A Qxx A'      (the dense code's  R' \\ A')
    apply(A)   Qxx A'
    diagonal   diag(Qxx)
    block(idx) Qxx(idx, idx)
    maxCorrelation, condLowerBound, downdate (snoop's one-row removal)

Run: python3 blockcofactor.py   -> every check prints its max error and PASS/FAIL."""
import numpy as np
from scipy.linalg import solve_triangular, qr

def local_half(Nll):
    """H with H'H = pinv(Nll): inv(L) from chol when PD, else the eigen pseudo-inverse."""
    nl = Nll.shape[0]
    if nl == 0:
        return np.zeros((0, 0))
    S = (Nll + Nll.T) / 2
    if 1 / np.linalg.cond(S, 1) >= np.finfo(float).eps:           # the accumulator's own switch (rcond < eps)
        try:
            L = np.linalg.cholesky(S)
            return solve_triangular(L, np.eye(nl), lower=True)      # inv(L): inv(L)'inv(L) = inv(LL')
        except np.linalg.LinAlgError:
            pass
    if True:
        dv, V = np.linalg.eigh((Nll + Nll.T) / 2)
        keep = (dv > nl * np.finfo(float).eps * dv.max()) & (dv > 0)
        H = np.zeros((nl, nl))
        H[:keep.sum(), :] = (V[:, keep] / np.sqrt(dv[keep])).T
        return H

class BlockCofactor:
    def __init__(self, Rc, systems, u):
        self.Rc = Rc; self.nc = Rc.shape[0]; self.u = u
        self.sys = []
        for (cols, Nll, Nlc) in systems:
            H = local_half(Nll)
            G = H.T @ (H @ Nlc)
            self.sys.append(dict(cols=np.asarray(cols, int), Nll=Nll.copy(), H=H, G=G))

    def _split(self, A):
        Z = A[:, :self.nc].T.copy()
        Bs = []
        for s in self.sys:
            B = A[:, s['cols']].T
            touched = B.size > 0 and np.any(B != 0)
            if touched:
                Z -= s['G'].T @ B
            Bs.append(B if touched else None)
        return Z, Bs

    def half(self, A):
        Z, Bs = self._split(A)
        W = np.zeros((self.u, A.shape[0]))
        W[:self.nc] = solve_triangular(self.Rc.T, Z, lower=True)    # Rc' \ Z
        for s, B in zip(self.sys, Bs):
            if B is not None:
                W[s['cols']] = s['H'] @ B
        return W

    def apply(self, A):
        Z, Bs = self._split(A)
        Y = solve_triangular(self.Rc, solve_triangular(self.Rc.T, Z, lower=True))   # Sxx Z
        out = np.zeros((self.u, A.shape[0]))
        out[:self.nc] = Y
        for s, B in zip(self.sys, Bs):
            out[s['cols']] = -s['G'] @ Y
            if B is not None:
                out[s['cols']] += s['H'].T @ (s['H'] @ B)
        return out

    def diagonal(self):
        Ri = solve_triangular(self.Rc, np.eye(self.nc))             # Sxx = Ri Ri'
        d = np.zeros(self.u)
        d[:self.nc] = np.sum(Ri**2, axis=1)
        for s in self.sys:
            d[s['cols']] = np.sum(s['H']**2, axis=0) + np.sum((s['G'] @ Ri)**2, axis=1)
        return d

    def block(self, idx):
        E = np.zeros((len(idx), self.u)); E[np.arange(len(idx)), idx] = 1
        W = self.half(E)
        return W.T @ W

    def maxCorrelation(self, panel=3):
        sig = np.sqrt(np.maximum(self.diagonal(), 0))
        best, pair = 0.0, (0, 0)
        for j0 in range(0, self.u, panel):
            J = np.arange(j0, min(j0 + panel, self.u))
            E = np.zeros((len(J), self.u)); E[np.arange(len(J)), J] = 1
            QJ = self.apply(E)                                       # Qxx(:, J)
            with np.errstate(divide='ignore', invalid='ignore'):
                C = QJ / np.outer(sig, sig[J])
            C[J, np.arange(len(J))] = 0
            C[~np.isfinite(C)] = 0
            k = np.argmax(np.abs(C))
            i, jj = np.unravel_index(k, C.shape)
            if abs(C[i, jj]) > best:
                best, pair = abs(C[i, jj]), (i, J[jj])
        return best, pair

    def times(self, V):
        """N_full V from the blocks (for the conditioning bound)."""
        Vc = V[:self.nc]
        out = np.zeros_like(V)
        out[:self.nc] = self.Rc.T @ (self.Rc @ Vc)
        for s in self.sys:
            T = s['Nll'] @ (s['G'] @ Vc + V[s['cols']])
            out[:self.nc] += s['G'].T @ T
            out[s['cols']] = T
        return out

    def condLowerBound(self, nIter=30):
        dQ = self.diagonal()
        dN = np.zeros(self.u)
        dN[:self.nc] = np.sum(self.Rc**2, axis=0)
        for s in self.sys:
            dN[:self.nc] += np.sum(s['G'] * (s['Nll'] @ s['G']), axis=0)
            dN[s['cols']] = np.diag(s['Nll'])
        lam, mu = dN.max(), dQ.max()
        for op, which in ((self.times, 'N'), (lambda v: self.apply(v.T), 'Q')):
            v = np.ones((self.u, 1)) / np.sqrt(self.u)
            for _ in range(nIter):
                w = op(v); rq = (v.T @ w).item()
                if which == 'N': lam = max(lam, rq)
                else:            mu = max(mu, rq)
                v = w / np.linalg.norm(w)
        return lam * mu

    def downdate(self, a, p):
        """Remove one row a (1 x u, at most ONE system's locals) of weight p."""
        a = a.reshape(1, -1)
        Z, Bs = self._split(a)
        z = Z[:, 0]
        touched = [k for k, B in enumerate(Bs) if B is not None]
        assert len(touched) <= 1
        ptil = p
        if touched:
            s = self.sys[touched[0]]; b = Bs[touched[0]][:, 0]
            c = s['H'].T @ (s['H'] @ b)
            d = 1 - p * b @ c
            if d <= 1e3 * np.finfo(float).eps: return False
            ptil = p / d
        Rnew, ok = cholupdate_minus(self.Rc, np.sqrt(ptil) * z)
        if not ok: return False
        self.Rc = Rnew
        if touched:
            s['Nll'] = s['Nll'] - p * np.outer(b, b)
            s['H'] = local_half(s['Nll'])
            s['G'] = s['G'] - ptil * np.outer(c, z)
        return True

def cholupdate_minus(R, x):
    """R'R - x x' (R upper), LINPACK-style hyperbolic rotations; ok=False when not PD."""
    R = R.copy(); x = x.copy(); n = len(x)
    for k in range(n):
        r2 = R[k, k]**2 - x[k]**2
        if r2 <= 0: return R, False
        r = np.sqrt(r2); c = r / R[k, k]; s = x[k] / R[k, k]
        R[k, k] = r
        if k + 1 < n:
            R[k, k+1:] = (R[k, k+1:] - s * x[k+1:]) / c
            x[k+1:] = c * x[k+1:] - s * R[k, k+1:]
    return R, True

# ------------------------------------------------------------------ the problem
rng = np.random.default_rng(86)
nc, nls = 5, [2, 3, 0, 4]
u = nc + sum(nls)
offs = np.cumsum([nc] + nls)[:-1]
rows, sysOf = [], []
for s, nl in enumerate(nls):
    for _ in range(12 + 3 * nl):
        r = np.zeros(u); r[:nc] = rng.normal(size=nc)
        r[offs[s]:offs[s] + nl] = rng.normal(size=nl)
        rows.append(r); sysOf.append(s)
A = np.array(rows)
K = rng.normal(size=(2, nc)) / 0.5                     # a common prior, K/sigma
Afull = np.vstack([A, np.c_[K, np.zeros((2, u - nc))]])
N = Afull.T @ Afull
Qd = np.linalg.inv(N)

def blocks_of(Arows, sysOf):
    Nred = K.T @ K
    syst = []
    for s, nl in enumerate(nls):
        idx = [i for i, t in enumerate(sysOf) if t == s]
        As = Arows[idx][:, :nc]; Bs = Arows[idx][:, offs[s]:offs[s] + nl]
        Nll, Nlc, Ncc = Bs.T @ Bs, Bs.T @ As, As.T @ As
        Nred += Ncc - (Nlc.T @ np.linalg.solve(Nll, Nlc) if nl else 0)
        syst.append((np.arange(offs[s], offs[s] + nl), Nll, Nlc))
    return Nred, syst

Nred, syst = blocks_of(A, sysOf)
# Rc as BlockAccumulator.commonFactor builds it on the QR path: [R1 ; K/sigma] -> qr
R1 = np.linalg.cholesky(Nred - K.T @ K).T
Rc = np.triu(qr(np.vstack([R1, K]), mode='economic')[1])
Q = BlockCofactor(Rc, syst, u)

def check(name, err, tol=1e-9):
    print(f'{name:<34s} {err:9.2e}  {"PASS" if err < tol else "FAIL"}')

check('Rc from QR fold == chol(N_red)', np.abs(Rc.T @ Rc - Nred).max() / np.abs(Nred).max())
check('half: W\'W = A Qxx A\'', np.abs(Q.half(A).T @ Q.half(A) - A @ Qd @ A.T).max())
check('apply: Qxx A\'', np.abs(Q.apply(A) - Qd @ A.T).max())
check('diagonal', np.abs(Q.diagonal() - np.diag(Qd)).max())
idx = [1, 6, 9, 13]
check('block across systems', np.abs(Q.block(idx) - Qd[np.ix_(idx, idx)]).max())
h = np.sum(Q.half(Afull)**2, axis=0)
check('sum h = u (incl. prior rows)', abs(h.sum() - u))
sd = np.sqrt(np.diag(Qd)); C = Qd / np.outer(sd, sd); np.fill_diagonal(C, 0)
cm, pr = Q.maxCorrelation()
check('maxCorrelation value', abs(cm - np.abs(C).max()))
print('   pair', pr, 'dense', np.unravel_index(np.argmax(np.abs(C)), C.shape))
ev = np.linalg.eigvalsh(N); kb = Q.condLowerBound()
print(f'cond(N) exact {ev[-1]/ev[0]:.6e}  lower bound {kb:.6e}  ratio {kb/(ev[-1]/ev[0]):.6f}'
      f'  {"PASS" if kb <= ev[-1]/ev[0]*(1+1e-12) else "FAIL"}')

# downdate: remove row i of system 1 (weight p = 1, whitened); then a common-only row (system 2)
l = A @ rng.normal(size=u) + rng.normal(size=A.shape[0])
lf = np.r_[l, K @ np.zeros(nc)]
x = np.linalg.solve(N, Afull.T @ lf)
for i in (20, [k for k, t in enumerate(sysOf) if t == 2][3]):
    g = A[i]; v = g @ x - lf[i]; p = 1.0
    ok = Q.downdate(g, p)
    N = N - p * np.outer(g, g); Qd = np.linalg.inv(N)
    keep = np.ones(len(lf), bool); keep[i] = False
    lf_i = lf.copy()
    x_new = x + p * v * Q.apply(g.reshape(1, -1))[:, 0]
    xd = np.linalg.solve(N, Afull.T @ (lf * keep))
    Afull = Afull * keep[:, None]; lf = lf * keep
    check(f'downdate row {i}: ok', 0 if ok else 1)
    check(f'downdate row {i}: apply', np.abs(Q.apply(A) - Qd @ A.T).max())
    check(f'downdate row {i}: diagonal', np.abs(Q.diagonal() - np.diag(Qd)).max())
    check(f'downdate row {i}: x update', np.abs(x_new - xd).max())
    x = xd

# pinv fallback: a singular local block (two equal local columns in system 0)
A2 = A.copy(); A2[:, offs[0] + 1] = A2[:, offs[0]]
Nred2, syst2 = None, []
Nred2 = K.T @ K
for s, nl in enumerate(nls):
    ids = [i for i, t in enumerate(sysOf) if t == s]
    As = A2[ids][:, :nc]; Bs = A2[ids][:, offs[s]:offs[s] + nl]
    Nll, Nlc = Bs.T @ Bs, Bs.T @ As
    Nred2 += As.T @ As - (Nlc.T @ np.linalg.pinv(Nll) @ Nlc if nl else 0)
    syst2.append((np.arange(offs[s], offs[s] + nl), Nll, Nlc))
Q2 = BlockCofactor(np.linalg.cholesky(Nred2).T, syst2, u)
H = Q2.sys[0]['H']
check('pinv fallback: H\'H = pinv(Nll)', np.abs(H.T @ H - np.linalg.pinv(syst2[0][1])).max(), 1e-6)
Af2 = np.vstack([A2, np.c_[K, np.zeros((2, u - nc))]])
check('pinv fallback: hat = A pinv(N) A\'',
      np.abs(Q2.half(A2).T @ Q2.half(A2) - A2 @ np.linalg.pinv(Af2.T @ Af2) @ A2.T).max(), 1e-6)
