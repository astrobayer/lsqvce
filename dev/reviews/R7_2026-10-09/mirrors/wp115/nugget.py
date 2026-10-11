"""WP-115 part 3b mirror: the additive nugget and the weighted redundancy.

The robust factor f of an epoch becomes an additive nugget c = Sigma_ii (1/f - 1)
on the temporal covariance (same marginal variance as the Koch fold, the
neighbours' covariances untouched): Sigma = T + diag(c), T the law (Toeplitz).

  1. whitener  X = G F A, F the causal law whitener (F'F = T^-1), G = I - Q (I - S) Q',
               Q R = qr(F U), S = (I + R diag(c) R')^-1/2: X'X == A' Sigma^-1 A
  2. dS        lag sums of Sigma^-1 = T^-1 - H K^-1 H' (Woodbury; H = T^-1 U,
               K = diag(1/c) + U' T^-1 U) == the dense inverse's
  3. s^w       E_T[Omega_n] = lagsum_n(Rbar - (Rbar U) diag(c) (Rbar U)')  (Rbar U = Rbar(:, z))
               == lagsum_n of Sigma^-1 P T P' Sigma^-1, the exact expectation
               under the law-only noise (clean data down-weighted)
  4. the robust joint loop of t3_robust_level.py with nugget + s^w: level ~ 1
"""
import sys, os
import numpy as np
from scipy.linalg import toeplitz, cholesky, solve_triangular, qr, eigh
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'dom_vce'))
from kern import dctI, redundancyLags, columnLagSums

def lagsums(M, nL):
    """d(n+1) = sum_i M(i+n, i)."""
    return np.array([np.trace(M, offset=-n) for n in range(nL)])

def dense_s_om(Ab, W, UiZ, eb, nL):
    """redundancyLags for a general (non-Toeplitz) W: the kernel's fallback path."""
    Wi = np.linalg.inv(W)
    Si = Wi @ Wi.T
    Rt = Wi @ (Ab @ UiZ)
    s = lagsums(Si, nL) - columnLagSums(Rt, nL); s[1:] *= 2
    om = None
    if eb is not None:
        om = columnLagSums(Wi @ eb, nL); om[1:] *= 2
    return s, om

def nugget_whiten(T, c, A):
    """X = G F A (mirror of covModel>nuggetRows)."""
    m = T.shape[0]
    L = cholesky(T, lower=True)
    FA = solve_triangular(L, A, lower=True)
    z = np.flatnonzero(c > 0)
    if z.size == 0:
        return FA
    E = np.zeros((m, z.size)); E[z, np.arange(z.size)] = 1
    FU = solve_triangular(L, E, lower=True)
    Q, Rq = qr(FU, mode='economic')
    K = np.diag(1 / c[z]) + Rq.T @ Rq                 # O(1) entries, whatever c is
    B = Rq @ np.linalg.solve(K, Rq.T)                 # eigenvalues beta in [0, 1)
    beta, V = eigh((B + B.T) / 2)
    g = 1 - np.sqrt(np.clip(1 - beta, 0, None))       # G = I - Q V diag(g) V' Q'
    return FA - Q @ (V @ (g[:, None] * (V.T @ (Q.T @ FA))))

def nugget_dS(T, c, nL):
    """Lag sums of (T + diag c)^-1 via the law's factor and a k x k Woodbury."""
    m = T.shape[0]
    dT = lagsums(np.linalg.inv(T), nL)            # (the kernel's Gohberg-Semencul path)
    z = np.flatnonzero(c > 0)
    WT = cholesky(T)                              # upper, T = WT' WT
    E = np.zeros((m, z.size)); E[z, np.arange(z.size)] = 1
    Y = solve_triangular(WT, E, trans='T')        # WT^-T U
    H = solve_triangular(WT, Y)                   # T^-1 U
    K = np.diag(1 / c[z]) + Y.T @ Y
    LK = cholesky(K, lower=True)
    G = solve_triangular(LK, H.T, lower=True).T   # H LK^-T
    return dT - columnLagSums(G, nL)

rng = np.random.default_rng(115)
ok = True
m, nL = 200, 40
t = np.arange(m) * 5.0
A = np.c_[np.ones(m), (t - t.mean()) / np.ptp(t), np.sin(2*np.pi*t/1500.), np.cos(2*np.pi*t/700.)]
for phi in (0.5, 0.9):
    T = 1.7 * toeplitz(phi ** np.arange(m))
    f = np.ones(m)
    z = rng.choice(m, 30, replace=False)
    f[z[:10]] = 1e-12                              # zeroed (the loader's floor)
    f[z[10:]] = rng.uniform(0.05, 0.95, 20)        # partial
    c = np.diag(T) * (1 / f - 1)
    Sig = T + np.diag(c)
    Si = np.linalg.inv(Sig)
    # 1. whitener
    X = nugget_whiten(T, c, A)
    e1 = np.abs(X.T @ X - A.T @ Si @ A).max() / np.abs(A.T @ Si @ A).max()
    # 2. dS
    e2 = np.abs(nugget_dS(T, c, nL) - lagsums(Si, nL)).max() / abs(np.trace(Si))
    # 3. weighted redundancy, exact hat (UiZ = U^-1, Z = P)
    N = A.T @ Si @ A
    U = cholesky(N)
    W = cholesky(Sig)
    Abar = solve_triangular(W, A, trans='T')
    s, _ = dense_s_om(Abar, W, np.linalg.inv(U), None, nL)
    Rbar = Si - Si @ A @ np.linalg.solve(N, A.T @ Si)
    RU = Rbar[:, c > 0] * np.sqrt(c[c > 0])[None, :]
    d = columnLagSums(RU, nL)
    sw = s.copy(); sw[0] -= d[0]; sw[1:] -= 2 * d[1:]
    P = np.eye(m) - A @ np.linalg.solve(N, A.T @ Si)
    Etrue = Si @ P @ T @ P.T @ Si
    ex = lagsums(Etrue, nL); ex[1:] *= 2
    e3 = np.abs(sw - ex).max() / abs(ex[0])
    e3m = np.abs(s - lagsums(Rbar, nL) * np.r_[1, 2*np.ones(nL-1)]).max() / abs(s[0])
    print(f'phi {phi}: whitener {e1:.1e}  dS {e2:.1e}  s(model) vs dense {e3m:.1e}  s^w vs E_T[Omega] {e3:.1e}'
          f'   (s_0 model {s[0]:.2f}, weighted {sw[0]:.2f})')
    ok &= max(e1, e2, e3, e3m) < 1e-9

# 4. the robust joint loop (t3_robust_level.py) with the nugget and s^w
def igg3(u, k0=1.5, k1=3.0):
    a = np.abs(u); w = np.ones_like(a)
    mid = (a > k0) & (a <= k1); w[mid] = (k0/a[mid])*((k1-a[mid])/(k1-k0))**2; w[a > k1] = 0
    return w
mA, nArc, nc, nLag = 256, 9, 4, 60
tt = np.arange(nArc*mA)*5.0; w0 = 2*np.pi/(0.5*86400)
Afull = np.c_[np.sin(w0*tt), np.cos(w0*tt), np.sin(2*w0*tt), np.cos(2*w0*tt)]
def run(seed, mode, phi=0.5, passes=6):
    acv = phi**np.arange(mA)/(1-phi**2)
    Rt = toeplitz(acv); Wt = cholesky(Rt)
    r = np.random.default_rng(seed)
    e = np.concatenate([Wt.T@r.normal(size=mA) for _ in range(nArc)])
    l = Afull@np.array([.8, -.45, .3, .15]) + e
    level = 1.0; w = np.ones(len(l))
    for it in range(passes):
        Sigs = []
        N = np.zeros((nc, nc)); nv = np.zeros(nc)
        for k in range(nArc):
            rows = slice(k*mA, (k+1)*mA)
            if mode == 'fold':                    # t3: descale by sqrt(w), law level*Rt
                d = np.sqrt(np.maximum(w[rows], 1e-12)); Sig = level*Rt
                cvec = np.zeros(mA)
            else:
                d = np.ones(mA)
                cvec = level*acv[0]*(1/np.maximum(w[rows], 1e-12) - 1); Sig = level*Rt + np.diag(cvec)
            Ws = cholesky(Sig)
            Ab = solve_triangular(Ws, Afull[rows]*d[:, None], trans='T'); lb = solve_triangular(Ws, l[rows]*d, trans='T')
            N += Ab.T@Ab; nv += Ab.T@lb; Sigs.append((Ws, cvec, rows, d))
        x = np.linalg.solve(N, nv)
        UiZ = np.linalg.inv(cholesky(N))
        a2 = []
        for (Ws, cvec, rows, d) in Sigs:
            if mode == 'hard':
                cvec = np.where(w[rows] == 0, cvec, 0.0)
                Ws = cholesky(level*Rt + np.diag(cvec))
            Ab = solve_triangular(Ws, Afull[rows]*d[:, None], trans='T')
            eb = solve_triangular(Ws, (l[rows] - Afull[rows]@x)*d, trans='T')
            s, om = dense_s_om(Ab, Ws, UiZ, eb, nLag)
            if mode == 'nugget+sw' and (cvec > 0).any():
                z = cvec > 0
                Si_z = solve_triangular(Ws, solve_triangular(Ws, np.eye(mA)[:, z], trans='T'))
                Rtl = solve_triangular(Ws, Ab @ UiZ)
                RU = (Si_z - Rtl @ (Rtl[z, :].T)) * np.sqrt(cvec[z])[None, :]
                dd = columnLagSums(RU, nLag); s = s.copy(); s[0] -= dd[0]; s[1:] -= 2*dd[1:]
            C = acv[:nLag]/acv[0]
            a2.append(dctI(level*C*om).sum()/dctI(level*C*s).sum())
        level = np.mean(np.array(a2)*level)
        v = l - Afull@x
        sc = 1.4826*np.median(np.abs(v - np.median(v)))
        w = igg3(v/sc)
    return level
MODES = tuple(sys.argv[1:]) or ('fold', 'nugget', 'nugget+sw', 'hard')
print('robust joint loop (IGG3), clean AR(1), law declared, 6 passes; level per seed (truth 1).')
print('IGG3 under N(0,1): E[f] = 0.932, E[f u^2] = 0.676, E[f u^2]/E[f] = 0.726 (the selection factor)')
for phi in (0.5, 0.0):
    for mode in MODES:
        L = np.array([run(sd, mode, phi) for sd in range(8)])
        print(f'  phi {phi}  {mode:10s} mean {L.mean():.4f}  se {L.std(ddof=1)/np.sqrt(len(L)):.4f}   {np.round(L, 3)}')
print('ALL OK' if ok else 'FAILED')
sys.exit(0 if ok else 1)
