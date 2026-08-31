#!/usr/bin/env python3
"""v6 validation: the VAR path (ROADMAP 3.1) + the spmd day-merge combine.

VAR = vector AR for TEMPORAL cross-axis correlation, the case lsq.crossCorr
flags as 'needs-var' and separableWhiten cannot represent. Everything here is
brute-forced against dense block-Toeplitz linear algebra before a line of
MATLAB is written:

  * companion-form cross-ACV of a VAR(p)  (lsq.varACV)
  * Whittle multichannel Levinson: exact Gamma -> generating (A_k, Sigma)
  * matrix-FIR whitening with exact block startup == chol(C)^-1  (lsq.varWhiten)
  * GLS identities under the true VAR model
  * batch continuation with p-epoch overlap == whole-arc
  * mergeRaw combine is grouping-independent (the spmd tree-reduce claim)

    python3 v6verify.py
"""
import numpy as np
from numpy.linalg import norm, solve, cholesky, inv, lstsq

R = []
def T(n, f):
    try:
        f(); R.append((n, 'PASS', ''))
    except AssertionError as e:
        R.append((n, 'FAIL', str(e)[:70]))
    except Exception as e:
        R.append((n, 'ERROR', f'{type(e).__name__}: {e}'[:70]))

rng = np.random.default_rng(61)


# ---------------------------------------------------------- VAR reference ops
def stable_var(d=3, p=2, seed=3, rho=0.5):
    """A stationary VAR(p) with genuine cross-axis coupling."""
    r = np.random.default_rng(seed)
    A = 0.35 * r.standard_normal((p, d, d))
    A[0] += 0.3 * np.eye(d)
    # scale to spectral radius rho of the companion matrix
    F = companion(A)
    A *= (rho / max(abs(np.linalg.eigvals(F))))
    L = 0.2 * r.standard_normal((d, d))
    Sig = L @ L.T + 0.5 * np.eye(d)
    return A, Sig

def companion(A):
    p, d, _ = A.shape
    F = np.zeros((p * d, p * d))
    F[:d, :] = np.concatenate([A[k] for k in range(p)], axis=1)
    F[d:, :-d] = np.eye((p - 1) * d)
    return F

def var_acv(A, Sig, K):
    """Gamma(0..K), Gamma(j) = E[z_t z_{t-j}'], via the companion Lyapunov
    equation solved by vectorisation, then the Yule-Walker recursion."""
    p, d, _ = A.shape
    F = companion(A)
    Q = np.zeros_like(F); Q[:d, :d] = Sig
    n = F.shape[0]
    GX = solve(np.eye(n * n) - np.kron(F, F), Q.reshape(-1)).reshape(n, n)
    G = [GX[:d, k * d:(k + 1) * d].copy() for k in range(p)]   # Gamma(0..p-1)
    for j in range(p, K + 1):
        G.append(sum(A[k] @ G[j - 1 - k] for k in range(p)))
    return G

def block_toeplitz(G, nE):
    d = G[0].shape[0]
    C = np.zeros((nE * d, nE * d))
    for i in range(nE):
        for j in range(nE):
            k = i - j
            B = G[k] if k >= 0 else G[-k].T
            C[i * d:(i + 1) * d, j * d:(j + 1) * d] = B
    return C

def whittle(G, p):
    """Whittle, the standard symmetric update: Sf -= Kf Sb Kf', Sb -= Kb Sf_old Kb'."""
    d = G[0].shape[0]
    A = np.zeros((0, d, d)); B = np.zeros((0, d, d))
    Sf = G[0].copy(); Sb = G[0].copy()
    for n in range(p):
        D = G[n + 1] - sum(A[k] @ G[n - k] for k in range(n))
        Kf = solve(Sb.T, D.T).T
        Kb = solve(Sf.T, D).T
        Sf_old = Sf.copy()
        An = np.zeros((n + 1, d, d)); Bn = np.zeros((n + 1, d, d))
        An[n] = Kf; Bn[n] = Kb
        for k in range(n):
            An[k] = A[k] - Kf @ B[n - 1 - k]
            Bn[k] = B[k] - Kb @ A[n - 1 - k]
        Sf = Sf - Kf @ Sb @ Kf.T
        Sb = Sb - Kb @ Sf_old @ Kb.T
        A, B = An, Bn
    return A, (Sf + Sf.T) / 2

def var_whiten(A, Sig, G, Z, first=True, state=None):
    """The matrix FIR with exact block startup, column-wise on Z (m x q rows
    epoch-major). Returns whitened Z and the trailing p-epoch state."""
    p, d, _ = A.shape
    m, q = Z.shape
    nE = m // d
    Lf = cholesky(Sig)                          # lower
    Ze = Z.reshape(nE, d, q)
    W = np.empty_like(Ze)
    if first:
        k0 = min(p, nE)
        Cs = block_toeplitz(G, k0)
        Fs = inv(cholesky(Cs))                  # exact stationary startup
        W[:k0] = (Fs @ Ze[:k0].reshape(k0 * d, q)).reshape(k0, d, q)
        start = k0
        hist = Ze[:p].copy() if nE >= p else None
    else:
        start = 0
    for t in range(start, nE):
        acc = Ze[t].copy()
        for k in range(p):
            zk = Ze[t - 1 - k] if t - 1 - k >= 0 else state[p - 1 - (k - t)]
            acc = acc - A[k] @ zk
        W[t] = solve(Lf, acc)
    newstate = Ze[-p:].copy() if nE >= p else None
    return W.reshape(m, q), newstate


# =================== 1. companion ACV: Yule-Walker consistency (deterministic)
def t_var_acv_yule_walker():
    A, Sig = stable_var()
    p = A.shape[0]
    G = var_acv(A, Sig, 12)
    # Gamma(0) identity: G0 = sum A_k G(-k)' ... use  G0 - sum A_k G(k)' = Sig
    G0res = G[0] - sum(A[k] @ G[k + 1].T for k in range(p))
    assert norm(G0res - Sig) / norm(Sig) < 1e-11, f'{norm(G0res-Sig)/norm(Sig):.2e}'
    for j in range(p, 12):
        Gj = sum(A[k] @ G[j - 1 - k] for k in range(p))
        assert norm(G[j] - Gj) / max(norm(G[j]), 1e-300) < 1e-11

def t_var_acv_matches_simulation():
    """Statistical sanity: long simulation reproduces Gamma(0..2) to a few %."""
    A, Sig = stable_var()
    p, d, _ = A.shape
    G = var_acv(A, Sig, 3)
    r = np.random.default_rng(8)
    L = cholesky(Sig); n = 400000
    z = np.zeros((n, d))
    for t in range(p, n):
        z[t] = sum(A[k] @ z[t - 1 - k] for k in range(p)) + L @ r.standard_normal(d)
    z = z[2000:]
    for j in range(3):
        Ge = (z[j:].T @ z[:len(z) - j]) / (len(z) - j)
        e = norm(Ge - G[j]) / norm(G[j])
        assert e < 0.05, f'lag {j}: rel {e:.2e}'


# ============================ 2. Whittle recovers the generator from exact ACV
def t_whittle_recovers_generator():
    A, Sig = stable_var()
    p = A.shape[0]
    G = var_acv(A, Sig, p + 1)
    Ah, Sh = whittle(G, p)
    ea = max(norm(Ah[k] - A[k]) / max(norm(A[k]), 1e-300) for k in range(p))
    es = norm(Sh - Sig) / norm(Sig)
    assert ea < 1e-9 and es < 1e-10, f'A {ea:.2e}  Sig {es:.2e}'

def t_whittle_overfit_gives_zero_extra_lags():
    """Fit order p+2 on an exact VAR(p) ACV: the extra matrices must vanish."""
    A, Sig = stable_var()
    p = A.shape[0]
    G = var_acv(A, Sig, p + 3)
    Ah, _ = whittle(G, p + 2)
    extra = max(norm(Ah[k]) for k in range(p, p + 2))
    ref = norm(A[0])
    assert extra / ref < 1e-8, f'{extra/ref:.2e}'


# ================ 3. the whitener IS the dense block-Cholesky inverse
def t_var_whiten_equals_dense_chol():
    A, Sig = stable_var()
    d = Sig.shape[0]; nE = 60; m = nE * d
    G = var_acv(A, Sig, nE)
    C = block_toeplitz(G, nE)
    Wd = inv(cholesky(C))                       # dense reference operator

    Z = rng.standard_normal((m, 5))
    Y, _ = var_whiten(A, Sig, G, Z, first=True)
    Yd = Wd @ Z
    # F'F = inv(C) exactly, and F == chol-inverse row for row
    e = norm(Y - Yd) / norm(Yd)
    assert e < 1e-10, f'rel {e:.2e}'

def t_var_gls_identities():
    """Whiten with the TRUE VAR model: sigma0 ~ 1, x within formal errors,
    equality with the dense GLS solve."""
    A, Sig = stable_var()
    p, d, _ = A.shape
    nE = 900; m = nE * d; nc = 5
    G = var_acv(A, Sig, nE)
    r = np.random.default_rng(17)
    L = cholesky(Sig)
    z = np.zeros((nE + 300, d))
    for t in range(p, nE + 300):
        z[t] = sum(A[k] @ z[t - 1 - k] for k in range(p)) + L @ r.standard_normal(d)
    e = z[300:].reshape(m)

    Ades = r.standard_normal((m, nc))
    xt = r.standard_normal(nc)
    l = Ades @ xt + e

    J = np.column_stack([Ades, l])
    Jw, _ = var_whiten(A, Sig, G, J, first=True)
    Aw, lw = Jw[:, :nc], Jw[:, nc]
    N = Aw.T @ Aw
    x = solve(N, Aw.T @ lw)
    v = lw - Aw @ x
    s0 = np.sqrt(v @ v / (m - nc))
    assert abs(s0 - 1) < 0.05, f'sigma0 {s0:.4f}'
    sig = np.sqrt(np.diag(inv(N)))
    assert np.all(np.abs(x - xt) < 5 * sig), 'x outside 5 formal sigmas'

    C = block_toeplitz(G, nE)
    xg = solve(Ades.T @ solve(C, Ades), Ades.T @ solve(C, l))
    assert norm(x - xg) / norm(xg) < 1e-9, f'{norm(x-xg)/norm(xg):.2e}'


# ================== 4. batch continuation: p-epoch overlap state == whole-arc
def t_var_whiten_batched_equals_whole():
    A, Sig = stable_var()
    p, d, _ = A.shape
    nE = 80; m = nE * d
    G = var_acv(A, Sig, nE)
    Z = rng.standard_normal((m, 4))

    Yw, _ = var_whiten(A, Sig, G, Z, first=True)

    cut = 33                                    # epochs in the first batch
    Y1, st = var_whiten(A, Sig, G, Z[:cut * d], first=True)
    Y2, _ = var_whiten(A, Sig, G, Z[cut * d:], first=False, state=st)
    Yb = np.vstack([Y1, Y2])
    e = norm(Yb - Yw) / norm(Yw)
    assert e < 1e-12, f'rel {e:.2e}'


# ====== 5. mergeRaw combine is grouping-independent (the spmd tree reduce)
def t_qr_merge_grouping_independent():
    """spmdReduce combines partials in a TREE whose shape depends on the
    worker count. Assert ((a+b)+(c+d)) == (((a+b)+c)+d) in every derived
    quantity for compressed QR partials."""
    d = 9
    r = np.random.default_rng(23)
    parts = []
    for _ in range(4):
        M = r.standard_normal((30, d + 1))
        Rc = np.triu(np.linalg.qr(M, mode='r'))[:d + 1]
        parts.append(Rc)

    def stack(a, b):
        Rn = np.triu(np.linalg.qr(np.vstack([a, b]), mode='r'))
        return Rn[:d + 1]

    seq = stack(stack(stack(parts[0], parts[1]), parts[2]), parts[3])
    tree = stack(stack(parts[0], parts[1]), stack(parts[2], parts[3]))

    def derived(Rc):
        R1 = Rc[:d, :d]; rr = Rc[:d, d]
        ltl = rr @ rr + (Rc[d, d] ** 2 if Rc.shape[0] > d else 0.0)
        return R1.T @ R1, R1.T @ rr, ltl

    for a, b in zip(derived(seq), derived(tree)):
        e = norm(np.atleast_1d(a - b)) / max(norm(np.atleast_1d(a)), 1e-300)
        assert e < 5e-13, f'rel {e:.2e}'


T('1. VAR companion ACV satisfies Yule-Walker exactly',
  t_var_acv_yule_walker)
T('2. VAR companion ACV matches a long simulation (statistical)',
  t_var_acv_matches_simulation)
T('3. Whittle on the exact ACV recovers the generating (A, Sigma)',
  t_whittle_recovers_generator)
T('4. Whittle overfit: extra coefficient matrices vanish',
  t_whittle_overfit_gives_zero_extra_lags)
T('5. matrix-FIR + exact block startup == dense chol(C)^-1, row for row',
  t_var_whiten_equals_dense_chol)
T('6. GLS identities under the true VAR model (sigma0, x, dense equality)',
  t_var_gls_identities)
T('7. batched whitening with p-epoch state == whole-arc',
  t_var_whiten_batched_equals_whole)
T('8. QR mergeRaw is grouping-independent (spmd tree reduce is safe)',
  t_qr_merge_grouping_independent)

n = sum(1 for _, s, _ in R if s == 'PASS')
for name, s, msg in R:
    print(f'  [{s:5s}] {name}' + (f'  -- {msg}' if msg else ''))
print(f'  {n}/{len(R)} passed')
