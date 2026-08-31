import numpy as np
from numpy.linalg import solve, inv, norm, qr, lstsq
from scipy.linalg import cholesky, solve_triangular, toeplitz, block_diag
from scipy.signal import lfilter
import sys
sys.path.insert(0, '/home/claude')
from v3core import *

R = []
def T(n, f):
    try:
        f(); R.append((n, 'PASS', ''))
    except AssertionError as e:
        R.append((n, 'FAIL', str(e)[:64]))
    except Exception as e:
        R.append((n, 'ERROR', f'{type(e).__name__}: {e}'[:64]))

rng = np.random.default_rng(11)

# ============================================ 1. joint leverage == brute force
def t_leverage_exact():
    nc, nl, nSys = 8, 4, 3
    blocks = []
    for s in range(nSys):
        m = 200 + 30*s
        blocks.append((rng.standard_normal((m, nc)), rng.standard_normal((m, nl))))

    # brute-force joint hat matrix
    ntot = nc + nSys*nl
    rows = []
    for s, (A, B) in enumerate(blocks):
        m = A.shape[0]
        blk = np.zeros((m, ntot))
        blk[:, :nc] = A
        blk[:, nc+s*nl:nc+(s+1)*nl] = B
        rows.append(blk)
    Aw = np.vstack(rows)
    N = Aw.T @ Aw
    H = Aw @ inv(N) @ Aw.T
    h_true = np.diag(H)

    # via the block-arrow route (what the toolbox will do)
    Ncc_red = np.zeros((nc, nc))
    sysinfo = []
    for A, B in blocks:
        Ncc = A.T@A; Nlc = B.T@A; Nll = B.T@B
        Ncc_red += Ncc - Nlc.T @ solve(Nll, Nlc)
        sysinfo.append((Nll, Nlc))

    h_est = []
    for (A, B), (Nll, Nlc) in zip(blocks, sysinfo):
        Lc, Ll, G, SG, GSG = leverage_blocks(Ncc_red, Nll, Nlc)
        h_est.append(joint_leverage(A, B, Lc, Ll, G, SG, GSG))
    h_est = np.concatenate(h_est)

    e = np.abs(h_est - h_true).max()
    assert e < 1e-9, f'max abs err {e:.2e}'
    # and the identity sum(h_ii) = u
    assert abs(h_est.sum() - ntot) < 1e-8, f'sum(h)={h_est.sum():.6f} != u={ntot}'

def t_redundancy_identity():
    """sum(r_i) = m - u  is the strongest global check."""
    nc, nl, nSys = 6, 3, 4
    blocks = [(rng.standard_normal((150, nc)), rng.standard_normal((150, nl)))
              for _ in range(nSys)]
    Ncc_red = np.zeros((nc, nc)); info = []
    for A, B in blocks:
        Nlc = B.T@A; Nll = B.T@B
        Ncc_red += A.T@A - Nlc.T @ solve(Nll, Nlc)
        info.append((Nll, Nlc))
    h = []
    for (A, B), (Nll, Nlc) in zip(blocks, info):
        Lc, Ll, G, SG, GSG = leverage_blocks(Ncc_red, Nll, Nlc)
        h.append(joint_leverage(A, B, Lc, Ll, G, SG, GSG))
    h = np.concatenate(h)
    m = sum(A.shape[0] for A, _ in blocks)
    u = nc + nSys*nl
    r = 1 - h
    assert abs(r.sum() - (m-u)) < 1e-7, f'sum(r)={r.sum():.4f} != m-u={m-u}'
    assert (h > -1e-10).all() and (h < 1+1e-10).all(), 'h_ii outside [0,1]'

# ============================================ 2. common-part regularisation
def t_prior_normals_and_solve():
    nc, nk = 10, 10
    K = np.diag(1.0/np.arange(1, nc+1)**2)          # Kaula-like
    x0 = np.zeros(nc)
    pri = CommonPrior(K, x0, sigma=0.5)
    Nk, nk_ = pri.normals()
    w2 = 1/0.5**2
    assert norm(Nk - w2*(K.T@K)) < 1e-12
    assert norm(nk_ - w2*(K.T@x0)) < 1e-12

def t_prior_redundancy():
    """r_K = nk - trace(Sxx K'K)/sig^2, and DoF absorbed in [0, nk]."""
    nc = 12
    A = rng.standard_normal((200, nc))
    K = np.eye(nc)
    for sig in [0.01, 1.0, 100.0]:
        pri = CommonPrior(K, sigma=sig)
        Nk, _ = pri.normals()
        N = A.T@A + Nk
        Lc = cholesky(N, lower=True)
        rK, absorbed = pri.redundancy(Lc)
        # brute force
        tr = np.trace(inv(N) @ (K.T@K)) / sig**2
        assert abs(absorbed - tr) < 1e-9, f'sig={sig}: {absorbed:.6f} vs {tr:.6f}'
        assert -1e-8 <= absorbed <= nc + 1e-8, f'absorbed={absorbed} out of [0,{nc}]'
    # strong prior absorbs almost everything; weak prior absorbs almost nothing
    a_strong = CommonPrior(K, sigma=0.01).redundancy(
        cholesky(A.T@A + CommonPrior(K, sigma=0.01).normals()[0], lower=True))[1]
    a_weak = CommonPrior(K, sigma=100.).redundancy(
        cholesky(A.T@A + CommonPrior(K, sigma=100.).normals()[0], lower=True))[1]
    assert a_strong > 0.9*nc, f'strong prior absorbed only {a_strong:.2f}/{nc}'
    assert a_weak < 0.1*nc, f'weak prior absorbed {a_weak:.2f}/{nc}'

def t_prior_regularises_rank_deficient():
    """Rank-deficient common block: the prior must make it solvable."""
    nc = 8
    A = rng.standard_normal((50, nc))
    A[:, -1] = A[:, -2]                             # exactly collinear -> singular
    N = A.T@A
    try:
        cholesky(N, lower=True); assert False, 'should be singular'
    except Exception:
        pass
    pri = CommonPrior(np.eye(nc), sigma=1.0)
    Nk, nk_ = pri.normals()
    Lc = cholesky(N + Nk, lower=True)               # must succeed
    assert np.all(np.isfinite(Lc))

def t_prior_basis_agnostic():
    """K can be anything: identity, Kaula, mascon Laplacian, subset."""
    nc = 16
    A = rng.standard_normal((100, nc))
    # mascon smoothness: 1-D graph Laplacian
    L = np.zeros((nc-1, nc))
    for i in range(nc-1):
        L[i, i] = -1; L[i, i+1] = 1
    # subset: regularise only the last 4
    Ksub = np.zeros((4, nc)); Ksub[:, -4:] = np.eye(4)
    for K, name in [(np.eye(nc), 'identity'),
                    (np.diag(1/np.arange(1, nc+1)**2), 'kaula'),
                    (L, 'laplacian'),
                    (Ksub, 'subset')]:
        pri = CommonPrior(K, sigma=1.0)
        Nk, nk_ = pri.normals()
        Lc = cholesky(A.T@A + Nk, lower=True)
        rK, absorbed = pri.redundancy(Lc)
        assert np.isfinite(rK) and np.isfinite(absorbed), name
        assert -1e-8 <= absorbed <= pri.nk + 1e-8, f'{name}: absorbed={absorbed}'

# ============================================ 3. separable multi-axis
def t_separable_equals_dense():
    """The separable whitener must equal the dense Cholesky of the full C."""
    d, nE = 3, 60
    m = nE*d
    a = np.array([1, -1.2, 0.4]); s = 0.8; p = len(a)-1
    r = s**2 * arma_acv(a, [1.], nE)
    Rt = toeplitz(r/r[0])                            # temporal CORRELATION (unit diag)
    sig = 0.5 + rng.random((nE, d))                  # per-epoch, per-axis sigma

    # full covariance:  C = D (Rt kron I_d) D ,  scaled so temporal var = r[0]
    Ckron = np.kron(Rt, np.eye(d)) * r[0]
    D = np.diag(sig.ravel())
    C = D @ Ckron @ D

    A = rng.standard_normal((m, 5)); B = rng.standard_normal((m, 2))
    l = rng.standard_normal(m)

    Wd = solve_triangular(cholesky(C, lower=True), np.eye(m), lower=True)
    Aw_d, Bw_d, lw_d = Wd@A, Wd@B, Wd@l

    Aw, Bw, lw = separable_whiten(A, B, l, a, s, r, sig, d, True)

    # normals must agree (the factors may differ by an orthogonal rotation)
    e1 = norm(Aw.T@Aw - Aw_d.T@Aw_d)/norm(Aw_d.T@Aw_d)
    e2 = norm(Aw.T@lw - Aw_d.T@lw_d)/norm(Aw_d.T@lw_d)
    assert e1 < 1e-8, f'A normals differ by {e1:.2e}'
    assert e2 < 1e-8, f'rhs differs by {e2:.2e}'

def t_separable_full_epoch_cov():
    """Per-epoch FULL dxd covariance (the realistic PPP case)."""
    d, nE = 3, 40
    m = nE*d
    a = np.array([1, -0.9]); s = 1.0
    r = s**2 * arma_acv(a, [1.], nE)
    Rt = toeplitz(r/r[0])

    # per-epoch 3x3 covariances: vertical worse, correlated with horizontal
    Ce = np.zeros((nE, d, d))
    for k in range(nE):
        Q = rng.standard_normal((d, d))*0.2 + np.diag([1., 1., 2.5])
        Ce[k] = Q@Q.T
    Dblk = block_diag(*[cholesky(Ce[k], lower=True) for k in range(nE)])
    C = Dblk @ (np.kron(Rt, np.eye(d))*r[0]) @ Dblk.T

    A = rng.standard_normal((m, 4)); B = rng.standard_normal((m, 2))
    l = rng.standard_normal(m)
    Wd = solve_triangular(cholesky(C, lower=True), np.eye(m), lower=True)

    Aw, Bw, lw = separable_whiten(A, B, l, a, s, r, Ce, d, True)
    e = norm(Aw.T@Aw - (Wd@A).T@(Wd@A))/norm((Wd@A).T@(Wd@A))
    assert e < 1e-8, f'{e:.2e}'

def t_separable_cost_is_scalar():
    """Multi-axis must cost the same as scalar: O(m*p), no d^2 penalty."""
    import time
    a = np.array([1, -1.5, 0.7]); s = 1.0
    for d in [1, 3]:
        nE = 20000
        m = nE*d
        r = s**2*arma_acv(a, [1.], 8)
        sig = np.ones((nE, d))
        A = rng.standard_normal((m, 10)); B = rng.standard_normal((m, 3))
        l = rng.standard_normal(m)
        t0 = time.perf_counter()
        separable_whiten(A, B, l, a, s, r, sig, d, True)
        dt = time.perf_counter()-t0
        per_row = dt/m
        if d == 1: base = per_row
    # per-row cost must not blow up with d
    assert per_row < 4*base, f'multi-axis per-row cost {per_row/base:.1f}x scalar'

# ============================================ 4. cross-correlation diagnostic
def t_crosscorr_detects_ok():
    d, nE = 3, 4000
    w = rng.standard_normal(nE*d)                    # perfectly white
    c = cross_corr(w, d)
    assert c['verdict'] == 'ok', f"{c['verdict']}: {c['advice']}"

def t_crosscorr_detects_cross_axis():
    d, nE = 3, 4000
    Z = rng.standard_normal((nE, d))
    M = np.array([[1, .6, .3], [.6, 1, .4], [.3, .4, 1.]])
    W = Z @ cholesky(M, lower=True).T                # cross-axis at zero lag
    c = cross_corr(W.ravel(), d)
    assert c['verdict'] == 'cross-axis-zero-lag', f"{c['verdict']}"

def t_crosscorr_detects_ar_too_low():
    d, nE = 3, 4000
    Z = rng.standard_normal((nE, d))
    W = lfilter([1.], [1, -0.6], Z, axis=0)          # leftover TEMPORAL colour
    W = W/W.std(axis=0)
    c = cross_corr(W.ravel(), d)
    assert c['verdict'] == 'ar-order-low', f"{c['verdict']}"

def t_crosscorr_detects_needs_var():
    d, nE = 3, 6000
    Z = rng.standard_normal((nE, d))
    W = Z.copy()
    W[1:, 0] += 0.5*Z[:-1, 1]                        # axis0(k) depends on axis1(k-1)
    W = W/W.std(axis=0)
    c = cross_corr(W.ravel(), d)
    assert c['verdict'] == 'needs-var', f"{c['verdict']}"

# ============================================ 5. AR order selection
def t_order_selection_finds_true_order():
    dt = 5.0
    a_true = np.array([1, -1.6, 0.9, -0.25])
    f = np.linspace(1e-9, 0.5/dt, 4096)
    H = 1/(np.exp(-1j*2*np.pi*np.outer(f, np.arange(4))*dt) @ a_true)
    S = 0.7**2*dt*2*np.abs(H)**2
    info = select_ar_order(f, S, dt, acvTol=1e-3)
    assert info['verdict'] == 'good', info['verdict']
    assert info['p'] >= 3, f"selected p={info['p']}, too low for an AR(3)"
    assert info['acvRelErr'] < 1e-3

def t_order_selection_flags_hopeless():
    """A spectrum with a razor-sharp line: AR should struggle and say so."""
    dt = 1.0
    f = np.linspace(1e-6, 0.5, 8192)
    S = 1e-4*np.ones_like(f)
    S += 5.0*np.exp(-((f-0.1)/2e-4)**2)             # very narrow spectral line
    info = select_ar_order(f, S, dt, grid=[2, 5, 10, 20], acvTol=1e-3)
    assert info['verdict'] in ('marginal', 'bad'), info['verdict']

def t_order_selection_reports_curve():
    dt = 5.0
    a = np.array([1, -1.6, 0.9, -0.25])
    f = np.linspace(1e-9, 0.5/dt, 2048)
    H = 1/(np.exp(-1j*2*np.pi*np.outer(f, np.arange(4))*dt) @ a)
    S = dt*2*np.abs(H)**2
    info = select_ar_order(f, S, dt, grid=[2, 5, 10, 20, 40])
    assert len(info['curve']) >= 2
    errs = [c['acvRelErr'] for c in info['curve']]
    assert errs[0] > errs[-1], 'error should fall with order'

# ============================================ 6. VCE end-to-end
def t_vce_recovers_variance_components():
    nc, nl, nSys, nGrp, mGrp = 30, 5, 2, 5, 400
    TRUE = np.array([1.0, 2.0, 0.5, 3.0, 1.5])
    r2 = np.random.default_rng(99)
    xt = r2.standard_normal(nc)
    yt = [r2.standard_normal(nl) for _ in range(nSys)]
    data = []
    for s in range(nSys):
        for g in range(nGrp):
            A = r2.standard_normal((mGrp, nc))
            B = r2.standard_normal((mGrp, nl))
            l = A@xt + B@yt[s] + TRUE[g]*r2.standard_normal(mGrp)
            data.append(dict(A=A, B=B, l=l, sys=s, grp=g, m=mGrp))

    sig = np.ones(nGrp)
    for it in range(6):
        Ncc_red = np.zeros((nc, nc)); ncc = np.zeros(nc); info = []
        for s in range(nSys):
            ks = [k for k, d in enumerate(data) if d['sys'] == s]
            Ncc = np.zeros((nc, nc)); Nlc = np.zeros((nl, nc)); Nll = np.zeros((nl, nl))
            nc_ = np.zeros(nc); nl_ = np.zeros(nl)
            for k in ks:
                d = data[k]; w = 1/sig[d['grp']]
                A, B, l = d['A']*w, d['B']*w, d['l']*w
                Ncc += A.T@A; Nlc += B.T@A; Nll += B.T@B
                nc_ += A.T@l; nl_ += B.T@l
            Ncc_red += Ncc - Nlc.T@solve(Nll, Nlc)
            ncc += nc_ - Nlc.T@solve(Nll, nl_)
            info.append((Nll, Nlc, nl_, ks))
        x = solve(Ncc_red, ncc)

        rg = np.zeros(nGrp); vtv = np.zeros(nGrp)
        for (Nll, Nlc, nl_, ks) in info:
            Lc, Ll, G, SG, GSG = leverage_blocks(Ncc_red, Nll, Nlc)
            y = solve(Nll, nl_ - Nlc@x)
            for k in ks:
                d = data[k]; g = d['grp']; w = 1/sig[g]
                A, B, l = d['A']*w, d['B']*w, d['l']*w
                h = joint_leverage(A, B, Lc, Ll, G, SG, GSG)
                rg[g] += d['m'] - h.sum()
                v = l - A@x - B@y
                vtv[g] += v@v
        sig = sig*np.sqrt(vtv/np.maximum(rg, 1))

    M = sum(d['m'] for d in data); u = nc + nSys*nl
    assert abs(rg.sum() - (M-u)) < 1e-6, f'sum(r_g)={rg.sum():.2f} != {M-u}'
    ratio = sig/TRUE
    assert np.abs(ratio-1).max() < 0.12, f'ratio {np.round(ratio,3)}'

def t_outlier_detection():
    nc, nl, m = 20, 4, 800
    r2 = np.random.default_rng(5)
    A = r2.standard_normal((m, nc)); B = r2.standard_normal((m, nl))
    xt = r2.standard_normal(nc); yt = r2.standard_normal(nl)
    l = A@xt + B@yt + r2.standard_normal(m)
    bad = [17, 200, 555, 700]
    l[bad] += 10.0

    Ncc = A.T@A; Nlc = B.T@A; Nll = B.T@B
    Ncc_red = Ncc - Nlc.T@solve(Nll, Nlc)
    ncc = A.T@l - Nlc.T@solve(Nll, B.T@l)
    x = solve(Ncc_red, ncc)
    y = solve(Nll, B.T@l - Nlc@x)

    Lc, Ll, G, SG, GSG = leverage_blocks(Ncc_red, Nll, Nlc)
    h = joint_leverage(A, B, Lc, Ll, G, SG, GSG)
    v = l - A@x - B@y
    ri = 1 - h
    wnorm = np.abs(v)/np.sqrt(np.maximum(ri, 1e-9))
    flagged = set(np.where(wnorm > 4.0)[0])
    assert set(bad) <= flagged, f'missed {set(bad)-flagged}'
    assert len(flagged - set(bad)) <= 2, f'{len(flagged-set(bad))} false alarms'

for nm, fn in [
    ('leverage == brute-force hat matrix',      t_leverage_exact),
    ('sum(r_i) = m - u  identity',              t_redundancy_identity),
    ('prior normals correct',                   t_prior_normals_and_solve),
    ('prior redundancy / DoF absorbed',         t_prior_redundancy),
    ('prior rescues rank-deficient common',     t_prior_regularises_rank_deficient),
    ('prior is basis-agnostic (I/Kaula/L/sub)', t_prior_basis_agnostic),
    ('separable whiten == dense Cholesky',      t_separable_equals_dense),
    ('separable with full per-epoch 3x3 cov',   t_separable_full_epoch_cov),
    ('multi-axis costs the same as scalar',     t_separable_cost_is_scalar),
    ('crossCorr: white -> ok',                  t_crosscorr_detects_ok),
    ('crossCorr: detects cross-axis',           t_crosscorr_detects_cross_axis),
    ('crossCorr: detects AR order too low',     t_crosscorr_detects_ar_too_low),
    ('crossCorr: detects needs-VAR',            t_crosscorr_detects_needs_var),
    ('AR order: finds a sufficient order',      t_order_selection_finds_true_order),
    ('AR order: flags a hopeless spectrum',     t_order_selection_flags_hopeless),
    ('AR order: reports the curve',             t_order_selection_reports_curve),
    ('VCE recovers variance components',        t_vce_recovers_variance_components),
    ('outlier detection from the same h_ii',    t_outlier_detection),
]:
    T(nm, fn)

print('='*74)
npass = sum(1 for r in R if r[1] == 'PASS')
for nm, st, msg in R:
    print(f"  [{'ok  ' if st=='PASS' else 'FAIL'}] {nm:<42s} {msg}")
print('='*74)
print(f'  {npass}/{len(R)} passed')
