"""Integrated v3: prior in solve(), block-diagonal Schur, epoch-aligned batching,
VCE + outliers through the full accumulator path. Mirrors the MATLAB semantics."""
import numpy as np
from numpy.linalg import solve, inv, norm, lstsq
from scipy.linalg import cholesky, solve_triangular, toeplitz, block_diag
from scipy.signal import lfilter
import sys
sys.path.insert(0, '/home/claude')
from v3core import (joint_leverage, leverage_blocks, CommonPrior,
                    separable_whiten, cross_corr, select_ar_order,
                    psd2acv, levinson, arma_acv)

R = []
def T(n, f):
    try:
        f(); R.append((n, 'PASS', ''))
    except AssertionError as e:
        R.append((n, 'FAIL', str(e)[:60]))
    except Exception as e:
        R.append((n, 'ERROR', f'{type(e).__name__}: {e}'[:60]))

rng = np.random.default_rng(21)


# ==================================================== accumulator (v3 semantics)
class Acc:
    def __init__(self, nc, prior=None, blkSize=None):
        self.nc = nc; self.prior = prior; self.blkSize = blkSize
        self.Ncc = np.zeros((nc, nc)); self.ncc = np.zeros(nc)
        self.ltl = 0.0; self.m = 0; self.store = {}

    def addSystem(self, sid, nl, arcs):
        """arcs: list of (Aw, Bw, lw) already whitened."""
        Ncc = np.zeros((self.nc, self.nc)); Nlc = np.zeros((nl, self.nc))
        Nll = np.zeros((nl, nl)); nc_ = np.zeros(self.nc); nl_ = np.zeros(nl)
        ltl = 0.0; m = 0
        for A, B, l in arcs:
            Ncc += A.T@A; Nlc += B.T@A; Nll += B.T@B
            nc_ += A.T@l; nl_ += B.T@l; ltl += l@l; m += len(l)

        RHS = np.column_stack([Nlc, nl_])
        # block-diagonal Schur when the locals are arc-wise
        if self.blkSize and nl % self.blkSize == 0 and nl > self.blkSize:
            nB = nl // self.blkSize
            Z = np.zeros((nl, self.nc + 1))
            for b in range(nB):
                ib = slice(b*self.blkSize, (b+1)*self.blkSize)
                Z[ib] = solve(Nll[ib, ib], RHS[ib])
        else:
            Z = solve(Nll, RHS)

        self.Ncc += Ncc - Nlc.T@Z[:, :self.nc]
        self.ncc += nc_ - Nlc.T@Z[:, -1]
        self.ltl += ltl - nl_@Z[:, -1]
        self.m += m - nl
        self.store[sid] = dict(Nll=Nll, Nlc=Nlc, nl=nl_, nlp=nl, arcs=arcs, m=m)

    def reducedNormalMatrix(self):
        N = self.Ncc.copy()
        if self.prior is not None:
            N = N + self.prior.normals()[0]
        return (N + N.T)/2

    def solve(self):
        if self.prior is None:
            N = (self.Ncc + self.Ncc.T)/2
            x = solve(N, self.ncc)
            vtv = self.ltl - self.ncc@x
            dof = self.m - self.nc
            absorbed, rK = 0.0, np.nan
        else:
            Nk, nk = self.prior.normals()
            N = (self.Ncc + Nk + (self.Ncc + Nk).T)/2
            n = self.ncc + nk
            Lc = cholesky(N, lower=True)
            x = solve_triangular(Lc.T, solve_triangular(Lc, n, lower=True))
            x0w = self.prior.x0 @ self.prior.x0 / self.prior.sigma**2
            vtv = max(self.ltl + x0w - n@x, 0.0)
            rK, absorbed = self.prior.redundancy(Lc)
            dof = self.m + self.prior.nk - self.nc
        s0sq = vtv/max(dof, 1)
        return dict(x=x, s0sq=s0sq, sigma0=np.sqrt(s0sq), dof=dof, vtv=vtv,
                    priorAbsorbed=absorbed, priorRedundancy=rK)


# ==================================================== tests
def t_prior_in_solve_matches_augmented():
    """solve() with a prior must equal an explicitly augmented LSQ."""
    nc, nl, m = 10, 4, 300
    A = rng.standard_normal((m, nc)); B = rng.standard_normal((m, nl))
    l = rng.standard_normal(m)
    K = np.diag(1/np.arange(1, nc+1)**1.5)
    x0 = 0.1*rng.standard_normal(nc)
    sig = 0.7
    pri = CommonPrior(K, x0, sig)

    acc = Acc(nc, prior=pri); acc.addSystem(1, nl, [(A, B, l)])
    sol = acc.solve()

    # explicit augmentation of the JOINT system, then extract x_c
    Aa = np.zeros((m + K.shape[0], nc + nl))
    Aa[:m, :nc] = A; Aa[:m, nc:] = B
    Aa[m:, :nc] = K/sig
    la = np.concatenate([l, x0/sig])
    xf = lstsq(Aa, la, rcond=None)[0][:nc]

    e = norm(sol['x'] - xf)/norm(xf)
    assert e < 1e-9, f'{e:.2e}'


def t_prior_absorbed_bounds():
    nc, nl, m = 12, 4, 200
    A = rng.standard_normal((m, nc)); B = rng.standard_normal((m, nl))
    l = rng.standard_normal(m)
    for sig, lo, hi in [(0.01, 0.9*nc, nc+1e-6), (100.0, -1e-6, 0.1*nc)]:
        pri = CommonPrior(np.eye(nc), None, sig)
        acc = Acc(nc, prior=pri); acc.addSystem(1, nl, [(A, B, l)])
        s = acc.solve()
        a = s['priorAbsorbed']
        assert lo <= a <= hi, f'sigma={sig}: absorbed={a:.3f} not in [{lo:.2f},{hi:.2f}]'


def t_prior_rescues_singular_common():
    nc, nl, m = 8, 3, 60
    A = rng.standard_normal((m, nc)); A[:, -1] = A[:, -2]      # singular
    B = rng.standard_normal((m, nl)); l = rng.standard_normal(m)
    acc0 = Acc(nc); acc0.addSystem(1, nl, [(A, B, l)])
    try:
        cholesky(acc0.reducedNormalMatrix(), lower=True)
        assert False, 'should be singular without a prior'
    except Exception:
        pass
    acc = Acc(nc, prior=CommonPrior(np.eye(nc), None, 1.0))
    acc.addSystem(1, nl, [(A, B, l)])
    s = acc.solve()
    assert np.all(np.isfinite(s['x']))


def t_block_schur_equals_dense():
    """Arc-wise locals: block-diagonal Schur must equal the dense Schur."""
    nc, nArc, nlA = 8, 20, 6
    nl = nArc*nlA
    arcs = []
    for a in range(nArc):
        mA = 60
        A = rng.standard_normal((mA, nc))
        B = np.zeros((mA, nl))
        B[:, a*nlA:(a+1)*nlA] = rng.standard_normal((mA, nlA))
        arcs.append((A, B, rng.standard_normal(mA)))

    a1 = Acc(nc);                a1.addSystem(1, nl, arcs)
    a2 = Acc(nc, blkSize=nlA);   a2.addSystem(1, nl, arcs)

    e = norm(a2.Ncc - a1.Ncc)/norm(a1.Ncc)
    en = norm(a2.ncc - a1.ncc)/norm(a1.ncc)
    assert e < 1e-10 and en < 1e-10, f'Ncc {e:.2e}  ncc {en:.2e}'
    s1, s2 = a1.solve(), a2.solve()
    assert norm(s2['x']-s1['x'])/norm(s1['x']) < 1e-9


def t_block_schur_is_faster():
    import time
    nc, nArc, nlA = 12, 60, 6
    nl = nArc*nlA
    arcs = []
    for a in range(nArc):
        A = rng.standard_normal((80, nc))
        B = np.zeros((80, nl)); B[:, a*nlA:(a+1)*nlA] = rng.standard_normal((80, nlA))
        arcs.append((A, B, rng.standard_normal(80)))

    def run(bs):
        t0 = time.perf_counter()
        a = Acc(nc, blkSize=bs); a.addSystem(1, nl, arcs)
        return time.perf_counter()-t0
    td = min(run(None) for _ in range(3))
    tb = min(run(nlA) for _ in range(3))
    assert tb < td, f'block {tb*1e3:.1f}ms not faster than dense {td*1e3:.1f}ms'


def t_vce_end_to_end():
    """Full VCE through the accumulator, with leverage-based redundancy."""
    nc, nl, nSys, nGrp, mG = 25, 4, 2, 4, 350
    TRUE = np.array([1.0, 2.5, 0.6, 1.5])
    r2 = np.random.default_rng(77)
    xt = r2.standard_normal(nc); yt = [r2.standard_normal(nl) for _ in range(nSys)]
    raw = []
    for s in range(nSys):
        for g in range(nGrp):
            A = r2.standard_normal((mG, nc)); B = r2.standard_normal((mG, nl))
            l = A@xt + B@yt[s] + TRUE[g]*r2.standard_normal(mG)
            raw.append((s, g, A, B, l))

    sig = np.ones(nGrp)
    for it in range(6):
        acc = Acc(nc)
        for s in range(nSys):
            arcs = [(A/sig[g], B/sig[g], l/sig[g])
                    for (ss, g, A, B, l) in raw if ss == s]
            acc.addSystem(s, nl, arcs)
        sol = acc.solve()
        Nred = acc.reducedNormalMatrix()

        rg = np.zeros(nGrp); vtv = np.zeros(nGrp)
        for s in range(nSys):
            st = acc.store[s]
            blk = leverage_blocks(Nred, st['Nll'], st['Nlc'])
            y = solve(st['Nll'], st['nl'] - st['Nlc']@sol['x'])
            for (ss, g, A, B, l) in raw:
                if ss != s: continue
                Aw, Bw, lw = A/sig[g], B/sig[g], l/sig[g]
                h = joint_leverage(Aw, Bw, *blk)
                rg[g] += len(lw) - h.sum()
                v = lw - Aw@sol['x'] - Bw@y
                vtv[g] += v@v
        sig = sig*np.sqrt(vtv/np.maximum(rg, 1))

    M = sum(len(l) for (_, _, _, _, l) in raw); u = nc + nSys*nl
    assert abs(rg.sum() - (M-u)) < 1e-6, f'sum(r)={rg.sum():.3f} != {M-u}'
    ratio = sig/TRUE
    assert np.abs(ratio-1).max() < 0.15, f'ratio={np.round(ratio,3)}'


def t_epoch_aligned_batching():
    """Multi-axis: batches must not split an epoch, and the result must be
    identical to whitening the whole arc."""
    d, nE = 3, 300
    m = nE*d
    a = np.array([1, -1.1, 0.35]); s = 0.9
    r = s**2*arma_acv(a, [1.], nE)
    sig = 0.5 + rng.random((nE, d))
    A = rng.standard_normal((m, 5)); B = rng.standard_normal((m, 3))
    l = rng.standard_normal(m)

    Aw, Bw, lw = separable_whiten(A, B, l, a, s, r, sig, d, True)
    Nref = np.hstack([Aw, Bw]).T @ np.hstack([Aw, Bw])

    # batched, overlap = p epochs = p*d rows, aligned to epochs
    p = len(a)-1
    bsE = 64                                  # epochs per batch
    bs, nov = bsE*d, p*d
    N = np.zeros((8, 8))
    for i0 in range(0, m, bs):
        i1 = min(i0+bs, m); j0 = max(0, i0-nov)
        idx = slice(j0, i1); keep = slice(i0-j0, i1-j0)
        e0, e1 = j0//d, i1//d
        Ab, Bb, lb = separable_whiten(A[idx], B[idx], l[idx], a, s, r,
                                      sig[e0:e1], d, j0 == 0)
        X = np.hstack([Ab, Bb])[keep]
        N += X.T@X
    e = norm(N - Nref)/norm(Nref)
    assert e < 1e-8, f'batched normals differ by {e:.2e}'


def t_outliers_via_leverage():
    nc, nl, m = 18, 4, 700
    r2 = np.random.default_rng(9)
    A = r2.standard_normal((m, nc)); B = r2.standard_normal((m, nl))
    l = A@r2.standard_normal(nc) + B@r2.standard_normal(nl) + r2.standard_normal(m)
    bad = [23, 111, 400, 650]
    l[bad] += 11.0

    acc = Acc(nc); acc.addSystem(1, nl, [(A, B, l)])
    sol = acc.solve()
    st = acc.store[1]
    blk = leverage_blocks(acc.reducedNormalMatrix(), st['Nll'], st['Nlc'])
    y = solve(st['Nll'], st['nl'] - st['Nlc']@sol['x'])
    h = joint_leverage(A, B, *blk)
    v = l - A@sol['x'] - B@y
    w = np.abs(v)/np.sqrt(np.maximum(1-h, 1e-9))
    flag = set(np.where(w > 4.0)[0])
    assert set(bad) <= flag, f'missed {set(bad)-flag}'
    assert len(flag - set(bad)) <= 2, f'{len(flag-set(bad))} false alarms'


def t_v2_unbroken_no_prior():
    """The v3 changes must not alter the v2 result when no prior is set."""
    nc, nl, m = 10, 4, 400
    A = rng.standard_normal((m, nc)); B = rng.standard_normal((m, nl))
    l = rng.standard_normal(m)
    acc = Acc(nc); acc.addSystem(1, nl, [(A, B, l)])
    sol = acc.solve()

    # monolithic reference
    X = np.hstack([A, B])
    xf = lstsq(X, l, rcond=None)[0][:nc]
    assert norm(sol['x']-xf)/norm(xf) < 1e-9
    assert sol['priorAbsorbed'] == 0


for nm, fn in [
    ('prior in solve == augmented LSQ',    t_prior_in_solve_matches_augmented),
    ('prior absorbed DoF in [0, nk]',      t_prior_absorbed_bounds),
    ('prior rescues singular common',      t_prior_rescues_singular_common),
    ('block Schur == dense Schur',         t_block_schur_equals_dense),
    ('block Schur is faster',              t_block_schur_is_faster),
    ('VCE end-to-end via accumulator',     t_vce_end_to_end),
    ('epoch-aligned batching is exact',    t_epoch_aligned_batching),
    ('outliers via joint leverage',        t_outliers_via_leverage),
    ('v2 result unchanged (no prior)',     t_v2_unbroken_no_prior),
]:
    T(nm, fn)

print('='*72)
npass = sum(1 for r in R if r[1] == 'PASS')
for nm, st, msg in R:
    print(f"  [{'ok  ' if st=='PASS' else 'FAIL'}] {nm:<38s} {msg}")
print('='*72)
print(f'  {npass}/{len(R)} passed')
