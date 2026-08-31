#!/usr/bin/env python3
"""v5 validation: day-distributed accumulation WITHIN one system (ROADMAP 3.2).

The normal-form day partition was already proven (v4verify test 1). What is NEW
here is the QR path: per-day COMPRESSED triangular factors ((dim+1) rows, the
last row carrying the residual mass) merged by stack + re-QR + truncate must
reproduce the fused serial QR accumulation exactly -- including the l'l
bookkeeping that sigma0 depends on. If that holds, the day-distributed path
keeps the QR conditioning advantage end to end and needs no normal-form detour.

    python3 v5verify.py
"""
import numpy as np
from numpy.linalg import norm, solve, lstsq
from scipy.signal import lfilter

R = []
def T(n, f):
    try:
        f(); R.append((n, 'PASS', ''))
    except AssertionError as e:
        R.append((n, 'FAIL', str(e)[:70]))
    except Exception as e:
        R.append((n, 'ERROR', f'{type(e).__name__}: {e}'[:70]))

rng = np.random.default_rng(51)


# ------------------------------------------------------------ problem builder
def ar_noise(a, s, m, r):
    return lfilter([1.0], a, s * r.standard_normal(m + 200))[200:]

def whiten_fir(a, s, X):
    if X.ndim == 1:
        return lfilter(a, [1.0], X) / s
    return np.column_stack([lfilter(a, [1.0], X[:, j])
                            for j in range(X.shape[1])]) / s

def build_days(nc=6, nl=4, nDay=5, arcsPerDay=2, mArc=200, seed=11, nSys=2):
    """Whitened per-day arc lists per system + the truth."""
    r = np.random.default_rng(seed)
    a = np.array([1.0, -1.2, 0.5]); s = 0.8
    xt = r.standard_normal(nc)
    sys = {}
    for sid in range(nSys):
        yt = r.standard_normal(nl)
        days = []
        for _ in range(nDay):
            arcs = []
            for _ in range(arcsPerDay):
                A = r.standard_normal((mArc, nc))
                B = r.standard_normal((mArc, nl))
                l = A @ xt + B @ yt + ar_noise(a, s, mArc, r)
                arcs.append((whiten_fir(a, s, A), whiten_fir(a, s, B),
                             whiten_fir(a, s, l)))
            days.append(arcs)
        sys[sid] = days
    return xt, sys, nc, nl


# ------------------------------------------------- compressed-QR reference ops
def qr_stack(Rc, rows, dim):
    """The QRSystem.add / mergeRaw operation: stack, re-triangularise,
    truncate to dim+1 rows (row dim+1 = accumulated residual norm)."""
    M = np.vstack([Rc, rows]) if Rc.size else np.atleast_2d(rows)
    Rn = np.linalg.qr(M, mode='r')
    Rn = np.triu(Rn)
    return Rn[:min(Rn.shape[0], dim + 1), :]

def qr_accum(arclist, dim):
    """Fused serial QRSystem over a list of whitened arcs."""
    Rc = np.zeros((0, dim + 1)); m = 0
    for A, B, l in arclist:
        Rc = qr_stack(Rc, np.column_stack([A, B, l]), dim)
        m += len(l)
    return Rc, m

def qr_blocks(Rc, dim, nc):
    """QRSystem.blocks: N = R1'R1, n = R1'r, ltl incl. the residual row."""
    R1 = np.zeros((dim, dim)); rr = np.zeros(dim)
    k = min(Rc.shape[0], dim)
    R1[:k, :] = Rc[:k, :dim]; rr[:k] = Rc[:k, dim]
    N = R1.T @ R1; n = R1.T @ rr
    ltl = rr @ rr + (Rc[dim, dim] ** 2 if Rc.shape[0] > dim else 0.0)
    ic = np.arange(nc); il = np.arange(nc, dim)
    return (N[np.ix_(il, il)], N[np.ix_(il, ic)], N[np.ix_(ic, ic)],
            n[il], n[ic], ltl)

def schur(Nll, Nlc, Ncc, nlr, ncr, ltl, m, nl):
    Z = solve(Nll, np.column_stack([Nlc, nlr]))
    return (Ncc - Nlc.T @ Z[:, :-1], ncr - Nlc.T @ Z[:, -1],
            ltl - nlr @ Z[:, -1], m - nl)


# ============ 1. compression is lossless (the property the design rests on)
def t_compressed_factor_is_lossless():
    """A day with m >> dim+1 rows, truncated to dim+1: N, n, l'l preserved
    EXACTLY. Column-norm preservation under orthogonal transformation -- but
    assert it numerically, including the residual row."""
    _, S, nc, nl = build_days(nDay=1, arcsPerDay=3, mArc=300)
    dim = nc + nl
    arcs = S[0][0]
    M = np.vstack([np.column_stack(a) for a in arcs])
    Rc, _ = qr_accum(arcs, dim)
    assert Rc.shape[0] == dim + 1

    Nll, Nlc, Ncc, nlr, ncr, ltl = qr_blocks(Rc, dim, nc)
    Nref = M[:, :dim].T @ M[:, :dim]
    nref = M[:, :dim].T @ M[:, dim]
    e1 = norm(np.block([[Ncc, Nlc.T], [Nlc, Nll]]) - Nref) / norm(Nref)
    e2 = norm(np.concatenate([ncr, nlr]) - nref) / norm(nref)
    e3 = abs(ltl - M[:, dim] @ M[:, dim]) / (M[:, dim] @ M[:, dim])
    assert e1 < 1e-13 and e2 < 1e-13 and e3 < 1e-13, f'{e1:.2e} {e2:.2e} {e3:.2e}'


# ================== 2. day-partial stacking merge == fused serial QR
def t_qr_day_merge_equals_fused():
    """Per-day compressed factors merged in day order == the fused loop, in
    every derived quantity (N, n, ltl). This is exactly QRSystem.mergeRaw."""
    _, S, nc, nl = build_days()
    dim = nc + nl
    days = S[0]

    Rf, mf = qr_accum([a for d in days for a in d], dim)

    Rm = np.zeros((0, dim + 1)); mm = 0
    for d in days:
        Rd, md = qr_accum(d, dim)
        Rm = qr_stack(Rm, Rd, dim)          # mergeRaw: partial R as pseudo-obs
        mm += md

    bf = qr_blocks(Rf, dim, nc); bm = qr_blocks(Rm, dim, nc)
    assert mm == mf
    for k, (a, b) in enumerate(zip(bf, bm)):
        e = norm(np.atleast_1d(a - b)) / max(norm(np.atleast_1d(a)), 1e-300)
        assert e < 5e-13, f'block {k}: rel {e:.2e}'

def t_qr_day_merge_order_irrelevant():
    """Scrambled day order: R itself differs (signs / rotations), but every
    derived quantity is identical."""
    _, S, nc, nl = build_days()
    dim = nc + nl
    days = S[0]

    def merged(order):
        Rm = np.zeros((0, dim + 1))
        for i in order:
            Rd, _ = qr_accum(days[i], dim)
            Rm = qr_stack(Rm, Rd, dim)
        return qr_blocks(Rm, dim, nc)

    ba = merged(range(len(days)))
    bb = merged(rng.permutation(len(days)))
    for k, (a, b) in enumerate(zip(ba, bb)):
        e = norm(np.atleast_1d(a - b)) / max(norm(np.atleast_1d(a)), 1e-300)
        assert e < 5e-13, f'block {k}: rel {e:.2e}'


# ============ 3. the sigma0 guard: residual-row mass survives the merge
def t_ltl_bookkeeping_survives_merge():
    """The known silent failure mode: dropping the residual row makes sigma0
    collapse while x stays plausible. Assert the merged compressed factor
    reproduces sum(l'l) exactly, and that WITHOUT the residual rows it does
    not (the guard actually guards)."""
    _, S, nc, nl = build_days(nDay=4)
    dim = nc + nl
    days = S[0]
    ltl_true = sum(l @ l for d in days for _, _, l in d)

    Rm = np.zeros((0, dim + 1))
    Rbad = np.zeros((0, dim + 1))
    for d in days:
        Rd, _ = qr_accum(d, dim)
        Rm = qr_stack(Rm, Rd, dim)
        Rbad = qr_stack(Rbad, Rd[:dim, :], dim)    # residual row dropped

    ltl = qr_blocks(Rm, dim, nc)[5]
    ltl_bad = qr_blocks(Rbad, dim, nc)[5]
    assert abs(ltl - ltl_true) / ltl_true < 1e-13, f'{abs(ltl-ltl_true)/ltl_true:.2e}'
    assert (ltl_true - ltl_bad) / ltl_true > 1e-3, 'guard test has no teeth'


# ====== 4. full pipeline, both backends, vs the monolithic dense reference
def t_full_pipeline_day_distributed_vs_monolithic():
    """Day-partitioned partials -> merge -> Schur -> cross-system merge -> x,
    sigma0. Both backends, against one joint dense lstsq on everything."""
    xt, S, nc, nl = build_days(nDay=4, arcsPerDay=2, mArc=150, nSys=2)
    dim = nc + nl
    nSys = len(S)

    # ---- monolithic reference -------------------------------------------
    rows = []; rhs = []
    for sid in range(nSys):
        for d in S[sid]:
            for A, B, l in d:
                blk = np.zeros((len(l), nc + nSys * nl))
                blk[:, :nc] = A
                blk[:, nc + sid * nl: nc + (sid + 1) * nl] = B
                rows.append(blk); rhs.append(l)
    Aj = np.vstack(rows); lj = np.concatenate(rhs)
    xj = lstsq(Aj, lj, rcond=None)[0]
    vj = lj - Aj @ xj
    s0j = np.sqrt(vj @ vj / (len(lj) - Aj.shape[1]))

    # ---- toolbox pipelines ----------------------------------------------
    def via(backend):
        Nt = np.zeros((nc, nc)); nt = np.zeros(nc); lt = 0.0; dof = 0
        for sid in range(nSys):
            if backend == 'qr':
                Rm = np.zeros((0, dim + 1)); m = 0
                for d in S[sid]:
                    Rd, md = qr_accum(d, dim)
                    Rm = qr_stack(Rm, Rd, dim); m += md
                Nll, Nlc, Ncc, nlr, ncr, ltl = qr_blocks(Rm, dim, nc)
            else:
                Nll = np.zeros((nl, nl)); Nlc = np.zeros((nl, nc))
                Ncc = np.zeros((nc, nc)); nlr = np.zeros(nl)
                ncr = np.zeros(nc); ltl = 0.0; m = 0
                for d in S[sid]:                    # per-day partials, summed
                    for A, B, l in d:
                        Nll += B.T @ B; Nlc += B.T @ A; Ncc += A.T @ A
                        nlr += B.T @ l; ncr += A.T @ l
                        ltl += l @ l; m += len(l)
            Nr, nr, lr, dr = schur(Nll, Nlc, Ncc, nlr, ncr, ltl, m, nl)
            Nt += Nr; nt += nr; lt += lr; dof += dr
        x = solve((Nt + Nt.T) / 2, nt)
        s0 = np.sqrt(max(lt - nt @ x, 0) / (dof - nc))
        return x, s0

    for be in ('qr', 'normal'):
        x, s0 = via(be)
        ex = norm(x - xj[:nc]) / norm(xj[:nc])
        es = abs(s0 - s0j) / s0j
        assert ex < 1e-10 and es < 1e-10, f'{be}: x {ex:.2e}  s0 {es:.2e}'


# ==== 5. partial-vs-running summation gap is the only serial/parallel diff
def t_day_partial_summation_gap_is_epsilon():
    """Normal path: per-day partial sums then reduced vs the fused running
    sum. The difference is pure summation-order rounding -- quantify it so the
    MATLAB regression tolerance (1e-12) is justified, not guessed."""
    _, S, nc, nl = build_days(nDay=6, arcsPerDay=3, mArc=250)
    days = S[0]

    Nf = np.zeros((nc + nl, nc + nl))
    for d in days:
        for A, B, l in d:
            J = np.column_stack([A, B]); Nf += J.T @ J

    parts = []
    for d in days:
        Nd = np.zeros_like(Nf)
        for A, B, l in d:
            J = np.column_stack([A, B]); Nd += J.T @ J
        parts.append(Nd)
    Np = sum(parts)

    e = norm(Np - Nf) / norm(Nf)
    assert e < 1e-14, f'rel {e:.2e}'


T('1. compressed (dim+1)-row factor is lossless (N, n, ltl exact)',
  t_compressed_factor_is_lossless)
T('2. QR day-partial stack-merge == fused serial accumulation',
  t_qr_day_merge_equals_fused)
T('3. QR day-merge is order-independent in all derived quantities',
  t_qr_day_merge_order_irrelevant)
T('4. residual-row mass survives the merge (sigma0 guard has teeth)',
  t_ltl_bookkeeping_survives_merge)
T('5. full day-distributed pipeline == monolithic, both backends',
  t_full_pipeline_day_distributed_vs_monolithic)
T('6. day-partial summation gap ~ eps (justifies the 1e-12 tolerance)',
  t_day_partial_summation_gap_is_epsilon)

n = sum(1 for _, s, _ in R if s == 'PASS')
for name, s, msg in R:
    print(f'  [{s:5s}] {name}' + (f'  -- {msg}' if msg else ''))
print(f'  {n}/{len(R)} passed')
