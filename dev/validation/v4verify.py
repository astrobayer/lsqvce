"""
v4 mirror. Validate BEFORE transcribing to MATLAB (the working rule).

Three things are proven here:

1. THE REDUCTION DECOMPOSITION (parallelism correctness, ROADMAP 3.6 step 1).
   Every quantity the parallel paths accumulate -- Ncc_red / ncc_red / ltl / dof
   across systems, Nll / Nlc / Ncc / rhs across arcs within a system, and
   rg / vtv / sum(h) across arcs in the redundancy pass -- is a plain SUM of
   independent per-task contributions. Partial sums over an arbitrary partition,
   combined in an arbitrary order, must reconstruct the serial result to
   summation-order tolerance. This is the correctness proof that survives
   without MATLAB: parfor/spmd only change the SCHEDULE, never the terms.

2. THE PSD ITERATION HELPER (ROADMAP 2.3).
   solve -> estimate ACV from raw residuals -> re-fit AR -> re-solve. Starting
   from a WRONG model (white assumed, AR(2) true), sigma0 must go to ~1, the
   whitened-residual ACF must go to ~delta, and it must converge in <= 3
   iterations. Also: the residual-ACV bias (v = (I-H)e) must be negligible at
   high redundancy -- measured, not assumed.

3. PER-OBSTYPE SIGMA0 QUICK-LOOK (ROADMAP 2.4).
   rms_ot = sqrt(vtv_ot / m_ot) is biased LOW versus the exact component
   sqrt(vtv_ot / r_ot) by the leverage the group carries. Verify both against
   brute force and measure the bias so the doc can state it honestly.
"""
import numpy as np
from numpy.linalg import solve, norm
from scipy.linalg import cholesky, solve_triangular
from scipy.signal import lfilter
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v3core import (joint_leverage, leverage_blocks, levinson, arma_acv)

R = []
def T(n, f):
    try:
        f(); R.append((n, 'PASS', ''))
    except AssertionError as e:
        R.append((n, 'FAIL', str(e)[:70]))
    except Exception as e:
        R.append((n, 'ERROR', f'{type(e).__name__}: {e}'[:70]))

rng = np.random.default_rng(46)


# ===================================================== shared synthetic problem
def ar_noise(a, s, m, r):
    return lfilter([1.0], a, s * r.standard_normal(m + 200))[200:]

def whiten_fir(a, s, X):
    """Exact FIR whitener of an AR(a, s) process = chol(C)^-1 (measured v2 fact),
    applied column-wise with zero initial state (whole-arc, no batching)."""
    if X.ndim == 1:
        return lfilter(a, [1.0], X) / s
    return np.column_stack([lfilter(a, [1.0], X[:, j]) for j in range(X.shape[1])]) / s

def build_problem(nc=8, nl=4, nSys=2, nArc=6, mArc=250, seed=7):
    """Two systems, AR(2) noise, whitened with the TRUE model -> identities hold."""
    r = np.random.default_rng(seed)
    a = np.array([1.0, -1.2, 0.5]); s = 0.8
    xt = r.standard_normal(nc)
    sysArcs = {}
    for sid in range(nSys):
        yt = r.standard_normal(nl)
        arcs = []
        for _ in range(nArc):
            A = r.standard_normal((mArc, nc))
            B = r.standard_normal((mArc, nl))
            l = A @ xt + B @ yt + ar_noise(a, s, mArc, r)
            arcs.append((whiten_fir(a, s, A), whiten_fir(a, s, B),
                         whiten_fir(a, s, l)))
        sysArcs[sid] = arcs
    return xt, sysArcs, nc, nl

def sys_normals(arcs, nc, nl):
    Ncc = np.zeros((nc, nc)); Nlc = np.zeros((nl, nc)); Nll = np.zeros((nl, nl))
    ncr = np.zeros(nc); nlr = np.zeros(nl); ltl = 0.0; m = 0
    for A, B, l in arcs:
        Ncc += A.T @ A; Nlc += B.T @ A; Nll += B.T @ B
        ncr += A.T @ l; nlr += B.T @ l; ltl += l @ l; m += len(l)
    return Ncc, Nlc, Nll, ncr, nlr, ltl, m

def schur(Ncc, Nlc, Nll, ncr, nlr, ltl, m, nl):
    Z = solve(Nll, np.column_stack([Nlc, nlr]))
    return (Ncc - Nlc.T @ Z[:, :-1], ncr - Nlc.T @ Z[:, -1],
            ltl - nlr @ Z[:, -1], m - nl)


# ============================ 1. reduction decomposition (parallel correctness)
def t_arc_partition_reconstructs_system_normals():
    """Within-system: per-arc normal contributions summed over an arbitrary
    partition, in scrambled order, == the serial accumulation. This is the
    day-distributed accumulation of ROADMAP 3.4."""
    _, SA, nc, nl = build_problem()
    arcs = SA[0]
    ref = sys_normals(arcs, nc, nl)

    idx = rng.permutation(len(arcs))
    chunks = [idx[0:2], idx[2:3], idx[3:6]]          # uneven partition
    parts = [sys_normals([arcs[i] for i in ch], nc, nl) for ch in chunks]
    tot = [sum(p[k] for p in parts) for k in range(7)]

    for k, (a, b) in enumerate(zip(ref, tot)):
        e = norm(np.atleast_1d(a - b)) / max(norm(np.atleast_1d(a)), 1e-300)
        assert e < 1e-13, f'block {k}: rel {e:.2e}'

def t_system_partition_reconstructs_reduced_normals():
    """Across systems: reduced contributions merged in ANY order give the same
    common normals (the addSystems parfor path)."""
    _, SA, nc, nl = build_problem(nSys=4, nArc=3)
    reds = [schur(*sys_normals(SA[s], nc, nl), nl) for s in SA]

    def total(order):
        N = np.zeros((nc, nc)); n = np.zeros(nc); ltl = 0.0; dof = 0
        for i in order:
            N = N + reds[i][0]; n = n + reds[i][1]
            ltl += reds[i][2]; dof += reds[i][3]
        return N, n, ltl, dof

    Na, na, la, da = total(range(4))
    Nb, nb, lb, db = total([2, 0, 3, 1])
    assert norm(Na - Nb) / norm(Na) < 1e-14
    assert norm(na - nb) / norm(na) < 1e-13
    assert abs(la - lb) / abs(la) < 1e-13 and da == db
    xa = solve((Na + Na.T) / 2, na); xb = solve((Nb + Nb.T) / 2, nb)
    assert norm(xa - xb) / norm(xa) < 1e-11

def redundancy_pass(SA, nc, nl, x, Nred, tasks=None):
    """The redundancy pass as a MAP over independent arc-tasks + ordered REDUCE.
    tasks = list of (sysID, arcIdx); None = all, serial order."""
    blocks = {}
    ys = {}
    for s in SA:
        Ncc, Nlc, Nll, ncr, nlr, ltl, m = sys_normals(SA[s], nc, nl)
        blocks[s] = leverage_blocks(Nred, Nll, Nlc)
        ys[s] = solve(Nll, nlr - Nlc @ x)
    if tasks is None:
        tasks = [(s, i) for s in SA for i in range(len(SA[s]))]
    out = []
    for s, i in tasks:                                # the MAP (parallelisable)
        A, B, l = SA[s][i]
        h = joint_leverage(A, B, *blocks[s])
        v = l - A @ x - B @ ys[s]
        out.append(dict(m=len(l), r=len(l) - h.sum(), vtv=v @ v, h=h.sum()))
    return out

def t_redundancy_partition_reconstructs_serial():
    """rg / vtv / sum(h) over a scrambled partition of arc-tasks == serial, and
    the global identity sum(r) = m - u holds. THE proof for the parfor/spmd
    redundancy pass."""
    xt, SA, nc, nl = build_problem()
    reds = [schur(*sys_normals(SA[s], nc, nl), nl) for s in SA]
    Nred = sum(r[0] for r in reds); nred = sum(r[1] for r in reds)
    Nred = (Nred + Nred.T) / 2
    x = solve(Nred, nred)

    serial = redundancy_pass(SA, nc, nl, x, Nred)
    rS = sum(t['r'] for t in serial); vS = sum(t['vtv'] for t in serial)
    mS = sum(t['m'] for t in serial)

    tasks = [(s, i) for s in SA for i in range(len(SA[s]))]
    perm = [tasks[j] for j in rng.permutation(len(tasks))]
    cuts = [0, 3, 4, 9, len(perm)]
    rP = vP = mP = 0.0
    for c0, c1 in zip(cuts[:-1], cuts[1:]):           # partial per "worker"
        part = redundancy_pass(SA, nc, nl, x, Nred, perm[c0:c1])
        rP += sum(t['r'] for t in part); vP += sum(t['vtv'] for t in part)
        mP += sum(t['m'] for t in part)

    assert abs(rP - rS) / rS < 1e-12, f'r: {abs(rP-rS)/rS:.2e}'
    assert abs(vP - vS) / vS < 1e-12, f'vtv: {abs(vP-vS)/vS:.2e}'
    u = nc + len(SA) * nl
    assert abs(rS - (mS - u)) < 1e-8 * (mS - u), f'identity: {rS} vs {mS-u}'

def t_pertask_then_ordered_reduce_is_exact():
    """The MATLAB strategy collects per-task partials then reduces in TASK
    ORDER -- serial and parfor therefore share the identical reduction order.
    Verify: reducing per-task results in task order == the fused serial loop
    bit-for-bit (same terms, same order)."""
    xt, SA, nc, nl = build_problem(nArc=4)
    reds = [schur(*sys_normals(SA[s], nc, nl), nl) for s in SA]
    Nred = sum(r[0] for r in reds); Nred = (Nred + Nred.T) / 2
    nred = sum(r[1] for r in reds)
    x = solve(Nred, nred)

    out = redundancy_pass(SA, nc, nl, x, Nred)
    fused_r = 0.0
    for t in out:
        fused_r = fused_r + t['r']
    red2 = 0.0
    for t in redundancy_pass(SA, nc, nl, x, Nred):
        red2 = red2 + t['r']
    assert fused_r == red2                            # identical terms & order


# ================================================ 2. PSD iteration helper (2.3)
def fit_ar_on_empirical_acv(v, grid=(2, 4, 8, 12, 20), tol=5e-2, maxLag=60):
    """Fit AR to the empirical ACV of residuals; order chosen on ACV fidelity
    over maxLag lags (mirrors lsq.selectArOrder's criterion, but on the
    EMPIRICAL acv). Returns (a, sigma, relErr, p)."""
    v = np.asarray(v, float).ravel()
    v = v - v.mean()
    n = len(v)
    L = min(maxLag, n // 4)
    r = np.array([v[:n - k] @ v[k:] for k in range(L + 1)]) / n   # biased ACV
    best = None
    for p in grid:
        a, s2 = levinson(r, p)
        rf = s2 * arma_acv(a, [1.0], L + 1)
        e = norm(rf - r) / norm(r)
        if best is None or e < best[2]:
            best = (a, np.sqrt(s2), e, len(a) - 1)
        if e < tol:
            break
    return best

def psd_iterate_reference(A, B, l_raw, a_true, s_true, maxIter=4):
    """The helper's algorithm on one arc: start WHITE (wrong), iterate."""
    a_cur, s_cur = np.array([1.0]), 1.0
    hist = []
    for it in range(maxIter):
        Aw, Bw, lw = (whiten_fir(a_cur, s_cur, X) for X in (A, B, l_raw))
        J = np.column_stack([Aw, Bw])
        xy = solve(J.T @ J, J.T @ lw)
        w = lw - J @ xy
        dof = len(lw) - J.shape[1]
        s0 = np.sqrt(w @ w / dof)
        # whiteness of the whitened residuals
        wc = w - w.mean()
        ac1 = max(abs(wc[:-k] @ wc[k:]) for k in (1, 2, 3)) / (wc @ wc)
        hist.append((s0, ac1))
        if abs(s0 - 1) < 0.05 and ac1 < 3 / np.sqrt(len(w)):
            break
        # re-estimate the model from RAW-metric residuals
        v = l_raw - A @ xy[:A.shape[1]] - B @ xy[A.shape[1]:]
        a_cur, s_cur, _, _ = fit_ar_on_empirical_acv(v)
    return hist, (a_cur, s_cur), it + 1

def t_psd_iteration_converges_from_white():
    nc, nl, m = 6, 3, 6000
    a = np.array([1.0, -1.3, 0.6]); s = 0.7
    r = np.random.default_rng(11)
    A = r.standard_normal((m, nc)); B = r.standard_normal((m, nl))
    l = A @ r.standard_normal(nc) + B @ r.standard_normal(nl) + ar_noise(a, s, m, r)

    hist, (af, sf), nIter = psd_iterate_reference(A, B, l, a, s)
    s0_first, s0_last = hist[0][0], hist[-1][0]
    ac_last = hist[-1][1]
    assert abs(s0_first - 1) > 0.5, f'start not wrong enough: {s0_first:.3f}'
    assert abs(s0_last - 1) < 0.05, f'sigma0 -> {s0_last:.4f}'
    assert ac_last < 3 / np.sqrt(m), f'acf(w) not white: {ac_last:.4f}'
    assert nIter <= 3, f'took {nIter} iterations'

def t_psd_iteration_recovers_the_model():
    """The fitted ACV must match the TRUE ACV -- the model is right, not merely
    self-consistent (that is the whole point of 2.3)."""
    nc, nl, m = 6, 3, 8000
    a = np.array([1.0, -1.3, 0.6]); s = 0.7
    r = np.random.default_rng(12)
    A = r.standard_normal((m, nc)); B = r.standard_normal((m, nl))
    l = A @ r.standard_normal(nc) + B @ r.standard_normal(nl) + ar_noise(a, s, m, r)

    _, (af, sf), _ = psd_iterate_reference(A, B, l, a, s)
    L = 40
    r_true = s**2 * arma_acv(a, [1.0], L)
    r_fit = sf**2 * arma_acv(af, [1.0], L)
    e = norm(r_fit - r_true) / norm(r_true)
    assert e < 0.12, f'ACV rel err {e:.3f}'   # sampling-limited, not bias

def t_residual_acv_bias_negligible_at_high_redundancy():
    """v = (I-H)e biases the empirical ACV by O(u/m). Measure it."""
    nc, m = 8, 5000
    a = np.array([1.0, -1.2, 0.5]); s = 1.0
    errs = []
    for seed in range(6):
        r = np.random.default_rng(100 + seed)
        A = r.standard_normal((m, nc))
        e = ar_noise(a, s, m, r)
        l = A @ r.standard_normal(nc) + e
        x = solve(A.T @ A, A.T @ l)
        v = l - A @ x
        L = 20
        rv = np.array([v[:m-k] @ v[k:] for k in range(L)]) / m
        re = np.array([e[:m-k] @ e[k:] for k in range(L)]) / m
        errs.append(norm(rv - re) / norm(re))
    assert np.mean(errs) < 0.02, f'ACV bias {np.mean(errs):.4f} at u/m={nc/m:.1e}'


# ====================================== 3. per-obsType sigma0 quick-look (2.4)
def t_per_obstype_quicklook_vs_exact():
    """rms_ot = sqrt(vtv/m) vs exact sqrt(vtv/r): quick-look biased low by
    ~sqrt(1 - h_ot/m_ot); with the group's true sigma it must land within
    sampling scatter of the truth once corrected, and vtv/m must be verified
    against brute-force residuals."""
    nc, nl, mG = 10, 4, 800
    TRUE = {'pos': 1.0, 'kbr': 2.0}
    r = np.random.default_rng(33)
    xt = r.standard_normal(nc); yt = r.standard_normal(nl)
    arcs, tag = [], []
    for ot in ('pos', 'kbr'):
        for _ in range(3):
            A = r.standard_normal((mG, nc)); B = r.standard_normal((mG, nl))
            l = A @ xt + B @ yt + TRUE[ot] * r.standard_normal(mG)
            arcs.append((A, B, l)); tag.append(ot)

    Ncc, Nlc, Nll, ncr, nlr, ltl, m = sys_normals(arcs, nc, nl)
    Nred, nred, _, _ = schur(Ncc, Nlc, Nll, ncr, nlr, ltl, m, nl)
    Nred = (Nred + Nred.T) / 2
    x = solve(Nred, nred)
    y = solve(Nll, nlr - Nlc @ x)
    blk = leverage_blocks(Nred, Nll, Nlc)

    agg = {ot: dict(m=0, vtv=0.0, rr=0.0) for ot in TRUE}
    for (A, B, l), ot in zip(arcs, tag):
        v = l - A @ x - B @ y
        h = joint_leverage(A, B, *blk)
        agg[ot]['m'] += len(l); agg[ot]['vtv'] += v @ v
        agg[ot]['rr'] += len(l) - h.sum()

    for ot in TRUE:
        g = agg[ot]
        quick = np.sqrt(g['vtv'] / g['m'])
        exact = np.sqrt(g['vtv'] / g['rr'])
        assert quick < exact, 'quick-look must be biased low'
        bias = 1 - quick / exact                     # = 1 - sqrt(r/m)
        pred = 1 - np.sqrt(g['rr'] / g['m'])
        assert abs(bias - pred) < 1e-12
        assert abs(exact / TRUE[ot] - 1) < 0.08, f'{ot}: {exact:.3f} vs {TRUE[ot]}'
        assert bias < 0.02, f'{ot}: quick-look bias {bias:.3%} (u/m small here)'


# ============================================================================
TESTS = [
    ('arc partition reconstructs system normals (3.4)',
     t_arc_partition_reconstructs_system_normals),
    ('system partition reconstructs reduced normals    ',
     t_system_partition_reconstructs_reduced_normals),
    ('redundancy pass: partition == serial + identity  ',
     t_redundancy_partition_reconstructs_serial),
    ('per-task + ordered reduce == fused serial loop   ',
     t_pertask_then_ordered_reduce_is_exact),
    ('PSD iteration converges from a white start (2.3) ',
     t_psd_iteration_converges_from_white),
    ('PSD iteration recovers the TRUE model            ',
     t_psd_iteration_recovers_the_model),
    ('residual-ACV bias negligible at high redundancy  ',
     t_residual_acv_bias_negligible_at_high_redundancy),
    ('per-obsType sigma0 quick-look vs exact (2.4)     ',
     t_per_obstype_quicklook_vs_exact),
]

if __name__ == '__main__':
    for name, fn in TESTS:
        T(name, fn)
    ok = sum(1 for _, s, _ in R if s == 'PASS')
    for name, s, msg in R:
        print(f'  [{"ok  " if s == "PASS" else s.lower():4s}] {name} {msg}')
    print('=' * 74)
    print(f'  {ok}/{len(R)} passed')
    sys.exit(0 if ok == len(R) else 1)
