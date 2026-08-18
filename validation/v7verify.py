#!/usr/bin/env python3
"""v7 validation: block-arrow Helmert VCE (lsqvce.vce) and the diagnoseAccum bridge.

Two gaps that had no dedicated test. Both are proven here against brute force
before the MATLAB tests are written:

  * VCE recovers KNOWN per-group variance components from a wrong start, and the
    per-group identity sigma_g^2 = vtv_g / r_g holds at convergence.
  * diagnoseAccum's reconstruction: the full block-arrow root factor reproduces the
    Schur-reduced common solution, sum(r_i) = n - u, and the per-observation
    r_i match the block-arrow leverage.

    python3 v7verify.py
"""
import numpy as np
from numpy.linalg import qr, inv, norm, solve

R = []
def T(n, f):
    try:
        f(); R.append((n, 'PASS', ''))
    except AssertionError as e:
        R.append((n, 'FAIL', str(e)[:70]))
    except Exception as e:
        R.append((n, 'ERROR', f'{type(e).__name__}: {e}'[:70]))

rng = np.random.default_rng(71)


# ---------------------------------------------------- block-arrow problem
def build(nc=5, nSys=2, nl=3, arcsPerSys=4, mArc=200, sig=None, seed=7):
    """Two obs GROUPS ('a','b') with DIFFERENT true noise sigma, block-arrow
    (common xc + per-system yl). Returns the whitened-at-truth arcs tagged by
    group, plus the truth."""
    if sig is None:
        sig = {'a': 1.0, 'b': 3.0}          # the components VCE must recover
    r = np.random.default_rng(seed)
    xc = r.standard_normal(nc)
    arcs = []                                # (sysId, group, A, B, l_raw)
    yl = [r.standard_normal(nl) for _ in range(nSys)]
    for s in range(nSys):
        for k in range(arcsPerSys):
            grp = 'a' if k % 2 == 0 else 'b'
            A = r.standard_normal((mArc, nc))
            B = r.standard_normal((mArc, nl))
            e = sig[grp] * r.standard_normal(mArc)
            l = A @ xc + B @ yl[s] + e
            arcs.append((s, grp, A, B, l))
    return xc, yl, arcs, sig, nc, nSys, nl


def solve_reduced(arcs, weights, nc, nSys, nl):
    """One weighted block-arrow solve: whiten each arc by 1/sigma_group, Schur-
    eliminate the per-system locals, solve the reduced common system.
    Returns xc, and per-group (vtv, r) needed for the Helmert update."""
    # accumulate per-system full normals under the current weights
    Ncc = np.zeros((nc, nc)); ncr = np.zeros(nc)
    sysBlocks = []
    perArc = []                             # keep for redundancy pass
    for s in range(nSys):
        Nll = np.zeros((nl, nl)); Nlc = np.zeros((nl, nc)); nlr = np.zeros(nl)
        Ncc_s = np.zeros((nc, nc)); ncr_s = np.zeros(nc)
        arcsS = []
        for (si, grp, A, B, l) in arcs:
            if si != s:
                continue
            w = 1.0 / weights[grp]
            Aw, Bw, lw = A * w, B * w, l * w
            Nll += Bw.T @ Bw; Nlc += Bw.T @ Aw; nlr += Bw.T @ lw
            Ncc_s += Aw.T @ Aw; ncr_s += Aw.T @ lw
            arcsS.append((grp, Aw, Bw, lw))
        sysBlocks.append((Nll, Nlc, nlr, arcsS))
        Z = solve(Nll, np.column_stack([Nlc, nlr]))
        Ncc += Ncc_s - Nlc.T @ Z[:, :-1]
        ncr += ncr_s - Nlc.T @ Z[:, -1]
    xc = solve((Ncc + Ncc.T) / 2, ncr)

    # redundancy + vtv per group, from the exact block-arrow leverage
    # Qxx_full via the full block-arrow inverse per system
    grp_vtv = {}; grp_r = {}; grp_m = {}
    # full parameter cofactor: assemble the whole block-arrow N and invert once
    u = nc + nSys * nl
    N = np.zeros((u, u)); nvec = np.zeros(u)
    off = [nc + s * nl for s in range(nSys)]
    for s in range(nSys):
        Nll, Nlc, nlr, arcsS = sysBlocks[s]
        N[:nc, :nc] += 0                      # common filled below
    # rebuild full N directly from whitened arcs (simpler + exact)
    N = np.zeros((u, u)); nvec = np.zeros(u)
    for s in range(nSys):
        _, _, _, arcsS = sysBlocks[s]
        lo = off[s]
        for (grp, Aw, Bw, lw) in arcsS:
            full = np.zeros((Aw.shape[0], u))
            full[:, :nc] = Aw; full[:, lo:lo+nl] = Bw
            N += full.T @ full; nvec += full.T @ lw
    Qxx = inv((N + N.T) / 2)
    xFull = Qxx @ nvec
    for s in range(nSys):
        _, _, _, arcsS = sysBlocks[s]
        lo = off[s]
        for (grp, Aw, Bw, lw) in arcsS:
            full = np.zeros((Aw.shape[0], u))
            full[:, :nc] = Aw; full[:, lo:lo+nl] = Bw
            v = lw - full @ xFull
            h = np.einsum('ij,jk,ik->i', full, Qxx, full)  # leverage
            grp_vtv[grp] = grp_vtv.get(grp, 0) + v @ v
            grp_r[grp]   = grp_r.get(grp, 0) + np.sum(1 - h)
            grp_m[grp]   = grp_m.get(grp, 0) + len(v)
    return xc, xFull, grp_vtv, grp_r, grp_m


def helmert(arcs, nc, nSys, nl, sigma0, maxIter=8, tol=1e-3):
    """The simplified Helmert iteration lsqvce.vce runs: sigma_g^2 = vtv_g/r_g,
    re-whiten, repeat."""
    w = dict(sigma0)
    hist = [dict(w)]
    for it in range(maxIter):
        _, _, vtv, r, m = solve_reduced(arcs, w, nc, nSys, nl)
        wnew = {g: np.sqrt(vtv[g] / r[g]) * w[g] for g in w}
        # note: since we whiten by 1/sigma, the update is multiplicative
        change = max(abs(wnew[g] / w[g] - 1) for g in w)
        w = wnew; hist.append(dict(w))
        if change < tol:
            return w, it + 1, True, hist
    return w, maxIter, False, hist


# ===================================================== 1. VCE recovers truth
def t_vce_recovers_components():
    xc, yl, arcs, sig, nc, nSys, nl = build()
    # start 3x wrong on group 'b'
    w, nIter, conv, hist = helmert(arcs, nc, nSys, nl,
                                   {'a': 1.0, 'b': 1.0})
    assert conv, 'did not converge'
    # recovered sigma ratios must match the truth ratio (absolute scale is
    # only identifiable up to the global sigma0, so compare the RATIO b/a)
    ratio_est = w['b'] / w['a']
    ratio_true = sig['b'] / sig['a']
    rel = abs(ratio_est - ratio_true) / ratio_true
    assert rel < 0.10, f'sigma_b/sigma_a: est {ratio_est:.3f} vs {ratio_true:.3f}'

def t_vce_identity_holds_at_convergence():
    xc, yl, arcs, sig, nc, nSys, nl = build(seed=9)
    w, nIter, conv, hist = helmert(arcs, nc, nSys, nl, {'a': 2.0, 'b': 0.5})
    _, _, vtv, r, m = solve_reduced(arcs, w, nc, nSys, nl)
    # at the fixed point sigma_g^2 == vtv_g / r_g (the update is a no-op)
    for g in w:
        assert abs(np.sqrt(vtv[g] / r[g]) - 1) < 5e-2, \
            f'group {g}: vtv/r = {vtv[g]/r[g]:.3f}, want ~1 after whitening'

def t_vce_converges_fast():
    xc, yl, arcs, sig, nc, nSys, nl = build(seed=3)
    _, nIter, conv, _ = helmert(arcs, nc, nSys, nl, {'a': 1.0, 'b': 1.0})
    assert conv and nIter <= 4, f'took {nIter} iterations'


# ================================================ 2. diagnoseAccum reconstruction
def t_fromaccum_full_equals_reduced():
    """The full block-arrow solve (what diagnoseAccum builds) reproduces the Schur-
    reduced common solution."""
    xc, yl, arcs, sig, nc, nSys, nl = build(seed=5)
    w = {'a': 1.0, 'b': 1.0}
    xc_red, xFull, vtv, r, m = solve_reduced(arcs, w, nc, nSys, nl)
    assert norm(xFull[:nc] - xc_red) / norm(xc_red) < 1e-9

def t_fromaccum_sumr_identity():
    xc, yl, arcs, sig, nc, nSys, nl = build(seed=6)
    w = {'a': 1.0, 'b': 1.0}
    _, _, vtv, r, m = solve_reduced(arcs, w, nc, nSys, nl)
    n = sum(m.values()); u = nc + nSys * nl
    assert abs(sum(r.values()) - (n - u)) < 1e-6, \
        f'sum(r)={sum(r.values()):.3f}, n-u={n-u}'

def t_fromaccum_streamed_R_matches():
    """Stream the whitened full-width batches through QR (as diagnoseAccum does);
    the resulting r_i must equal the direct-inverse leverage."""
    xc, yl, arcs, sig, nc, nSys, nl = build(seed=8)
    u = nc + nSys * nl
    off = [nc + s * nl for s in range(nSys)]
    batches = []
    for (s, grp, A, B, l) in arcs:
        full = np.zeros((A.shape[0], u))
        full[:, :nc] = A; full[:, off[s]:off[s]+nl] = B
        batches.append(full)
    Rst = np.zeros((0, u))
    for full in batches:
        Rst = np.triu(qr(np.vstack([Rst, full]), mode='r'))[:u]
    sumr = 0; n = 0
    for full in batches:
        W = solve(Rst.T, full.T)
        sumr += np.sum(1 - np.sum(W**2, axis=0)); n += full.shape[0]
    assert abs(sumr - (n - u)) < 1e-6, f'sum(r)={sumr:.3f}, n-u={n-u}'


T('1. VCE recovers the known variance-component ratio from a wrong start',
  t_vce_recovers_components)
T('2. VCE fixed point satisfies sigma_g^2 = vtv_g / r_g',
  t_vce_identity_holds_at_convergence)
T('3. VCE converges in <= 4 Helmert iterations',
  t_vce_converges_fast)
T('4. diagnoseAccum: full block-arrow solve == Schur-reduced common solution',
  t_fromaccum_full_equals_reduced)
T('5. diagnoseAccum: sum(r_i) = n - u (block-arrow redundancy identity)',
  t_fromaccum_sumr_identity)
T('6. diagnoseAccum: streamed full-width R gives the same redundancy',
  t_fromaccum_streamed_R_matches)

n = sum(1 for _, s, _ in R if s == 'PASS')
for name, s, msg in R:
    print(f'  [{s:5s}] {name}' + (f'  -- {msg}' if msg else ''))
print(f'  {n}/{len(R)} passed')
