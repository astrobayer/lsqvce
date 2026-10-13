#!/usr/bin/env python3
"""WP-132 mirror: the joint estimator's arc factor with the WHOLE law as its
component (lsqvce.kernel.redundancyLags .omegaLaw / .sLaw, review low 23).

The arc factor (eq. 6.5.27) is Helmert's ratio q / t with the arc's whole law
covariance T as the component. sweepArcs formed it from the lag sums,
q = r' Omega, t = r' s with r = T(:, 1): exact for a Toeplitz T, but for any
other T it is the ratio for a DIFFERENT component (toeplitz(r)) -- still
unbiased (E[Omega_V] = sigma^2 s_V for every V), not efficient.

Closed form (the probe realisation redundancyLags uses: Rt = W^-1 Abar UiZ,
K = W^-1 Bbar Ll^-T, Rbar = Sigma^-1 - Rt Rt' - K K', etilde = W^-1 ebar):
  no nugget (T = Sigma = W'W):  q = ebar'ebar,  t = N - |Abar UiZ|^2 - |Bbar Ll^-T|^2
  nugget (Sigma = T + diag(cN), T = WL'WL):
      q = |WL etilde|^2,  t = N - sum cN_i (Sigma^-1)_ii - |WL Rt|^2 - |WL K|^2
Checked against the explicit etilde' T etilde and trace(Rbar T), and against
r'Omega / r's on a Toeplitz T; then the efficiency on a per-row sigma 1..10.

Run:  python3 dev/validation/wp132_arc_component_mirror.py
Created 2026-10-13 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-132.
"""
import numpy as np


def lagsums(M):
    """s_0 = trace M, s_k = 2 sum of the k-th sub-diagonal (redundancyLags' V_n)."""
    return np.array([np.trace(M) if k == 0 else 2 * np.trace(M, -k) for k in range(M.shape[0])])


def closed_form(Abar, Bbar, ebar, W, UiZ, Ll, cN=None, WL=None):
    """lsqvce.kernel.redundancyLags .omegaLaw / .sLaw, line for line."""
    N = ebar.size
    P = Abar @ UiZ
    Kb = np.zeros((N, 0)) if Bbar.shape[1] == 0 else np.linalg.solve(Ll, Bbar.T).T   # Bbar Ll^-T
    if cN is None:
        return ebar @ ebar, N - np.sum(P ** 2) - np.sum(Kb ** 2)
    et = WL @ np.linalg.solve(W, ebar)
    rej = np.flatnonzero(cN > 0)
    dSi = np.array([np.sum(np.linalg.solve(W.T, np.eye(N)[:, i]) ** 2) for i in rej])
    Pr = WL @ np.linalg.solve(W, P)
    Kr = WL @ np.linalg.solve(W, Kb)
    return et @ et, N - np.sum(cN[rej] * dSi) - np.sum(Pr ** 2) - np.sum(Kr ** 2)


def explicit(Abar, Bbar, ebar, W, UiZ, Ll, T):
    """etilde' T etilde and trace(Rbar T) from the dense Rbar."""
    Wi = np.linalg.inv(W)
    Rt = Wi @ Abar @ UiZ
    K = Wi @ np.linalg.solve(Ll, Bbar.T).T
    Rbar = Wi @ Wi.T - Rt @ Rt.T - K @ K.T
    et = Wi @ ebar
    return et @ T @ et, np.trace(Rbar @ T), Rbar, et


def main():
    rng = np.random.default_rng(132)
    N, P, nl, Z = 50, 4, 2, 100
    t = np.arange(N)
    A = rng.standard_normal((N, P))
    B = np.c_[np.ones(N), t / N]
    UiZ = rng.standard_normal((P, Z)) / np.sqrt(Z)        # any probe block: the identity is algebraic
    #                                                      (not a projector here, so t may be < 0)
    Ll = np.linalg.cholesky(B.T @ B + np.eye(nl))
    e = rng.standard_normal(N)
    worst = 0.0
    # 1) a Toeplitz law: closed form == explicit == r'Omega / r's
    T = 2.0 * 0.8 ** np.abs(t[:, None] - t[None, :])
    W = np.linalg.cholesky(T).T
    Abar, Bbar, ebar = (np.linalg.solve(W.T, M) for M in (A, B, e))
    q, s = closed_form(Abar, Bbar, ebar, W, UiZ, Ll)
    qe, se, Rbar, et = explicit(Abar, Bbar, ebar, W, UiZ, Ll, T)
    r = T[:, 0]
    qr, sr = r @ lagsums(np.outer(et, et)), r @ lagsums(Rbar)
    print(f'Toeplitz T:     q {q:.12g} | explicit {qe:.12g} | r\'Omega {qr:.12g}')
    print(f'                t {s:.12g} | explicit {se:.12g} | r\'s     {sr:.12g}')
    worst = max(worst, abs(q / qe - 1), abs(s / se - 1), abs(qr / qe - 1), abs(sr / se - 1))
    # 2) a non-Toeplitz law (per-row sigma): closed form == explicit, r'Omega is not
    sig = np.linspace(1, 3, N)
    Tn = np.outer(sig, sig) * 0.6 ** np.abs(t[:, None] - t[None, :])
    W = np.linalg.cholesky(Tn).T
    Abar, Bbar, ebar = (np.linalg.solve(W.T, M) for M in (A, B, e))
    q, s = closed_form(Abar, Bbar, ebar, W, UiZ, Ll)
    qe, se, Rbar, et = explicit(Abar, Bbar, ebar, W, UiZ, Ll, Tn)
    qr = Tn[:, 0] @ lagsums(np.outer(et, et))
    print(f'non-Toeplitz T: q {q:.12g} | explicit {qe:.12g} | r\'Omega {qr:.6g} (a different component)')
    worst = max(worst, abs(q / qe - 1), abs(s / se - 1))
    # 3) a nugget on three rows (the robust pass's rejected epochs): T is the component
    cN = np.zeros(N)
    cN[[5, 17, 40]] = [1e3, 50.0, 7.0]
    Sig = Tn + np.diag(cN)
    W, WL = np.linalg.cholesky(Sig).T, np.linalg.cholesky(Tn).T
    Abar, Bbar, ebar = (np.linalg.solve(W.T, M) for M in (A, B, e))
    q, s = closed_form(Abar, Bbar, ebar, W, UiZ, Ll, cN, WL)
    qe, se, _, _ = explicit(Abar, Bbar, ebar, W, UiZ, Ll, Tn)
    print(f'nugget:         q {q:.12g} | explicit {qe:.12g};  t {s:.12g} | explicit {se:.12g}')
    worst = max(worst, abs(q / qe - 1), abs(s / se - 1))
    print(f'worst relative disagreement with the explicit forms: {worst:.1e}')
    # 4) efficiency: white noise, per-row sigma 1..10 over 60 epochs, factor 1.7
    m, p = 60, 4
    tt = np.arange(m)
    Ad = np.c_[np.ones(m), tt / m, np.sin(tt / 7), np.cos(tt / 7)]
    sg = np.linspace(1, 10, m)
    Td = np.diag(sg ** 2)
    W = np.diag(sg)
    Abar = Ad / sg[:, None]
    L = np.linalg.cholesky(Abar.T @ Abar)
    UiZx = np.linalg.inv(L.T)                              # exact projector (Z -> inf)
    fo, fn = [], []
    for _ in range(4000):
        l = Ad @ rng.standard_normal(p) + np.sqrt(1.7) * sg * rng.standard_normal(m)
        lb = l / sg
        eb = lb - Abar @ np.linalg.lstsq(Abar, lb, rcond=None)[0]
        q, s = closed_form(Abar, np.zeros((m, 0)), eb, W, UiZx, None)
        _, _, Rb, et = explicit(Abar, np.zeros((m, 0)), eb, W, UiZx, np.eye(0), Td)
        r = Td[:, 0]
        fo.append((r @ lagsums(np.outer(et, et))) / (r @ lagsums(Rb)))
        fn.append(q / s)
    fo, fn = np.array(fo), np.array(fn)
    bique = 1.7 * np.sqrt(2 / (m - p))
    print(f'per-row sigma 1..10, factor 1.7: r\'Omega mean {fo.mean():.3f} sd {fo.std():.3f} | '
          f'whole law mean {fn.mean():.3f} sd {fn.std():.3f} (BIQUE sd {bique:.3f})')
    assert worst < 1e-10
    assert abs(fn.std() / bique - 1) < 0.1 and fo.std() > 1.5 * fn.std()
    assert abs(fo.mean() / 1.7 - 1) < 0.05 and abs(fn.mean() / 1.7 - 1) < 0.05


if __name__ == '__main__':
    main()
